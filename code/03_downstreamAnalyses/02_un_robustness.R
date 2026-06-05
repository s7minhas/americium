# un voting robustness checks

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))
source(paste0(path_funcs, 'project_metadata.R'))

load_pkg(c('data.table', 'lame'))
set.seed(seed)

# weighted lijphart agreement
build_agreement = function(pos_mat, vote_weight){
	actors = rownames(pos_mat)
	num = matrix(0, length(actors), length(actors), dimnames = list(actors, actors))
	den = matrix(0, length(actors), length(actors), dimnames = list(actors, actors))
	for(v in seq_len(ncol(pos_mat))){
		present = which(!is.na(pos_mat[, v]))
		if(length(present) < 2){ next }
		positions = pos_mat[present, v]
		agree = 1 - abs(outer(positions, positions, '-')) / 2
		num[present, present] = num[present, present] + vote_weight[v] * agree
		den[present, present] = den[present, present] + vote_weight[v]
	}
	out = num / den
	out[den == 0] = NA
	diag(out) = NA_real_
	out
}
place_on_master = function(mat, actors){
	out = matrix(NA_real_, length(actors), length(actors), dimnames = list(actors, actors))
	out[rownames(mat), colnames(mat)] = mat
	out
}

# recorded positions since 1989
load(paste0(path_in, 'un_votes_long.rda'))
setDT(un_votes_long)
votes = un_votes_long[year >= 1989 & ms_vote %in% c('Y', 'N', 'A')]
pos_map = c(Y = 1, A = 0, N = -1)
votes[, pos := pos_map[ms_vote]]

# resolution competitiveness
res_info = unique(votes[, .(resolution, year, total_yes, total_no)])
res_info[, competitiveness := fifelse(
	total_yes + total_no > 0,
	1 - abs(total_yes - total_no) / (total_yes + total_no),
	0
)]
weight_floor = 0.1
res_info[, weight := weight_floor + (1 - weight_floor) * competitiveness]

# drop top-ten-percent divisive votes
res_info[, drop_divisive := competitiveness >= quantile(competitiveness, 0.9, na.rm = TRUE), by = year]
drop_set = res_info[drop_divisive == TRUE, resolution]

master_actors = sort(unique(votes$ms_code))
all_years = sort(unique(votes$year))

# build unweighted agreement networks
build_year_list = function(keep_votes, weighted = FALSE){
	out = vector('list', length(all_years))
	names(out) = all_years
	for(yr in all_years){
		year_votes = keep_votes[year == yr]
		pos_wide = dcast(year_votes, ms_code ~ resolution, value.var = 'pos')
		pos_mat = as.matrix(pos_wide[, -1])
		rownames(pos_mat) = pos_wide$ms_code
		vote_weight = if(weighted) {
			res_info[match(colnames(pos_mat), resolution), weight]
		} else {
			rep(1, ncol(pos_mat))
		}
		out[[as.character(yr)]] = place_on_master(
			build_agreement(pos_mat, vote_weight), master_actors
		)
	}
	out
}

variants = list(
	unweighted = build_year_list(votes, weighted = FALSE),
	drop_divisive_weighted = build_year_list(votes[!resolution %in% drop_set], weighted = TRUE),
	drop_divisive_unweighted = build_year_list(votes[!resolution %in% drop_set], weighted = FALSE)
)

# refit each variant
allies = amer_actor_cats$allies
rivals = amer_actor_cats$rivals
hedgers = amer_actor_cats$hedgers
order_actors = amer_focal_actors
un_robustness = list()
for(v in names(variants)){
	fit = lame(
		variants[[v]], R = 2, family = 'normal', symmetric = TRUE,
		dynamic_uv = TRUE, dynamic_uv_kind = 'snap',
		prior = amer_snap_prior,
		nscan = 3000, burn = 1500, odens = 10, keep_snap_draws = 'draws',
		seed = seed, gof = FALSE, verbose = FALSE
	)
	sp = fit$snap_prob
	idx = amer_snap_index(sp, order_actors)
	un_robustness[[v]] = list(
		us_2025 = sp['USA', '2025'],
		allies_2025 = mean(sp[intersect(allies, rownames(sp)), '2025'], na.rm = TRUE),
		rivals_2025 = mean(sp[intersect(rivals, rownames(sp)), '2025'], na.rm = TRUE),
		hedgers_2025 = mean(sp[intersect(hedgers, rownames(sp)), '2025'], na.rm = TRUE),
		index_2025 = idx['2025'],
		rank_2025 = amer_rank_year(idx, '2025'),
		top_2025 = sort(sp[, '2025'], decreasing = TRUE)[1:8]
	)
}

un_robustness$americium_build = amer_build_metadata()

save(un_robustness, file = paste0(path_out, 'un_robustness.rda'))
