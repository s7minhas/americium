# build un voting networks

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))

load_pkg(c('data.table'))

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

# align one year
place_on_master = function(mat, actors){
	out = matrix(NA_real_, length(actors), length(actors), dimnames = list(actors, actors))
	out[rownames(mat), colnames(mat)] = mat
	out
}

# keep recorded positions since 1989
load(paste0(path_in, 'un_votes_long.rda'))
setDT(un_votes_long)
votes = un_votes_long[year >= 1989 & ms_vote %in% c('Y', 'N', 'A')]

# encode vote positions
pos_map = c(Y = 1, A = 0, N = -1)
votes[, pos := pos_map[ms_vote]]

# resolution salience
res_info = unique(votes[, .(
	resolution, year, date, title, agenda_title, subjects,
	total_yes, total_no, total_abstentions, total_non_voting, total_ms
)])
res_info[, present := total_yes + total_no + total_abstentions]
res_info[, turnout := present / total_ms]
res_info[, competitiveness := fifelse(
	total_yes + total_no > 0,
	1 - abs(total_yes - total_no) / (total_yes + total_no),
	0
)]
res_info[, is_emergency := grepl('ES-', resolution)]

# weight contested votes
weight_floor = 0.1
res_info[, weight := weight_floor + (1 - weight_floor) * competitiveness]

# actor-years
master_actors = sort(unique(votes$ms_code))
all_years = sort(unique(votes$year))

# build yearly networks
un_agree_w = vector('list', length(all_years))
un_agree_u = vector('list', length(all_years))
names(un_agree_w) = all_years
names(un_agree_u) = all_years
for(yr in all_years){
	year_votes = votes[year == yr]
	pos_wide = dcast(year_votes, ms_code ~ resolution, value.var = 'pos')
	pos_mat = as.matrix(pos_wide[, -1])
	rownames(pos_mat) = pos_wide$ms_code
	vote_weight = res_info[match(colnames(pos_mat), resolution), weight]
	un_agree_w[[as.character(yr)]] = place_on_master(
		build_agreement(pos_mat, vote_weight), master_actors
	)
	un_agree_u[[as.character(yr)]] = place_on_master(
		build_agreement(pos_mat, rep(1, ncol(pos_mat))), master_actors
	)
}

un_networks = list(
	weighted = un_agree_w,
	unweighted = un_agree_u,
	actors = master_actors,
	years = all_years,
	vote_info = res_info[]
)

save(un_networks, file = paste0(path_out, 'un_networks.rda'))
