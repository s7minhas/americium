# decompose the us shift

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))
source(paste0(path_funcs, 'project_metadata.R'))

load_pkg('data.table')

load(paste0(path_in, 'un_votes_long.rda'))
setDT(un_votes_long)

pos_map = c(Y = 1, A = 0, N = -1)

# pairwise lijphart agreement
agree = function(yr, a, b){
	d = un_votes_long[year == yr & ms_code %in% c(a, b) & ms_vote %in% c('Y', 'N', 'A')]
	w = dcast(d, resolution ~ ms_code, value.var = 'ms_vote')
	if(!all(c(a, b) %in% names(w))){ return(NA_real_) }
	w = w[!is.na(get(a)) & !is.na(get(b))]
	if(nrow(w) == 0){ return(NA_real_) }
	1 - mean(abs(pos_map[w[[a]]] - pos_map[w[[b]]]) / 2)
}

allies = amer_actor_cats$allies
rivals = amer_actor_cats$rivals

# us voting profile by year
us_profile = un_votes_long[ms_code == 'USA' & ms_vote %in% c('Y', 'N', 'A'),
	.(n = .N, yes = sum(ms_vote == 'Y'), no = sum(ms_vote == 'N'),
		abstain = sum(ms_vote == 'A')), by = year][order(year)]

# us agreement by group
years = 2018:2025
us_alignment = data.table(
	year = years,
	us_allies = sapply(years, function(y) mean(sapply(allies, function(x) agree(y, 'USA', x)), na.rm = TRUE)),
	us_rivals = sapply(years, function(y) mean(sapply(rivals, function(x) agree(y, 'USA', x)), na.rm = TRUE))
)

us_decomposition = list(
	us_profile = us_profile,
	us_alignment = us_alignment
)
save(us_decomposition, file = paste0(path_out, 'us_decomposition.rda'))
