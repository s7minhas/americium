# align model arrays

rm(list = ls())
setup_candidates = c(
	file.path(here::here(), 'setup.R'),
	file.path(here::here(), 'code', 'setup.R'),
	file.path(getwd(), 'setup.R'),
	file.path(getwd(), 'code', 'setup.R'),
	file.path(dirname(getwd()), 'code', 'setup.R')
)
setup_candidates = setup_candidates[file.exists(setup_candidates)]
if(length(setup_candidates) == 0){ stop('Could not locate setup.R.') }
source(setup_candidates[1])

load(paste0(path_out, 'un_networks.rda'))
load(paste0(path_out, 'trade_networks.rda'))
load(paste0(path_out, 'sanctions_networks.rda'))

# align networks to the actor grid
align_relation = function(net_list, actors, years, fill = NA_real_){
	out = vector('list', length(years))
	names(out) = as.character(years)
	for(yr in as.character(years)){
		grid = matrix(
			fill, length(actors), length(actors),
			dimnames = list(actors, actors)
		)
		mat = net_list[[yr]]
		if(!is.null(mat)){
			shared_r = intersect(rownames(mat), actors)
			shared_c = intersect(colnames(mat), actors)
			grid[shared_r, shared_c] = mat[shared_r, shared_c]
		}
		diag(grid) = NA_real_
		out[[yr]] = grid
	}
	out
}

# use un actors as the grid
system_actors = sort(un_networks$actors)
common_years = sort(Reduce(
	intersect,
	list(un_networks$years, trade_networks$years, sanctions_networks$years)
))

# align each relation
arr_list_lame = list(
	un_voting = align_relation(un_networks$weighted, system_actors, common_years),
	trade = align_relation(trade_networks$total_usd, system_actors, common_years),
	sanctions = align_relation(sanctions_networks$new_designations, system_actors, common_years, fill = 0),
	actors = system_actors,
	years = common_years
)

save(arr_list_lame, file = paste0(path_in, 'arr_list_lame.rda'))
