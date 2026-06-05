# build sanctions networks

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))

load_pkg('data.table')

# keep dated designations since 1989
load(paste0(path_in, 'sanctions_long.rda'))
setDT(sanctions_long)
sanc = sanctions_long[!is.na(year) & year >= 1989]

# count new designations
flow = sanc[, .(n = .N), by = .(sender, target, year)]

# actor-years
master_actors = sort(unique(c(flow$sender, flow$target)))
all_years = sort(unique(flow$year))

# build yearly sanctions grids
sanctions_new = vector('list', length(all_years))
names(sanctions_new) = all_years
for(yr in all_years){
	year_flow = flow[year == yr]
	mat = matrix(
		0, length(master_actors), length(master_actors),
		dimnames = list(master_actors, master_actors)
	)
	mat[cbind(year_flow$sender, year_flow$target)] = year_flow$n
	diag(mat) = NA
	sanctions_new[[as.character(yr)]] = mat
}

sanctions_networks = list(
	new_designations = sanctions_new,
	actors = master_actors,
	years = all_years
)

save(sanctions_networks, file = paste0(path_out, 'sanctions_networks.rda'))
