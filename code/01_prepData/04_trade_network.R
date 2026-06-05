# build trade networks

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))

load_pkg('data.table')

# keep reported positive flows since 1989
load(paste0(path_in, 'trade_long.rda'))
setDT(trade_long)
trade = trade_long[year >= 1989 & exports_usd > 0]

# actor-years
master_actors = sort(unique(c(trade$reporter, trade$partner)))
all_years = sort(unique(trade$year))

# build yearly trade grids
trade_exports = vector('list', length(all_years))
trade_total = vector('list', length(all_years))
names(trade_exports) = all_years
names(trade_total) = all_years
for(yr in all_years){
	year_trade = trade[year == yr]
	exp_mat = matrix(
		NA_real_, length(master_actors), length(master_actors),
		dimnames = list(master_actors, master_actors)
	)
	exp_mat[cbind(year_trade$reporter, year_trade$partner)] = year_trade$exports_usd
	diag(exp_mat) = NA

	# combine both export directions
	reported = !is.na(exp_mat) | !is.na(t(exp_mat))
	exp_zero = exp_mat
	exp_zero[is.na(exp_zero)] = 0
	tot_mat = exp_zero + t(exp_zero)
	tot_mat[!reported] = NA
	diag(tot_mat) = NA

	trade_exports[[as.character(yr)]] = exp_mat
	trade_total[[as.character(yr)]] = tot_mat
}

trade_networks = list(
	exports_usd = trade_exports,
	total_usd = trade_total,
	actors = master_actors,
	years = all_years
)

save(trade_networks, file = paste0(path_out, 'trade_networks.rda'))
