# 2025 snap correlates

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))

load_pkg(c('vdemdata', 'peacesciencer', 'countrycode'))

# outcome
load(paste0(path_out, 'snap_fit_un_voting.rda'))
snap_2025 = fit$snap_prob[, '2025']
actors = names(snap_2025)

# regime type
vd = vdem[vdem$year == 2024, c('COWcode', 'v2x_polyarchy')]
vd$iso3 = suppressWarnings(countrycode(vd$COWcode, 'cown', 'iso3c'))
polyarchy = setNames(vd$v2x_polyarchy, vd$iso3)

# defense-pact ally
us_dyads = create_dyadyears(system = 'cow', directed = TRUE)
us_dyads = us_dyads[us_dyads$year == 2018 & us_dyads$ccode1 == 2, ]
us_dyads = add_atop_alliance(us_dyads)
us_dyads$iso3 = suppressWarnings(countrycode(us_dyads$ccode2, 'cown', 'iso3c'))
us_ally = setNames(us_dyads$atop_defense, us_dyads$iso3)

# exposure covariates
load(paste0(path_in, 'arr_list_lame.rda'))
un_2024 = arr_list_lame$un_voting[['2024']]
trade_2024 = arr_list_lame$trade[['2024']]
sanctions_pre = Reduce('+', lapply(as.character(1989:2024), function(y){
	m = arr_list_lame$sanctions[[y]]
	m[is.na(m)] = 0
	colSums(m)
}))
trade_total = rowSums(trade_2024, na.rm = TRUE)

covars = data.frame(
	actor = actors,
	snap = snap_2025[actors],
	us_alignment = un_2024[actors, 'USA'],
	rival_alignment = rowMeans(un_2024[actors, c('RUS', 'CHN')], na.rm = TRUE),
	us_trade_dep = trade_2024[actors, 'USA'] / trade_total[actors],
	chn_trade_dep = trade_2024[actors, 'CHN'] / trade_total[actors],
	trade_centrality = log1p(trade_total[actors]),
	sanctions_exposure = log1p(sanctions_pre[actors]),
	polyarchy = polyarchy[actors],
	us_ally = us_ally[actors],
	stringsAsFactors = FALSE
)
covars$us_ally[is.na(covars$us_ally)] = 0
covars = covars[is.finite(covars$snap) & complete.cases(covars), ]

# fractional logit
covars$snap_b = pmin(pmax(covars$snap, 1e-4), 1 - 1e-4)
continuous_covars = c(
	'us_alignment', 'rival_alignment', 'us_trade_dep', 'chn_trade_dep',
	'trade_centrality', 'sanctions_exposure', 'polyarchy'
)
for(v in continuous_covars){
	covars[[paste0(v, '_s')]] = (covars[[v]] - mean(covars[[v]], na.rm = TRUE)) /
		(2 * sd(covars[[v]], na.rm = TRUE))
}

who_model = glm(
	snap_b ~ us_alignment_s + rival_alignment_s + us_trade_dep_s +
		chn_trade_dep_s + trade_centrality_s + sanctions_exposure_s +
		polyarchy_s + us_ally,
	data = covars,
	family = quasibinomial()
)
coef_table = as.data.frame(summary(who_model)$coefficients)
coef_table$term = rownames(coef_table)
coef_table = coef_table[, c('term', 'Estimate', 'Std. Error', 't value', 'Pr(>|t|)')]
rownames(coef_table) = NULL

who_snaps_full = list(
	coefficients = summary(who_model)$coefficients,
	coefficient_table = coef_table,
	n = nrow(covars),
	mean_snap_by_ally = tapply(covars$snap, covars$us_ally, mean),
	data = covars
)
save(who_snaps_full, file = paste0(path_out, 'who_snaps_full.rda'))
write.csv(coef_table, paste0(path_out, 'who_snaps_full_coefficients.csv'), row.names = FALSE)
