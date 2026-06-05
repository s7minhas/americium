# 2025 snap regression

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))

# outcome
load(paste0(path_out, 'snap_fit_un_voting.rda'))
snap_2025 = fit$snap_prob[, '2025']

# exposure covariates
load(paste0(path_in, 'arr_list_lame.rda'))
actors = arr_list_lame$actors
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
	us_alignment = un_2024[, 'USA'],
	rival_alignment = rowMeans(un_2024[, c('RUS','CHN')], na.rm = TRUE),
	us_trade_dep = trade_2024[, 'USA'] / trade_total,
	chn_trade_dep = trade_2024[, 'CHN'] / trade_total,
	trade_centrality = log1p(trade_total),
	sanctions_exposure = log1p(sanctions_pre[actors]),
	stringsAsFactors = FALSE
)
covars = covars[is.finite(covars$snap) & complete.cases(covars), ]

# fractional logit
covars$snap_b = pmin(pmax(covars$snap, 1e-4), 1 - 1e-4)
who_model = glm(
	snap_b ~ us_alignment + rival_alignment + us_trade_dep + chn_trade_dep +
		trade_centrality + sanctions_exposure,
	data = covars,
	family = quasibinomial()
)

who_snaps = list(
	coefficients = summary(who_model)$coefficients,
	n = nrow(covars),
	data = covars
)
save(who_snaps, file = paste0(path_out, 'who_snaps.rda'))
