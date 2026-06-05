# 2025 held-out surprise

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))

load_pkg('lame')
set.seed(seed)

load(paste0(path_in, 'arr_list_lame.rda'))
un = arr_list_lame$un_voting
years = as.integer(names(un))

# train on 1989-2024
train = as.character(years[years <= 2024])
fit = lame(
	un[train], R = 2, family = 'normal', symmetric = TRUE,
	dynamic_uv = TRUE, dynamic_uv_kind = 'ar1',
	nscan = 4000, burn = 2000, odens = 10, seed = seed, gof = FALSE, verbose = FALSE
)

u = fit$U
fit_years = dimnames(u)[[3]]
actors = rownames(u[, , 1])
add_eff = setNames(as.numeric(fit$APM), actors)
rho = mean(fit$rho_uv)
s2 = mean(fit$VC[, 've'])

# one-step drift forecast
forecast_eta = function(u_prev){
	u_hat = rho * u_prev
	outer(add_eff, add_eff, '+') + u_hat %*% t(u_hat)
}

# per-state log density
state_logdens = function(obs, eta){
	ld = dnorm(obs, mean = eta, sd = sqrt(s2), log = TRUE)
	diag(ld) = NA
	rowMeans(ld, na.rm = TRUE)
}

# training-period baseline
null_ld = matrix(NA_real_, length(actors), length(fit_years) - 1,
	dimnames = list(actors, fit_years[-1]))
for(t in 2:length(fit_years)){
	eta = forecast_eta(u[, , t-1])
	obs = un[[fit_years[t]]][actors, actors]
	null_ld[, t-1] = state_logdens(obs, eta)
}

# 2025 forecast score
eta_2025 = forecast_eta(u[, , length(fit_years)])
obs_2025 = un[['2025']][actors, actors]
ld_2025 = state_logdens(obs_2025, eta_2025)

# per-state 2025 surprise
state_surprise = vapply(actors, function(i){
	hist = null_ld[i, ]
	hist = hist[is.finite(hist)]
	if(length(hist) < 5 || !is.finite(ld_2025[i])){ return(NA_real_) }
	(median(hist) - ld_2025[i]) / (mad(hist) + 1e-8)
}, numeric(1))

oos_surprise = list(
	state_surprise_2025 = state_surprise,
	ld_2025 = ld_2025,
	null_ld = null_ld,
	rho = rho,
	s2 = s2
)
save(oos_surprise, file = paste0(path_out, 'oos_surprise.rda'))
