# simulate planted snaps

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
source(paste0(path_funcs, 'project_metadata.R'))

load_pkg('lame')
set.seed(seed)

# simulate latent paths
n = 30
n_time = 8
rank = 2
rho_true = 0.9
sigma_true = 0.2
kappa_true = 2

# plant snap years
snap_events = list(c(2, 3), c(5, 5), c(9, 6), c(14, 4), c(20, 7), c(24, 8))
delta_true = matrix(0, n, n_time)
for(ev in snap_events){
	delta_true[ev[1], ev[2]] = 1
}

u_true = array(0, dim = c(n, rank, n_time))
u_true[, , 1] = matrix(rnorm(n * rank), n, rank)
for(t in 2:n_time){
	for(i in 1:n){
		if(delta_true[i, t] == 1){
			u_true[i, , t] = rnorm(rank, 0, kappa_true)
		} else {
			u_true[i, , t] = rho_true * u_true[i, , t-1] + rnorm(rank, 0, sigma_true)
		}
	}
}

# build symmetric networks
y_list = lapply(1:n_time, function(t){
	mat = u_true[, , t] %*% t(u_true[, , t]) + matrix(rnorm(n * n, 0, 0.3), n, n)
	mat = (mat + t(mat)) / 2
	diag(mat) = NA
	rownames(mat) = colnames(mat) = paste0('a', 1:n)
	mat
})

# fit the snap model
fit = lame(
	y_list,
	R = rank,
	family = 'normal',
	symmetric = TRUE,
	dynamic_uv = TRUE,
	dynamic_uv_kind = 'snap',
	prior = c(list(snap_kappa = kappa_true), amer_snap_prior),
	nscan = 4000,
	burn = 2000,
	odens = 10,
	keep_snap_draws = 'draws',
	seed = seed,
	gof = FALSE,
	verbose = FALSE
)

# align actors by name
rownames(delta_true) = paste0('a', 1:n)
snap_prob = fit$snap_prob[rownames(delta_true), ]

# rank snap recovery
eligible = col(snap_prob) > 1
pos = snap_prob[eligible & delta_true == 1]
neg = snap_prob[eligible & delta_true == 0]
auc = mean(outer(pos, neg, '>')) + 0.5 * mean(outer(pos, neg, '=='))

# check drift-only false snaps
u_drift = array(0, dim = c(n, rank, n_time))
u_drift[, , 1] = matrix(rnorm(n * rank), n, rank)
for(t in 2:n_time){
	u_drift[, , t] = rho_true * u_drift[, , t-1] + matrix(rnorm(n * rank, 0, sigma_true), n, rank)
}
y_drift = lapply(1:n_time, function(t){
	mat = u_drift[, , t] %*% t(u_drift[, , t]) + matrix(rnorm(n * n, 0, 0.3), n, n)
	mat = (mat + t(mat)) / 2
	diag(mat) = NA
	rownames(mat) = colnames(mat) = paste0('a', 1:n)
	mat
})
fit_drift = lame(
	y_drift,
	R = rank,
	family = 'normal',
	symmetric = TRUE,
	dynamic_uv = TRUE,
	dynamic_uv_kind = 'snap',
	prior = c(list(snap_kappa = kappa_true), amer_snap_prior),
	nscan = 4000,
	burn = 2000,
	odens = 10,
	keep_snap_draws = 'draws',
	seed = seed,
	gof = FALSE,
	verbose = FALSE
)

snap_recovery = list(
	auc = auc,
	pi_snap_post = mean(fit$pi_snap),
	rho_true = rho_true,
	sigma_true = sigma_true,
	kappa_true = kappa_true,
	planted = snap_events,
	detected_prob = vapply(snap_events, function(ev) snap_prob[ev[1], ev[2]], numeric(1)),
	max_false_prob = max(neg, na.rm = TRUE),
	drift_only_max_snap = max(fit_drift$snap_prob[, -1], na.rm = TRUE),
	drift_only_terminal_max_snap = max(fit_drift$snap_prob[, n_time], na.rm = TRUE),
	drift_only_terminal_mean_snap = mean(fit_drift$snap_prob[, n_time], na.rm = TRUE),
	drift_only_mean_snap = mean(fit_drift$snap_prob[, -1], na.rm = TRUE),
	drift_only_pi_snap = mean(fit_drift$pi_snap),
	terminal_planted_prob = snap_prob['a24', n_time],
	delta_true = delta_true,
	snap_prob_recovered = snap_prob,
	snap_prob_drift_only = fit_drift$snap_prob,
	pos_probs = pos,
	neg_probs = neg
)

save(snap_recovery, file = paste0(path_out, 'snap_recovery.rda'))
