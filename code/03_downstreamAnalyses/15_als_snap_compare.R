# compare als and mcmc snaps

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
source(file.path(path_funcs, 'project_metadata.R'))

lame_dev_path = path.expand('~/Research/netify_dev/lame')
if(requireNamespace('pkgload', quietly = TRUE) && dir.exists(lame_dev_path)){
	pkgload::load_all(lame_dev_path, quiet = TRUE)
} else {
	load_pkg('lame')
}

load(paste0(path_in, 'arr_list_lame.rda'))
load(paste0(path_out, 'snap_fit_un_voting.rda'))
fit_mcmc = fit
rm(fit)

y_un = lapply(arr_list_lame$un_voting, function(m){ m })
mcmc_rho = mean(fit_mcmc$rho_uv, na.rm = TRUE)
mcmc_sigma = mean(fit_mcmc$sigma_uv, na.rm = TRUE)
mcmc_pi = mean(fit_mcmc$pi_snap, na.rm = TRUE)
or_na = function(x, default = NA){
	if(is.null(x)){ default } else { x }
}

als_source_file = file.path(lame_dev_path, 'R', 'lame_snap_als.R')
als_source_md5 = if(file.exists(als_source_file)){
	paste(names(tools::md5sum(als_source_file)), tools::md5sum(als_source_file),
		sep = '=')
} else {
	NA_character_
}

variant_config = list(
	package_default_estimated_hyper = list(
		max_iter = 30
	),
	global_estimated_hyper = list(
		max_iter = 30,
		estimate_rho_uv = TRUE,
		estimate_sigma_uv = TRUE,
		align = 'global'
	),
	default_estimated_hyper = list(
		max_iter = 30,
		estimate_rho_uv = TRUE,
		estimate_sigma_uv = TRUE,
		snap_pi_prior = c(1, 999),
		align = 'global'
	),
	fixed_mcmc_hyper_rare = list(
		max_iter = 20,
		estimate_rho_uv = FALSE,
		estimate_sigma_uv = FALSE,
		rho_uv = mcmc_rho,
		sigma_uv = mcmc_sigma,
		snap_pi_prior = c(1, 999),
		align = 'global'
	),
	fixed_mcmc_hyper_permissive = list(
		max_iter = 20,
		estimate_rho_uv = FALSE,
		estimate_sigma_uv = FALSE,
		rho_uv = mcmc_rho,
		sigma_uv = mcmc_sigma,
		snap_pi_prior = c(1, 1),
		align = 'global'
	),
	fixed_mcmc_hyper_pi_mean = list(
		max_iter = 20,
		estimate_rho_uv = FALSE,
		estimate_sigma_uv = FALSE,
		rho_uv = mcmc_rho,
		sigma_uv = mcmc_sigma,
		snap_pi_prior = c(mcmc_pi * 1000, (1 - mcmc_pi) * 1000),
		align = 'global'
	),
	fixed_mcmc_hyper_pi_mean_long = list(
		max_iter = 80,
		estimate_rho_uv = FALSE,
		estimate_sigma_uv = FALSE,
		rho_uv = mcmc_rho,
		sigma_uv = mcmc_sigma,
		snap_pi_prior = c(mcmc_pi * 1000, (1 - mcmc_pi) * 1000),
		align = 'global'
	),
	fixed_mcmc_hyper_pi_mean_sequential = list(
		max_iter = 20,
		estimate_rho_uv = FALSE,
		estimate_sigma_uv = FALSE,
		rho_uv = mcmc_rho,
		sigma_uv = mcmc_sigma,
		snap_pi_prior = c(mcmc_pi * 1000, (1 - mcmc_pi) * 1000),
		align = 'sequential'
	)
)

force_refit = identical(Sys.getenv('ALS_FORCE_REFIT', unset = '0'), '1')

fit_one = function(label, cfg){
	out_file = paste0(path_out, 'snap_als_fit_un_voting_', label, '.rda')
	if(file.exists(out_file) && !force_refit){
		load(out_file)
		return(fit_als)
	}
	args = c(
		list(
			Y = y_un,
			R = 2,
			family = 'normal',
			symmetric = TRUE,
			snap_kappa = amer_snap_prior$snap_kappa,
			tol = 1e-5,
			verbose = FALSE,
			seed = seed
		),
		cfg
	)
	tm = system.time({
		fit_als = do.call(lame_snap_als, args)
	})
	fit_als$americium_build = amer_build_metadata(list(
		relation = 'un_voting',
		estimator = 'lame_snap_als',
		variant = label,
		lame_snap_als_source_md5 = als_source_md5,
		elapsed_seconds = unname(tm[['elapsed']])
	))
	save(fit_als, file = out_file)
	fit_als
}

align_mcmc = function(sp){
	fit_mcmc$snap_prob[rownames(sp), colnames(sp), drop = FALSE]
}

threshold_share = function(x, threshold){
	mean(x > threshold, na.rm = TRUE)
}

snap_index_rank = function(sp, year){
	idx = amer_snap_index(sp, amer_focal_actors)
	unname(rank(-idx[!is.na(idx)], ties.method = 'min')[as.character(year)])
}

summarize_fit = function(label, fit_als){
	sp_als = fit_als$snap_prob
	sp_mcmc = align_mcmc(sp_als)
	focal = intersect(amer_focal_actors, rownames(sp_als))
	year_2025 = '2025'
	elapsed_seconds = if(is.null(fit_als$americium_build$elapsed_seconds)){
		NA_real_
	} else {
		fit_als$americium_build$elapsed_seconds
	}
	common = is.finite(sp_als) & is.finite(sp_mcmc)
	common_2025 = is.finite(sp_als[, year_2025]) & is.finite(sp_mcmc[, year_2025])
	als_hi = sp_als > 0.5
	mcmc_hi = sp_mcmc > 0.5
	data.frame(
		variant = label,
		iterations = fit_als$iterations,
		converged = fit_als$converged,
		elapsed_seconds = elapsed_seconds,
		als_rho_uv = fit_als$rho_uv,
		als_sigma_uv = fit_als$sigma_uv,
		als_pi_snap = fit_als$pi_snap,
		als_align = or_na(fit_als$meta$align),
		als_hyper_update = or_na(fit_als$meta$hyper_update),
		als_sigma_update_floor = or_na(fit_als$meta$sigma_update_floor,
			NA_real_),
		als_final_max_snap_delta = or_na(
			fit_als$convergence$final_max_snap_delta, NA_real_),
		als_final_snap_delta_q95 = or_na(
			fit_als$convergence$final_snap_delta_q95, NA_real_),
		als_final_snap_delta_q99 = or_na(
			fit_als$convergence$final_snap_delta_q99, NA_real_),
		als_final_snap_delta_conv = or_na(
			fit_als$convergence$final_snap_delta_conv, NA_real_),
		als_final_n_snap_delta_gt_tol = or_na(
			fit_als$convergence$final_n_snap_delta_gt_tol, NA_integer_),
		als_final_n_snap_class_changed = or_na(
			fit_als$convergence$final_n_snap_class_changed, NA_integer_),
		als_final_share_snap_class_changed = or_na(
			fit_als$convergence$final_share_snap_class_changed, NA_real_),
		als_snap_stable = or_na(fit_als$convergence$snap_stable, NA),
		als_max_snap_stable = or_na(fit_als$convergence$max_snap_stable,
			NA),
		als_objective_increases = or_na(
			fit_als$convergence$objective_increases, NA_integer_),
		als_max_objective_increase = or_na(
			fit_als$convergence$max_objective_increase, NA_real_),
		mcmc_rho_uv = mcmc_rho,
		mcmc_sigma_uv = mcmc_sigma,
		mcmc_pi_snap = mcmc_pi,
		cor_all_actor_years = cor(as.vector(sp_als), as.vector(sp_mcmc),
			use = 'complete.obs'),
		cor_2025 = cor(sp_als[common_2025, year_2025], sp_mcmc[common_2025, year_2025]),
		agreement_gt_50_all = mean(als_hi[common] == mcmc_hi[common], na.rm = TRUE),
		als_all_2025_mean = mean(sp_als[, year_2025], na.rm = TRUE),
		mcmc_all_2025_mean = mean(sp_mcmc[, year_2025], na.rm = TRUE),
		als_all_2025_gt_50 = threshold_share(sp_als[, year_2025], 0.5),
		mcmc_all_2025_gt_50 = threshold_share(sp_mcmc[, year_2025], 0.5),
		als_focal_2025_mean = mean(sp_als[focal, year_2025], na.rm = TRUE),
		mcmc_focal_2025_mean = mean(sp_mcmc[focal, year_2025], na.rm = TRUE),
		als_focal_2025_gt_50 = threshold_share(sp_als[focal, year_2025], 0.5),
		mcmc_focal_2025_gt_50 = threshold_share(sp_mcmc[focal, year_2025], 0.5),
		als_focal_2025_gt_99 = threshold_share(sp_als[focal, year_2025], 0.99),
		mcmc_focal_2025_gt_99 = threshold_share(sp_mcmc[focal, year_2025], 0.99),
		als_2025_rank = snap_index_rank(sp_als, year_2025),
		mcmc_2025_rank = snap_index_rank(sp_mcmc, year_2025),
		stringsAsFactors = FALSE
	)
}

category_2025 = function(label, fit_als){
	sp_als = fit_als$snap_prob
	sp_mcmc = align_mcmc(sp_als)
	do.call(rbind, lapply(names(amer_actor_cats), function(cat_name){
		actors = intersect(amer_actor_cats[[cat_name]], rownames(sp_als))
		data.frame(
			variant = label,
			category = amer_category_labels[[cat_name]],
			n = length(actors),
			als_mean_2025 = mean(sp_als[actors, '2025'], na.rm = TRUE),
			mcmc_mean_2025 = mean(sp_mcmc[actors, '2025'], na.rm = TRUE),
			als_gt_50_2025 = threshold_share(sp_als[actors, '2025'], 0.5),
			mcmc_gt_50_2025 = threshold_share(sp_mcmc[actors, '2025'], 0.5),
			stringsAsFactors = FALSE
		)
	}))
}

top_years = function(label, fit_als, n_top = 8){
	sp_als = fit_als$snap_prob
	sp_mcmc = align_mcmc(sp_als)
	idx_als = amer_snap_index(sp_als, amer_focal_actors)
	idx_mcmc = amer_snap_index(sp_mcmc, amer_focal_actors)
	make = function(estimator, idx){
		ord = order(idx, decreasing = TRUE, na.last = NA)
		ord = ord[seq_len(min(n_top, length(ord)))]
		data.frame(
			variant = label,
			estimator = estimator,
			rank = seq_along(ord),
			year = names(idx)[ord],
			focal_index = unname(idx[ord]),
			stringsAsFactors = FALSE
		)
	}
	rbind(make('ALS', idx_als), make('MCMC', idx_mcmc))
}

fits = lapply(names(variant_config), function(label){
	message('Fitting/loading ', label)
	fit_one(label, variant_config[[label]])
})
names(fits) = names(variant_config)

comparison_summary = do.call(rbind, Map(summarize_fit, names(fits), fits))
comparison_category_2025 = do.call(rbind, Map(category_2025, names(fits), fits))
comparison_top_years = do.call(rbind, Map(top_years, names(fits), fits))

write.csv(comparison_summary,
	file = paste0(path_out, 'snap_als_comparison_summary.csv'),
	row.names = FALSE)
write.csv(comparison_category_2025,
	file = paste0(path_out, 'snap_als_comparison_category_2025.csv'),
	row.names = FALSE)
write.csv(comparison_top_years,
	file = paste0(path_out, 'snap_als_comparison_top_years.csv'),
	row.names = FALSE)

snap_als_comparison = list(
	summary = comparison_summary,
	category_2025 = comparison_category_2025,
	top_years = comparison_top_years,
	mcmc_hyperparameters = c(rho_uv = mcmc_rho, sigma_uv = mcmc_sigma,
		pi_snap = mcmc_pi),
	als_source_md5 = als_source_md5,
	note = 'ALS snap scores are approximate point-estimator scores, not posterior probabilities.'
)
save(snap_als_comparison, file = paste0(path_out, 'snap_als_comparison.rda'))

print(comparison_summary)
