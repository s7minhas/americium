# compare transition models

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))
source(paste0(path_funcs, 'project_metadata.R'))

load_pkg('lame')
set.seed(seed)

load(paste0(path_in, 'arr_list_lame.rda'))
un = arr_list_lame$un_voting
status_file = paste0(path_out, 'model_comparison_status.log')
write_status = function(...){
	write(
		paste(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), paste(..., collapse = ' ')),
		file = status_file,
		append = TRUE
	)
}

# fit and score each model
model_comparison = list()
for(k in c('ar1', 'snap', 't')){
	write_status('fitting comparison model', k)
	ll_path = paste0(path_out, 'loglik_', k, '.bin')
	prior = if(k == 'snap') amer_snap_prior else list()
	fit = lame(
		un, R = 2, family = 'normal', symmetric = TRUE,
		dynamic_uv = TRUE, dynamic_uv_kind = k,
		prior = prior,
		nscan = 4000, burn = 2000, odens = 40,
		seed = seed, gof = FALSE, verbose = FALSE,
		save_log_lik = TRUE, log_lik_path = ll_path
	)
	model_comparison[[k]] = list(
		waic = tryCatch(waic(fit), error = function(e) conditionMessage(e)),
		pi_snap = if(k == 'snap') mean(fit$pi_snap, na.rm = TRUE) else NA_real_,
		build = amer_build_metadata(list(comparison_model = k))
	)
	write_status('finished comparison model', k)
}

save(model_comparison, file = paste0(path_out, 'model_comparison.rda'))
