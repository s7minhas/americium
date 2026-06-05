# rupture index by domain

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))
source(paste0(path_funcs, 'project_metadata.R'))

# focal actors
order_actors = amer_focal_actors
shock_years = c('1991','2001','2003','2008','2014','2016','2020','2022','2025')

rupture_index = list()
for(rel in c('un_voting', 'trade', 'sanctions')){
	load(paste0(path_out, 'snap_fit_', rel, '.rda'))
	sp = fit$snap_prob
	sp_receiver = if(!is.null(fit$snap_prob_v)) fit$snap_prob_v else NULL
	sp_combined = if(!is.null(sp_receiver)) pmax(sp, sp_receiver) else NULL

	# mean snap by year
	idx = amer_snap_index(sp, order_actors)
	idx_receiver = if(!is.null(sp_receiver)) amer_snap_index(sp_receiver, order_actors) else NULL
	idx_combined = if(!is.null(sp_combined)) amer_snap_index(sp_combined, order_actors) else NULL

	rupture_index[[rel]] = list(
		index_by_year = round(idx, 4),
		receiver_index_by_year = if(!is.null(idx_receiver)) round(idx_receiver, 4) else NULL,
		combined_index_by_year = if(!is.null(idx_combined)) round(idx_combined, 4) else NULL,
		shocks = round(idx[intersect(shock_years, names(idx))], 4),
		rank_2025 = amer_rank_year(idx, '2025'),
		receiver_rank_2025 = if(!is.null(idx_receiver)) amer_rank_year(idx_receiver, '2025') else NULL,
		combined_rank_2025 = if(!is.null(idx_combined)) amer_rank_year(idx_combined, '2025') else NULL,
		n_years = sum(!is.na(idx)),
		actors = intersect(order_actors, rownames(sp)),
		estimand = if(!is.null(sp_receiver)) 'sender/designator-side; receiver and combined indices stored separately' else 'actor position'
	)
}

save(rupture_index, file = paste0(path_out, 'rupture_index.rda'))
