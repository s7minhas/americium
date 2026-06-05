# project metadata

amer_actor_cats = list(
	hegemon = 'USA',
	allies = c('FRA','DEU','GBR','ITA','CAN','JPN','AUS','NLD','ESP','KOR','BEL','NOR','DNK','POL'),
	rivals = c('RUS','CHN','PRK','IRN','BLR','SYR','VEN','CUB','NIC'),
	hedgers = c('IND','TUR','BRA','ZAF','IDN','SAU','MEX','EGY','ARE','VNM','THA')
)

amer_category_labels = c(
	hegemon = 'Hegemon',
	allies = 'Core Allies',
	rivals = 'Strategic Rivals',
	hedgers = 'Pivotal Hedgers'
)

amer_category_colors = c(
	'Hegemon' = '#d62728',
	'Core Allies' = '#1f77b4',
	'Strategic Rivals' = '#2ca02c',
	'Pivotal Hedgers' = 'grey50'
)

amer_focal_actors = unique(unlist(amer_actor_cats, use.names = FALSE))

amer_actor_type = function(codes){
	out = rep(NA_character_, length(codes))
	for(cat_name in names(amer_actor_cats)){
		out[codes %in% amer_actor_cats[[cat_name]]] = amer_category_labels[[cat_name]]
	}
	out
}

amer_actor_category_table = do.call(rbind, lapply(names(amer_actor_cats), function(cat_name){
	data.frame(
		actor = amer_actor_cats[[cat_name]],
		category = amer_category_labels[[cat_name]],
		category_key = cat_name,
		focal_index = TRUE,
		stringsAsFactors = FALSE
	)
}))
rownames(amer_actor_category_table) = NULL

amer_current_build = list(
	un_diagonal = 'NA',
	lame_snap_hyperparameters = 'drift_only_transitions',
	lame_sigma_uv_prior = 'configured_ig_drift_only',
	lame_snap_transition = 'forward_gated_by_next_delta',
	lame_snap_pi_prior = 'Beta(1,999)',
	lame_keep_snap_draws = 'draws',
	snap_first_year = 'NA'
)

amer_snap_prior = list(
	snap_kappa = 2,
	snap_pi_a = 1,
	snap_pi_b = 999
)

amer_lame_source_md5 = function(){
	src = c(
		'/home/s7m/Research/netify_dev/lame/R/lame.R',
		'/home/s7m/Research/netify_dev/lame/src/dynamic_latent_factors.cpp'
	)
	src = src[file.exists(src)]
	if(length(src) == 0){ return(NA_character_) }
	paste(names(tools::md5sum(src)), tools::md5sum(src), sep = '=', collapse = ';')
}

amer_fit_is_current = function(fit){
	if(is.null(fit$americium_build)){ return(FALSE) }
	labels_current = all(vapply(names(amer_current_build), function(nm){
		identical(fit$americium_build[[nm]], amer_current_build[[nm]])
	}, logical(1)))
	src_current = identical(fit$americium_build$lame_source_md5, amer_lame_source_md5())
	labels_current && src_current
}

amer_build_metadata = function(extra = list()){
	base = c(
		amer_current_build,
		list(
			lame_version = as.character(utils::packageVersion('lame')),
			lame_source_md5 = amer_lame_source_md5(),
			built_at = as.character(Sys.time())
		)
	)
	if(length(extra) > 0){
		base[names(extra)] = extra
	}
	base
}

amer_snap_index = function(sp, actors){
	keep = intersect(actors, rownames(sp))
	x = sp[keep, , drop = FALSE]
	idx = colMeans(x, na.rm = TRUE)
	idx[colSums(!is.na(x)) == 0] = NA_real_
	idx
}

amer_rank_year = function(idx, year){
	eligible = !is.na(idx)
	ranks = rank(-idx[eligible], ties.method = 'min')
	unname(ranks[as.character(year)])
}
