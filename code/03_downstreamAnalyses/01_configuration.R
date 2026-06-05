# 2025 actor configuration

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))
source(paste0(path_funcs, 'project_metadata.R'))

# actor groups
actor_cats = amer_actor_cats

# latent distance summary
intra_dist = function(pos, set){
	set = intersect(set, rownames(pos))
	if(length(set) < 2){ return(NA_real_) }
	d = as.matrix(dist(pos[set, , drop = FALSE]))
	mean(d[upper.tri(d)])
}
cross_dist = function(pos, set_a, set_b){
	set_a = intersect(set_a, rownames(pos))
	set_b = intersect(set_b, rownames(pos))
	if(length(set_a) < 1 || length(set_b) < 1){ return(NA_real_) }
	full = as.matrix(dist(rbind(pos[set_a, , drop = FALSE], pos[set_b, , drop = FALSE])))
	mean(full[1:length(set_a), (length(set_a) + 1):(length(set_a) + length(set_b)), drop = FALSE])
}

relations = c('un_voting', 'trade', 'sanctions')
configuration = list()
for(rel in relations){
	load(paste0(path_out, 'snap_fit_', rel, '.rda'))
	u_array = fit$U
	years = dimnames(u_array)[[3]]

	# 2025 snaps by category
	actors = rownames(fit$snap_prob)
	snap_tab = data.frame(
		actor = actors,
		snap_prob_2025 = fit$snap_prob[, '2025'],
		snap_prob_2025_receiver = if(!is.null(fit$snap_prob_v)) fit$snap_prob_v[, '2025'] else NA_real_,
		snap_prob_2025_combined = if(!is.null(fit$snap_prob_v)) {
			pmax(fit$snap_prob[, '2025'], fit$snap_prob_v[, '2025'])
		} else {
			fit$snap_prob[, '2025']
		},
		category = NA_character_,
		stringsAsFactors = FALSE
	)
	for(cat_name in names(actor_cats)){
		snap_tab$category[snap_tab$actor %in% actor_cats[[cat_name]]] = cat_name
	}

	# bloc distance paths
	bloc_dist = t(sapply(years, function(y){
		pos = u_array[, , y]
		c(
			intra_allies = intra_dist(pos, actor_cats$allies),
			intra_rivals = intra_dist(pos, actor_cats$rivals),
			allies_to_us = cross_dist(pos, actor_cats$allies, actor_cats$hegemon),
			us_to_rivals = cross_dist(pos, actor_cats$hegemon, actor_cats$rivals),
			allies_to_rivals = cross_dist(pos, actor_cats$allies, actor_cats$rivals)
		)
	}))
	rownames(bloc_dist) = years

	configuration[[rel]] = list(
		snap_table = snap_tab,
		bloc_dist = bloc_dist,
		pi_snap = mean(fit$pi_snap)
	)
}

save(configuration, file = paste0(path_out, 'configuration.rda'))
