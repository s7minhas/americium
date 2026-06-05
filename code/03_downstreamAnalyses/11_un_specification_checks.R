# un voting spec checks

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))
source(paste0(path_funcs, 'project_metadata.R'))

load_pkg(c('lame', 'data.table'))
set.seed(seed)

load(paste0(path_in, 'arr_list_lame.rda'))
status_file = paste0(path_out, 'un_specification_checks_status.log')
write_status = function(...){
	write(
		paste(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), paste(..., collapse = ' ')),
		file = status_file,
		append = TRUE
	)
}

actor_cats = amer_actor_cats
order_actors = amer_focal_actors

# chain settings
robust_nscan = 3000
robust_burn = 1500
robust_odens = 10
run_heavy = Sys.getenv('RUN_HEAVY_ROBUSTNESS', unset = '0') == '1'

actor_observed = function(net_list){
	obs = sapply(net_list, function(m){
		rowSums(!is.na(m), na.rm = TRUE) + colSums(!is.na(m), na.rm = TRUE) > 0
	})
	rowSums(obs) == ncol(obs)
}

balanced_actors = names(which(actor_observed(arr_list_lame$un_voting)))

subset_net = function(net_list, actors){
	lapply(net_list, function(m){ m[actors, actors] })
}

fit_path = function(name){
	paste0(path_out, 'snap_fit_un_voting_', name, '.rda')
}

load_fit_file = function(file){
	load(file)
	fit
}

fit_un_variant = function(name, y, rank = 2, prior = list(), force = FALSE){
	fit_file = paste0(path_out, 'snap_fit_un_voting_', name, '.rda')
	if(file.exists(fit_file) && !force){
		cached_fit = load_fit_file(fit_file)
		if(amer_fit_is_current(cached_fit) &&
		   identical(cached_fit$americium_build$specification, name)){
			write_status('using cached', name)
			return(cached_fit)
		}
		if(!run_heavy){
			return(NULL)
		}
	}
	if(!run_heavy){
		return(NULL)
	}
	fit_prior = utils::modifyList(amer_snap_prior, prior)
	write_status('fitting', name)
	fit = lame(
		y,
		R = rank,
		family = 'normal',
		symmetric = TRUE,
		dynamic_uv = TRUE,
		dynamic_uv_kind = 'snap',
		prior = fit_prior,
		nscan = robust_nscan,
		burn = robust_burn,
		odens = robust_odens,
		keep_snap_draws = 'draws',
		seed = seed,
		gof = FALSE,
		verbose = FALSE
	)
	fit$americium_build = amer_build_metadata(list(specification = name))
	save(fit, file = fit_file)
	write_status('saved', name)
	fit
}

posthoc_subset_fit = function(fit, actors){
	actors = intersect(actors, rownames(fit$snap_prob))
	x = fit
	x$snap_prob = fit$snap_prob[actors, , drop = FALSE]
	x$U = fit$U[actors, , , drop = FALSE]
	x
}

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

summarize_fit = function(fit, spec, latent_rank, prior_label, refit_type){
	sp = fit$snap_prob
	actors = rownames(sp)
	idx = amer_snap_index(sp, order_actors)
	pos_2024 = fit$U[, , '2024']
	pos_2025 = fit$U[, , '2025']
	data.frame(
		spec = spec,
		refit_type = refit_type,
		latent_rank = latent_rank,
		prior_label = prior_label,
		n_actors = length(actors),
		n_draws = if(length(dim(fit$BETA)) >= 2) dim(fit$BETA)[1] else length(fit$BETA),
		pi_snap = mean(fit$pi_snap, na.rm = TRUE),
		us_2025 = sp['USA', '2025'],
		allies_2025 = mean(sp[intersect(actor_cats$allies, actors), '2025'], na.rm = TRUE),
		rivals_2025 = mean(sp[intersect(actor_cats$rivals, actors), '2025'], na.rm = TRUE),
		hedgers_2025 = mean(sp[intersect(actor_cats$hedgers, actors), '2025'], na.rm = TRUE),
		index_2024 = idx['2024'],
		index_2025 = idx['2025'],
		rank_2025 = amer_rank_year(idx, '2025'),
		intra_allies_2024 = intra_dist(pos_2024, actor_cats$allies),
		intra_allies_2025 = intra_dist(pos_2025, actor_cats$allies),
		allies_to_us_2024 = cross_dist(pos_2024, actor_cats$allies, actor_cats$hegemon),
		allies_to_us_2025 = cross_dist(pos_2025, actor_cats$allies, actor_cats$hegemon),
		stringsAsFactors = FALSE,
		row.names = NULL
	)
}

# load main fit
load(paste0(path_out, 'snap_fit_un_voting.rda'))
main_fit = fit
if(!amer_fit_is_current(main_fit)){
	stop('Main UN voting fit is stale; rerun 02_simMeasure/01_run_snap_lame.R first.')
}

base_specs = list(
	main = list(
		fit = main_fit,
		rank = 2,
		prior_label = 'main',
		refit_type = 'main_analysis'
	),
	balanced_posthoc = list(
		fit = posthoc_subset_fit(main_fit, balanced_actors),
		rank = 2,
		prior_label = 'main',
		refit_type = 'posthoc_restriction'
	)
)

heavy_specs = list(
	r3 = list(
		y = arr_list_lame$un_voting,
		rank = 3,
		prior = list(),
		prior_label = 'main',
		refit_type = 'short_refit'
	),
	kappa_1p5 = list(
		y = arr_list_lame$un_voting,
		rank = 2,
		prior = list(snap_kappa = 1.5),
		prior_label = 'kappa=1.5',
		refit_type = 'short_refit'
	),
	kappa_3 = list(
		y = arr_list_lame$un_voting,
		rank = 2,
		prior = list(snap_kappa = 3),
		prior_label = 'kappa=3',
		refit_type = 'short_refit'
	),
	rare_pi = list(
		y = arr_list_lame$un_voting,
		rank = 2,
		prior = list(snap_pi_a = 1, snap_pi_b = 1999),
		prior_label = 'pi~Beta(1,1999)',
		refit_type = 'short_refit'
	),
	balanced_panel = list(
		y = subset_net(arr_list_lame$un_voting, balanced_actors),
		rank = 2,
		prior = list(),
		prior_label = 'main',
		refit_type = 'balanced_refit'
	)
)

summaries = lapply(names(base_specs), function(nm){
	x = base_specs[[nm]]
	summarize_fit(x$fit, nm, x$rank, x$prior_label, x$refit_type)
})

heavy_status = rbindlist(lapply(names(heavy_specs), function(nm){
	x = heavy_specs[[nm]]
	fit_file = fit_path(nm)
	available_before = file.exists(fit_file)
	fit = fit_un_variant(nm, x$y, rank = x$rank, prior = x$prior)
	available_after = !is.null(fit)
	if(available_after){
		summaries[[length(summaries) + 1]] <<- summarize_fit(
			fit, nm, x$rank, x$prior_label, x$refit_type
		)
	}
	data.frame(
		spec = nm,
		file = fit_file,
		available_before = available_before,
		available_after = available_after,
		run_heavy = run_heavy,
		stringsAsFactors = FALSE
	)
}))

spec_summary = rbindlist(summaries)

un_specification_checks = list(
	summary = as.data.frame(spec_summary),
	balanced_actors = balanced_actors,
	heavy_status = as.data.frame(heavy_status),
	run_heavy = run_heavy,
	robust_chain = list(nscan = robust_nscan, burn = robust_burn, odens = robust_odens)
)

save(un_specification_checks, file = paste0(path_out, 'un_specification_checks.rda'))
