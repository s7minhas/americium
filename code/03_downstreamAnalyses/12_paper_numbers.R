# refresh paper numbers

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
load_pkg('data.table')

actor_cats = amer_actor_cats
order_actors = amer_focal_actors
shock_years = c('1991','2001','2003','2008','2014','2016','2020','2022','2025')

group_mean = function(sp, actors, year = '2025'){
	keep = intersect(actors, rownames(sp))
	if(length(keep) == 0){ return(NA_real_) }
	mean(sp[keep, year], na.rm = TRUE)
}

fit_has_current_build_labels = function(fit){
	if(is.null(fit$americium_build)){ return(FALSE) }
	all(vapply(names(amer_current_build), function(nm){
		identical(fit$americium_build[[nm]], amer_current_build[[nm]])
	}, logical(1)))
}

fit_has_recorded_source = function(fit){
	!is.null(fit$americium_build$lame_source_md5) &&
		nzchar(fit$americium_build$lame_source_md5)
}

summarize_fit = function(rel){
	load(paste0(path_out, 'snap_fit_', rel, '.rda'))
	if(!fit_has_current_build_labels(fit)){
		stop(paste0('Fit is missing current build labels: ', rel))
	}
	sp = fit$snap_prob
	sp_receiver = if(!is.null(fit$snap_prob_v)) fit$snap_prob_v else NULL
	sp_combined = if(!is.null(sp_receiver)) pmax(sp, sp_receiver) else NULL
	idx = amer_snap_index(sp, order_actors)
	idx_receiver = if(!is.null(sp_receiver)) amer_snap_index(sp_receiver, order_actors) else NULL
	idx_combined = if(!is.null(sp_combined)) amer_snap_index(sp_combined, order_actors) else NULL
	data.frame(
		relation = rel,
		n_actors = nrow(sp),
		n_years = ncol(sp),
		n_draws = if(length(dim(fit$BETA)) >= 2) dim(fit$BETA)[1] else length(fit$BETA),
		pi_snap = mean(fit$pi_snap, na.rm = TRUE),
		index_2024 = idx['2024'],
		index_2025 = idx['2025'],
		rank_2025 = amer_rank_year(idx, '2025'),
		receiver_index_2025 = if(!is.null(idx_receiver)) idx_receiver['2025'] else NA_real_,
		receiver_rank_2025 = if(!is.null(idx_receiver)) amer_rank_year(idx_receiver, '2025') else NA_real_,
		combined_index_2025 = if(!is.null(idx_combined)) idx_combined['2025'] else NA_real_,
		combined_rank_2025 = if(!is.null(idx_combined)) amer_rank_year(idx_combined, '2025') else NA_real_,
		us_2025 = sp['USA', '2025'],
		us_2025_receiver = if(!is.null(sp_receiver)) sp_receiver['USA', '2025'] else NA_real_,
		us_2025_combined = if(!is.null(sp_combined)) sp_combined['USA', '2025'] else NA_real_,
		allies_2025 = group_mean(sp, actor_cats$allies),
		rivals_2025 = group_mean(sp, actor_cats$rivals),
		hedgers_2025 = group_mean(sp, actor_cats$hedgers),
		stringsAsFactors = FALSE,
		row.names = NULL
	)
}

relations = c('un_voting', 'trade', 'sanctions')
relation_summary = do.call(rbind, lapply(relations, summarize_fit))

draw_summary = function(x){
	x = as.numeric(x)
	x = x[is.finite(x)]
	if(length(x) == 0){
		return(data.frame(
			mean = NA_real_,
			q025 = NA_real_,
			q50 = NA_real_,
			q975 = NA_real_,
			stringsAsFactors = FALSE
		))
	}
	data.frame(
		mean = mean(x),
		q025 = unname(stats::quantile(x, 0.025, names = FALSE)),
		q50 = stats::median(x),
		q975 = unname(stats::quantile(x, 0.975, names = FALSE)),
		stringsAsFactors = FALSE
	)
}

snap_draw_matrix = function(fit, actors, side = c('sender','receiver','combined')){
	side = match.arg(side)
	if(side == 'sender'){
		arr = fit$snap_draws
	} else if(side == 'receiver'){
		arr = fit$snap_draws_v
	} else {
		arr = if(!is.null(fit$snap_draws_v)) pmax(fit$snap_draws, fit$snap_draws_v) else NULL
	}
	if(is.null(arr)){ return(NULL) }
	actors = intersect(actors, dimnames(arr)[[2]])
	years = dimnames(arr)[[3]]
	out = sapply(seq_along(years), function(j){
		rowMeans(arr[, actors, j, drop = FALSE], na.rm = TRUE)
	})
	colnames(out) = years
	out
}

snap_rank_draws = function(mat, target = '2025'){
	keep = colSums(is.finite(mat)) > 0
	mat = mat[, keep, drop = FALSE]
	if(!(target %in% colnames(mat))){
		stop(paste0('Target year not found in draw matrix: ', target))
	}
	apply(mat, 1, function(x){
		rank(-x, ties.method = 'min')[target]
	})
}

relation_draw_row = function(rel, side = c('sender','receiver','combined')){
	side = match.arg(side)
	load(paste0(path_out, 'snap_fit_', rel, '.rda'))
	mat = snap_draw_matrix(fit, order_actors, side = side)
	if(is.null(mat)){ return(NULL) }
	summ = draw_summary(mat[, '2025'])
	ranks = snap_rank_draws(mat, target = '2025')
	bench = c('1991','2001','2016','2022')
	bench = intersect(bench, colnames(mat))
	data.frame(
		relation = rel,
		index_side = side,
		index_2025_mean = summ$mean,
		index_2025_q025 = summ$q025,
		index_2025_median = summ$q50,
		index_2025_q975 = summ$q975,
		rank_2025_mean = mean(ranks, na.rm = TRUE),
		rank_2025_median = stats::median(ranks, na.rm = TRUE),
		pr_rank_1 = mean(ranks == 1, na.rm = TRUE),
		pr_rank_le_3 = mean(ranks <= 3, na.rm = TRUE),
		pr_2025_gt_selected_benchmarks = mean(vapply(
			bench,
			function(yr){ mean(mat[, '2025'] > mat[, yr], na.rm = TRUE) },
			numeric(1)
		), na.rm = TRUE),
		stringsAsFactors = FALSE,
		row.names = NULL
	)
}

relation_draw_summary = do.call(rbind, Filter(
	Negate(is.null),
	list(
		relation_draw_row('un_voting', 'sender'),
		relation_draw_row('trade', 'sender'),
		relation_draw_row('sanctions', 'sender'),
		relation_draw_row('sanctions', 'receiver'),
		relation_draw_row('sanctions', 'combined')
	)
))

chain_summary = function(x){
	x = as.numeric(x)
	x = x[is.finite(x)]
	if(length(x) == 0){
		return(c(mean = NA_real_, sd = NA_real_, ess = NA_real_))
	}
	x_sd = stats::sd(x)
	if(is.finite(x_sd) && x_sd < .Machine$double.eps){
		return(c(mean = mean(x), sd = x_sd, ess = length(x)))
	}
	ess = NA_real_
	if(requireNamespace('coda', quietly = TRUE)){
		ess = tryCatch(
			as.numeric(coda::effectiveSize(x)),
			error = function(e){ NA_real_ }
		)
	}
	c(mean = mean(x), sd = x_sd, ess = ess)
}

mcmc_row = function(rel){
	load(paste0(path_out, 'snap_fit_', rel, '.rda'))
	call_args = as.list(fit$call)
	call_integer = function(name, fallback = NA_integer_){
		if(!is.null(fit$americium_build[[paste0('snap_', name)]])){
			return(as.integer(fit$americium_build[[paste0('snap_', name)]]))
		}
		x = call_args[[name]]
		if(is.null(x)){ return(fallback) }
		if(is.numeric(x) || is.integer(x)){ return(as.integer(x)) }
		fallback
	}
	pi_sum = chain_summary(fit$pi_snap)
	rho_sum = chain_summary(fit$rho_uv)
	sigma_sum = chain_summary(fit$sigma_uv)
	data.frame(
		relation = rel,
		family = fit$family,
		nscan = call_integer('nscan'),
		burn = call_integer('burn'),
		thin = call_integer('odens'),
		retained_draws = length(fit$pi_snap),
		pi_mean = pi_sum[['mean']],
		pi_sd = pi_sum[['sd']],
		rho_uv_mean = rho_sum[['mean']],
		rho_uv_sd = rho_sum[['sd']],
		sigma_uv_mean = sigma_sum[['mean']],
		sigma_uv_sd = sigma_sum[['sd']],
		min_ess = min(c(pi_sum[['ess']], rho_sum[['ess']], sigma_sum[['ess']]), na.rm = TRUE),
		terminated_early = isTRUE(fit$terminated_early),
		source_hash_recorded = fit_has_recorded_source(fit),
		local_source_match = identical(fit$americium_build$lame_source_md5, amer_lame_source_md5()),
		stringsAsFactors = FALSE,
		row.names = NULL
	)
}

mcmc_diagnostics = do.call(rbind, lapply(relations, mcmc_row))

load(paste0(path_out, 'snap_fit_un_voting.rda'))
un_sp = fit$snap_prob
un_year_index = amer_snap_index(un_sp, order_actors)

un_category_draw_sets = list(
	'Hegemon' = 'USA',
	'Core allies' = actor_cats$allies,
	'Strategic rivals' = actor_cats$rivals,
	'Pivotal hedgers' = actor_cats$hedgers
)
un_category_draw_summary = do.call(rbind, lapply(names(un_category_draw_sets), function(category){
	actors = intersect(un_category_draw_sets[[category]], dimnames(fit$snap_draws)[[2]])
	vals = rowMeans(fit$snap_draws[, actors, '2025', drop = FALSE], na.rm = TRUE)
	summ = draw_summary(vals)
	data.frame(
		category = category,
		n_actors = length(actors),
		index_2025_mean = summ$mean,
		index_2025_q025 = summ$q025,
		index_2025_median = summ$q50,
		index_2025_q975 = summ$q975,
		stringsAsFactors = FALSE,
		row.names = NULL
	)
}))

un_shocks = data.frame(
	year = intersect(shock_years, names(un_year_index)),
	rupture_index = as.numeric(un_year_index[intersect(shock_years, names(un_year_index))]),
	rank = as.numeric(rank(-un_year_index, ties.method = 'min')[intersect(shock_years, names(un_year_index))]),
	row.names = NULL
)
un_top_2025 = data.frame(
	actor = names(sort(un_sp[, '2025'], decreasing = TRUE)[1:20]),
	snap_prob_2025 = as.numeric(sort(un_sp[, '2025'], decreasing = TRUE)[1:20]),
	row.names = NULL
)

observed_actors = function(fit, year){
	y_yr = fit$Y[, , year]
	rownames(y_yr)[rowSums(!is.na(y_yr)) + colSums(!is.na(y_yr)) > 0]
}

index_row = function(label, actors){
	idx = amer_snap_index(un_sp, actors)
	actors_2025 = intersect(actors, observed_actors(fit, '2025'))
	actors_index = intersect(actors, rownames(un_sp))
	data.frame(
		contrast = label,
		n_index = length(actors_index),
		n_observed_2025 = length(actors_2025),
		index_2024 = idx['2024'],
		index_2025 = idx['2025'],
		rank_2025 = amer_rank_year(idx, '2025'),
		stringsAsFactors = FALSE,
		row.names = NULL
	)
}

observed_2025 = observed_actors(fit, '2025')
non_focal_2025 = setdiff(observed_2025, order_actors)
low_exposure_actors = character(0)
low_exposure_note = NA_character_
who_snaps_file = paste0(path_out, 'who_snaps_full.rda')
if(file.exists(who_snaps_file)){
	load(who_snaps_file)
	ws_dat = who_snaps_full$data
	ws_dat = ws_dat[!(ws_dat$actor %in% order_actors), ]
	us_trade_cut = quantile(ws_dat$us_trade_dep, 0.25, na.rm = TRUE)
	low_exposure_actors = ws_dat$actor[ws_dat$us_trade_dep <= us_trade_cut]
	low_exposure_note = paste0('non-focal bottom quartile of 2024 U.S. trade dependence, cutoff=', round(us_trade_cut, 4))
}
if(length(low_exposure_actors) == 0){
	y_2024 = fit$Y[, , '2024']
	obs_2024 = observed_actors(fit, '2024')
	non_focal_2024 = setdiff(obs_2024, order_actors)
	us_agreement = y_2024[non_focal_2024, 'USA']
	us_agree_cut = quantile(us_agreement, 0.25, na.rm = TRUE)
	low_exposure_actors = names(us_agreement)[us_agreement <= us_agree_cut]
	low_exposure_note = paste0('non-focal bottom quartile of 2024 U.S. agreement, cutoff=', round(us_agree_cut, 4))
}

un_denominator_contrasts = rbind(
	index_row('Theory-selected focal actors', order_actors),
	index_row('All observed states', observed_2025),
	index_row('Observed non-focal states', non_focal_2025),
	index_row('Low-exposure non-focal contrast', low_exposure_actors)
)
un_denominator_contrasts$definition = c(
	'hegemon, core allies, strategic rivals, and pivotal hedgers',
	'all states observed in the 2025 UN-voting network',
	'observed states not in the theory-selected focal set',
	low_exposure_note
)

if(nrow(un_denominator_contrasts) < 1 ||
		un_denominator_contrasts$n_index[un_denominator_contrasts$contrast == 'Theory-selected focal actors'] != 35){
	stop('Focal denominator must contain exactly 35 theory-selected actors.')
}
if(!('VEN' %in% rownames(un_sp)) || is.na(un_sp['VEN', '2025']) ||
		round(un_sp['VEN', '2025'], 3) != 0.068){
	stop('Expected Venezuela 2025 convention is not present in the UN focal index.')
}

load(paste0(path_in, 'un_votes_long.rda'))
setDT(un_votes_long)
pos_map = c(Y = 1, A = 0, N = -1)
un_votes_clean = un_votes_long[year >= 1989 & ms_vote %in% c('Y', 'N', 'A')]

agreement_pair = function(yr, a, b){
	d = un_votes_clean[year == yr & ms_code %in% c(a, b)]
	w = dcast(d, resolution ~ ms_code, value.var = 'ms_vote')
	if(!all(c(a, b) %in% names(w))){ return(NA_real_) }
	w = w[!is.na(get(a)) & !is.na(get(b))]
	if(nrow(w) == 0){ return(NA_real_) }
	1 - mean(abs(pos_map[w[[a]]] - pos_map[w[[b]]]) / 2)
}

mean_agreement_between = function(yr, group_a, group_b, paired = FALSE){
	group_a = unique(group_a)
	group_b = unique(group_b)
	pairs = expand.grid(a = group_a, b = group_b, stringsAsFactors = FALSE)
	if(paired){
		pairs = pairs[pairs$a < pairs$b, , drop = FALSE]
	} else {
		pairs = pairs[pairs$a != pairs$b, , drop = FALSE]
	}
	if(nrow(pairs) == 0){ return(NA_real_) }
	mean(mapply(function(a, b){ agreement_pair(yr, a, b) }, pairs$a, pairs$b), na.rm = TRUE)
}

vote_profile = function(yr, actor){
	x = un_votes_clean[year == yr & ms_code == actor]
	if(nrow(x) == 0){
		return(data.frame(n_votes = NA_integer_, yes = NA_real_, no = NA_real_, abstain = NA_real_))
	}
	data.frame(
		n_votes = nrow(x),
		yes = mean(x$ms_vote == 'Y'),
		no = mean(x$ms_vote == 'N'),
		abstain = mean(x$ms_vote == 'A')
	)
}

europe_allies = c('FRA','DEU','GBR','ITA','NLD','ESP','BEL','NOR','DNK','POL')
indo_pacific_allies = c('JPN','KOR','AUS')
raw_metric = function(label, fun){
	val_2024 = fun(2024)
	val_2025 = fun(2025)
	data.frame(
		metric = label,
		value_2024 = val_2024,
		value_2025 = val_2025,
		change = val_2025 - val_2024,
		stringsAsFactors = FALSE,
		row.names = NULL
	)
}
un_raw_vote_diagnostics = rbind(
	raw_metric('U.S. share voting yes', function(y){ vote_profile(y, 'USA')$yes }),
	raw_metric('U.S. share voting no', function(y){ vote_profile(y, 'USA')$no }),
	raw_metric('U.S. share abstaining', function(y){ vote_profile(y, 'USA')$abstain }),
	raw_metric('U.S. agreement with core allies', function(y){ mean_agreement_between(y, 'USA', actor_cats$allies) }),
	raw_metric('U.S. agreement with strategic rivals', function(y){ mean_agreement_between(y, 'USA', actor_cats$rivals) }),
	raw_metric('U.S. agreement with pivotal hedgers', function(y){ mean_agreement_between(y, 'USA', actor_cats$hedgers) }),
	raw_metric('Agreement within core allies', function(y){ mean_agreement_between(y, actor_cats$allies, actor_cats$allies, paired = TRUE) }),
	raw_metric('Europe--Indo-Pacific allied agreement', function(y){ mean_agreement_between(y, europe_allies, indo_pacific_allies) }),
	raw_metric('Core-allied agreement with strategic rivals', function(y){ mean_agreement_between(y, actor_cats$allies, actor_cats$rivals) }),
	raw_metric('Core-allied agreement with pivotal hedgers', function(y){ mean_agreement_between(y, actor_cats$allies, actor_cats$hedgers) })
)
un_raw_vote_diagnostics$interpretation = c(
	'U.S. vote profile in the annual UN agenda',
	'U.S. vote profile in the annual UN agenda',
	'U.S. vote profile in the annual UN agenda',
	'direct U.S.-ally divergence',
	'U.S. does not simply join the rival voting profile',
	'U.S. diverges from pivotal hedgers as well',
	'allied cohesion excluding U.S. ties',
	'Europe and Indo-Pacific allies remain close without U.S. ties',
	'cross-bloc agreement also changes on the annual agenda',
	'allied relation to hedgers changes on the annual agenda'
)

load(paste0(path_out, 'leave_hegemon_out.rda'))
un_leaveout_summary = leave_hegemon_out[leave_hegemon_out$year %in% c(2024, 2025), ]
un_leaveout_summary$full_change = un_leaveout_summary$full - un_leaveout_summary$full[un_leaveout_summary$year == 2024]
un_leaveout_summary$excl_us_change = un_leaveout_summary$excl_us - un_leaveout_summary$excl_us[un_leaveout_summary$year == 2024]
un_leaveout_summary$excl_us_rus_chn_change = un_leaveout_summary$excl_us_rus_chn - un_leaveout_summary$excl_us_rus_chn[un_leaveout_summary$year == 2024]
un_leaveout_summary$excl_us_rivals_change = un_leaveout_summary$excl_us_rivals - un_leaveout_summary$excl_us_rivals[un_leaveout_summary$year == 2024]

benchmark_windows = list(
	'Post-Cold War transition' = 1990:1992,
	'9/11, Afghanistan, and Iraq' = 2001:2003,
	'Global financial crisis' = 2008:2009,
	'Crimea annexation' = 2014:2015,
	'First Trump administration' = 2017:2020,
	'Russia full-scale invasion' = 2022:2023,
	'Second Trump administration' = 2025
)
un_historical_calibration = do.call(rbind, lapply(names(benchmark_windows), function(label){
	yrs = intersect(as.character(benchmark_windows[[label]]), names(un_year_index))
	vals = un_year_index[yrs]
	max_year = yrs[which.max(vals)]
	data.frame(
		benchmark = label,
		transition_years = paste(range(as.integer(yrs)), collapse = '--'),
		mean_index = mean(vals, na.rm = TRUE),
		max_index = max(vals, na.rm = TRUE),
		max_year = max_year,
		max_rank = amer_rank_year(un_year_index, max_year),
		ratio_to_2025 = max(vals, na.rm = TRUE) / un_year_index['2025'],
		stringsAsFactors = FALSE,
		row.names = NULL
	)
}))

format_top_years = function(vals, n = 4){
	vals = vals[is.finite(vals)]
	vals = sort(vals, decreasing = TRUE)
	vals = vals[seq_len(min(n, length(vals)))]
	paste(paste0(names(vals), ' (', sprintf('%.2f', vals), ')'), collapse = '; ')
}

un_country_snap_history = do.call(rbind, lapply(order_actors, function(actor){
	vals = un_sp[actor, ]
	vals = vals[is.finite(vals)]
	top_vals = sort(vals, decreasing = TRUE)
	top_year = names(top_vals)[1]
	data.frame(
		actor = actor,
		category = amer_actor_type(actor),
		max_year = top_year,
		max_snap_prob = unname(top_vals[1]),
		snap_prob_2025 = un_sp[actor, '2025'],
		is_2025_max = top_year == '2025',
		major_snap_years = paste(names(vals)[vals >= 0.75], collapse = ', '),
		top_snap_years = format_top_years(vals, n = 4),
		stringsAsFactors = FALSE,
		row.names = NULL
	)
}))

period_windows = data.frame(
	period = c(
		'Post-Cold War',
		'Late 1990s',
		'9/11 and Bush',
		'Financial crisis',
		'Crimea aftermath',
		'Trump I',
		'Ukraine invasion',
		'2025 rupture'
	),
	years = I(list(
		1990:1992,
		1996:1998,
		2001:2003,
		2008:2009,
		2014:2015,
		2017:2020,
		2022:2023,
		2025
	)),
	stringsAsFactors = FALSE
)

category_mean = function(category_key, year){
	actors = intersect(actor_cats[[category_key]], rownames(un_sp))
	mean(un_sp[actors, as.character(year)], na.rm = TRUE)
}

un_period_category = do.call(rbind, lapply(seq_len(nrow(period_windows)), function(k){
	period = period_windows$period[k]
	yrs = intersect(as.character(period_windows$years[[k]]), colnames(un_sp))
	do.call(rbind, lapply(names(actor_cats), function(category_key){
		vals = vapply(yrs, function(year){ category_mean(category_key, year) }, numeric(1))
		max_year = names(vals)[which.max(vals)]
		data.frame(
			period = period,
			category_key = category_key,
			category = amer_category_labels[[category_key]],
			years = paste(range(as.integer(yrs)), collapse = '--'),
			mean_snap = mean(vals, na.rm = TRUE),
			max_snap = max(vals, na.rm = TRUE),
			max_year = max_year,
			stringsAsFactors = FALSE,
			row.names = NULL
		)
	}))
}))
un_period_category$period = factor(un_period_category$period, levels = period_windows$period)
un_period_category$category = factor(un_period_category$category, levels = amer_category_labels)

episode_years = c('1991','1998','2001','2009','2013','2019','2023','2025')
episode_signatures = c(
	'1991' = 'post-Cold War adjustment; allies, rivals, and hedgers move more than the United States',
	'1998' = 'late-1990s rival and hedger re-coordination, with Turkey, China, Belarus, and Japan prominent',
	'2001' = 'post-9/11 and Bush-era re-coordination; rivals and hedgers are most exposed',
	'2009' = 'financial-crisis-era diplomatic re-coordination; broad but less synchronized than 2025',
	'2013' = 'pre-Crimea and Arab Spring-era bloc repositioning; allies and rivals are high',
	'2019' = 'Trump I-era re-coordination; broad but moderate relative to 2025',
	'2023' = 'Russia-invasion aftermath; United States, allies, and hedgers are high',
	'2025' = 'bloc-preserving hegemonic exit; all focal categories snap together'
)
un_episode_typology = do.call(rbind, lapply(episode_years, function(year){
	cat_vals = vapply(names(actor_cats), function(category_key){
		category_mean(category_key, year)
	}, numeric(1))
	names(cat_vals) = amer_category_labels[names(cat_vals)]
	dominant = names(cat_vals)[cat_vals >= 0.5]
	data.frame(
		year = year,
		focal_index = un_year_index[year],
		us_snap = un_sp['USA', year],
		allies = cat_vals['Core Allies'],
		rivals = cat_vals['Strategic Rivals'],
		hedgers = cat_vals['Pivotal Hedgers'],
		focal_actors_above_50 = sum(un_sp[intersect(order_actors, rownames(un_sp)), year] >= 0.5, na.rm = TRUE),
		dominant_categories = paste(dominant, collapse = ', '),
		signature = episode_signatures[[year]],
		stringsAsFactors = FALSE,
		row.names = NULL
	)
}))

source_md5 = function(path){
	path = path.expand(path)
	if(!file.exists(path)){ return(NA_character_) }
	unname(tools::md5sum(path))
}
cache_stamp = function(path){
	path = path.expand(path)
	if(!file.exists(path)){ return(NA_character_) }
	format(file.info(path)$mtime, '%Y-%m-%d %H:%M:%S %Z')
}
load(paste0(path_out, 'un_votes_meta.rda'))
load(paste0(path_out, 'trade_meta.rda'))
load(paste0(path_out, 'sanctions_meta.rda'))
fit_source_hash = function(rel){
	load(paste0(path_out, 'snap_fit_', rel, '.rda'))
	if(is.null(fit$americium_build$lame_source_md5)){
		return(NA_character_)
	}
	fit$americium_build$lame_source_md5
}
saved_fit_source_hashes = unique(vapply(relations, fit_source_hash, character(1)))
saved_fit_source_hashes = saved_fit_source_hashes[!is.na(saved_fit_source_hashes) & nzchar(saved_fit_source_hashes)]
source_manifest = data.frame(
	layer = c('UN voting', 'Trade', 'Sanctions', 'Snap-shift package'),
	source = c(
		'UN Digital Library bulk UNGA voting file',
		'IMF IMTS API annual FOB exports',
		'OpenSanctions sanctions targets nested JSON',
		'local lame branch source recorded in saved fits'
	),
	file_or_endpoint = c(
		basename(un_votes_meta$csv_url),
		'IMTS API, indicator XG_FOB_USD',
		basename(sanctions_meta$source),
		'lame source files'
	),
	retrieval_or_build = c(
		paste0('bulk file vintage 2026_02_06; cache ', cache_stamp(paste0(path_in, '2026_02_06_ga_voting.csv'))),
		paste0('current API response cached as trade_long.rda; cache ', cache_stamp(paste0(path_in, 'trade_long.rda'))),
		paste0('latest endpoint cached as opensanctions_targets_nested.json; cache ', cache_stamp(paste0(path_in, 'opensanctions_targets_nested.json'))),
		'fit build time recorded in result objects'
	),
	hash_or_identifier = c(
		source_md5(paste0(path_in, '2026_02_06_ga_voting.csv')),
		source_md5(paste0(path_in, 'trade_long.rda')),
		source_md5(paste0(path_in, 'opensanctions_targets_nested.json')),
		paste(saved_fit_source_hashes, collapse = '; ')
	),
	notes = c(
		paste0('years ', paste(un_votes_meta$year_range, collapse = '--'), '; ES-11/7 present'),
		paste0('years ', paste(trade_meta$year_range, collapse = '--'), '; ', trade_meta$n_obs, ' dyadic observations'),
		paste0('dated records ', sanctions_meta$n_dated, '; source year range ', paste(sanctions_meta$year_range, collapse = '--')),
		paste0('current local source match in regenerated numbers: ', ifelse(all(mcmc_diagnostics$local_source_match), 'yes', 'no'))
	),
	stringsAsFactors = FALSE
)

paper_numbers = list(
	relation_summary = relation_summary,
	relation_draw_summary = relation_draw_summary,
	un_shocks = un_shocks,
	un_top_2025 = un_top_2025,
	un_category_draw_summary = un_category_draw_summary,
	un_denominator_contrasts = un_denominator_contrasts,
	un_raw_vote_diagnostics = un_raw_vote_diagnostics,
	un_leaveout_summary = un_leaveout_summary,
	un_historical_calibration = un_historical_calibration,
	un_country_snap_history = un_country_snap_history,
	un_period_category = un_period_category,
	un_episode_typology = un_episode_typology,
	mcmc_diagnostics = mcmc_diagnostics,
	source_manifest = source_manifest,
	actor_categories = amer_actor_category_table,
	build = amer_build_metadata()
)

optional_results = c(
	configuration = paste0(path_out, 'configuration.rda'),
	rupture_index = paste0(path_out, 'rupture_index.rda'),
	leave_hegemon_out = paste0(path_out, 'leave_hegemon_out.rda'),
	pivotal_hedgers = paste0(path_out, 'pivotal_hedgers.rda'),
	who_snaps_full = paste0(path_out, 'who_snaps_full.rda'),
	un_specification_checks = paste0(path_out, 'un_specification_checks.rda'),
	un_robustness = paste0(path_out, 'un_robustness.rda')
)
for(nm in names(optional_results)){
	if(file.exists(optional_results[[nm]])){
		load(optional_results[[nm]])
		paper_numbers[[nm]] = get(nm)
	}
}

save(paper_numbers, file = paste0(path_out, 'paper_numbers.rda'))
write.csv(relation_summary, paste0(path_out, 'paper_numbers_relation_summary.csv'), row.names = FALSE)
write.csv(relation_draw_summary, paste0(path_out, 'paper_numbers_relation_draw_summary.csv'), row.names = FALSE)
write.csv(un_shocks, paste0(path_out, 'paper_numbers_un_shocks.csv'), row.names = FALSE)
write.csv(un_top_2025, paste0(path_out, 'paper_numbers_un_top_2025.csv'), row.names = FALSE)
write.csv(un_category_draw_summary, paste0(path_out, 'paper_numbers_un_category_draw_summary.csv'), row.names = FALSE)
write.csv(un_denominator_contrasts, paste0(path_out, 'paper_numbers_un_denominator_contrasts.csv'), row.names = FALSE)
write.csv(un_raw_vote_diagnostics, paste0(path_out, 'paper_numbers_un_raw_vote_diagnostics.csv'), row.names = FALSE)
write.csv(un_leaveout_summary, paste0(path_out, 'paper_numbers_un_leaveout_summary.csv'), row.names = FALSE)
write.csv(un_historical_calibration, paste0(path_out, 'paper_numbers_un_historical_calibration.csv'), row.names = FALSE)
write.csv(un_country_snap_history, paste0(path_out, 'paper_numbers_un_country_snap_history.csv'), row.names = FALSE)
write.csv(un_period_category, paste0(path_out, 'paper_numbers_un_period_category.csv'), row.names = FALSE)
write.csv(un_episode_typology, paste0(path_out, 'paper_numbers_un_episode_typology.csv'), row.names = FALSE)
write.csv(mcmc_diagnostics, paste0(path_out, 'paper_numbers_mcmc_diagnostics.csv'), row.names = FALSE)
write.csv(source_manifest, paste0(path_out, 'paper_numbers_source_manifest.csv'), row.names = FALSE)
write.csv(amer_actor_category_table, paste0(path_out, 'paper_actor_categories.csv'), row.names = FALSE)

path_generated_tables = paste0(path_paper, 'generated_tables/')
dir.create(path_generated_tables, recursive = TRUE, showWarnings = FALSE)

tex_escape = function(x){
	x = gsub('\\\\', '\\\\textbackslash{}', x)
	x = gsub('&', '\\\\&', x)
	x = gsub('%', '\\\\%', x)
	x = gsub('_', '\\\\_', x)
	x
}
tex_hash = function(x){
	if(grepl('^[0-9a-f]{32}$', x)){
		return(paste0('\\texttt{\\seqsplit{', x, '}}'))
	}
	if(grepl('R/lame\\.R=', x)){
		lame_hash = sub('.*R/lame\\.R=([0-9a-f]{32}).*', '\\1', x)
		cpp_hash = sub('.*dynamic_latent_factors\\.cpp=([0-9a-f]{32}).*', '\\1', x)
		return(paste0(
			'\\texttt{R/lame.R}=\\texttt{\\seqsplit{', lame_hash, '}}; ',
			'\\texttt{cpp}=\\texttt{\\seqsplit{', cpp_hash, '}}'
		))
	}
	tex_escape(x)
}

country_rows = apply(un_country_snap_history, 1, function(row){
	paste(
		tex_escape(row[['actor']]),
		tex_escape(row[['category']]),
		tex_escape(row[['max_year']]),
		sprintf('%.3f', as.numeric(row[['max_snap_prob']])),
		sprintf('%.3f', as.numeric(row[['snap_prob_2025']])),
		tex_escape(row[['top_snap_years']]),
		sep = ' & '
	)
})
country_table = c(
	'\\begin{sidewaystable}[ht]',
	'\\centering',
	'\\caption{Country-level snap histories for theory-selected focal actors in the UN-voting layer. The final column lists the four highest posterior snap-probability years for each actor.}',
	'\\label{tab:country-snap-history}',
	'\\scriptsize',
	'\\begin{tabularx}{\\textwidth}{@{}l l r r r X@{}}',
	'\\toprule',
	'\\textbf{Actor} & \\textbf{Category} & \\textbf{Max. year} & \\textbf{Max. prob.} & \\textbf{2025 prob.} & \\textbf{Top snap years} \\\\',
	'\\midrule',
	paste0(country_rows, ' \\\\'),
	'\\bottomrule',
	'\\end{tabularx}',
	'\\end{sidewaystable}'
)
writeLines(country_table, paste0(path_generated_tables, 'country_snap_history.tex'))

episode_rows = apply(un_episode_typology, 1, function(row){
	paste(
		tex_escape(row[['year']]),
		sprintf('%.3f', as.numeric(row[['focal_index']])),
		sprintf('%.3f', as.numeric(row[['us_snap']])),
		sprintf('%.3f', as.numeric(row[['allies']])),
		sprintf('%.3f', as.numeric(row[['rivals']])),
		sprintf('%.3f', as.numeric(row[['hedgers']])),
		tex_escape(row[['focal_actors_above_50']]),
		tex_escape(row[['signature']]),
		sep = ' & '
	)
})
episode_table = c(
	'\\begin{table}[ht]',
	'\\centering',
	'\\caption{Selected diplomatic snap episodes in the UN-voting layer. Category entries are posterior mean snap probabilities for the focal actors in that category.}',
	'\\label{tab:episode-typology}',
	'\\scriptsize',
	'\\begin{tabularx}{\\textwidth}{@{}r r r r r r r X@{}}',
	'\\toprule',
	'\\textbf{Year} & \\textbf{Index} & \\textbf{U.S.} & \\textbf{Allies} & \\textbf{Rivals} & \\textbf{Hedgers} & \\textbf{$>0.5$} & \\textbf{Signature} \\\\',
	'\\midrule',
	paste0(episode_rows, ' \\\\'),
	'\\bottomrule',
	'\\end{tabularx}',
	'\\end{table}'
)
writeLines(episode_table, paste0(path_generated_tables, 'episode_typology.tex'))

raw_rows = apply(un_raw_vote_diagnostics, 1, function(row){
	paste(
		tex_escape(row[['metric']]),
		sprintf('%.3f', as.numeric(row[['value_2024']])),
		sprintf('%.3f', as.numeric(row[['value_2025']])),
		sprintf('%+.3f', as.numeric(row[['change']])),
		tex_escape(row[['interpretation']]),
		sep = ' & '
	)
})
raw_table = c(
	'\\begin{table}[ht]',
	'\\centering',
	'\\caption{Raw UN-voting diagnostics around 2025. Entries are descriptive annual-agenda summaries, not model estimates and not same-resolution before--after comparisons. Agreement is Lijphart three-category agreement on yes, no, and abstain votes.}',
	'\\label{tab:raw-un-diagnostics}',
	'\\scriptsize',
	'\\begin{tabularx}{\\textwidth}{@{}X r r r X@{}}',
	'\\toprule',
	'\\textbf{Raw diagnostic} & \\textbf{2024} & \\textbf{2025} & \\textbf{Change} & \\textbf{Interpretation} \\\\',
	'\\midrule',
	paste0(raw_rows, ' \\\\'),
	'\\bottomrule',
	'\\end{tabularx}',
	'\\end{table}'
)
writeLines(raw_table, paste0(path_generated_tables, 'raw_un_diagnostics.tex'))

leaveout_wide = data.frame(
	contrast = c(
		'Full focal index',
		'Excluding United States',
		'Excluding United States, Russia, and China',
		'Excluding United States and all strategic rivals'
	),
	value_2024 = c(
		un_leaveout_summary$full[un_leaveout_summary$year == 2024],
		un_leaveout_summary$excl_us[un_leaveout_summary$year == 2024],
		un_leaveout_summary$excl_us_rus_chn[un_leaveout_summary$year == 2024],
		un_leaveout_summary$excl_us_rivals[un_leaveout_summary$year == 2024]
	),
	value_2025 = c(
		un_leaveout_summary$full[un_leaveout_summary$year == 2025],
		un_leaveout_summary$excl_us[un_leaveout_summary$year == 2025],
		un_leaveout_summary$excl_us_rus_chn[un_leaveout_summary$year == 2025],
		un_leaveout_summary$excl_us_rivals[un_leaveout_summary$year == 2025]
	),
	stringsAsFactors = FALSE
)
leaveout_wide$change = leaveout_wide$value_2025 - leaveout_wide$value_2024
leaveout_rows = apply(leaveout_wide, 1, function(row){
	paste(
		tex_escape(row[['contrast']]),
		sprintf('%.3f', as.numeric(row[['value_2024']])),
		sprintf('%.3f', as.numeric(row[['value_2025']])),
		sprintf('%+.3f', as.numeric(row[['change']])),
		sep = ' & '
	)
})
leaveout_table = c(
	'\\begin{table}[ht]',
	'\\centering',
	'\\caption{Leave-out focal rupture indices in the main UN-voting fit. Entries recompute the reported index over different actor sets using the same posterior snap probabilities; they are aggregation checks, not separate no-U.S. refits.}',
	'\\label{tab:leaveout-index}',
	'\\small',
	'\\begin{tabularx}{\\textwidth}{@{}X r r r@{}}',
	'\\toprule',
	'\\textbf{Index denominator} & \\textbf{2024} & \\textbf{2025} & \\textbf{Change} \\\\',
	'\\midrule',
	paste0(leaveout_rows, ' \\\\'),
	'\\bottomrule',
	'\\end{tabularx}',
	'\\end{table}'
)
writeLines(leaveout_table, paste0(path_generated_tables, 'leaveout_index.tex'))

relation_labels = c(
	un_voting = 'UN voting',
	trade = 'Trade',
	sanctions = 'Sanctions'
)
index_side_labels = c(
	sender = 'State/sender',
	receiver = 'Receiver/target',
	combined = 'Sender or receiver'
)
uncertainty_rows = apply(relation_draw_summary, 1, function(row){
	paste(
		tex_escape(relation_labels[[row[['relation']]]]),
		tex_escape(index_side_labels[[row[['index_side']]]]),
		sprintf('%.3f', as.numeric(row[['index_2025_mean']])),
		paste0(
			'[',
			sprintf('%.3f', as.numeric(row[['index_2025_q025']])),
			', ',
			sprintf('%.3f', as.numeric(row[['index_2025_q975']])),
			']'
		),
		sprintf('%.2f', as.numeric(row[['pr_rank_1']])),
		sprintf('%.2f', as.numeric(row[['pr_rank_le_3']])),
		sep = ' & '
	)
})
uncertainty_table = c(
	'\\begin{table}[ht]',
	'\\centering',
	'\\caption{Draw-level uncertainty for the 2025 focal rupture indices. Entries average retained snap-indicator draws within the theory-selected focal actor set before summarizing the index. Rank probabilities compare 2025 with all eligible transition years in the same relation and index side.}',
	'\\label{tab:snap-index-uncertainty}',
	'\\small',
	'\\begin{tabularx}{\\textwidth}{@{}l l r r r r@{}}',
	'\\toprule',
	'\\textbf{Relation} & \\textbf{Index side} & \\textbf{Mean} & \\textbf{95\\% interval} & \\textbf{Pr(rank 1)} & \\textbf{Pr(rank $\\leq$ 3)} \\\\',
	'\\midrule',
	paste0(uncertainty_rows, ' \\\\'),
	'\\bottomrule',
	'\\end{tabularx}',
	'\\end{table}'
)
writeLines(uncertainty_table, paste0(path_generated_tables, 'snap_index_uncertainty.tex'))

diag_labels = c(
	un_voting = 'UN voting',
	trade = 'Trade',
	sanctions = 'Sanctions'
)
diag_rows = apply(mcmc_diagnostics, 1, function(row){
	paste(
		diag_labels[[row[['relation']]]],
		tex_escape(row[['family']]),
		tex_escape(paste0(row[['nscan']], '/', row[['burn']], '/', row[['thin']])),
		tex_escape(row[['retained_draws']]),
		sprintf('%.3f (%.3f)', as.numeric(row[['pi_mean']]), as.numeric(row[['pi_sd']])),
		sprintf('%.3f (%.3f)', as.numeric(row[['rho_uv_mean']]), as.numeric(row[['rho_uv_sd']])),
		sprintf('%.3f (%.3f)', as.numeric(row[['sigma_uv_mean']]), as.numeric(row[['sigma_uv_sd']])),
		sprintf('%.0f', as.numeric(row[['min_ess']])),
		tex_escape(ifelse(as.logical(row[['source_hash_recorded']]), 'yes', 'no')),
		tex_escape(ifelse(as.logical(row[['local_source_match']]), 'yes', 'no')),
		sep = ' & '
	)
})
diag_table = c(
	'\\begin{table}[ht]',
	'\\centering',
	'\\caption{Scalar MCMC diagnostics for the three main snap-shift fits. The scan column reports total scans, burn-in, and thinning interval; retained draws are posterior draws after burn-in and thinning. Parentheses report posterior standard deviations. ESS is the minimum effective sample size across the saved scalar chains for $\\pi$, $\\rho_{uv}$, and $\\sigma_u$. Recorded source indicates that the saved fit contains the snap-shift source hash; local match reports whether the currently mutable local \\pkg{lame} branch matches that saved hash.}',
	'\\label{tab:mcmc-diagnostics}',
	'\\scriptsize',
	'\\resizebox{\\textwidth}{!}{%',
	'\\begin{tabular}{@{}l l l r r r r r l l@{}}',
	'\\toprule',
	'\\textbf{Relation} & \\textbf{Family} & \\textbf{Scans/burn/thin} & \\textbf{Draws} & $\\bar{\\pi}$ & $\\rho_{uv}$ & $\\sigma_u$ & \\textbf{Min. ESS} & \\textbf{Recorded source} & \\textbf{Local match} \\\\',
	'\\midrule',
	paste0(diag_rows, ' \\\\'),
	'\\bottomrule',
	'\\end{tabular}',
	'}',
	'\\end{table}'
)
writeLines(diag_table, paste0(path_generated_tables, 'mcmc_diagnostics.tex'))
