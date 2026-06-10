# rolling prediction check

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

load_pkg(c('lame', 'ggplot2', 'ggrepel'))
set.seed(seed)

load(paste0(path_in, 'arr_list_lame.rda'))
un = arr_list_lame$un_voting
years = as.integer(names(un))
targets = c(2024, 2025)
high_surprise_threshold = 2

tex_escape = function(x){
	x = gsub('\\\\', '\\\\textbackslash{}', x)
	x = gsub('&', '\\\\&', x)
	x = gsub('%', '\\\\%', x)
	x = gsub('_', '\\\\_', x)
	x
}

fit_ar = function(train_years){
	train_label = max(train_years)
	out_file = paste0(path_out, 'oos_fit_un_voting_train_through_', train_label, '.rda')
	if(file.exists(out_file)){
		load(out_file)
		return(fit)
	}
	fit = lame(
		un[as.character(train_years)],
		R = 2,
		family = 'normal',
		symmetric = TRUE,
		dynamic_uv = TRUE,
		dynamic_uv_kind = 'ar1',
		nscan = 4000,
		burn = 2000,
		odens = 10,
		seed = seed,
		gof = FALSE,
		verbose = FALSE
	)
	fit$americium_build = amer_build_metadata(list(
		relation = 'un_voting',
		specification = 'rolling_oos_ar1',
		train_through = train_label,
		oos_nscan = 4000,
		oos_burn = 2000,
		oos_odens = 10
	))
	save(fit, file = out_file)
	fit
}

forecast_eta = function(u_prev, rho, add_eff){
	u_hat = rho * u_prev
	outer(add_eff, add_eff, '+') + u_hat %*% t(u_hat)
}

state_logdens = function(obs, eta, s2){
	ld = dnorm(obs, mean = eta, sd = sqrt(s2), log = TRUE)
	diag(ld) = NA
	rowMeans(ld, na.rm = TRUE)
}

category_mean = function(x, actors){
	keep = intersect(actors, names(x))
	if(length(keep) == 0){ return(NA_real_) }
	mean(x[keep], na.rm = TRUE)
}

run_target = function(target_year){
	train_years = years[years < target_year]
	fit = fit_ar(train_years)
	u = fit$U
	fit_years = dimnames(u)[[3]]
	actors = rownames(u[, , 1])
	add_eff = setNames(as.numeric(fit$APM), actors)
	rho = mean(fit$rho_uv)
	s2 = mean(fit$VC[, 've'])

	null_ld = matrix(
		NA_real_,
		length(actors),
		length(fit_years) - 1,
		dimnames = list(actors, fit_years[-1])
	)
	for(t in 2:length(fit_years)){
		eta = forecast_eta(u[, , t-1], rho, add_eff)
		obs = un[[fit_years[t]]][actors, actors]
		null_ld[, t-1] = state_logdens(obs, eta, s2)
	}

	eta_target = forecast_eta(u[, , length(fit_years)], rho, add_eff)
	obs_target = un[[as.character(target_year)]][actors, actors]
	ld_target = state_logdens(obs_target, eta_target, s2)

	state_surprise = vapply(actors, function(i){
		hist = null_ld[i, ]
		hist = hist[is.finite(hist)]
		if(length(hist) < 5 || !is.finite(ld_target[i])){ return(NA_real_) }
		(median(hist) - ld_target[i]) / (mad(hist) + 1e-8)
	}, numeric(1))

	observed = names(state_surprise)[is.finite(state_surprise)]
	focal = intersect(amer_focal_actors, observed)
	agg_hist = colMeans(null_ld[observed, , drop = FALSE], na.rm = TRUE)
	agg_target = mean(ld_target[observed], na.rm = TRUE)
	focal_hist = colMeans(null_ld[focal, , drop = FALSE], na.rm = TRUE)
	focal_target = mean(ld_target[focal], na.rm = TRUE)

	summary = data.frame(
		target_year = target_year,
		train_start = min(train_years),
		train_end = max(train_years),
		n_train_years = length(train_years),
		n_observed_states = length(observed),
		n_focal_states = length(focal),
		all_mean_surprise = mean(state_surprise[observed], na.rm = TRUE),
		all_median_surprise = median(state_surprise[observed], na.rm = TRUE),
		focal_mean_surprise = mean(state_surprise[focal], na.rm = TRUE),
		focal_median_surprise = median(state_surprise[focal], na.rm = TRUE),
		us_surprise = state_surprise['USA'],
		allies_mean_surprise = category_mean(state_surprise, amer_actor_cats$allies),
		rivals_mean_surprise = category_mean(state_surprise, amer_actor_cats$rivals),
		hedgers_mean_surprise = category_mean(state_surprise, amer_actor_cats$hedgers),
		focal_share_high_surprise = mean(state_surprise[focal] > high_surprise_threshold, na.rm = TRUE),
		all_target_logdens = agg_target,
		all_target_worst_rank = rank(c(agg_target, agg_hist), ties.method = 'min')[1],
		all_comparison_years = length(agg_hist) + 1,
		all_target_worse_percentile = mean(agg_hist >= agg_target, na.rm = TRUE),
		focal_target_logdens = focal_target,
		focal_target_worst_rank = rank(c(focal_target, focal_hist), ties.method = 'min')[1],
		focal_comparison_years = length(focal_hist) + 1,
		focal_target_worse_percentile = mean(focal_hist >= focal_target, na.rm = TRUE),
		stringsAsFactors = FALSE,
		row.names = NULL
	)

	actor_summary = data.frame(
		target_year = target_year,
		actor = actors,
		category = amer_actor_type(actors),
		focal_actor = actors %in% amer_focal_actors,
		log_density = as.numeric(ld_target[actors]),
		surprise = as.numeric(state_surprise[actors]),
		stringsAsFactors = FALSE
	)
	actor_summary$category[is.na(actor_summary$category)] = 'Other States'

	list(
		target_year = target_year,
		train_years = train_years,
		fit_file = paste0(path_out, 'oos_fit_un_voting_train_through_', max(train_years), '.rda'),
		null_ld = null_ld,
		ld_target = ld_target,
		state_surprise = state_surprise,
		summary = summary,
		actor_summary = actor_summary
	)
}

rolling = lapply(targets, run_target)
names(rolling) = paste0('predict_', targets)
oos_summary = do.call(rbind, lapply(rolling, `[[`, 'summary'))
oos_actor_summary = do.call(rbind, lapply(rolling, `[[`, 'actor_summary'))

oos_surprise = list(
	rolling = rolling,
	summary = oos_summary,
	actor_summary = oos_actor_summary,
	high_surprise_threshold = high_surprise_threshold,
	build = amer_build_metadata(list(
		relation = 'un_voting',
		specification = 'rolling_oos_ar1',
		target_years = paste(targets, collapse = ',')
	))
)
save(oos_surprise, file = paste0(path_out, 'oos_surprise.rda'))
write.csv(oos_summary, paste0(path_out, 'oos_surprise_summary.csv'), row.names = FALSE)

plot_df = oos_actor_summary[is.finite(oos_actor_summary$surprise), ]
plot_df$target_label = factor(
	paste0('predict ', plot_df$target_year),
	levels = paste0('predict ', targets)
)
plot_df$category = factor(plot_df$category, levels = c(amer_category_labels, 'Other States'))
highlight_df = plot_df[plot_df$focal_actor | plot_df$actor == 'USA', ]
label_df = highlight_df[highlight_df$actor %in% c('USA','DEU','FRA','GBR','JPN','KOR','RUS','CHN','IND','TUR','BRA'), ]

fig_oos = ggplot(plot_df, aes(x = target_label, y = surprise)) +
	geom_hline(yintercept = 0, color = 'grey70', linewidth = 0.3) +
	geom_hline(yintercept = high_surprise_threshold, linetype = 'dashed', color = 'grey55', linewidth = 0.3) +
	geom_violin(fill = 'grey90', color = 'grey65', linewidth = 0.3, width = 0.7, trim = FALSE) +
	geom_jitter(data = plot_df[!plot_df$focal_actor, ], width = 0.08, height = 0, size = 0.5, alpha = 0.25, color = 'grey45') +
	geom_point(data = highlight_df, aes(color = category), position = position_jitter(width = 0.06, height = 0), size = 1.5, alpha = 0.9) +
	geom_point(data = plot_df[plot_df$actor == 'USA', ], shape = 23, fill = 'white', color = 'black', size = 2.7, stroke = 0.7) +
	ggrepel::geom_text_repel(
		data = label_df,
		aes(label = actor, color = category),
		size = 2.4,
		min.segment.length = 0,
		max.overlaps = 30,
		show.legend = FALSE
	) +
	scale_color_manual(values = c(amer_category_colors, 'Other States' = 'grey45'), name = 'Focal actor category') +
	labs(x = NULL, y = 'Standardized Forecast Surprise (Higher = Worse)') +
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.ticks = element_blank(),
		legend.position = 'top',
		legend.title = element_text(size = 9)
	)
ggsave(paste0(path_graphics, 'oos_prediction_surprise.pdf'), fig_oos, width = 6.6, height = 4.4)
ggsave(paste0(path_graphics, 'oos_prediction_surprise.png'), fig_oos, width = 6.6, height = 4.4, dpi = 600)

fmt = function(x){ sprintf('%.2f', as.numeric(x)) }
pct = function(x){ paste0(sprintf('%.0f', 100 * as.numeric(x)), '\\%') }
oos_rows = apply(oos_summary, 1, function(row){
	paste(
		tex_escape(paste0(row[['train_start']], '--', row[['train_end']], ' -> ', row[['target_year']])),
		row[['n_observed_states']],
		row[['n_focal_states']],
		fmt(row[['all_median_surprise']]),
		fmt(row[['focal_median_surprise']]),
		fmt(row[['us_surprise']]),
		fmt(row[['allies_mean_surprise']]),
		fmt(row[['rivals_mean_surprise']]),
		fmt(row[['hedgers_mean_surprise']]),
		pct(row[['focal_share_high_surprise']]),
		pct(row[['focal_target_worse_percentile']]),
		sep = ' & '
	)
})
oos_table = c(
	'\\begin{table}[!htbp]',
	'\\centering',
	"\\caption{Rolling one-year-ahead prediction diagnostics for the UN-voting layer. Each row fits a smooth AR(1) latent-space model through the training end year and forecasts the next observed annual network. Surprise is the per-state loss in predictive log density relative to that state's in-sample one-step forecast history, scaled by the in-sample median absolute deviation. Higher values mean worse prediction.}",
	'\\label{tab:oos-prediction}',
	'\\scriptsize',
	'\\resizebox{\\textwidth}{!}{%',
	'\\begin{tabular}{@{}l r r r r r r r r r r@{}}',
	'\\toprule',
	'\\textbf{Train $\\rightarrow$ target} & \\textbf{States} & \\textbf{Focal} & \\textbf{All med.} & \\textbf{Focal med.} & \\textbf{U.S.} & \\textbf{Allies} & \\textbf{Rivals} & \\textbf{Hedgers} & \\textbf{Focal $>2$} & \\textbf{Focal worse pct.} \\\\',
	'\\midrule',
	paste0(oos_rows, ' \\\\'),
	'\\bottomrule',
	'\\end{tabular}',
	'}',
	'\\end{table}'
)
dir.create(paste0(path_paper, 'generated_tables'), recursive = TRUE, showWarnings = FALSE)
writeLines(oos_table, paste0(path_paper, 'generated_tables/oos_prediction.tex'))
