# terminal-year stress tests

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

load_pkg(c('lame', 'data.table'))
set.seed(seed)

load(paste0(path_in, 'arr_list_lame.rda'))
load(paste0(path_out, 'snap_fit_un_voting.rda'))
main_fit = fit

un = arr_list_lame$un_voting
all_years = names(un)

observed_actors_matrix = function(net_list){
	obs = sapply(net_list, function(m){
		rowSums(!is.na(m), na.rm = TRUE) + colSums(!is.na(m), na.rm = TRUE) > 0
	})
	obs
}

actors_observed_in_year = function(fit_obj, yr){
	y = fit_obj$Y[, , yr]
	rownames(y)[rowSums(!is.na(y), na.rm = TRUE) + colSums(!is.na(y), na.rm = TRUE) > 0]
}

subset_net = function(net_list, actors){
	lapply(net_list, function(m){ m[actors, actors] })
}

fit_subset = function(net_list, label){
	out_file = paste0(path_out, 'snap_stress_fit_', label, '.rda')
	if(file.exists(out_file)){
		load(out_file)
		return(fit)
	}
	fit = lame(
		net_list,
		R = 2,
		family = 'normal',
		symmetric = TRUE,
		dynamic_uv = TRUE,
		dynamic_uv_kind = 'snap',
		prior = amer_snap_prior,
		nscan = 1500,
		burn = 750,
		odens = 5,
		keep_snap_draws = 'draws',
		seed = seed,
		gof = FALSE,
		verbose = FALSE
	)
	fit$americium_build = amer_build_metadata(list(
		relation = 'un_voting',
		specification = paste0('stress_', label),
		stress_nscan = 1500,
		stress_burn = 750,
		stress_odens = 5
	))
	save(fit, file = out_file)
	fit
}

summarize_terminal = function(fit_obj, label, terminal_year){
	obs = actors_observed_in_year(fit_obj, terminal_year)
	focal_obs = intersect(amer_focal_actors, obs)
	sp = fit_obj$snap_prob
	data.frame(
		test = label,
		terminal_year = terminal_year,
		n_actors = nrow(sp),
		n_observed_terminal = length(obs),
		focal_observed_terminal = length(focal_obs),
		all_mean = mean(sp[obs, terminal_year], na.rm = TRUE),
		all_ge50 = mean(sp[obs, terminal_year] >= 0.5, na.rm = TRUE),
		all_ge99 = mean(sp[obs, terminal_year] >= 0.99, na.rm = TRUE),
		focal_mean = mean(sp[focal_obs, terminal_year], na.rm = TRUE),
		focal_ge50 = mean(sp[focal_obs, terminal_year] >= 0.5, na.rm = TRUE),
		focal_ge99 = mean(sp[focal_obs, terminal_year] >= 0.99, na.rm = TRUE),
		stringsAsFactors = FALSE
	)
}

# use focal and balanced actors
obs_mat = observed_actors_matrix(un)
balanced_actors = rownames(obs_mat)[rowSums(obs_mat) == ncol(obs_mat)]
focal_available = intersect(amer_focal_actors, rownames(un[[1]]))
nonfocal_balanced = sort(setdiff(balanced_actors, focal_available))
diagnostic_actors = unique(c(focal_available, head(nonfocal_balanced, 45)))

un_diag = subset_net(un, diagnostic_actors)
un_diag_2025_copy_2024 = un_diag
un_diag_2025_copy_2024[['2025']] = un_diag_2025_copy_2024[['2024']]

tests = list(
	main_full_2025 = summarize_terminal(main_fit, 'main_full_2025', '2025')
)

stress_specs = list(
	subset_actual_2025 = list(
		y = un_diag[all_years <= '2025'],
		terminal = '2025'
	),
	subset_copy_2024_as_2025 = list(
		y = un_diag_2025_copy_2024[all_years <= '2025'],
		terminal = '2025'
	),
	subset_terminal_2024 = list(
		y = un_diag[all_years <= '2024'],
		terminal = '2024'
	),
	subset_terminal_2023 = list(
		y = un_diag[all_years <= '2023'],
		terminal = '2023'
	)
)

for(nm in names(stress_specs)){
	spec = stress_specs[[nm]]
	fit_stress = fit_subset(spec$y, nm)
	tests[[nm]] = summarize_terminal(fit_stress, nm, spec$terminal)
}

snap_stress_tests = list(
	diagnostic_actors = diagnostic_actors,
	stress_chain = list(nscan = 1500, burn = 750, odens = 5),
	summary = rbindlist(tests, fill = TRUE),
	note = paste(
		'Short MCMC diagnostic fits on focal actors plus 45 balanced non-focal actors.',
		'Use for model-pathology triage, not final inference.'
	)
)

save(snap_stress_tests, file = paste0(path_out, 'snap_stress_tests.rda'))
write.csv(snap_stress_tests$summary,
	file = paste0(path_out, 'snap_stress_tests_summary.csv'),
	row.names = FALSE)

tex_escape = function(x){
	x = gsub('\\\\', '\\\\textbackslash{}', x)
	x = gsub('_', '\\\\_', x)
	x
}

stress_tex = copy(snap_stress_tests$summary)
stress_tex$test = c(
	'Main full panel, actual 2025',
	'Diagnostic subset, actual 2025',
	'Diagnostic subset, 2024 copied as 2025',
	'Diagnostic subset, terminal 2024',
	'Diagnostic subset, terminal 2023'
)
stress_tex = stress_tex[, .(
	Test = test,
	`Terminal year` = terminal_year,
	`Actors` = n_actors,
	`Observed focal` = focal_observed_terminal,
	`Focal mean` = sprintf('%.3f', focal_mean),
	`Focal > .50` = sprintf('%.3f', focal_ge50),
	`Focal > .99` = sprintf('%.3f', focal_ge99),
	`All mean` = sprintf('%.3f', all_mean)
)]

stress_lines = c(
	'\\begin{table}[!htbp]',
	'\\centering',
	'\\begingroup',
	'\\scriptsize',
	'\\caption{Stress tests for terminal-year snap saturation in the UN-voting layer.}',
	'\\label{tab:snap-stress-tests}',
	'\\begin{tabularx}{\\textwidth}{@{}>{\\raggedright\\arraybackslash}Xrrrrrrr@{}}',
	'\\toprule',
	'Test & Year & Actors & Focal obs. & Focal mean & Focal $>.50$ & Focal $>.99$ & All mean \\\\',
	'\\midrule'
)
for(i in seq_len(nrow(stress_tex))){
	stress_lines = c(stress_lines, paste(
		tex_escape(stress_tex$Test[i]),
		stress_tex$`Terminal year`[i],
		stress_tex$Actors[i],
		stress_tex$`Observed focal`[i],
		stress_tex$`Focal mean`[i],
		stress_tex$`Focal > .50`[i],
		stress_tex$`Focal > .99`[i],
		stress_tex$`All mean`[i],
		sep = ' & '
	) |> paste0(' \\\\'))
}
stress_lines = c(
	stress_lines,
	'\\bottomrule',
	'\\end{tabularx}',
	'\\endgroup',
	'\\end{table}'
)
writeLines(stress_lines, paste0(path_paper, 'generated_tables/snap_stress_tests.tex'))
