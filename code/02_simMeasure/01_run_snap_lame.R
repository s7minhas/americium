# fit snap models

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
run_lame_version = as.character(utils::packageVersion('lame'))
run_lame_source_md5 = amer_lame_source_md5()
set.seed(seed)

load(paste0(path_in, 'arr_list_lame.rda'))
status_file = paste0(path_out, 'snap_fit_status.log')
write_status = function(...){
	write(
		paste(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), paste(..., collapse = ' ')),
		file = status_file,
		append = TRUE
	)
}

# relation settings
relation_config = list(
	un_voting = list(
		transform = function(m) m,
		family = 'normal',
		symmetric = TRUE
	),
	trade = list(
		transform = function(m){ m[m <= 0] = NA; log(m) },
		family = 'normal',
		symmetric = TRUE
	),
	sanctions = list(
		transform = function(m) m,
		family = 'poisson',
		symmetric = FALSE
	)
)

requested_relations = Sys.getenv('RELATIONS', unset = '')
if(nchar(requested_relations) > 0){
	requested_relations = trimws(strsplit(requested_relations, ',')[[1]])
	relation_config = relation_config[intersect(requested_relations, names(relation_config))]
}
if(length(relation_config) == 0){ stop('No valid relations selected.') }

snap_nscan = as.integer(Sys.getenv('SNAP_NSCAN', unset = '6000'))
snap_burn = as.integer(Sys.getenv('SNAP_BURN', unset = '3000'))
snap_odens = as.integer(Sys.getenv('SNAP_ODENS', unset = '10'))
snap_output_suffix = Sys.getenv('SNAP_OUTPUT_SUFFIX', unset = '')
if(!grepl('^[A-Za-z0-9_\\-]*$', snap_output_suffix)){
	stop('SNAP_OUTPUT_SUFFIX may contain only letters, numbers, underscores, and hyphens.')
}
if(any(!is.finite(c(snap_nscan, snap_burn, snap_odens))) ||
		snap_nscan <= 0 || snap_burn < 0 || snap_odens <= 0 || snap_burn >= snap_nscan){
	stop('Invalid SNAP_NSCAN, SNAP_BURN, or SNAP_ODENS settings.')
}

# fit and save each relation
for(rel in names(relation_config)){
	write_status('fitting', rel, 'nscan', snap_nscan, 'burn', snap_burn, 'odens', snap_odens)
	cfg = relation_config[[rel]]
	y = lapply(arr_list_lame[[rel]], cfg$transform)
	fit = lame(
		y,
		R = 2,
		family = cfg$family,
		symmetric = cfg$symmetric,
		dynamic_uv = TRUE,
		dynamic_uv_kind = 'snap',
		prior = amer_snap_prior,
		nscan = snap_nscan,
		burn = snap_burn,
		odens = snap_odens,
		keep_snap_draws = 'draws',
		seed = seed,
		gof = FALSE,
		verbose = FALSE
	)
	fit$americium_build = amer_build_metadata(list(
		relation = rel,
		lame_version = run_lame_version,
		lame_source_md5 = run_lame_source_md5,
		snap_nscan = snap_nscan,
		snap_burn = snap_burn,
		snap_odens = snap_odens,
		snap_output_suffix = snap_output_suffix
	))
	out_file = paste0(path_out, 'snap_fit_', rel, snap_output_suffix, '.rda')
	save(fit, file = out_file)
	write_status('saved', rel, out_file)
}
