# replay and verify the trade model chains

helper_file = if (sys.nframe() > 0L && !is.null(sys.frame(1)$ofile)) {
	sys.frame(1)$ofile
} else {
	sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
}

source(
	file.path(dirname(normalizePath(helper_file)), "00_data_helpers.R"),
	local = TRUE
)

args = commandArgs(trailingOnly = TRUE)
argument = sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
root = dirname(dirname(normalizePath(argument)))
source(file.path(root, "code", "01_observed_votes.R"))
source(file.path(root, "code", "00_sampler.R"))
source(file.path(root, "code", "10_models.R"))
replay_trade_space = function(specification = "current_baseline", workers = 2L) {
	stopifnot(
		specification %in% c(
			"current_baseline",
			"both_directions",
			"fresh_terminal",
			"twelve_month_terminal",
			"partner_reallocation_terminal"
		),
		is.finite(workers),
		workers >= 1L,
		as.character(packageVersion("netify")) == "1.5.3"
	)
	output = file.path(root, "results", "models", "trade_space_replay", specification)
	dir.create(output, recursive = TRUE, showWarnings = FALSE)
	context = model_load(file.path(root, "code"))
	inputs = readRDS(file.path(root, "results/models/inputs.rds"))$trade
	starts = readRDS(file.path(root, "data/model_initial_states.rds"))$trade
	plan = read.csv(file.path(root, "data/model_chains.csv"))
	plan = plan[plan$domain == "trade" & plan$specification == specification, ]
	stopifnot(nrow(plan) == 2L, !anyDuplicated(plan$seed))
	y = inputs[[specification]]
	prior = model_priors()$trade
	jobs = parallel::mclapply(
		seq_len(nrow(plan)),
		function(i) {
			job = plan[i, ]
			initial = starts[[specification]]$chains[[as.character(job$seed)]]
			shifts = starts[[specification]]$shifts
			stopifnot(
				identical(digest::digest(y, algo = "sha256"), job$input_sha256),
				identical(digest::digest(initial, algo = "sha256"), job$initial_state_sha256),
				identical(digest::digest(shifts, algo = "sha256"), job$shifts_sha256)
			)
			fit = fit_model(
				y,
				paste(specification, job$seed, sep = "_"),
				job$seed,
				prior,
				context,
				output,
				nscan = job$sampling_sweeps,
				burn = job$warmup_sweeps,
				odens = job$thin,
				initial_state = initial,
				frozen_shifts = shifts,
				indicator_every = job$indicator_every,
				temporal_tau = job$temporal_tau
			)
			stopifnot(isTRUE(fit$complete))
			data.frame(
				domain = "trade",
				specification = specification,
				seed = job$seed,
				complete = TRUE,
				fit_file = fit$revision_build$fit_file,
				checkpoint_file = fit$revision_build$checkpoint_file
			)
		},
		mc.cores = min(workers, nrow(plan)),
		mc.preschedule = FALSE
	)
	stopifnot(all(vapply(jobs, is.data.frame, TRUE)))
	manifest = do.call(rbind, jobs)
	write.csv(manifest, file.path(output, "fit_manifest.csv"), row.names = FALSE)
	invisible(manifest)
}

verify_trade_space = function(specification) {
	output = file.path(root, "results/models/trade_space_replay", specification)
	manifest = file.path(output, "fit_manifest.csv")
	if (!file.exists(manifest)) {
		runs = read.csv(file.path(root, "results/models/trade_replay/fit_manifest.csv"))
		runs = runs[runs$domain == "trade" & runs$specification == specification, ]
		stopifnot(nrow(runs) == 2L, all(runs$complete))
		dir.create(output, recursive = TRUE, showWarnings = FALSE)
		write.csv(runs, manifest, row.names = FALSE)
	}
	runs = read.csv(manifest)
	plan = read.csv(file.path(root, "data/model_chains.csv"))
	plan = plan[plan$domain == "trade" & plan$specification == specification, ]
	stopifnot(
		nrow(runs) == 2L,
		all(runs$complete),
		setequal(runs$seed, plan$seed)
	)
	original = readRDS(file.path(root, "results/models/trade/chain_features.rds"))[[specification]]
	backend = summary_backend(root)
	comparisons = list()
	recovered = list()
	fits = list()
	for (i in seq_len(nrow(runs))) {
		f = readRDS(runs$fit_file[i])
		job = plan[plan$seed == runs$seed[i], ]
		cfg = f$revision_build$config
		stopifnot(
			isTRUE(f$complete),
			cfg$seed == job$seed,
			cfg$nscan == job$sampling_sweeps,
			cfg$burn == job$warmup_sweeps,
			cfg$odens == job$thin,
			cfg$indicator_every == job$indicator_every,
			cfg$temporal_tau == job$temporal_tau,
			identical(cfg$prior, model_priors()$trade),
			identical(digest::digest(f$Y, algo = "sha256"), job$input_sha256),
			identical(
				digest::digest(cfg$initial_state, algo = "sha256"),
				job$initial_state_sha256
			),
			identical(
				digest::digest(cfg$frozen_shifts, algo = "sha256"),
				job$shifts_sha256
			),
			identical(unname(tools::md5sum(names(cfg$code))), unname(cfg$code))
		)
		expected = original[[which(vapply(original, function(x) x$seed == job$seed, TRUE))]]
		actual = summarize_chain(f, "trade", backend)
		stopifnot(setequal(names(expected), names(actual)))
		for (field in names(expected)) {
			equal = isTRUE(all.equal(actual[[field]], expected[[field]], tolerance = 1e-12))
			error = if (is.numeric(actual[[field]]) && identical(dim(actual[[field]]), dim(expected[[field]]))) {
				max(abs(actual[[field]] - expected[[field]]), na.rm = TRUE)
			} else {
				NA_real_
			}
			comparisons[[length(comparisons) + 1L]] = data.frame(
				seed = job$seed,
				field = field,
				matches_release = equal,
				identical_to_release = identical(actual[[field]], expected[[field]]),
				maximum_absolute_difference = error
			)
		}
		recovered[[i]] = actual
		fits[[i]] = list(
			seed = job$seed,
			path = runs$fit_file[i],
			sha256 = digest::digest(file = runs$fit_file[i], algo = "sha256"),
			checkpoint = runs$checkpoint_file[i],
			draws = dim(f$U_draws)[4]
		)
		rm(f, actual)
		gc(FALSE)
	}
	comparisons = do.call(rbind, comparisons)
	write.csv(
		comparisons,
		file.path(output, "release_comparison.csv"),
		row.names = FALSE
	)
	diagnostic = summarize_features(recovered, specification)
	passed = all(comparisons$matches_release) && diagnostic$complete
	jsonlite::write_json(
		list(
			passed = passed,
			specification = specification,
			fits = fits,
			tolerance = 1e-12,
			release_diagnostics_pass = diagnostic$complete,
			release_features_sha256 = digest::digest(
				file = file.path(root, "results/models/trade/chain_features.rds"),
				algo = "sha256"
			),
			source_hashes = setNames(
				lapply(
					c(
						"code/00_sampler.R",
						"code/sampler.cpp",
						"code/10_models.R",
						"data/model_chains.csv",
						"data/model_initial_states.rds"
					),
					function(p) digest::digest(file = file.path(root, p), algo = "sha256")
				),
				c(
					"code/00_sampler.R",
					"code/sampler.cpp",
					"code/10_models.R",
					"data/model_chains.csv",
					"data/model_initial_states.rds"
				)
			)
		),
		file.path(output, "verification.json"),
		pretty = TRUE,
		auto_unbox = TRUE
	)
	stopifnot(passed)
	invisible(TRUE)
}

mode = if (length(args)) args[1] else "verify"
specification = if (length(args) >= 2L) args[2] else "current_baseline"
workers = if (length(args) >= 3L) as.integer(args[3]) else 2L
stopifnot(mode %in% c("fit", "verify"))
switch(
	mode,
	fit = replay_trade_space(specification, workers),
	verify = verify_trade_space(specification)
)
