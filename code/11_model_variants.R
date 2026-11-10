# fit alternative specifications and simulated controls

helper_file = if (sys.nframe() > 0L && !is.null(sys.frame(1)$ofile)) {
	sys.frame(1)$ofile
} else {
	sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
}

source(
	file.path(dirname(normalizePath(helper_file)), "00_data_helpers.R"),
	local = TRUE
)

suppressPackageStartupMessages(
	{
		library(data.table)
		library(jsonlite)
		library(digest)
	}
)
argument = sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
root = dirname(dirname(normalizePath(argument)))
source(file.path(root, "code", "01_observed_votes.R"))
source(file.path(root, "code", "00_sampler.R"))
source(file.path(root, "code", "10_models.R"))
variant_root = file.path(root, "results", "models", "un_variants")
chain_seeds = c(6886L, 19421L)
allies = actor_groups$allies
regular_share = .75
source_files = c(
	"code/11_model_variants.R",
	"code/00_sampler.R",
	"code/sampler.cpp",
	"code/10_models.R",
	"code/01_observed_votes.R"
)
file_hashes = function(paths) {
	setNames(
		lapply(paths, function(path) digest(file = path, algo = "sha256")),
		sub(paste0("^", root, "/"), "", paths)
	)
}

peak_rss_gb = function() {
	status = readLines("/proc/self/status")
	line = status[grepl("^VmHWM:", status)]
	as.numeric(gsub("[^0-9]", "", line)) / 1024^2
}

load_compiled = function(loader) {
	for (attempt in seq_len(6L)) {
		result = tryCatch(loader(), error = function(e) NULL)
		if (!is.null(result)) {
			return(result)
		}
		Sys.sleep(10)
	}
	stop("The compiled backend could not be loaded.")
}

variant_directory = function(variant, output) {
	if (output == "main") file.path(variant_root, variant) else file.path(variant_root, output, variant)
}

participation_shares = function(votes) {
	votes[year %in% 1989:2025, .(share_cast = mean(recorded_vote %in% c("Y", "N", "A"))), by = .(actor = voting_country_code, year)]
}

weighted_from_votes = function(votes) {
	directory = tempfile("votes_")
	dir.create(directory)
	un_votes_long = votes
	save(un_votes_long, file = file.path(directory, "un_votes_long.rda"))
	networks = build_weighted_networks(directory)
	unlink(directory, recursive = TRUE)
	simplify2array(networks$original)
}

# unweighted agreement in one year from yes, abstain and no scored 1, 0 and -1
unweighted_agreement = function(votes, year_value, actors) {
	rows = votes[year == year_value & recorded_vote %in% c("Y", "N", "A")]
	ids = sort(unique(rows$resolution))
	position = matrix(NA_real_, length(actors), length(ids), dimnames = list(actors, ids))
	position[cbind(match(rows$voting_country_code, actors), match(rows$resolution, ids))] = c(Y = 1, A = 0, N = -1)[rows$recorded_vote]
	m = 1 * is.finite(position)
	p = position
	p[!is.finite(position)] = 0
	a = 1 * (is.finite(position) & position == 0)
	support = tcrossprod(m)
	y = .5 * (support + tcrossprod(p) + tcrossprod(a)) / support
	y[support == 0] = NA_real_
	diag(y) = NA_real_
	y
}

observed_country_years = function(y) {
	sum(
		vapply(
			seq_len(dim(y)[3]),
			function(t) sum(rowSums(is.finite(y[, , t])) > 0),
			integer(1)
		)
	)
}

plant_following = function(y, weight, pair_rule = "independent") {
	stopifnot(pair_rule %in% c("independent", "joint"))
	actors = dimnames(y)[[1]]
	terminal = y[, , "2025"]
	planted = terminal
	us = terminal["USA", ]
	others = setdiff(actors, c("USA", allies))
	for (ally in allies) {
		for (j in others) {
			if (is.finite(terminal[ally, j]) && is.finite(us[j])) {
				planted[ally, j] = weight * terminal[ally, j] + (1 - weight) * us[j]
				planted[j, ally] = planted[ally, j]
			}
		}
		if (is.finite(terminal[ally, "USA"])) {
			planted[ally, "USA"] = weight * terminal[ally, "USA"] + (1 - weight)
			planted["USA", ally] = planted[ally, "USA"]
		}
	}
	for (i in seq_along(allies)) {
		for (j in seq_along(allies)) {
			if (i >= j) next
			a = allies[i]
			b = allies[j]
			if (!is.finite(terminal[a, b])) next
			if (pair_rule == "independent" && is.finite(us[a]) && is.finite(us[b])) {
				planted[a, b] = weight^2 * terminal[a, b] + weight * (1 - weight) * (us[a] + us[b]) + (1 - weight)^2
				planted[b, a] = planted[a, b]
			}
			if (pair_rule == "joint") {
				planted[a, b] = weight * terminal[a, b] + (1 - weight)
				planted[b, a] = planted[a, b]
			}
		}
	}
	out = y
	out[, , "2025"] = planted
	stopifnot(
		identical(is.finite(out), is.finite(y)),
		identical(out, aperm(out, c(2, 1, 3))),
		all(out[is.finite(out)] >= 0 & out[is.finite(out)] <= 1)
	)
	out
}

# agreement between allies in 2025 falls by a fixed amount, other ties unchanged
plant_division = function(y, drop) {
	terminal = y[, , "2025"]
	for (a in allies) {
		for (b in allies) {
			if (a != b && is.finite(terminal[a, b])) terminal[a, b] = max(terminal[a, b] - drop, 0)
		}
	}
	out = y
	out[, , "2025"] = terminal
	stopifnot(
		identical(is.finite(out), is.finite(y)),
		identical(out, aperm(out, c(2, 1, 3)))
	)
	out
}

allied_tie_means = function(y) {
	terminal = y[, , "2025"]
	others = setdiff(dimnames(y)[[1]], c("USA", allies))
	within = terminal[allies, allies]
	diag(within) = NA_real_
	to_us = terminal[allies, "USA"]
	to_others = terminal[allies, others]
	c(
		within_allied = mean(within, na.rm = TRUE),
		ally_us = mean(to_us, na.rm = TRUE),
		ally_other = mean(to_others, na.rm = TRUE),
		within_allied_ties = sum(is.finite(within)) / 2,
		ally_us_ties = sum(is.finite(to_us)),
		ally_other_ties = sum(is.finite(to_others))
	)
}

write_atomic = function(path, writer) {
	temporary = tempfile(paste0(basename(path), "_"), tmpdir = dirname(path))
	writer(temporary)
	stopifnot(file.rename(temporary, path))
	invisible(path)
}

prepare_variants = function() {
	dir.create(variant_root, recursive = TRUE, showWarnings = FALSE)
	vote_env = new.env(parent = emptyenv())
	load_vote_archive(
		file.path(root, "data", "un_votes_long.rda"),
		envir = vote_env
	)
	votes = as.data.table(vote_env$un_votes_long)
	released = read_analysis_table(
		file.path(root, "results", "models", "input_manifest.csv")
	)
	full = simplify2array(build_weighted_networks(file.path(root, "data"))$original)
	stopifnot(
		identical(
			digest(full, algo = "sha256"),
			released[domain == "un" & specification == "full", input_sha256]
		)
	)
	actors = dimnames(full)[[1]]
	years = dimnames(full)[[3]]

	matches = match_recorded_titles(
		file.path(root, "data"),
		years = 2024:2025,
		horizons = 1
	)$matches
	stopifnot(nrow(matches) == 59L)
	recurring_votes = votes[!year %in% 2024:2025 | (year == 2024L & resolution %in% matches$resolution_baseline) |
		(year == 2025L & resolution %in% matches$resolution_target)]
	recurring = weighted_from_votes(recurring_votes)
	stopifnot(
		identical(dimnames(recurring), dimnames(full)),
		identical(recurring[, , 1:35], full[, , 1:35])
	)
	# every resolution weighted one: equal yes and no totals give weight 0.1 + 0.9
	equal_votes = copy(votes)
	equal_votes[, `:=`(total_yes = 1L, total_no = 1L)]
	equal = weighted_from_votes(equal_votes)
	stopifnot(
		identical(dimnames(equal), dimnames(full)),
		identical(is.finite(equal), is.finite(full))
	)
	direct = unweighted_agreement(votes, 2025L, actors)
	stopifnot(max(abs(direct - equal[, , "2025"]), na.rm = TRUE) < 1e-12)

	common = is.finite(full[, , "2024"]) & is.finite(full[, , "2025"])
	terminal = full[, , "2024"]
	terminal[!common] = NA_real_
	us_ties = recurring["USA", , "2025"]
	terminal["USA", ] = us_ties
	terminal[, "USA"] = us_ties
	terminal[!common] = NA_real_
	us_recurring = full
	us_recurring[, , "2025"] = terminal
	released_inputs = readRDS(file.path(root, "results", "models", "inputs.rds"))$un
	other = setdiff(actors, "USA")
	stopifnot(
		identical(
			us_recurring[other, other, "2025"],
			released_inputs$us_only_change[other, other, "2025"]
		)
	)

	participation = participation_shares(votes)
	sparse = participation[share_cast < regular_share & actor %in% actors]
	sparse[, had_ties := vapply(
		seq_len(.N),
		function(k) any(is.finite(full[actor[k], , as.character(year[k])])),
		logical(1)
	)]
	regular = full
	for (k in seq_len(nrow(sparse))) {
		regular[sparse$actor[k], , as.character(sparse$year[k])] = NA_real_
		regular[, sparse$actor[k], as.character(sparse$year[k])] = NA_real_
	}
	fwrite(
		sparse[order(year, actor)],
		file.path(variant_root, "regular_voters_removed.csv")
	)
	truncated = full[, , years[as.integer(years) <= 2020L]]
	base_prior = model_priors()$full
	inputs = list(
		recurring_both_years = recurring,
		equal_weights = equal,
		us_recurring_only_change = us_recurring,
		rho_random_walk = full,
		rho_free_bound = full,
		cold_full = full,
		regular_voters_only = regular,
		truncate_2020 = truncated,
		partial_following_75 = plant_following(full, .75),
		partial_following_50 = plant_following(full, .50),
		allies_divide = plant_division(full, .20),
		partial_following_75_joint = plant_following(full, .75, "joint"),
		partial_following_50_joint = plant_following(full, .50, "joint"),
		partial_following_75_joint_controlled = plant_following(released_inputs$us_only_change, .75, "joint")
	)
	priors = lapply(inputs, function(y) base_prior)
	priors$rho_random_walk = c(base_prior, list(rho_fixed = .999, rho_bound = .9999))
	priors$rho_free_bound = c(base_prior, list(rho_bound = .9999))
	settings = data.table(
		variant = names(inputs),
		tier = c("A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B", "B", "B"),
		domain = "un",
		temporal_tau = .01,
		indicator_every = 0L,
		description = c(
			"2024 and 2025 networks built from the 59 titles recorded in both years; earlier years and weights as released",
			"every resolution weighted 1 in every year",
			"2025 is the 2024 network on dyads observed in both years, with the American 2025 ties computed from the 59 recurring titles",
			"released full network; persistence held at 0.999 with its stationary initial prior; truncation bound 0.9999",
			"released full network; persistence estimated under a prior truncated at 0.9999 instead of 0.99",
			"released full network; cold starts",
			"released full network with every country-year casting fewer than three-quarters of the year's recorded votes removed from the likelihood",
			"released full network ending in 2020",
			"released full network; each ally's 2025 ties replaced by the expected agreement when it casts 75 percent of its votes as it did and 25 percent as the United States did in 2025",
			"released full network; each ally's 2025 ties replaced by the expected agreement when it casts 50 percent of its votes as it did and 50 percent as the United States did in 2025",
			"released full network; 2025 agreement between allies reduced by 0.20",
			"released full network; each ally's 2025 ties replaced by the expected agreement when it casts 75 percent of its votes as it did and 25 percent as the United States did in 2025, every ally copying the same votes, so a tie between two allies is 0.75 times its observed value plus 0.25",
			"released full network; each ally's 2025 ties replaced by the expected agreement when it casts 50 percent of its votes as it did and 50 percent as the United States did in 2025, every ally copying the same votes, so a tie between two allies is 0.50 times its observed value plus 0.50",
			"released us_only_change network (2025 is the 2024 network on dyads observed in both years with the observed American 2025 ties); each ally's 2025 ties replaced by the expected agreement when it casts 75 percent of its votes as in 2024 and 25 percent as the United States did in 2025, every ally copying the same votes"
		)
	)
	manifest = settings[, .(variant, tier, domain, description)]
	manifest[, actors := vapply(inputs, function(y) dim(y)[1], integer(1))]
	manifest[, periods := vapply(inputs, function(y) dim(y)[3], integer(1))]
	manifest[, terminal_year := vapply(inputs, function(y) as.integer(tail(dimnames(y)[[3]], 1)), integer(1))]
	manifest[, observed_country_years := vapply(inputs, observed_country_years, integer(1))]
	manifest[, observed_dyads := vapply(inputs, function(y) sum(is.finite(y)) / 2, numeric(1))]
	manifest[, rho_bound := vapply(
		priors,
		function(p) if (is.null(p$rho_bound)) .99 else p$rho_bound,
		numeric(1)
	)]
	manifest[, rho_fixed := vapply(
		priors,
		function(p) if (is.null(p$rho_fixed)) NA_real_ else p$rho_fixed,
		numeric(1)
	)]
	manifest[, input_sha256 := vapply(inputs, function(y) digest(y, algo = "sha256"), character(1))]
	manifest[, prior_sha256 := vapply(priors, function(p) digest(p, algo = "sha256"), character(1))]
	write_atomic(
		file.path(variant_root, "inputs.rds"),
		function(path) saveRDS(list(inputs = inputs, priors = priors, settings = settings), path)
	)
	write_atomic(
		file.path(variant_root, "input_manifest.csv"),
		function(path) fwrite(manifest, path)
	)
	# mean allied ties in 2025 before and after planting, for every planted variant
	bases = list(full = full, us_only_change = released_inputs$us_only_change)
	planting = data.table(
		variant = c(
			"partial_following_75",
			"partial_following_50",
			"partial_following_75_joint",
			"partial_following_50_joint",
			"partial_following_75_joint_controlled",
			"allies_divide"
		),
		base = c("full", "full", "full", "full", "us_only_change", "full"),
		own_share = c(.75, .50, .75, .50, .75, NA_real_),
		pair_rule = c("independent", "independent", "joint", "joint", "joint", "division")
	)
	means = rbindlist(
		lapply(
			seq_len(nrow(planting)),
			function(k) {
				before = allied_tie_means(bases[[planting$base[k]]])
				after = allied_tie_means(inputs[[planting$variant[k]]])
				data.table(
					within_allied_ties = before[["within_allied_ties"]],
					ally_us_ties = before[["ally_us_ties"]],
					ally_other_ties = before[["ally_other_ties"]],
					within_allied_before = before[["within_allied"]],
					within_allied_after = after[["within_allied"]],
					ally_us_before = before[["ally_us"]],
					ally_us_after = after[["ally_us"]],
					ally_other_before = before[["ally_other"]],
					ally_other_after = after[["ally_other"]]
				)
			}
		)
	)
	fwrite(cbind(planting, means), file.path(variant_root, "planting_check.csv"))
	us_dyads = data.table(
		specification = c("us_only_change", "us_recurring_only_change"),
		us_ties_2025 = c(
			sum(is.finite(released_inputs$us_only_change["USA", , "2025"])),
			sum(is.finite(us_recurring["USA", , "2025"]))
		)
	)
	fwrite(us_dyads, file.path(variant_root, "us_ties_2025.csv"))
	write_json(
		list(
			script = "code/11_model_variants.R",
			mode = "prepare",
			inputs = c(
				file_hashes(
					file.path(
						root,
						c(
							"data/un_votes_long.rda",
							"results/models/inputs.rds",
							"results/models/input_manifest.csv"
						)
					)
				)
			),
			sources = file_hashes(file.path(root, source_files)),
			recurring_titles_2024_2025 = nrow(matches),
			regular_share = regular_share,
			sparse_country_years_removed = nrow(sparse),
			sparse_country_years_with_ties = sum(sparse$had_ties),
			chain_seeds = chain_seeds,
			R = R.version.string,
			completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")
		),
		file.path(variant_root, "input_provenance.json"),
		pretty = TRUE,
		auto_unbox = TRUE
	)
	invisible(manifest)
}

simulate_level_change = function() {
	output = file.path(variant_root, "level_change_us_size")
	dir.create(output, recursive = TRUE, showWarnings = FALSE)
	grid = expand.grid(
		seed = c(6886L, 19421L, 36583L),
		sigma = c(.05, .10),
		drop = c(.30, .45),
		stringsAsFactors = FALSE
	)
	inputs = list()
	truth = list()
	rows = list()
	for (i in seq_len(nrow(grid))) {
		base = simulate_control(grid$seed[i], grid$sigma[i], 12L, "drift")
		y = base$y
		periods = dim(y)[3]
		y[1, -1, periods] = y[1, -1, periods] - grid$drop[i]
		y[-1, 1, periods] = y[1, -1, periods]
		stopifnot(
			identical(is.finite(y), base$mask),
			identical(y, aperm(y, c(2, 1, 3))),
			identical(y[, , -periods], base$y[, , -periods]),
			sum(base$truth$reset_indicator) == 0L
		)
		hash = digest(y, algo = "sha256")
		key = paste0("level_", substr(hash, 1, 12))
		inputs[[key]] = y
		truth[[key]] = list(
			specification = base$specification,
			level_change_actor = dimnames(y)[[1]][1],
			level_drop = grid$drop[i],
			event_period = periods,
			reset_indicator = base$truth$reset_indicator,
			positions = base$truth$u
		)
		rows[[i]] = data.frame(
			case_id = i,
			grid[i, ],
			event_period = periods,
			specification = key,
			input_sha256 = hash
		)
	}
	cases = do.call(rbind, rows)
	stopifnot(!anyDuplicated(cases$specification))
	saveRDS(
		list(
			inputs = inputs,
			prior = model_priors()$controls,
			domain = "controls",
			temporal_tau = 0,
			indicator_every = 0L
		),
		file.path(output, "inputs.rds")
	)
	saveRDS(truth, file.path(output, "truth.rds"))
	fwrite(cases, file.path(output, "cases.csv"))
	write_json(
		list(
			script = "code/11_model_variants.R",
			mode = "simulate",
			case = "level_change_us_size",
			sources = file_hashes(file.path(root, source_files)),
			arrays = nrow(cases),
			generator = "simulate_control drift case with the first actor's final-period ties lowered by the drop",
			R = R.version.string,
			completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")
		),
		file.path(output, "simulate_provenance.json"),
		pretty = TRUE,
		auto_unbox = TRUE
	)
	invisible(cases)
}

load_variant = function(variant) {
	un = readRDS(file.path(variant_root, "inputs.rds"))
	if (variant %in% names(un$inputs)) {
		settings = as.list(as.data.frame(un$settings)[match(variant, un$settings$variant), ])
		return(
			list(
				domain = settings$domain,
				inputs = un$inputs[variant],
				prior = un$priors[[variant]],
				temporal_tau = settings$temporal_tau,
				indicator_every = settings$indicator_every
			)
		)
	}
	simulated = file.path(variant_root, variant, "inputs.rds")
	if (!file.exists(simulated)) stop("Unknown variant.")
	s = readRDS(simulated)
	list(
		domain = s$domain,
		inputs = s$inputs,
		prior = s$prior,
		temporal_tau = s$temporal_tau,
		indicator_every = s$indicator_every
	)
}

fit_variant = function(
	variant,
	workers,
	warmup = 3000L,
	sampling = 12000L,
	thin = 5L,
	output = "main",
	chains = 2L
) {
	spec = load_variant(variant)
	context = load_compiled(function() model_load(file.path(root, "code")))
	output_dir = variant_directory(variant, output)
	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
	jobs = expand.grid(
		specification = names(spec$inputs),
		seed = chain_seeds[seq_len(chains)],
		stringsAsFactors = FALSE
	)
	manifest = parallel::mclapply(
		seq_len(nrow(jobs)),
		function(i) {
			job = jobs[i, ]
			y = spec$inputs[[job$specification]]
			fit = fit_model(
				y,
				paste(job$specification, job$seed, sep = "_"),
				job$seed,
				spec$prior,
				context,
				output_dir,
				nscan = sampling,
				burn = warmup,
				odens = thin,
				initial_state = NULL,
				frozen_shifts = NULL,
				indicator_every = spec$indicator_every,
				temporal_tau = spec$temporal_tau
			)
			data.frame(
				variant = variant,
				domain = spec$domain,
				specification = job$specification,
				seed = job$seed,
				complete = isTRUE(fit$complete),
				warmup_sweeps = warmup,
				sampling_sweeps = sampling,
				thin = thin,
				stored_draws = if (isTRUE(fit$complete)) nrow(fit$trace) else NA_integer_,
				input_sha256 = digest(y, algo = "sha256"),
				prior_sha256 = digest(spec$prior, algo = "sha256"),
				fit_file = fit$revision_build$fit_file,
				checkpoint_file = fit$revision_build$checkpoint_file,
				elapsed_seconds = fit$revision_build$elapsed_seconds,
				peak_rss_gb = peak_rss_gb()
			)
		},
		mc.cores = workers,
		mc.preschedule = FALSE
	)
	stopifnot(all(vapply(manifest, is.data.frame, TRUE)))
	manifest = do.call(rbind, manifest)
	fwrite(manifest, file.path(output_dir, "fit_manifest.csv"))
	write_json(
		list(
			script = "code/11_model_variants.R",
			mode = "fit",
			variant = variant,
			output = output,
			domain = spec$domain,
			warmup = warmup,
			sampling = sampling,
			thin = thin,
			chains = chains,
			seeds = chain_seeds[seq_len(chains)],
			cold_start = TRUE,
			temporal_tau = spec$temporal_tau,
			indicator_every = spec$indicator_every,
			prior = spec$prior,
			sources = file_hashes(file.path(root, source_files)),
			fits = file_hashes(manifest$fit_file),
			R = R.version.string,
			completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")
		),
		file.path(output_dir, "fit_provenance.json"),
		pretty = TRUE,
		auto_unbox = TRUE
	)
	invisible(manifest)
}

break_threshold = function(sigma, pi, kappa, rank, norm) {
	bracket = norm^2 / (2 * kappa^2) - qlogis(pi) - rank * log(sigma / kappa)
	if (bracket <= 0) {
		return(NA_real_)
	}
	sigma * sqrt(2 * bracket)
}

reset_probability_table = function(features, participation) {
	by_chain = rbindlist(
		lapply(
			seq_along(features),
			function(chain) {
				log_odds = features[[chain]]$actor_log_odds
				data.table(
					chain = chain,
					id = colnames(log_odds),
					probability = colMeans(plogis(log_odds))
				)
			}
		)
	)
	by_chain[, c("actor", "year") := tstrsplit(id, "_")]
	by_chain[, year := as.integer(year)]
	by_chain = merge(by_chain, participation, by = c("actor", "year"), all.x = TRUE)
	probs = by_chain[, .(
		probability = mean(probability),
		chain_min = min(probability),
		chain_max = max(probability)
	), by = .(actor, year, share_cast)]
	probs[, group := fifelse(
		actor == "USA",
		"United States",
		fifelse(actor %in% allies, "Fourteen allies", "Other countries")
	)]
	probs[, regular := share_cast >= regular_share]
	probs[order(year, -probability)]
}

summarize_variant = function(variant, output = "main") {
	spec = load_variant(variant)
	output_dir = variant_directory(variant, output)
	manifest = read_analysis_table(file.path(output_dir, "fit_manifest.csv"))
	stopifnot(all(manifest$complete))
	backend = load_compiled(function() summary_backend(root))
	results = lapply(
		unique(manifest$specification),
		function(name) {
			rows = manifest[manifest$specification == name, ]
			fits = lapply(rows$fit_file, readRDS)
			stopifnot(
				all(
					vapply(
						fits,
						function(f) isTRUE(f$complete) && identical(f$Y, fits[[1]]$Y),
						TRUE
					)
				)
			)
			features = lapply(fits, summarize_chain, domain = spec$domain, backend = backend)
			summarize_features(features, name)
		}
	)
	domain_path = substring(output_dir, nchar(file.path(root, "results", "models")) + 2L)
	write_model_summaries(root, domain_path, results)
	diagnostics = read.csv(file.path(output_dir, "diagnostics.csv"))
	actors = read_analysis_table(file.path(output_dir, "actors.csv"))
	vote_env = new.env(parent = emptyenv())
	load_vote_archive(
		file.path(root, "data", "un_votes_long.rda"),
		envir = vote_env
	)
	participation = participation_shares(as.data.table(vote_env$un_votes_long))
	summary_rows = list()
	for (result in results) {
		continuous = do.call(rbind, lapply(result$features, `[[`, "continuous"))
		kappa = result$features[[1]]$prior$snap_kappa
		rank = 2L
		for (parameter in c("sigma_uv", "rho_uv", "pi_snap", "ve", "va", "intercept")) {
			x = continuous[, parameter]
			d = diagnostics[diagnostics$specification == result$specification & diagnostics$family == "continuous" & diagnostics$parameter == parameter, ]
			summary_rows[[length(summary_rows) + 1L]] = data.table(
				specification = result$specification,
				quantity = parameter,
				mean = mean(x),
				sd = sd(x),
				lower = unname(quantile(x, .025)),
				upper = unname(quantile(x, .975)),
				rhat = d$rhat,
				ess_bulk = d$ess_bulk,
				ess_tail = d$ess_tail
			)
		}
		sigma = mean(continuous[, "sigma_uv"])
		pi = mean(continuous[, "pi_snap"])
		for (norm in c(0, 1)) {
			summary_rows[[length(summary_rows) + 1L]] = data.table(
				specification = result$specification,
				quantity = paste0("threshold_norm_", norm),
				mean = break_threshold(sigma, pi, kappa, rank, norm),
				sd = NA_real_,
				lower = NA_real_,
				upper = NA_real_,
				rhat = NA_real_,
				ess_bulk = NA_real_,
				ess_tail = NA_real_
			)
		}
		d = diagnostics[diagnostics$specification == result$specification & !diagnostics$constant, ]
		for (quantity in c("max_rhat", "min_ess_bulk", "min_ess_tail")) {
			value = switch(
	quantity,
				max_rhat = max(d$rhat, na.rm = TRUE),
				min_ess_bulk = min(d$ess_bulk, na.rm = TRUE),
				min_ess_tail = min(d$ess_tail, na.rm = TRUE)
			)
			summary_rows[[length(summary_rows) + 1L]] = data.table(
				specification = result$specification,
				quantity = quantity,
				mean = value,
				sd = NA_real_,
				lower = NA_real_,
				upper = NA_real_,
				rhat = NA_real_,
				ess_bulk = NA_real_,
				ess_tail = NA_real_
			)
		}
	}
	fwrite(rbindlist(summary_rows), file.path(output_dir, "variant_summary.csv"))
	if (spec$domain == "un") {
		stopifnot(length(results) == 1L)
		probs = reset_probability_table(results[[1]]$features, participation)
		fwrite(probs, file.path(output_dir, "reset_probabilities.csv"))
		terminal_year = max(probs$year)
		terminal = probs[year == terminal_year & actor %in% c("USA", allies)]
		sampled = actors[actors$year == terminal_year, .(sampled_probability = mean(probability)), by = actor]
		terminal = merge(terminal, sampled, by = "actor", all.x = TRUE)
		terminal[, specification := results[[1]]$specification]
		fwrite(
			terminal[order(-probability)],
			file.path(output_dir, "terminal_probabilities.csv")
		)
	} else {
		cases = read_analysis_table(file.path(variant_root, variant, "cases.csv"))
		detection = rbindlist(
			lapply(
				results,
				function(result) {
					terminal = as.data.table(result$actors)[, .(
						probability = mean(probability),
						conditional_probability = mean(conditional_probability)
					), by = .(actor, year)]
					first = dimnames(spec$inputs[[result$specification]])[[1]][1]
					data.table(
						specification = result$specification,
						year = unique(terminal$year),
						level_change_probability = terminal[actor == first, probability],
						level_change_conditional_probability = terminal[actor == first, conditional_probability],
						unplanted_actors = terminal[actor != first, .N],
						unplanted_above_half = terminal[actor != first, sum(probability > .5)],
						unplanted_max = terminal[actor != first, max(probability)],
						unplanted_mean = terminal[actor != first, mean(probability)]
					)
				}
			)
		)
		detection = merge(cases, detection, by = "specification")
		fwrite(
			detection[order(case_id)],
			file.path(output_dir, "final_period_detection.csv")
		)
	}
	write_json(
		list(
			script = "code/11_model_variants.R",
			mode = "summarize",
			variant = variant,
			output = output,
			domain = spec$domain,
			chains = length(results[[1]]$features),
			draws_per_chain = results[[1]]$draws,
			fits = file_hashes(manifest$fit_file),
			sources = file_hashes(file.path(root, source_files)),
			inputs = file_hashes(file.path(root, "data", "un_votes_long.rda")),
			R = R.version.string,
			completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")
		),
		file.path(output_dir, "provenance.json"),
		pretty = TRUE,
		auto_unbox = TRUE
	)
	invisible(results)
}

time_variant = function(variant, warmup = 200L, sampling = 400L) {
	manifest = fit_variant(variant, 1L, warmup, sampling, 5L, "timing", chains = 1L)
	sweeps = warmup + sampling
	per_sweep = manifest$elapsed_seconds / sweeps
	memory = readLines("/proc/meminfo")
	total_gb = as.numeric(gsub("[^0-9]", "", memory[grepl("^MemTotal:", memory)])) / 1024^2
	timing = data.table(
		variant = variant,
		specification = manifest$specification,
		seed = manifest$seed,
		actors = dim(load_variant(variant)$inputs[[1]])[1],
		periods = dim(load_variant(variant)$inputs[[1]])[3],
		warmup_sweeps = warmup,
		sampling_sweeps = sampling,
		sweeps = sweeps,
		elapsed_seconds = manifest$elapsed_seconds,
		elapsed_minutes = manifest$elapsed_seconds / 60,
		seconds_per_sweep = per_sweep,
		peak_rss_gb = manifest$peak_rss_gb,
		full_budget_sweeps = 15000L,
		full_budget_minutes = per_sweep * 15000 / 60,
		reduced_budget_sweeps = 10000L,
		reduced_budget_minutes = per_sweep * 10000 / 60,
		threads = parallel::detectCores(),
		memory_gb = total_gb,
		R = R.version.string,
		completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")
	)
	fwrite(timing, file.path(variant_root, "timing", "timing.csv"))
	invisible(timing)
}

if (sys.nframe() == 0L) {
	args = commandArgs(trailingOnly = TRUE)
	mode = if (length(args)) args[1] else "prepare"
	stopifnot(mode %in% c("prepare", "simulate", "fit", "summarize", "timing"))
	if (mode == "prepare") prepare_variants()
	if (mode == "simulate") {
		stopifnot(length(args) >= 2L, args[2] == "level_change_us_size")
		simulate_level_change()
	}
	if (mode == "fit") {
		stopifnot(length(args) >= 3L)
		fit_variant(
			args[2],
			as.integer(args[3]),
			warmup = if (length(args) >= 4L) as.integer(args[4]) else 3000L,
			sampling = if (length(args) >= 5L) as.integer(args[5]) else 12000L,
			thin = if (length(args) >= 6L) as.integer(args[6]) else 5L,
			output = if (length(args) >= 7L) args[7] else "main"
		)
	}
	if (mode == "summarize") {
		stopifnot(length(args) >= 2L)
		summarize_variant(
			args[2],
			output = if (length(args) >= 3L) args[3] else "main"
		)
	}
	if (mode == "timing") {
		stopifnot(length(args) >= 2L)
		time_variant(
			args[2],
			warmup = if (length(args) >= 3L) as.integer(args[3]) else 200L,
			sampling = if (length(args) >= 4L) as.integer(args[4]) else 400L
		)
	}
}
