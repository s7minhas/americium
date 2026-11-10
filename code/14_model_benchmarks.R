# calculate model benchmarks and simulation summaries

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
set.seed(20261031)
allies = c(
	"AUS",
	"BEL",
	"CAN",
	"DNK",
	"FRA",
	"DEU",
	"ITA",
	"JPN",
	"NLD",
	"NOR",
	"POL",
	"KOR",
	"ESP",
	"GBR"
)
focal = c("USA", allies)
minimum_transitions = 3L
final_period = 12L
final_year = "2012"
signature_cases = c("drift", "common_level", "agenda_only")
input_path = file.path(root, "results", "models", "inputs.rds")
manifest_path = file.path(root, "results", "models", "input_manifest.csv")
probability_path = file.path(root, "results", "models", "un", "reset_probabilities.csv")
case_path = file.path(root, "results", "models", "control_cases.csv")
truth_path = file.path(root, "results", "models", "control_truth.rds")
control_actor_path = file.path(root, "results", "models", "controls", "actors.csv")
control_status_path = file.path(root, "results", "models", "controls", "status.csv")
benchmark_root = file.path(root, "results", "models", "un_variants", "benchmarks")
signature_root = file.path(root, "results", "models", "un_variants", "signatures")
dir.create(benchmark_root, recursive = TRUE, showWarnings = FALSE)
dir.create(signature_root, recursive = TRUE, showWarnings = FALSE)
file_record = function(path) {
	list(
		path = sub(paste0("^", root, "/"), "", path),
		sha256 = digest(file = path, algo = "sha256")
	)
}

# the released full weighted network, checked against the input manifest
network = readRDS(input_path)$un$full
manifest = read_analysis_table(manifest_path)
stopifnot(
	identical(
		digest(network, algo = "sha256"),
		manifest[domain == "un" & specification == "full", input_sha256]
	)
)
actors = dimnames(network)[[1]]
years = as.integer(dimnames(network)[[3]])
stopifnot(all(focal %in% actors), length(dim(network)) == 3L)

transition_rows = list()
for (t in seq_along(years)[-1]) {
	previous = network[, , t - 1L]
	current = network[, , t]
	for (i in seq_along(actors)) {
		common = which(is.finite(previous[i, ]) & is.finite(current[i, ]))
		if (!length(common)) {
			next
		}
		difference = 100 * (current[i, common] - previous[i, common])
		transition_rows[[length(transition_rows) + 1L]] = data.table(
			actor = actors[i],
			year = years[t],
			partners = length(common),
			mean_agreement_previous = 100 * mean(previous[i, common]),
			mean_agreement_current = 100 * mean(current[i, common]),
			mean_change = mean(difference),
			rms_change = sqrt(mean(difference^2))
		)
	}
}

changes = rbindlist(transition_rows)
setorder(changes, actor, year)

standardize = function(x) {
	n = length(x)
	full = if (n >= minimum_transitions) (x - mean(x)) / sd(x) else rep(NA_real_, n)
	leave_one_out = vapply(
		seq_len(n),
		function(k) {
			other = x[-k]
			if (length(other) >= minimum_transitions) (x[k] - mean(other)) / sd(other) else NA_real_
		},
		numeric(1)
	)
	list(full, leave_one_out)
}

changes[, own_transitions := .N, by = actor]
changes[, c("z_mean_change", "z_mean_change_loo") := standardize(mean_change), by = actor]
changes[, c("z_rms_change", "z_rms_change_loo") := standardize(rms_change), by = actor]
changes[, own_rank_abs_mean := frank(-abs(mean_change), ties.method = "min"), by = actor]
changes[, own_rank_rms := frank(-rms_change, ties.method = "min"), by = actor]

probabilities = read_analysis_table(probability_path)
changes = merge(
	changes,
	probabilities[, .(actor, year, share_cast, regular, group, model_probability = probability)],
	by = c("actor", "year"),
	all.x = TRUE
)
model_only = probabilities[!changes, on = c("actor", "year")]

rank_within = function(value, keep) {
	out = rep(NA_integer_, length(value))
	eligible = keep & is.finite(value)
	out[eligible] = frank(-value[eligible], ties.method = "min")
	out
}

changes[, rank_abs_z_mean := rank_within(abs(z_mean_change), rep(TRUE, .N)), by = year]
changes[, rank_abs_z_mean_loo := rank_within(abs(z_mean_change_loo), rep(TRUE, .N)), by = year]
changes[, rank_abs_mean := rank_within(abs(mean_change), rep(TRUE, .N)), by = year]
changes[, rank_rms := rank_within(rms_change, rep(TRUE, .N)), by = year]
changes[, rank_z_rms := rank_within(z_rms_change, rep(TRUE, .N)), by = year]
changes[, rank_abs_z_mean_regular := rank_within(abs(z_mean_change), regular %in% TRUE), by = year]
changes[, rank_rms_regular := rank_within(rms_change, regular %in% TRUE), by = year]
changes[, countries_z := sum(is.finite(z_mean_change)), by = year]
changes[, countries_rms := sum(is.finite(rms_change)), by = year]
changes[, countries_z_regular := sum(is.finite(z_mean_change) & regular %in% TRUE), by = year]
changes[, countries_rms_regular := sum(is.finite(rms_change) & regular %in% TRUE), by = year]
# the united states and the allies in 2025
ranks_2025 = changes[year == 2025L & actor %in% focal]
ranks_2025 = ranks_2025[order(actor != "USA", rank_abs_z_mean)]
rank_columns = c(
	"actor",
	"group",
	"partners",
	"mean_agreement_previous",
	"mean_agreement_current",
	"mean_change",
	"z_mean_change",
	"z_mean_change_loo",
	"own_rank_abs_mean",
	"own_transitions",
	"rank_abs_z_mean",
	"countries_z",
	"rank_abs_z_mean_regular",
	"countries_z_regular",
	"rms_change",
	"z_rms_change",
	"own_rank_rms",
	"rank_rms",
	"countries_rms",
	"rank_rms_regular",
	"countries_rms_regular",
	"share_cast",
	"model_probability"
)
ranks_2025 = ranks_2025[, ..rank_columns]

leaders = function(value, actor, keep) {
	ordered = order(-value, na.last = NA)
	ordered = ordered[keep[ordered]]
	list(
		actor[ordered[1]],
		value[ordered[1]],
		actor[ordered[2]],
		value[ordered[2]]
	)
}

top = changes[, c(
	leaders(abs(z_mean_change), actor, rep(TRUE, .N)),
	leaders(rms_change, actor, rep(TRUE, .N)),
	leaders(abs(z_mean_change), actor, regular %in% TRUE),
	leaders(rms_change, actor, regular %in% TRUE)
), by = year]
top_names = c(
	"year",
	"top_abs_z_actor",
	"top_abs_z",
	"second_abs_z_actor",
	"second_abs_z",
	"top_rms_actor",
	"top_rms",
	"second_rms_actor",
	"second_rms",
	"top_abs_z_regular_actor",
	"top_abs_z_regular",
	"second_abs_z_regular_actor",
	"second_abs_z_regular",
	"top_rms_regular_actor",
	"top_rms_regular",
	"second_rms_regular_actor",
	"second_rms_regular"
)
setnames(top, top_names)
us_columns = c(
	"year",
	"partners",
	"mean_agreement_previous",
	"mean_agreement_current",
	"mean_change",
	"z_mean_change",
	"z_mean_change_loo",
	"rank_abs_z_mean",
	"countries_z",
	"rank_abs_z_mean_regular",
	"countries_z_regular",
	"rms_change",
	"z_rms_change",
	"rank_rms",
	"countries_rms",
	"rank_rms_regular",
	"countries_rms_regular",
	"model_probability"
)
us_years = merge(
	changes[actor == "USA", ..us_columns],
	top,
	by = "year"
)
us_years[, us_first_by_abs_z := rank_abs_z_mean == 1L]
us_years[, us_first_by_rms := rank_rms == 1L]
us_years[, us_first_by_abs_z_regular := rank_abs_z_mean_regular == 1L]
us_years[, us_first_by_rms_regular := rank_rms_regular == 1L]

year_summary = changes[, .(
	countries = .N,
	regular_countries = sum(regular %in% TRUE),
	rose = sum(mean_change > 0),
	fell = sum(mean_change < 0),
	median_mean_change = median(mean_change),
	own_largest_mean = sum(own_rank_abs_mean == 1L & own_transitions >= minimum_transitions),
	own_largest_rms = sum(own_rank_rms == 1L & own_transitions >= minimum_transitions),
	own_largest_mean_regular = sum(
		own_rank_abs_mean == 1L & own_transitions >= minimum_transitions & regular %in% TRUE
	),
	own_largest_rms_regular = sum(
		own_rank_rms == 1L & own_transitions >= minimum_transitions & regular %in% TRUE
	),
	abs_z_at_least_3 = sum(abs(z_mean_change) >= 3, na.rm = TRUE),
	us_rank_abs_z_mean = rank_abs_z_mean[actor == "USA"],
	us_rank_rms = rank_rms[actor == "USA"],
	us_rank_abs_z_mean_regular = rank_abs_z_mean_regular[actor == "USA"],
	us_rank_rms_regular = rank_rms_regular[actor == "USA"]
), by = year][order(year)]

selected_breaks = data.table(
	actor = c(
		"BGR",
		"CSK",
		"HUN",
		"POL",
		"ROU",
		"SUN",
		"UKR",
		"BLR",
		"ARG",
		"RUS",
		"BLR",
		"KAZ",
		"SRB",
		"ARG",
		"USA"
	),
	year = c(
		1990L,
		1990L,
		1990L,
		1990L,
		1990L,
		1991L,
		1991L,
		1991L,
		1991L,
		2000L,
		2000L,
		2000L,
		2008L,
		2024L,
		2025L
	)
)
selected_columns = c(
	"actor",
	"year",
	"partners",
	"mean_change",
	"z_mean_change",
	"own_rank_abs_mean",
	"own_transitions",
	"rank_abs_z_mean",
	"countries_z",
	"rank_abs_z_mean_regular",
	"rms_change",
	"own_rank_rms",
	"rank_rms",
	"countries_rms",
	"rank_rms_regular",
	"share_cast",
	"regular",
	"model_probability"
)
selected_breaks = merge(
	selected_breaks,
	changes[, ..selected_columns],
	by = c("actor", "year"),
	all.x = TRUE,
	sort = FALSE
)
setorder(changes, year, actor)
fwrite(changes, file.path(benchmark_root, "country_year_changes.csv"))
fwrite(ranks_2025, file.path(benchmark_root, "ranks_2025.csv"))
fwrite(
	us_years[order(year)],
	file.path(benchmark_root, "us_ranks_by_year.csv")
)
fwrite(year_summary, file.path(benchmark_root, "year_summary.csv"))
fwrite(
	selected_breaks,
	file.path(benchmark_root, "selected_breaks_ranks.csv")
)
network_record = c(
	file_record(input_path),
	list(element = "un$full", input_sha256 = digest(network, algo = "sha256"))
)
benchmark_provenance = list(
	script = "code/14_model_benchmarks.R",
	completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
	inputs = list(
		network = network_record,
		manifest = file_record(manifest_path),
		probabilities = file_record(probability_path)
	),
	settings = list(
		scale = "agreement in percentage points, weighted as in the model input",
		partners = "other members with observed ties to the country in both adjacent years",
		mean_change = "mean over partners of the change in agreement",
		rms_change = "square root of the mean over partners of the squared change in agreement",
		standardization = "the change less the mean of the country's own changes, divided by their standard deviation",
		leave_one_out = "the same standardization with the transition itself excluded from the mean and standard deviation",
		minimum_transitions = minimum_transitions,
		ranks = "rank 1 is the largest value in the year; absolute values are ranked for the standardized mean change",
		regular = "countries casting at least three quarters of the year's recorded votes, from reset_probabilities.csv"
	),
	coverage = list(
		country_years = nrow(changes),
		countries = uniqueN(changes$actor),
		years = range(changes$year),
		country_years_without_model_probability = sum(!is.finite(changes$model_probability)),
		model_country_years_without_benchmark = nrow(model_only)
	),
	allies = allies,
	R = R.version.string
)
write_json(
	benchmark_provenance,
	file.path(benchmark_root, "provenance.json"),
	pretty = TRUE,
	auto_unbox = TRUE
)

cases = read_analysis_table(case_path)
truth = readRDS(truth_path)
control_actors = read_analysis_table(control_actor_path)
control_status = read_analysis_table(control_status_path)
selected = cases[event_period == final_period & case %in% signature_cases][order(case_id)]
selected[, distinct_input := !duplicated(specification)]
count_rows = list()
for (row in seq_len(nrow(selected))) {
	case_row = selected[row]
	known = truth[[as.character(case_row$case_id)]]
	planted = known$truth$reset_indicator[, final_year]
	fitted = control_actors[specification == case_row$specification]
	stopifnot(
		nrow(fitted) > 0L,
		all(fitted$year == as.integer(final_year)),
		uniqueN(fitted$seed) == 2L
	)
	fitted[, planted := planted[actor]]
	pooled = fitted[, .(
		probability = mean(probability),
		conditional_probability = mean(conditional_probability),
		planted = planted[1],
		chains = .N
	), by = actor]
	by_chain = fitted[planted == 0, .(above_half = sum(probability > .5)), by = seed]
	unplanted = pooled[planted == 0]
	count_rows[[row]] = data.table(
		case_id = case_row$case_id,
		case = case_row$case,
		generation_seed = case_row$seed,
		sigma = case_row$sigma,
		specification = case_row$specification,
		distinct_input = case_row$distinct_input,
		level_change_final = isTRUE(known$truth$level_change[final_period]),
		agenda_change_final = isTRUE(known$truth$agenda_change[final_period]),
		countries = nrow(pooled),
		planted = sum(pooled$planted),
		unplanted_above_half = sum(unplanted$probability > .5),
		unplanted_above_half_chain_min = min(by_chain$above_half),
		unplanted_above_half_chain_max = max(by_chain$above_half),
		unplanted_above_half_conditional = sum(unplanted$conditional_probability > .5),
		unplanted_above_nine_tenths = sum(unplanted$probability > .9),
		unplanted_max = max(unplanted$probability),
		unplanted_mean = mean(unplanted$probability),
		first_country_probability = pooled[actor == "USA", probability],
		flagged = paste(unplanted[probability > .5][order(-probability), actor], collapse = ";"),
		all_diagnostics_pass = control_status[specification == case_row$specification, all_diagnostics_pass]
	)
}

counts = rbindlist(count_rows)
case_summary = merge(
	counts[, .(case_rows = .N), by = case],
	counts[distinct_input == TRUE, .(
		distinct_arrays = .N,
		counts_above_half = paste(unplanted_above_half, collapse = ";"),
		minimum_above_half = min(unplanted_above_half),
		median_above_half = as.numeric(median(unplanted_above_half)),
		maximum_above_half = max(unplanted_above_half),
		arrays_with_any_above_half = sum(unplanted_above_half > 0L),
		arrays_at_least_tenth_above_half = sum(unplanted_above_half >= ceiling(.1 * countries)),
		arrays_at_least_quarter_above_half = sum(unplanted_above_half >= ceiling(.25 * countries)),
		highest_unplanted_probability = max(unplanted_max),
		all_diagnostics_pass = all(all_diagnostics_pass)
	), by = case],
	by = "case"
)
case_summary = case_summary[match(signature_cases, case)]
fwrite(counts, file.path(signature_root, "final_period_counts.csv"))
fwrite(
	case_summary,
	file.path(signature_root, "final_period_case_summary.csv")
)
signature_provenance = list(
	script = "code/14_model_benchmarks.R",
	completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
	inputs = list(
		cases = file_record(case_path),
		truth = file_record(truth_path),
		actors = file_record(control_actor_path),
		status = file_record(control_status_path)
	),
	settings = list(
		cases = signature_cases,
		event_period = final_period,
		probability = "posterior mean of the sampled break indicator in the final period, averaged over the two chains of each array",
		unplanted = "simulated countries without a planted break in the final period; none is planted in these cases",
		distinct_input = "the agenda-change arrays do not depend on the drift standard deviation, so each pair of case rows shares one input"
	),
	R = R.version.string
)
write_json(
	signature_provenance,
	file.path(signature_root, "provenance.json"),
	pretty = TRUE,
	auto_unbox = TRUE
)
