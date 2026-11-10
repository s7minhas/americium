# calculate observed voting comparisons

helper_file = if (sys.nframe() > 0L && !is.null(sys.frame(1)$ofile)) {
	sys.frame(1)$ofile
} else {
	sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
}

source(
	file.path(dirname(normalizePath(helper_file)), "00_data_helpers.R"),
	local = TRUE
)

suppressPackageStartupMessages(library(data.table))
Sys.setenv(TZ = "UTC", OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
data.table::setDTthreads(1L)
pipeline_paths = function() {
	args = commandArgs(trailingOnly = TRUE)
	list(
		input = normalizePath(if (length(args)) args[1] else "data"),
		output = normalizePath(
			if (length(args) >= 2L) args[2] else "results",
			mustWork = FALSE
		)
	)
}

record_execution = function(paths, name, sources, settings = list()) {
	dir.create(paths$output, recursive = TRUE, showWarnings = FALSE)
	argument = sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
	scripts = file.path(
		dirname(normalizePath(argument)),
		c("01_observed_votes.R", "02_adoptions.R", "19_figures.R")
	)
	sources = unique(c(sources, scripts))
	jsonlite::write_json(
		list(
			completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"),
			sources = lapply(
				sources,
				function(p) list(path = p, sha256 = digest::digest(file = p, algo = "sha256"))
			),
			settings = settings,
			R = R.version.string
		),
		file.path(paths$output, paste0(name, "_execution.json")),
		pretty = TRUE,
		auto_unbox = TRUE
	)
}

actor_groups = list(
	hegemon = "USA",
	allies = c(
		"FRA",
		"DEU",
		"GBR",
		"ITA",
		"CAN",
		"JPN",
		"AUS",
		"NLD",
		"ESP",
		"KOR",
		"BEL",
		"NOR",
		"DNK",
		"POL"
	),
	rivals = c("RUS", "CHN", "PRK", "IRN", "BLR", "SYR", "VEN", "CUB", "NIC"),
	hedgers = c("IND", "TUR", "BRA", "ZAF", "IDN", "SAU", "MEX", "EGY", "ARE", "VNM", "THA")
)
# prepare the fixed recurring-title comparison panel
normalize_title = function(title) {
	out = trimws(gsub("[[:space:]]+", " ", tolower(title)))
	out = sub(
		"[[:space:]]*:?[[:space:]]*resolution[[:space:]]*/[[:space:]]*adopted by (the )?general assembly[[:space:]]*$",
		"",
		out
	)
	trimws(out)
}

match_recorded_titles = function(input_root, years = 1989:2025, horizons = 1:2) {
	path = file.path(input_root, "un_votes_long.rda")
	env = new.env(parent = emptyenv())
	load_vote_archive(path, envir = env)
	raw = data.table::copy(env$un_votes_long[year %in% years])
	stopifnot(
		anyDuplicated(raw[, c("year", "resolution", "voting_country_code")]) == 0L
	)
	titles = unique(raw[, c("year", "resolution", "date", "title")])
	stopifnot(anyDuplicated(titles[, c("year", "resolution")]) == 0L)
	titles[, normalized_title := normalize_title(title)]
	titles[, within_year_resolution_count := .N, by = c("year", "normalized_title")]
	titles[, eligible_unique_title := within_year_resolution_count == 1L &
		!is.na(normalized_title) & nzchar(normalized_title)]
	votes = raw[, c(
		"year",
		"resolution",
		"voting_country_code",
		"voting_country_name",
		"recorded_vote"
	)]
	votes[, actor := voting_country_code]
	votes[year == 1989L & voting_country_code == "GER", actor := "DEU"]
	votes[, position := unname(c(Y = 1, A = 0, N = -1)[recorded_vote])]
	stopifnot(
		anyDuplicated(votes[, c("year", "resolution", "actor")]) == 0L,
		all(votes$recorded_vote %in% c("Y", "N", "A", "X"))
	)
	transitions = data.table::rbindlist(
		lapply(
			horizons,
			function(horizon) {
				baseline = years[(years + horizon) %in% years]
				data.table::data.table(
					horizon = horizon,
					baseline = baseline,
					target = baseline + horizon
				)
			}
		)
	)
	matched = support = vector("list", nrow(transitions))
	for (index in seq_len(nrow(transitions))) {
		transition = transitions[index]
		before = titles[year == transition$baseline & eligible_unique_title]
		after = titles[year == transition$target & eligible_unique_title]
		m = merge(
			before,
			after,
			by = "normalized_title",
			suffixes = c("_baseline", "_target")
		)
		data.table::setorder(m, normalized_title)
		m[, `:=`(
			horizon = transition$horizon,
			baseline = transition$baseline,
			target = transition$target
		)]
		stopifnot(
			anyDuplicated(m$resolution_baseline) == 0L,
			anyDuplicated(m$resolution_target) == 0L
		)
		matched[[index]] = m
		support[[index]] = data.table::data.table(
			horizon = transition$horizon,
			baseline = transition$baseline,
			target = transition$target,
			baseline_resolutions = nrow(titles[year == transition$baseline]),
			target_resolutions = nrow(titles[year == transition$target]),
			baseline_ambiguous = nrow(titles[year == transition$baseline & !eligible_unique_title]),
			target_ambiguous = nrow(titles[year == transition$target & !eligible_unique_title]),
			matched_titles = nrow(m)
		)
	}
	list(
		votes = votes,
		titles = titles,
		matches = data.table::rbindlist(matched),
		support = data.table::rbindlist(support),
		transitions = transitions,
		source = path
	)
}

classify_agenda = function(title, agenda_title, subjects) {
	text = paste(title, agenda_title, subjects)
	out = rep("Other/institutional", length(text))
	out[grepl(
		"DEVELOPMENT|POVERTY|FINANC|ECONOMIC|TRADE|TECHNOLOGY|ENVIRONMENT|CLIMATE|SUSTAINABLE|ISLAND|MOUNTAIN",
		text,
		ignore.case = TRUE
	)] = "Development/economic"
	out[grepl(
		"NUCLEAR|DISARMAMENT|WEAPON|ARMS|TERRORISM|SECURITY",
		text,
		ignore.case = TRUE
	)] = "Security/disarmament"
	out[grepl(
		"HUMAN RIGHTS|RIGHTS--|RACISM|MINORITIES|SELF-DETERMINATION|WOMEN",
		text,
		ignore.case = TRUE
	)] = "Human rights/self-determination"
	out[grepl(
		"PALESTIN|ISRAEL|MIDDLE EAST|UNRWA|JERUSALEM|GAZA",
		text,
		ignore.case = TRUE
	)] = "Palestine/Israel"
	out[grepl("UKRAINE|RUSSIAN|CRIMEA", text, ignore.case = TRUE)] = "Ukraine/Russia"
	out
}

build_weighted_networks = function(input_root) {
	env = new.env(parent = emptyenv())
	load_vote_archive(file.path(input_root, "un_votes_long.rda"), envir = env)
	raw = copy(env$un_votes_long)
	votes = raw[year >= 1989L & year <= 2025L & recorded_vote %in% c("Y", "N", "A")]
	votes[, pos := unname(c(Y = 1, A = 0, N = -1)[recorded_vote])]
	actors = sort(unique(votes$voting_country_code))
	years = 1989:2025
	original = without_us = setNames(vector("list", length(years)), years)
	build = function(mat, weight) {
		numerator = denominator = matrix(0, nrow(mat), nrow(mat), dimnames = list(rownames(mat), rownames(mat)))
		for (v in seq_len(ncol(mat))) {
			present = which(!is.na(mat[, v]))
			if (length(present) < 2L) {
				next
			}
			agreement = 1 - abs(outer(mat[present, v], mat[present, v], "-")) / 2
			numerator[present, present] = numerator[present, present] + weight[v] * agreement
			denominator[present, present] = denominator[present, present] + weight[v]
		}
		y = numerator / denominator
		y[denominator == 0] = NA_real_
		diag(y) = NA_real_
		y
	}
	for (year_value in years) {
		d = votes[year == year_value]
		wide = dcast(d, voting_country_code ~ resolution, value.var = "pos")
		mat = as.matrix(wide[, -1])
		rownames(mat) = wide$voting_country_code
		res = unique(d[, .(resolution, total_yes, total_no)])
		res = res[match(colnames(mat), resolution)]
		weight = 0.1 + 0.9 * fifelse(
			res$total_yes + res$total_no > 0,
			1 - abs(res$total_yes - res$total_no) / (res$total_yes + res$total_no),
			0
		)
		grid = matrix(
			NA_real_,
			length(actors),
			length(actors),
			dimnames = list(actors, actors)
		)
		grid[rownames(mat), rownames(mat)] = build(mat, weight)
		original[[as.character(year_value)]] = grid
		# the original no-us reconstruction accumulates resolutions in source order
		res = unique(as.data.frame(d)[c("resolution", "total_yes", "total_no")])
		us = d[voting_country_code == "USA", .(resolution, recorded_vote)]
		us_vote = us$recorded_vote[match(res$resolution, us$resolution)]
		yes = res$total_yes - as.integer(!is.na(us_vote) & us_vote == "Y")
		no = res$total_no - as.integer(!is.na(us_vote) & us_vote == "N")
		weight = 0.1 + 0.9 * ifelse(yes + no > 0, 1 - abs(yes - no) / (yes + no), 0)
		other = setdiff(actors, "USA")
		mat = matrix(
			NA_real_,
			length(other),
			nrow(res),
			dimnames = list(other, res$resolution)
		)
		d = d[voting_country_code %in% other]
		mat[cbind(
			match(d$voting_country_code, other),
			match(d$resolution, res$resolution)
		)] = d$pos
		without_us[[as.character(year_value)]] = build(mat, weight)
	}
	list(original = original, without_us = without_us)
}

allied_histories = function(paths) {
	output_root = file.path(paths$output, "observed/allies")
	dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
	input_root = paths$input
	replicates = 1999L
	seed = 2026091354L
	inputs = match_recorded_titles(input_root)
	fwrite(inputs$matches, file.path(output_root, "matched_title_pairs.csv"))
	fwrite(inputs$support, file.path(output_root, "matching_support.csv"))
	fwrite(inputs$titles, file.path(output_root, "title_audit.csv"))
	allies = actor_groups$allies
	europe = c("FRA", "DEU", "GBR", "ITA", "NLD", "ESP", "BEL", "NOR", "DNK", "POL")
	indopacific = c("JPN", "KOR", "AUS")
	actors = c("USA", allies)
	pairs = as.data.table(t(combn(actors, 2)))
	setnames(pairs, c("a", "b"))
	group_membership = list(
		us_allies = pairs$a == "USA",
		us_europe = pairs$a == "USA" & pairs$b %in% europe,
		us_indopacific = pairs$a == "USA" & pairs$b %in% indopacific,
		within_allies = pairs$a %in% allies & pairs$b %in% allies,
		within_europe = pairs$a %in% europe & pairs$b %in% europe,
		within_indopacific = pairs$a %in% indopacific & pairs$b %in% indopacific,
		europe_indopacific = (pairs$a %in% europe & pairs$b %in% indopacific) |
			(pairs$b %in% europe & pairs$a %in% indopacific)
	)
	for (actor in allies) {
		group_membership[[paste0("peer_", actor)]] = pairs$a != "USA" &
			(pairs$a == actor | pairs$b == actor)
	}
	group_matrix = do.call(rbind, group_membership) * 1
	group_labels = c(
		us_allies = "US–allies",
		us_europe = "US–Europe",
		us_indopacific = "US–Indo-Pacific",
		within_allies = "Within allies",
		within_europe = "Within Europe",
		within_indopacific = "Within Indo-Pacific",
		europe_indopacific = "Europe–Indo-Pacific"
	)
	vote_coverage = inputs$votes[actor %in% actors, .(
		records = .N,
		yes = sum(recorded_vote == "Y"),
		abstain = sum(recorded_vote == "A"),
		no = sum(recorded_vote == "N"),
		explicit_nonvote = sum(recorded_vote == "X"),
		original_codes = paste(sort(unique(voting_country_code)), collapse = ";")
	), by = c("year", "actor")]
	vote_coverage = merge(
		CJ(year = 1989:2025, actor = actors),
		vote_coverage,
		by = c("year", "actor"),
		all.x = TRUE
	)
	vote_coverage = merge(
		vote_coverage,
		inputs$titles[, .(resolutions = .N), by = "year"],
		by = "year"
	)
	vote_coverage[, absent_records := resolutions - fifelse(is.na(records), 0L, records)]
	fwrite(vote_coverage, file.path(output_root, "actor_vote_coverage.csv"))
	families = sort(unique(inputs$matches$normalized_title))
	set.seed(seed)
	family_weights = matrix(
		rexp(length(families) * replicates),
		length(families),
		replicates,
		dimnames = list(families, NULL)
	)
	build_codes = function(year_value, resolutions) {
		out = matrix(
			NA_character_,
			length(actors),
			length(resolutions),
			dimnames = list(actors, resolutions)
		)
		d = inputs$votes[year == year_value & resolution %in% resolutions & actor %in% actors]
		out[cbind(match(d$actor, actors), match(d$resolution, resolutions))] = d$recorded_vote
		out
	}
	positions = function(codes) {
		out = matrix(unname(c(Y = 1, A = 0, N = -1)[codes]), nrow(codes), ncol(codes))
		dimnames(out) = dimnames(codes)
		out
	}
	interval = function(draws) {
		t(apply(draws, 1, quantile, probs = c(.025, .975), names = FALSE))
	}
	agreement_bounds = function(left, right) {
		known = is.finite(left) & is.finite(right)
		lower = matrix(0, nrow(left), ncol(left))
		upper = matrix(1, nrow(left), ncol(left))
		lower[known] = upper[known] = (1 - abs(left - right) / 2)[known]
		lower[!known & ((is.finite(left) & left == 0) | (is.finite(right) & right == 0))] = .5
		list(lower = lower, upper = upper)
	}
	summary_rows = dyad_rows = bound_rows = support_rows = list()
	bootstrap = list()
	checks = list(
		data.table(
			check = "all_historical_transitions",
			passed = nrow(inputs$transitions) == 71L
		)
	)
	for (index in seq_len(nrow(inputs$transitions))) {
		transition = inputs$transitions[index]
		matched = inputs$matches[horizon == transition$horizon & target == transition$target]
		before_codes = build_codes(transition$baseline, matched$resolution_baseline)
		after_codes = build_codes(transition$target, matched$resolution_target)
		before = positions(before_codes)
		after = positions(after_codes)
		available = rowSums(!is.na(before_codes)) > 0 & rowSums(!is.na(after_codes)) > 0
		common = colSums(
			is.finite(before[available, , drop = FALSE]) &
				is.finite(after[available, , drop = FALSE])
		) == sum(available)
		weight = family_weights[matched$normalized_title, , drop = FALSE]
		original_before = 1 - abs(before[pairs$a, , drop = FALSE] - before[pairs$b, , drop = FALSE]) / 2
		original_after = 1 - abs(after[pairs$a, , drop = FALSE] - after[pairs$b, , drop = FALSE]) / 2
		original_shared = is.finite(original_before) & is.finite(original_after)
		support_rows[[index]] = data.table(
			horizon = transition$horizon,
			baseline = transition$baseline,
			target = transition$target,
			matched_titles = nrow(matched),
			common_recorded_titles = sum(common),
			available_actors = sum(available),
			available_allies = sum(available[-1]),
			unavailable_actors = paste(actors[!available], collapse = ";"),
			baseline_explicit_nonvotes = sum(before_codes == "X", na.rm = TRUE),
			target_explicit_nonvotes = sum(after_codes == "X", na.rm = TRUE),
			baseline_absent_records = sum(is.na(before_codes)),
			target_absent_records = sum(is.na(after_codes))
		)
		for (spec in c(
			"ordinal",
			"exact_category",
			"yes_no_only",
			"common_agenda",
			"support_80",
			"support_100"
		)) {
			before_agreement = original_before
			after_agreement = original_after
			shared = original_shared
			if (spec == "exact_category") {
				before_agreement = (before[pairs$a, , drop = FALSE] == before[pairs$b, , drop = FALSE]) * 1
				after_agreement = (after[pairs$a, , drop = FALSE] == after[pairs$b, , drop = FALSE]) * 1
			}
			if (spec == "yes_no_only") {
				shared = shared & abs(before[pairs$a, , drop = FALSE]) == 1 &
					abs(before[pairs$b, , drop = FALSE]) == 1 & abs(after[pairs$a, , drop = FALSE]) == 1 &
					abs(after[pairs$b, , drop = FALSE]) == 1
				shared[is.na(shared)] = FALSE
			}
			if (spec == "common_agenda") {
				shared = sweep(shared, 2, common, "&")
			}
			minimum_fraction = switch(
	spec,
				support_80 = .8,
				support_100 = 1,
				0
			)
			keep = rowSums(shared) > 0 & rowMeans(shared) >= minimum_fraction
			counts = rowSums(shared)
			before_agreement[!shared] = 0
			after_agreement[!shared] = 0
			before_mean = rowSums(before_agreement) / counts
			after_mean = rowSums(after_agreement) / counts
			difference = after_agreement - before_agreement
			draws = (difference %*% weight) / ((shared * 1) %*% weight)
			draws[!keep, ] = 0
			group_weights = sweep(group_matrix, 2, keep, "*")
			group_counts = rowSums(group_weights)
			group_weights = group_weights / pmax(group_counts, 1)
			before_mean[!keep] = after_mean[!keep] = 0
			change = after_mean - before_mean
			group_draws = group_weights %*% draws
			group_draws[group_counts == 0, ] = NA_real_
			group_intervals = t(
				vapply(
					seq_len(nrow(group_draws)),
					function(row) {
						if (group_counts[row] == 0) {
							return(c(NA_real_, NA_real_))
						}
						quantile(group_draws[row, ], c(.025, .975), names = FALSE)
					},
					numeric(2)
				)
			)
			rows = data.table(
				horizon = transition$horizon,
				baseline = transition$baseline,
				target = transition$target,
				specification = spec,
				group = names(group_membership),
				matched_titles = nrow(matched),
				n_dyads = group_counts,
				expected_dyads = rowSums(group_matrix),
				complete_group = group_counts == rowSums(group_matrix),
				baseline_agreement = drop(group_weights %*% before_mean),
				target_agreement = drop(group_weights %*% after_mean),
				change = drop(group_weights %*% change),
				mean_absolute_dyad_change = drop(group_weights %*% abs(change)),
				rms_dyad_change = sqrt(drop(group_weights %*% change^2)),
				lower_agenda_95 = group_intervals[, 1],
				upper_agenda_95 = group_intervals[, 2],
				minimum_titles_per_dyad = vapply(
					group_membership,
					function(g) {
						if (!any(g & keep)) {
							return(NA_real_)
						}
						min(counts[g & keep])
					},
					numeric(1)
				),
				mean_titles_per_dyad = drop(group_weights %*% counts)
			)
			rows[n_dyads == 0, c(
				"baseline_agreement",
				"target_agreement",
				"change",
				"mean_absolute_dyad_change",
				"rms_dyad_change",
				"mean_titles_per_dyad"
			) := NA_real_]
			key = paste(spec, transition$horizon, transition$target, sep = "__")
			bootstrap[[key]] = group_draws
			summary_rows[[length(summary_rows) + 1L]] = rows
			if (spec == "ordinal") {
				dyad_rows[[index]] = data.table(
					horizon = transition$horizon,
					baseline = transition$baseline,
					target = transition$target,
					a = pairs$a,
					b = pairs$b,
					shared_titles = counts,
					included = keep,
					baseline_agreement = ifelse(keep, before_mean, NA_real_),
					target_agreement = ifelse(keep, after_mean, NA_real_),
					change = ifelse(keep, change, NA_real_)
				)
				checks[[length(checks) + 1L]] = data.table(
					check = paste("unit_weights", transition$horizon, transition$target, sep = "_"),
					passed = max(abs(rowSums(difference)[keep] / counts[keep] - change[keep])) < 1e-12
				)
			}
		}
		before_bounds = agreement_bounds(
			before[pairs$a, , drop = FALSE],
			before[pairs$b, , drop = FALSE]
		)
		after_bounds = agreement_bounds(
			after[pairs$a, , drop = FALSE],
			after[pairs$b, , drop = FALSE]
		)
		present = !is.na(before_codes[pairs$a, , drop = FALSE]) & !is.na(before_codes[pairs$b, , drop = FALSE]) &
			!is.na(after_codes[pairs$a, , drop = FALSE]) & !is.na(after_codes[pairs$b, , drop = FALSE])
		lower = after_bounds$lower - before_bounds$upper
		upper = after_bounds$upper - before_bounds$lower
		lower[!present] = upper[!present] = 0
		keep = rowSums(present) > 0
		lower_mean = rowSums(lower) / pmax(rowSums(present), 1)
		upper_mean = rowSums(upper) / pmax(rowSums(present), 1)
		group_weights = sweep(group_matrix, 2, keep, "*")
		group_counts = rowSums(group_weights)
		group_weights = group_weights / pmax(group_counts, 1)
		bound_rows[[index]] = data.table(
			horizon = transition$horizon,
			baseline = transition$baseline,
			target = transition$target,
			group = names(group_membership),
			n_dyads = group_counts,
			lower_unknown_nonvote = drop(group_weights %*% lower_mean),
			upper_unknown_nonvote = drop(group_weights %*% upper_mean)
		)
	}
	summaries = rbindlist(summary_rows)
	dyads = rbindlist(dyad_rows)
	bounds = rbindlist(bound_rows)
	fwrite(summaries, file.path(output_root, "group_history.csv"))
	fwrite(dyads, file.path(output_root, "dyad_history.csv"))
	fwrite(bounds, file.path(output_root, "explicit_nonvote_bounds.csv"))
	fwrite(
		rbindlist(support_rows),
		file.path(output_root, "matched_actor_support.csv")
	)
	# preserve the joint contrast and historical extreme comparison
	contrast_rows = list()
	for (spec in unique(summaries$specification)) {
		for (index in seq_len(nrow(inputs$transitions))) {
			transition = inputs$transitions[index]
			key = paste(spec, transition$horizon, transition$target, sep = "__")
			d = summaries[specification == spec & horizon == transition$horizon & target == transition$target]
			for (group_name in c(
				"within_allies",
				"within_europe",
				"within_indopacific",
				"europe_indopacific"
			)) {
				group = paste0(group_name, "_minus_us_allies")
				value = d[group == group_name, change] - d[group == "us_allies", change]
				draws = bootstrap[[key]][group_name, ] - bootstrap[[key]]["us_allies", ]
				bootstrap[[key]] = rbind(bootstrap[[key]], draws)
				rownames(bootstrap[[key]])[nrow(bootstrap[[key]])] = group
				contrast_rows[[length(contrast_rows) + 1L]] = data.table(
					horizon = transition$horizon,
					baseline = transition$baseline,
					target = transition$target,
					specification = spec,
					group = group,
					change = value,
					complete_group = all(d[group %in% c(group_name, "us_allies"), complete_group]),
					lower_agenda_95 = if (all(is.finite(draws))) quantile(draws, .025, names = FALSE) else NA_real_,
					upper_agenda_95 = if (all(is.finite(draws))) quantile(draws, .975, names = FALSE) else NA_real_
				)
			}
		}
	}
	contrasts = rbindlist(contrast_rows)
	fwrite(contrasts, file.path(output_root, "allied_us_contrast_history.csv"))
	rank_input = rbindlist(
		list(
			summaries[, .(
				horizon,
				baseline,
				target,
				specification,
				group,
				change,
				complete_group
			)],
			contrasts[, .(
				horizon,
				baseline,
				target,
				specification,
				group,
				change,
				complete_group
			)]
		)
	)
	rank_rows = list()
	for (spec in unique(rank_input$specification)) {
		for (lag in 1:2) {
			for (group_name in unique(rank_input$group)) {
				for (history_sample in c("all_supported", "complete_group")) {
					d = rank_input[specification == spec & horizon == lag & group == group_name & is.finite(change)]
					if (history_sample == "complete_group") {
						d = d[complete_group == TRUE]
					}
					current = d[target == 2025]
					past = d[target < 2025]
					if (nrow(current) != 1L || nrow(past) == 0L) {
						next
					}
					key = paste(spec, lag, 2025, sep = "__")
					current_draws = bootstrap[[key]][group_name, ]
					past_draws = do.call(
						rbind,
						lapply(
							past$target,
							function(year) {
								bootstrap[[paste(spec, lag, year, sep = "__")]][group_name, ]
							}
						)
					)
					is_contrast = grepl("_minus_", group_name, fixed = TRUE)
					extreme = if (is_contrast) max else min
					gap_draws = current_draws - apply(past_draws, 2, extreme)
					absolute_gap_draws = abs(current_draws) - apply(abs(past_draws), 2, max)
					rank_rows[[length(rank_rows) + 1L]] = data.table(
						specification = spec,
						horizon = lag,
						group = group_name,
						history_sample = history_sample,
						baseline = current$baseline,
						target = 2025,
						historical_transitions = nrow(past),
						first_target = min(past$target),
						change_2025 = current$change,
						decline_rank = 1L + sum(past$change < current$change - 1e-12),
						increase_rank = 1L + sum(past$change > current$change + 1e-12),
						absolute_change_rank = 1L + sum(abs(past$change) > abs(current$change) + 1e-12),
						comparison_extreme = if (is_contrast) "maximum" else "minimum",
						prior_extreme = extreme(past$change),
						prior_extreme_target = past$target[which.min(if (is_contrast) -past$change else past$change)],
						gap_to_prior_extreme = current$change - extreme(past$change),
						gap_lower_agenda_95 = quantile(gap_draws, .025, names = FALSE),
						gap_upper_agenda_95 = quantile(gap_draws, .975, names = FALSE),
						fraction_reweights_more_extreme = mean(if (is_contrast) gap_draws > 0 else gap_draws < 0),
						absolute_gap_lower_agenda_95 = quantile(absolute_gap_draws, .025, names = FALSE),
						absolute_gap_upper_agenda_95 = quantile(absolute_gap_draws, .975, names = FALSE)
					)
				}
			}
		}
	}
	ranks = rbindlist(rank_rows)
	fwrite(ranks, file.path(output_root, "historical_ranks.csv"))
	fwrite(
		summaries[target == 2025],
		file.path(output_root, "terminal_group_results.csv")
	)
	fwrite(
		summaries[target == 2025 & grepl("^peer_", group)],
		file.path(output_root, "country_peer_results.csv")
	)
	for (left in c(-1, 0, 1, NA_real_)) {
		for (right in c(-1, 0, 1, NA_real_)) {
			candidate_left = if (is.na(left)) -1:1 else left
			candidate_right = if (is.na(right)) -1:1 else right
			possible = as.vector(1 - abs(outer(candidate_left, candidate_right, "-")) / 2)
			actual = agreement_bounds(matrix(left, 1, 1), matrix(right, 1, 1))
			checks[[length(checks) + 1L]] = data.table(
				check = paste("nonvote_bounds", left, right, sep = "_"),
				passed = actual$lower == min(possible) && actual$upper == max(possible)
			)
		}
	}
	checks[[length(checks) + 1L]] = data.table(
		check = "agreement_range",
		passed = all(abs(summaries$change[is.finite(summaries$change)]) <= 1 + 1e-12)
	)
	checks[[length(checks) + 1L]] = data.table(
		check = "two_terminal_match_counts",
		passed = identical(
			inputs$support[target == 2025][order(horizon), matched_titles],
			c(59L, 66L)
		)
	)
	validation = rbindlist(checks)
	fwrite(validation, file.path(output_root, "validation.csv"))
	stopifnot(all(validation$passed))
	invisible(NULL)
}

network_histories = function(paths) {
	output_root = file.path(paths$output, "observed/network")
	dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
	input_root = paths$input
	panel = match_recorded_titles(input_root)
	votes = panel$votes
	years = 1989:2025
	actors = sort(unique(votes$actor))
	positions = setNames(
		lapply(
			years,
			function(year_value) {
				rows = votes[year == year_value]
				resolutions = unique(rows$resolution)
				out = matrix(
					NA_real_,
					length(actors),
					length(resolutions),
					dimnames = list(actors, resolutions)
				)
				out[cbind(match(rows$actor, actors), match(rows$resolution, resolutions))] = rows$position
				out
			}
		),
		years
	)
	membership = data.table(actor = actors, category = "Other")
	membership[actor %in% actor_groups$allies, category := "Allies"]
	membership[actor %in% actor_groups$rivals, category := "Rivals"]
	membership[actor %in% actor_groups$hedgers, category := "Hedgers"]
	membership[actor == "USA", category := "USA"]
	fwrite(membership, file.path(output_root, "actor_groups.csv"))
	fwrite(
		unique(
			votes[, c("year", "actor", "voting_country_code", "voting_country_name")]
		),
		file.path(output_root, "actor_lineage.csv")
	)
	# use equal resolution weights and keep missing votes out of both sums
	agreement_network = function(position, availability = is.finite(position)) {
		stopifnot(
			identical(dim(position), dim(availability)),
			all(is.finite(position[availability]))
		)
		m = 1 * availability
		p = position
		p[!availability] = 0
		a = 1 * (availability & position == 0)
		a[!availability] = 0
		denominator = tcrossprod(m)
		numerator = 0.5 * (denominator + tcrossprod(p) + tcrossprod(a))
		y = numerator / denominator
		y[denominator == 0] = NA_real_
		diag(y) = NA_real_
		dimnames(y) = list(rownames(position), rownames(position))
		dimnames(denominator) = dimnames(y)
		list(y = y, support = denominator)
	}
	annual = lapply(positions, agreement_network)
	# fit the additive edge subspace; the constant belongs to this subspace
	projection_cache = new.env(parent = emptyenv())
	edge_projection = function(delta, from, to) {
		nodes = sort(unique(c(from, to)))
		i = match(from, nodes)
		j = match(to, nodes)
		n = length(nodes)
		m = length(delta)
		stopifnot(n >= 3L, m >= 3L, all(is.finite(delta)), all(i != j))
		key = digest::digest(list(nodes, i, j))
		if (exists(key, projection_cache, inherits = FALSE)) {
			prepared = get(key, projection_cache, inherits = FALSE)
		} else {
			gram = matrix(0, n, n)
			gram[cbind(i, j)] = 1
			gram[cbind(j, i)] = 1
			degree = tabulate(c(i, j), nbins = n)
			diag(gram) = degree
			factor = tryCatch(
				chol(gram),
				error = function(e) {
					NULL
				}
			)
			if (is.null(factor)) {
				spectral = eigen(gram, symmetric = TRUE)
				keep = spectral$values > max(spectral$values) * 1e-10
				inverse = tcrossprod(
					sweep(
						spectral$vectors[, keep, drop = FALSE],
						2,
						sqrt(spectral$values[keep]),
						"/"
					)
				)
			} else {
				inverse = NULL
				keep = rep(TRUE, n)
			}
			prepared = list(factor = factor, inverse = inverse, rank = sum(keep), degree = degree)
			assign(key, prepared, projection_cache)
		}
		rhs = rowsum(c(delta, delta), c(i, j), reorder = TRUE)[as.character(seq_len(n)), 1]
		coef = if (is.null(prepared$factor)) {
			as.vector(prepared$inverse %*% rhs)
		} else {
			backsolve(prepared$factor, forwardsolve(t(prepared$factor), rhs))
		}
		mu = mean(delta)
		actor = coef[i] + coef[j] - mu
		residual = delta - mu - actor
		energy = mean(delta^2)
		parts = c(level = mu^2, actor = mean(actor^2), residual = mean(residual^2))
		checks = c(
			mean_actor = abs(mean(actor)),
			mean_residual = abs(mean(residual)),
			actor_residual_inner_product = abs(mean(actor * residual)),
			energy_identity = abs(energy - sum(parts))
		)
		stopifnot(max(checks) < 1e-10)
		list(
			mean = mu,
			actor = actor,
			residual = residual,
			energy = energy,
			parts = parts,
			nodes = nodes,
			effects = coef - mean(coef),
			degree = prepared$degree,
			design_rank = prepared$rank,
			checks = checks
		)
	}
	# validate the projection against an independent explicit design matrix
	validation = list()
	record_check = function(check, error, tolerance = 1e-10) {
		assign(
			"validation",
			append(
				validation,
				list(
					data.table(
						check = check,
						error = error,
						tolerance = tolerance,
						passed = is.finite(error) && error <= tolerance
					)
				)
			),
			envir = parent.env(environment())
		)
	}
	for (kind in c("complete", "irregular")) {
		pairs = t(combn(letters[1:7], 2))
		if (kind == "irregular") {
			pairs = pairs[-c(2, 5, 8, 13), , drop = FALSE]
		}
		delta = sin(seq_len(nrow(pairs))) / 7 + 0.02
		actual = edge_projection(delta, pairs[, 1], pairs[, 2])
		design = outer(pairs[, 1], letters[1:7], "==") + outer(pairs[, 2], letters[1:7], "==")
		reference = as.vector(design %*% qr.solve(design, delta))
		record_check(
			paste0(kind, "_explicit_least_squares"),
			max(abs(reference - actual$mean - actual$actor))
		)
		record_check(
			paste0(kind, "_residual"),
			max(abs(delta - reference - actual$residual))
		)
	}
	toy = matrix(
		c(1, 0, -1, NA, 1, -1, 0, 1, NA, 0, 0, -1),
		nrow = 3,
		dimnames = list(letters[1:3], NULL)
	)
	toy_network = agreement_network(toy)$y
	for (i in 1:2) {
		for (j in (i + 1):3) {
			record_check(
				paste0("agreement_direct_", i, "_", j),
				abs(toy_network[i, j] - mean(1 - abs(toy[i, ] - toy[j, ]) / 2, na.rm = TRUE))
			)
		}
	}
	# keep each comparison rule fixed across the entire history
	summaries = blocks = effects = decomposed_2025 = graph_membership = list()
	category_order = c("USA", "Allies", "Rivals", "Hedgers", "Other")
	analyze_change = function(
		before,
		after,
		key,
		support = NULL,
		minimum_support = 1L,
		fixed_mask = NULL
	) {
		stopifnot(identical(dimnames(before), dimnames(after)))
		valid = is.finite(before) & is.finite(after) & upper.tri(before)
		if (!is.null(support)) {
			valid = valid & support >= minimum_support
		}
		if (!is.null(fixed_mask)) {
			valid = valid & fixed_mask
		}
		scopes = if ("USA" %in% rownames(before) && any(is.finite(before["USA", ]) & is.finite(after["USA", ]))) {
			c("full", "no_us")
		} else {
			"no_us"
		}
		for (scope in scopes) {
			selected = valid
			if (scope == "no_us") {
				selected[rownames(selected) == "USA", ] = FALSE
				selected[, colnames(selected) == "USA"] = FALSE
			}
			idx = which(selected, arr.ind = TRUE)
			if (nrow(idx) < 3L) {
				stop("A declared network comparison has fewer than three dyads.")
			}
			from = rownames(before)[idx[, 1]]
			to = colnames(before)[idx[, 2]]
			delta = after[idx] - before[idx]
			p = edge_projection(delta, from, to)
			k = c(key, list(scope = scope))
			us = from == "USA" | to == "USA"
			assign(
				"summaries",
				append(
					summaries,
					list(
						as.data.table(
							c(
								k,
								list(
									n_actors = length(p$nodes),
									n_dyads = length(delta),
									design_rank = p$design_rank,
									baseline_agreement = mean(before[idx]),
									target_agreement = mean(after[idx]),
									mean_change = p$mean,
									total_rms = sqrt(p$energy),
									level_rms = abs(p$mean),
									actor_rms = sqrt(p$parts["actor"]),
									residual_rms = sqrt(p$parts["residual"]),
									mean_absolute_change = mean(abs(delta)),
									fraction_positive = mean(delta > 0),
									level_energy_share = p$parts["level"] / p$energy,
									actor_energy_share = p$parts["actor"] / p$energy,
									residual_energy_share = p$parts["residual"] / p$energy,
									us_edge_fraction = mean(us),
									us_squared_change_share = sum(delta[us]^2) / sum(delta^2),
									minimum_dyad_support = if (is.null(support)) NA_real_ else min(support[idx]),
									median_dyad_support = if (is.null(support)) NA_real_ else median(support[idx]),
									max_projection_error = max(p$checks)
								)
							)
						)
					)
				),
				envir = parent.env(environment())
			)
			category_from = membership$category[match(from, membership$actor)]
			category_to = membership$category[match(to, membership$actor)]
			left = pmin(match(category_from, category_order), match(category_to, category_order))
			right = pmax(match(category_from, category_order), match(category_to, category_order))
			edges = data.table(
				from = from,
				to = to,
				group = paste(category_order[left], category_order[right], sep = "__"),
				before = before[idx],
				after = after[idx],
				change = delta,
				actor_component = p$actor,
				residual = p$residual
			)
			group_rows = edges[, .(
				n_dyads = .N,
				baseline_agreement = mean(before),
				target_agreement = mean(after),
				mean_change = mean(change),
				total_rms = sqrt(mean(change^2)),
				actor_component_mean = mean(actor_component),
				residual_mean = mean(residual),
				residual_rms = sqrt(mean(residual^2))
			), by = group]
			assign(
				"blocks",
				append(blocks, list(cbind(as.data.table(k), group_rows))),
				envir = parent.env(environment())
			)
			assign(
				"effects",
				append(
					effects,
					list(
						cbind(
							as.data.table(k),
							data.table(
								actor = p$nodes,
								centered_actor_coefficient = p$effects,
								n_partners = p$degree
							)
						)
					)
				),
				envir = parent.env(environment())
			)
			if (key$target == 2025L) {
				assign(
					"decomposed_2025",
					append(decomposed_2025, list(cbind(as.data.table(k), edges))),
					envir = parent.env(environment())
				)
			}
			assign(
				"graph_membership",
				append(
					graph_membership,
					list(cbind(as.data.table(k), data.table(actor = p$nodes)))
				),
				envir = parent.env(environment())
			)
		}
	}
	annual_methods = list(unweighted = lapply(annual, `[[`, "y"))
	weighted = build_weighted_networks(input_root)
	harmonize_network = function(y, year_value) {
		labels = rownames(y)
		if (year_value == 1989L && "GER" %in% labels) {
			if ("DEU" %in% labels) {
				stopifnot(all(!is.finite(y["DEU", ])))
				y = y[labels != "DEU", labels != "DEU", drop = FALSE]
				labels = rownames(y)
			}
			labels[labels == "GER"] = "DEU"
		}
		out = matrix(
			NA_real_,
			length(actors),
			length(actors),
			dimnames = list(actors, actors)
		)
		keep = labels %in% actors
		out[labels[keep], labels[keep]] = y[keep, keep, drop = FALSE]
		diag(out) = NA_real_
		out
	}
	annual_methods$archived_weighted = setNames(Map(harmonize_network, weighted$original, years), years)
	annual_methods$us_free_weights = setNames(Map(harmonize_network, weighted$without_us, years), years)
	fixed_panels = list(transition = NULL)
	for (start in c(1989L, 1992L)) {
		year_names = as.character(start:2025)
		fixed_panels[[paste0("balanced_", start)]] = Reduce(`&`, lapply(annual_methods$unweighted[year_names], is.finite))
	}
	fixed_pair_rows = list()
	for (panel_name in setdiff(names(fixed_panels), "transition")) {
		idx = which(
			fixed_panels[[panel_name]] & upper.tri(fixed_panels[[panel_name]]),
			arr.ind = TRUE
		)
		fixed_pair_rows[[panel_name]] = data.table(panel = panel_name, from = actors[idx[, 1]], to = actors[idx[, 2]])
	}
	fwrite(
		rbindlist(fixed_pair_rows),
		file.path(output_root, "balanced_dyad_membership.csv")
	)
	for (method in names(annual_methods)) {
		for (panel_name in names(fixed_panels)) {
			for (horizon in 1:2) {
				start_target = if (panel_name == "balanced_1992") 1992L + horizon else 1989L + horizon
				for (target in start_target:2025) {
					key = list(
						dataset = method,
						panel = panel_name,
						horizon = horizon,
						baseline = target - horizon,
						target = target,
						support_fraction = 0,
						matched_titles = NA_integer_
					)
					analyze_change(
						annual_methods[[method]][[as.character(target - horizon)]],
						annual_methods[[method]][[as.character(target)]],
						key,
						fixed_mask = fixed_panels[[panel_name]]
					)
				}
			}
		}
	}
	agenda_rows = agenda_dyads_2025 = resolution_contrasts = contrast_draws = list()
	families = sort(unique(panel$matches$normalized_title))
	set.seed(2026091354L)
	replicates = 1999L
	family_weights = matrix(
		rexp(length(families) * replicates),
		length(families),
		replicates,
		dimnames = list(families, NULL)
	)
	ally_pairs = t(combn(intersect(actors, actor_groups$allies), 2))
	outside_pairs = as.matrix(
		expand.grid(
			from = intersect(actors, actor_groups$allies),
			to = setdiff(actors, c("USA", actor_groups$allies)),
			stringsAsFactors = FALSE
		)
	)
	resolution_means = function(before, after, pairs) {
		a0 = 1 - abs(before[pairs[, 1], , drop = FALSE] - before[pairs[, 2], , drop = FALSE]) / 2
		a1 = 1 - abs(after[pairs[, 1], , drop = FALSE] - after[pairs[, 2], , drop = FALSE]) / 2
		shared = is.finite(a0) & is.finite(a1)
		a0[!shared] = NA_real_
		a1[!shared] = NA_real_
		list(
			before = colMeans(a0, na.rm = TRUE),
			after = colMeans(a1, na.rm = TRUE),
			n_dyads = colSums(shared)
		)
	}
	for (transition_index in seq_len(nrow(panel$transitions))) {
		transition = panel$transitions[transition_index]
		matching = panel$matches[horizon == transition$horizon & baseline == transition$baseline & target == transition$target]
		before = positions[[as.character(transition$baseline)]][, matching$resolution_baseline, drop = FALSE]
		after = positions[[as.character(transition$target)]][, matching$resolution_target, drop = FALSE]
		within = resolution_means(before, after, ally_pairs)
		between = resolution_means(before, after, outside_pairs)
		common_titles = within$n_dyads > 0L & between$n_dyads > 0L
		before_gap = within$before[common_titles] - between$before[common_titles]
		after_gap = within$after[common_titles] - between$after[common_titles]
		change_gap = after_gap - before_gap
		weight = family_weights[matching$normalized_title[common_titles], , drop = FALSE]
		resampled = as.vector(change_gap %*% weight) / colSums(weight)
		draw_key = paste(transition$horizon, transition$target, sep = "_")
		contrast_draws[[draw_key]] = resampled
		resolution_contrasts[[length(resolution_contrasts) + 1L]] = cbind(
			transition[, c("horizon", "baseline", "target")],
			data.table(
				weighting = "equal_title_after_equal_dyads_within_title",
				n_titles = sum(common_titles),
				within_baseline = mean(within$before[common_titles]),
				within_target = mean(within$after[common_titles]),
				between_baseline = mean(between$before[common_titles]),
				between_target = mean(between$after[common_titles]),
				contrast_baseline = mean(before_gap),
				contrast_target = mean(after_gap),
				contrast_change = mean(change_gap),
				lower_title_reweight_95 = unname(quantile(resampled, 0.025)),
				upper_title_reweight_95 = unname(quantile(resampled, 0.975)),
				minimum_within_dyads_per_title = min(within$n_dyads[common_titles]),
				minimum_between_dyads_per_title = min(between$n_dyads[common_titles]),
				mean_within_dyads_per_title = mean(within$n_dyads[common_titles]),
				mean_between_dyads_per_title = mean(between$n_dyads[common_titles])
			)
		)
		shared = is.finite(before) & is.finite(after)
		before_network = agreement_network(before, shared)
		after_network = agreement_network(after, shared)
		stopifnot(identical(before_network$support, after_network$support))
		literal_before = agreement_network(before)
		literal_after = agreement_network(after)
		annual_before = annual[[as.character(transition$baseline)]]
		annual_after = annual[[as.character(transition$target)]]
		idx = which(
			is.finite(before_network$y) & is.finite(after_network$y) & upper.tri(before_network$y),
			arr.ind = TRUE
		)
		account = data.table(
			from = actors[idx[, 1]],
			to = actors[idx[, 2]],
			annual_before = annual_before$y[idx],
			annual_after = annual_after$y[idx],
			title_before = literal_before$y[idx],
			title_after = literal_after$y[idx],
			shared_before = before_network$y[idx],
			shared_after = after_network$y[idx],
			annual_support_before = annual_before$support[idx],
			annual_support_after = annual_after$support[idx],
			title_support_before = literal_before$support[idx],
			title_support_after = literal_after$support[idx],
			shared_support = before_network$support[idx]
		)
		account[, `:=`(
			annual_change = annual_after - annual_before,
			matched_change = shared_after - shared_before,
			availability_component = (title_after - shared_after) - (title_before - shared_before),
			title_selection_component = (annual_after - title_after) - (annual_before - title_before)
		)]
		account[, decomposition_error := abs(
			annual_change - matched_change - availability_component - title_selection_component
		)]
		stopifnot(
			max(account$decomposition_error) < 1e-12,
			all(account$annual_support_before >= account$title_support_before),
			all(account$annual_support_after >= account$title_support_after),
			all(account$title_support_before >= account$shared_support),
			all(account$title_support_after >= account$shared_support)
		)
		for (scope in c("full", "no_us", "within_allies", "us_allies")) {
			selected = switch(
	scope,
				full = rep(TRUE, nrow(account)),
				no_us = account$from != "USA" & account$to != "USA",
				within_allies = account$from %in% actor_groups$allies & account$to %in% actor_groups$allies,
				us_allies = (account$from == "USA" & account$to %in% actor_groups$allies) |
					(account$to == "USA" & account$from %in% actor_groups$allies)
			)
			d = account[selected]
			agenda_rows[[length(agenda_rows) + 1L]] = cbind(
				transition[, c("horizon", "baseline", "target")],
				d[, .(
					scope = scope,
					n_dyads = .N,
					annual_change = mean(annual_change),
					matched_change = mean(matched_change),
					availability_component = mean(availability_component),
					title_selection_component = mean(title_selection_component),
					matched_ratio = mean(matched_change) / mean(annual_change),
					availability_ratio = mean(availability_component) / mean(annual_change),
					title_selection_ratio = mean(title_selection_component) / mean(annual_change),
					mean_annual_support_before = mean(annual_support_before),
					mean_annual_support_after = mean(annual_support_after),
					mean_title_support_before = mean(title_support_before),
					mean_title_support_after = mean(title_support_after),
					mean_shared_support = mean(shared_support),
					max_identity_error = max(decomposition_error)
				)]
			)
		}
		if (transition$target == 2025L) {
			agenda_dyads_2025[[length(agenda_dyads_2025) + 1L]] = cbind(transition[, c("horizon", "baseline", "target")], account)
		}
		for (fraction in c(0, 0.5, 0.8, 1)) {
			key = list(
				dataset = "matched_unweighted",
				panel = "transition",
				horizon = transition$horizon,
				baseline = transition$baseline,
				target = transition$target,
				support_fraction = fraction,
				matched_titles = nrow(matching)
			)
			analyze_change(
				before_network$y,
				after_network$y,
				key,
				support = before_network$support,
				minimum_support = max(1L, ceiling(fraction * nrow(matching)))
			)
		}
		if (transition$target == 2025L) {
			pairs = list(c("USA", "FRA"), c("FRA", "JPN"), c("CHN", "IND"))
			for (pair in pairs) {
				keep = shared[pair[1], ] & shared[pair[2], ]
				direct = mean(
					(abs(before[pair[1], keep] - before[pair[2], keep]) -
						abs(after[pair[1], keep] - after[pair[2], keep])) / 2
				)
				record_check(
					paste0(
						"matched_direct_",
						transition$baseline,
						"_",
						paste(pair, collapse = "_")
					),
					abs(
						direct - after_network$y[pair[1], pair[2]] + before_network$y[pair[1], pair[2]]
					)
				)
			}
		}
	}
	summary_table = rbindlist(summaries)
	block_table = rbindlist(blocks)
	effect_table = rbindlist(effects)
	ranking_keys = c("dataset", "panel", "horizon", "support_fraction", "scope")
	metrics = c("total_rms", "level_rms", "actor_rms", "residual_rms", "mean_change")
	rankings = melt(
		summary_table,
		id.vars = c(ranking_keys, "baseline", "target", "n_actors", "n_dyads"),
		measure.vars = metrics,
		variable.name = "metric",
		value.name = "value"
	)
	rankings[, `:=`(
		rank_highest = frank(-value, ties.method = "min"),
		rank_lowest = frank(value, ties.method = "min"),
		n_transitions = .N
	), by = c(ranking_keys, "metric")]
	rankings[, highest_historical_target := target[which.max(value)], by = c(ranking_keys, "metric")]
	block_rankings = melt(
		block_table,
		id.vars = c(ranking_keys, "baseline", "target", "group", "n_dyads"),
		measure.vars = c("mean_change", "residual_mean", "residual_rms"),
		variable.name = "metric",
		value.name = "value"
	)
	block_rankings[, `:=`(
		rank_highest = frank(-value, ties.method = "min"),
		rank_lowest = frank(value, ties.method = "min"),
		n_transitions = .N
	), by = c(ranking_keys, "group", "metric")]
	coalition_keys = c(ranking_keys, "baseline", "target")
	within = block_table[group == "Allies__Allies", c(
		coalition_keys,
		"n_dyads",
		"baseline_agreement",
		"target_agreement",
		"residual_mean"
	), with = FALSE]
	setnames(
		within,
		c("n_dyads", "baseline_agreement", "target_agreement", "residual_mean"),
		c(
			"within_dyads",
			"within_baseline",
			"within_target",
			"within_residual_change"
		)
	)
	between = block_table[group %in% c("Allies__Rivals", "Allies__Hedgers", "Allies__Other"),
		.(
			between_dyads = sum(n_dyads),
			between_baseline = weighted.mean(baseline_agreement, n_dyads),
			between_target = weighted.mean(target_agreement, n_dyads),
			between_residual_change = weighted.mean(residual_mean, n_dyads)
		),
		by = coalition_keys
	]
	coalition = merge(within, between, by = coalition_keys)
	coalition[, `:=`(
		contrast_baseline = within_baseline - between_baseline,
		contrast_target = within_target - between_target,
		contrast_change = (within_target - within_baseline) - (between_target - between_baseline),
		residual_contrast_change = within_residual_change - between_residual_change
	)]
	coalition[, `:=`(
		change_rank_highest = frank(-contrast_change, ties.method = "min"),
		change_rank_lowest = frank(contrast_change, ties.method = "min"),
		target_gap_rank_highest = frank(-contrast_target, ties.method = "min"),
		target_gap_rank_lowest = frank(contrast_target, ties.method = "min"),
		n_transitions = .N
	), by = ranking_keys]
	resolution_coalition = rbindlist(resolution_contrasts)
	resolution_coalition[, `:=`(
		change_rank_highest = frank(-contrast_change, ties.method = "min"),
		change_rank_lowest = frank(contrast_change, ties.method = "min"),
		target_gap_rank_highest = frank(-contrast_target, ties.method = "min"),
		target_gap_rank_lowest = frank(contrast_target, ties.method = "min"),
		n_transitions = .N
	), by = horizon]
	for (horizon_value in 1:2) {
		targets = resolution_coalition[horizon == horizon_value & target < 2025, target]
		previous = do.call(rbind, contrast_draws[paste(horizon_value, targets, sep = "_")])
		current = contrast_draws[[paste(horizon_value, 2025, sep = "_")]]
		resolution_coalition[
			horizon == horizon_value & target == 2025,
			fraction_reweights_new_largest_decline := mean(current < apply(previous, 2, min))
		]
	}
	checks = rbindlist(validation)
	checks = rbind(
		checks,
		data.table(
			check = "all_observed_energy_and_orthogonality_identities",
			error = max(summary_table$max_projection_error),
			tolerance = 1e-10,
			passed = max(summary_table$max_projection_error) < 1e-10
		)
	)
	stopifnot(
		all(checks$passed),
		all(is.finite(summary_table$total_rms)),
		all(
			abs(
				summary_table$level_energy_share + summary_table$actor_energy_share + summary_table$residual_energy_share - 1
			) < 1e-8
		)
	)
	outputs = list(
		network_history = summary_table,
		network_rankings = rankings,
		network_2025_ranks = rankings[target == 2025],
		group_history = block_table,
		group_2025_ranks = block_rankings[target == 2025],
		actor_effects = effect_table,
		graph_membership = rbindlist(graph_membership),
		agenda_accounting_history = rbindlist(agenda_rows),
		coalition_distinctiveness_history = coalition,
		coalition_distinctiveness_by_title = resolution_coalition,
		validation = checks
	)
	for (name in names(outputs)) {
		fwrite(outputs[[name]], file.path(output_root, paste0(name, ".csv")))
	}
	fwrite(panel$support, file.path(output_root, "matching_support.csv"))
	invisible(NULL)
}

vote_accounting = function(paths) {
	output_root = file.path(paths$output, "observed/vote_changes")
	dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
	input_root = paths$input
	set.seed(20260913)
	votes_path = file.path(input_root, "un_votes_long.rda")
	votes_env = new.env(parent = emptyenv())
	load_vote_archive(votes_path, envir = votes_env)
	votes = as.data.frame(votes_env$un_votes_long)
	votes = votes[votes$year %in% 1989:2025, ]
	votes$source_actor = votes$voting_country_code
	votes$voting_country_code[votes$voting_country_code == "GER" & votes$year == 1989] = "DEU"
	stopifnot(
		anyDuplicated(votes[c("year", "resolution", "voting_country_code")]) == 0
	)
	titles = unique(votes[c("year", "resolution", "date", "title")])
	stopifnot(anyDuplicated(titles[c("year", "resolution")]) == 0)
	titles$normalized_title = normalize_title(titles$title)
	title_counts = table(paste(titles$year, titles$normalized_title))
	titles$title_count = as.integer(title_counts[paste(titles$year, titles$normalized_title)])
	titles$eligible = titles$title_count == 1 & !is.na(titles$normalized_title) & nzchar(titles$normalized_title)
	position_map = c(Y = 1, A = 0, N = -1)
	votes$position = unname(position_map[votes$recorded_vote])
	actors = c("USA", actor_groups$allies)
	votes = votes[votes$voting_country_code %in% actors, ]
	agreement = function(a, b, scoring) {
		if (scoring == "Partial credit for abstention") 1 - abs(a - b) / 2 else 1 * (a == b)
	}
	contributions = function(u_before, u_after, a_before, a_after, scoring) {
		g_before = agreement(u_before, a_before, scoring)
		g_after = agreement(u_after, a_after, scoring)
		g_us_first = agreement(u_after, a_before, scoring)
		g_ally_first = agreement(u_before, a_after, scoring)
		list(
			before = g_before,
			after = g_after,
			us = ((g_us_first - g_before) + (g_after - g_ally_first)) / 2,
			ally = ((g_ally_first - g_before) + (g_after - g_us_first)) / 2
		)
	}
	# check every combination of recorded votes and both scoring rules
	test_grid = expand.grid(u_before = -1:1, u_after = -1:1, a_before = -1:1, a_after = -1:1)
	validation = list()
	for (scoring in c("Partial credit for abstention", "Exact vote equality")) {
		z = do.call(contributions, c(as.list(test_grid), list(scoring = scoring)))
		z_swapped = contributions(
			test_grid$a_before,
			test_grid$a_after,
			test_grid$u_before,
			test_grid$u_after,
			scoring
		)
		z_reverse = contributions(
			test_grid$u_after,
			test_grid$u_before,
			test_grid$a_after,
			test_grid$a_before,
			scoring
		)
		checks = c(
			add_up = max(abs(z$us + z$ally - z$after + z$before)),
			symmetry_us = max(abs(z$us - z_swapped$ally)),
			symmetry_ally = max(abs(z$ally - z_swapped$us)),
			reverse_us = max(abs(z$us + z_reverse$us)),
			reverse_ally = max(abs(z$ally + z_reverse$ally)),
			unchanged_us = max(abs(z$us[test_grid$u_before == test_grid$u_after])),
			unchanged_ally = max(abs(z$ally[test_grid$a_before == test_grid$a_after]))
		)
		stopifnot(all(checks < 1e-12))
		validation[[scoring]] = data.frame(scoring = scoring, check = names(checks), max_error = checks)
	}
	build_positions = function(year, ids) {
		x = matrix(NA_real_, length(actors), length(ids), dimnames = list(actors, ids))
		d = votes[votes$year == year & votes$resolution %in% ids, ]
		x[cbind(match(d$voting_country_code, actors), match(d$resolution, ids))] = d$position
		x
	}
	weighted_summary = function(x, support, weights) {
		x[!support] = 0
		colMeans((x %*% weights) / ((1 * support) %*% weights), na.rm = TRUE)
	}
	rows = list()
	detail = list()
	response = list()
	us_changes = list()
	matched_rows = list()
	source_rows = list()
	replicates = 2000L
	comparisons = rbind(
		data.frame(baseline = 1989:2024, target = 1990:2025),
		data.frame(baseline = 2023L, target = 2025L)
	)
	for (comparison in seq_len(nrow(comparisons))) {
		baseline = comparisons$baseline[comparison]
		target = comparisons$target[comparison]
		before_titles = titles[titles$year == baseline & titles$eligible, ]
		after_titles = titles[titles$year == target & titles$eligible, ]
		matched = merge(
			before_titles,
			after_titles,
			by = "normalized_title",
			suffixes = c("_baseline", "_target")
		)
		matched = matched[order(matched$normalized_title), ]
		stopifnot(
			nrow(matched) > 1,
			!anyDuplicated(matched$resolution_baseline),
			!anyDuplicated(matched$resolution_target)
		)
		matched_rows[[comparison]] = matched
		before = build_positions(baseline, matched$resolution_baseline)
		after = build_positions(target, matched$resolution_target)
		u_before = matrix(
			before["USA", ],
			nrow = length(actors) - 1,
			ncol = nrow(matched),
			byrow = TRUE
		)
		u_after = matrix(
			after["USA", ],
			nrow = length(actors) - 1,
			ncol = nrow(matched),
			byrow = TRUE
		)
		a_before = before[-1, , drop = FALSE]
		a_after = after[-1, , drop = FALSE]
		shared = is.finite(u_before) & is.finite(u_after) & is.finite(a_before) & is.finite(a_after)
		keep = rowSums(shared) > 0
		u_before = u_before[keep, , drop = FALSE]
		u_after = u_after[keep, , drop = FALSE]
		a_before = a_before[keep, , drop = FALSE]
		a_after = a_after[keep, , drop = FALSE]
		shared = shared[keep, , drop = FALSE]
		weights = matrix(rexp(nrow(matched) * replicates), nrow(matched), replicates)
		for (scoring in c("Partial credit for abstention", "Exact vote equality")) {
			z = contributions(u_before, u_after, a_before, a_after, scoring)
			for (name in names(z)) z[[name]][!shared] = NA_real_
			stopifnot(max(abs(z$us + z$ally - z$after + z$before), na.rm = TRUE) < 1e-12)
			for (quantity in c("before", "after", "us", "ally", "change")) {
				x = if (quantity == "change") z$after - z$before else z[[quantity]]
				resampled = weighted_summary(x, shared, weights)
				rows[[length(rows) + 1]] = data.frame(
					baseline = baseline,
					target = target,
					scoring = scoring,
					quantity = quantity,
					estimate = mean(rowMeans(x, na.rm = TRUE)),
					lower_reweighting_95 = unname(quantile(resampled, .025)),
					upper_reweighting_95 = unname(quantile(resampled, .975)),
					matched_titles = nrow(matched),
					allies = sum(keep),
					paired_votes = sum(shared),
					minimum_pairs_per_ally = min(rowSums(shared)),
					replicates = replicates,
					stringsAsFactors = FALSE
				)
			}
			indices = which(shared, arr.ind = TRUE)
			detail[[length(detail) + 1]] = data.frame(
				baseline = baseline,
				target = target,
				scoring = scoring,
				ally = rownames(a_before)[indices[, 1]],
				normalized_title = matched$normalized_title[indices[, 2]],
				u_before = u_before[shared],
				u_after = u_after[shared],
				a_before = a_before[shared],
				a_after = a_after[shared],
				agreement_before = z$before[shared],
				agreement_after = z$after[shared],
				us_contribution = z$us[shared],
				ally_contribution = z$ally[shared],
				stringsAsFactors = FALSE
			)
		}
		us_switch = shared & u_before != u_after
		room_to_follow = us_switch & a_before != u_after
		closer = abs(a_after - u_after) < abs(a_before - u_after)
		further = abs(a_after - u_after) > abs(a_before - u_after)
		for (i in seq_len(nrow(shared))) {
			response[[length(response) + 1]] = data.frame(
				baseline = baseline,
				target = target,
				ally = rownames(shared)[i],
				paired_votes = sum(shared[i, ]),
				ally_changes = sum(shared[i, ] & a_before[i, ] != a_after[i, ], na.rm = TRUE),
				us_switch_pairs = sum(us_switch[i, ], na.rm = TRUE),
				ally_unchanged_when_us_switches = sum(us_switch[i, ] & a_before[i, ] == a_after[i, ], na.rm = TRUE),
				room_to_follow = sum(room_to_follow[i, ], na.rm = TRUE),
				closer_to_new_us_vote = sum(room_to_follow[i, ] & closer[i, ], na.rm = TRUE),
				further_from_new_us_vote = sum(room_to_follow[i, ] & further[i, ], na.rm = TRUE),
				stringsAsFactors = FALSE
			)
		}
		u_shared = is.finite(before["USA", ]) & is.finite(after["USA", ])
		us_changes[[comparison]] = data.frame(
			baseline = baseline,
			target = target,
			matched_titles = nrow(matched),
			paired_us_votes = sum(u_shared),
			us_changes = sum(before["USA", u_shared] != after["USA", u_shared]),
			us_changes_to_no = sum(before["USA", u_shared] != -1 & after["USA", u_shared] == -1),
			us_changes_from_no = sum(before["USA", u_shared] == -1 & after["USA", u_shared] != -1),
			baseline_resolution_count = sum(titles$year == baseline),
			target_resolution_count = sum(titles$year == target)
		)
	}
	summary = do.call(rbind, rows)
	detail = do.call(rbind, detail)
	response = do.call(rbind, response)
	response_totals = aggregate(
		response[setdiff(names(response), c("baseline", "target", "ally"))],
		response[c("baseline", "target")],
		sum
	)
	response_totals$unchanged_share_when_us_switches = with(response_totals, ally_unchanged_when_us_switches / us_switch_pairs)
	response_totals$closer_share_when_room = with(response_totals, closer_to_new_us_vote / room_to_follow)
	response_totals$ally_switch_share = with(response_totals, ally_changes / paired_votes)
	us_changes = do.call(rbind, us_changes)
	us_changes$us_switch_share = with(us_changes, us_changes / paired_us_votes)
	rank_rows = list()
	for (scoring in unique(summary$scoring)) {
		for (quantity in c("change", "us", "ally")) {
			d = summary[summary$scoring == scoring & summary$quantity == quantity & summary$target == summary$baseline + 1, ]
			d$rank_most_negative = rank(d$estimate, ties.method = "min")
			d$rank_largest_absolute = rank(-abs(d$estimate), ties.method = "min")
			rank_rows[[length(rank_rows) + 1]] = d
		}
	}
	tables = list(
		agreement_attribution = summary,
		individual_vote_accounting = detail,
		ally_responses = response,
		response_totals = response_totals,
		us_vote_changes = us_changes,
		historical_ranks = do.call(rbind, rank_rows),
		matched_titles = do.call(rbind, matched_rows),
		identity_checks = do.call(rbind, validation),
		actor_code_mapping = unique(votes[c("year", "source_actor", "voting_country_code")])
	)
	for (name in names(tables)) {
		write.csv(
			tables[[name]],
			file.path(output_root, paste0(name, ".csv")),
			row.names = FALSE
		)
	}
	invisible(NULL)
}

annual_agreement_levels = function(paths) {
	output = file.path(paths$output, "observed/annual")
	dir.create(output, recursive = TRUE, showWarnings = FALSE)
	environment = new.env(parent = emptyenv())
	load_vote_archive(
		file.path(paths$input, "un_votes_long.rda"),
		envir = environment
	)
	votes = copy(environment$un_votes_long[year %in% 1989:2025])
	votes[, source_actor := voting_country_code]
	votes[year == 1989L & voting_country_code == "GER", voting_country_code := "DEU"]
	stopifnot(!anyDuplicated(votes[, .(year, resolution, voting_country_code)]))
	votes[, position := unname(c(Y = 1, A = 0, N = -1)[recorded_vote])]
	american = votes[voting_country_code == "USA", .(year, resolution, us_position = position)]
	paired = merge(
		votes[voting_country_code != "USA"],
		american,
		by = c("year", "resolution")
	)
	paired = paired[is.finite(position) & is.finite(us_position)]
	paired[, agreement := 1 - abs(position - us_position) / 2]
	countries = paired[, .(agreement = mean(agreement), jointly_recorded_items = .N), by = .(year, voting_country_code)]
	groups = list(
		"Core Allies" = actor_groups$allies,
		"Strategic Rivals" = actor_groups$rivals
	)
	annual = rbindlist(
		lapply(
			names(groups),
			function(group_name) {
				countries[voting_country_code %in% groups[[group_name]], .(group = group_name, agreement = mean(agreement), partners = .N), by = year]
			}
		)
	)
	setorder(annual, group, year)
	context = rbindlist(
		lapply(
			names(groups),
			function(group_name) {
				x = annual[group == group_name]
				previous = x[year < 2025]
				data.table(
					group = group_name,
					historical_years = nrow(previous),
					historical_mean = mean(previous$agreement),
					historical_minimum = min(previous$agreement),
					historical_maximum = max(previous$agreement),
					agreement_2024 = x[year == 2024, agreement],
					agreement_2025 = x[year == 2025, agreement],
					minimum_partners = min(x$partners),
					maximum_partners = max(x$partners)
				)
			}
		)
	)
	profile = votes[voting_country_code == "USA", .(
		resolutions = .N,
		recorded_votes = sum(recorded_vote %in% c("Y", "A", "N")),
		yes = sum(recorded_vote == "Y"),
		no = sum(recorded_vote == "N"),
		abstain = sum(recorded_vote == "A"),
		nonvote = sum(recorded_vote == "X")
	), by = year]
	profile[, `:=`(
		yes_share_recorded = yes / recorded_votes,
		no_share_recorded = no / recorded_votes,
		abstain_share_recorded = abstain / recorded_votes
	)]
	fwrite(annual, file.path(output, "annual_group_agreement.csv"))
	fwrite(context, file.path(output, "group_context.csv"))
	fwrite(profile, file.path(output, "us_vote_profile.csv"))
	fwrite(
		unique(votes[, .(year, source_actor, voting_country_code)]),
		file.path(output, "actor_code_mapping.csv")
	)
	invisible(context)
}

simple_forecasts = function(paths) {
	networks = build_weighted_networks(paths$input)$original
	years = as.integer(names(networks))
	mask = Reduce("&", lapply(networks[as.character(2009:2025)], is.finite))
	diag(mask) = FALSE
	actors = rownames(mask)[rowSums(mask) > 0]
	mask = mask[actors, actors, drop = FALSE]
	edges = which(upper.tri(mask) & mask, arr.ind = TRUE)
	values = vapply(
		networks,
		function(m) {
			m[actors, actors][edges]
		},
		numeric(nrow(edges))
	)
	groups = c(
		list(all = actors, focal = unique(unlist(actor_groups, use.names = FALSE))),
		actor_groups
	)
	names(groups)[names(groups) == "hegemon"] = "us"
	groups = lapply(groups, intersect, y = actors)
	actor_average = function(x) {
		m = matrix(
			NA_real_,
			length(actors),
			length(actors),
			dimnames = list(actors, actors)
		)
		m[edges] = x
		m[cbind(edges[, 2], edges[, 1])] = x
		rowMeans(m, na.rm = TRUE)
	}
	actor_rows = summary_rows = list()
	for (window in c(1L, 3L)) {
		method = if (window == 1L) "persistence" else "previous_three_year_mean"
		for (year in 2010:2025) {
			t = match(year, years)
			previous = seq.int(t - window, t - 1L)
			stopifnot(all(years[previous] < year))
			error = values[, t] - rowMeans(values[, previous, drop = FALSE], na.rm = TRUE)
			stopifnot(all(is.finite(error)))
			by_actor = data.frame(
				actor = actors,
				mae = actor_average(abs(error)),
				rmse = sqrt(actor_average(error^2)),
				method = method,
				target_year = year
			)
			actor_rows[[paste(method, year)]] = by_actor
			for (group in names(groups)) {
				x = by_actor[by_actor$actor %in% groups[[group]], ]
				summary_rows[[paste(method, year, group)]] = data.frame(
					method = method,
					target_year = year,
					group = group,
					n_actors = nrow(x),
					n_common_dyads = nrow(edges),
					mae = mean(x$mae),
					rmse = sqrt(mean(x$rmse^2))
				)
			}
		}
	}
	summary = do.call(rbind, summary_rows)
	series = interaction(summary$method, summary$group, drop = TRUE)
	summary$worst_mae_rank = ave(-summary$mae, series, FUN = rank)
	summary$worst_rmse_rank = ave(-summary$rmse, series, FUN = rank)
	summary$comparison_years = ave(summary$target_year, series, FUN = length)
	output = file.path(paths$output, "observed", "forecasts")
	dir.create(output, recursive = TRUE, showWarnings = FALSE)
	write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
	write.csv(
		do.call(rbind, actor_rows),
		file.path(output, "actors.csv"),
		row.names = FALSE
	)
	jsonlite::write_json(
		list(
			target_years = 2010:2025,
			pair_selection_years = 2009:2025,
			actors = actors,
			pairs = nrow(edges),
			groups = groups,
			definition = "Each actor receives the mean absolute or squared error across its eligible partners, then actors receive equal group weight. RMSE is the square root after averaging squared errors. The three-year prediction averages available prior values. Future availability selects the evaluation panel; future values do not enter predictions."
		),
		file.path(output, "definition.json"),
		pretty = TRUE,
		auto_unbox = TRUE
	)
	invisible(summary)
}

if (sys.nframe() == 0L) {
	paths = pipeline_paths()
	RNGkind("Mersenne-Twister", "Inversion", "Rejection")
	allied_histories(paths)
	network_histories(paths)
	vote_accounting(paths)
	annual_agreement_levels(paths)
	simple_forecasts(paths)
	record_execution(
		paths,
		"observed",
		file.path(paths$input, "un_votes_long.rda"),
		list(
			years = 1989:2025,
			horizons = 1:2,
			title_replicates = 1999L,
			title_seed = 2026091354L,
			attribution_replicates = 2000L,
			attribution_seed = 20260913L
		)
	)
}
