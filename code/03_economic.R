# calculate trade shares and import components

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
setDTthreads(1L)
source_record = function(path) {
	list(
		path = normalizePath(path),
		bytes = file.info(path)$size,
		modified_utc = format(file.info(path)$mtime, tz = "UTC", usetz = TRUE),
		sha256 = digest(file = path, algo = "sha256")
	)
}

trade_analysis = function(trade, actors, output_dir, source_files, label, code_file) {
	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
	trade = copy(trade)
	stopifnot(
		!anyDuplicated(trade[, .(reporter, partner, year)]),
		all(is.finite(trade$exports_usd)),
		all(trade$exports_usd > 0)
	)
	fwrite(
		trade[year %in% 2024:2025],
		file.path(output_dir, "cached_terminal.csv")
	)
	fwrite(
		trade[, .(
			records = .N,
			reporters = uniqueN(reporter),
			partners = uniqueN(partner)
		), by = year],
		file.path(output_dir, "source_coverage.csv")
	)
	years = 1989:2025
	exports = setNames(
		lapply(
			years,
			function(current_year) {
				m = matrix(
					NA_real_,
					length(actors),
					length(actors),
					dimnames = list(actors, actors)
				)
				x = trade[year == current_year & reporter %in% actors & partner %in% actors]
				m[cbind(match(x$reporter, actors), match(x$partner, actors))] = x$exports_usd
				diag(m) = NA_real_
				m
			}
		),
		years
	)
	sum_directions = function(m) {
		observed = is.finite(m) | is.finite(t(m))
		m[!is.finite(m)] = 0
		out = m + t(m)
		out[!observed] = NA_real_
		diag(out) = NA_real_
		out
	}
	raw_networks = lapply(exports, sum_directions)
	model_inputs = lapply(raw_networks, log)
	saveRDS(model_inputs, file.path(output_dir, "model_inputs.rds"))
	history = list()
	actor_changes = list()
	for (index in seq.int(2, length(years))) {
		previous = exports[[index - 1L]]
		current = exports[[index]]
		common = is.finite(previous) & is.finite(current)
		p = previous
		c = current
		p[!common] = NA_real_
		c[!common] = NA_real_
		paired_previous = sum_directions(p)
		paired_current = sum_directions(c)
		delta = log(paired_current) - log(paired_previous)
		eligible = which(upper.tri(delta) & is.finite(delta), arr.ind = TRUE)
		y = delta[eligible]
		x = matrix(0, nrow(eligible), length(actors))
		x[cbind(seq_len(nrow(x)), eligible[, 1])] = 1
		x[cbind(seq_len(nrow(x)), eligible[, 2])] = 1
		x = sweep(x, 2, colMeans(x))
		fit = lm.fit(cbind(1, x), y)
		mean_change = mean(y)
		actor_fitted = fit$fitted.values - mean_change
		parts = c(
			mean = length(y) * mean_change^2,
			actors = sum(actor_fitted^2),
			pairs = sum(fit$residuals^2)
		)
		stopifnot(abs(sum(parts) - sum(y^2)) < 1e-7 * max(1, sum(y^2)))
		raw_delta = log(raw_networks[[index]]) - log(raw_networks[[index - 1L]])
		prior_directions = is.finite(previous) + t(is.finite(previous))
		next_directions = is.finite(current) + t(is.finite(current))
		history[[index - 1L]] = data.table(
			year = years[index],
			common_directions = sum(common),
			common_pairs = length(y),
			raw_pairs = sum(is.finite(raw_networks[[index]]) & upper.tri(delta)),
			lost_pairs = sum(prior_directions > 0 & next_directions == 0 & upper.tri(delta)),
			new_pairs = sum(prior_directions == 0 & next_directions > 0 & upper.tri(delta)),
			two_to_one = sum(prior_directions == 2 & next_directions == 1 & upper.tri(delta)),
			total_ratio = sum(c, na.rm = TRUE) / sum(p, na.rm = TRUE),
			mean_log_change = mean_change,
			median_log_change = median(y),
			mean_abs_log_change = mean(abs(y)),
			rms_log_change = sqrt(mean(y^2)),
			raw_mean_abs_log_change = mean(abs(raw_delta[upper.tri(delta)]), na.rm = TRUE),
			mean_squared_share = unname(parts["mean"] / sum(parts)),
			actor_squared_share = unname(parts["actors"] / sum(parts)),
			pair_squared_share = unname(parts["pairs"] / sum(parts)),
			pair_residual_rms = sqrt(mean(fit$residuals^2)),
			identical_directions = sum(previous[common] == current[common]),
			directed_previous_total = sum(p, na.rm = TRUE),
			directed_current_total = sum(c, na.rm = TRUE)
		)
		actor_changes[[index - 1L]] = rbindlist(
			lapply(
				seq_along(actors),
				function(i) {
					keep = common[i, ]
					if (sum(keep) < 10L) {
						return(NULL)
					}
					before = previous[i, keep]
					after = current[i, keep]
					log_change = log(after / before)
					data.table(
						year = years[index],
						actor = actors[i],
						partners = sum(keep),
						before_total = sum(before),
						after_total = sum(after),
						share_reallocated = sum(abs(after / sum(after) - before / sum(before))) / 2,
						log_change_sd = sd(log_change),
						identical_share = mean(before == after),
						coverage_previous = sum(before) / sum(previous[i, ], na.rm = TRUE),
						coverage_current = sum(after) / sum(current[i, ], na.rm = TRUE)
					)
				}
			)
		)
	}
	history = rbindlist(history)
	actor_changes = rbindlist(actor_changes)
	composition = actor_changes[, .(
		reporters = .N,
		mean_share_reallocated = mean(share_reallocated),
		median_share_reallocated = median(share_reallocated),
		prior_export_weighted_share_reallocated = weighted.mean(share_reallocated, before_total),
		median_coverage_previous = median(coverage_previous),
		median_coverage_current = median(coverage_current),
		usa_share_reallocated = share_reallocated[actor == "USA"][1],
		constant_growth_reporters = sum(log_change_sd < 1e-8)
	), by = year]
	history = merge(history, composition, by = "year")
	for (variable in c(
		"mean_abs_log_change",
		"rms_log_change",
		"pair_residual_rms",
		"mean_share_reallocated",
		"prior_export_weighted_share_reallocated",
		"usa_share_reallocated"
	)) {
		set(
			history,
			j = paste0(variable, "_rank"),
			value = frank(-history[[variable]], ties.method = "min")
		)
	}
	fwrite(history, file.path(output_dir, "trade_history.csv"))
	fwrite(actor_changes, file.path(output_dir, "trade_actor_history.csv"))
	# assess a fixed set of directions throughout the historical comparison
	balanced = Reduce(`&`, lapply(exports, is.finite))
	balanced_history = rbindlist(
		lapply(
			seq.int(2, length(years)),
			function(index) {
				p = exports[[index - 1L]]
				c = exports[[index]]
				p[!balanced] = NA_real_
				c[!balanced] = NA_real_
				change = log(sum_directions(c)) - log(sum_directions(p))
				change = change[upper.tri(change) & is.finite(change)]
				data.table(
					year = years[index],
					directions = sum(balanced),
					pairs = length(change),
					total_ratio = sum(c, na.rm = TRUE) / sum(p, na.rm = TRUE),
					mean_abs_log_change = mean(abs(change)),
					rms_log_change = sqrt(mean(change^2))
				)
			}
		)
	)
	balanced_history[, mean_abs_log_change_rank := frank(-mean_abs_log_change, ties.method = "min")]
	fwrite(balanced_history, file.path(output_dir, "trade_balanced_history.csv"))
	write_json(
		list(
			label = label,
			inputs = lapply(source_files, source_record),
			code = source_record(code_file),
			definition = "Exports in nominal US dollars; common observed directions within each adjacent-year comparison. Partner-share changes retain reporters with at least ten common partners. The fixed-history comparison retains directions observed in all 37 years. Components are an orthogonal descriptive projection of pairwise log changes, not causal effects.",
			created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
		),
		file.path(output_dir, "trade_provenance.json"),
		auto_unbox = TRUE,
		pretty = TRUE
	)
	invisible(
		list(
			model_inputs = model_inputs,
			history = history,
			actor_changes = actor_changes,
			balanced_history = balanced_history
		)
	)
}

trade_partner_sensitivity = function(archived, refreshed, actors, output_dir, source_files, code_file) {
	started = Sys.time()
	years = 1989:2025
	usa = match("USA", actors)
	stopifnot(!is.na(usa), !anyDuplicated(actors))
	make_exports = function(trade) {
		setNames(
			lapply(
				years,
				function(year_value) {
					m = matrix(
						NA_real_,
						length(actors),
						length(actors),
						dimnames = list(actors, actors)
					)
					x = trade[year == year_value & reporter %in% actors & partner %in% actors]
					m[cbind(match(x$reporter, actors), match(x$partner, actors))] = x$exports_usd
					diag(m) = NA_real_
					m
				}
			),
			years
		)
	}
	panels = list(
		archived = make_exports(archived),
		refreshed_terminal = make_exports(refreshed)
	)
	us = lapply(
		panels,
		function(panel) do.call(rbind, lapply(panel, function(m) m[usa, ]))
	)
	high_rank = function(values) vapply(values, function(value) 1L + sum(values > value + 1e-14), integer(1))
	share_change = function(before, after) sum(abs(after / sum(after) - before / sum(before))) / 2
	us_rows = list()
	global_rows = list()
	for (vintage in names(panels)) {
		m = us[[vintage]]
		for (panel_name in c("adjacent", "fixed_1989_2025", "fixed_1995_2025")) {
			start_year = if (panel_name == "fixed_1995_2025") 1995L else 1989L
			indices = which(years > start_year)
			fixed = colSums(is.finite(m[years >= start_year, , drop = FALSE])) == sum(years >= start_year)
			destinations = if (panel_name == "adjacent") colSums(is.finite(m)) > 0 else fixed
			for (dropped in c("none", actors[destinations])) {
				series = rbindlist(
					lapply(
						indices,
						function(index) {
							keep = if (panel_name == "adjacent") is.finite(m[index - 1L, ]) & is.finite(m[index, ]) else fixed
							if (dropped != "none") keep[match(dropped, actors)] = FALSE
							before = m[index - 1L, keep]
							after = m[index, keep]
							stopifnot(
								sum(keep) >= 2L,
								all(is.finite(before)),
								all(is.finite(after)),
								all(before > 0),
								all(after > 0)
							)
							value = share_change(before, after)
							# retain the earlier audit's cached-history comparison explicitly
							reference = value
							if (vintage == "refreshed_terminal" && years[index] < 2025) {
								old = us$archived
								reference_keep = if (panel_name == "adjacent") is.finite(old[index - 1L, ]) & is.finite(old[index, ]) else fixed
								if (dropped != "none") reference_keep[match(dropped, actors)] = FALSE
								stopifnot(
									all(is.finite(old[index - 1L, reference_keep])),
									all(is.finite(old[index, reference_keep]))
								)
								reference = share_change(old[index - 1L, reference_keep], old[index, reference_keep])
							}
							data.table(
								vintage = vintage,
								panel = panel_name,
								dropped_partner = dropped,
								year = years[index],
								partners = sum(keep),
								before_total = sum(before),
								after_total = sum(after),
								share_reallocated = value,
								coverage_previous = sum(before) / sum(m[index - 1L, ], na.rm = TRUE),
								coverage_current = sum(after) / sum(m[index, ], na.rm = TRUE),
								cached_history_reference_change = reference
							)
						}
					)
				)
				series[, rank_descending := high_rank(share_reallocated)]
				series[, rank_against_cached_history := high_rank(cached_history_reference_change)]
				series[, historical_transitions := .N]
				us_rows[[length(us_rows) + 1L]] = series
			}
		}
		exports = panels[[vintage]]
		fixed_directions = Reduce(`&`, lapply(exports, is.finite))
		fixed_reporters = Reduce(
			`&`,
			lapply(
				seq.int(2L, length(years)),
				function(index) rowSums(is.finite(exports[[index - 1L]]) & is.finite(exports[[index]])) >= 10L
			)
		)
		for (rule in c("same_reporters_adjacent_partners", "same_reporters_fixed_directions")) {
			selected = which(
				if (rule == "same_reporters_adjacent_partners") fixed_reporters else rowSums(fixed_directions) >= 10L
			)
			series = rbindlist(
				lapply(
					seq.int(2L, length(years)),
					function(index) {
						before = exports[[index - 1L]]
						after = exports[[index]]
						by_reporter = rbindlist(
							lapply(
								selected,
								function(actor) {
									keep = if (rule == "same_reporters_adjacent_partners") is.finite(before[actor, ]) & is.finite(after[actor, ]) else fixed_directions[actor, ]
									b = before[actor, keep]
									a = after[actor, keep]
									data.table(
										directions = sum(keep),
										before_total = sum(b),
										after_total = sum(a),
										share_reallocated = share_change(b, a)
									)
								}
							)
						)
						data.table(
							vintage = vintage,
							panel = rule,
							year = years[index],
							reporters = length(selected),
							directions = sum(by_reporter$directions),
							mean_tv = mean(by_reporter$share_reallocated),
							median_tv = median(by_reporter$share_reallocated),
							weighted_tv = weighted.mean(by_reporter$share_reallocated, by_reporter$before_total),
							retained_previous_value_share = sum(by_reporter$before_total) / sum(before, na.rm = TRUE),
							retained_current_value_share = sum(by_reporter$after_total) / sum(after, na.rm = TRUE)
						)
					}
				)
			)
			for (metric in c("mean_tv", "median_tv", "weighted_tv")) set(series, j = paste0(metric, "_rank"), value = high_rank(series[[metric]]))
			global_rows[[length(global_rows) + 1L]] = series
		}
	}
	us_result = rbindlist(us_rows)
	global_result = rbindlist(global_rows)
	checks = data.table(
		check = c(
			"unique_us_history_keys",
			"share_changes_in_unit_interval",
			"coverage_in_unit_interval",
			"rank_in_history_range",
			"one_deletion_per_destination",
			"fixed_panel_partner_count_constant",
			"unique_global_history_keys",
			"fixed_reporter_count_constant"
		),
		passed = c(
			!anyDuplicated(us_result[, .(vintage, panel, dropped_partner, year)]),
			all(us_result$share_reallocated >= 0 & us_result$share_reallocated <= 1),
			all(
				us_result$coverage_previous >= 0 & us_result$coverage_previous <= 1 + 1e-12 & us_result$coverage_current >= 0 & us_result$coverage_current <= 1 + 1e-12
			),
			all(
				us_result$rank_descending >= 1 & us_result$rank_descending <= us_result$historical_transitions
			),
			!any(us_result$dropped_partner == "USA"),
			all(
				us_result[panel != "adjacent", uniqueN(partners), by = .(vintage, panel, dropped_partner)]$V1 == 1L
			),
			!anyDuplicated(global_result[, .(vintage, panel, year)]),
			all(global_result[, uniqueN(reporters), by = .(vintage, panel)]$V1 == 1L)
		)
	)
	stopifnot(all(checks$passed))
	fwrite(us_result, file.path(output_dir, "us_partner_sensitivity.csv"))
	fwrite(
		global_result,
		file.path(output_dir, "reporter_composition_sensitivity.csv")
	)
	fwrite(checks, file.path(output_dir, "partner_sensitivity_checks.csv"))
	write_json(
		list(
			created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
			elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
			code = source_record(code_file),
			input_files = lapply(source_files, source_record),
			definitions = list(
				measure = "Half the sum of absolute changes in export destination shares. This is the minimum redistribution at a fixed total; it does not measure realized physical diversion or import sourcing.",
				adjacent = "Partners with positive observed US exports in both adjacent years; a destination is deleted from every year where it is observed.",
				fixed = "The same partners have positive observations in every year from the named start through 2025 within that input vintage. Every retained destination is deleted in turn, with none preserving the full panel.",
				refreshed = "1989-2023 retain the archived values; 2024 and 2025 use the refreshed annual values. The main series coherently revises both the 2023-2024 and 2024-2025 comparisons.",
				cached_history_reference = "An additional reference retains the cached earlier-transition values and substitutes only the refreshed final comparison, using the refreshed fixed partner support where applicable. Its rank is labeled separately and reproduces the earlier independent audit; it is not a uniformly refreshed history.",
				ranks = "Largest first, counting values exceeding the candidate by more than 1e-14. Rank comparisons share the same panel and destination deletion.",
				fixed_reporters = "Reporters with at least ten common positive partners in every adjacent-year comparison; partners may change by year.",
				fixed_directions = "Reporters with at least ten directions observed throughout 1989-2025; the same directions enter every comparison. Export-value coverage is relative to all observed directed values on the modeled country grid.",
				weights = "Country averages are unweighted, median, or weighted by each retained reporter previous-year exports. These select durable reporting relationships and are not world-trade aggregates."
			),
			mathematical_function_sha256 = digest(body(trade_partner_sensitivity), algo = "sha256")
		),
		file.path(output_dir, "partner_sensitivity_provenance.json"),
		pretty = TRUE,
		auto_unbox = TRUE
	)
	invisible(list(us = us_result, reporters = global_result, checks = checks))
}

trade_actor_coverage = function(trade, actors, output_dir) {
	incident = rbind(
		trade[, .(actor = reporter, partner, year)],
		trade[, .(actor = partner, partner = reporter, year)]
	)
	counts = incident[, .(
		source_records = .N,
		model_partners = uniqueN(partner[partner %in% actors & partner != actor])
	), by = .(actor, year)]
	coverage = merge(
		CJ(actor = actors, year = 1989:2025),
		counts,
		by = c("actor", "year"),
		all.x = TRUE
	)
	coverage[is.na(source_records), source_records := 0L]
	coverage[is.na(model_partners), model_partners := 0L]
	coverage[, eligible_reset := model_partners > 0 & shift(model_partners, fill = 0L) > 0, by = actor]
	coverage[, first_source_year := if (any(source_records > 0)) min(year[source_records > 0]) else NA_integer_, by = actor]
	coverage[, first_model_year := if (any(model_partners > 0)) min(year[model_partners > 0]) else NA_integer_, by = actor]
	fwrite(coverage, file.path(output_dir, "historical_actor_coverage.csv"))
	invisible(coverage)
}

trade_political_pattern = function(
	archived,
	refreshed,
	twelve_support,
	actors,
	output_dir,
	source_files,
	code_file,
	economy_scope = "political_grid"
) {
	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
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
	)
	years = 1989:2025
	stopifnot(
		length(allies) == 14L,
		all(c("USA", allies, "CHE", "CHN") %in% actors)
	)
	make_panel = function(trade) {
		setNames(
			lapply(
				years,
				function(year_value) {
					m = matrix(
						NA_real_,
						length(actors),
						length(actors),
						dimnames = list(actors, actors)
					)
					x = trade[year == year_value & reporter %in% actors & partner %in% actors]
					stopifnot(
						!anyDuplicated(x[, .(reporter, partner)]),
						all(is.finite(x$exports_usd) & x$exports_usd > 0)
					)
					m[cbind(match(x$reporter, actors), match(x$partner, actors))] = x$exports_usd
					diag(m) = NA_real_
					m
				}
			),
			years
		)
	}
	panels = list(archived = make_panel(archived), refreshed = make_panel(refreshed))
	monthly_mask = matrix(FALSE, length(actors), length(actors), dimnames = list(actors, actors))
	monthly = twelve_support[reporter %in% actors & partner %in% actors]
	monthly_mask[cbind(match(monthly$reporter, actors), match(monthly$partner, actors))] = TRUE
	diag(monthly_mask) = FALSE
	panels$refreshed_twelve_month = lapply(
		panels$refreshed,
		function(m) {
			m[!monthly_mask] = NA_real_
			m
		}
	)
	add_directions = function(m) {
		seen = is.finite(m) | is.finite(t(m))
		m[!is.finite(m)] = 0
		x = m + t(m)
		x[!seen] = NA_real_
		x
	}
	settings = list(
		paired_available_1989 = list(
			start = 1989L,
			fixed = FALSE,
			members = allies,
			excluded = character(),
			complete = FALSE
		),
		paired_core13_1989 = list(
			start = 1989L,
			fixed = FALSE,
			members = setdiff(allies, "BEL"),
			excluded = "BEL",
			complete = TRUE
		),
		paired_core14_1997 = list(
			start = 1997L,
			fixed = FALSE,
			members = allies,
			excluded = character(),
			complete = TRUE
		),
		fixed_core14_1997_2025 = list(
			start = 1997L,
			fixed = TRUE,
			members = allies,
			excluded = character(),
			complete = TRUE
		)
	)
	history = list()
	countries = list()
	us_partners = list()
	checks = list()
	check = function(name, passed) {
		assign(
			"checks",
			append(checks, list(data.table(check = name, passed = isTRUE(passed)))),
			envir = parent.env(environment())
		)
	}
	coverage = rbindlist(
		lapply(
			names(panels),
			function(vintage) {
				rbindlist(
					lapply(
						names(panels[[vintage]]),
						function(year) {
							m = panels[[vintage]][[year]]
							data.table(
								vintage = vintage,
								year = as.integer(year),
								actor = c("USA", allies),
								export_partners = rowSums(is.finite(m))[c("USA", allies)],
								mirror_import_partners = colSums(is.finite(m))[c("USA", allies)]
							)
						}
					)
				)
			}
		)
	)
	fwrite(coverage, file.path(output_dir, "political_trade_coverage.csv"))
	us_coverage = rbindlist(
		lapply(
			names(panels),
			function(vintage) {
				rbindlist(
					lapply(
						names(panels[[vintage]]),
						function(year) {
							m = panels[[vintage]][[year]]
							data.table(
								economy_scope = economy_scope,
								vintage = vintage,
								year = as.integer(year),
								partner = actors,
								us_exports_usd = m["USA", actors],
								mirror_exports_to_us_usd = m[actors, "USA"],
								exports_observed = is.finite(m["USA", actors]),
								mirror_incoming_observed = is.finite(m[actors, "USA"])
							)
						}
					)
				)
			}
		)
	)
	fwrite(us_coverage, file.path(output_dir, "political_trade_us_coverage.csv"))
	for (vintage in names(panels)) {
		panel = panels[[vintage]]
		fixed_mask = Reduce(`&`, lapply(panel[as.character(1997:2025)], is.finite))
		for (support in names(settings)) {
			setting = settings[[support]]
			members = setting$members
			for (comparison in c("annual", "two_year", "annual_timing")) {
				horizon = if (comparison == "two_year") 2L else 1L
				lookback = if (comparison == "annual") 1L else 2L
				for (target in seq.int(setting$start + lookback, 2025L)) {
					baseline = target - horizon
					before = panel[[as.character(baseline)]]
					after = panel[[as.character(target)]]
					earlier = if (comparison == "annual_timing") panel[[as.character(target - 2L)]] else NULL
					common = is.finite(before) & is.finite(after)
					if (!is.null(earlier)) common = common & is.finite(earlier)
					if (setting$fixed) common = common & fixed_mask
					for (omitted in c("none", "CHE", "CHN")) {
						keep = common
						excluded = c(setting$excluded, if (omitted != "none") omitted)
						if (length(excluded)) {
							keep[excluded, ] = FALSE
							keep[, excluded] = FALSE
						}
						b = before
						a = after
						e = earlier
						b[!keep] = NA_real_
						a[!keep] = NA_real_
						if (!is.null(e)) e[!keep] = NA_real_
						check(
							paste(vintage, support, comparison, target, omitted, "fixed_directions"),
							identical(is.finite(b), is.finite(a)) && (is.null(e) || identical(is.finite(b), is.finite(e)))
						)
						for (direction in c("exports", "mirror_imports", "symmetric")) {
							orient = function(m) {
								switch(
	direction,
									exports = m,
									mirror_imports = t(m),
									symmetric = add_directions(m)
								)
							}
							previous = orient(b)
							current = orient(a)
							first = if (!is.null(e)) orient(e) else NULL
							direction_count = switch(
	direction,
								exports = keep * 1L,
								mirror_imports = t(keep) * 1L,
								symmetric = keep + t(keep)
							)
							for (group in c("us_allies", "allied_peers_without_us")) {
								focal = if (group == "us_allies") "USA" else members
								denominator_partners = if (group == "us_allies") actors else setdiff(actors, "USA")
								by_actor = rbindlist(
									lapply(
										focal,
										function(actor) {
											partners = setdiff(denominator_partners, actor)
											peers = intersect(partners, members)
											data.table(
												actor = actor,
												denominator_before = sum(previous[actor, partners], na.rm = TRUE),
												denominator_after = sum(current[actor, partners], na.rm = TRUE),
												numerator_before = sum(previous[actor, peers], na.rm = TRUE),
												numerator_after = sum(current[actor, peers], na.rm = TRUE),
												denominator_earlier = if (!is.null(first)) sum(first[actor, partners], na.rm = TRUE) else NA_real_,
												numerator_earlier = if (!is.null(first)) sum(first[actor, peers], na.rm = TRUE) else NA_real_,
												numerator_directions = sum(direction_count[actor, peers]),
												denominator_directions = sum(direction_count[actor, partners]),
												allied_partners = sum(is.finite(previous[actor, peers])),
												all_partners = sum(is.finite(previous[actor, partners]))
											)
										}
									)
								)
								valid = by_actor$denominator_before > 0 & by_actor$denominator_after > 0 & by_actor$numerator_directions > 0 & by_actor$denominator_directions > by_actor$numerator_directions
								if (setting$complete) {
									check(
										paste(
											vintage,
											support,
											comparison,
											target,
											omitted,
											direction,
											group,
											"fixed_focal_group"
										),
										all(valid)
									)
									if (!all(valid)) {
										stop(
											paste(
												"Insufficient fixed focal coverage:",
												vintage,
												support,
												comparison,
												target,
												omitted,
												direction,
												group,
												paste(by_actor$actor[!valid], collapse = ",")
											)
										)
									}
								}
								unobserved = paste(by_actor$actor[!valid], collapse = ";")
								by_actor = by_actor[valid]
								stopifnot(nrow(by_actor) > 0L)
								by_actor[, `:=`(
									share_before = numerator_before / denominator_before,
									share_after = numerator_after / denominator_after,
									share_earlier = numerator_earlier / denominator_earlier
								)]
								by_actor[, `:=`(
									share_change_pp = 100 * (share_after - share_before),
									prior_share_change_pp = 100 * (share_before - share_earlier)
								)]
								by_actor[, acceleration_pp := share_change_pp - prior_share_change_pp]
								if (group == "us_allies") {
									partners = members[is.finite(previous["USA", members]) & is.finite(current["USA", members])]
									part = data.table(
										partner = partners,
										previous_flow_usd = previous["USA", partners],
										current_flow_usd = current["USA", partners],
										earlier_flow_usd = if (!is.null(first)) first["USA", partners] else NA_real_,
										us_total_before_usd = by_actor$denominator_before,
										us_total_after_usd = by_actor$denominator_after,
										us_total_earlier_usd = by_actor$denominator_earlier,
										reporting_directions = direction_count["USA", partners]
									)
									part[, `:=`(
										share_before = previous_flow_usd / us_total_before_usd,
										share_after = current_flow_usd / us_total_after_usd,
										share_earlier = earlier_flow_usd / us_total_earlier_usd
									)]
									part[, `:=`(
										share_change_pp = 100 * (share_after - share_before),
										prior_share_change_pp = 100 * (share_before - share_earlier),
										nominal_growth_percent = 100 * (current_flow_usd / previous_flow_usd - 1)
									)]
									part[, acceleration_pp := share_change_pp - prior_share_change_pp]
									part[, `:=`(
										vintage = vintage,
										support = support,
										omitted_partner = omitted,
										comparison = comparison,
										horizon = horizon,
										baseline = baseline,
										target = target,
										direction = direction
									)]
									check(
										paste(
											vintage,
											support,
											comparison,
											target,
											omitted,
											direction,
											"partner_contributions_sum"
										),
										abs(sum(part$share_change_pp) - by_actor$share_change_pp) < 1e-10
									)
									us_partners[[length(us_partners) + 1L]] = part
								}
								for (weighting in if (group == "us_allies") "value_weighted" else c("value_weighted", "equal_country")) {
									num_before = sum(by_actor$numerator_before)
									num_after = sum(by_actor$numerator_after)
									den_before = sum(by_actor$denominator_before)
									den_after = sum(by_actor$denominator_after)
									share_before = if (weighting == "value_weighted") num_before / den_before else mean(by_actor$share_before)
									share_after = if (weighting == "value_weighted") num_after / den_after else mean(by_actor$share_after)
									share_earlier = if (weighting == "value_weighted") sum(by_actor$numerator_earlier) / sum(by_actor$denominator_earlier) else mean(by_actor$share_earlier)
									history[[length(history) + 1L]] = data.table(
										vintage = vintage,
										support = support,
										omitted_partner = omitted,
										comparison = comparison,
										horizon = horizon,
										baseline = baseline,
										target = target,
										direction = direction,
										group = group,
										weighting = weighting,
										allied_group_size = length(members),
										focal_countries = nrow(by_actor),
										unobserved_focal = unobserved,
										share_before = share_before,
										share_after = share_after,
										share_change_pp = 100 * (share_after - share_before),
										share_earlier = share_earlier,
										prior_share_change_pp = 100 * (share_before - share_earlier),
										acceleration_pp = 100 * (share_after - 2 * share_before + share_earlier),
										allied_flow_before_usd = num_before,
										allied_flow_after_usd = num_after,
										comparison_flow_before_usd = den_before - num_before,
										comparison_flow_after_usd = den_after - num_after,
										allied_flow_growth_percent = 100 * (num_after / num_before - 1),
										comparison_flow_growth_percent = 100 * ((den_after - num_after) / (den_before - num_before) - 1),
										relative_log_growth = log(num_after / num_before) - log((den_after - num_after) / (den_before - num_before)),
										numerator_directions = sum(by_actor$numerator_directions),
										denominator_directions = sum(by_actor$denominator_directions),
										min_allied_partners = min(by_actor$allied_partners),
										min_all_partners = min(by_actor$all_partners),
										countries_share_declines = sum(by_actor$share_change_pp < 0)
									)
								}
								by_actor[, `:=`(
									vintage = vintage,
									support = support,
									omitted_partner = omitted,
									comparison = comparison,
									horizon = horizon,
									baseline = baseline,
									target = target,
									direction = direction,
									group = group
								)]
								countries[[length(countries) + 1L]] = by_actor
							}
						}
					}
				}
			}
		}
	}
	result = rbindlist(history)
	keys = c(
		"vintage",
		"support",
		"omitted_partner",
		"comparison",
		"direction",
		"group",
		"weighting"
	)
	result[, `:=`(
		rank_largest_decline = frank(share_change_pp, ties.method = "min"),
		rank_largest_increase = frank(-share_change_pp, ties.method = "min"),
		rank_largest_slowdown = if (all(is.finite(acceleration_pp))) frank(acceleration_pp, ties.method = "min") else NA_integer_,
		historical_transitions = .N
	), by = keys]
	check(
		"unique_history_keys",
		!anyDuplicated(result[, c(keys, "target"), with = FALSE])
	)
	check(
		"share_bounds",
		all(
			result$share_before >= 0 & result$share_before <= 1 & result$share_after >= 0 & result$share_after <= 1
		)
	)
	check(
		"complete_fourteen_allies_after_coverage_start",
		all(
			result[group == "allied_peers_without_us" & support %in% c("paired_core14_1997", "fixed_core14_1997_2025"), focal_countries] == 14L
		)
	)
	check(
		"complete_thirteen_allies_in_fixed_membership_history",
		all(
			result[group == "allied_peers_without_us" & support == "paired_core13_1989", focal_countries] == 13L
		)
	)
	check(
		"fixed_history_direction_counts",
		all(
			result[support == "fixed_core14_1997_2025", uniqueN(denominator_directions), by = keys]$V1 == 1L
		)
	)
	check(
		"symmetric_internal_growth_equals_directed_growth",
		max(
			result[group == "allied_peers_without_us", .(
				difference = max(allied_flow_growth_percent) - min(allied_flow_growth_percent)
			), by = .(vintage, support, omitted_partner, comparison, target)]$difference
		) < 1e-10
	)
	validation = rbindlist(checks)
	stopifnot(all(validation$passed))
	country_history = rbindlist(countries)
	country_keys = c(setdiff(keys, "weighting"), "actor")
	country_history[, `:=`(
		rank_largest_decline = frank(share_change_pp, ties.method = "min"),
		rank_largest_increase = frank(-share_change_pp, ties.method = "min"),
		rank_largest_slowdown = if (all(is.finite(acceleration_pp))) frank(acceleration_pp, ties.method = "min") else NA_integer_,
		historical_transitions = .N
	), by = country_keys]
	stopifnot(
		!anyDuplicated(country_history[, c(country_keys, "target"), with = FALSE])
	)
	fwrite(result, file.path(output_dir, "political_trade_history.csv"))
	fwrite(
		result[target == 2025],
		file.path(output_dir, "political_trade_terminal.csv")
	)
	fwrite(
		country_history[target == 2025],
		file.path(output_dir, "political_trade_countries_2025.csv")
	)
	partner_history = rbindlist(us_partners)
	partner_keys = c(setdiff(keys, c("weighting", "group")), "partner")
	partner_history[, `:=`(
		rank_largest_decline = frank(share_change_pp, ties.method = "min"),
		rank_largest_increase = frank(-share_change_pp, ties.method = "min"),
		historical_transitions = .N
	), by = partner_keys]
	stopifnot(
		!anyDuplicated(partner_history[, c(partner_keys, "target"), with = FALSE])
	)
	fwrite(
		partner_history[target == 2025],
		file.path(output_dir, "political_trade_us_partners_2025.csv")
	)
	fwrite(validation, file.path(output_dir, "political_trade_checks.csv"))
	write_json(
		list(
			created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
			code = source_record(code_file),
			input_files = lapply(source_files, source_record),
			economy_scope = economy_scope,
			economies = actors,
			allies = allies,
			prespecified_design = list(
				comparisons = c("annual", "two_year", "annual_timing"),
				source_vintages = names(panels),
				supports = names(settings),
				omitted_partners = c("none", "CHE", "CHN"),
				directions = c("exports", "mirror_imports", "symmetric"),
				metrics = c("us_allies", "allied_peers_without_us"),
				weights = c("value_weighted", "equal_country")
			),
			definitions = list(
				us_allies = "Allied-group share of US recorded exports to all destinations, partners recorded exports to the US, or their combined bilateral sum.",
				allied_cohesion = "Other group members share of each allied country trade with all non-US partners. The US is excluded from numerator and denominator, so US movement cannot mechanically raise this measure.",
				imports = "Mirror imports are partners recorded exports to the destination, using nominal FOB export observations. They are not independently reported US imports and are incomplete when a reporter is absent; notably Taiwan has no retained exports-to-US series. Missing directions remain unknown.",
				other_economies = "The comparison comprises economies outside the fourteen-country analytical group, not a comprehensive category of nonallies. Taiwan and other US partners can lie outside this group.",
				symmetric = "Sum of available common export directions. Group totals aggregate country endpoints: an internal directed allied flow enters twice in symmetric allied exposure, once at each endpoint. This is not a share of distinct world trade flows.",
				weights = "Primary group shares divide summed allied flows by summed retained flows. Equal-country shares average the observed country-specific shares; aggregate nominal-growth fields remain value-weighted under both labels.",
				coverage = "Belgium first appears in the retained trade source in 1997; the other thirteen allies appear in 1989. The initially proposed fixed 1992 direction panel could not represent all fourteen allies and was replaced before examining results. The longest available-membership history reports observed-country counts. The fixed thirteen-country history excludes Belgium from every numerator, denominator and reporting direction. Complete fourteen-country comparisons begin in 1997; their first annual change ends 1998.",
				support = "Every comparison retains identical positive reporting directions at its endpoints. The fixed fourteen-country panel additionally requires each direction to be observed in every year 1997-2025. Unobserved countries and missing directions are not assigned zero trade.",
				timing = "Annual timing uses the same directions in target-2, target-1 and target. It compares the prior annual share change with the current annual change and reports current minus prior change. This is descriptive acceleration, not a causal or significance test.",
				twelve_month = "The positive annual 2024/2025 directions having twelve populated 2025 monthly values are retained throughout this historical sensitivity, using refreshed annual values. This selects reporting relationships by 2025 availability; it does not establish historical monthly completeness or unestimated observations.",
				vintage = "Refreshed histories replace 2024 and 2025 coherently, keeping archived 1989-2023 values. Thus comparisons involving 2023 cross the original and refreshed source vintages. A uniformly revised history is unavailable.",
				omissions = "Switzerland and China are removed separately from every year and both ends of every direction. The declared allied membership is otherwise unchanged.",
				ranks = "Largest-decline rank one is the smallest percentage-point share change. Largest-slowdown rank one is the smallest current-minus-prior annual change. Largest-increase ranks are separate. The support and horizon determine the historical comparison count.",
				growth = "Nominal change on retained directions, not real trade or total national trade. Relative log growth compares flows with group members versus other retained flows.",
				country_contributions = "Every core partner uses the same US total and retained directions within a comparison. The partner share changes sum exactly to the total US allied-group share change. Each partner also retains its nominal growth and historical rank; no country is selected for omission based on the observed result.",
				interpretation = "Descriptive concentration of nominal trade in a declared partner group. Changes in shares do not directly measure political support, causal policy responses or comprehensive decoupling."
			),
			raw_source_network_fetch = FALSE
		),
		file.path(output_dir, "political_trade_definition.json"),
		pretty = TRUE,
		auto_unbox = TRUE
	)
	invisible(result)
}

trade_census_pattern = function(root, output_dir, code_file) {
	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
	workbook = file.path(root, "data", "census", "country.xlsx")
	manifest_file = file.path(root, "data", "census", "source_manifest.json")
	acquisition = fromJSON(manifest_file)
	stopifnot(
		identical(digest(file = workbook, algo = "sha256"), acquisition$sha256)
	)
	raw = as.data.table(read_census_table(workbook, sheet = "country"))
	raw[, `:=`(
		year = as.integer(year),
		census_code = sprintf("%04d", as.integer(source_country_code))
	)]
	stopifnot(!anyDuplicated(raw[, .(year, census_code)]))
	mapping = data.table(
		actor = c(
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
			"GBR",
			"CHE",
			"CHN",
			"TWN",
			"HKG",
			"MAC"
		),
		census_code = c(
			"6021",
			"4231",
			"1220",
			"4099",
			"4279",
			"4280",
			"4759",
			"5880",
			"4210",
			"4039",
			"4550",
			"5800",
			"4700",
			"4120",
			"4419",
			"5700",
			"5830",
			"5820",
			"5660"
		)
	)
	mapping[, core_ally := seq_len(.N) <= 14L]
	names_in_source = unique(raw[, .(census_code, source_name = source_country_name)])
	stopifnot(!anyDuplicated(names_in_source$census_code))
	mapping = merge(
		mapping,
		names_in_source,
		by = "census_code",
		all.x = TRUE,
		sort = FALSE
	)
	stopifnot(!anyNA(mapping$source_name))
	fwrite(mapping, file.path(output_dir, "country_mapping.csv"))
	values = rbindlist(
		list(
			raw[, .(
				year,
				census_code,
				source_name = source_country_name,
				flow = "exports",
				value_usd = annual_exports_millions * 1e6
			)],
			raw[, .(
				year,
				census_code,
				source_name = source_country_name,
				flow = "imports",
				value_usd = annual_imports_millions * 1e6
			)]
		)
	)
	stopifnot(all(is.finite(values$value_usd)), all(values$value_usd >= 0))
	values[, complete_calendar_year := year <= 2025L]
	coverage = values[, .(
		source_rows = .N,
		country_rows = sum(as.integer(census_code) >= 1000L),
		core_allies = sum(census_code %in% mapping[core_ally == TRUE, census_code]),
		world_record_present = any(census_code == "0015"),
		complete_calendar_year = all(complete_calendar_year)
	), by = .(year, flow)]
	fwrite(coverage, file.path(output_dir, "coverage.csv"))
	world = values[census_code == "0015" & complete_calendar_year == TRUE, .(year, flow, world_total_usd = value_usd)]
	country_sum = values[as.integer(census_code) >= 1000L & complete_calendar_year == TRUE, .(country_sum_usd = sum(value_usd), country_rows = .N), by = .(year, flow)]
	reconciliation = merge(world, country_sum, by = c("year", "flow"), all.x = TRUE)
	reconciliation[, `:=`(
		difference_usd = country_sum_usd - world_total_usd,
		relative_difference = (country_sum_usd - world_total_usd) / world_total_usd
	)]
	reconciliation[, within_one_part_per_million := is.finite(relative_difference) & abs(relative_difference) <= 1e-6]
	reconciliation[, primary_history := year >= 1992L]
	stopifnot(
		all(reconciliation[primary_history == TRUE, within_one_part_per_million])
	)
	fwrite(reconciliation, file.path(output_dir, "reconciliation.csv"))
	years = 1992:2025
	selected = merge(
		values[year %in% years],
		mapping[, .(census_code, actor, core_ally)],
		by = "census_code"
	)
	stopifnot(
		nrow(selected) == length(years) * nrow(mapping) * 2L,
		!anyNA(selected$value_usd)
	)
	two_way = values[year %in% years, .(value_usd = sum(value_usd), flow_count = .N), by = .(year, census_code, source_name)]
	stopifnot(all(two_way$flow_count == 2L))
	two_way[, flow := "two_way"]
	all_values = rbindlist(
		list(
			values[year %in% years, .(year, census_code, source_name, flow, value_usd)],
			two_way[, .(year, census_code, source_name, flow, value_usd)]
		)
	)
	selected = merge(
		all_values,
		mapping[, .(census_code, actor, core_ally)],
		by = "census_code"
	)
	fwrite(selected, file.path(output_dir, "selected_source_levels.csv"))
	omissions = list(none = character(), CHE = "4419", CHN = "5700")
	levels = list()
	country_levels = list()
	index = 0L
	for (variant in names(omissions)) {
		for (direction in c("exports", "imports", "two_way")) {
			for (current_year in years) {
				x = all_values[year == current_year & flow == direction]
				w = x[census_code == "0015", value_usd]
				omitted = x[census_code %in% omissions[[variant]], value_usd]
				stopifnot(length(w) == 1L, length(omitted) == length(omissions[[variant]]))
				denominator = w - sum(omitted)
				core = merge(
					x,
					mapping[core_ally == TRUE, .(census_code, actor)],
					by = "census_code"
				)
				stopifnot(
					nrow(core) == 14L,
					denominator > 0,
					all(core$value_usd > 0),
					sum(core$value_usd) <= denominator
				)
				index = index + 1L
				levels[[index]] = data.table(
					variant = variant,
					flow = direction,
					year = current_year,
					allied_value_usd = sum(core$value_usd),
					world_total_usd = w,
					omitted_value_usd = sum(omitted),
					denominator_usd = denominator,
					allied_share = sum(core$value_usd) / denominator,
					allies = nrow(core)
				)
				country_levels[[index]] = core[, .(
					variant = variant,
					flow = direction,
					year = current_year,
					actor,
					census_code,
					source_name,
					value_usd,
					denominator_usd = denominator,
					share = value_usd / denominator
				)]
			}
		}
	}
	levels = rbindlist(levels)
	country_levels = rbindlist(country_levels)
	fwrite(levels, file.path(output_dir, "levels.csv"))
	history = list()
	country_changes = list()
	index = 0L
	for (variant_name in names(omissions)) {
		for (direction in c("exports", "imports", "two_way")) {
			for (h in 1:2) {
				for (target in years[years >= min(years) + h]) {
					before = levels[variant == variant_name & flow == direction & year == target - h]
					after = levels[variant == variant_name & flow == direction & year == target]
					earlier = levels[variant == variant_name & flow == direction & year == target - 2L * h]
					previous_share = if (nrow(earlier)) earlier$allied_share else NA_real_
					previous_change = 100 * (before$allied_share - previous_share)
					change = 100 * (after$allied_share - before$allied_share)
					country_before = country_levels[variant == variant_name & flow == direction & year == target - h]
					country_after = country_levels[variant == variant_name & flow == direction & year == target]
					country_earlier = country_levels[variant == variant_name & flow == direction & year == target - 2L * h]
					partners = merge(
						country_before[, .(actor, previous_value_usd = value_usd, previous_share = share)],
						country_after[, .(actor, census_code, source_name, value_usd, share)],
						by = "actor"
					)
					partners = merge(
						partners,
						country_earlier[, .(actor, earlier_share = share)],
						by = "actor",
						all.x = TRUE
					)
					partners[, `:=`(
						variant = variant_name,
						flow = direction,
						horizon = h,
						previous_year = target - h,
						year = target,
						earlier_year = if (nrow(earlier)) target - 2L * h else NA_integer_,
						change_pp = 100 * (share - previous_share),
						prior_change_pp = 100 * (previous_share - earlier_share),
						value_growth_percent = 100 * (value_usd / previous_value_usd - 1)
					)]
					partners[, acceleration_pp := change_pp - prior_change_pp]
					index = index + 1L
					country_changes[[index]] = partners
					history[[index]] = data.table(
						variant = variant_name,
						flow = direction,
						horizon = h,
						previous_year = target - h,
						year = target,
						earlier_year = if (nrow(earlier)) target - 2L * h else NA_integer_,
						previous_allied_share = before$allied_share,
						allied_share = after$allied_share,
						change_pp = change,
						prior_change_pp = previous_change,
						acceleration_pp = change - previous_change,
						allies = 14L,
						allies_share_declines = sum(partners$change_pp < 0),
						allies_nominal_declines = sum(partners$value_growth_percent < 0),
						allied_nominal_growth_percent = 100 * (after$allied_value_usd / before$allied_value_usd - 1),
						denominator_growth_percent = 100 * (after$denominator_usd / before$denominator_usd - 1),
						previous_allied_value_usd = before$allied_value_usd,
						allied_value_usd = after$allied_value_usd,
						previous_denominator_usd = before$denominator_usd,
						denominator_usd = after$denominator_usd
					)
				}
			}
		}
	}
	history = rbindlist(history)
	country_changes = rbindlist(country_changes)
	history[, `:=`(
		decline_rank = vapply(change_pp, function(v) 1L + sum(change_pp < v - 1e-12), integer(1)),
		change_history_n = .N,
		acceleration_rank_most_negative = vapply(
			acceleration_pp,
			function(v) if (is.finite(v)) 1L + sum(acceleration_pp < v - 1e-12, na.rm = TRUE) else NA_integer_,
			integer(1)
		),
		acceleration_history_n = sum(is.finite(acceleration_pp))
	), by = .(variant, flow, horizon)]
	country_changes[, `:=`(
		decline_rank = vapply(change_pp, function(v) 1L + sum(change_pp < v - 1e-12), integer(1)),
		change_history_n = .N
	), by = .(variant, flow, horizon, actor)]
	identities = country_changes[, .(
		country_sum_change_pp = sum(change_pp),
		country_sum_acceleration_pp = if (all(is.finite(acceleration_pp))) sum(acceleration_pp) else NA_real_
	), by = .(variant, flow, horizon, year)]
	checks = merge(history, identities, by = c("variant", "flow", "horizon", "year"))
	checks = checks[, .(
		variant,
		flow,
		horizon,
		year,
		contribution_error_pp = abs(change_pp - country_sum_change_pp),
		acceleration_contribution_error_pp = abs(acceleration_pp - country_sum_acceleration_pp),
		all_fourteen_present = allies == 14L,
		share_is_probability = previous_allied_share >= 0 & previous_allied_share <= 1 & allied_share >= 0 & allied_share <= 1,
		historical_count_correct = change_history_n == 34L - horizon & acceleration_history_n == 34L - 2L * horizon
	)]
	checks[, passed := contribution_error_pp < 1e-10 & (is.na(acceleration_contribution_error_pp) | acceleration_contribution_error_pp < 1e-10) & all_fourteen_present & share_is_probability & historical_count_correct]
	stopifnot(all(checks$passed))
	fwrite(history, file.path(output_dir, "history.csv"))
	fwrite(history[year == 2025L], file.path(output_dir, "terminal.csv"))
	fwrite(country_changes, file.path(output_dir, "country_changes.csv"))
	fwrite(
		country_changes[year == 2025L],
		file.path(output_dir, "country_changes_2025.csv")
	)
	fwrite(checks, file.path(output_dir, "checks.csv"))
	write_json(
		list(
			code = source_record(code_file),
			inputs = lapply(c(workbook, manifest_file), source_record),
			acquisition = acquisition,
			sheet = "country",
			first_level_year = 1992L,
			last_complete_year = 2025L,
			world_code = "0015",
			allies = mapping[core_ally == TRUE, actor],
			omitted_partner_codes = omissions,
			definitions = list(
				values = "US Census Bureau nominal merchandise trade, millions of US dollars in the workbook converted here to dollars. Total exports are FAS values; general imports are Customs values. Census basis, not seasonally adjusted. Both are US-published directions, unlike IMF partner-reported exports into the US.",
				denominator = "Published world annual total, minus the named omitted economy where applicable. Individual partner totals are never substituted for the published world total; all fourteen numerator series are present in every primary year.",
				coverage = "Source rows before 1992 are retained in coverage and country/world reconciliation diagnostics. Ten direction-year totals in 1987-1991 do not reconcile within one part per million. From 1992 both directions reconcile within that tolerance; small residuals, including imports in 2004-2007, remain recorded rather than forced to zero. Partial 2026 is excluded.",
				horizons = "For horizon h, current change is S[t]-S[t-h], prior change is S[t-h]-S[t-2h], and acceleration is their difference, all in percentage points. Annual timing therefore compares three consecutive years. The two-year timing statistic compares nonoverlapping two-year intervals.",
				ranks = "Within each fixed direction, horizon and partner-omission history; rank one is the largest decline or most negative acceleration. The primary annual history has 33 changes and 32 accelerations, and the two-year history has 32 changes and 30 accelerations.",
				country_contributions = "Each core ally flow is divided by the identical official-world or omission-adjusted denominator. Country percentage-point changes sum exactly to the group change. Nominal declines count reported dollar decreases, not share decreases.",
				interpretation = "Descriptive US merchandise-trade composition. Source publication is distinct from IMF but the underlying observations need not be statistically independent; US exports to Canada are compiled using Canadian import records. These data do not measure bilateral flows among the other allies, price-adjusted volume, services, or causal policy effects."
			),
			network_fetch = FALSE
		),
		file.path(output_dir, "definition.json"),
		pretty = TRUE,
		auto_unbox = TRUE
	)
	invisible(history)
}

trade_pipeline = function(root, code_file) {
	started = Sys.time()
	data_dir = file.path(root, "data")
	output_dir = file.path(root, "results", "trade")
	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
	source_file = file.path(data_dir, "trade.rda")
	actor_file = file.path(data_dir, "actors.txt")
	actors = readLines(actor_file)
	stopifnot(length(actors) > 1L, !anyDuplicated(actors), all(nzchar(actors)))
	input_environment = new.env(parent = emptyenv())
	load(source_file, envir = input_environment)
	trade = as.data.table(input_environment$trade_long)
	baseline = trade_analysis(
		trade,
		actors,
		file.path(output_dir, "baseline"),
		c(source_file, actor_file),
		"cached annual observations",
		code_file
	)
	trade_actor_coverage(trade, actors, output_dir)
	# the parser verifies archived bytes and never fetches in this call
	python = Sys.getenv("PYTHON", "python3")
	status = system2(
		python,
		c(shQuote(file.path(root, "code", "03_prepare_imf_trade.py")), "--root", shQuote(root))
	)
	stopifnot(identical(status, 0L))
	annual = read_analysis_table(file.path(output_dir, "imf", "imf_current_annual.csv"))
	coverage = read_analysis_table(file.path(output_dir, "imf", "imf_monthly_coverage.csv"))
	columns = c("reporter", "partner", "year", "exports_usd")
	stopifnot(
		all(annual$year %in% 2024:2025),
		!anyDuplicated(annual[, .(reporter, partner, year)]),
		all(is.finite(annual$exports_usd)),
		all(annual$exports_usd > 0)
	)
	fresh_trade = rbindlist(list(trade[year < 2024], annual[, ..columns]), use.names = TRUE)
	terminal_support = merge(
		annual[year == 2024, .(reporter, partner)],
		annual[year == 2025, .(reporter, partner)],
		by = c("reporter", "partner")
	)
	twelve_support = merge(
		terminal_support,
		coverage[months == 12, .(reporter, partner)],
		by = c("reporter", "partner")
	)
	stopifnot(!anyDuplicated(twelve_support))
	twelve_annual = merge(annual[, ..columns], twelve_support, by = c("reporter", "partner"))
	twelve_trade = rbindlist(
		list(trade[year < 2024], twelve_annual[, ..columns]),
		use.names = TRUE
	)
	refresh_sources = c(
		source_file,
		actor_file,
		file.path(
			data_dir,
			"imf",
			c("imf_all_annual.json", "imf_all_monthly_2025.json", "source_manifest.json")
		)
	)
	fresh = trade_analysis(
		fresh_trade,
		actors,
		file.path(output_dir, "fresh"),
		refresh_sources,
		"refreshed 2024 and 2025 annual observations",
		code_file
	)
	twelve = trade_analysis(
		twelve_trade,
		actors,
		file.path(output_dir, "twelve_month"),
		refresh_sources,
		"refreshed annual observations on directions with twelve populated months in 2025 and annual values in both terminal years",
		code_file
	)
	inputs = list(
		fresh_terminal = fresh$model_inputs,
		twelve_month_terminal = twelve$model_inputs
	)
	checks = rbindlist(
		lapply(
			names(inputs),
			function(label) {
				rbindlist(
					lapply(
						names(inputs[[label]]),
						function(year) {
							m = inputs[[label]][[year]]
							data.table(
								input = label,
								year = as.integer(year),
								actor_order = identical(rownames(m), actors) && identical(colnames(m), actors),
								symmetric = identical(m, t(m)),
								diagonal_missing = all(is.na(diag(m))),
								earlier_year_unchanged = if (as.integer(year) <= 2023) identical(m, baseline$model_inputs[[year]]) else NA
							)
						}
					)
				)
			}
		)
	)
	stopifnot(
		all(checks$actor_order),
		all(checks$symmetric),
		all(checks$diagonal_missing),
		all(checks$earlier_year_unchanged, na.rm = TRUE)
	)
	saveRDS(inputs, file.path(output_dir, "trade_refresh_inputs.rds"))
	fwrite(checks, file.path(output_dir, "input_checks.csv"))
	annual_acquisition = fromJSON(file.path(data_dir, "imf", "imf_all_annual.manifest.json"))
	source_date = as.Date(substr(annual_acquisition$retrieved_utc, 1L, 10L))
	manifest = rbindlist(
		lapply(
			names(inputs),
			function(label) {
				data.table(
					input = label,
					input_hash = digest(inputs[[label]], algo = "sha256"),
					earlier_years_identical = TRUE,
					years_replaced = "2024,2025",
					source_vintage = as.character(source_date),
					years_original = "1989-2023"
				)
			}
		)
	)
	fwrite(manifest, file.path(output_dir, "trade_refresh_manifest.csv"))
	trade_partner_sensitivity(
		trade,
		fresh_trade,
		actors,
		output_dir,
		refresh_sources,
		code_file
	)
	trade_political_pattern(
		trade,
		fresh_trade,
		twelve_support,
		actors,
		output_dir,
		refresh_sources,
		code_file
	)
	all_annual = read_analysis_table(
		file.path(output_dir, "imf", "imf_all_economies_annual.csv")
	)
	all_monthly = read_analysis_table(
		file.path(output_dir, "imf", "imf_all_economies_monthly_coverage.csv")
	)
	all_codes = sort(
		unique(
			c(trade$reporter, trade$partner, all_annual$reporter, all_annual$partner)
		)
	)
	economy_codes = data.table(
		code = all_codes,
		included = grepl("^[A-Z]{3}$", all_codes),
		in_political_grid = all_codes %in% actors
	)
	economy_codes[, rule := ifelse(
		included,
		"three-letter source economy; includes territories and historical economy codes",
		"source aggregate or special non-economy code"
	)]
	fwrite(
		economy_codes,
		file.path(output_dir, "political_trade_economy_codes.csv")
	)
	economies = economy_codes[included == TRUE, code]
	all_annual = all_annual[reporter %in% economies & partner %in% economies & reporter != partner, ..columns]
	stopifnot(
		all(is.finite(all_annual$exports_usd) & all_annual$exports_usd > 0),
		!anyDuplicated(all_annual[, .(reporter, partner, year)])
	)
	all_fresh = rbindlist(list(trade[year < 2024], all_annual), use.names = TRUE)
	all_terminal = merge(
		all_annual[year == 2024, .(reporter, partner)],
		all_annual[year == 2025, .(reporter, partner)],
		by = c("reporter", "partner")
	)
	all_twelve = merge(
		all_terminal,
		all_monthly[months == 12, .(reporter, partner)],
		by = c("reporter", "partner")
	)
	trade_political_pattern(
		trade,
		all_fresh,
		all_twelve,
		economies,
		file.path(output_dir, "all_economies"),
		refresh_sources,
		code_file,
		economy_scope = "available_source_economies"
	)
	trade_census_pattern(root, file.path(output_dir, "census"), code_file)
	selected = rbindlist(
		list(
			baseline$history[year %in% c(1991, 2008, 2009, 2020, 2024, 2025)][, vintage := "cached"],
			fresh$history[year == 2025][, vintage := "fresh"],
			twelve$history[year == 2025][, vintage := "twelve_month"]
		),
		use.names = TRUE
	)
	fwrite(selected, file.path(output_dir, "selected_history.csv"))
	comparison = rbindlist(
		list(
			baseline$history[year %in% c(2009, 2020, 2025)][, source := "Cached May 2026"],
			fresh$history[year == 2025][, source := paste("Fresh", format(source_date, "%B %Y"))],
			twelve$history[year == 2025][, source := "Fresh; twelve populated months in 2025"]
		)
	)
	comparison = comparison[, .(
		year,
		source,
		common_directions,
		common_bilateral_pairs = common_pairs,
		nominal_change_percent = 100 * (total_ratio - 1),
		mean_absolute_log_change = mean_abs_log_change,
		previous_total_usd = directed_previous_total,
		current_total_usd = directed_current_total
	)]
	fwrite(comparison, file.path(output_dir, "economic_comparison.csv"))
	write_json(
		list(
			title = "Nominal exports and bilateral changes on common reporting directions",
			note = "Each row uses export directions with positive annual values in both adjacent years. Bilateral values sum only those common directions. The final row additionally requires twelve populated monthly observations in 2025; it does not certify 2024 monthly coverage or impute zero for excluded months. Nominal growth is a matched-panel total, not world trade growth, and mean absolute log change is an unweighted average across supported bilateral pairs."
		),
		file.path(output_dir, "economic_comparison_note.json"),
		pretty = TRUE,
		auto_unbox = TRUE
	)
	write_json(
		list(
			created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
			elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
			code = source_record(code_file),
			parser = source_record(file.path(root, "code", "03_prepare_imf_trade.py")),
			input_files = lapply(refresh_sources, source_record),
			archived_first_year = min(trade$year),
			archived_last_year = max(trade$year),
			refresh_years = c(2024, 2025),
			common_twelve_month_directions = nrow(twelve_support),
			definitions = list(
				values = "Original IMF dollar observation values; SCALE metadata is not an instruction to multiply these returned observations again.",
				twelve_month = "Twelve distinct nonmissing monthly observations in 2025, including explicit numerical zeros if any, and positive annual values in both 2024 and 2025. This does not establish twelve populated months in 2024. Missing monthly observations are not set to zero.",
				flags = "STATUS values are retained literally: the fetched scheme provides no authoritative meaning for e or r. Blank STATUS does not establish complete direct reporting. DERIVATION_TYPE M denotes Mixed-type data. Twelve populated monthly values do not imply unestimated observations.",
				vintage = "The cached panel and refreshed terminal years represent different source vintages. Current derivation/status flags cannot reconstruct metadata discarded from the older cache.",
				model_inputs = "Log of the sum of observed export directions for each unordered country pair; one direction is sufficient, both missing gives NA, diagonal remains NA. These are reconstructed data, not fitted model results."
			),
			raw_source_network_fetch = FALSE
		),
		file.path(output_dir, "provenance.json"),
		pretty = TRUE,
		auto_unbox = TRUE
	)
	invisible(inputs)
}

trade_census_products = function(input_xlsx, country_xlsx, mapping_csv, out) {
	data.table::setDTthreads(1L)
	dir.create(out, recursive = TRUE, showWarnings = FALSE)
	d = data.table::as.data.table(
		read_census_table(
			input_xlsx,
			col_types = c(rep("text", 4), rep("numeric", 10))
		)
	)
	stopifnot(
		!anyDuplicated(d[, .(source_country_code, enduse_category_code)]),
		all(grepl("^[0-5][0-9]{4}$", d$enduse_category_code))
	)
	value_columns = paste0("value_", 16:25)
	stopifnot(
		identical(
			names(d),
			c(
				"source_country_code",
				"source_country_name",
				"enduse_category_code",
				"enduse_category_description",
				value_columns
			)
		),
		length(unique(d$enduse_category_code)) == 142L,
		length(unique(d$source_country_code)) == 235L
	)
	stopifnot(
		!anyNA(d[, ..value_columns]),
		all(as.matrix(d[, ..value_columns]) >= 0),
		all(as.matrix(d[, ..value_columns]) %% 1 == 0)
	)
	t = data.table::melt(
		d,
		id.vars = c(
			"source_country_code",
			"source_country_name",
			"enduse_category_code",
			"enduse_category_description"
		),
		measure.vars = value_columns,
		variable.name = "year",
		value.name = "value_usd"
	)
	data.table::setnames(
		t,
		c(
			"source_country_code",
			"source_country_name",
			"enduse_category_code",
			"enduse_category_description"
		),
		c("country_code", "country_name", "enduse_code", "description")
	)
	t[, year := 2000L + as.integer(sub("value_", "", year))]
	t[, observation_status := ifelse(value_usd == 0, "reported_zero", "reported_positive")]
	data.table::fwrite(t, file.path(out, "tidy.csv"))
	world = t[country_code == "0000", .(world_usd = sum(value_usd)), by = year]
	ws = t[country_code == "0000", .(year, enduse_code, world_product = value_usd)]
	ps = t[country_code != "0000", .(partner_sum = sum(value_usd)), by = .(year, enduse_code)]
	test = merge(ws, ps, by = c("year", "enduse_code"))
	stopifnot(nrow(test) == 1420L, all(test$world_product == test$partner_sum))
	cc = data.table::as.data.table(read_census_table(country_xlsx))
	cc[, `:=`(
		source_country_code = sprintf("%04d", as.integer(source_country_code)),
		year = as.integer(year),
		annual_imports_millions = as.numeric(annual_imports_millions)
	)]
	stopifnot(!anyDuplicated(cc[, .(source_country_code, year)]))
	totals = t[, .(product_total = sum(value_usd)), by = .(country_code, year)]
	totals[, source_country_code := ifelse(country_code == "0000", "0015", country_code)]
	reconciliation = merge(
		totals,
		cc[year %in% 2016:2025, .(
			source_country_code,
			year,
			country_total_usd = annual_imports_millions * 1e6
		)],
		by = c("source_country_code", "year"),
		all.x = TRUE
	)
	reconciliation[, `:=`(
		difference_usd = product_total - country_total_usd,
		rounding_bound_usd = 0.501,
		comparison_status = ifelse(
			is.na(country_total_usd),
			"country_total_absent_product_sum_reported_zero",
			"available"
		)
	)]
	stopifnot(
		all(reconciliation[is.na(country_total_usd), product_total] == 0),
		all(abs(reconciliation[!is.na(country_total_usd), difference_usd]) <= 0.501)
	)
	data.table::fwrite(reconciliation, file.path(out, "reconciliation.csv"))
	m = read_analysis_table(mapping_csv)
	m[, census_code := sprintf("%04d", as.integer(census_code))]
	allies = m[core_ally == TRUE, census_code]
	stopifnot(
		length(allies) == 14L,
		!anyDuplicated(allies),
		setequal(
			allies,
			c(
				"6021",
				"4231",
				"1220",
				"4099",
				"4279",
				"4280",
				"4759",
				"5880",
				"4210",
				"4039",
				"4550",
				"5800",
				"4700",
				"4120"
			)
		)
	)
	groups = list(
		fourteen = allies,
		canada = "1220",
		other_thirteen = setdiff(allies, "1220")
	)
	a = t[country_code %in% allies]
	a[, sector := substr(enduse_code, 1, 1)]
	contributions = list()
	components = list()
	history = list()
	leaveout = list()
	n = 0L
	k = 0L
	j = 0L
	for (h in 1:2) {
		for (y in seq.int(2016 + h, 2025)) {
			b = y - h
			w0 = world[year == b, world_usd]
			w1 = world[year == y, world_usd]
			z = merge(
				a[year == b, .(
					country_code,
					country_name,
					enduse_code,
					description,
					sector,
					value_start = value_usd
				)],
				a[year == y, .(country_code, enduse_code, value_end = value_usd)],
				by = c("country_code", "enduse_code"),
				all = TRUE
			)
			stopifnot(!anyNA(z))
			z[, `:=`(
				start_year = b,
				end_year = y,
				horizon = h,
				contribution_pp = 100 * (value_end / w1 - value_start / w0)
			)]
			n = n + 1L
			contributions[[n]] = z
			loo = z[, .(omitted_contribution_pp = sum(contribution_pp)), by = .(country_code, country_name)]
			loo[, `:=`(
				start_year = b,
				end_year = y,
				horizon = h,
				remaining_change_pp = sum(z$contribution_pp) - omitted_contribution_pp,
				denominator = "full world unchanged"
			)]
			leaveout[[n]] = loo
			for (g in names(groups)) {
				cs = groups[[g]]
				q = z[country_code %in% cs]
				ct = q[, .(v0 = sum(value_start), v1 = sum(value_end), delta = sum(contribution_pp)), by = country_code]
				j = j + 1L
				history[[j]] = data.table::data.table(
					group = g,
					start_year = b,
					end_year = y,
					horizon = h,
					share_start = 100 * sum(q$value_start) / w0,
					share_end = 100 * sum(q$value_end) / w1,
					change_pp = sum(q$contribution_pp),
					nominal_change_pct = 100 * (sum(q$value_end) / sum(q$value_start) - 1),
					share_decline_countries = sum(ct$delta < 0),
					dollar_decline_countries = sum(ct$v1 < ct$v0)
				)
				for (level in c("principal", "five_digit")) {
					at = data.table::copy(a[country_code %in% cs & year %in% c(b, y)])
					wt = data.table::copy(t[country_code == "0000" & year %in% c(b, y)])
					at[, code := if (level == "principal") substr(enduse_code, 1, 1) else enduse_code]
					wt[, code := if (level == "principal") substr(enduse_code, 1, 1) else enduse_code]
					ax = at[, .(allied_imports_usd = sum(value_usd)), by = .(year, code)]
					wx = wt[, .(world_imports_usd = sum(value_usd)), by = .(year, code)]
					v = merge(wx, ax, by = c("year", "code"), all.y = TRUE)
					stopifnot(!anyNA(v))
					p = merge(
						v[year == b, .(
							sector = code,
							world_start = world_imports_usd,
							allied_start = allied_imports_usd
						)],
						v[year == y, .(
							sector = code,
							world_end = world_imports_usd,
							allied_end = allied_imports_usd
						)],
						by = "sector",
						all = TRUE
					)
					stopifnot(!anyNA(p), all(p$world_start > 0), all(p$world_end > 0))
					p[, `:=`(
						within_sector_pp = 100 * (world_start / w0 + world_end / w1) / 2 * (allied_end / world_end - allied_start / world_start),
						sector_composition_pp = 100 * (allied_start / world_start + allied_end / world_end) / 2 * (world_end / w1 - world_start / w0),
						contribution_pp = 100 * (allied_end / w1 - allied_start / w0),
						entry_exit_pp = 0,
						group = g,
						start_year = b,
						end_year = y,
						horizon = h,
						classification = level
					)]
					p[, residual_pp := contribution_pp - within_sector_pp - sector_composition_pp]
					stopifnot(
						max(abs(p$residual_pp)) < 1e-8,
						abs(sum(p$contribution_pp) - sum(q$contribution_pp)) < 1e-8
					)
					k = k + 1L
					components[[k]] = p
				}
			}
		}
	}
	hh = data.table::rbindlist(history)
	data.table::setorder(hh, group, horizon, end_year)
	hh[, `:=`(
		decline_rank = data.table::frank(change_pp, ties.method = "min"),
		comparison_count = .N,
		previous_change_pp = data.table::shift(change_pp, n = horizon[1])
	), by = .(group, horizon)]
	hh[, change_from_previous_pp := change_pp - previous_change_pp]
	data.table::fwrite(hh, file.path(out, "history.csv"))
	data.table::fwrite(
		data.table::rbindlist(contributions),
		file.path(out, "contributions.csv")
	)
	data.table::fwrite(
		data.table::rbindlist(components),
		file.path(out, "sector_components.csv")
	)
	data.table::fwrite(
		data.table::rbindlist(leaveout),
		file.path(out, "leave_one_out.csv")
	)
	invisible(hh)
}

trade_census_timing = function(country_xlsx, mapping_csv, out) {
	data.table::setDTthreads(1L)
	dir.create(out, recursive = TRUE, showWarnings = FALSE)
	c = data.table::as.data.table(read_census_table(country_xlsx))
	c[, `:=`(
		source_country_code = sprintf("%04d", as.integer(source_country_code)),
		year = as.integer(year)
	)]
	m = read_analysis_table(mapping_csv)
	allies = sprintf("%04d", as.integer(m[core_ally == TRUE, census_code]))
	columns = paste0("imports_", tolower(month.abb), "_millions")
	stopifnot(
		length(allies) == 14L,
		!anyDuplicated(allies),
		!anyDuplicated(c[, .(source_country_code, year)])
	)
	c = c[year %in% 2023:2025 & source_country_code %in% c(allies, "0015")]
	c[, (c(columns, "annual_imports_millions")) := lapply(.SD, as.numeric), .SDcols = c(columns, "annual_imports_millions")]
	stopifnot(nrow(c) == 45L, !anyNA(c[, ..columns]))
	rec = c[, .(
		year,
		source_country_code,
		annual_imports_millions,
		monthly_sum_millions = rowSums(.SD)
	), .SDcols = columns]
	rec[, `:=`(
		difference_usd = (monthly_sum_millions - annual_imports_millions) * 1e6,
		rounding_bound_usd = 6.501
	)]
	stopifnot(all(abs(rec$difference_usd) <= rec$rounding_bound_usd))
	data.table::fwrite(rec, file.path(out, "monthly_reconciliation.csv"))
	t = data.table::melt(
		c,
		id.vars = c("year", "source_country_code"),
		measure.vars = columns,
		variable.name = "month",
		value.name = "value_millions"
	)
	t[, month := match(month, columns)]
	t[, date := data.table::as.IDate(sprintf("%04d-%02d-01", year, month))]
	stopifnot(
		!anyDuplicated(t[, .(source_country_code, year, month)]),
		all(is.finite(t$value_millions))
	)
	w = t[source_country_code == "0015", .(date, world_usd = value_millions * 1e6)]
	stopifnot(all(is.finite(w$world_usd)), all(w$world_usd > 0))
	groups = list(
		fourteen = allies,
		canada = "1220",
		other_thirteen = setdiff(allies, "1220")
	)
	result = list()
	for (g in names(groups)) {
		a = t[source_country_code %in% groups[[g]], .(numerator_usd = sum(value_millions) * 1e6), by = date]
		z = merge(a, w, by = "date")
		data.table::setorder(z, date)
		z[, `:=`(
			group = g,
			monthly_share_pct = 100 * numerator_usd / world_usd,
			trailing12_share_pct = 100 * data.table::frollsum(numerator_usd, 12) / data.table::frollsum(world_usd, 12)
		)]
		z[, `:=`(
			same_month_change_pp = monthly_share_pct - data.table::shift(monthly_share_pct, 12),
			trailing12_change_pp = trailing12_share_pct - data.table::shift(trailing12_share_pct, 12)
		)]
		result[[g]] = z
	}
	all = data.table::rbindlist(result)
	data.table::fwrite(all, file.path(out, "timing.csv"))
	all[, `:=`(
		year = as.integer(format(date, "%Y")),
		quarter = (as.integer(format(date, "%m")) - 1L) %/% 3L + 1L
	)]
	qq = all[, .(numerator_usd = sum(numerator_usd), world_usd = sum(world_usd)), by = .(group, year, quarter)]
	qq[, share_pct := 100 * numerator_usd / world_usd]
	data.table::setorder(qq, group, quarter, year)
	qq[, same_quarter_change_pp := share_pct - data.table::shift(share_pct), by = .(group, quarter)]
	data.table::fwrite(qq, file.path(out, "quarterly_timing.csv"))
	invisible(result)
}

trade_price_coverage = function(raw_dir, out) {
	dir.create(out, recursive = TRUE, showWarnings = FALSE)
	rows = list()
	k = 0L
	label_checks = list()
	labels = c(
		"Foods, feeds, and beverages",
		"Industrial supplies and materials",
		"Capital goods",
		"Automotive vehicles, parts and engines",
		"Consumer goods excluding automotives"
	)
	for (i in 0:4) {
		id = paste0("EIUIR", i)
		response = jsonlite::fromJSON(
			file.path(raw_dir, paste0("bls_", id, ".json")),
			simplifyVector = FALSE
		)
		stopifnot(
			identical(response$status, "REQUEST_SUCCEEDED"),
			length(response$Results$series) == 1L,
			identical(response$Results$series[[1]]$seriesID, id)
		)
		a = response$Results$series[[1]]$data
		v = data.table::rbindlist(
			lapply(
				a,
				function(x) {
					data.table::data.table(
						year = as.integer(x$year),
						period = x$period,
						value = as.numeric(x$value)
					)
				}
			)
		)
		stopifnot(
			!anyDuplicated(v[, .(year, period)]),
			all(is.finite(v$value)),
			all(v$value > 0)
		)
		fred_id = paste0("IR", i)
		mirror = read_analysis_table(file.path(raw_dir, paste0("fred_", fred_id, ".csv")))
		setnames(mirror, fred_id, "mirror_value")
		mirror[, mirror_value := suppressWarnings(as.numeric(mirror_value))]
		mirror[, date := as.character(observation_date)]
		bv = v[period %in% sprintf("M%02d", 1:12), .(date = paste0(year, "-", substring(period, 2), "-01"), value)]
		both = merge(bv, mirror[, .(date, mirror_value)], by = "date", all = TRUE)
		stopifnot(
			!anyDuplicated(mirror$date),
			all(is.na(both$value) == is.na(both$mirror_value)),
			all(both[!is.na(value), value] == both[!is.na(value), mirror_value])
		)
		html = paste(
			readLines(file.path(raw_dir, paste0("fred_", fred_id, ".html")), warn = FALSE),
			collapse = " "
		)
		title = regmatches(html, regexpr("<title>.*?</title>", html, perl = TRUE))
		stopifnot(length(title) == 1L)
		title = gsub("</?title>", "", title)
		label_checks[[id]] = data.table(
			bls_series = id,
			fred_series = fred_id,
			fred_page_title = title,
			matched_months = sum(!is.na(both$value)),
			observed_values_equal = TRUE,
			missing_months_equal = TRUE
		)
		for (y in 2023:2025) {
			q = v[year == y & period %in% sprintf("M%02d", 1:12)]
			stopifnot(!anyDuplicated(q$period))
			k = k + 1L
			rows[[k]] = data.table::data.table(
				series_id = id,
				description = labels[i + 1],
				enduse_sector = as.character(i),
				year = y,
				observed_months = nrow(q),
				missing_months = paste(setdiff(sprintf("M%02d", 1:12), q$period), collapse = ";"),
				annual_average = if (nrow(q) == 12) mean(q$value) else NA_real_,
				annual_status = if (nrow(q) == 12) "complete" else "unavailable_no_imputation"
			)
		}
	}
	f = read_analysis_table(file.path(raw_dir, "fred_CANTOT.csv"))
	f[, `:=`(
		year = as.integer(substr(observation_date, 1, 4)),
		month = as.integer(substr(observation_date, 6, 7)),
		canada_import_price_index = suppressWarnings(as.numeric(canada_import_price_index))
	)]
	stopifnot(
		!anyDuplicated(f[, .(year, month)]),
		all(
			is.finite(f[!is.na(canada_import_price_index), canada_import_price_index])
		),
		all(f[!is.na(canada_import_price_index), canada_import_price_index] > 0)
	)
	for (y in 2023:2025) {
		q = f[year == y & !is.na(canada_import_price_index)]
		k = k + 1L
		rows[[k]] = data.table::data.table(
			series_id = "CANTOT",
			description = "Import Price Index by Origin (NAICS): All Industries for Canada",
			enduse_sector = "Canada all",
			year = y,
			observed_months = nrow(q),
			missing_months = paste(sprintf("M%02d", setdiff(1:12, q$month)), collapse = ";"),
			annual_average = if (nrow(q) == 12) mean(q$canada_import_price_index) else NA_real_,
			annual_status = if (nrow(q) == 12) "complete" else "unavailable_no_imputation"
		)
	}
	data.table::fwrite(
		data.table::rbindlist(label_checks),
		file.path(out, "sector_label_verification.csv")
	)
	data.table::fwrite(
		data.table::rbindlist(rows),
		file.path(out, "price_coverage.csv")
	)
	invisible(rows)
}

trade_census_detail_pipeline = function(root, code_file) {
	out = file.path(root, "results", "trade", "census", "detail")
	raw = file.path(root, "data", "census", "detail")
	prices = file.path(root, "data", "census", "prices")
	mapping = file.path(root, "results", "trade", "census", "country_mapping.csv")
	country = file.path(root, "data", "census", "country.xlsx")
	trade_census_products(
		file.path(raw, "enduse_imports.xlsx"),
		country,
		mapping,
		out
	)
	trade_census_timing(country, mapping, out)
	trade_price_coverage(prices, out)
	# compact summaries are derived from exhaustive retained accounting tables.
	components = read_analysis_table(file.path(out, "sector_components.csv"))
	summary = components[, .(
		within_sector_pp = sum(within_sector_pp),
		sector_composition_pp = sum(sector_composition_pp),
		entry_exit_pp = sum(entry_exit_pp),
		change_pp = sum(contribution_pp),
		algebraic_residual_pp = sum(residual_pp),
		reported_categories = .N
	), by = .(group, start_year, end_year, horizon, classification)]
	fwrite(summary, file.path(out, "component_summary.csv"))
	contributions = read_analysis_table(
		file.path(out, "contributions.csv"),
		colClasses = c(country_code = "character", enduse_code = "character")
	)
	concentration = contributions[, .(contribution_pp = sum(contribution_pp)), by = .(start_year, end_year, horizon, enduse_code, description)]
	fwrite(concentration, file.path(out, "product_concentration.csv"))
	major = concentration[, .(
		passenger_cars_pp = sum(contribution_pp[enduse_code == "30000"]),
		crude_oil_pp = sum(contribution_pp[enduse_code == "10000"]),
		total_change_pp = sum(contribution_pp)
	), by = .(start_year, end_year, horizon)]
	major[, `:=`(
		cars_and_oil_pp = passenger_cars_pp + crude_oil_pp,
		cars_and_oil_fraction_net_change = (passenger_cars_pp + crude_oil_pp) / total_change_pp
	)]
	fwrite(major, file.path(out, "cars_oil_summary.csv"))
	countries = contributions[, .(
		value_start = sum(value_start),
		value_end = sum(value_end),
		contribution_pp = sum(contribution_pp)
	), by = .(start_year, end_year, horizon, country_code, country_name)]
	countries[, nominal_change_pct := 100 * (value_end / value_start - 1)]
	fwrite(countries, file.path(out, "country_contributions.csv"))
	history = read_analysis_table(file.path(out, "history.csv"))
	fwrite(
		history[end_year %in% c(2024L, 2025L)],
		file.path(out, "terminal_summary.csv")
	)
	tidy = read_analysis_table(
		file.path(out, "tidy.csv"),
		colClasses = c(country_code = "character", enduse_code = "character")
	)
	scope = tidy[country_code %in% c("5830", "4190", "4419"), .(value_usd = sum(value_usd)), by = .(country_code, country_name, year)]
	world = tidy[country_code == "0000", .(world_usd = sum(value_usd)), by = year]
	scope = merge(scope, world, by = "year")
	scope[, share_pct := 100 * value_usd / world_usd]
	scope_rows = list()
	k = 0L
	for (h in 1:2) {
		z = copy(scope)
		setorder(z, country_code, year)
		z[, `:=`(
			share_change_pp = share_pct - shift(share_pct, h),
			nominal_change_pct = 100 * (value_usd / shift(value_usd, h) - 1)
		), by = country_code]
		z[, `:=`(start_year = year - h, end_year = year, horizon = h)]
		k = k + 1L
		scope_rows[[k]] = z[year >= 2016 + h]
	}
	fwrite(rbindlist(scope_rows), file.path(out, "scope_comparison.csv"))
	price = read_analysis_table(file.path(out, "price_coverage.csv"))
	tidy[, sector := substr(enduse_code, 1, 1)]
	mapping_data = read_analysis_table(mapping)
	allies = sprintf("%04d", as.integer(mapping_data[core_ally == TRUE, census_code]))
	coverage = list()
	k = 0L
	for (y in 2023:2025) {
		for (g in c("fourteen", "world")) {
			codes = if (g == "fourteen") allies else "0000"
			z = tidy[year == y & country_code %in% codes]
			available = price[year == y & annual_status == "complete" & enduse_sector %in% as.character(0:5), enduse_sector]
			k = k + 1L
			coverage[[k]] = data.table(
				year = y,
				population = g,
				nominal_total_usd = sum(z$value_usd),
				sectors_with_complete_prices = paste(sort(available), collapse = ";"),
				complete_price_sector_value_usd = sum(z[sector %in% available, value_usd]),
				coverage_pct = 100 * sum(z[sector %in% available, value_usd]) / sum(z$value_usd)
			)
		}
	}
	fwrite(rbindlist(coverage), file.path(out, "price_trade_coverage.csv"))
	feasibility = read_analysis_table(
		file.path(prices, "quantity_feasibility.csv"),
		colClasses = c(country_code = "character", enduse_code = "character")
	)
	fwrite(feasibility, file.path(out, "quantity_feasibility.csv"))
	# check numerical identities after writing summaries.
	stopifnot(
		max(
			abs(
				summary$change_pp - summary$within_sector_pp - summary$sector_composition_pp - summary$entry_exit_pp
			)
		) < 1e-8
	)
	inputs = c(
		file.path(
			prices,
			as.vector(outer(paste0("fred_IR", 0:4), c(".csv", ".html"), paste0))
		),
		file.path(raw, "enduse_imports.xlsx"),
		file.path(raw, "sources.json"),
		country,
		mapping,
		file.path(prices, paste0("bls_EIUIR", 0:4, ".json")),
		file.path(prices, "fred_CANTOT.csv"),
		file.path(prices, "sources.json"),
		file.path(prices, "quantity_feasibility.csv")
	)
	jsonlite::write_json(
		list(
			created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
			code = source_record(code_file),
			inputs = lapply(inputs, source_record),
			definitions = list(
				values = "Nominal US general merchandise imports at customs value; source end-use cells are dollar integers. Country workbook units are millions converted to dollars.",
				group = "The same fourteen pre-existing allies, with Canada and the remaining thirteen compared using the same published world denominator.",
				classification = "Six principal sectors and142five-digit categories. Other thirteen have141reported categories: no electric-energy10400rows at any date. No missing observations are filled. Reported categories recover every group total exactly.",
				source_reconciliation = "All1420world-product totals equal individual-country sums. Country annual precision is one dollar: bound0.5dollar plus0.001floating-point allowance. Fourteen nonallied country-year totals absent in country workbook remain missing comparisons against product-reported zeros.",
				decomposition = "Exact symmetric accounting using mean category weights and mean allied category shares. All retained world categories are positive in both years. Components depend on classification; neither is causal or fully price adjusted.",
				history = "2016-2025 only:9annual and8two-year changes. Acceleration uses nonoverlapping previous intervals, shift(horizon).",
				timing = "January2023-December2025, not seasonally adjusted. Monthly same-month comparisons, ratios of summed trailing-twelve-month flows, and fixed calendar-quarter ratios. No breakpoint selection.",
				prices = "Industrial-supplies,capital-goods,automotive,consumer-goods,andCanadaorigin price series lackOctober2025. Annual means require12months. No imputation,no11monthannualaverage,no complete annual price adjustment.",
				quantity = "Frozen Canada crude-oil and Korea passenger-car cells. EIA quantities are separate context without verified Census crosswalk; no unit values calculated. Korea end-use source has no quantity fields. See retained feasibility records."
			),
			raw_source_network_fetch = FALSE
		),
		file.path(out, "definition.json"),
		pretty = TRUE,
		auto_unbox = TRUE
	)
	invisible(out)
}

if (sys.nframe() == 0L) {
	args = commandArgs(trailingOnly = TRUE)
	mode = if (length(args) >= 1L) args[1L] else "trade"
	stopifnot(mode %in% c("trade", "detail"))
	script_argument = commandArgs()[grepl("^--file=", commandArgs())][1L]
	code_file = normalizePath(sub("^--file=", "", script_argument))
	root = if (length(args) >= 2L) normalizePath(args[2L]) else dirname(dirname(code_file))
	if (identical(mode, "trade")) trade_pipeline(root, code_file) else trade_census_detail_pipeline(root, code_file)
}
