# calculate request continuation across adopted decisions

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
suppressPackageStartupMessages(library(digest))
setDTthreads(1L)
main = function() {
	args = commandArgs(trailingOnly = TRUE)
	archive = if (length(args) >= 1L) args[1L] else "data/proposal_comparison_sandbox.zip"
	output_dir = if (length(args) >= 2L) args[2L] else "results/proposal_comparison"
	stopifnot(file.exists(archive))
	sandbox_root = tempfile("americium_comparison_")
	dir.create(sandbox_root)
	on.exit(unlink(sandbox_root, recursive = TRUE), add = TRUE)
	members = unzip(archive, list = TRUE)$Name
	stopifnot(!any(grepl("(^/|(^|/)\\.\\.(/|$))", members)))
	unzip(archive, exdir = sandbox_root)
	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
	# verify every archived file against the separate preservation manifest.
	manifest_path = sub("[.]zip$", ".sha256", archive)
	stopifnot(file.exists(manifest_path))
	lines = readLines(manifest_path, warn = FALSE)
	expected = substr(lines, 1L, 64L)
	relative = sub("^[0-9a-f]{64}  ", "", lines)
	stopifnot(
		length(lines) > 0L,
		!anyDuplicated(relative),
		!any(grepl("(^/|(^|/)\\.\\.(/|$))", relative))
	)
	actual = vapply(
		file.path(sandbox_root, relative),
		function(path) {
			digest(file = path, algo = "sha256")
		},
		character(1L)
	)
	stopifnot(all(actual == expected))
	# accepted registry pins the 18 batches; no superseded readings enter by glob.
	input_receipt = read_analysis_table(
		file.path(sandbox_root, "work/review/final_accepted_inputs.csv")
	)
	stopifnot(nrow(input_receipt) == 72L, !anyDuplicated(input_receipt$path))
	stopifnot(
		all(
			vapply(
				file.path(sandbox_root, input_receipt$path),
				function(path) {
					digest(file = path, algo = "sha256")
				},
				character(1L)
			) == input_receipt$sha256
		)
	)
	read_table = function(relative) {
		read_analysis_table(
			file.path(sandbox_root, relative),
			colClasses = "character",
			na.strings = NULL
		)
	}
	write_table = function(x, filename) {
		stopifnot(!grepl("[/\\\\]", filename))
		fwrite(x, file.path(output_dir, filename), na = "")
	}
	truth = function(x) {
		tolower(x) %in% c("true", "1")
	}
	feature_present = function(x) {
		!is.na(x) & !tolower(trimws(x)) %in% c("", "false", "0", "no", "none", "not applicable")
	}
	validation = list()
	record_check = function(name, passed, detail) {
		assign(
			"validation",
			append(
				validation,
				list(
					data.table(
						check = name,
						result = if (passed) "pass" else "fail",
						detail = detail
					)
				)
			),
			envir = parent.env(environment())
		)
	}
	accepted_batches = read_table("work/review/accepted_batches.csv")
	content_dirs = accepted_batches[status == "accepted_ai_reconciliation", path]
	content_names = c(
		"adopted_requests",
		"request_correspondence",
		"draft_stage_comparison",
		"paragraph_accounting"
	)
	integrated = list()
	for (name in content_names) {
		paths = file.path(content_dirs, paste0(name, ".csv"))
		paths = paths[file.exists(file.path(sandbox_root, paths))]
		stopifnot(length(paths) == 18L)
		pieces = lapply(
			paths,
			function(path) {
				x = read_table(path)
				x[, working_input := path]
				x
			}
		)
		x = rbindlist(pieces, fill = TRUE, use.names = TRUE)
		integrated[[name]] = x
		write_table(x, paste0(name, ".csv"))
	}
	requests = integrated$adopted_requests
	requests[, core_request := count_status %in% c("included", "count") & category %in% c("substantive", "substantive_request")]
	write_table(requests, "adopted_requests.csv")
	record_check(
		"included_substantive_categories_accounted",
		all(
			requests[count_status %in% c("included", "count") & grepl("^substantive", category), core_request]
		),
		"Both accepted substantive category aliases enter the core denominator."
	)
	correspondence = integrated$request_correspondence
	stages = integrated$draft_stage_comparison
	record_check(
		"core_evidence_fields_present",
		all(
			requests[core_request == TRUE, nzchar(quote) & nzchar(source_id) & nzchar(paragraph) & nzchar(requested_actor) & nzchar(action)]
		),
		"Core units retain a source passage, locator, actor and requested activity."
	)
	record_check(
		"request_ids_unique",
		!anyDuplicated(requests$request_id),
		"Unit inventory including excluded and uncertain candidates; uniqueness is not semantic accuracy."
	)
	record_check(
		"correspondence_keys_unique",
		!anyDuplicated(
			correspondence[, .(pair_id, prior_request_id, later_request_ids)]
		),
		"Each explicit correspondence tuple occurs once."
	)
	prior_ids = correspondence$prior_request_id[nzchar(correspondence$prior_request_id)]
	record_check(
		"prior_request_ids_resolve",
		all(prior_ids %in% requests$request_id),
		"Every nonblank prior request ID resolves in integrated inventory."
	)
	later_ids = trimws(
		unlist(
			strsplit(
				correspondence$later_request_ids[nzchar(correspondence$later_request_ids)],
				"[|;]"
			)
		)
	)
	record_check(
		"later_request_ids_resolve",
		all(later_ids %in% requests$request_id),
		"Every destination, including splits, resolves to an inventoried request."
	)
	record_check(
		"one_prior_row_per_pair",
		!anyDuplicated(
			correspondence[nzchar(prior_request_id), .(pair_id, prior_request_id)]
		),
		"Prior requests count once even if destinations split."
	)
	for (name in c("family_pairs", "pilot_selection")) {
		write_table(
			read_table(file.path("work/content", paste0(name, ".csv"))),
			paste0(name, ".csv")
		)
	}
	pairs = read_table("work/content/family_pairs.csv")
	pairs[, `:=`(prior_year = as.integer(prior_year), later_year = as.integer(later_year))]
	record_check(
		"pair_ids_unique",
		!anyDuplicated(pairs$pair_id),
		"Each family transition appears once in the reviewed family inventory."
	)
	record_check(
		"coded_pairs_are_eligible",
		all(correspondence$pair_id %in% pairs[truth(eligible_one_to_one), pair_id]),
		"No unmatched or off-cycle absence is treated as a coded comparison."
	)
	prior_scope = merge(
		correspondence[nzchar(prior_request_id), .(pair_id, prior_request_id)],
		requests[, .(prior_request_id = request_id, actual_root = root_id)],
		by = "prior_request_id"
	)
	prior_scope = merge(prior_scope, pairs[, .(pair_id, prior_root)], by = "pair_id")
	record_check(
		"prior_units_belong_to_pair_root",
		all(prior_scope$actual_root == prior_scope$prior_root),
		"The prior request is from the earlier adopted document assigned to that transition."
	)
	later_scope = correspondence[nzchar(later_request_ids), .(request_id = trimws(unlist(strsplit(later_request_ids, "[|;]")))), by = pair_id]
	later_scope = merge(
		later_scope,
		requests[, .(request_id, actual_root = root_id)],
		by = "request_id"
	)
	later_scope = merge(later_scope, pairs[, .(pair_id, later_root)], by = "pair_id")
	record_check(
		"later_units_belong_to_pair_root",
		all(later_scope$actual_root == later_scope$later_root),
		"Destinations resolve to the assigned later adopted document, including split mappings."
	)
	expected_later = merge(
		pairs[pair_id %in% correspondence$pair_id, .(pair_id, later_root)],
		requests[core_request == TRUE, .(later_root = root_id, request_id)],
		by = "later_root",
		allow.cartesian = TRUE
	)
	record_check(
		"every_later_core_has_mapping_or_addition",
		nrow(
			fsetdiff(
				expected_later[, .(pair_id, request_id)],
				later_scope[, .(pair_id, request_id)]
			)
		) == 0L,
		"Every later core request is mapped from an earlier request or explicitly classified as an addition."
	)
	pair_summary = pairs[, .(families = uniqueN(family_id)), by = .(prior_year, later_year, status)]
	write_table(pair_summary, "family_pair_summary.csv")
	inventory = requests[, .(candidate_rows = .N), by = .(family_id, year, count_status, coding_status)]
	write_table(inventory, "request_inventory_summary.csv")
	summary = correspondence[, .(rows = .N), by = .(pair_id, family_id, count_status, outcome, coding_status)]
	summary[, interpretation := "Reconciled AI descriptive inventory; complete-population status is recorded separately in comparison_coverage; not human validation or a causal effect."]
	write_table(summary, "adopted_comparison_summary.csv")
	stage_summary = stages[, .(rows = .N), by = .(pair_id, visible_change_timing, coding_status)]
	write_table(stage_summary, "draft_stage_summary.csv")
	prior_inventory = requests[, .(
		prior_request_id = request_id,
		core_request,
		category,
		prior_inventory_count_status = count_status,
		continuing_or_time_limited,
		named_institution_mandate,
		specific_resources_appointments,
		reporting_monitoring,
		deadline,
		access_participation,
		explicit_scope
	)]
	comparisons = merge(
		correspondence[nzchar(prior_request_id)],
		prior_inventory,
		by = "prior_request_id",
		all.x = TRUE,
		sort = FALSE
	)
	comparisons = merge(
		comparisons,
		pairs[, .(
			pair_id,
			prior_year,
			later_year,
			common_three_years,
			internal_decision
		)],
		by = "pair_id",
		all.x = TRUE,
		sort = FALSE
	)
	valid_core_outcomes = c(
		"continued_without_substantive_change",
		"continued_stronger_specified_requirements",
		"continued_weaker_specified_requirements",
		"activity_changed_not_orderable",
		"completed_or_expired",
		"no_corresponding_request_found",
		"unresolved"
	)
	record_check(
		"core_outcomes_use_reviewed_categories",
		all(comparisons[core_request == TRUE, outcome %in% valid_core_outcomes]),
		"Core earlier requests use the seven interpretable outcomes; additions and exclusions cannot enter their denominator."
	)
	subsets = list(
		all_core = comparisons[core_request == TRUE],
		continuing_core = comparisons[core_request == TRUE & continuing_or_time_limited %in% c("continuing", "ongoing_or_recurring", "recurrent")],
		boundary_candidates = comparisons[grepl("boundary|actor", prior_inventory_count_status)]
	)
	for (field in c(
		"named_institution_mandate",
		"specific_resources_appointments",
		"reporting_monitoring",
		"deadline",
		"access_participation",
		"explicit_scope"
	)) {
		subsets[[field]] = comparisons[core_request == TRUE & feature_present(get(field))]
	}
	family_results = rbindlist(
		lapply(
			names(subsets),
			function(label) {
				x = subsets[[label]]
				if (!nrow(x)) {
					return(NULL)
				}
				x[, .(numerator = .N, request_ids = paste(prior_request_id, collapse = "|")),
					by = .(
						pair_id,
						family_id,
						prior_year,
						later_year,
						common_three_years,
						internal_decision,
						outcome
					)
				][, subset := label]
			}
		),
		fill = TRUE
	)
	family_results[, denominator := sum(numerator), by = .(pair_id, subset)]
	family_results[, fraction := numerator / denominator]
	family_results[, interpretation := "Reconciled AI descriptive inventory; boundary directives separate; clauses not independent observations of policy change."]
	write_table(family_results, "family_outcome_fractions.csv")
	subset_grid = unique(comparisons[, .(pair_id, family_id)])[, .(subset = names(subsets)), by = .(pair_id, family_id)]
	subset_grid = merge(
		subset_grid,
		unique(family_results[, .(pair_id, subset, denominator)]),
		by = c("pair_id", "subset"),
		all.x = TRUE
	)
	subset_grid[is.na(denominator), denominator := 0L]
	subset_grid[, fraction_defined := denominator > 0L]
	write_table(subset_grid, "subset_denominators.csv")
	exclusions = comparisons[, .(
		inventory_rows = .N,
		core_rows = sum(core_request),
		boundary_rows = sum(grepl("boundary|actor", prior_inventory_count_status)),
		other_excluded_rows = sum(!core_request & !grepl("boundary|actor", prior_inventory_count_status))
	),
	by = .(pair_id, family_id)
	]
	write_table(exclusions, "denominator_accounting.csv")

	all_outcomes = sort(unique(comparisons$outcome))
	family_denominators = unique(
		family_results[, .(
			pair_id,
			family_id,
			prior_year,
			later_year,
			common_three_years,
			internal_decision,
			subset,
			denominator
		)]
	)
	full_grid = family_denominators[, .(outcome = all_outcomes),
		by = .(
			pair_id,
			family_id,
			prior_year,
			later_year,
			common_three_years,
			internal_decision,
			subset,
			denominator
		)
	]
	full_grid = merge(
		full_grid,
		family_results[, .(pair_id, subset, outcome, numerator)],
		by = c("pair_id", "subset", "outcome"),
		all.x = TRUE
	)
	full_grid[is.na(numerator), numerator := 0L]
	full_grid[, fraction := numerator / denominator]
	weighted = rbindlist(
		lapply(
			c("coded_pairs", "common_three_year_families"),
			function(population) {
				x = full_grid[!truth(internal_decision)]
				if (population == "common_three_year_families") x = x[truth(common_three_years)]
				x[, .(
					families = uniqueN(family_id),
					mean_family_fraction = mean(fraction),
					numerator_sum = sum(numerator),
					denominator_sum = sum(denominator)
				),
				by = .(prior_year, later_year, subset, outcome)
				][, population := population]
			}
		)
	)
	write_table(weighted, "equal_family_summary.csv")
	headline = dcast(
		weighted[subset == "all_core" & outcome != "excluded_from_substantive_comparison"],
		population + prior_year + later_year + families + denominator_sum ~ outcome,
		value.var = "mean_family_fraction"
	)
	headline[, interpretation := "Equal-family fractions for substantive requests with identified or inherited actors; denominator_sum is pooled prior requests, not the averaging weight."]
	write_table(headline, "headline_comparison.csv")
	leave_one = full_grid[!truth(internal_decision) & subset == "all_core",
		{
			lapply(
				unique(family_id),
				function(omitted) {
					keep = family_id != omitted
					data.table(
						omitted_family = omitted,
						remaining_families = uniqueN(family_id[keep]),
						mean_family_fraction = if (any(keep)) mean(fraction[keep]) else NA_real_
					)
				}
			) |> rbindlist()
		},
		by = .(prior_year, later_year, outcome)
	]
	write_table(leave_one, "leave_one_family_out.csv")
	write_table(
		comparisons[outcome %in% c("completed_or_expired", "unresolved", "no_corresponding_request_found")],
		"completion_omission_uncertainty.csv"
	)
	write_table(correspondence[!nzchar(prior_request_id)], "additions.csv")
	coverage = pairs[truth(eligible_one_to_one), .(
		pair_id,
		family_id,
		prior_year,
		later_year,
		internal_decision
	)]
	coverage[, reconciled_rows_present := pair_id %in% comparisons$pair_id]
	write_table(coverage, "comparison_coverage.csv")
	if (TRUE) {
		record_check(
			"all_eligible_pairs_reconciled",
			all(coverage$reconciled_rows_present),
			paste(
				sum(coverage$reconciled_rows_present),
				"of",
				nrow(coverage),
				"eligible pairs have reconciled rows."
			)
		)
	}
	expected_keys = pairs[
		truth(eligible_one_to_one) & pair_id %in% comparisons$pair_id,
		.(pair_id, prior_root)
	]
	expected_keys = merge(
		expected_keys,
		requests[, .(prior_root = root_id, prior_request_id = request_id)],
		by = "prior_root",
		allow.cartesian = TRUE
	)
	observed_keys = correspondence[nzchar(prior_request_id), .(pair_id, prior_request_id)]
	record_check(
		"every_prior_candidate_accounted",
		nrow(fsetdiff(expected_keys[, .(pair_id, prior_request_id)], observed_keys)) == 0L,
		"All inventoried prior candidates, including excluded ones, have a status in each coded eligible pair."
	)
	record_check(
		"draft_stage_keys_unique",
		!anyDuplicated(stages[nzchar(prior_request_id), .(pair_id, prior_request_id)]),
		"Each prior request has one three-stage comparison row."
	)
	record_check(
		"every_prior_has_draft_stage",
		nrow(fsetdiff(observed_keys, stages[, .(pair_id, prior_request_id)])) == 0L,
		"Draft-stage table accounts for every prior correspondence, including explicit unavailable stages."
	)
	record_check(
		"core_temporal_status_explicit",
		all(
			requests[core_request == TRUE, continuing_or_time_limited %in%
				c(
					"continuing",
					"ongoing_or_recurring",
					"recurrent",
					"time_limited",
					"time_limited_task",
					"one_time",
					"unresolved"
				)]
		),
		"No blank or placeholder temporal classes silently excluded from continuing-only summary."
	)
	record_check(
		"family_fraction_accounting",
		all(
			family_results[, abs(sum(fraction) - 1) < 1e-10,
				by = .(pair_id, subset)
			]$V1
		),
		"Outcome counts exhaust each defined subset; zero-core pairs have no invented fraction."
	)
	pilot = read_table("work/content/pilot_selection.csv")
	record_check(
		"six_pilot_families_present",
		all(pilot$family_id %in% requests$family_id),
		"All selected families must be represented; this does not certify complete text or independent review."
	)
	for (name in c(
		"process_coverage",
		"coverage_summary",
		"earlier_demand_comparability",
		"statement_inventory"
	)) {
		write_table(
			read_table(file.path("work/coverage", paste0(name, ".csv"))),
			paste0(name, ".csv")
		)
	}
	for (name in c("process_events", "actor_actions")) {
		write_table(
			read_table(file.path("work/process", paste0(name, ".csv"))),
			paste0(name, ".csv")
		)
	}
	events = read_table("work/process/process_events.csv")
	statement_audit = read_table("work/process/general_statement_classification_addendum.csv")
	write_table(statement_audit, "general_statement_classification_addendum.csv")
	events[, broader_explicit_institutional_direction := explicit_current_text_demand]
	events[, classification_addendum := ""]
	for (i in seq_len(nrow(statement_audit))) {
		row = statement_audit[i]
		events[event_id == row$event_id, `:=`(
			broader_explicit_institutional_direction = row$broader_explicit_institutional_direction,
			classification_addendum = row$rationale,
			exact_quote = gsub("\\\\n", "\n", row$exact_quote)
		)]
	}
	write_table(events, "process_events.csv")
	record_check(
		"process_event_ids_unique",
		!anyDuplicated(events$event_id),
		"Stable records distinguish shared actions and multiple references."
	)
	events[, us_demand := truth(US_involved) & truth(explicit_current_text_demand)]
	event_summary = events[, .(
		records = .N,
		distinct_actions = uniqueN(action_id),
		root_count = uniqueN(root_id[nzchar(root_id)])
	), by = .(year, scope, action_type, us_demand, relative_to_adoption, outcome)]
	write_table(event_summary, "process_summary.csv")
	write_table(
		read_table("work/process/process_summary.csv"),
		"process_audit_summary.csv"
	)
	manifest_files = list.files(
		file.path(sandbox_root, "work"),
		pattern = "^source_entries[.]csv$",
		recursive = TRUE,
		full.names = TRUE
	)
	manifest_paths = c(
		substring(manifest_files, nchar(sandbox_root) + 2L),
		"work/review/additional_sources.csv"
	)
	manifest = rbindlist(
		lapply(
			manifest_paths,
			function(path) {
				x = read_table(path)
				x[, work_package := dirname(path)]
				x
			}
		),
		fill = TRUE,
		use.names = TRUE
	)
	manifest[, semantic_review_note := ifelse(
		work_package == "work/content" & (stage == "adopted" | as.integer(year) > 2023L),
		"Whole main texts read; source-specific footnotes and incorporated scope qualified in reconciliation records.",
		""
	)]
	coalesce_text = function(x, fields) {
		answer = rep("", nrow(x))
		for (field in intersect(fields, names(x))) {
			take = !nzchar(answer) & !is.na(x[[field]]) & nzchar(x[[field]])
			answer[take] = x[[field]][take]
		}
		answer
	}
	manifest[, local_path := coalesce_text(manifest, c("local_path", "path"))]
	manifest[, archive_member := coalesce_text(manifest, c("archive_member", "archive_path"))]
	manifest[, original_url := coalesce_text(manifest, c("original_url", "url"))]
	manifest[, manifest_entry := seq_len(.N)]
	record_check(
		"request_source_ids_resolve",
		all(requests$source_id %in% manifest$source_id),
		"All request inventory source IDs have a provenance entry; passage accuracy rests on saved source readings."
	)
	fwrite(manifest, file.path(output_dir, "source_manifest.csv"), na = "")
	variants = manifest[nzchar(source_id), .(
		distinct_hashes = uniqueN(sha256[nzchar(sha256)]),
		work_packages = paste(unique(work_package), collapse = "|")
	), by = source_id]
	write_table(variants[distinct_hashes > 1L], "source_identity_variants.csv")
	record_check(
		"coverage_root_count",
		uniqueN(read_table("work/coverage/process_coverage.csv")$root_id) == 128L,
		"Source availability inventory covers 128 roots, not 128 verified complete negotiations."
	)
	record_check(
		"new_statement_not_multiplied_into_roots",
		all(events[grepl("^NEW_2025_US_LONGFORM_", event_id), !nzchar(root_id)]),
		"Shared statement instructions remain separate from root-specific applicability."
	)
	links = read_table("work/review/process_request_links.csv")
	write_table(links, "process_request_links.csv")
	record_check(
		"selected_process_links_resolve",
		all(links$event_id %in% events$event_id) & all(links$request_id %in% requests$request_id),
		"Selected ICT process-to-request links; not a complete new census of actor-response chains."
	)
	# these are reviewed qualitative assessments, not estimates inferred by code.
	for (name in c("claim_evidence", "explanation_assessment")) {
		write_table(
			read_table(file.path("work/review", paste0(name, ".csv"))),
			paste0(name, ".csv")
		)
	}
	checks = rbindlist(validation)
	fwrite(checks, file.path(output_dir, "validation.csv"))
	stopifnot(all(checks$result == "pass"))
	# compare every generated analytical table with the frozen accepted table.
	# character reads compare identifiers and labels as well as numeric output.
	table_files = list.files(file.path(sandbox_root, "tables"), pattern = "[.]csv$")
	reproduction = rbindlist(
		lapply(
			table_files,
			function(name) {
				generated = file.path(output_dir, name)
				accepted_file = file.path(sandbox_root, "tables", name)
				stopifnot(file.exists(generated))
				a = read_analysis_table(generated, colClasses = "character", na.strings = NULL)
				b = read_analysis_table(
					accepted_file,
					colClasses = "character",
					na.strings = NULL
				)
				# identical values and row order; ignore data.table's internal pointer.
				same = identical(names(a), names(b)) && identical(as.list(a), as.list(b))
				data.table(table = name, rows = nrow(a), matches_accepted = same)
			}
		)
	)
	fwrite(reproduction, file.path(output_dir, "accepted_table_reproduction.csv"))
	stopifnot(all(reproduction$matches_accepted))
	out_dir = output_dir
	family = read_analysis_table(file.path(output_dir, "family_outcome_fractions.csv"))
	outcome_labels = c(
		continued_without_substantive_change = "Continued without substantive change",
		continued_stronger_specified_requirements = "Continued with stronger requirements",
		continued_weaker_specified_requirements = "Continued with weaker requirements",
		activity_changed_not_orderable = "Changed, not comparable in strength",
		completed_or_expired = "Completed or expired",
		no_corresponding_request_found = "No corresponding request",
		unresolved = "Unresolved"
	)
	family = family[internal_decision == FALSE & outcome %in% names(outcome_labels)]
	# mean of family fractions, zero counts included
	summarize_families = function(x) {
		keys = c("prior_year", "later_year", "subset")
		units = unique(x[, c(keys, "family_id", "denominator"), with = FALSE])
		totals = units[, .(families = .N, denominator_sum = sum(denominator)), by = keys]
		sums = x[, .(fraction_sum = sum(numerator / denominator), numerator_sum = sum(numerator)), by = c(keys, "outcome")]
		grid = totals[, .(outcome = names(outcome_labels)), by = c(keys, "families", "denominator_sum")]
		out = merge(grid, sums, by = c(keys, "outcome"), all.x = TRUE)
		out[is.na(fraction_sum), `:=`(fraction_sum = 0, numerator_sum = 0)]
		out[, mean_family_fraction := fraction_sum / families]
		out[, fraction_sum := NULL]
		out[]
	}
	populations = list(
		coded_pairs = family,
		common_three_year_families = family[common_three_years == TRUE]
	)
	summary_all = rbindlist(
		lapply(
			names(populations),
			function(name) {
				cbind(population = name, summarize_families(populations[[name]]))
			}
		)
	)
	summary_all[, outcome_label := outcome_labels[outcome]]
	setorder(summary_all, population, prior_year, later_year, subset, outcome)
	# full equal_family_summary.csv was already reconstructed above.
	# headline rows used in the paper
	headline = summary_all[subset == "all_core" & (
		(population == "common_three_year_families") |
			(population == "coded_pairs" & prior_year == 2023 & later_year == 2025)
	)]
	headline[, percent := round(100 * mean_family_fraction, 1)]
	fwrite(headline, file.path(out_dir, "headline.csv"))
	# unchanged share by request feature, common families, annual comparisons
	unchanged = "continued_without_substantive_change"
	subsets = summary_all[
		population == "common_three_year_families" & outcome == unchanged &
			later_year - prior_year == 1 & subset != "boundary_candidates"
	]
	subsets[, percent := round(100 * mean_family_fraction, 2)]
	fwrite(subsets, file.path(out_dir, "subsets.csv"))
	# unchanged share with each family left out in turn
	core = family[subset == "all_core"]
	comparisons = list(
		core[common_three_years == TRUE & prior_year == 2023 & later_year == 2024],
		core[common_three_years == TRUE & prior_year == 2024 & later_year == 2025],
		core[prior_year == 2023 & later_year == 2025]
	)
	leave_one = rbindlist(
		lapply(
			comparisons,
			function(x) {
				rbindlist(
					lapply(
						unique(x$family_id),
						function(omitted) {
							kept = summarize_families(x[family_id != omitted])[outcome == unchanged]
							kept[, .(
								prior_year,
								later_year,
								omitted_family = omitted,
								remaining_families = families,
								mean_family_fraction
							)]
						}
					)
				)
			}
		)
	)
	fwrite(leave_one, file.path(out_dir, "unchanged_leave_one_family_out.csv"))
	leave_one_range = leave_one[, .(
		minimum_percent = round(100 * min(mean_family_fraction), 2),
		maximum_percent = round(100 * max(mean_family_fraction), 2),
		families = .N
	), by = .(prior_year, later_year)]
	fwrite(leave_one_range, file.path(out_dir, "leave_one_family_out_range.csv"))
}

main()
