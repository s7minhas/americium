# compare adopted texts and considered proposals

helper_file = if (sys.nframe() > 0L && !is.null(sys.frame(1)$ofile)) {
	sys.frame(1)$ofile
} else {
	sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
}

source(
	file.path(dirname(normalizePath(helper_file)), "00_data_helpers.R"),
	local = TRUE
)

script_file = sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
source(file.path(dirname(normalizePath(script_file)), "01_observed_votes.R"))
suppressPackageStartupMessages(library(xml2))
suppressPackageStartupMessages(library(rvest))
recent_adoptions = function(paths) {
	output_root = file.path(paths$output, "adoptions/recent")
	dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
	input_root = paths$input
	source_root = file.path(paths$input, "official")
	read_session = function(session) {
		filename = paste0("ga", session, "_resolutions.html")
		path = file.path(source_root, filename)
		manifest_path = file.path(source_root, paste0("ga", session, "_fetch_manifest.json"))
		url = paste0("https://public.e-delegate.un.org/reports/", filename)
		if (!file.exists(path) || !file.exists(manifest_path)) {
			stop("Required cached official source is missing: ", path)
		}
		manifest = jsonlite::read_json(manifest_path, simplifyVector = TRUE)
		stopifnot(
			identical(manifest$sha256, digest::digest(file = path, algo = "sha256")),
			manifest$status == 200L,
			isTRUE(manifest$accepted)
		)
		list(path = path, manifest = manifest)
	}
	sources = lapply(78:80, read_session)
	parse_session = function(index) {
		session = (78:80)[index]
		document = read_html(sources[[index]]$path)
		rows = xml_find_all(document, "//table//tr")
		out = rbindlist(
			lapply(
				rows,
				function(row) {
					cells = trimws(gsub("[[:space:]]+", " ", xml_text(xml_find_all(row, "./th|./td"))))
					if (length(cells) != 3L || !grepl("^[0-9]+/", cells[1])) {
						return(NULL)
					}
					data.table(
						number = cells[1],
						title = cells[2],
						action = cells[3],
						session = session
					)
				}
			)
		)
		out[, source_url := sources[[index]]$manifest$url]
		out[, source_fetch_utc := sources[[index]]$manifest$fetch_utc]
		out
	}
	adoptions = rbindlist(lapply(seq_along(sources), parse_session))
	adoptions[, resolution := paste0("A/RES/", number)]
	adoptions[, date_text := sub(".*adopted on ([0-9]+ [A-Za-z]+ [0-9]{4}).*", "\\1", action)]
	adoptions[, date := as.Date(date_text, "%d %B %Y")]
	adoptions[, year := as.integer(format(date, "%Y"))]
	adoptions[, method := fcase(
		grepl("without a vote", action, fixed = TRUE),
		"without_vote",
		grepl("with a vote", action, fixed = TRUE),
		"with_vote",
		default = "unclassified"
	)]
	adoptions[, normalized_title := normalize_title(title)]
	adoptions[, title_multiplicity := .N, by = c("year", "normalized_title")]
	adoptions[, unambiguous_title := title_multiplicity == 1L & nzchar(normalized_title)]
	stopifnot(
		!anyNA(adoptions$date),
		all(adoptions$method != "unclassified"),
		anyDuplicated(adoptions$resolution) == 0L
	)
	fwrite(adoptions, file.path(output_root, "official_session_adoptions.csv"))
	calendar = adoptions[year %in% 2024:2025]
	coverage = calendar[, .(adoption_entries = .N), by = c("year", "method")]
	coverage[, total_adoption_entries := sum(adoption_entries), by = "year"]
	coverage[, fraction := adoption_entries / total_adoption_entries]
	fwrite(coverage, file.path(output_root, "calendar_adoption_coverage.csv"))
	source_env = new.env(parent = emptyenv())
	load_vote_archive(
		file.path(input_root, "un_votes_long.rda"),
		envir = source_env
	)
	raw_votes = copy(source_env$un_votes_long[year %in% 2024:2025])
	recorded = unique(
		raw_votes[, .(
			year,
			resolution,
			date,
			title,
			total_yes,
			total_no,
			total_abstentions,
			total_non_voting,
			total_ms
		)]
	)
	stopifnot(anyDuplicated(recorded[, .(year, resolution)]) == 0L)
	join_audit = merge(
		recorded,
		calendar[, .(
			year,
			resolution,
			date,
			title,
			action,
			method,
			normalized_title
		)],
		by = c("year", "resolution"),
		all = TRUE,
		suffixes = c("_voting_source", "_adoption_source")
	)
	join_audit[, has_member_votes := !is.na(title_voting_source)]
	join_audit[, regular_session := !grepl("A/RES/ES-", resolution, fixed = TRUE)]
	join_audit[, regular_resolution := regular_session & grepl("A/RES/", resolution, fixed = TRUE)]
	join_audit[, title_same_after_normalization :=
		normalize_title(title_voting_source) == normalized_title]
	with_votes = join_audit[method == "with_vote" & has_member_votes]
	pattern = ".*with a vote \\(([0-9]+)-([0-9]+)-([0-9]+)\\).*"
	with_votes[, `:=`(
		official_yes = as.integer(sub(pattern, "\\1", action)),
		official_no = as.integer(sub(pattern, "\\2", action)),
		official_abstain = as.integer(sub(pattern, "\\3", action))
	)]
	with_votes[, totals_agree := total_yes == official_yes & total_no == official_no &
		total_abstentions == official_abstain]
	fwrite(join_audit, file.path(output_root, "member_vote_join_audit.csv"))
	fwrite(with_votes, file.path(output_root, "recorded_total_validation.csv"))
	fwrite(
		with_votes[totals_agree == FALSE],
		file.path(output_root, "official_total_discrepancies.csv")
	)
	fwrite(
		join_audit[has_member_votes & is.na(method)],
		file.path(output_root, "votes_outside_regular_resolution_scope.csv")
	)
	before = calendar[year == 2024 & unambiguous_title == TRUE]
	after = calendar[year == 2025 & unambiguous_title == TRUE]
	matched = merge(before, after, by = "normalized_title", suffixes = c("_2024", "_2025"))
	matched[, transition := paste(method_2024, method_2025, sep = "_to_")]
	fwrite(matched, file.path(output_root, "matched_adoption_methods.csv"))
	method_transitions = merge(
		CJ(
			method_2024 = c("with_vote", "without_vote"),
			method_2025 = c("with_vote", "without_vote")
		),
		matched[, .(matched_titles = .N), by = c("method_2024", "method_2025")],
		by = c("method_2024", "method_2025"),
		all.x = TRUE
	)
	method_transitions[is.na(matched_titles), matched_titles := 0L]
	fwrite(
		method_transitions,
		file.path(output_root, "method_transition_counts.csv")
	)
	matched_ids = matched$resolution_2025
	calendar[, match_status := fcase(
		!unambiguous_title,
		"ambiguous_within_year",
		year == 2025 & resolution %in% matched_ids,
		"matched",
		year == 2024 & resolution %in% matched$resolution_2024,
		"matched",
		default = "no_unique_exact_title_match"
	)]
	fwrite(calendar, file.path(output_root, "all_adoption_title_match_audit.csv"))
	fwrite(
		calendar[, .(adoption_entries = .N), by = c("year", "method", "match_status")],
		file.path(output_root, "adoption_title_match_support.csv")
	)
	# retain differences between the two official title representations
	inputs = match_recorded_titles(input_root)
	existing = inputs$matches[horizon == 1L & target == 2025]
	existing_audit = merge(
		existing[, .(normalized_title, resolution_baseline, resolution_target)],
		matched[, .(
			resolution_baseline = resolution_2024,
			resolution_target = resolution_2025,
			adoption_normalized_title = normalized_title
		)],
		by = c("resolution_baseline", "resolution_target"),
		all = TRUE
	)
	existing_audit[, `:=`(
		matched_in_voting_source = !is.na(normalized_title),
		matched_in_adoption_source = !is.na(adoption_normalized_title)
	)]
	fwrite(
		existing_audit,
		file.path(output_root, "matching_representation_audit.csv")
	)
	cross_source = copy(recorded[year == 2025])
	cross_source[, normalized_title := normalize_title(title)]
	cross_source[, title_multiplicity := .N, by = "normalized_title"]
	cross_source = merge(
		before,
		cross_source[title_multiplicity == 1L],
		by = "normalized_title",
		suffixes = c("_adoption_2024", "_voting_2025")
	)
	fwrite(
		cross_source,
		file.path(output_root, "cross_source_2024_adoption_2025_vote_matches.csv")
	)
	# use the archived voting title uniformly wherever a recorded vote exists
	harmonized = merge(
		calendar,
		recorded[, .(year, resolution, archived_voting_title = title)],
		by = c("year", "resolution"),
		all.x = TRUE
	)
	harmonized[, native_title := title]
	harmonized[, title_source := fifelse(
		is.na(archived_voting_title),
		"session_adoption_list",
		"archived_un_voting_dataset"
	)]
	harmonized[, title := fifelse(is.na(archived_voting_title), native_title, archived_voting_title)]
	harmonized[, normalized_title := normalize_title(title)]
	harmonized[, title_multiplicity := .N, by = c("year", "normalized_title")]
	harmonized[, unambiguous_title := title_multiplicity == 1L & nzchar(normalized_title)]
	harmonized_matched = merge(
		harmonized[year == 2024 & unambiguous_title == TRUE],
		harmonized[year == 2025 & unambiguous_title == TRUE],
		by = "normalized_title",
		suffixes = c("_2024", "_2025")
	)
	fwrite(
		harmonized,
		file.path(output_root, "harmonized_adoption_title_audit.csv")
	)
	fwrite(
		harmonized_matched,
		file.path(output_root, "harmonized_matched_adoption_methods.csv")
	)
	harmonized_counts = merge(
		CJ(
			method_2024 = c("with_vote", "without_vote"),
			method_2025 = c("with_vote", "without_vote")
		),
		harmonized_matched[, .(matched_titles = .N), by = c("method_2024", "method_2025")],
		by = c("method_2024", "method_2025"),
		all.x = TRUE
	)
	harmonized_counts[is.na(matched_titles), matched_titles := 0L]
	fwrite(
		harmonized_counts,
		file.path(output_root, "harmonized_method_transition_counts.csv")
	)
	harmonized_new = harmonized_matched[method_2024 == "without_vote" & method_2025 == "with_vote"]
	harmonized_old = harmonized_matched[method_2024 == "with_vote" & method_2025 == "with_vote"]
	harmonized_vote_support = merge(
		existing[, .(normalized_title, resolution_baseline, resolution_target)],
		harmonized_matched[, .(
			resolution_baseline = resolution_2024,
			resolution_target = resolution_2025,
			harmonized_normalized_title = normalized_title
		)],
		by = c("resolution_baseline", "resolution_target"),
		all = TRUE
	)
	fwrite(
		harmonized_vote_support,
		file.path(output_root, "harmonized_voting_match_support.csv")
	)
	# summarize actual 2025 member votes without imputing consensus votes
	newly_recorded = matched[method_2024 == "without_vote" & method_2025 == "with_vote"]
	old_recorded = matched[method_2024 == "with_vote" & method_2025 == "with_vote"]
	groups = list(
		consensus_to_recorded = newly_recorded$resolution_2025,
		recorded_to_recorded = old_recorded$resolution_2025,
		harmonized_consensus_to_recorded = harmonized_new$resolution_2025,
		harmonized_recorded_to_recorded = harmonized_old$resolution_2025,
		all_2025_regular_recorded = calendar[year == 2025 & method == "with_vote", resolution]
	)
	profile_rows = agreement_rows = list()
	for (group_name in names(groups)) {
		ids = groups[[group_name]]
		v = raw_votes[year == 2025 & resolution %in% ids]
		for (actor_group in c("USA", "allies", "all_other_members")) {
			selected = switch(
	actor_group,
				USA = v[voting_country_code == "USA"],
				allies = v[voting_country_code %in% actor_groups$allies],
				all_other_members = v[!voting_country_code %in% c("USA", actor_groups$allies)]
			)
			profile_rows[[length(profile_rows) + 1L]] = data.table(
				resolution_group = group_name,
				actor_group = actor_group,
				resolutions = length(ids),
				records = nrow(selected),
				yes = sum(selected$recorded_vote == "Y"),
				no = sum(selected$recorded_vote == "N"),
				abstain = sum(selected$recorded_vote == "A"),
				explicit_nonvote = sum(selected$recorded_vote == "X")
			)
		}
		actors = c("USA", actor_groups$allies)
		positions = matrix(NA_real_, length(actors), length(ids), dimnames = list(actors, ids))
		s = v[voting_country_code %in% actors]
		positions[cbind(match(s$voting_country_code, actors), match(s$resolution, ids))] =
			unname(c(Y = 1, A = 0, N = -1)[s$recorded_vote])
		pair_groups = list(
			us_allies = data.table(a = "USA", b = actor_groups$allies),
			within_allies = as.data.table(t(combn(actor_groups$allies, 2)))
		)
		setnames(pair_groups$within_allies, c("a", "b"))
		for (pair_group in names(pair_groups)) {
			p = pair_groups[[pair_group]]
			agreement = 1 - abs(positions[p$a, , drop = FALSE] - positions[p$b, , drop = FALSE]) / 2
			pair_mean = rowMeans(agreement, na.rm = TRUE)
			agreement_rows[[length(agreement_rows) + 1L]] = data.table(
				resolution_group = group_name,
				pair_group = pair_group,
				resolutions = length(ids),
				dyads = sum(is.finite(pair_mean)),
				mean_2025_recorded_agreement = mean(pair_mean[is.finite(pair_mean)]),
				minimum_recorded_votes_per_dyad = min(rowSums(is.finite(agreement)))
			)
		}
	}
	new_vote_sets = lapply(
		list(native = newly_recorded, harmonized = harmonized_new),
		function(d) {
			merge(
				d[, .(
					normalized_title,
					resolution_2024,
					resolution = resolution_2025,
					title_2024,
					title_2025,
					action_2024,
					action_2025
				)],
				raw_votes[year == 2025],
				by = "resolution"
			)
		}
	)
	for (representation in names(new_vote_sets)) {
		new_votes = new_vote_sets[[representation]]
		prefix = if (representation == "native") "" else "harmonized_"
		fwrite(
			new_votes,
			file.path(output_root, paste0(prefix, "newly_recorded_2025_member_votes.csv"))
		)
		new_summary = new_votes[, .(
			normalized_title = normalized_title[1],
			resolution_2024 = resolution_2024[1],
			title_2025 = title_2025[1],
			action_2024 = action_2024[1],
			action_2025 = action_2025[1],
			us_vote = recorded_vote[voting_country_code == "USA"],
			member_yes = sum(recorded_vote == "Y"),
			member_no = sum(recorded_vote == "N"),
			member_abstain = sum(recorded_vote == "A"),
			member_nonvote = sum(recorded_vote == "X"),
			nonus_yes = sum(voting_country_code != "USA" & recorded_vote == "Y"),
			nonus_no = sum(voting_country_code != "USA" & recorded_vote == "N"),
			nonus_abstain = sum(voting_country_code != "USA" & recorded_vote == "A"),
			ally_yes = sum(voting_country_code %in% actor_groups$allies & recorded_vote == "Y"),
			ally_no = sum(voting_country_code %in% actor_groups$allies & recorded_vote == "N"),
			ally_abstain = sum(voting_country_code %in% actor_groups$allies & recorded_vote == "A"),
			cached_total_yes = total_yes[1],
			cached_total_no = total_no[1],
			cached_total_abstentions = total_abstentions[1]
		), by = "resolution"]
		fwrite(
			new_summary,
			file.path(output_root, paste0(prefix, "newly_recorded_title_summary.csv"))
		)
	}
	fwrite(
		rbindlist(profile_rows),
		file.path(output_root, "recorded_vote_profiles.csv")
	)
	fwrite(
		rbindlist(agreement_rows),
		file.path(output_root, "recorded_group_agreement.csv")
	)
	# isolate the one decision included in the archived voting dataset
	decision = raw_votes[resolution == "A/DEC/80/506"]
	fwrite(decision, file.path(output_root, "decision_80_506_member_votes.csv"))
	annual_rows = list()
	annual_samples = list(
		all_2024 = raw_votes[year == 2024],
		all_2025 = raw_votes[year == 2025],
		resolutions_only_2025 = raw_votes[year == 2025 & resolution != "A/DEC/80/506"]
	)
	annual_actors = sort(unique(raw_votes$voting_country_code))
	annual_pairs = list(
		us_allies = data.table(a = "USA", b = actor_groups$allies),
		within_allies = as.data.table(t(combn(actor_groups$allies, 2))),
		all_nonus = as.data.table(t(combn(setdiff(annual_actors, "USA"), 2)))
	)
	setnames(annual_pairs$within_allies, c("a", "b"))
	setnames(annual_pairs$all_nonus, c("a", "b"))
	for (sample_name in names(annual_samples)) {
		v = annual_samples[[sample_name]]
		ids = sort(unique(v$resolution))
		positions = matrix(
			NA_real_,
			length(annual_actors),
			length(ids),
			dimnames = list(annual_actors, ids)
		)
		positions[cbind(match(v$voting_country_code, annual_actors), match(v$resolution, ids))] =
			unname(c(Y = 1, A = 0, N = -1)[v$recorded_vote])
		for (pair_group in names(annual_pairs)) {
			p = annual_pairs[[pair_group]]
			agreement = 1 - abs(positions[p$a, , drop = FALSE] - positions[p$b, , drop = FALSE]) / 2
			pair_mean = rowMeans(agreement, na.rm = TRUE)
			annual_rows[[length(annual_rows) + 1L]] = data.table(
				sample = sample_name,
				pair_group = pair_group,
				recorded_items = length(ids),
				dyads = sum(is.finite(pair_mean)),
				agreement = mean(pair_mean[is.finite(pair_mean)])
			)
		}
	}
	annual_sensitivity = rbindlist(annual_rows)
	annual_sensitivity = merge(
		annual_sensitivity,
		annual_sensitivity[sample == "all_2024", .(pair_group, baseline_agreement = agreement)],
		by = "pair_group"
	)
	annual_sensitivity[, change_from_2024 := agreement - baseline_agreement]
	fwrite(
		annual_sensitivity,
		file.path(output_root, "decision_exclusion_sensitivity.csv")
	)
	# check the recorded subset after full-adoption title ambiguity is applied
	matched_positions = function(year_value, resolution_ids) {
		out = matrix(
			NA_real_,
			length(annual_actors),
			length(resolution_ids),
			dimnames = list(annual_actors, resolution_ids)
		)
		d = raw_votes[year == year_value & resolution %in% resolution_ids]
		out[cbind(
			match(d$voting_country_code, annual_actors),
			match(d$resolution, resolution_ids)
		)] =
			unname(c(Y = 1, A = 0, N = -1)[d$recorded_vote])
		out
	}
	before_positions = matched_positions(2024, harmonized_old$resolution_2024)
	after_positions = matched_positions(2025, harmonized_old$resolution_2025)
	harmonized_agreement_rows = lapply(
		names(annual_pairs),
		function(pair_group) {
			p = annual_pairs[[pair_group]]
			before_agreement = 1 - abs(
				before_positions[p$a, , drop = FALSE] - before_positions[p$b, , drop = FALSE]
			) / 2
			after_agreement = 1 - abs(
				after_positions[p$a, , drop = FALSE] - after_positions[p$b, , drop = FALSE]
			) / 2
			shared = is.finite(before_agreement) & is.finite(after_agreement)
			before_agreement[!shared] = after_agreement[!shared] = NA_real_
			before_mean = rowMeans(before_agreement, na.rm = TRUE)
			after_mean = rowMeans(after_agreement, na.rm = TRUE)
			keep = is.finite(before_mean) & is.finite(after_mean)
			data.table(
				pair_group = pair_group,
				matched_titles = nrow(harmonized_old),
				dyads = sum(keep),
				baseline_agreement = mean(before_mean[keep]),
				target_agreement = mean(after_mean[keep]),
				change = mean(after_mean[keep] - before_mean[keep])
			)
		}
	)
	fwrite(
		rbindlist(harmonized_agreement_rows),
		file.path(output_root, "harmonized_recorded_agreement_change.csv")
	)
	checks = data.table(
		check = c(
			"all_adoption_dates_parsed",
			"all_adoption_methods_classified",
			"unique_adoption_symbols",
			"all_regular_recorded_entries_have_member_votes",
			"all_cached_regular_resolutions_have_adoption_entries",
			"no_without_vote_entry_has_member_votes",
			"recorded_dates_match",
			"recorded_totals_parse",
			"one_to_one_matched_titles",
			"all_newly_recorded_members_are_observed_votes"
		),
		passed = c(
			!anyNA(adoptions$date),
			all(adoptions$method != "unclassified"),
			anyDuplicated(adoptions$resolution) == 0L,
			all(join_audit[method == "with_vote", has_member_votes]),
			all(!is.na(join_audit[regular_resolution & has_member_votes, method])),
			all(!join_audit[method == "without_vote", has_member_votes]),
			all(with_votes$date_voting_source == with_votes$date_adoption_source),
			all(
				is.finite(with_votes$official_yes) & is.finite(with_votes$official_no) & is.finite(with_votes$official_abstain)
			),
			anyDuplicated(matched$resolution_2024) == 0L && anyDuplicated(matched$resolution_2025) == 0L,
			all(
				vapply(
					new_vote_sets,
					function(d) all(d$recorded_vote %in% c("Y", "N", "A", "X")),
					logical(1)
				)
			)
		)
	)
	fwrite(checks, file.path(output_root, "validation.csv"))
	stopifnot(all(checks$passed))
	invisible(NULL)
}

adoption_histories = function(paths) {
	output_root = file.path(paths$output, "adoptions/history")
	dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
	input_root = paths$input
	source_root = file.path(paths$input, "official")
	canonical_symbol = function(x) {
		gsub("[[:space:]\\[\\]]", "", x, perl = TRUE)
	}
	clean_space = function(x) {
		trimws(gsub("[[:space:]]+", " ", x))
	}
	source_rows = parsed_rows = list()
	for (session in 64:80) {
		keys = c(
			paste0("ga", session),
			paste0("static_ga", session),
			paste0("www_ga", session)
		)
		chosen = NULL
		for (key in keys) {
			manifest_path = file.path(source_root, paste0(key, "_fetch_manifest.json"))
			if (!file.exists(manifest_path)) {
				next
			}
			manifest = jsonlite::read_json(manifest_path, simplifyVector = TRUE)
			if (isTRUE(manifest$accepted)) {
				chosen = list(key = key, manifest = manifest)
				break
			}
		}
		if (is.null(chosen)) {
			source_rows[[length(source_rows) + 1L]] = data.table(
				session = session,
				available = FALSE,
				url = NA_character_,
				fetch_utc = NA_character_,
				sha256 = NA_character_,
				schema = NA_character_
			)
			next
		}
		path = file.path(source_root, paste0(chosen$key, "_resolutions.html"))
		stopifnot(
			identical(
				digest::digest(file = path, algo = "sha256"),
				chosen$manifest$sha256
			)
		)
		document_text = as.character(read_html(path))
		document_text = gsub("(?i)<br[[:space:]]*/?>", "<br/> ", document_text, perl = TRUE)
		document = read_html(document_text)
		link_rows = rbindlist(
			lapply(
				xml_find_all(document, "//table//tr"),
				function(row) {
					cell = xml_find_first(row, "./td|./th")
					label = clean_space(xml_text(cell))
					href = xml_attr(xml_find_first(cell, ".//a[@href]"), "href")
					if (is.na(href) || !grepl("A/RES/", href, fixed = TRUE)) {
						return(NULL)
					}
					href = URLdecode(href)
					pattern = "A/RES/[0-9]+/[0-9]+([[:space:]]*[A-Z](-[A-Z])?)?"
					value = regmatches(href, regexpr(pattern, href))
					if (length(value) != 1L) {
						return(NULL)
					}
					data.table(label = label, document_link_symbol = value)
				}
			)
		)
		tables = html_table(document, fill = TRUE)
		table = as.data.table(tables[[which.max(vapply(tables, nrow, integer(1)))]])
		setnames(table, paste0("v", seq_len(ncol(table))))
		if (ncol(table) == 3L) {
			out = table[grepl("^[0-9]+/", v1), .(resolution_display = paste0("A/RES/", v1), title = v2, action = v3)]
			schema = "e_delegate_three_columns"
		} else if (ncol(table) == 6L) {
			out = table[grepl("A/RES/", v1, fixed = TRUE), .(resolution_display = v1, title = v6, action = v4)]
			schema = "official_session_six_columns"
		} else {
			stop("Unexpected official table schema for session ", session)
		}
		out[, c("resolution_display", "title", "action") := lapply(.SD, clean_space),
			.SDcols = c("resolution_display", "title", "action")
		]
		out[, document_link_symbol := link_rows$document_link_symbol[match(resolution_display, link_rows$label)]]
		out[, resolution_parse := resolution_display]
		out[!grepl("^A/RES/", resolution_parse), resolution_parse := document_link_symbol]
		out[, identifier_repaired_from_link := resolution_parse != resolution_display]
		stopifnot(!anyNA(out$resolution_parse))
		out[, `:=`(
			session = session,
			source_url = chosen$manifest$url,
			source_fetch_utc = chosen$manifest$fetch_utc,
			source_schema = schema
		)]
		stopifnot(nrow(out) > 100L)
		parsed_rows[[length(parsed_rows) + 1L]] = out
		source_rows[[length(source_rows) + 1L]] = data.table(
			session = session,
			available = TRUE,
			url = chosen$manifest$url,
			fetch_utc = chosen$manifest$fetch_utc,
			sha256 = chosen$manifest$sha256,
			schema = schema
		)
	}
	source_coverage = rbindlist(source_rows)
	displayed = rbindlist(parsed_rows)
	displayed[, source_row := .I]
	displayed[, method := fcase(
		grepl("silence procedure", action, ignore.case = TRUE),
		"silence_procedure",
		grepl("without (a )?vote", action, ignore.case = TRUE),
		"without_vote",
		grepl(
			"[0-9]+[[:space:]]*[-–][[:space:]]*[0-9]+[[:space:]]*[-–][[:space:]]*[0-9]+",
			action
		),
		"with_vote",
		default = "unclassified"
	)]
	date_pattern = "([0-3]?[0-9])[[:space:]]*(Jan(?:uary)?|Feb(?:ruary)?|Mar(?:ch)?|Apr(?:il)?|May|Jun(?:e)?|Jul(?:y)?|Aug(?:ust)?|Sep(?:t(?:ember)?)?|Oct(?:ober)?|Nov(?:ember)?|Dec(?:ember)?)[.]?[[:space:]]*([12][0-9]{3})"
	date_fields = regmatches(
		displayed$action,
		regexec(date_pattern, displayed$action, perl = TRUE, ignore.case = TRUE)
	)
	displayed[, date := as.Date(
		vapply(
			date_fields,
			function(x) {
				if (length(x) != 4L) {
					return(NA_character_)
				}
				month = match(tolower(substr(x[3], 1, 3)), tolower(month.abb))
				sprintf("%s-%02d-%02d", x[4], month, as.integer(x[2]))
			},
			character(1)
		)
	)]
	displayed[, year := as.integer(format(date, "%Y"))]
	displayed[, canonical_display := canonical_symbol(resolution_parse)]
	displayed[, source_issue := is.na(year) | method == "unclassified" |
		!(year == session + 1945L | year == session + 1946L)]
	fwrite(displayed, file.path(output_root, "all_displayed_source_rows.csv"))
	fwrite(
		displayed[identifier_repaired_from_link == TRUE],
		file.path(output_root, "document_link_identifier_repairs.csv")
	)
	fwrite(source_coverage, file.path(output_root, "source_coverage.csv"))
	# expand explicitly bundled lettered parts without inventing an adoption method
	expanded = rbindlist(
		lapply(
			seq_len(nrow(displayed)),
			function(index) {
				d = displayed[index]
				fields = regmatches(
					d$canonical_display,
					regexec("^(.*/[0-9]+)([A-Z])-([A-Z])$", d$canonical_display)
				)[[1]]
				if (length(fields) == 4L) {
					parts = LETTERS[seq(match(fields[3], LETTERS), match(fields[4], LETTERS))]
					out = d[rep(1L, length(parts))]
					out[, canonical_resolution := paste0(fields[2], parts)]
					out[, bundled_parts := length(parts)]
					return(out)
				}
				d[, `:=`(canonical_resolution = canonical_display, bundled_parts = 1L)]
				d
			}
		)
	)
	duplicate_audit = expanded[duplicated(canonical_resolution) | duplicated(canonical_resolution, fromLast = TRUE)]
	fwrite(
		duplicate_audit,
		file.path(output_root, "duplicate_identifier_audit.csv")
	)
	units = expanded[, .(
		resolution_display = paste(sort(unique(resolution_display)), collapse = ";"),
		native_title = title[1],
		title_candidates = paste(sort(unique(title)), collapse = ";"),
		action_candidates = paste(sort(unique(action)), collapse = " | "),
		date = if (uniqueN(date) == 1L) date[1] else as.Date(NA),
		date_candidates = paste(sort(unique(as.character(date))), collapse = ";"),
		year = if (uniqueN(year) == 1L) year[1] else NA_integer_,
		method = if (uniqueN(method) == 1L) method[1] else "conflicting_method",
		session = session[1],
		source_url = source_url[1],
		source_schema = source_schema[1],
		source_rows = .N,
		bundled_parts = max(bundled_parts),
		identifier_repaired_from_link = any(identifier_repaired_from_link),
		date_conflict = uniqueN(date) > 1L,
		title_conflict = uniqueN(title) > 1L,
		classification_conflict = uniqueN(year) != 1L | uniqueN(method) != 1L | any(source_issue)
	), by = "canonical_resolution"]
	fwrite(units, file.path(output_root, "all_adoption_units.csv"))
	calendar = units[year %in% 2010:2025]
	year_coverage = data.table(year = 2010:2025)
	year_coverage[, sessions_required := paste(year - 1946L, year - 1945L, sep = ";")]
	year_coverage[, available_sessions := vapply(
		year,
		function(y) {
			sum(source_coverage[session %in% c(y - 1946L, y - 1945L), available])
		},
		integer(1)
	)]
	year_coverage[, pages_complete := available_sessions == 2L]
	year_coverage[, source_issue_entries := vapply(
		year,
		function(y) {
			sum(units[year == y, classification_conflict]) +
				sum(
					units[is.na(year) & session %in% c(y - 1946L, y - 1945L), classification_conflict]
				)
		},
		integer(1)
	)]
	year_coverage[, classification_complete := pages_complete & source_issue_entries == 0L]
	annual = merge(
		year_coverage,
		calendar[, .(
			adoption_units = .N,
			with_vote = sum(method == "with_vote"),
			without_vote = sum(method == "without_vote"),
			silence_procedure = sum(method == "silence_procedure"),
			date_conflicts = sum(date_conflict)
		), by = "year"],
		by = "year",
		all.x = TRUE
	)
	annual[, recorded_fraction := fifelse(classification_complete, with_vote / adoption_units, NA_real_)]
	displayed_year = unique(displayed[, .(canonical_display, year, method)])
	annual = merge(
		annual,
		displayed_year[year %in% 2010:2025, .(
			displayed_identifiers = .N,
			displayed_with_vote = sum(method == "with_vote"),
			displayed_without_vote = sum(method == "without_vote"),
			displayed_silence_procedure = sum(method == "silence_procedure")
		), by = "year"],
		by = "year",
		all.x = TRUE
	)
	annual[, displayed_recorded_fraction := fifelse(
		classification_complete,
		displayed_with_vote / displayed_identifiers,
		NA_real_
	)]
	setorder(annual, year)
	annual[, recorded_fraction_change := recorded_fraction - shift(recorded_fraction)]
	fwrite(annual, file.path(output_root, "annual_adoption_methods.csv"))
	# link only existing recorded member votes and retain source disagreements
	source_env = new.env(parent = emptyenv())
	load_vote_archive(
		file.path(input_root, "un_votes_long.rda"),
		envir = source_env
	)
	raw_votes = copy(source_env$un_votes_long[year %in% 2010:2025])
	raw_votes[, canonical_resolution := canonical_symbol(resolution)]
	recorded = unique(raw_votes[, .(year, canonical_resolution, resolution, date, title)])
	stopifnot(anyDuplicated(recorded[, .(year, canonical_resolution)]) == 0L)
	join_audit = merge(
		calendar,
		recorded[, .(
			year,
			canonical_resolution,
			archived_resolution = resolution,
			archived_date = date,
			archived_title = title
		)],
		by = c("year", "canonical_resolution"),
		all = TRUE
	)
	join_audit[, has_archived_votes := !is.na(archived_resolution)]
	join_audit[, recorded_date_agrees := date == archived_date]
	join_audit[, regular_resolution := grepl("^A/RES/[0-9]+/", canonical_resolution)]
	join_audit[, composite_archived_symbol := grepl("[A-Z]-[A-Z]$", canonical_resolution)]
	fwrite(join_audit, file.path(output_root, "archived_vote_join_audit.csv"))
	fwrite(
		join_audit[has_archived_votes & (is.na(method) | method != "with_vote")],
		file.path(output_root, "archived_vote_source_disagreements.csv")
	)
	fwrite(
		join_audit[method == "with_vote" & (!has_archived_votes | recorded_date_agrees == FALSE)],
		file.path(output_root, "recorded_source_coverage_disagreements.csv")
	)
	join_coverage = join_audit[, .(
		listed_recorded_units = sum(method == "with_vote", na.rm = TRUE),
		listed_recorded_with_archive = sum(method == "with_vote" & has_archived_votes, na.rm = TRUE),
		listed_unrecorded_with_archive = sum(method %in% c("without_vote", "silence_procedure") & has_archived_votes),
		archived_regular_unlisted = sum(
			has_archived_votes & regular_resolution & !composite_archived_symbol & is.na(method)
		),
		archived_composites_unlisted = sum(has_archived_votes & composite_archived_symbol & is.na(method)),
		recorded_date_disagreements = sum(
			method == "with_vote" & has_archived_votes & recorded_date_agrees == FALSE,
			na.rm = TRUE
		)
	), by = "year"]
	fwrite(
		join_coverage,
		file.path(output_root, "annual_recorded_source_coverage.csv")
	)
	archive_annual = merge(
		annual,
		join_coverage[, .(
			year,
			listed_unrecorded_with_archive,
			archived_regular_unlisted
		)],
		by = "year"
	)
	archive_annual[, `:=`(
		archive_confirmed_with_vote = with_vote + listed_unrecorded_with_archive + archived_regular_unlisted,
		archive_confirmed_adoption_units = adoption_units + archived_regular_unlisted
	)]
	archive_annual[, archive_confirmed_recorded_fraction := archive_confirmed_with_vote / archive_confirmed_adoption_units]
	fwrite(
		archive_annual,
		file.path(output_root, "archive_confirmed_annual_sensitivity.csv")
	)
	# match each adjacent year under both unchanged title representations
	method_levels = c("with_vote", "without_vote", "silence_procedure")
	match_rows = transition_rows = vote_rows = title_rows = list()
	for (representation in c("native", "harmonized", "harmonized_archive_confirmed")) {
		d = copy(join_audit[!is.na(method)])
		d[, title_used := native_title]
		d[, title_source := "official_session_list"]
		if (grepl("^harmonized", representation)) {
			d[has_archived_votes == TRUE, `:=`(title_used = archived_title, title_source = "archived_un_voting_dataset")]
		}
		if (representation == "harmonized_archive_confirmed") {
			d[has_archived_votes == TRUE, method := "with_vote"]
		}
		d[, normalized_title := normalize_title(title_used)]
		d[, title_multiplicity := .N, by = c("year", "normalized_title")]
		d[, eligible_title := title_multiplicity == 1L & nzchar(normalized_title) & !classification_conflict & !title_conflict]
		d[, representation := representation]
		title_rows[[representation]] = d
		for (target in 2011:2025) {
			baseline = target - 1L
			complete = all(annual[year %in% c(baseline, target), classification_complete])
			before = d[year == baseline & eligible_title == TRUE]
			after = d[year == target & eligible_title == TRUE]
			m = merge(
				before,
				after,
				by = "normalized_title",
				suffixes = c("_baseline", "_target")
			)
			m[, `:=`(baseline = baseline, target = target, representation = representation)]
			stopifnot(
				anyDuplicated(m$canonical_resolution_baseline) == 0L,
				anyDuplicated(m$canonical_resolution_target) == 0L
			)
			match_rows[[length(match_rows) + 1L]] = m
			for (scope in c("without_vote_only", "including_silence")) {
				unrecorded = if (scope == "without_vote_only") "without_vote" else c("without_vote", "silence_procedure")
				new = m[method_baseline %in% unrecorded & method_target == "with_vote"]
				reverse = m[method_baseline == "with_vote" & method_target %in% unrecorded]
				baseline_unrecorded = sum(m$method_baseline %in% unrecorded)
				baseline_recorded = sum(m$method_baseline == "with_vote")
				new_votes = raw_votes[year == target & canonical_resolution %in% new$canonical_resolution_target]
				us_votes = new_votes[voting_country_code == "USA"]
				ally_votes = new_votes[voting_country_code %in% actor_groups$allies]
				transition_rows[[length(transition_rows) + 1L]] = data.table(
					representation = representation,
					scope = scope,
					baseline = baseline,
					target = target,
					calendar_coverage_complete = complete,
					matched_titles = nrow(m),
					baseline_eligible_titles = nrow(before),
					target_eligible_titles = nrow(after),
					baseline_matched_unrecorded = baseline_unrecorded,
					baseline_matched_recorded = baseline_recorded,
					unrecorded_to_recorded = if (complete) nrow(new) else NA_integer_,
					recorded_to_unrecorded = if (complete) nrow(reverse) else NA_integer_,
					new_recorded_rate = if (complete && baseline_unrecorded > 0L) nrow(new) / baseline_unrecorded else NA_real_,
					reverse_rate = if (complete && baseline_recorded > 0L) nrow(reverse) / baseline_recorded else NA_real_,
					net_recorded_switches = if (complete) nrow(new) - nrow(reverse) else NA_integer_,
					us_yes = sum(us_votes$recorded_vote == "Y"),
					us_no = sum(us_votes$recorded_vote == "N"),
					us_abstain = sum(us_votes$recorded_vote == "A"),
					us_explicit_nonvote = sum(us_votes$recorded_vote == "X"),
					us_missing_records = nrow(new) - nrow(us_votes),
					us_no_fraction_substantive = if (any(us_votes$recorded_vote %in% c("Y", "N", "A"))) {
						mean(us_votes[recorded_vote %in% c("Y", "N", "A"), recorded_vote] == "N")
					} else {
						NA_real_
					},
					ally_yes = sum(ally_votes$recorded_vote == "Y"),
					ally_no = sum(ally_votes$recorded_vote == "N"),
					ally_abstain = sum(ally_votes$recorded_vote == "A"),
					ally_explicit_nonvote = sum(ally_votes$recorded_vote == "X"),
					ally_missing_records = nrow(new) * length(actor_groups$allies) - nrow(ally_votes)
				)
				if (nrow(new_votes) > 0L) {
					profiles = new_votes[, .(
						yes = sum(recorded_vote == "Y"),
						no = sum(recorded_vote == "N"),
						abstain = sum(recorded_vote == "A"),
						nonvote = sum(recorded_vote == "X"),
						us_vote = if (any(voting_country_code == "USA")) recorded_vote[voting_country_code == "USA"] else NA_character_
					),
					by = c("canonical_resolution", "resolution")
					]
					profiles = merge(
						new[, .(
							canonical_resolution = canonical_resolution_target,
							baseline_resolution = canonical_resolution_baseline,
							normalized_title,
							baseline_title = title_used_baseline,
							target_title = title_used_target,
							baseline_method = method_baseline,
							target_method = method_target
						)],
						profiles,
						by = "canonical_resolution",
						all.x = TRUE
					)
					profiles[, `:=`(
						representation = representation,
						scope = scope,
						baseline = baseline,
						target = target
					)]
					vote_rows[[length(vote_rows) + 1L]] = profiles
				}
			}
		}
	}
	matches = rbindlist(match_rows, fill = TRUE)
	transitions = rbindlist(transition_rows)
	fwrite(
		rbindlist(title_rows, fill = TRUE),
		file.path(output_root, "title_representation_audit.csv")
	)
	fwrite(matches, file.path(output_root, "historical_title_matches.csv"))
	fwrite(
		transitions,
		file.path(output_root, "historical_method_transitions.csv")
	)
	fwrite(
		rbindlist(vote_rows, fill = TRUE),
		file.path(output_root, "newly_recorded_title_vote_profiles.csv")
	)
	fwrite(
		matches[method_baseline == "with_vote" & method_target %in% c("without_vote", "silence_procedure")],
		file.path(output_root, "historical_reverse_transition_titles.csv")
	)
	# apply the unchanged coarse issue classifier
	issue_metadata = unique(
		raw_votes[, .(year, canonical_resolution, title, agenda_title, subjects)]
	)
	stopifnot(
		anyDuplicated(issue_metadata[, .(year, canonical_resolution)]) == 0L
	)
	issue_metadata[, metadata_available := rowSums(
		vapply(
			.SD,
			function(x) {
				!is.na(x) & nzchar(trimws(x))
			},
			logical(.N)
		)
	) > 0L, .SDcols = c("title", "agenda_title", "subjects")]
	issue_metadata[, issue := classify_agenda(title, agenda_title, subjects)]
	issue_metadata[metadata_available == FALSE, issue := "Uncoded metadata"]
	issue_profiles = merge(
		rbindlist(vote_rows, fill = TRUE),
		issue_metadata,
		by.x = c("target", "canonical_resolution"),
		by.y = c("year", "canonical_resolution"),
		all.x = TRUE
	)
	issue_profiles[is.na(issue), issue := "Uncoded metadata"]
	issue_counts = issue_profiles[, .(
		titles = .N,
		us_yes = sum(us_vote == "Y", na.rm = TRUE),
		us_no = sum(us_vote == "N", na.rm = TRUE),
		us_abstain = sum(us_vote == "A", na.rm = TRUE),
		us_nonvote = sum(us_vote == "X", na.rm = TRUE),
		us_missing = sum(is.na(us_vote))
	),
	by = c("representation", "scope", "baseline", "target", "issue")
	]
	fwrite(
		issue_profiles,
		file.path(output_root, "newly_recorded_title_issue_profiles.csv")
	)
	fwrite(
		issue_counts,
		file.path(output_root, "newly_recorded_issue_counts.csv")
	)
	method_counts = matches[, .(matched_titles = .N),
		by = c("representation", "baseline", "target", "method_baseline", "method_target")
	]
	fwrite(
		method_counts,
		file.path(output_root, "historical_three_method_counts.csv")
	)
	source_sensitivity = merge(
		transitions[representation == "harmonized"],
		transitions[representation == "harmonized_archive_confirmed"],
		by = c("scope", "baseline", "target"),
		suffixes = c("_listed", "_archive_confirmed")
	)
	for (metric in c(
		"unrecorded_to_recorded",
		"recorded_to_unrecorded",
		"new_recorded_rate",
		"net_recorded_switches",
		"us_no"
	)) {
		source_sensitivity[, (paste0(metric, "_difference")) :=
			get(paste0(metric, "_archive_confirmed")) - get(paste0(metric, "_listed"))]
	}
	fwrite(
		source_sensitivity,
		file.path(output_root, "method_disagreement_sensitivity.csv")
	)
	# report the complete historical ranking without selecting comparison years
	rank_rows = list()
	for (metric in c(
		"recorded_fraction",
		"recorded_fraction_change",
		"displayed_recorded_fraction"
	)) {
		past = annual[year < 2025 & is.finite(get(metric))]
		current = annual[year == 2025, get(metric)]
		rank_rows[[length(rank_rows) + 1L]] = data.table(
			representation = "adoption_inventory",
			scope = "regular_resolution_units",
			metric = metric,
			value_2025 = current,
			earlier_comparators = nrow(past),
			rank_2025_descending = 1L + sum(past[[metric]] > current),
			previous_maximum = max(past[[metric]]),
			previous_maximum_target = past$year[which.max(past[[metric]])]
		)
	}
	for (representation_name in unique(transitions$representation)) {
		for (scope_name in unique(transitions$scope)) {
			for (metric in c(
				"unrecorded_to_recorded",
				"recorded_to_unrecorded",
				"new_recorded_rate",
				"net_recorded_switches",
				"us_no",
				"us_no_fraction_substantive"
			)) {
				d = transitions[representation == representation_name & scope == scope_name & calendar_coverage_complete]
				past = d[target < 2025 & is.finite(get(metric))]
				current = d[target == 2025, get(metric)]
				rank_rows[[length(rank_rows) + 1L]] = data.table(
					representation = representation_name,
					scope = scope_name,
					metric = metric,
					value_2025 = current,
					earlier_comparators = nrow(past),
					rank_2025_descending = 1L + sum(past[[metric]] > current),
					previous_maximum = max(past[[metric]]),
					previous_maximum_target = past$target[which.max(past[[metric]])]
				)
			}
		}
	}
	fwrite(rbindlist(rank_rows), file.path(output_root, "historical_ranks.csv"))
	# check the extension against the completed terminal audit
	old_coverage = read_analysis_table(
		file.path(
			paths$output,
			"adoptions",
			"recent",
			"calendar_adoption_coverage.csv"
		)
	)
	checks = list(
		data.table(
			check = "all17requested_session_sources_available",
			passed = all(source_coverage$available)
		)
	)
	checks[[length(checks) + 1L]] = data.table(
		check = "all_source_dates_and_methods_classified",
		passed = !any(displayed$source_issue)
	)
	checks[[length(checks) + 1L]] = data.table(
		check = "all_calendar_year_classifications_complete",
		passed = all(annual$classification_complete)
	)
	checks[[length(checks) + 1L]] = data.table(
		check = "canonical_unit_identifiers_unique",
		passed = anyDuplicated(units$canonical_resolution) == 0L
	)
	checks[[length(checks) + 1L]] = data.table(
		check = "identifier_repairs_match_document_links",
		passed = all(
			displayed[identifier_repaired_from_link == TRUE, resolution_parse == document_link_symbol]
		)
	)
	checks[[length(checks) + 1L]] = data.table(
		check = "all_listed_recorded_units_have_archived_votes",
		passed = all(
			join_coverage$listed_recorded_units == join_coverage$listed_recorded_with_archive
		)
	)
	checks[[length(checks) + 1L]] = data.table(
		check = "all_archived_individual_regular_units_are_listed",
		passed = all(join_coverage$archived_regular_unlisted == 0L)
	)
	checks[[length(checks) + 1L]] = data.table(
		check = "issue_metadata_join_preserves_all_new_title_profiles",
		passed = nrow(issue_profiles) == nrow(rbindlist(vote_rows, fill = TRUE))
	)
	for (year_value in 2024:2025) {
		for (method_name in c("with_vote", "without_vote")) {
			a = old_coverage[year == year_value & method == method_name, adoption_entries]
			b = calendar[year == year_value & method == method_name, .N]
			checks[[length(checks) + 1L]] = data.table(
				check = paste("recent_coverage", year_value, method_name, sep = "_"),
				passed = identical(a, b)
			)
		}
	}
	for (representation_name in c("native", "harmonized")) {
		old_name = if (representation_name == "native") "method_transition_counts.csv" else "harmonized_method_transition_counts.csv"
		old = read_analysis_table(file.path(paths$output, "adoptions", "recent", old_name))
		current = transitions[representation == representation_name & scope == "without_vote_only" & target == 2025]
		checks[[length(checks) + 1L]] = data.table(
			check = paste("recent_transitions", representation_name, sep = "_"),
			passed = current$unrecorded_to_recorded == old[method_2024 == "without_vote" & method_2025 == "with_vote", matched_titles] &&
				current$recorded_to_unrecorded == old[method_2024 == "with_vote" & method_2025 == "without_vote", matched_titles]
		)
	}
	validation = rbindlist(checks)
	fwrite(validation, file.path(output_root, "validation.csv"))
	stopifnot(all(validation$passed))
	invisible(NULL)
}

resolution_content = function(paths) {
	# keep every requested pair, including updates and unavailable evidence.
	review = read_analysis_table(file.path(paths$input, "resolution_review.csv"))
	output_root = file.path(paths$output, "adoptions/content")
	dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
	stopifnot(
		anyDuplicated(
			review[, c("review_set", "baseline", "target", "symbol_before", "symbol_after")]
		) == 0L,
		all(review$strict_comparable == (review$classification == "routine_updates")),
		all(
			review$classification %in% c(
				"routine_updates",
				"substantive_updates_same_question",
				"material_change_limiting_vote_interpretation",
				"insufficient_evidence"
			)
		)
	)
	counts = review[, .(reviewed_pairs = .N), by = c("review_set", "baseline", "target", "classification")]
	fwrite(counts, file.path(output_root, "review_counts.csv"))
	accounting = read_analysis_table(
		file.path(
			paths$output,
			"observed/vote_changes/individual_vote_accounting.csv"
		)
	)[
		target == 2025 & baseline %in% 2023:2024 & scoring == "Partial credit for abstention"
	]
	accounting[, original_pairs_per_ally := .N, by = c("baseline", "target", "ally")]
	accounting[, number_allies := uniqueN(ally), by = c("baseline", "target")]
	stopifnot(all(accounting$number_allies == 14L))

	accounting[, original_weight := 100 / original_pairs_per_ally / number_allies]
	match_ids = read_analysis_table(
		file.path(paths$output, "observed/vote_changes/matched_titles.csv")
	)[
		year_target == 2025 & year_baseline %in% 2023:2024,
		.(
			baseline = year_baseline,
			target = year_target,
			normalized_title,
			symbol_before = resolution_baseline,
			symbol_after = resolution_target
		)
	]
	switches = unique(accounting[u_before != u_after, .(baseline, target, normalized_title)])
	switches = merge(switches, match_ids, by = c("baseline", "target", "normalized_title"))
	coding = review[review_set == "us_switches"]
	key = c("baseline", "target", "symbol_before", "symbol_after")
	stopifnot(
		nrow(switches) == 29L,
		nrow(coding) == 29L,
		nrow(merge(switches, coding, by = key)) == 29L
	)
	coded = merge(switches, coding[, c(key, "classification"), with = FALSE], by = key)
	votes = merge(
		accounting[u_before != u_after],
		coded,
		by = c("baseline", "target", "normalized_title")
	)
	stopifnot(nrow(votes) == nrow(accounting[u_before != u_after]))
	parts = votes[, .(
		reviewed_pairs = uniqueN(normalized_title),
		paired_votes = .N,
		allied_votes_retained = sum(a_before == a_after),
		us_contribution_pp = sum(us_contribution * original_weight)
	),
	by = c("baseline", "target", "classification")
	]
	totals = accounting[, .(
		full_us_contribution_pp = sum(us_contribution * original_weight),
		full_agreement_change_pp = sum((agreement_after - agreement_before) * original_weight)
	),
	by = c("baseline", "target")
	]
	parts = merge(parts, totals, by = c("baseline", "target"))
	fwrite(parts, file.path(output_root, "us_switch_contributions.csv"))
	stopifnot(
		all(
			abs(
				parts[, sum(us_contribution_pp) - unique(full_us_contribution_pp),
					by = c("baseline", "target")
				]$V1
			) < 1e-10
		)
	)
	new_titles = read_analysis_table(
		file.path(
			paths$output,
			"adoptions/recent/harmonized_newly_recorded_title_summary.csv"
		)
	)
	new_coding = review[review_set == "newly_recorded"]
	joined = merge(
		new_titles,
		new_coding,
		by.x = c("resolution_2024", "resolution"),
		by.y = c("symbol_before", "symbol_after")
	)
	stopifnot(
		nrow(new_titles) == 33L,
		nrow(new_coding) == 33L,
		nrow(joined) == 33L
	)
	support = joined[, .(
		reviewed_pairs = .N,
		us_yes = sum(us_vote == "Y"),
		us_no = sum(us_vote == "N"),
		ally_yes = sum(ally_yes),
		ally_no = sum(ally_no),
		ally_abstain = sum(ally_abstain),
		allied_votes = sum(ally_yes + ally_no + ally_abstain)
	),
	by = "classification"
	]
	fwrite(support, file.path(output_root, "newly_recorded_support.csv"))
	invisible(NULL)
}

documentary_event_summary = function(events_path, output_root) {
	e = read_analysis_table(events_path, na.strings = c(""))
	required = c(
		"event_id",
		"family_id",
		"year",
		"source_id",
		"stage",
		"observation_or_stated_reason",
		"page_or_timestamp",
		"uncertainty"
	)
	stopifnot(
		all(required %in% names(e)),
		!anyDuplicated(e$event_id),
		all(e$year %in% 2023:2025),
		!anyNA(e[, .(family_id, year, source_id, page_or_timestamp)])
	)
	panel = unique(e[, .(family_id, year)])
	stopifnot(
		nrow(panel) == 39L,
		data.table::uniqueN(panel$family_id) == 13L,
		all(panel[, .N, by = family_id]$N == 3L)
	)
	out = e[, .(
		families = data.table::uniqueN(family_id),
		event_records = .N,
		sources = data.table::uniqueN(source_id),
		plenary_actions = sum(stage == "plenary_action")
	),
	by = year
	][order(year)]
	stopifnot(all(out$plenary_actions == 13L))
	dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
	cases = data.table::as.data.table(
		read.csv(
			unz("data/resolution_sources.zip", "documentary_panel/design/cases.csv"),
			stringsAsFactors = FALSE
		)
	)
	cases[, `:=`(
		principal = tolower(principal) == "true",
		selected_comparator = tolower(selected_comparator) == "true"
	)]
	stopifnot(
		!anyDuplicated(cases$family_id),
		setequal(cases$family_id, panel$family_id)
	)
	design = cases[, .(
		families = .N,
		principal_families = sum(principal),
		comparison_families = sum(selected_comparator),
		first = sum(committee == "first" & selected_comparator),
		fourth = sum(committee == "fourth" & selected_comparator),
		second = sum(committee == "second" & selected_comparator),
		plenary = sum(committee == "plenary" & selected_comparator),
		case_years = nrow(panel)
	)]
	stopifnot(design$principal_families == 7L, design$comparison_families == 6L)
	data.table::fwrite(design, file.path(output_root, "documentary_design.csv"))
	data.table::fwrite(out, file.path(output_root, "documentary_inventory.csv"))
	invisible(out)
}

# complete both proposal comparisons from the retained source coding.
# inputs are adjudicated records; this file never classifies source text.
proposal_complete_stage_comparisons = function(commitments, coverage, stage_inventory) {
	c = data.table::copy(data.table::as.data.table(commitments))
	v = data.table::copy(data.table::as.data.table(coverage))
	s = data.table::copy(data.table::as.data.table(stage_inventory))
	stopifnot(
		all(
			c(
				"root_id",
				"earliest_issued_version",
				"earliest_documented_considered_version"
			) %in% names(s)
		),
		!anyDuplicated(s$root_id),
		!anyNA(s$earliest_documented_considered_version)
	)
	absent_issued = is.na(s$earliest_issued_version) | s$earliest_issued_version == ""
	stopifnot(
		all(s$proposal_type[absent_issued] == "substantive_internal_decision")
	)
	s[, same_initial_version := (is.na(earliest_issued_version) | earliest_issued_version == "") | earliest_issued_version == earliest_documented_considered_version]
	same = s[same_initial_version == TRUE, root_id]
	different = s[same_initial_version == FALSE, root_id]
	# the sensitivity includes all roots. only changed initial versions require new
	# coding; identical source comparisons are inherited explicitly, not omitted.
	stopifnot(
		setequal(v[comparison_basis == "formally_considered", root_id], s$root_id),
		setequal(v[comparison_basis == "earliest_issued", root_id], different),
		!any(c[comparison_basis == "earliest_issued", root_id] %in% same)
	)
	c[, comparison_origin := "separately_coded_source_comparison"]
	v[, comparison_origin := "separately_coded_source_comparison"]
	carry_c = data.table::copy(
		c[comparison_basis == "formally_considered" & root_id %in% same]
	)
	carry_v = data.table::copy(
		v[comparison_basis == "formally_considered" & root_id %in% same]
	)
	carry_c[, `:=`(
		comparison_basis = "earliest_issued",
		comparison_origin = "same_initial_and_final_sources_as_primary"
	)]
	carry_v[, `:=`(
		comparison_basis = "earliest_issued",
		comparison_origin = "same_initial_and_final_sources_as_primary"
	)]
	internal_ids = s[absent_issued, root_id]
	carry_c[root_id %in% internal_ids, comparison_origin := "unnumbered_internal_decision; no_separate_issued_draft; source_record_as_primary"]
	carry_v[root_id %in% internal_ids, comparison_origin := "unnumbered_internal_decision; no_separate_issued_draft; source_record_as_primary"]
	list(
		commitments = data.table::rbindlist(list(c, carry_c), use.names = TRUE),
		coverage = data.table::rbindlist(list(v, carry_v), use.names = TRUE)
	)
}

proposal_content_summary = function(roots, commitments, coverage) {
	r = data.table::copy(data.table::as.data.table(roots))
	c = data.table::copy(data.table::as.data.table(commitments))
	v = data.table::copy(data.table::as.data.table(coverage))
	require_fields = function(x, fields) {
		if (!all(fields %in% names(x))) stop("Missing fields: ", paste(setdiff(fields, names(x)), collapse = ", "))
	}
	require_fields(
		r,
		c("root_id", "year", "proposal_type", "adopted_continuation")
	)
	require_fields(
		c,
		c(
			"root_id",
			"comparison_basis",
			"commitment_id",
			"material_request",
			"retention",
			"review_status",
			"request_origin",
			"count_in_denominator",
			"duplicate_of_commitment_ids"
		)
	)
	require_fields(
		v,
		c(
			"root_id",
			"comparison_basis",
			"coding_complete",
			"material_inventory_complete",
			"review_status"
		)
	)
	stopifnot(
		!anyDuplicated(r$root_id),
		!anyNA(r$adopted_continuation),
		is.logical(r$adopted_continuation),
		!anyDuplicated(c[, .(root_id, comparison_basis, commitment_id)]),
		!anyDuplicated(v[, .(root_id, comparison_basis)]),
		all(c$root_id %in% r$root_id),
		all(v$root_id %in% r$root_id),
		all(v$comparison_basis %in% c("formally_considered", "earliest_issued")),
		all(c$comparison_basis %in% c("formally_considered", "earliest_issued")),
		all(c$request_origin %in% c("initial", "addition", "challenged_preamble")),
		all(c$review_status == "adjudicated"),
		all(v$review_status == "adjudicated"),
		is.logical(v$coding_complete),
		all(v$coding_complete),
		is.logical(v$material_inventory_complete),
		!anyNA(v$material_inventory_complete),
		is.logical(c$material_request),
		is.logical(c$count_in_denominator),
		!anyNA(c$count_in_denominator),
		setequal(v[comparison_basis == "formally_considered", root_id], r$root_id)
	)
	# repeated main-text/annex duties remain visible without counting the same
	# request twice. every omitted material row must point to retained canonical
	# requests within the same root, source comparison and request origin.
	duplicates = c[material_request == TRUE & count_in_denominator == FALSE]
	if (nrow(duplicates)) {
		for (i in seq_len(nrow(duplicates))) {
			row = duplicates[i]
			ids = strsplit(row$duplicate_of_commitment_ids, "|", fixed = TRUE)[[1L]]
			stopifnot(
				length(ids) > 0L,
				all(nzchar(ids)),
				!anyNA(ids),
				!anyDuplicated(ids)
			)
			canonical = c[root_id == row$root_id & comparison_basis == row$comparison_basis &
				request_origin == row$request_origin & commitment_id %in% ids]
			stopifnot(
				nrow(canonical) == length(ids),
				all(canonical$count_in_denominator),
				all(canonical$material_request),
				!anyNA(canonical$material_request)
			)
		}
	}

	levels = c(
		"substantively_retained",
		"narrowed_weakened",
		"removed",
		"substantively_replaced",
		"unknown"
	)
	m = c[request_origin == "initial" & material_request & count_in_denominator]
	stopifnot(all(m$retention %in% levels), !anyNA(m$retention))
	key = c("root_id", "comparison_basis")
	stopifnot(
		nrow(merge(unique(c[, ..key]), v[, ..key], by = key)) == nrow(unique(c[, ..key]))
	)
	ambiguous_material = unique(c[request_origin == "initial" & is.na(material_request), ..key])
	if (nrow(ambiguous_material)) {
		check = merge(ambiguous_material, v, by = key)
		stopifnot(!any(check$material_inventory_complete))
	}
	s = m[, .(
		material_requests = .N,
		retained = sum(retention == "substantively_retained"),
		weakened = sum(retention == "narrowed_weakened"),
		removed = sum(retention == "removed"),
		replaced = sum(retention == "substantively_replaced"),
		unknown = sum(retention == "unknown")
	), by = key]
	out = merge(
		v[, c(key, "material_inventory_complete"), with = FALSE],
		s,
		by = key,
		all.x = TRUE,
		sort = FALSE
	)
	counts = c(
		"material_requests",
		"retained",
		"weakened",
		"removed",
		"replaced",
		"unknown"
	)
	for (field in counts) data.table::set(out, which(is.na(out[[field]])), field, 0L)
	out = merge(out, r, by = "root_id", all.x = TRUE, sort = FALSE)
	a = c[request_origin == "addition" & material_request & count_in_denominator,
		.(added_requests = .N),
		by = key
	]
	out = merge(out, a, by = key, all.x = TRUE, sort = FALSE)
	out[is.na(added_requests), added_requests := 0L]
	au = c[request_origin == "addition" & is.na(material_request) & count_in_denominator,
		.(unknown_added_requests = .N),
		by = key
	]
	out = merge(out, au, by = key, all.x = TRUE, sort = FALSE)
	out[is.na(unknown_added_requests), unknown_added_requests := 0L]
	out[, added_requests_upper := added_requests + unknown_added_requests]
	out[, `:=`(
		any_weakening = weakened > 0L,
		any_removal = removed > 0L,
		any_replacement = replaced > 0L,
		any_unknown = unknown > 0L | !material_inventory_complete
	)]
	out[, known_change := any_weakening | any_removal | any_replacement]
	out[, content_category := data.table::fcase(
		!adopted_continuation,
		"no_adopted_continuation",
		material_requests == 0L & material_inventory_complete,
		"no_classifiable_initial_material_requests",
		known_change,
		"at_least_one_known_substantive_change",
		any_unknown,
		"unresolved_without_known_change",
		retained == material_requests,
		"fully_retained",
		default = "invalid"
	)]
	stopifnot(
		!any(out$content_category == "invalid"),
		all(
			out$material_requests == out$retained + out$weakened + out$removed + out$replaced + out$unknown
		)
	)
	# unknown requests are not deleted or retained in a point estimate.
	out[, `:=`(retained_fraction_lower = NA_real_, retained_fraction_upper = NA_real_)]
	out[
		adopted_continuation & material_requests > 0L,
		`:=`(
			retained_fraction_lower = retained / material_requests,
			retained_fraction_upper = (retained + unknown) / material_requests
		)
	]
	out[
		adopted_continuation & !material_inventory_complete,
		`:=`(retained_fraction_lower = 0, retained_fraction_upper = 1)
	]
	out[, fully_retained_lower := as.integer(content_category == "fully_retained")]
	out[, fully_retained_upper := as.integer(
		content_category %in% c("fully_retained", "unresolved_without_known_change")
	)]
	data.table::setorder(out, comparison_basis, year, root_id)
	out[]
}

proposal_content_tables = function(root_summary) {
	x = data.table::copy(root_summary)
	pooled = data.table::copy(x)
	pooled[, proposal_type := "all_eligible_types"]
	x = data.table::rbindlist(list(x, pooled), use.names = TRUE)
	counts = x[, .(roots = .N), by = .(comparison_basis, year, proposal_type, content_category)]
	totals = x[, .(
		considered_roots = .N,
		adopted_paths = sum(adopted_continuation),
		roots_with_initial_material_requests = sum(adopted_continuation & material_requests > 0L),
		roots_with_incomplete_material_inventory = sum(adopted_continuation & !material_inventory_complete),
		roots_with_weakening = sum(adopted_continuation & any_weakening),
		roots_with_removal = sum(adopted_continuation & any_removal),
		roots_with_replacement = sum(adopted_continuation & any_replacement),
		roots_with_unknown = sum(adopted_continuation & any_unknown),
		roots_with_additions = sum(adopted_continuation & added_requests > 0L),
		roots_with_additions_upper = sum(adopted_continuation & added_requests_upper > 0L),
		distinct_additions_lower = sum(added_requests),
		distinct_additions_upper = sum(added_requests_upper)
	),
	by = .(comparison_basis, year, proposal_type)
	]
	bounds = x[adopted_continuation & (material_requests > 0L | !material_inventory_complete),
		.(
			eligible_roots = .N,
			fully_retained_lower = mean(fully_retained_lower),
			fully_retained_upper = mean(fully_retained_upper),
			equal_root_request_retention_lower = mean(retained_fraction_lower),
			equal_root_request_retention_upper = mean(retained_fraction_upper)
		),
		by = .(comparison_basis, year, proposal_type)
	]
	bounds = merge(
		totals[, .(comparison_basis, year, proposal_type)],
		bounds,
		by = c("comparison_basis", "year", "proposal_type"),
		all.x = TRUE
	)
	bounds[is.na(eligible_roots), eligible_roots := 0L]
	list(counts = counts, totals = totals, bounds = bounds)
}

proposal_adoption_bounds = function(n, a, u, k) {
	stopifnot(
		length(n) == length(a),
		length(n) == length(u),
		length(n) == length(k),
		all(is.finite(n)),
		all(is.finite(a)),
		all(is.finite(u)),
		all(is.finite(k)),
		all(n >= 0),
		all(a >= 0),
		all(u >= 0),
		all(k >= 0),
		all(n == floor(n)),
		all(a == floor(a)),
		all(u == floor(u)),
		all(k == floor(k)),
		all(a + u <= n)
	)
	denominator = n + k
	data.table::data.table(
		n = n,
		adopted = a,
		unresolved_fates = u,
		potential_additional_roots = k,
		lower = ifelse(denominator > 0, a / denominator, NA_real_),
		upper = ifelse(denominator > 0, (a + u + k) / denominator, NA_real_)
	)
}

proposal_fate_summary = function(roots) {
	r = data.table::copy(data.table::as.data.table(roots))
	required = c(
		"root_id",
		"year",
		"proposal_type",
		"authority_route",
		"committee_status",
		"plenary_yearend_status",
		"plenary_primary_status",
		"adopted_successors",
		"plenary_final_symbol"
	)
	stopifnot(all(required %in% names(r)), !anyDuplicated(r$root_id))
	normalize_status = function(x) {
		y = data.table::fcase(
			x %in% c(
				"adopted",
				"adopted_successor",
				"committee_adopted",
				"adopted_internal_decision"
			),
			"adopted",
			x %in% c(
				"rejected",
				"withdrawn_without_successor",
				"pending",
				"deferred",
				"unknown",
				"not_applicable"
			),
			x,
			default = "invalid"
		)
		if (any(y == "invalid")) {
			stop(
				"Unrecognized disposition: ",
				paste(unique(x[y == "invalid"]), collapse = ", ")
			)
		}
		y
	}
	pieces = lapply(
		c("committee_yearend", "plenary_yearend", "plenary_primary"),
		function(stage) {
			field = switch(
	stage,
				committee_yearend = "committee_status",
				plenary_yearend = "plenary_yearend_status",
				plenary_primary = "plenary_primary_status"
			)
			z = r[, .(root_id, year, proposal_type, authority_route)]
			z[, `:=`(
				stage = stage,
				status = normalize_status(r[[field]]),
				successors = if (stage == "committee_yearend") r$adopted_successors else r$plenary_final_symbol
			)]
			z
		}
	)
	paths = data.table::rbindlist(pieces)
	internal = grepl("internal", r$authority_route, fixed = TRUE)
	ids = r$root_id[internal]
	paths[root_id %in% ids & stage == "committee_yearend" & status == "adopted" &
		(is.na(successors) | successors == ""), successors := paste0("internal:", root_id)]
	stopifnot(
		!any(paths[status == "adopted", is.na(successors) | successors == ""])
	)
	primary = data.table::copy(paths[stage == "plenary_primary" & !root_id %in% ids])
	primary = data.table::rbindlist(
		list(primary, paths[stage == "committee_yearend" & root_id %in% ids])
	)
	primary[, stage := "competent_authority_primary"]
	stopifnot(
		nrow(primary) == nrow(r),
		!anyDuplicated(primary$root_id),
		!any(primary$status == "not_applicable")
	)
	paths = data.table::rbindlist(list(paths, primary))
	pooled = data.table::copy(paths)
	pooled[, proposal_type := "all_eligible_types"]
	x = data.table::rbindlist(list(paths, pooled))
	tables = x[,
		{
			adopted_ids = successors[status == "adopted"]
			finals = unique(unlist(strsplit(adopted_ids, "|", fixed = TRUE), use.names = FALSE))
			list(
				considered_roots = .N,
				applicable_roots = sum(status != "not_applicable"),
				adopted_paths = sum(status == "adopted"),
				unique_adopted_successors = length(finals),
				rejected = sum(status == "rejected"),
				withdrawn_without_successor = sum(status == "withdrawn_without_successor"),
				pending = sum(status == "pending"),
				deferred = sum(status == "deferred"),
				unresolved_fates = sum(status == "unknown"),
				not_applicable = sum(status == "not_applicable")
			)
		},
		by = .(year, proposal_type, stage)
	]
	tables[, adoption_path_share := ifelse(applicable_roots > 0, adopted_paths / applicable_roots, NA_real_)]
	list(paths = paths, counts = tables)
}

proposal_year_comparisons = function(bounds, baselines = c(2023L, 2024L), target = 2025L) {
	b = data.table::copy(data.table::as.data.table(bounds))
	keys = c("comparison_basis", "proposal_type")
	stopifnot(
		all(
			c(
				keys,
				"year",
				"eligible_roots",
				"fully_retained_lower",
				"fully_retained_upper",
				"equal_root_request_retention_lower",
				"equal_root_request_retention_upper"
			) %in% names(b)
		),
		!anyDuplicated(b[, c(keys, "year"), with = FALSE])
	)
	result = lapply(
		baselines,
		function(base) {
			before = b[year == base]
			after = b[year == target]
			pair = merge(before, after, by = keys, all = TRUE, suffixes = c("_before", "_after"))
			pair[, `:=`(baseline = base, target = target)]
			pair[, .(
				comparison_basis,
				proposal_type,
				baseline,
				target,
				roots_before = eligible_roots_before,
				roots_after = eligible_roots_after,
				fully_retained_change_lower_pp = 100 * (fully_retained_lower_after - fully_retained_upper_before),
				fully_retained_change_upper_pp = 100 * (fully_retained_upper_after - fully_retained_lower_before),
				request_retention_change_lower_pp = 100 *
					(equal_root_request_retention_lower_after - equal_root_request_retention_upper_before),
				request_retention_change_upper_pp = 100 *
					(equal_root_request_retention_upper_after - equal_root_request_retention_lower_before)
			)]
		}
	)
	data.table::rbindlist(result, use.names = TRUE)
}

proposal_us_votes = function(roots, member_votes) {
	r = data.table::copy(data.table::as.data.table(roots))
	votes = data.table::as.data.table(member_votes)
	stopifnot(
		all(
			c("root_id", "year", "plenary_final_symbol", "plenary_adoption_date") %in% names(r)
		),
		all(
			c("resolution", "voting_country_code", "recorded_vote", "date", "undl_id") %in% names(votes)
		)
	)
	us = votes[
		resolution %in% r$plenary_final_symbol & voting_country_code == "USA",
		.(
			resolution,
			us_recorded_vote = recorded_vote,
			vote_date = as.character(date),
			undl_id
		)
	]
	stopifnot(
		!anyDuplicated(us$resolution),
		all(us$us_recorded_vote %in% c("Y", "N", "A", "X"))
	)
	x = merge(
		r[, .(
			root_id,
			cohort_year = year,
			plenary_final_symbol,
			plenary_adoption_date = as.character(plenary_adoption_date)
		)],
		us,
		by.x = "plenary_final_symbol",
		by.y = "resolution",
		all.x = TRUE
	)
	stopifnot(
		nrow(x) == nrow(r),
		!anyDuplicated(x$root_id),
		all(x[!is.na(us_recorded_vote), vote_date == plenary_adoption_date])
	)
	x[, position_basis := ifelse(
		is.na(us_recorded_vote),
		"no_member_vote_imputed; consult_actual_action_and_statement_records",
		"recorded_plenary_member_vote"
	)]
	x[]
}

proposal_matched_families = function(
	root_summary,
	family_map,
	baselines = c(2023L, 2024L),
	target = 2025L
) {
	x = data.table::copy(data.table::as.data.table(root_summary))
	f = data.table::copy(data.table::as.data.table(family_map))
	stopifnot(
		all(c("root_id", "family_id", "year", "review_status") %in% names(f)),
		!anyDuplicated(f$root_id),
		!anyNA(f$family_id),
		all(nzchar(f$family_id)),
		all(f$review_status == "adjudicated"),
		setequal(f$root_id, unique(x$root_id))
	)
	source_years = unique(x[, .(root_id, year)])
	check = merge(source_years, f[, .(root_id, family_year = year)], by = "root_id")
	stopifnot(nrow(check) == nrow(f), all(check$year == check$family_year))
	membership = list()
	content = list()
	counts = list()
	differences = list()
	for (base in baselines) {
		pair = f[year %in% c(base, target)]
		before = pair[year == base, .(roots_before = .N), by = family_id]
		after = pair[year == target, .(roots_after = .N), by = family_id]
		frame = merge(before, after, by = "family_id", all = TRUE)
		frame[is.na(roots_before), roots_before := 0L]
		frame[is.na(roots_after), roots_after := 0L]
		frame[, match_status := data.table::fcase(
			roots_before > 1L | roots_after > 1L,
			"multiple_roots_in_at_least_one_cohort",
			roots_before == 1L & roots_after == 1L,
			"one_root_in_each_cohort",
			roots_before == 0L,
			"absent_from_baseline_cohort",
			roots_after == 0L,
			"absent_from_target_cohort",
			default = "invalid"
		)]
		frame[, `:=`(baseline = base, target = target)]
		member = merge(pair, frame, by = "family_id", all.x = TRUE)
		member[, interpretation := "Family absence is not a failed proposal or evidence of substantive novelty."]
		membership[[as.character(base)]] = member
		counts[[as.character(base)]] = frame[, .(
			families = .N,
			roots_before = sum(roots_before),
			roots_after = sum(roots_after)
		),
		by = .(baseline, target, match_status)
		]
		ids = member[match_status == "one_root_in_each_cohort", root_id]
		if (!length(ids)) next
		selected = data.table::copy(x[root_id %in% ids])
		tab = proposal_content_tables(selected)
		tab$bounds[, `:=`(comparison_baseline = base, comparison_target = target)]
		content[[as.character(base)]] = tab$bounds
		differences[[as.character(base)]] = proposal_year_comparisons(tab$bounds, baselines = base, target = target)
	}
	list(
		membership = data.table::rbindlist(membership, use.names = TRUE),
		coverage = data.table::rbindlist(counts, use.names = TRUE),
		content_bounds = data.table::rbindlist(content, use.names = TRUE),
		comparisons = data.table::rbindlist(differences, use.names = TRUE)
	)
}

# summarize documented demands and amendments within the proposal analysis.
resolution_process = function(paths, read_member) {
	demands = read_analysis_table(
		file.path(paths$input, "resolution_demands.csv"),
		na.strings = c("", "NA")
	)
	strata = read_member(
		"proposal_study/adjudication/process/resolution_process_roots.csv"
	)
	conditions = read_member(
		"proposal_study/adjudication/process/resolution_condition_observability.csv"
	)
	out = file.path(paths$output, "adoptions/process")
	dir.create(out, recursive = TRUE, showWarnings = FALSE)
	truth = function(x) !is.na(x) & tolower(as.character(x)) == "true"
	stopifnot(
		!anyDuplicated(demands$canonical_demand_id),
		!anyDuplicated(strata$root_id),
		all(is.na(demands$root_id) | demands$root_id %in% strata$root_id),
		!anyDuplicated(conditions[, c("root_id", "expectation")])
	)
	for (v in c(
		"formal_amendment",
		"count_as_distinct_action",
		"us_explicitly_involved",
		"explicit_us_targeted_substantive_request",
		"request_observed_before_root_disposition"
	)) {
		data.table::set(demands, j = v, value = truth(demands[[v]]))
	}
	formal = demands[formal_amendment == TRUE & count_as_distinct_action == TRUE]
	stopifnot(
		!anyDuplicated(formal[, c("year", "amendment_id")]),
		all(formal$formal_disposition %in% c("adopted", "rejected", "withdrawn"))
	)
	amendment_summary = formal[, .(amendments = .N, roots = data.table::uniqueN(root_id)),
		by = c("year", "formal_disposition", "us_explicitly_involved")
	]
	demand_summary = demands[count_as_distinct_action == TRUE,
		.(records = .N, roots = data.table::uniqueN(root_id, na.rm = TRUE)),
		by = c(
			"year",
			"action_kind",
			"explicit_us_targeted_substantive_request",
			"request_observed_before_root_disposition",
			"us_demand_role"
		)
	]
	stratum_summary = strata[, .(roots = .N), by = c("year", "stratum")]
	observability = conditions[, .(roots = data.table::uniqueN(root_id)),
		by = c(
			"expectation",
			"eligible_2025_US_substantive_challenge",
			"prior_condition_value",
			"test_status"
		)
	]
	data.table::fwrite(
		formal,
		file.path(out, "formal_amendment_dispositions.csv")
	)
	data.table::fwrite(
		amendment_summary,
		file.path(out, "formal_amendment_summary.csv")
	)
	data.table::fwrite(
		demand_summary,
		file.path(out, "demand_record_summary.csv")
	)
	data.table::fwrite(stratum_summary, file.path(out, "us_position_strata.csv"))
	data.table::fwrite(
		observability,
		file.path(out, "conditional_expectation_observability.csv")
	)
	invisible(
		list(
			formal = formal,
			demands = demands,
			strata = strata,
			conditions = conditions
		)
	)
}

# reproduce the accepted proposal and process summaries from local inputs.
# no network requests, source classification or model estimation occurs here.
proposal_pipeline = function(input_root = "data", output_root = "results") {
	data.table::setDTthreads(1L)
	archive = file.path(input_root, "resolution_sources.zip")
	manifest_path = file.path(input_root, "resolution_sources.json")
	manifest = jsonlite::read_json(manifest_path, simplifyVector = FALSE)
	stopifnot(
		identical(
			digest::digest(file = archive, algo = "sha256"),
			manifest$archive_sha256
		)
	)
	members = utils::unzip(archive, list = TRUE)$Name
	stopifnot(!anyDuplicated(members))
	read_member = function(name) {
		stopifnot(name %in% members)
		data.table::as.data.table(
			read.csv(
				unz(archive, name),
				stringsAsFactors = FALSE,
				check.names = FALSE,
				na.strings = c("", "NA")
			)
		)
	}
	read_json_member = function(name) {
		stopifnot(name %in% members)
		connection = unz(archive, name)
		on.exit(close(connection))
		jsonlite::fromJSON(
			paste(readLines(connection, warn = FALSE), collapse = "\n"),
			simplifyVector = FALSE
		)
	}
	accepted = read_json_member("proposal_study/adjudication/acceptance.json")
	stopifnot(identical(accepted$status, "coordinator_accepted_for_derivation"))
	for (relative in names(accepted$input_hashes)) {
		stopifnot(
			!grepl("(^/|(^|/)\\.\\.(/|$))", relative),
			identical(
				digest::digest(file = relative, algo = "sha256"),
				accepted$input_hashes[[relative]]
			)
		)
	}
	roots = read_member("proposal_study/inventory/roots.csv")
	stage = read_member("proposal_study/inventory/content_packet_inventory.csv")
	coverage = read_member("proposal_study/adjudication/coverage.csv")
	families = read_member("proposal_study/adjudication/family_map.csv")
	commitments = read_analysis_table(
		file.path(input_root, "resolution_commitments.csv"),
		na.strings = c("", "NA")
	)
	proposals = read_analysis_table(
		file.path(input_root, "resolution_proposals.csv"),
		na.strings = c("", "NA")
	)
	links = read_analysis_table(
		file.path(input_root, "resolution_proposal_links.csv"),
		na.strings = c("", "NA")
	)
	demands = read_analysis_table(
		file.path(input_root, "resolution_demands.csv"),
		na.strings = c("", "NA")
	)
	as_flag = function(x, missing = FALSE) {
		value = tolower(as.character(x))
		stopifnot(all(is.na(value) | value %in% c("true", "false")))
		if (!missing) stopifnot(!anyNA(value))
		ifelse(is.na(value), NA, value == "true")
	}
	commitments[, material_request := as_flag(material_request, TRUE)]
	commitments[, count_in_denominator := as_flag(count_in_denominator)]
	coverage[, coding_complete := as_flag(coding_complete)]
	coverage[, material_inventory_complete := as_flag(material_inventory_complete)]
	stopifnot(
		nrow(roots) == 128L,
		!anyDuplicated(roots$root_id),
		identical(as.integer(table(roots$year)), c(45L, 41L, 42L)),
		!anyDuplicated(proposals$version_id),
		all(commitments$root_id %in% roots$root_id),
		all(demands$root_id[!is.na(demands$root_id)] %in% roots$root_id),
		!anyDuplicated(demands$canonical_demand_id)
	)
	fate = proposal_fate_summary(roots)
	fate_primary = fate$paths[stage == "competent_authority_primary", .(root_id, adopted_continuation = status == "adopted")]
	roots = merge(roots, fate_primary, by = "root_id", all.x = TRUE)
	complete = proposal_complete_stage_comparisons(commitments, coverage, stage)
	root_summary = proposal_content_summary(roots, complete$commitments, complete$coverage)
	tables = proposal_content_tables(root_summary)
	comparisons = proposal_year_comparisons(tables$bounds)
	matched = proposal_matched_families(root_summary, families)
	vote_environment = new.env(parent = emptyenv())
	load_vote_archive(
		file.path(input_root, "un_votes_long.rda"),
		envir = vote_environment
	)
	votes = proposal_us_votes(roots, vote_environment$un_votes_long)
	destination = file.path(output_root, "adoptions/proposals")
	dir.create(destination, recursive = TRUE, showWarnings = FALSE)
	outputs = list(
		roots = roots,
		fate_paths = fate$paths,
		fates = fate$counts,
		content_by_root = root_summary,
		content_categories = tables$counts,
		content_totals = tables$totals,
		content_bounds = tables$bounds,
		year_comparisons = comparisons,
		family_membership = matched$membership,
		family_coverage = matched$coverage,
		matched_content_bounds = matched$content_bounds,
		matched_year_comparisons = matched$comparisons,
		us_votes = votes
	)
	for (name in names(outputs)) {
		data.table::fwrite(
			outputs[[name]],
			file.path(destination, paste0(name, ".csv")),
			na = "NA"
		)
	}
	record_execution(
		list(input = input_root, output = destination),
		"proposals",
		file.path(
			input_root,
			c(
				"resolution_sources.zip",
				"resolution_sources.json",
				"resolution_proposals.csv",
				"resolution_proposal_links.csv",
				"resolution_commitments.csv",
				"resolution_demands.csv",
				"un_votes_long.rda"
			)
		),
		list(
			network_access = FALSE,
			models_refitted = FALSE,
			cohorts = 2023:2025,
			initial_basis = "earliest_documented_considered",
			sensitivity = "earliest_issued_all_roots",
			weighting = "equal_proposal; average_request_retention_within_proposal_first",
			missingness = "explicit_outcome_and_material_inventory_bounds",
			national_votes_imputed = FALSE,
			adoption_cutoff = "30June_following_cohort",
			source_acceptance = "proposal_study/adjudication/acceptance.json"
		)
	)
	resolution_process(
		list(input = input_root, output = output_root),
		read_member
	)
	invisible(outputs)
}

if (sys.nframe() == 0L && identical(commandArgs(trailingOnly = TRUE)[1L], "proposals")) {
	proposal_pipeline()
	quit(save = "no", status = 0L)
}

if (sys.nframe() == 0L && identical(commandArgs(trailingOnly = TRUE)[1L], "events")) {
	documentary_event_summary(
		"data/resolution_events.csv",
		"results/adoptions/content"
	)
	quit(save = "no", status = 0L)
}

if (sys.nframe() == 0L) {
	paths = pipeline_paths()
	Sys.setlocale("LC_TIME", "C")
	recent_adoptions(paths)
	adoption_histories(paths)
	resolution_content(paths)
	sources = c(
		file.path(paths$input, "un_votes_long.rda"),
		list.files(file.path(paths$input, "official"), full.names = TRUE),
		file.path(
			paths$input,
			c(
				"resolution_review.csv",
				"resolution_sources.zip",
				"resolution_sources.json"
			)
		),
		file.path(
			paths$output,
			"observed/vote_changes",
			c("individual_vote_accounting.csv", "matched_titles.csv")
		)
	)
	record_execution(
		paths,
		"adoptions",
		sources,
		list(
			calendar_years = 2010:2025,
			sessions = 64:80,
			consensus_member_votes_imputed = FALSE,
			network_access = FALSE
		)
	)
}
