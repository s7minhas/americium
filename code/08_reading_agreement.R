# compare request codings and continuity bounds

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
args = commandArgs(trailingOnly = TRUE)
archive = if (length(args) >= 1L) args[1L] else "data/proposal_comparison_sandbox.zip"
output_dir = if (length(args) >= 2L) args[2L] else "results/proposal_comparison/reading_agreement"
manifest_path = sub("[.]zip$", ".sha256", archive)
stopifnot(file.exists(archive), file.exists(manifest_path))
`%||%` = function(a, b) if (is.null(a)) b else a

sandbox_root = tempfile("americium_reading_agreement_")
dir.create(sandbox_root)
members = unzip(archive, list = TRUE)$Name
stopifnot(!any(grepl("(^/|(^|/)\\.\\.(/|$))", members)))
unzip(archive, exdir = sandbox_root)
manifest_lines = readLines(manifest_path, warn = FALSE)
manifest = data.table(
	sha256 = substr(manifest_lines, 1L, 64L),
	path = sub("^[0-9a-f]{64}  ", "", manifest_lines)
)
manifest[, actual := vapply(
	file.path(sandbox_root, path),
	function(p) digest(file = p, algo = "sha256"),
	character(1L)
)]
stopifnot(nrow(manifest) > 0L, all(manifest$actual == manifest$sha256))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
tracker = new.env()
tracker$members_read = character(0)
read_table = function(relative) {
	tracker$members_read = c(tracker$members_read, relative)
	read_analysis_table(
		file.path(sandbox_root, relative),
		colClasses = "character",
		na.strings = NULL,
		encoding = "UTF-8"
	)
}

split_ids = function(x) {
	lapply(
		strsplit(x, "[;|]"),
		function(v) {
			v = trimws(v)
			v[nzchar(v)]
		}
	)
}

norm_paragraph = function(x) {
	y = toupper(gsub("\\s+", "", x))
	hit = grepl("^(OP|PP)[0-9]+", y)
	y[hit] = regmatches(y, regexpr("^(OP|PP)[0-9]+", y))
	y
}

root_from_id = function(x) {
	m = regexpr("^C2_[0-9]+_L[0-9]+", x)
	out = rep(NA_character_, length(x))
	out[m > 0L] = regmatches(x, m)
	out
}

seven = c(
	"continued_without_substantive_change",
	"continued_stronger_specified_requirements",
	"continued_weaker_specified_requirements",
	"activity_changed_not_orderable",
	"completed_or_expired",
	"no_corresponding_request_found",
	"unresolved"
)
alias = c(
	unchanged = seven[1],
	continued = seven[1],
	continued_via_verified_incorporation = seven[1],
	stronger = seven[2],
	continued_with_stronger_specified_requirements = seven[2],
	weaker = seven[3],
	continued_weaker = seven[3],
	continued_with_weaker_specified_requirements = seven[3],
	postponed_pending_task = seven[3],
	not_orderable = seven[4],
	changed_not_orderable = seven[4],
	activity_changed_in_way_cannot_be_ordered = seven[4],
	completed = seven[5],
	completed_task = seven[5],
	expired_opportunity = seven[5],
	no_corresponding_request = seven[6],
	pending_task_not_reissued = seven[7]
)
harmonize_outcome = function(x) {
	out = unname(alias[x])
	out[x %in% seven] = x[x %in% seven]
	out
}

noncore_labels = c(
	"boundary_only",
	"excluded_boundary_or_container",
	"administrative_separate",
	"duplicate",
	"excluded_from_substantive_comparison",
	"administrative_renewal"
)
addition_labels = c(
	"addition",
	"boundary_addition",
	"new_request",
	"new_boundary_or_context",
	"addition_or_newly_explicit_candidate",
	"unresolved_addition",
	"addition_or_noncorecandidate"
)
first_core_status = c("included", "count", "candidate_count")
independent_core_status = c("core", "count", "included", "request", "candidate", "time_limited")
pairs = read_table("work/content/family_pairs.csv")
pairs[, `:=`(
	prior_year = as.integer(prior_year),
	later_year = as.integer(later_year),
	eligible = tolower(eligible_one_to_one) == "true",
	internal = tolower(internal_decision) == "true",
	common = tolower(common_three_years) == "true"
)]
substantive_pairs = pairs[eligible & !internal]
root_lookup = substantive_pairs[, setNames(pair_id, paste(prior_root, later_root))]
pair_from_roots = function(prior_root, later_root) {
	unname(root_lookup[paste(prior_root, later_root)])
}

batches = read_table("work/review/accepted_batches.csv")
batches[, dir := dirname(path)]
assignments = read_table("work/review/batch_assignments.csv")
batch_names = batches$batch
first_sources = list(
	pilot_a = list(
		correspondence = c(
			"work/content/first_reader_frozen_request_correspondence.csv",
			"work/content_002_010/request_correspondence.csv"
		),
		units = c(
			"work/content/first_reader_frozen_adopted_requests.csv",
			"work/content_002_010/adopted_requests.csv"
		)
	),
	pilot_b = list(
		correspondence = c(
			"work/content/ict_first_reader_frozen_correspondence.csv",
			"work/content/poverty_first_reader_frozen_correspondence.csv"
		),
		units = c(
			"work/content/ict_first_reader_frozen_requests.csv",
			"work/content/poverty_first_reader_frozen_requests.csv"
		)
	),
	pilot_climate = list(
		correspondence = "work/content_048/request_correspondence.csv",
		units = "work/content_048/adopted_requests.csv"
	)
)
for (b in setdiff(batch_names, names(first_sources))) {
	first_sources[[b]] = list(
		correspondence = paste0("work/expansion_", b, "/request_correspondence.csv"),
		units = paste0("work/expansion_", b, "/adopted_requests.csv")
	)
}

load_accepted = function(batch, dir) {
	units = read_table(file.path(dir, "accepted", "adopted_requests.csv"))
	units[, core := count_status %in% c("included", "count") & category %in% c("substantive", "substantive_request")]
	rows = read_table(file.path(dir, "accepted", "request_correspondence.csv"))[nzchar(prior_request_id)]
	rows = merge(
		rows[, .(pair_id, accepted_id = prior_request_id, accepted_outcome = outcome)],
		units[, .(
			accepted_id = request_id,
			root_id,
			paragraph = norm_paragraph(paragraph),
			core
		)],
		by = "accepted_id",
		all.x = TRUE
	)
	rows[, batch := batch]
	rows[pair_id %in% substantive_pairs$pair_id]
}

load_first = function(batch) {
	src = first_sources[[batch]]
	units = rbindlist(lapply(src$units, read_table), fill = TRUE, use.names = TRUE)
	units = data.table(
		batch = batch,
		unit_id = units$request_id,
		root_id = units$root_id,
		paragraph = norm_paragraph(units$paragraph),
		status = units$count_status,
		core = units$count_status %in% first_core_status
	)
	rows = rbindlist(
		lapply(src$correspondence, read_table),
		fill = TRUE,
		use.names = TRUE
	)[nzchar(prior_request_id)]
	rows = data.table(
		batch = batch,
		pair_id = rows$pair_id,
		prior_id = rows$prior_request_id,
		outcome_raw = rows$outcome,
		row_status = rows$count_status
	)
	rows = merge(
		rows,
		units[, .(prior_id = unit_id, root_id, paragraph, unit_status = status)],
		by = "prior_id",
		all.x = TRUE
	)
	rows[, core := row_status %in% first_core_status]
	rows[, outcome := harmonize_outcome(outcome_raw)]
	list(units = units, rows = rows)
}

load_independent = function(batch, dir, family_ids) {
	x = read_table(file.path(dir, "independent_units.csv"))
	id_col = intersect(c("independent_id", "unit_id", "request_id"), names(x))[1L]
	root = if ("root_id" %in% names(x)) x$root_id else sub("_adopted$", "", x$source)
	paragraph = if ("paragraph" %in% names(x)) x$paragraph else x$locator
	status = if ("status" %in% names(x)) x$status else x$count_status
	units = data.table(
		batch = batch,
		unit_id = x[[id_col]],
		root_id = root,
		paragraph = norm_paragraph(paragraph),
		status = status,
		core = status %in% independent_core_status
	)
	y = read_table(file.path(dir, "independent_comparisons.csv"))
	prior_col = intersect(
		c("prior_id", "prior_request_id", "prior_unit", "prior_unit_id"),
		names(y)
	)[1L]
	outcome_col = intersect(c("outcome", "classification", "label"), names(y))[1L]
	status_col = intersect(c("prior_status", "count_status", "status"), names(y))
	if ("pair_id" %in% names(y)) {
		pair = y$pair_id
		if (any(grepl("__", pair))) {
			parts = tstrsplit(pair, "__")
			pair = pair_from_roots(parts[[1L]], parts[[2L]])
		}
	} else if (all(c("prior_source", "later_source") %in% names(y))) {
		pair = pair_from_roots(
			sub("_adopted$", "", y$prior_source),
			sub("_adopted$", "", y$later_source)
		)
	} else if (all(c("prior_root", "later_root") %in% names(y))) {
		pair = pair_from_roots(y$prior_root, y$later_root)
	} else {
		family = if ("family_id" %in% names(y)) y$family_id else family_ids[1L]
		pair = paste0(family, "_", gsub("-", "_", y$pair))
	}
	rows = data.table(
		batch = batch,
		pair_id = pair,
		prior_id = y[[prior_col]],
		outcome_raw = y[[outcome_col]],
		row_status = if (length(status_col)) y[[status_col[1L]]] else NA_character_
	)[nzchar(prior_id)]
	rows = merge(
		rows,
		units[, .(
			prior_id = unit_id,
			root_id,
			paragraph,
			unit_status = status,
			unit_core = core
		)],
		by = "prior_id",
		all.x = TRUE
	)
	rows[is.na(root_id), root_id := root_from_id(prior_id)]
	rows[, core := fifelse(!is.na(unit_core), unit_core, row_status %in% independent_core_status)]
	rows[, outcome := harmonize_outcome(outcome_raw)]
	list(units = units, rows = rows)
}

link = function(reader, original, accepted) {
	reader = rep_len(reader, length(original))
	o = split_ids(original)
	a = split_ids(accepted)
	out = rbindlist(
		lapply(
			seq_along(o),
			function(i) {
				if (!length(o[[i]]) || !length(a[[i]])) {
					return(NULL)
				}
				CJ(original_id = o[[i]], accepted_id = a[[i]])[, reader := reader[i]]
			}
		)
	)
	if (!nrow(out)) {
		return(
			data.table(
				original_id = character(0),
				accepted_id = character(0),
				reader = character(0)
			)
		)
	}
	unique(out)
}

load_links = function(batch, dir) {
	acc = file.path(dir, "accepted")
	out = switch(
	batch,
		a01 = {
			x = read_table(file.path(acc, "id_crosswalk.csv"))
			link("first", x$first_id, x$accepted_id)
		},
		a02 = {
			x = read_table(file.path(acc, "crosswalk.csv"))
			link("independent", x$independent_id, x$accepted_paragraph_candidates)
		},
		a03 = ,
		a04 = {
			x = read_table(file.path(acc, "id_crosswalk.csv"))
			link(x$reading, x$original_id, x$accepted_ids)
		},
		b01 = {
			x = read_table(file.path(acc, "id_crosswalk.csv"))
			link(x$reader, x$original_id, x$accepted_ids)
		},
		b02 = {
			x = read_table(file.path(acc, "id_crosswalk.csv"))
			link(x$reader, x$original_id, x$accepted_id)
		},
		b03 = {
			x = read_table(file.path(acc, "id_crosswalk.csv"))
			rbind(
				link("first", x$first_reader_ids, x$accepted_ids),
				link("independent", x$independent_id, x$accepted_ids)
			)
		},
		b04 = {
			x = read_table(file.path(acc, "id_crosswalk.csv"))
			link("independent", x$independent_id, x$accepted_ids)
		},
		c01 = ,
		c02 = ,
		c03 = {
			x = read_table(file.path(acc, "id_crosswalk.csv"))
			link("first", x$first_id, x$accepted_id)
		},
		c04 = {
			x = read_table(file.path(acc, "id_crosswalk.csv"))
			rbind(
				link(
					"first",
					x[relation %in% c("retained_id", "split"), original_request_id],
					x[relation %in% c("retained_id", "split"), accepted_request_id]
				),
				link(
					"independent",
					x[relation == "independent_to_accepted_action", original_request_id],
					x[relation == "independent_to_accepted_action", accepted_request_id]
				)
			)
		},
		root01 = {
			x = read_table(file.path(dir, "reconciliation.csv"))
			rbind(
				link("first", x$first_id, x$accepted_ids),
				link("independent", x$independent_ids, x$accepted_ids)
			)
		},
		root02 = {
			x = read_table(file.path(acc, "crosswalk.csv"))
			link("independent", x$independent_id, x$accepted_ids)
		},
		root03 = {
			x = read_table(file.path(acc, "crosswalk.csv"))
			link("independent", x$independent_id, x$accepted_id)
		},
		pilot_a = {
			x = read_table(file.path(acc, "id_crosswalk.csv"))
			rbind(
				link("first", x$first_reader_id, x$accepted_ids),
				link("independent", x$independent_same_paragraph_ids, x$accepted_ids)
			)
		},
		pilot_b = {
			x = read_table(file.path(dir, "reconciliation.csv"))
			rbind(
				link("first", x$first_reader_id, x$accepted_id),
				link("independent", x$independent_reader_ids, x$accepted_id)
			)
		},
		pilot_climate = {
			x = read_table(file.path(acc, "id_crosswalk.csv"))
			link(
				fifelse(x$reader == "second", "independent", "first"),
				x$original_id,
				x$accepted_ids
			)
		}
	)
	out[, batch := batch]
	out
}

load_direct_first = function(batch, dir) {
	if (!batch %in% c("b03", "b04")) {
		return(NULL)
	}
	x = read_table(file.path(dir, "accepted", "reconciliation.csv"))
	x = x[!first_reader_outcome %in% c("newsplit", "new_split_or_candidate")]
	data.table(
		batch = batch,
		pair_id = x$pair_id,
		accepted_id = x$request_id,
		outcome_raw = x$first_reader_outcome,
		outcome = harmonize_outcome(x$first_reader_outcome)
	)
}

accepted = list()
first_units = list()
first_rows = list()
indep_units = list()
indep_rows = list()
links = list()
direct_first = list()
for (i in seq_len(nrow(batches))) {
	b = batches$batch[i]
	d = batches$dir[i]
	families = strsplit(batches$families[i], "[|]")[[1L]]
	accepted[[b]] = load_accepted(b, d)
	f = load_first(b)
	first_units[[b]] = f$units
	first_rows[[b]] = f$rows
	s = load_independent(b, d, families)
	indep_units[[b]] = s$units
	indep_rows[[b]] = s$rows
	links[[b]] = load_links(b, d)
	direct_first[[b]] = load_direct_first(b, d)
}

accepted = rbindlist(accepted)
first_units = rbindlist(first_units)
first_rows = rbindlist(first_rows)
indep_units = rbindlist(indep_units)
indep_rows = rbindlist(indep_rows)
links = rbindlist(links)
direct_first = rbindlist(direct_first)
identity_links = rbind(
	unique(
		first_units[unit_id %in% accepted$accepted_id, .(original_id = unit_id, accepted_id = unit_id, reader = "first", batch)]
	),
	unique(
		indep_units[unit_id %in% accepted$accepted_id, .(original_id = unit_id, accepted_id = unit_id, reader = "independent", batch)]
	)
)
links = unique(rbind(links, identity_links))
accepted = merge(
	accepted,
	substantive_pairs[, .(pair_id, family_id, prior_year, later_year, common)],
	by = "pair_id"
)
stopifnot(
	!anyNA(accepted$core),
	!anyDuplicated(accepted[, .(batch, pair_id, accepted_id)])
)
core_accepted = accepted[core == TRUE]

align_reading = function(rows, units, reader_links, direct = NULL) {
	by_id = merge(
		core_accepted[, .(batch, pair_id, accepted_id)],
		reader_links[, .(batch, accepted_id, prior_id = original_id)],
		by = c("batch", "accepted_id"),
		allow.cartesian = TRUE
	)
	by_id = merge(
		by_id,
		rows[, .(batch, pair_id, prior_id, core, outcome)],
		by = c("batch", "pair_id", "prior_id")
	)
	by_id[, route := "identifier"]
	by_paragraph = merge(
		core_accepted[, .(batch, pair_id, accepted_id, root_id, paragraph)],
		rows[!is.na(paragraph), .(batch, pair_id, prior_id, root_id, paragraph, core, outcome)],
		by = c("batch", "pair_id", "root_id", "paragraph"),
		allow.cartesian = TRUE
	)
	by_paragraph[, route := "paragraph"]
	candidates = rbind(
		by_id,
		by_paragraph[, .(batch, pair_id, accepted_id, prior_id, core, outcome, route)]
	)
	if (!is.null(direct) && nrow(direct)) {
		candidates = rbind(
			candidates,
			direct[, .(
				batch,
				pair_id,
				accepted_id,
				prior_id = NA_character_,
				core = outcome %in% seven,
				outcome,
				route = "record"
			)]
		)
	}
	core_candidates = candidates[core == TRUE & !is.na(outcome)]
	summary = core_candidates[, .(
		identified = TRUE,
		route = paste(sort(unique(route)), collapse = "+"),
		outcome = {
			all_outcomes = unique(outcome)
			linked = unique(outcome[route %in% c("identifier", "record")])
			if (length(all_outcomes) == 1L) all_outcomes else if (length(linked) == 1L) linked else "ambiguous"
		},
		matched_units = uniqueN(prior_id[!is.na(prior_id)])
	),
	by = .(batch, pair_id, accepted_id)
	]
	# a unit of the reading at the paragraph that it did not count as a request
	present = rbind(
		merge(
			core_accepted[, .(batch, pair_id, accepted_id, root_id, paragraph)],
			units[, .(batch, root_id, paragraph, unit_id)],
			by = c("batch", "root_id", "paragraph")
		)[, .(batch, pair_id, accepted_id)],
		merge(
			core_accepted[, .(batch, pair_id, accepted_id)],
			reader_links[, .(batch, accepted_id, unit_id = original_id)],
			by = c("batch", "accepted_id")
		)[, .(batch, pair_id, accepted_id)]
	)
	present = unique(present)[, unit_present := TRUE]
	out = merge(
		core_accepted[, .(batch, pair_id, accepted_id)],
		summary,
		by = c("batch", "pair_id", "accepted_id"),
		all.x = TRUE
	)
	out = merge(out, present, by = c("batch", "pair_id", "accepted_id"), all.x = TRUE)
	out[is.na(identified), identified := FALSE]
	out[is.na(unit_present), unit_present := FALSE]
	out[identified == FALSE, route := "none"]
	matched_units = unique(core_candidates[!is.na(prior_id), .(batch, pair_id, prior_id)])
	list(alignment = out, matched_units = matched_units)
}

first_alignment = align_reading(first_rows, first_units, links[reader == "first"], direct_first)
independent_alignment = align_reading(indep_rows, indep_units, links[reader == "independent"])
setnames(
	first_alignment$alignment,
	c("identified", "route", "outcome", "matched_units", "unit_present"),
	c(
		"first_identified",
		"first_route",
		"first_outcome",
		"first_matched_units",
		"first_unit_present"
	)
)
setnames(
	independent_alignment$alignment,
	c("identified", "route", "outcome", "matched_units", "unit_present"),
	c(
		"independent_identified",
		"independent_route",
		"independent_outcome",
		"independent_matched_units",
		"independent_unit_present"
	)
)
requests = merge(
	core_accepted,
	first_alignment$alignment,
	by = c("batch", "pair_id", "accepted_id")
)
requests = merge(
	requests,
	independent_alignment$alignment,
	by = c("batch", "pair_id", "accepted_id")
)
requests[, both_identified := first_identified & independent_identified]
requests[, outcomes_compared := both_identified & first_outcome != "ambiguous" & independent_outcome != "ambiguous"]
requests[, readings_agree := fifelse(outcomes_compared, first_outcome == independent_outcome, NA)]
requests[, continued_accepted := accepted_outcome == seven[1]]
requests[, continued_both_readings := outcomes_compared & first_outcome == seven[1] & independent_outcome == seven[1]]

reading_only = function(rows, matched_units, label) {
	core_rows = rows[core == TRUE & !is.na(outcome) & pair_id %in% substantive_pairs$pair_id]
	core_rows = merge(
		core_rows,
		matched_units[, .(batch, pair_id, prior_id, matched = TRUE)],
		by = c("batch", "pair_id", "prior_id"),
		all.x = TRUE
	)
	core_rows[is.na(matched), matched := FALSE]
	core_rows[, reading := label]
	core_rows
}

first_only = reading_only(first_rows, first_alignment$matched_units, "first")
independent_only = reading_only(indep_rows, independent_alignment$matched_units, "independent")

first_links = merge(
	links[reader == "first"],
	first_units[, .(batch, original_id = unit_id, first_core = core)],
	by = c("batch", "original_id")
)
first_links = first_links[accepted_id %in% core_accepted$accepted_id]
first_links[, children := uniqueN(accepted_id), by = .(batch, original_id)]
origin = first_links[, .(origins = .N, children = max(children), any_core = any(first_core)), by = .(batch, accepted_id)]
origin[, unit_origin := fifelse(
	origins >= 2L,
	"merged from first-reading units",
	fifelse(
		children >= 2L,
		"split from a first-reading unit",
		fifelse(
			any_core,
			"first-reading unit retained",
			"first-reading unit promoted to request"
		)
	)
)]
requests = merge(
	requests,
	origin[, .(batch, accepted_id, unit_origin)],
	by = c("batch", "accepted_id"),
	all.x = TRUE
)
requests[is.na(unit_origin), unit_origin := "no first-reading identifier"]
requests[, independent_presence := fifelse(
	independent_identified,
	"independent request at the same paragraph",
	fifelse(
		independent_unit_present,
		"independent unit at the paragraph not counted as a request",
		"no independent unit at the paragraph"
	)
)]
requests[, outcome_adoption := fifelse(
	!outcomes_compared,
	"not compared",
	fifelse(
		readings_agree & accepted_outcome == first_outcome,
		"readings agreed; accepted follows them",
		fifelse(
			readings_agree,
			"readings agreed; accepted differs from both",
			fifelse(
				accepted_outcome == first_outcome,
				"readings disagreed; accepted follows the first reading",
				fifelse(
					accepted_outcome == independent_outcome,
					"readings disagreed; accepted follows the independent reading",
					"readings disagreed; accepted follows neither"
				)
			)
		)
	)
)]
cohen_kappa = function(a, b) {
	if (!length(a)) {
		return(NA_real_)
	}
	cats = sort(unique(c(a, b)))
	tab = table(factor(a, cats), factor(b, cats))
	n = sum(tab)
	po = sum(diag(tab)) / n
	pe = sum(rowSums(tab) * colSums(tab)) / n^2
	if (pe >= 1) {
		return(NA_real_)
	}
	(po - pe) / (1 - pe)
}

summarize_scope = function(d, scope, batch = NA_character_) {
	compared = d[outcomes_compared == TRUE]
	fo = first_only[pair_id %in% unique(d$pair_id) & batch %in% unique(d$batch)]
	io = independent_only[pair_id %in% unique(d$pair_id) & batch %in% unique(d$batch)]
	data.table(
		scope = scope,
		batch = batch,
		families = uniqueN(d$family_id),
		pairs = uniqueN(d$pair_id),
		accepted_core_prior_requests = nrow(d),
		first_identified = sum(d$first_identified),
		independent_identified = sum(d$independent_identified),
		both_identified = sum(d$both_identified),
		share_first_identified = mean(d$first_identified),
		share_independent_identified = mean(d$independent_identified),
		share_both_identified = mean(d$both_identified),
		first_request_units_without_accepted_request = sum(!fo$matched),
		independent_request_units_without_accepted_request = sum(!io$matched),
		independent_request_units = nrow(io),
		ambiguous_outcomes = sum(d$both_identified & !d$outcomes_compared),
		outcomes_compared = nrow(compared),
		outcomes_agree = sum(compared$readings_agree),
		percent_agreement = 100 * mean(compared$readings_agree),
		kappa = cohen_kappa(compared$first_outcome, compared$independent_outcome),
		continued_agree = sum(
			(compared$first_outcome == seven[1]) == (compared$independent_outcome == seven[1])
		),
		percent_agreement_continued = 100 * mean(
			(compared$first_outcome == seven[1]) == (compared$independent_outcome == seven[1])
		),
		kappa_continued = cohen_kappa(
			compared$first_outcome == seven[1],
			compared$independent_outcome == seven[1]
		)
	)
}

scopes = list()
for (b in batch_names) scopes[[b]] = summarize_scope(requests[batch == b], "batch", b)
common = requests[common == TRUE]
for (y in list(c(2023L, 2024L), c(2024L, 2025L), c(2023L, 2025L))) {
	label = paste0("common_families_", y[1L], "_", y[2L])
	scopes[[label]] = summarize_scope(common[prior_year == y[1L] & later_year == y[2L]], label)
}

scopes[["all_substantive_pairs"]] = summarize_scope(requests, "all_substantive_pairs")
batch_agreement = rbindlist(scopes)
# confusion matrices of the two readings over the seven categories
confusion = function(d, scope, batch = NA_character_) {
	d[outcomes_compared == TRUE, .(n = .N), by = .(first_outcome, independent_outcome)][, `:=`(scope = scope, batch = batch)]
}

outcome_agreement = rbind(
	rbindlist(
		lapply(batch_names, function(b) confusion(requests[batch == b], "batch", b))
	),
	rbindlist(
		lapply(
			list(c(2023L, 2024L), c(2024L, 2025L), c(2023L, 2025L)),
			function(y) {
				confusion(
					common[prior_year == y[1L] & later_year == y[2L]],
					paste0("common_families_", y[1L], "_", y[2L])
				)
			}
		)
	),
	confusion(requests, "all_substantive_pairs")
)
setcolorder(
	outcome_agreement,
	c("scope", "batch", "first_outcome", "independent_outcome", "n")
)
setorder(outcome_agreement, scope, batch, first_outcome, independent_outcome)

decision_table = function(d, scope, batch = NA_character_) {
	unit_level = unique(d[, .(batch, accepted_id, unit_origin, independent_presence)])
	rbind(
		unit_level[, .(level = "unit origin relative to the first reading", n = .N), by = .(type = unit_origin)][, .(level, type, n)],
		unit_level[, .(level = "independent reading at the accepted request", n = .N), by = .(type = independent_presence)][, .(level, type, n)],
		d[, .(level = "outcome adoption in the accepted coding", n = .N), by = .(type = outcome_adoption)][, .(level, type, n)],
		data.table(
			level = "first-reading request units without an accepted core request",
			type = "dropped or demoted in reconciliation",
			n = sum(
				!first_only[batch %in% unique(d$batch) & pair_id %in% unique(d$pair_id), matched]
			)
		),
		data.table(
			level = "independent request units without an accepted core request",
			type = "not adopted as a core request",
			n = sum(
				!independent_only[batch %in% unique(d$batch) & pair_id %in% unique(d$pair_id), matched]
			)
		)
	)[, `:=`(scope = scope, batch = batch)]
}

reconciliation_decisions = rbind(
	rbindlist(
		lapply(
			batch_names,
			function(b) decision_table(requests[batch == b], "batch", b)
		)
	),
	decision_table(requests, "all_substantive_pairs")
)
reconciliation_decisions[, share := n / sum(n), by = .(scope, batch, level)]
setcolorder(
	reconciliation_decisions,
	c("scope", "batch", "level", "type", "n", "share")
)

bound_family = common[, .(
	requests = .N,
	both_identified = sum(both_identified),
	outcomes_compared = sum(outcomes_compared),
	disagreements = sum(outcomes_compared & !readings_agree),
	continued_accepted = sum(continued_accepted),
	continued_disagreement_bound = sum(continued_accepted & !(outcomes_compared & !readings_agree)),
	continued_confirmed_bound = sum(continued_accepted & continued_both_readings)
),
by = .(prior_year, later_year, family_id)
]
bound_family[, `:=`(
	accepted_share = continued_accepted / requests,
	disagreement_bound_share = continued_disagreement_bound / requests,
	confirmed_bound_share = continued_confirmed_bound / requests
)]
bound_summary = bound_family[, .(
	family_id = "equal_family_mean",
	requests = sum(requests),
	both_identified = sum(both_identified),
	outcomes_compared = sum(outcomes_compared),
	disagreements = sum(disagreements),
	continued_accepted = sum(continued_accepted),
	continued_disagreement_bound = sum(continued_disagreement_bound),
	continued_confirmed_bound = sum(continued_confirmed_bound),
	accepted_share = mean(accepted_share),
	disagreement_bound_share = mean(disagreement_bound_share),
	confirmed_bound_share = mean(confirmed_bound_share),
	families = .N
), by = .(prior_year, later_year)]
bound_family[, families := 1L]
continuity_bound = rbind(
	bound_summary[, scope := "equal_family_mean"],
	bound_family[, scope := "family"],
	use.names = TRUE
)
continuity_bound[, `:=`(
	accepted_percent = round(100 * accepted_share, 2),
	disagreement_bound_percent = round(100 * disagreement_bound_share, 2),
	confirmed_bound_percent = round(100 * confirmed_bound_share, 2)
)]
setcolorder(
	continuity_bound,
	c("scope", "prior_year", "later_year", "family_id", "families", "requests")
)
setorder(continuity_bound, -scope, prior_year, later_year, family_id)
headline_check = list(available = FALSE)
headline_path = file.path(dirname(output_dir), "headline.csv")
if (file.exists(headline_path)) {
	headline = read_analysis_table(headline_path)[population == "common_three_year_families" & subset == "all_core" & outcome == seven[1]]
	check = merge(
		headline[, .(prior_year, later_year, headline_fraction = mean_family_fraction)],
		bound_summary[, .(prior_year, later_year, accepted_share)],
		by = c("prior_year", "later_year")
	)
	headline_check = list(
		available = TRUE,
		path = headline_path,
		max_absolute_difference = max(abs(check$headline_fraction - check$accepted_share))
	)
	stopifnot(headline_check$max_absolute_difference < 1e-9)
}

freeze_files = c(
	a01 = "INDEPENDENT_FREEZE.sha256",
	a02 = "freeze.sha256",
	a03 = "independent_freeze.json",
	a04 = "independent_freeze.json",
	b01 = "independent_freeze.json",
	b02 = "independent_freeze.json",
	b03 = "independent_freeze.json",
	b04 = "independent_freeze.json",
	c01 = "INDEPENDENT_FREEZE.sha256",
	c02 = "INDEPENDENT_FREEZE.sha256",
	c03 = "INDEPENDENT_FREEZE.sha256",
	c04 = "independent_freeze.json",
	root01 = "freeze.sha256",
	root02 = "freeze.sha256",
	root03 = "freeze.sha256",
	pilot_a = "INDEPENDENT_FREEZE.sha256",
	pilot_b = "freeze.sha256",
	pilot_climate = "freeze_receipt.json"
)
read_freeze = function(relative) {
	tracker$members_read = c(tracker$members_read, relative)
	path = file.path(sandbox_root, relative)
	if (grepl("[.]json$", path)) {
		j = fromJSON(path)
		hashes = j[["files"]] %||% j[["sha256"]] %||% j[["hashes"]]
		if (is.null(hashes)) {
			hashes = Filter(
				function(v) is.character(v) && length(v) == 1L && grepl("^[0-9a-f]{64}$", v),
				j
			)
		}
		stamp = j[["time_utc"]] %||% j[["utc"]] %||% j[["frozen_utc"]]
		data.table(
			file = names(hashes),
			sha256 = unlist(hashes),
			time_utc = stamp %||% NA_character_
		)
	} else if (grepl("[.]csv$", path)) {
		x = read_analysis_table(path, colClasses = "character")
		data.table(
			file = basename(x$file),
			sha256 = x$sha256,
			time_utc = NA_character_
		)
	} else {
		lines = readLines(path, warn = FALSE)
		hit = grepl("^[0-9a-f]{64}  ", lines)
		stamp = sub("^UTC ", "", lines[grepl("^UTC ", lines)])
		data.table(
			file = sub("^[0-9a-f]{64}  ", "", lines[hit]),
			sha256 = substr(lines[hit], 1L, 64L),
			time_utc = if (length(stamp)) stamp[1L] else NA_character_
		)
	}
}

first_freeze_files = c(
	a01 = "work/expansion_a01/first_reader_freeze.json",
	a02 = "work/expansion_a02/first_reader_freeze.json",
	a03 = "work/expansion_a03/first_reader_freeze.json",
	a04 = "work/expansion_a04/first_reader_freeze.json",
	b01 = "work/expansion_b01/FIRST_READER_FREEZE.sha256",
	b02 = "work/expansion_b02/FIRST_READER_FREEZE.sha256",
	b03 = "work/expansion_b03/FIRST_READER_FREEZE.sha256",
	b04 = "work/expansion_b04/FIRST_READER_FREEZE.sha256",
	c01 = "work/expansion_c01/first_reader_freeze.sha256",
	c02 = "work/expansion_c02/first_reader_freeze.sha256",
	c03 = "work/expansion_c03/first_reader_freeze.sha256",
	c04 = "work/expansion_c04/first_reader_freeze.json",
	root01 = "work/expansion_root01/first_reader_freeze.json",
	root02 = "work/expansion_root02/tax_first_reader_freeze.json",
	root03 = "work/expansion_root03/first_reader_freeze.json",
	pilot_climate = "work/content_048/freeze.csv"
)
first_reading_freeze = rbindlist(
	lapply(
		batch_names,
		function(b) {
			src = first_sources[[b]]
			files = data.table(
				batch = b,
				role = rep(
					c("correspondence", "units"),
					c(length(src$correspondence), length(src$units))
				),
				file = c(src$correspondence, src$units)
			)
			files[, current_sha256 := vapply(
				file.path(sandbox_root, file),
				function(p) digest(file = p, algo = "sha256"),
				character(1L)
			)]
			files[, `:=`(freeze_record = NA_character_, recorded_sha256 = NA_character_)]
			if (b %in% names(first_freeze_files)) {
				frozen = read_freeze(first_freeze_files[[b]])
				files[, recorded_sha256 := frozen$sha256[match(basename(file), frozen$file)]]
				files[!is.na(recorded_sha256), freeze_record := first_freeze_files[[b]]]
			}
			files[, match_freeze := fifelse(is.na(recorded_sha256), NA, recorded_sha256 == current_sha256)]
			files
		}
	)
)
first_freeze_summary = first_reading_freeze[, .(
	first_reading_freeze_record = paste(unique(freeze_record[!is.na(freeze_record)]), collapse = "|"),
	first_reading_files = .N,
	first_reading_files_with_record = sum(!is.na(recorded_sha256)),
	first_reading_files_match_record = sum(match_freeze %in% TRUE),
	first_reading_files_without_record = paste(file[is.na(recorded_sha256)], collapse = "|")
), by = batch]
findings_notes = data.table(batch = batch_names)
findings_notes[, first_labels_withheld_until_freeze := TRUE]
findings_notes[, fresh_context := FALSE]
findings_notes[, context_reused := c(
	pilot_a = "codebook and source metadata read before coding; first labels withheld until the freeze",
	pilot_b = "same agent had read the ICT bargaining chronology; first labels withheld until the freeze",
	pilot_climate = "agent context that had coded other families; first labels withheld until the freeze",
	root01 = "reused model-agent context disclosed in context.txt; first labels withheld until the freeze",
	c01 = "prior climate, education, energy and disaster-risk context reused; first labels withheld until the freeze",
	b01 = "prior-family context reused; a generic policy-container rule was forwarded before the freeze; first labels withheld",
	a01 = "earlier-family and codebook context reused; first labels withheld until the freeze",
	root02 = "general process context reused; first labels withheld until the freeze",
	a02 = "first labels withheld until the freeze; no fresh-context claim",
	c02 = "shared codebook and earlier-family context; one external text read before the freeze; first labels withheld",
	b02 = "context reused from other families; first labels withheld until the freeze",
	root03 = "internal decisions only; features retrospective; first labels withheld until the freeze",
	c03 = "shared earlier-family and codebook context; first labels withheld until the freeze",
	a03 = "earlier process, tax and Sevilla familiarity retained; first labels withheld until the freeze",
	b03 = "context reused from other families; first labels withheld until the freeze",
	c04 = "shared codebook and previous family context; first labels withheld until the freeze",
	b04 = "context reused from other families; first labels withheld until the freeze; no prospective feature-freeze claim for the 2023 SIDS text",
	a04 = "prior process, ICT and Sevilla context disclosed; first labels withheld until the freeze"
)[batch]]
independence = rbindlist(
	lapply(
		seq_len(nrow(batches)),
		function(i) {
			b = batches$batch[i]
			d = batches$dir[i]
			frozen = read_freeze(file.path(d, freeze_files[[b]]))
			current = function(name) digest(file = file.path(sandbox_root, d, name), algo = "sha256")
			matches = function(name) if (name %in% frozen$file) frozen[file == name, sha256] == current(name) else NA
			readers = assignments[batch == b]
			data.table(
				batch = b,
				families = batches$families[i],
				first_reader = paste(unique(readers$first_reader), collapse = "|"),
				independent_reader = paste(unique(readers$independent_reader), collapse = "|"),
				freeze_record = file.path(d, freeze_files[[b]]),
				freeze_time_utc = frozen$time_utc[1L],
				independent_units_match_freeze = matches("independent_units.csv"),
				independent_comparisons_match_freeze = matches("independent_comparisons.csv"),
				findings_file = if (file.exists(file.path(sandbox_root, d, "review_findings.txt"))) file.path(d, "review_findings.txt") else file.path(d, "accepted", "review_findings.txt")
			)
		}
	)
)
independence = merge(independence, findings_notes, by = "batch", sort = FALSE)
independence = merge(independence, first_freeze_summary, by = "batch", sort = FALSE)
independence[, classification := fifelse(
	independent_units_match_freeze & independent_comparisons_match_freeze & first_labels_withheld_until_freeze,
	"independent of the first reading's labels; context not fresh",
	"see notes"
)]
setcolorder(
	independence,
	c(
		"batch",
		"families",
		"first_reader",
		"independent_reader",
		"freeze_record",
		"freeze_time_utc",
		"independent_units_match_freeze",
		"independent_comparisons_match_freeze",
		names(first_freeze_summary)[-1L]
	)
)
write_table = function(x, name) {
	fwrite(x, file.path(output_dir, name), na = "")
}

write_table(batch_agreement, "batch_agreement.csv")
write_table(outcome_agreement, "outcome_agreement.csv")
write_table(reconciliation_decisions, "reconciliation_decisions.csv")
write_table(continuity_bound, "continuity_bound.csv")
write_table(independence, "batch_independence.csv")
write_table(first_reading_freeze, "first_reading_freeze.csv")
write_table(
	requests[, .(
		batch,
		family_id,
		pair_id,
		prior_year,
		later_year,
		common,
		accepted_id,
		root_id,
		paragraph,
		accepted_outcome,
		first_identified,
		first_route,
		first_outcome,
		independent_identified,
		independent_route,
		independent_outcome,
		both_identified,
		outcomes_compared,
		readings_agree,
		unit_origin,
		independent_presence,
		outcome_adoption
	)],
	"request_alignment.csv"
)
write_table(
	rbind(first_only, independent_only, fill = TRUE)[
		matched == FALSE,
		.(reading, batch, pair_id, prior_id, root_id, paragraph, outcome_raw, outcome)
	],
	"reading_only_units.csv"
)
output_files = c(
	"batch_agreement.csv",
	"outcome_agreement.csv",
	"reconciliation_decisions.csv",
	"continuity_bound.csv",
	"batch_independence.csv",
	"first_reading_freeze.csv",
	"request_alignment.csv",
	"reading_only_units.csv"
)
members_used = manifest[sub("^[.]/", "", path) %in% unique(tracker$members_read), .(path, sha256)]
write_json(
	list(
		script = "code/08_reading_agreement.R",
		run_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
		R = R.version.string,
		archive = list(path = archive, sha256 = digest(file = archive, algo = "sha256")),
		manifest = list(
			path = manifest_path,
			sha256 = digest(file = manifest_path, algo = "sha256"),
			members_verified = nrow(manifest)
		),
		members_used = members_used,
		batches = batch_names,
		headline_check = headline_check,
		rules = list(
			accepted_core_request = "count_status in {included, count} and category in {substantive, substantive_request}, as in code/07_proposal_comparison.R",
			first_reading_request_status = first_core_status,
			independent_reading_request_status = independent_core_status,
			identification = "a reading identified an accepted request when one of its units counted as a request is linked to that request by an identifier crosswalk or anchored to the same adopted text and paragraph, and has an outcome in the same year pair",
			outcome = "the outcome of the linked units; several distinct outcomes are resolved by the identifier link or the reconciliation record, otherwise the outcome is ambiguous and excluded from the agreement statistics",
			outcome_harmonization = as.list(alias),
			noncore_labels = noncore_labels,
			addition_labels = addition_labels,
			kappa = "Cohen's kappa over the seven outcome categories; kappa_continued over continued_without_substantive_change against the rest",
			disagreement_bound = "accepted continued_without_substantive_change counts as continuing unless both readings assigned outcomes and they differ",
			confirmed_bound = "a request counts as continuing only when the accepted coding and both readings assign continued_without_substantive_change",
			hand_coded_columns = "first_labels_withheld_until_freeze, fresh_context and context_reused in batch_independence.csv are coded from each batch's review_findings.txt",
			first_reading_freeze = paste(
				"each first-reading file the script reads is compared with the hash in the batch's",
				"first-reader freeze record; files with no record are listed in first_reading_freeze.csv"
			)
		),
		outputs = lapply(
			setNames(output_files, output_files),
			function(f) digest(file = file.path(output_dir, f), algo = "sha256")
		)
	),
	file.path(output_dir, "provenance.json"),
	pretty = TRUE,
	auto_unbox = TRUE
)
unlink(sandbox_root, recursive = TRUE)
