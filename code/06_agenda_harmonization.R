# match resolution titles and predecessors

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
set.seed(20261031)
allies = actor_groups$allies
window_years = 2022:2024
full_years = 2010:2024
minimum_part_tokens = 4L
inventory_path = file.path(
	root,
	"results",
	"adoptions",
	"history",
	"archived_vote_join_audit.csv"
)
strict_path = file.path(
	root,
	"results",
	"adoptions",
	"recent",
	"harmonized_matched_adoption_methods.csv"
)
recorded_path = file.path(
	root,
	"results",
	"adoptions",
	"recent",
	"recorded_total_validation.csv"
)
vote_path = file.path(root, "data", "un_votes_long.rda")
agenda_root = file.path(root, "results", "observed", "allied_agenda")
output_root = file.path(root, "results", "adoptions", "harmonized")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
# ordinal words for one to ninety-nine, compound forms written without a hyphen
ordinal_words = function() {
	units = c(
		"first",
		"second",
		"third",
		"fourth",
		"fifth",
		"sixth",
		"seventh",
		"eighth",
		"ninth"
	)
	teens = c(
		"tenth",
		"eleventh",
		"twelfth",
		"thirteenth",
		"fourteenth",
		"fifteenth",
		"sixteenth",
		"seventeenth",
		"eighteenth",
		"nineteenth"
	)
	tens = c(
		"twenty",
		"thirty",
		"forty",
		"fifty",
		"sixty",
		"seventy",
		"eighty",
		"ninety"
	)
	tens_ordinal = c(
		"twentieth",
		"thirtieth",
		"fortieth",
		"fiftieth",
		"sixtieth",
		"seventieth",
		"eightieth",
		"ninetieth"
	)
	words = c(units, teens)
	for (i in seq_along(tens)) {
		words = c(words, tens_ordinal[i], paste0(tens[i], units))
	}
	words
}

ordinals = ordinal_words()

# the colon is kept
clean_title = function(title) {
	out = normalize_title(title)
	out = gsub("[‐‑‒–—]", "-", out)
	out = gsub("[‘’“”]", "'", out)
	out = gsub("\u00a0", " ", out, fixed = TRUE)
	out = sub(
		"^[0-9]{1,3}[[:space:]]*(\\([a-z]+\\))?[[:space:]]*-?[[:space:]]+(?=[a-z])",
		"",
		out,
		perl = TRUE
	)
	for (n in seq_along(ordinals)) {
		out = gsub(paste0("\\b", n, "(st|nd|rd|th)\\b"), ordinals[n], out, perl = TRUE)
	}
	compound = paste0(
		"\\b(twenty|thirty|forty|fifty|sixty|seventy|eighty|ninety)[- ]+",
		"(first|second|third|fourth|fifth|sixth|seventh|eighth|ninth)\\b"
	)
	out = gsub(compound, "\\1\\2", out, perl = TRUE)
	out = gsub("[0-9]+", " ", out)
	trimws(gsub("[[:space:]]+", " ", out))
}

key_tokens = function(x) {
	out = gsub("(?<=[a-z])['-](?=[a-z])", "", x, perl = TRUE)
	out = gsub("[^a-z]", " ", out)
	out = gsub("\\b(a|an|the|st|nd|rd|th)\\b", " ", out, perl = TRUE)
	trimws(gsub("[[:space:]]+", " ", out))
}

token_count = function(key) lengths(strsplit(key, " ", fixed = TRUE)) * nzchar(key)

title_keys = function(title) {
	cleaned = clean_title(title)
	has_colon = grepl(":", cleaned, fixed = TRUE)
	before = ifelse(has_colon, sub(":.*$", "", cleaned), "")
	after = ifelse(has_colon, sub("^[^:]*:", "", cleaned), "")
	keys = data.table(
		key_full = key_tokens(cleaned),
		key_before = key_tokens(before),
		key_after = key_tokens(after)
	)
	keys[, `:=`(
		tokens_full = token_count(key_full),
		tokens_before = token_count(key_before),
		tokens_after = token_count(key_after)
	)]
	keys
}

# adoption inventory from the official session lists, 2010 to 2025
inventory = read_analysis_table(inventory_path)
inventory = inventory[!is.na(method) & nzchar(method)]
inventory[, title := fifelse(has_archived_votes == TRUE, archived_title, native_title)]
inventory[, title_source := fifelse(
	has_archived_votes == TRUE,
	"archived_un_voting_dataset",
	"official_session_list"
)]
inventory[, date := as.Date(date)]
inventory[, adopted_without_vote := method %in% c("without_vote", "silence_procedure")]
inventory = cbind(inventory, title_keys(inventory$title))
inventory[, key_multiplicity := .N, by = c("year", "key_full")]
title_key_columns = c(
	"year",
	"date",
	"canonical_resolution",
	"archived_resolution",
	"method",
	"title_source",
	"title",
	"key_full",
	"key_before",
	"key_after",
	"key_multiplicity"
)
fwrite(
	inventory[, title_key_columns, with = FALSE],
	file.path(output_root, "title_keys.csv")
)
targets = inventory[year == 2025 & method == "with_vote" & regular_resolution == TRUE]
recorded_2025 = read_analysis_table(recorded_path)[year == 2025, resolution]
stopifnot(
	nrow(targets) == 187L,
	!anyNA(targets$archived_resolution),
	setequal(targets$archived_resolution, recorded_2025)
)

strict = read_analysis_table(strict_path)
strict_label = strict[, .(
	archived_resolution = resolution_2025,
	strict_resolution_2024 = sub(" ([A-Z])$", "\\1", resolution_2024),
	strict_method_2024 = method_2024
)]
targets = merge(targets, strict_label, by = "archived_resolution", all.x = TRUE)
targets[, strict_group := fcase(
	strict_method_2024 == "without_vote",
	"newly_recorded",
	strict_method_2024 == "with_vote",
	"recurring_recorded",
	default = "unmatched"
)]
stopifnot(
	targets[, sum(strict_group == "newly_recorded")] == 33L,
	targets[, sum(strict_group == "recurring_recorded")] == 58L
)

candidates = inventory[year < 2025]
part_match = function(target, pool) {
	enough = minimum_part_tokens
	full_usable = nzchar(target$key_full) & target$tokens_full >= enough
	before_usable = nzchar(target$key_before) & target$tokens_before >= enough
	pool[
		(full_usable & key_after == target$key_full) |
			(before_usable & key_before == target$key_before) |
			(before_usable & key_full == target$key_before) |
			(full_usable & key_before == target$key_full)
	]
}

find_predecessor = function(target, years) {
	for (year_value in sort(years, decreasing = TRUE)) {
		pool = candidates[year == year_value]
		hit = pool[key_full == target$key_full]
		step = "harmonized_title"
		if (nrow(hit) == 0L) {
			hit = part_match(target, pool)
			step = "title_part"
		}
		if (nrow(hit) == 0L) {
			next
		}
		hit = hit[date == max(date)]
		method_value = if (uniqueN(hit$adopted_without_vote) == 1L) hit$method[1] else "ambiguous"
		return(
			data.table(
				predecessor = hit$canonical_resolution[1],
				predecessor_year = year_value,
				predecessor_date = hit$date[1],
				predecessor_method = method_value,
				predecessor_title = hit$title[1],
				predecessor_candidates = nrow(hit),
				match_step = step,
				lag_years = 2025L - year_value
			)
		)
	}
	data.table(
		predecessor = NA_character_,
		predecessor_year = NA_integer_,
		predecessor_date = as.Date(NA),
		predecessor_method = "none",
		predecessor_title = NA_character_,
		predecessor_candidates = 0L,
		match_step = "none",
		lag_years = NA_integer_
	)
}

predecessor_table = function(years, prefix) {
	out = rbindlist(
		lapply(
			seq_len(nrow(targets)),
			function(i) find_predecessor(targets[i], years)
		)
	)
	setnames(out, paste0(prefix, names(out)))
	out
}

targets = cbind(
	targets,
	predecessor_table(window_years, "window_"),
	predecessor_table(full_years, "full_")
)

jaccard = function(a, b) {
	a = unique(a)
	b = unique(b)
	length(intersect(a, b)) / length(union(a, b))
}

candidate_tokens = strsplit(candidates$key_full, " ", fixed = TRUE)
closest_title = function(i) {
	target_tokens = strsplit(targets$key_full[i], " ", fixed = TRUE)[[1]]
	similarity = vapply(candidate_tokens, jaccard, numeric(1), b = target_tokens)
	ordering = order(-similarity, -candidates$year)
	best = candidates[ordering[1]]
	data.table(
		closest_resolution = best$canonical_resolution,
		closest_year = best$year,
		closest_method = best$method,
		closest_similarity = similarity[ordering[1]],
		closest_title = best$title
	)
}

closest = rbindlist(lapply(seq_len(nrow(targets)), closest_title))
targets = cbind(targets, closest)

reviewed = data.table(
	archived_resolution = c("A/RES/80/53", "A/RES/80/136", "A/RES/80/111"),
	reviewed_predecessor = c("A/RES/79/240", "A/RES/79/203", "A/RES/77/288"),
	reviewed_reason = c(
		"the 2025 title adds the prefix Establishment of the to the 2024 title",
		"the annual follow-up item on small island developing states changes from the SAMOA Pathway to the Antigua and Barbuda Agenda",
		"the agenda item on the Chernobyl disaster is considered every third year; the 2016, 2019 and 2023 resolutions are titled Persistent legacy of the Chernobyl disaster"
	)
)
stopifnot(
	all(reviewed$archived_resolution %in% targets$archived_resolution),
	all(reviewed$reviewed_predecessor %in% candidates$canonical_resolution),
	all(
		targets[archived_resolution %in% reviewed$archived_resolution, window_predecessor_method] == "none"
	)
)
reviewed_candidates = candidates[, .(
	reviewed_predecessor = canonical_resolution,
	reviewed_year = year,
	reviewed_method = method,
	reviewed_title = title
)]
reviewed = merge(reviewed, reviewed_candidates, by = "reviewed_predecessor")
stopifnot(all(reviewed$reviewed_year %in% window_years))
targets = merge(targets, reviewed, by = "archived_resolution", all.x = TRUE)
# group labels under each rule set
group_of = function(method) {
	fcase(
		method == "none",
		"unmatched_regular",
		method == "ambiguous",
		"ambiguous_regular",
		method %in% c("without_vote", "silence_procedure"),
		"newly_recorded",
		method == "with_vote",
		"recurring_recorded"
	)
}

targets[, harmonized_group := group_of(window_predecessor_method)]
targets[, full_window_group := group_of(full_predecessor_method)]
targets[, reviewed_group := fifelse(
	!is.na(reviewed_predecessor),
	group_of(reviewed_method),
	harmonized_group
)]
targets[, reviewed_lag_years := fifelse(!is.na(reviewed_predecessor), 2025L - reviewed_year, window_lag_years)]
# a reading of the titles that the window rule leaves unmatched
unmatched_reading = data.table(
	archived_resolution = c(
		"A/RES/79/266",
		"A/RES/79/269",
		"A/RES/79/270",
		"A/RES/79/291",
		"A/RES/79/293",
		"A/RES/80/10",
		"A/RES/80/106",
		"A/RES/79/283",
		"A/RES/80/132",
		"A/RES/79/323",
		"A/RES/80/117",
		"A/RES/79/278",
		"A/RES/80/56",
		"A/RES/79/315",
		"A/RES/80/1",
		"A/RES/80/23",
		"A/RES/80/116",
		"A/RES/80/119",
		"A/RES/80/152",
		"A/RES/79/314",
		"A/RES/79/329",
		"A/RES/80/53",
		"A/RES/80/111",
		"A/RES/80/136"
	),
	reading = c(
		rep("new international day", 7L),
		rep("new United Nations decade", 2L),
		rep("new declaration or outcome document", 2L),
		rep("modalities of a new meeting or conference", 2L),
		rep("other new resolution", 6L),
		rep("recurring item with a predecessor before 2022", 2L),
		rep("renamed recurring item, reviewed link", 3L)
	),
	reading_note = c(
		rep("", 7L),
		rep("", 2L),
		"",
		"political declarations of the earlier high-level meetings on noncommunicable diseases were adopted without a vote in 2011 and 2018",
		"",
		"follows the 2023 and 2024 recorded resolutions on victim assistance and environmental remediation",
		"",
		"",
		"",
		"follows the 2024 recorded request for the advisory opinion",
		"",
		"",
		"the 2017 declaration of the first United Nations Ocean Conference was adopted without a vote",
		"the most recent adoption under this title in the lists is in 2013, without a vote",
		rep("see the reviewed link", 3L)
	)
)
unmatched_targets = targets[harmonized_group == "unmatched_regular", archived_resolution]
stopifnot(setequal(unmatched_reading$archived_resolution, unmatched_targets))
targets = merge(targets, unmatched_reading, by = "archived_resolution", all.x = TRUE)
vote_env = new.env(parent = emptyenv())
load_vote_archive(vote_path, envir = vote_env)
votes = as.data.table(vote_env$un_votes_long)[year %in% 2024:2025 & recorded_vote %in% c("Y", "A", "N")]
votes[, position := c(Y = 1, A = 0, N = -1)[recorded_vote]]
actors = sort(unique(votes$voting_country_code))
resolutions_2024 = unique(votes[year == 2024, resolution])
resolutions_2025 = unique(votes[year == 2025, resolution])
outside_regular = setdiff(resolutions_2025, targets$archived_resolution)
stopifnot(length(resolutions_2025) == 192L, length(outside_regular) == 5L)
us_vote = votes[year == 2025 & voting_country_code == "USA", .(archived_resolution = resolution, us_vote_2025 = recorded_vote)]
targets = merge(targets, us_vote, by = "archived_resolution", all.x = TRUE)
predecessor_columns = c(
	"archived_resolution",
	"canonical_resolution",
	"date",
	"title",
	"key_full",
	"us_vote_2025",
	"strict_group",
	"strict_resolution_2024",
	"strict_method_2024",
	"harmonized_group",
	"window_predecessor",
	"window_predecessor_year",
	"window_predecessor_date",
	"window_predecessor_method",
	"window_predecessor_title",
	"window_predecessor_candidates",
	"window_match_step",
	"window_lag_years",
	"full_window_group",
	"full_predecessor",
	"full_predecessor_year",
	"full_predecessor_method",
	"full_predecessor_title",
	"full_match_step",
	"full_lag_years",
	"reviewed_group",
	"reviewed_predecessor",
	"reviewed_year",
	"reviewed_method",
	"reviewed_title",
	"reviewed_reason",
	"reviewed_lag_years",
	"reading",
	"reading_note",
	"closest_resolution",
	"closest_year",
	"closest_method",
	"closest_similarity",
	"closest_title"
)
setorder(targets, date, archived_resolution)
fwrite(
	targets[, predecessor_columns, with = FALSE],
	file.path(output_root, "predecessors_2025.csv")
)
# counts of the 187 by rule set, group and lag, with the american vote
rule_sets = c(
	strict = "strict_group",
	harmonized = "harmonized_group",
	harmonized_reviewed = "reviewed_group",
	harmonized_full_window = "full_window_group"
)
lag_columns = c(
	strict = NA_character_,
	harmonized = "window_lag_years",
	harmonized_reviewed = "reviewed_lag_years",
	harmonized_full_window = "full_lag_years"
)
count_rows = list()
for (rule_name in names(rule_sets)) {
	group_column = rule_sets[[rule_name]]
	lag_column = lag_columns[[rule_name]]
	d = copy(targets)
	d[, group := get(group_column)]
	d[, lag := if (is.na(lag_column)) fifelse(group == "unmatched", NA_integer_, 1L) else get(lag_column)]
	count_rows[[rule_name]] = d[, .(
		rule_set = rule_name,
		resolutions = .N,
		us_no = sum(us_vote_2025 == "N", na.rm = TRUE),
		us_yes = sum(us_vote_2025 == "Y", na.rm = TRUE),
		us_abstain = sum(us_vote_2025 == "A", na.rm = TRUE),
		us_nonvote = sum(is.na(us_vote_2025))
	), by = c("group", "lag")]
}

predecessor_counts = rbindlist(count_rows)[order(rule_set, group, lag)]
group_totals = predecessor_counts[, .(
	resolutions = sum(resolutions),
	us_no = sum(us_no),
	us_yes = sum(us_yes),
	us_abstain = sum(us_abstain),
	us_nonvote = sum(us_nonvote)
), by = c("rule_set", "group")]
fwrite(
	predecessor_counts,
	file.path(output_root, "predecessor_counts_by_lag.csv")
)
fwrite(group_totals, file.path(output_root, "predecessor_counts.csv"))
# resolution groups over the 192 archive items for each rule set
group_lists = lapply(
	rule_sets,
	function(group_column) {
		d = targets[, .(archived_resolution, group = get(group_column))]
		out = split(d$archived_resolution, d$group)
		if (group_column == "strict_group") {
			out = list(
				newly_recorded = out$newly_recorded,
				recurring_recorded = out$recurring_recorded,
				other_recorded = c(out$unmatched, outside_regular)
			)
		} else {
			out$outside_regular = outside_regular
		}
		out
	}
)

position_matrix = function(year_value, ids) {
	rows = votes[year == year_value & resolution %in% ids]
	out = matrix(NA_real_, length(actors), length(ids), dimnames = list(actors, ids))
	out[cbind(match(rows$voting_country_code, actors), match(rows$resolution, ids))] = rows$position
	out
}

agreement_network = function(position) {
	m = 1 * is.finite(position)
	p = position
	p[!is.finite(position)] = 0
	a = 1 * (is.finite(position) & position == 0)
	support = tcrossprod(m)
	y = 0.5 * (support + tcrossprod(p) + tcrossprod(a)) / support
	y[support == 0] = NA_real_
	list(y = y, support = support)
}

full_2024 = agreement_network(position_matrix(2024, resolutions_2024))
full_2025 = agreement_network(position_matrix(2025, resolutions_2025))
others = setdiff(actors, c("USA", allies))
ally_pairs = t(combn(allies, 2))
pairs = rbind(
	data.table(
		a = ally_pairs[, 1],
		b = ally_pairs[, 2],
		pair_group = "within_allies"
	),
	CJ(a = allies, b = others)[, pair_group := "allies_other"]
)
index = cbind(match(pairs$a, actors), match(pairs$b, actors))
pairs[, `:=`(agreement_2024 = full_2024$y[index], agreement_2025 = full_2025$y[index])]
pairs = pairs[is.finite(agreement_2024) & is.finite(agreement_2025)]
index = cbind(match(pairs$a, actors), match(pairs$b, actors))
stopifnot(
	pairs[, sum(pair_group == "within_allies")] == 91L,
	pairs[, sum(pair_group == "allies_other")] == 2464L
)
group_level = function(ids, pair_name) {
	y = agreement_network(position_matrix(2025, ids))$y[index]
	mean(y[pairs$pair_group == pair_name], na.rm = TRUE)
}

us_no_count = function(ids) votes[year == 2025 & voting_country_code == "USA" & resolution %in% ids, sum(recorded_vote == "N")]
us_no_2024 = votes[year == 2024 & voting_country_code == "USA", sum(recorded_vote == "N")]
gap_rows = share_rows = list()
for (rule_name in names(group_lists)) {
	groups = group_lists[[rule_name]]
	measures = c(
		list(full_2024 = resolutions_2024, full_2025 = resolutions_2025),
		groups,
		list(
			without_newly_recorded = setdiff(resolutions_2025, groups$newly_recorded)
		)
	)
	gap_rows[[rule_name]] = data.table(
		rule_set = rule_name,
		measure = names(measures),
		resolutions = lengths(measures),
		us_no = c(us_no_2024, vapply(measures[-1], us_no_count, integer(1))),
		within_allies = c(
			mean(pairs[pair_group == "within_allies", agreement_2024]),
			vapply(measures[-1], group_level, numeric(1), pair_name = "within_allies")
		),
		allies_other = c(
			mean(pairs[pair_group == "allies_other", agreement_2024]),
			vapply(measures[-1], group_level, numeric(1), pair_name = "allies_other")
		)
	)
	group_of_resolution = function(r) which(vapply(groups, function(g) r %in% g, logical(1)))
	group_votes = votes[year == 2025]
	group_votes[, resolution_group := names(groups)[vapply(resolution, group_of_resolution, integer(1))]]
	group_votes[, actor_group := fifelse(
		voting_country_code == "USA",
		"USA",
		fifelse(voting_country_code %in% allies, "allies", "other_members")
	)]
	share_rows[[rule_name]] = group_votes[, .(
		rule_set = rule_name,
		resolutions = uniqueN(resolution),
		votes = .N,
		yes = sum(recorded_vote == "Y"),
		abstain = sum(recorded_vote == "A"),
		no = sum(recorded_vote == "N"),
		yes_share = mean(recorded_vote == "Y"),
		no_share = mean(recorded_vote == "N")
	), by = .(resolution_group, actor_group)][order(resolution_group, actor_group)]
}

gap_table = rbindlist(gap_rows)
gap_table[, gap := within_allies - allies_other]
share_table = rbindlist(share_rows)
setcolorder(share_table, c("rule_set", "resolution_group", "actor_group"))
fwrite(gap_table, file.path(output_root, "allied_gap_by_rule_set.csv"))
fwrite(share_table, file.path(output_root, "vote_shares_by_rule_set.csv"))

fuller = targets[harmonized_group == "newly_recorded"]
fuller_votes = votes[year == 2025 & resolution %in% fuller$archived_resolution]
fuller_summary = fuller_votes[, .(
	us_vote = if (any(voting_country_code == "USA")) recorded_vote[voting_country_code == "USA"] else NA_character_,
	member_yes = sum(recorded_vote == "Y"),
	member_no = sum(recorded_vote == "N"),
	member_abstain = sum(recorded_vote == "A"),
	ally_yes = sum(voting_country_code %in% allies & recorded_vote == "Y"),
	ally_no = sum(voting_country_code %in% allies & recorded_vote == "N"),
	ally_abstain = sum(voting_country_code %in% allies & recorded_vote == "A"),
	ally_no_codes = paste(
		sort(
			voting_country_code[voting_country_code %in% allies & recorded_vote == "N"]
		),
		collapse = ";"
	),
	ally_abstain_codes = paste(
		sort(
			voting_country_code[voting_country_code %in% allies & recorded_vote == "A"]
		),
		collapse = ";"
	)
), by = .(archived_resolution = resolution)]
fuller_columns = c(
	"archived_resolution",
	"date",
	"title",
	"strict_group",
	"window_predecessor",
	"window_predecessor_year",
	"window_predecessor_method",
	"window_predecessor_title",
	"window_match_step",
	"window_lag_years"
)
fuller_summary = merge(
	fuller[, fuller_columns, with = FALSE],
	fuller_summary,
	by = "archived_resolution"
)
setorder(fuller_summary, date, archived_resolution)
fwrite(fuller_summary, file.path(output_root, "newly_recorded_summary.csv"))

residual_columns = c(
	"archived_resolution",
	"date",
	"title",
	"us_vote_2025",
	"harmonized_group",
	"window_predecessor",
	"window_predecessor_year",
	"window_predecessor_method",
	"window_predecessor_title",
	"window_match_step",
	"window_lag_years",
	"reviewed_group",
	"reviewed_predecessor",
	"reviewed_lag_years",
	"full_window_group",
	"full_predecessor",
	"full_predecessor_year",
	"full_predecessor_method",
	"reading",
	"reading_note",
	"closest_resolution",
	"closest_year",
	"closest_method",
	"closest_similarity",
	"closest_title"
)
residual = targets[strict_group == "unmatched", residual_columns, with = FALSE]
outside = votes[year == 2025 & resolution %in% outside_regular, .(
	archived_resolution = resolution[1],
	date = date[1],
	title = title[1],
	us_vote_2025 = if (any(voting_country_code == "USA")) recorded_vote[voting_country_code == "USA"] else NA_character_
), by = "resolution"][, resolution := NULL]
outside[, `:=`(
	harmonized_group = "outside_regular",
	reviewed_group = "outside_regular",
	full_window_group = "outside_regular"
)]
residual = rbind(residual, outside, fill = TRUE)
composition = function(group, lag) {
	fcase(
		group == "outside_regular",
		"outside the regular-session adoption lists",
		group == "newly_recorded",
		"recurring item newly put to a vote",
		group == "recurring_recorded" & lag == 1L,
		"recorded item with a 2024 predecessor",
		group == "recurring_recorded" & lag > 1L,
		"recorded item with a 2022 or 2023 predecessor",
		group == "ambiguous_regular",
		"ambiguous predecessor",
		default = "no predecessor found under the rule"
	)
}

residual[, composition_harmonized := composition(harmonized_group, window_lag_years)]
residual[, composition_reviewed := composition(reviewed_group, reviewed_lag_years)]
setorder(residual, date, archived_resolution)
fwrite(residual, file.path(output_root, "residual_items.csv"))
residual_counts = rbind(
	residual[, .(
		rule_set = "harmonized",
		items = .N,
		us_no = sum(us_vote_2025 == "N", na.rm = TRUE)
	), by = .(composition = composition_harmonized)],
	residual[, .(
		rule_set = "harmonized_reviewed",
		items = .N,
		us_no = sum(us_vote_2025 == "N", na.rm = TRUE)
	), by = .(composition = composition_reviewed)]
)
setcolorder(residual_counts, c("rule_set", "composition"))
fwrite(
	residual_counts[order(rule_set, -items)],
	file.path(output_root, "residual_counts.csv")
)
reading_counts = residual[harmonized_group == "unmatched_regular", .(
	items = .N,
	us_no = sum(us_vote_2025 == "N", na.rm = TRUE)
), by = "reading"]
fwrite(
	reading_counts[order(-items)],
	file.path(output_root, "unmatched_reading_counts.csv")
)

existing_gap = read_analysis_table(
	file.path(agenda_root, "allied_gap_by_resolution_group.csv")
)
existing_shares = read_analysis_table(
	file.path(agenda_root, "vote_shares_by_resolution_group.csv")
)
strict_gap = merge(
	existing_gap,
	gap_table[rule_set == "strict"],
	by = "measure",
	suffixes = c("_existing", "_new")
)
strict_shares = merge(
	existing_shares,
	share_table[rule_set == "strict"],
	by = c("resolution_group", "actor_group"),
	suffixes = c("_existing", "_new")
)
strict_matched = targets[strict_group != "unmatched"]
covers_archive = function(g) setequal(unlist(g), resolutions_2025) && !anyDuplicated(unlist(g))
checks = data.table(
	check = c(
		"targets_are_the_187_recorded_regular_resolutions",
		"strict_groups_reproduce_33_and_58",
		"residual_has_101_items",
		"strict_gap_reproduces_allied_agenda",
		"strict_shares_reproduce_allied_agenda",
		"no_predecessor_in_2025",
		"window_predecessors_inside_window",
		"reviewed_links_only_for_unmatched_items",
		"every_archive_item_in_one_group_per_rule_set",
		"strict_matches_keep_their_2024_predecessor"
	),
	passed = c(
		nrow(targets) == 187L,
		targets[, sum(strict_group == "newly_recorded") == 33L & sum(strict_group == "recurring_recorded") == 58L],
		nrow(residual) == 101L,
		nrow(strict_gap) == nrow(existing_gap) &&
			max(abs(strict_gap$gap_existing - strict_gap$gap_new)) < 1e-9 &&
			all(strict_gap$resolutions_existing == strict_gap$resolutions_new) &&
			all(strict_gap$us_no_existing == strict_gap$us_no_new),
		nrow(strict_shares) == nrow(existing_shares) &&
			all(strict_shares$yes_existing == strict_shares$yes_new) &&
			all(strict_shares$no_existing == strict_shares$no_new) &&
			all(strict_shares$votes_existing == strict_shares$votes_new),
		all(
			is.na(targets$window_predecessor_year) | targets$window_predecessor_year < 2025L
		),
		all(
			is.na(targets$window_predecessor_year) | targets$window_predecessor_year %in% window_years
		),
		all(
			targets[!is.na(reviewed_predecessor), window_predecessor_method == "none"]
		),
		all(vapply(group_lists, covers_archive, logical(1))),
		strict_matched[, all(
			window_predecessor == strict_resolution_2024 &
				window_predecessor_year == 2024L &
				harmonized_group == strict_group
		)]
	)
)
fwrite(checks, file.path(output_root, "validation.csv"))
stopifnot(all(checks$passed))
input_record = function(relative_path, path) list(path = relative_path, sha256 = digest(file = path, algo = "sha256"))
provenance = list(
	script = "code/06_agenda_harmonization.R",
	inputs = list(
		inventory = input_record(
			"results/adoptions/history/archived_vote_join_audit.csv",
			inventory_path
		),
		strict_matches = input_record(
			"results/adoptions/recent/harmonized_matched_adoption_methods.csv",
			strict_path
		),
		recorded = input_record(
			"results/adoptions/recent/recorded_total_validation.csv",
			recorded_path
		),
		votes = input_record("data/un_votes_long.rda", vote_path)
	),
	settings = list(
		window_years = window_years,
		full_window_years = full_years,
		minimum_part_tokens = minimum_part_tokens,
		silence_procedure_counted_as = "without_vote",
		title_rules = c(
			"archive suffix removed",
			"typographic dashes and quotes replaced",
			"non-breaking spaces replaced",
			"agenda-item number prefix removed",
			"numeral ordinals written as words, compound ordinals joined",
			"digits removed",
			"hyphenated and possessive forms joined",
			"articles and ordinal suffixes dropped",
			"parts before and after a colon matched when the full key fails and the part has at least four tokens"
		),
		reviewed_links = reviewed[, .(archived_resolution, reviewed_predecessor, reviewed_reason)]
	),
	resolution_groups = lapply(group_lists, function(g) lapply(g, length)),
	allies = allies,
	R = R.version.string
)
write_json(
	provenance,
	file.path(output_root, "provenance.json"),
	pretty = TRUE,
	auto_unbox = TRUE
)
