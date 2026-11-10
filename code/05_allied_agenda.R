# compare allied agreement across recorded agendas

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
vote_path = file.path(root, "data", "un_votes_long.rda")
match_path = file.path(
	root,
	"results",
	"adoptions",
	"recent",
	"harmonized_matched_adoption_methods.csv"
)
output_root = file.path(root, "results", "observed", "allied_agenda")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
vote_env = new.env(parent = emptyenv())
load_vote_archive(vote_path, envir = vote_env)
votes = as.data.table(vote_env$un_votes_long)[year %in% 2024:2025 & recorded_vote %in% c("Y", "A", "N")]
votes[, position := c(Y = 1, A = 0, N = -1)[recorded_vote]]
actors = sort(unique(votes$voting_country_code))
# yes, abstain and no as 1, 0 and -1; nonvotes stay missing
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

matched = read_analysis_table(match_path)
resolutions_2025 = unique(votes[year == 2025, resolution])
groups = list(
	newly_recorded = matched[method_2024 == "without_vote" & method_2025 == "with_vote", resolution_2025],
	recurring_recorded = matched[method_2024 == "with_vote" & method_2025 == "with_vote", resolution_2025]
)
stopifnot(
	length(groups$newly_recorded) == 33L,
	length(groups$recurring_recorded) == 58L,
	all(unlist(groups) %in% resolutions_2025)
)
groups$other_recorded = setdiff(resolutions_2025, unlist(groups))
full_2024 = agreement_network(
	position_matrix(2024, unique(votes[year == 2024, resolution]))
)
full_2025 = agreement_network(position_matrix(2025, resolutions_2025))
parts = lapply(groups, function(ids) agreement_network(position_matrix(2025, ids)))

others = setdiff(actors, c("USA", allies))
pairs = rbind(
	data.table(
		a = t(combn(allies, 2))[, 1],
		b = t(combn(allies, 2))[, 2],
		pair_group = "within_allies"
	),
	CJ(a = allies, b = others)[, pair_group := "allies_other"]
)
index = cbind(match(pairs$a, actors), match(pairs$b, actors))
pairs[, `:=`(
	agreement_2024 = full_2024$y[index],
	agreement_2025 = full_2025$y[index],
	support_2025 = full_2025$support[index]
)]
pairs = pairs[is.finite(agreement_2024) & is.finite(agreement_2025)]
index = cbind(match(pairs$a, actors), match(pairs$b, actors))

for (group_name in names(groups)) {
	y = parts[[group_name]]$y[index]
	s = parts[[group_name]]$support[index]
	pairs[, (paste0("agreement_", group_name)) := y]
	pairs[, (paste0("weight_", group_name)) := s / support_2025]
	pairs[, (paste0("component_", group_name)) := fifelse(s > 0, y * s / support_2025, 0)]
}

pairs[, identity_error := abs(
	agreement_2025 - component_newly_recorded - component_recurring_recorded - component_other_recorded
)]
stopifnot(max(pairs$identity_error) < 1e-12)
summary_rows = list()
for (pair_name in c("within_allies", "allies_other")) {
	d = pairs[pair_group == pair_name]
	row = data.table(
		pair_group = pair_name,
		dyads = nrow(d),
		full_2024 = mean(d$agreement_2024),
		full_2025 = mean(d$agreement_2025)
	)
	for (group_name in names(groups)) {
		row[, (paste0(group_name, "_resolutions")) := length(groups[[group_name]])]
		row[, (paste0(group_name, "_agreement")) := mean(d[[paste0("agreement_", group_name)]], na.rm = TRUE)]
		row[, (paste0(group_name, "_mean_weight")) := mean(d[[paste0("weight_", group_name)]])]
		row[, (paste0(group_name, "_component")) := mean(d[[paste0("component_", group_name)]])]
	}
	summary_rows[[pair_name]] = row
}

summary_table = rbindlist(summary_rows)
# the 2025 gap without the newly recorded resolutions, from the same dyads
without_new = agreement_network(
	position_matrix(2025, c(groups$recurring_recorded, groups$other_recorded))
)
pairs[, agreement_without_newly_recorded := without_new$y[index]]
level = function(column, pair_name) pairs[pair_group == pair_name, mean(get(column), na.rm = TRUE)]
columns = c(
	full_2024 = "agreement_2024",
	full_2025 = "agreement_2025",
	newly_recorded = "agreement_newly_recorded",
	recurring_recorded = "agreement_recurring_recorded",
	other_recorded = "agreement_other_recorded",
	without_newly_recorded = "agreement_without_newly_recorded"
)
gap_table = data.table(
	measure = names(columns),
	resolutions = c(
		length(unique(votes[year == 2024, resolution])),
		length(resolutions_2025),
		lengths(groups),
		length(groups$recurring_recorded) + length(groups$other_recorded)
	),
	us_no = c(
		votes[year == 2024 & voting_country_code == "USA", sum(recorded_vote == "N")],
		votes[year == 2025 & voting_country_code == "USA", sum(recorded_vote == "N")],
		vapply(
			groups,
			function(ids) votes[year == 2025 & voting_country_code == "USA" & resolution %in% ids, sum(recorded_vote == "N")],
			integer(1)
		),
		votes[year == 2025 & voting_country_code == "USA" & !resolution %in% groups$newly_recorded, sum(recorded_vote == "N")]
	),
	within_allies = vapply(columns, level, numeric(1), pair_name = "within_allies"),
	allies_other = vapply(columns, level, numeric(1), pair_name = "allies_other")
)
gap_table[, gap := within_allies - allies_other]
# vote shares in each resolution group
group_votes = votes[year == 2025]
group_votes[, resolution_group := names(groups)[vapply(
	resolution,
	function(r) which(vapply(groups, function(g) r %in% g, logical(1))),
	integer(1)
)]]
group_votes[, actor_group := fifelse(
	voting_country_code == "USA",
	"USA",
	fifelse(voting_country_code %in% allies, "allies", "other_members")
)]
vote_table = group_votes[, .(
	resolutions = uniqueN(resolution),
	votes = .N,
	yes = sum(recorded_vote == "Y"),
	abstain = sum(recorded_vote == "A"),
	no = sum(recorded_vote == "N"),
	yes_share = mean(recorded_vote == "Y"),
	no_share = mean(recorded_vote == "N")
),
by = .(resolution_group, actor_group)
][order(resolution_group, actor_group)]
fwrite(pairs, file.path(output_root, "dyad_agreement_by_group.csv"))
fwrite(
	summary_table,
	file.path(output_root, "group_agreement_by_resolution_group.csv")
)
fwrite(
	gap_table,
	file.path(output_root, "allied_gap_by_resolution_group.csv")
)
fwrite(
	vote_table,
	file.path(output_root, "vote_shares_by_resolution_group.csv")
)
write_json(
	list(
		script = "code/05_allied_agenda.R",
		inputs = list(
			votes = list(
				path = "data/un_votes_long.rda",
				sha256 = digest(file = vote_path, algo = "sha256")
			),
			matched = list(
				path = "results/adoptions/recent/harmonized_matched_adoption_methods.csv",
				sha256 = digest(file = match_path, algo = "sha256")
			)
		),
		resolution_groups = lapply(groups, length),
		allies = allies,
		R = R.version.string
	),
	file.path(output_root, "provenance.json"),
	pretty = TRUE,
	auto_unbox = TRUE
)
