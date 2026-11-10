# compare national responses to american vote changes

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
set.seed(20261031)
argument = sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
root = dirname(dirname(normalizePath(argument)))
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
rivals = c("RUS", "CHN", "PRK", "IRN", "BLR", "SYR", "VEN", "CUB", "NIC")
plausible_followers = c("HUN", "ISR", "ARG", "PRY", "CZE", "TUR")
sparse_voters = c("BOL", "DMA", "MDG", "VCT", "STP", "VEN")
named_resolutions = c("A/RES/ES-10/19", "A/RES/ES-11/7", "A/RES/ES-11/8", "A/RES/79/323")
tax_pair = c("A/RES/79/235", "A/RES/80/128")
inputs = c(
	votes = "data/un_votes_long.rda",
	member_groups = "data/member_groups.csv",
	matched_titles = "results/observed/vote_changes/matched_titles.csv",
	us_vote_changes = "results/observed/vote_changes/us_vote_changes.csv",
	ally_responses = "results/observed/vote_changes/ally_responses.csv",
	country_peer_results = "results/observed/allies/country_peer_results.csv",
	reset_probabilities = "results/models/un/reset_probabilities.csv",
	us_votes = "results/adoptions/proposals/us_votes.csv",
	newly_recorded_votes = "results/adoptions/recent/harmonized_newly_recorded_2025_member_votes.csv",
	matched_methods = "results/adoptions/recent/harmonized_matched_adoption_methods.csv",
	vote_shares = "results/observed/allied_agenda/vote_shares_by_resolution_group.csv",
	recorded_totals = "results/adoptions/recent/recorded_total_validation.csv"
)
displacement_path = file.path(
	root,
	"results",
	"models",
	"un_variants",
	"threshold",
	"displacements.csv"
)
output_root = file.path(root, "results", "observed", "followers")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
read_input = function(name) read_analysis_table(file.path(root, inputs[[name]]))
checks = list()
add_check = function(name, value, expected) {
	assign(
		"checks",
		append(
			checks,
			list(
				data.table(
					check = name,
					value = paste(value, collapse = " "),
					expected = paste(expected, collapse = " "),
					passed = isTRUE(all.equal(as.numeric(value), as.numeric(expected)))
				)
			)
		),
		envir = parent.env(environment())
	)
}

# recorded votes with yes, abstain and no as 1, 0 and -1; nonvotes stay missing
vote_env = new.env(parent = emptyenv())
load_vote_archive(file.path(root, inputs[["votes"]]), envir = vote_env)
votes = as.data.table(vote_env$un_votes_long)[year %in% c(2017, 2024, 2025)]
votes[, position := unname(c(Y = 1, A = 0, N = -1)[recorded_vote])]
stopifnot(
	all(votes$recorded_vote %in% c("Y", "A", "N", "X")),
	anyDuplicated(votes[, .(year, resolution, voting_country_code)]) == 0L
)
members = sort(unique(votes[year == 2025, voting_country_code]))
stopifnot(setequal(members, unique(votes[year == 2024, voting_country_code])))

groups = read_input("member_groups")
stopifnot(setequal(groups$iso3, members), anyDuplicated(groups$iso3) == 0L)
groups[, treaty_group_2025 := fifelse(
	iso3 == "USA",
	"United States",
	fifelse(us_treaty_in_force_2025 %in% TRUE, us_treaty_group, "none")
)]
groups[, paper_group := fifelse(
	iso3 == "USA",
	"United States",
	fifelse(
		iso3 %in% allies,
		"allies",
		fifelse(iso3 %in% rivals, "rivals", "other members")
	)
)]
group_columns = groups[, .(
	iso3,
	member_name,
	paper_group,
	eu_member,
	treaty_group_2025,
	us_treaty,
	mnna
)]

matched = read_input("matched_titles")[year_baseline == 2024 & year_target == 2025]
us_changes = read_input("us_vote_changes")[baseline == 2024 & target == 2025]
add_check(
	"matched titles 2024 to 2025",
	nrow(matched),
	us_changes$matched_titles
)
us_pairs = merge(
	votes[year == 2024 & voting_country_code == "USA", .(
		resolution_2024 = resolution,
		us_vote_2024 = recorded_vote,
		us_position_2024 = position
	)],
	matched[, .(
		normalized_title,
		resolution_2024 = resolution_baseline,
		date_2024 = date_baseline,
		title_2024 = title_baseline,
		resolution_2025 = resolution_target,
		date_2025 = date_target,
		title_2025 = title_target
	)],
	by = "resolution_2024"
)
us_pairs = merge(
	us_pairs,
	votes[year == 2025 & voting_country_code == "USA", .(
		resolution_2025 = resolution,
		us_vote_2025 = recorded_vote,
		us_position_2025 = position,
		total_yes_2025 = total_yes,
		total_no_2025 = total_no,
		total_abstentions_2025 = total_abstentions
	)],
	by = "resolution_2025"
)
us_pairs = merge(
	us_pairs,
	unique(
		votes[year == 2024, .(
			resolution_2024 = resolution,
			total_yes_2024 = total_yes,
			total_no_2024 = total_no,
			total_abstentions_2024 = total_abstentions
		)]
	),
	by = "resolution_2024"
)
switches = us_pairs[is.finite(us_position_2024) & is.finite(us_position_2025) & us_position_2024 != us_position_2025]
setorder(switches, normalized_title)
add_check(
	"american switches 2024 to 2025",
	nrow(switches),
	us_changes$us_changes
)
add_check(
	"american switches to no",
	switches[us_vote_2025 == "N" & us_vote_2024 != "N", .N],
	us_changes$us_changes_to_no
)
switch_votes = rbindlist(
	lapply(
		seq_len(nrow(switches)),
		function(index) {
			before = votes[year == 2024 & resolution == switches$resolution_2024[index], .(
				iso3 = voting_country_code,
				vote_2024 = recorded_vote,
				position_2024 = position
			)]
			after = votes[year == 2025 & resolution == switches$resolution_2025[index], .(
				iso3 = voting_country_code,
				vote_2025 = recorded_vote,
				position_2025 = position
			)]
			out = merge(before, after, by = "iso3", all = TRUE)
			out[, `:=`(
				normalized_title = switches$normalized_title[index],
				resolution_2024 = switches$resolution_2024[index],
				resolution_2025 = switches$resolution_2025[index],
				us_vote_2024 = switches$us_vote_2024[index],
				us_vote_2025 = switches$us_vote_2025[index],
				us_position_2025 = switches$us_position_2025[index]
			)]
			out
		}
	)
)
switch_votes = switch_votes[iso3 != "USA"]
switch_votes[, paired := is.finite(position_2024) & is.finite(position_2025)]
switch_votes[, distance_2024 := abs(position_2024 - us_position_2025)]
switch_votes[, distance_2025 := abs(position_2025 - us_position_2025)]
switch_votes[, response := fcase(
	!paired,
	"not paired",
	position_2024 == position_2025,
	"unchanged",
	distance_2025 < distance_2024,
	"toward new US vote",
	distance_2025 > distance_2024,
	"away from new US vote",
	default = "changed, same distance"
)]
switch_votes[, already_at_new_us_vote := paired & position_2024 == us_position_2025]
switch_votes[, room_to_follow := paired & position_2024 != us_position_2025]
switch_votes = merge(switch_votes, group_columns, by = "iso3")
setcolorder(
	switch_votes,
	c(
		"normalized_title",
		"resolution_2024",
		"resolution_2025",
		"us_vote_2024",
		"us_vote_2025",
		"iso3",
		"member_name",
		"vote_2024",
		"vote_2025",
		"response"
	)
)
setorder(switch_votes, normalized_title, iso3)
member_responses = switch_votes[, .(
	switches = .N,
	paired = sum(paired),
	room_to_follow = sum(room_to_follow),
	toward = sum(response == "toward new US vote"),
	away = sum(response == "away from new US vote"),
	unchanged = sum(response == "unchanged"),
	changed_same_distance = sum(response == "changed, same distance"),
	not_paired = sum(!paired),
	already_at_new_us_vote = sum(already_at_new_us_vote),
	titles_toward = paste(normalized_title[response == "toward new US vote"], collapse = "; "),
	titles_away = paste(normalized_title[response == "away from new US vote"], collapse = "; ")
), by = iso3]
member_responses = merge(member_responses, group_columns, by = "iso3")
setorder(member_responses, -toward, -away, iso3)
ally_responses = read_input("ally_responses")[baseline == 2024 & target == 2025]
ally_check = merge(
	member_responses[, .(iso3, paired, unchanged, toward, away)],
	ally_responses[, .(
		iso3 = ally,
		us_switch_pairs,
		ally_unchanged_when_us_switches,
		closer_to_new_us_vote,
		further_from_new_us_vote
	)],
	by = "iso3"
)
add_check(
	"allied switch pairs agree with ally_responses.csv",
	ally_check$paired,
	ally_check$us_switch_pairs
)
add_check(
	"allied unchanged votes agree with ally_responses.csv",
	ally_check$unchanged,
	ally_check$ally_unchanged_when_us_switches
)
add_check(
	"allied moves toward agree with ally_responses.csv",
	ally_check$toward,
	ally_check$closer_to_new_us_vote
)
add_check(
	"allied moves away agree with ally_responses.csv",
	ally_check$away,
	ally_check$further_from_new_us_vote
)
summarize_responses = function(column) {
	member_responses[, .(
		grouping = column,
		members = .N,
		members_with_paired_votes = sum(paired > 0),
		switch_pairs = sum(paired),
		room_to_follow = sum(room_to_follow),
		toward = sum(toward),
		away = sum(away),
		unchanged = sum(unchanged),
		changed_same_distance = sum(changed_same_distance),
		not_paired = sum(not_paired),
		unchanged_share_of_paired = sum(unchanged) / sum(paired),
		toward_share_of_room = sum(toward) / sum(room_to_follow),
		members_toward_any = sum(toward >= 1),
		members_toward_two_or_more = sum(toward >= 2),
		members_away_any = sum(away >= 1)
	), by = .(level = as.character(get(column)))]
}

response_by_group = rbindlist(
	c(
		list(
			member_responses[, .(
				grouping = "all",
				level = "all members except the United States",
				members = .N,
				members_with_paired_votes = sum(paired > 0),
				switch_pairs = sum(paired),
				room_to_follow = sum(room_to_follow),
				toward = sum(toward),
				away = sum(away),
				unchanged = sum(unchanged),
				changed_same_distance = sum(changed_same_distance),
				not_paired = sum(not_paired),
				unchanged_share_of_paired = sum(unchanged) / sum(paired),
				toward_share_of_room = sum(toward) / sum(room_to_follow),
				members_toward_any = sum(toward >= 1),
				members_toward_two_or_more = sum(toward >= 2),
				members_away_any = sum(away >= 1)
			)]
		),
		lapply(
			c("paper_group", "eu_member", "treaty_group_2025", "mnna"),
			summarize_responses
		)
	),
	use.names = TRUE
)
setcolorder(response_by_group, c("grouping", "level"))
followers_two_or_more = member_responses[toward >= 2]
response_by_switch = switch_votes[, .(
	us_vote_2024 = us_vote_2024[1],
	us_vote_2025 = us_vote_2025[1],
	paired = sum(paired),
	room_to_follow = sum(room_to_follow),
	toward = sum(response == "toward new US vote"),
	away = sum(response == "away from new US vote"),
	unchanged = sum(response == "unchanged"),
	not_paired = sum(!paired),
	members_toward = paste(sort(iso3[response == "toward new US vote"]), collapse = " "),
	members_away = paste(sort(iso3[response == "away from new US vote"]), collapse = " ")
), by = .(normalized_title, resolution_2024, resolution_2025)]
response_by_switch = merge(
	response_by_switch,
	switches[, .(
		resolution_2025,
		total_yes_2024,
		total_no_2024,
		total_abstentions_2024,
		total_yes_2025,
		total_no_2025,
		total_abstentions_2025
	)],
	by = "resolution_2025"
)
setorder(response_by_switch, normalized_title)

methods = read_input("matched_methods")
resolutions_2025 = sort(unique(votes[year == 2025, resolution]))
agenda_groups = list(
	newly_recorded = methods[method_2024 == "without_vote" & method_2025 == "with_vote", resolution_2025],
	recurring_recorded = methods[method_2024 == "with_vote" & method_2025 == "with_vote", resolution_2025]
)
agenda_groups$other_recorded = setdiff(resolutions_2025, unlist(agenda_groups))
add_check(
	"newly recorded resolutions",
	length(agenda_groups$newly_recorded),
	33L
)
add_check(
	"recurring recorded resolutions",
	length(agenda_groups$recurring_recorded),
	58L
)
add_check(
	"other recorded resolutions",
	length(agenda_groups$other_recorded),
	101L
)
newly = read_input("newly_recorded_votes")
stopifnot(setequal(unique(newly$resolution), agenda_groups$newly_recorded))
archive_newly = votes[year == 2025 & resolution %in% agenda_groups$newly_recorded, .(
	resolution,
	voting_country_code,
	archive_vote = recorded_vote
)]
newly_check = merge(
	newly[, .(resolution, voting_country_code, recorded_vote)],
	archive_newly,
	by = c("resolution", "voting_country_code")
)
add_check(
	"newly recorded member votes agree with the archive",
	sum(newly_check$recorded_vote == newly_check$archive_vote),
	nrow(newly)
)
us_no_2025 = votes[year == 2025 & voting_country_code == "USA" & recorded_vote == "N", resolution]
count_group_votes = function(ids, label) {
	votes[year == 2025 & resolution %in% ids, .(
		resolution_group = label,
		resolutions = length(ids),
		us_no_resolutions = sum(ids %in% us_no_2025),
		votes_cast = sum(recorded_vote != "X"),
		yes = sum(recorded_vote == "Y"),
		abstain = sum(recorded_vote == "A"),
		no = sum(recorded_vote == "N"),
		absent = sum(recorded_vote == "X"),
		no_with_us_no = sum(recorded_vote == "N" & resolution %in% us_no_2025),
		no_resolutions = paste(sort(resolution[recorded_vote == "N"]), collapse = " ")
	), by = .(iso3 = voting_country_code)]
}

no_votes_long = rbind(
	count_group_votes(agenda_groups$newly_recorded, "newly_recorded"),
	count_group_votes(agenda_groups$other_recorded, "other_recorded"),
	count_group_votes(agenda_groups$recurring_recorded, "recurring_recorded")
)
no_votes_long = merge(no_votes_long, group_columns, by = "iso3")
setorder(no_votes_long, resolution_group, -no, iso3)
no_votes = dcast(
	no_votes_long[resolution_group != "recurring_recorded"],
	iso3 ~ resolution_group,
	value.var = c("votes_cast", "yes", "abstain", "no", "absent", "no_with_us_no")
)
no_votes = merge(no_votes, group_columns, by = "iso3")
no_votes[, no_total := no_newly_recorded + no_other_recorded]
no_votes[, no_with_us_no_total := no_with_us_no_newly_recorded + no_with_us_no_other_recorded]
setorder(no_votes, -no_newly_recorded, -no_other_recorded, iso3)
members_three_or_more = no_votes[no_newly_recorded >= 3 | no_other_recorded >= 3]
shares = read_input("vote_shares")
add_check(
	"allied no votes on the newly recorded items agree with vote_shares_by_resolution_group.csv",
	no_votes[iso3 %in% allies, sum(no_newly_recorded)],
	shares[resolution_group == "newly_recorded" & actor_group == "allies", no]
)
add_check(
	"allied votes cast on the newly recorded items agree with vote_shares_by_resolution_group.csv",
	no_votes[iso3 %in% allies, sum(votes_cast_newly_recorded)],
	shares[resolution_group == "newly_recorded" & actor_group == "allies", votes]
)
add_check(
	"allied no votes on the other recorded items agree with vote_shares_by_resolution_group.csv",
	no_votes[iso3 %in% allies, sum(no_other_recorded)],
	shares[resolution_group == "other_recorded" & actor_group == "allies", no]
)
add_check(
	"allied votes cast on the other recorded items agree with vote_shares_by_resolution_group.csv",
	no_votes[iso3 %in% allies, sum(votes_cast_other_recorded)],
	shares[resolution_group == "other_recorded" & actor_group == "allies", votes]
)
add_check(
	"other members no votes on the newly recorded items agree with vote_shares_by_resolution_group.csv",
	no_votes[paper_group %in% c("rivals", "other members"), sum(no_newly_recorded)],
	shares[resolution_group == "newly_recorded" & actor_group == "other_members", no]
)
no_votes_by_group = rbindlist(
	lapply(
		c("paper_group", "eu_member", "treaty_group_2025", "mnna"),
		function(column) {
			no_votes[, .(
				grouping = column,
				members = .N,
				votes_cast_newly_recorded = sum(votes_cast_newly_recorded),
				no_newly_recorded = sum(no_newly_recorded),
				no_with_us_no_newly_recorded = sum(no_with_us_no_newly_recorded),
				no_share_newly_recorded = sum(no_newly_recorded) / sum(votes_cast_newly_recorded),
				members_no_newly_recorded = sum(no_newly_recorded >= 1),
				votes_cast_other_recorded = sum(votes_cast_other_recorded),
				no_other_recorded = sum(no_other_recorded),
				no_with_us_no_other_recorded = sum(no_with_us_no_other_recorded),
				no_share_other_recorded = sum(no_other_recorded) / sum(votes_cast_other_recorded),
				members_no_other_recorded = sum(no_other_recorded >= 1),
				members_three_or_more_either = sum(no_newly_recorded >= 3 | no_other_recorded >= 3),
				members_three_or_more_with_us_no = sum(no_with_us_no_newly_recorded >= 3 | no_with_us_no_other_recorded >= 3)
			), by = .(level = as.character(get(column)))]
		}
	)
)
setcolorder(no_votes_by_group, c("grouping", "level"))
newly_or_other = c(agenda_groups$newly_recorded, agenda_groups$other_recorded)
no_voters_by_resolution = votes[year == 2025 & resolution %in% newly_or_other, .(
	resolution_group = fifelse(
		resolution[1] %in% agenda_groups$newly_recorded,
		"newly_recorded",
		"other_recorded"
	),
	date = date[1],
	title = title[1],
	total_yes = total_yes[1],
	total_no = total_no[1],
	total_abstentions = total_abstentions[1],
	us_vote = recorded_vote[voting_country_code == "USA"],
	no_votes = sum(recorded_vote == "N"),
	allied_no = sum(voting_country_code %in% allies & recorded_vote == "N"),
	no_voters = paste(sort(voting_country_code[recorded_vote == "N"]), collapse = " ")
), by = resolution]
setorder(no_voters_by_resolution, resolution_group, date, resolution)

us_votes = read_input("us_votes")
second_2025 = us_votes[cohort_year == 2025 & !is.na(us_recorded_vote)]
second_2024 = us_votes[cohort_year == 2024 & !is.na(us_recorded_vote)]
add_check(
	"second committee resolutions with an american no in 2025",
	second_2025[us_recorded_vote == "N", .N],
	41L
)
add_check(
	"second committee resolutions recorded in 2024",
	nrow(second_2024),
	8L
)
second_ids = rbind(
	second_2025[, .(cohort_year, resolution = plenary_final_symbol, root_id, us_recorded_vote)],
	second_2024[, .(cohort_year, resolution = plenary_final_symbol, root_id, us_recorded_vote)]
)
stopifnot(
	all(second_ids$resolution %in% votes[year %in% 2024:2025, resolution])
)
second_long = merge(
	votes[year %in% 2024:2025 & voting_country_code %in% c("USA", allies), .(
		resolution,
		year,
		date,
		title,
		total_yes,
		total_no,
		total_abstentions,
		iso3 = voting_country_code,
		vote = recorded_vote
	)],
	second_ids,
	by = "resolution"
)
stopifnot(all(second_long$year == second_long$cohort_year))
add_check(
	"archived american votes agree with us_votes.csv on the second committee items",
	second_long[iso3 == "USA", sum(vote == us_recorded_vote)],
	nrow(second_ids)
)
second_wide = dcast(
	second_long,
	cohort_year + resolution + root_id + date + title + total_yes + total_no + total_abstentions ~ iso3,
	value.var = "vote"
)
setcolorder(
	second_wide,
	c(
		"cohort_year",
		"resolution",
		"root_id",
		"date",
		"title",
		"total_yes",
		"total_no",
		"total_abstentions",
		"USA",
		allies
	)
)
second_by_resolution = second_long[iso3 != "USA", .(
	us_vote = us_recorded_vote[1],
	allies_yes = sum(vote == "Y"),
	allies_abstain = sum(vote == "A"),
	allies_no = sum(vote == "N"),
	allies_absent = sum(vote == "X"),
	allies_abstaining = paste(sort(iso3[vote == "A"]), collapse = " "),
	allies_no_voters = paste(sort(iso3[vote == "N"]), collapse = " ")
), by = .(
	cohort_year,
	resolution,
	date,
	title,
	total_yes,
	total_no,
	total_abstentions
)]
second_by_resolution[, agenda_group := fcase(
	cohort_year == 2024,
	"recorded in 2024",
	resolution %in% agenda_groups$newly_recorded,
	"newly_recorded",
	resolution %in% agenda_groups$recurring_recorded,
	"recurring_recorded",
	resolution %in% agenda_groups$other_recorded,
	"other_recorded"
)]
setorder(second_by_resolution, cohort_year, resolution)
second_by_agenda_group = second_by_resolution[, .(
	resolutions = .N,
	allied_votes = 14L * .N,
	allies_yes = sum(allies_yes),
	allies_abstain = sum(allies_abstain),
	allies_no = sum(allies_no),
	allies_absent = sum(allies_absent),
	all_fourteen_yes = sum(allies_yes == 14),
	majority_yes = sum(allies_yes > 7),
	any_allied_no = sum(allies_no > 0),
	all_fourteen_abstain = sum(allies_abstain == 14),
	all_fourteen_no = sum(allies_no == 14)
), by = .(cohort_year, us_vote, agenda_group)]
setorder(second_by_agenda_group, cohort_year, us_vote, agenda_group)
second_wide = merge(
	second_wide,
	second_by_resolution[, .(
		resolution,
		agenda_group,
		allies_yes,
		allies_abstain,
		allies_no,
		allies_absent
	)],
	by = "resolution"
)
setorder(second_wide, cohort_year, resolution)
second_by_ally = second_long[iso3 != "USA", .(
	resolutions = .N,
	yes = sum(vote == "Y"),
	abstain = sum(vote == "A"),
	no = sum(vote == "N"),
	absent = sum(vote == "X")
), by = .(cohort_year, iso3)]
second_by_ally = merge(
	second_by_ally,
	group_columns[, .(iso3, member_name, eu_member)],
	by = "iso3"
)
setorder(second_by_ally, cohort_year, iso3)
second_totals = second_long[iso3 != "USA", .(
	resolutions = uniqueN(resolution),
	allied_votes = .N,
	yes = sum(vote == "Y"),
	abstain = sum(vote == "A"),
	no = sum(vote == "N"),
	absent = sum(vote == "X"),
	yes_share = mean(vote == "Y"),
	no_share = mean(vote == "N")
), by = .(cohort_year, us_recorded_vote)]
second_support = second_by_resolution[, .(
	all_fourteen_yes = sum(allies_yes == 14),
	majority_yes = sum(allies_yes > 7),
	any_allied_no = sum(allies_no > 0),
	no_allied_no = sum(allies_no == 0),
	all_fourteen_abstain = sum(allies_abstain == 14),
	all_fourteen_no = sum(allies_no == 14)
), by = .(cohort_year, us_vote)]
second_totals = merge(
	second_totals,
	second_support,
	by.x = c("cohort_year", "us_recorded_vote"),
	by.y = c("cohort_year", "us_vote")
)
setorder(second_totals, cohort_year, us_recorded_vote)
recorded_totals = read_input("recorded_totals")
tax_votes = votes[resolution %in% tax_pair & voting_country_code %in% c("USA", allies), .(
	resolution,
	year,
	date,
	title,
	total_yes,
	total_no,
	total_abstentions,
	iso3 = voting_country_code,
	vote = recorded_vote
)]
tax_archive_totals = unique(
	tax_votes[, .(resolution, year, total_yes, total_no, total_abstentions)]
)
tax_validation = recorded_totals[resolution %in% tax_pair, .(
	resolution,
	official_yes,
	official_no,
	official_abstain,
	totals_agree
)]
tax_totals = merge(tax_archive_totals, tax_validation, by = "resolution")
add_check(
	"tax pair totals agree with recorded_total_validation.csv",
	c(tax_totals$total_yes, tax_totals$total_no, tax_totals$total_abstentions),
	c(
		tax_totals$official_yes,
		tax_totals$official_no,
		tax_totals$official_abstain
	)
)
tax_wide = dcast(tax_votes, iso3 ~ year, value.var = "vote")
setnames(tax_wide, c("2024", "2025"), c("vote_2024", "vote_2025"))
tax_wide = merge(tax_wide, group_columns[, .(iso3, member_name, eu_member)], by = "iso3")
tax_wide[, change := fifelse(
	vote_2024 == vote_2025,
	"unchanged",
	paste(vote_2024, "to", vote_2025)
)]
setorder(tax_wide, -vote_2024, iso3)
tax_summary = tax_votes[iso3 != "USA", .(
	us_vote = tax_votes[iso3 == "USA" & resolution == .BY$resolution, vote],
	allies_yes = sum(vote == "Y"),
	allies_abstain = sum(vote == "A"),
	allies_no = sum(vote == "N"),
	allies_absent = sum(vote == "X"),
	allies_no_voters = paste(sort(iso3[vote == "N"]), collapse = " ")
), by = .(resolution, year, date, total_yes, total_no, total_abstentions)]
tax_summary = merge(
	tax_summary,
	votes[resolution %in% tax_pair, .(
		members_no = sum(recorded_vote == "N"),
		no_voters = paste(sort(voting_country_code[recorded_vote == "N"]), collapse = " ")
	), by = resolution],
	by = "resolution"
)
setorder(tax_summary, year)
# the fourteen allies one by one
peers = read_input("country_peer_results")[
	horizon == 1 & baseline == 2024 & target == 2025 & specification == "ordinal" & grepl("^peer_", group)
]
peers[, iso3 := sub("^peer_", "", group)]
resets = read_input("reset_probabilities")
ally_table = merge(
	ally_responses[, .(
		iso3 = ally,
		switch_pairs = us_switch_pairs,
		votes_kept = ally_unchanged_when_us_switches,
		toward = closer_to_new_us_vote,
		away = further_from_new_us_vote
	)],
	peers[, .(
		iso3,
		peer_agreement_2024 = baseline_agreement,
		peer_agreement_2025 = target_agreement,
		peer_change = change
	)],
	by = "iso3"
)
ally_table = merge(
	ally_table,
	resets[year == 2025, .(
		iso3 = actor,
		break_probability_2025 = probability,
		share_cast_2025 = share_cast
	)],
	by = "iso3"
)
ally_table = merge(
	ally_table,
	groups[, .(
		iso3,
		member_name,
		eu_member,
		eu_accession_year,
		us_treaty,
		us_treaty_in_force_year,
		mnna
	)],
	by = "iso3"
)
ally_table[, votes_kept_share := votes_kept / switch_pairs]
ally_table[, ally_group := fifelse(eu_member, "EU member", "non-EU ally")]
displacement_source = NA_character_
if (file.exists(displacement_path)) {
	displacements = read_analysis_table(displacement_path)
	displacement_column = grep("sigma", names(displacements), value = TRUE)[1]
	if (all(c("actor", "year") %in% names(displacements)) && !is.na(displacement_column)) {
		ally_table = merge(
			ally_table,
			displacements[year == 2025, .(iso3 = actor, displacement_sigma_2025 = get(displacement_column))],
			by = "iso3",
			all.x = TRUE
		)
		displacement_source = paste(
			"results/models/un_variants/threshold/displacements.csv",
			displacement_column
		)
	}
}

if (!"displacement_sigma_2025" %in% names(ally_table)) ally_table[, displacement_sigma_2025 := NA_real_]
setcolorder(
	ally_table,
	c(
		"iso3",
		"member_name",
		"ally_group",
		"eu_member",
		"eu_accession_year",
		"us_treaty",
		"us_treaty_in_force_year",
		"mnna",
		"switch_pairs",
		"votes_kept",
		"votes_kept_share",
		"toward",
		"away",
		"peer_agreement_2024",
		"peer_agreement_2025",
		"peer_change",
		"break_probability_2025",
		"share_cast_2025",
		"displacement_sigma_2025"
	)
)
setorder(ally_table, ally_group, iso3)
summarize_allies = function(d, label) {
	d[, .(
		ally_group = label,
		allies = .N,
		switch_pairs = sum(switch_pairs),
		votes_kept = sum(votes_kept),
		votes_kept_share = sum(votes_kept) / sum(switch_pairs),
		toward = sum(toward),
		away = sum(away),
		mean_peer_change = mean(peer_change),
		min_peer_change = min(peer_change),
		max_peer_change = max(peer_change),
		mean_break_probability = mean(break_probability_2025),
		min_break_probability = min(break_probability_2025),
		max_break_probability = max(break_probability_2025),
		mean_displacement_sigma = mean(displacement_sigma_2025)
	)]
}

ally_group_means = rbind(
	summarize_allies(ally_table[ally_group == "EU member"], "EU member"),
	summarize_allies(ally_table[ally_group == "non-EU ally"], "non-EU ally"),
	summarize_allies(ally_table, "all fourteen")
)
# the most plausible followers outside the group
follower_table = resets[year == 2025 & actor %in% plausible_followers, .(
	iso3 = actor,
	break_probability_2025 = probability,
	chain_min = chain_min,
	chain_max = chain_max,
	share_cast_2025 = share_cast
)]
follower_table = merge(
	follower_table,
	member_responses[, .(
		iso3,
		member_name,
		paired,
		toward,
		away,
		unchanged,
		not_paired,
		titles_toward,
		titles_away
	)],
	by = "iso3",
	all.x = TRUE
)
follower_table = merge(
	follower_table,
	no_votes[, .(
		iso3,
		no_newly_recorded,
		no_with_us_no_newly_recorded,
		votes_cast_newly_recorded,
		no_other_recorded,
		no_with_us_no_other_recorded,
		votes_cast_other_recorded
	)],
	by = "iso3",
	all.x = TRUE
)
follower_table = merge(
	follower_table,
	group_columns[, .(iso3, eu_member, treaty_group_2025, us_treaty, mnna)],
	by = "iso3",
	all.x = TRUE
)
earlier_breaks = resets[actor %in% plausible_followers & probability >= 0.5 & year < 2025, .(
	years_with_break_probability_at_least_one_half = paste(year, collapse = " ")
), by = .(iso3 = actor)]
follower_table = merge(follower_table, earlier_breaks, by = "iso3", all.x = TRUE)
follower_table[, reweighted_break_probability_2025 := NA_real_]
setorder(follower_table, -break_probability_2025)
follower_switch_votes = switch_votes[iso3 %in% plausible_followers, .(
	normalized_title,
	resolution_2024,
	resolution_2025,
	us_vote_2024,
	us_vote_2025,
	iso3,
	vote_2024,
	vote_2025,
	response
)]
setorder(follower_switch_votes, iso3, normalized_title)
follower_history = resets[actor %in% plausible_followers, .(iso3 = actor, year, share_cast, probability)]
setorder(follower_history, iso3, year)
# four named votes
named = votes[resolution %in% named_resolutions]
stopifnot(uniqueN(named[, .(year, resolution)]) == length(named_resolutions))
named_summary = named[, .(
	year = year[1],
	date = date[1],
	title = title[1],
	total_yes = total_yes[1],
	total_no = total_no[1],
	total_abstentions = total_abstentions[1],
	us_vote = recorded_vote[voting_country_code == "USA"],
	allies_yes = sum(voting_country_code %in% allies & recorded_vote == "Y"),
	allies_abstain = sum(voting_country_code %in% allies & recorded_vote == "A"),
	allies_no = sum(voting_country_code %in% allies & recorded_vote == "N"),
	allies_absent = sum(voting_country_code %in% allies & recorded_vote == "X"),
	allies_abstaining = paste(
		sort(
			voting_country_code[voting_country_code %in% allies & recorded_vote == "A"]
		),
		collapse = " "
	),
	allies_no_voters = paste(
		sort(
			voting_country_code[voting_country_code %in% allies & recorded_vote == "N"]
		),
		collapse = " "
	),
	rivals_yes = sum(voting_country_code %in% rivals & recorded_vote == "Y"),
	rivals_abstain = sum(voting_country_code %in% rivals & recorded_vote == "A"),
	rivals_no = sum(voting_country_code %in% rivals & recorded_vote == "N"),
	rivals_absent = sum(voting_country_code %in% rivals & recorded_vote == "X"),
	no_voters = paste(sort(voting_country_code[recorded_vote == "N"]), collapse = " ")
), by = resolution]
setorder(named_summary, date)
named_votes = dcast(
	named[voting_country_code %in% c("USA", allies), .(resolution, year, iso3 = voting_country_code, vote = recorded_vote)],
	resolution + year ~ iso3,
	value.var = "vote"
)
setcolorder(named_votes, c("resolution", "year", "USA", allies))
setorder(named_votes, year, resolution)
# recorded votes cast by the sparse voters in 2025
sparse = votes[year == 2025 & voting_country_code %in% c(sparse_voters, "AFG"), .(
	member_name = voting_country_name[1],
	resolutions_2025 = .N,
	recorded_votes_2025 = sum(recorded_vote != "X"),
	yes = sum(recorded_vote == "Y"),
	abstain = sum(recorded_vote == "A"),
	no = sum(recorded_vote == "N"),
	absent = sum(recorded_vote == "X"),
	share_recorded_2025 = mean(recorded_vote != "X"),
	first_recorded_vote_2025 = if (any(recorded_vote != "X")) as.character(min(date[recorded_vote != "X"])) else NA_character_,
	last_recorded_vote_2025 = if (any(recorded_vote != "X")) as.character(max(date[recorded_vote != "X"])) else NA_character_
), by = .(iso3 = voting_country_code)]
sparse = merge(
	sparse,
	votes[year == 2024 & voting_country_code %in% c(sparse_voters, "AFG"), .(
		resolutions_2024 = .N,
		recorded_votes_2024 = sum(recorded_vote != "X")
	), by = .(iso3 = voting_country_code)],
	by = "iso3"
)
sparse = merge(
	sparse,
	resets[year == 2025, .(
		iso3 = actor,
		model_share_cast_2025 = share_cast,
		break_probability_2025 = probability,
		regular
	)],
	by = "iso3",
	all.x = TRUE
)
setorder(sparse, recorded_votes_2025)
# outputs
outputs = list(
	us_switches = switches,
	switch_votes = switch_votes,
	member_responses = member_responses,
	followers_two_or_more = followers_two_or_more,
	response_by_group = response_by_group,
	response_by_switch = response_by_switch,
	no_votes_by_member = no_votes,
	no_votes_by_member_long = no_votes_long,
	members_three_or_more_no = members_three_or_more,
	no_votes_by_group = no_votes_by_group,
	no_voters_by_resolution = no_voters_by_resolution,
	second_committee_votes = second_wide,
	second_committee_by_resolution = second_by_resolution,
	second_committee_by_ally = second_by_ally,
	second_committee_totals = second_totals,
	second_committee_by_agenda_group = second_by_agenda_group,
	tax_pair_votes = tax_wide,
	tax_pair_summary = tax_summary,
	ally_table = ally_table,
	ally_group_means = ally_group_means,
	plausible_followers = follower_table,
	plausible_follower_switch_votes = follower_switch_votes,
	plausible_follower_history = follower_history,
	named_votes_summary = named_summary,
	named_votes = named_votes,
	sparse_voters_2025 = sparse
)
checks_table = rbindlist(checks)
stopifnot(all(checks_table$passed))
outputs$checks = checks_table
for (name in names(outputs)) {
	fwrite(outputs[[name]], file.path(output_root, paste0(name, ".csv")))
}

write_json(
	list(
		script = "code/16_followers.R",
		completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
		inputs = lapply(
			inputs,
			function(path) {
				list(
					path = path,
					sha256 = digest(file = file.path(root, path), algo = "sha256")
				)
			}
		),
		displacement_source = displacement_source,
		counts = list(
			members = length(members),
			american_switches = nrow(switches),
			newly_recorded = length(agenda_groups$newly_recorded),
			recurring_recorded = length(agenda_groups$recurring_recorded),
			other_recorded = length(agenda_groups$other_recorded),
			second_committee_2025 = nrow(second_2025),
			second_committee_2024 = nrow(second_2024)
		),
		allies = allies,
		rivals = rivals,
		plausible_followers = plausible_followers,
		sparse_voters = sparse_voters,
		named_resolutions = named_resolutions,
		tax_pair = tax_pair,
		outputs = paste0(names(outputs), ".csv"),
		R = R.version.string
	),
	file.path(output_root, "provenance.json"),
	pretty = TRUE,
	auto_unbox = TRUE
)
