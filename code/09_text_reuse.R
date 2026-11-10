# summarize wording reuse and save its figure

helper_file = if (sys.nframe() > 0L && !is.null(sys.frame(1)$ofile)) {
	sys.frame(1)$ofile
} else {
	sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
}

source(
	file.path(dirname(normalizePath(helper_file)), "00_data_helpers.R"),
	local = TRUE
)

library(data.table)
library(ggplot2)
source("code/00_plot_style.R")
result_dir = "results/text_reuse"
figure_dir = "results/figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
pairs_all = read_analysis_table(file.path(result_dir, "pairs.csv"))
pairs = pairs_all[link == "same_title"]
first_year = 1994
# year summaries with each pair weighted equally
year_share = function(x, by_columns) {
	x[, .(
		pairs = .N,
		carried_forward = 100 * mean(carried_80),
		new_paragraphs = 100 * mean(new_80)
	), by = by_columns]
}

series = year_share(
	pairs[!is.na(adoption_year) & adoption_year >= first_year & adoption_year <= 2025],
	"adoption_year"
)
setorder(series, adoption_year)
series[, rank_highest := frank(-carried_forward)]
series[, change := carried_forward - shift(carried_forward)]
series[, change_rank := frank(-abs(change), na.last = "keep")]
fwrite(series, file.path(result_dir, "series_main.csv"))
# how far 2024 and 2025 move when one committee or one pair is left out
leave_out = rbindlist(
	lapply(
		c(2024, 2025),
		function(year) {
			x = pairs[adoption_year == year]
			by_committee = sapply(
				unique(x$committee),
				function(k) 100 * mean(x[committee != k, carried_80])
			)
			by_pair = (sum(x$carried_80) - x$carried_80) / (nrow(x) - 1) * 100
			data.table(
				adoption_year = year,
				pairs = nrow(x),
				carried_forward = 100 * mean(x$carried_80),
				leave_one_committee_min = min(by_committee),
				leave_one_committee_max = max(by_committee),
				leave_one_pair_min = min(by_pair),
				leave_one_pair_max = max(by_pair)
			)
		}
	)
)
fwrite(leave_out, file.path(result_dir, "leave_out_ranges.csv"))
# committees, each against its own history
committee_labels = c(
	"0" = "Plenary, No Main Committee",
	"1" = "First Committee",
	"2" = "Second Committee",
	"3" = "Third Committee",
	"4" = "Fourth Committee",
	"5" = "Fifth Committee",
	"6" = "Sixth Committee"
)
by_committee = year_share(
	pairs[!is.na(adoption_year) & adoption_year >= first_year & adoption_year <= 2025],
	c("committee", "adoption_year")
)
by_committee = by_committee[pairs >= 5]
by_committee[, committee_label := factor(committee_labels[as.character(committee)], levels = committee_labels)]
by_committee[, rank_highest := frank(-carried_forward), by = committee]
by_committee[, years_ranked := .N, by = committee]
setorder(by_committee, committee, adoption_year)
fwrite(by_committee, file.path(result_dir, "series_committee.csv"))
# american votes on the earlier and later resolution
load_vote_archive("data/un_votes_long.rda")
us_votes = unique(
	un_votes_long[voting_country_code == "USA", .(resolution, us_vote = recorded_vote)],
	by = "resolution"
)
voted = merge(
	pairs,
	us_votes[, .(prior = resolution, prior_vote = us_vote)],
	by = "prior",
	all.x = TRUE
)
voted = merge(
	voted,
	us_votes[, .(later = resolution, later_vote = us_vote)],
	by = "later",
	all.x = TRUE
)
voted[, vote_group := fcase(
	is.na(later_vote),
	"Later resolution adopted without a vote",
	is.na(prior_vote) & later_vote == "N",
	"Newly put to a vote, United States voted No",
	is.na(prior_vote),
	"Newly put to a vote, other American vote",
	prior_vote == "N" & later_vote == "N",
	"United States voted No both times",
	prior_vote == later_vote,
	"Same American vote both times, not No",
	default = "American vote changed"
)]
vote_groups = voted[!is.na(adoption_year) & adoption_year >= 2010, .(
	pairs = .N,
	carried_forward = 100 * mean(carried_80),
	new_paragraphs = 100 * mean(new_80)
), by = .(adoption_year, vote_group)]
setorder(vote_groups, adoption_year, -pairs)
fwrite(vote_groups, file.path(result_dir, "us_vote_groups.csv"))

adoption_list = read_analysis_table("results/adoptions/history/all_adoption_units.csv")[, .(prior = canonical_resolution, prior_method = method)]
newly_voted = merge(
	voted[adoption_year == 2025 & vote_group == "Newly put to a vote, United States voted No"],
	unique(adoption_list, by = "prior"),
	by = "prior",
	all.x = TRUE
)
fwrite(
	newly_voted[, .(
		resolutions = .N,
		predecessor_without_vote = sum(prior_method == "without_vote", na.rm = TRUE)
	), by = .(sessions_since_predecessor = gap)],
	file.path(result_dir, "newly_voted_predecessors.csv")
)
# the hand review of resolutions newly put to a vote in 2025
review = read_analysis_table("data/resolution_review.csv")[review_set == "newly_recorded"]
review = merge(
	review[, .(prior = trimws(symbol_before), later = trimws(symbol_after), classification)],
	pairs[, .(prior, later, carried_80, new_80)],
	by = c("prior", "later"),
	all.x = TRUE
)
review_labels = c(
	routine_updates = "Routine updates",
	substantive_updates_same_question = "Substantive updates, same question",
	material_change_limiting_vote_interpretation = "Material change"
)
review[, classification_label := review_labels[classification]]
fwrite(review, file.path(result_dir, "newly_recorded_review_pairs.csv"))
fwrite(
	rbind(
		review[!is.na(carried_80), .(
			resolutions = .N,
			carried_forward = 100 * mean(carried_80),
			new_paragraphs = 100 * mean(new_80)
		), by = classification_label],
		review[!is.na(carried_80), .(
			classification_label = "All matched",
			resolutions = .N,
			carried_forward = 100 * mean(carried_80),
			new_paragraphs = 100 * mean(new_80)
		)]
	),
	file.path(result_dir, "newly_recorded_review.csv")
)
# agreement with the request coding of the second committee
archive = "data/proposal_comparison_sandbox.zip"
work_dir = tempfile("text_reuse_")
unzip(
	archive,
	files = c("tables/family_pairs.csv", "tables/family_outcome_fractions.csv"),
	exdir = work_dir
)
family_pairs = read_analysis_table(file.path(work_dir, "tables/family_pairs.csv"))
fractions = read_analysis_table(
	file.path(work_dir, "tables/family_outcome_fractions.csv")
)
fractions = fractions[subset == "all_core" & internal_decision == FALSE]
coded = fractions[, .(
	coded_unchanged = sum(fraction[outcome == "continued_without_substantive_change"])
), by = .(pair_id, prior_year, later_year)]
coded = merge(
	coded,
	family_pairs[, .(pair_id, prior = prior_adopted, later = later_adopted)],
	by = "pair_id"
)
coded = merge(coded, pairs[, .(prior, later, carried_80)], by = c("prior", "later"))
fwrite(coded, file.path(result_dir, "request_coding_pairs.csv"))
fwrite(
	rbind(
		coded[, .(
			pairs = .N,
			coded_unchanged = 100 * mean(coded_unchanged),
			carried_forward = 100 * mean(carried_80),
			correlation = cor(coded_unchanged, carried_80)
		), by = .(prior_year, later_year)],
		coded[, .(
			prior_year = NA,
			later_year = NA,
			pairs = .N,
			coded_unchanged = 100 * mean(coded_unchanged),
			carried_forward = 100 * mean(carried_80),
			correlation = cor(coded_unchanged, carried_80)
		)]
	),
	file.path(result_dir, "request_coding_comparison.csv")
)
# figures
plot_theme = paper_plot_theme()
marked = data.table(start = c(2016.5, 2024.5), end = c(2020.5, 2025.5))
figure_main = ggplot(series, aes(x = adoption_year, y = carried_forward))
figure_main = figure_main + geom_rect(
	data = marked,
	aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
	inherit.aes = FALSE,
	alpha = 0.12
)
figure_main = figure_main + geom_line()
figure_main = figure_main + geom_point()
figure_main = figure_main + scale_x_continuous(breaks = seq(1995, 2025, 5))
figure_main = figure_main + labs(x = "Year of adoption", y = "Operative paragraphs carried forward (%)")
figure_main = figure_main + expand_limits(y = c(0, 100))
figure_main = figure_main + plot_theme
ggsave(
	filename = file.path(figure_dir, "text_carried_forward.pdf"),
	plot = figure_main,
	width = 6.8,
	height = 3.8,
	device = cairo_pdf
)
