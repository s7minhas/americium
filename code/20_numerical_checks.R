# check reported numerical summaries

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
set.seed(6886)
setDTthreads(1L)
output = "results/revision"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
check_store = new.env(parent = emptyenv())
check_store$rows = list()
add_check = function(name, actual, expected, tolerance = 1e-10) {
	check_store$rows[[length(check_store$rows) + 1L]] = data.table(
		check = name,
		actual = actual,
		expected = expected,
		passed = isTRUE(
			all.equal(actual, expected, tolerance = tolerance, check.attributes = FALSE)
		)
	)
}

variants = read_analysis_table(file.path(output, "variant_comparison.csv"))
for (variant in variants$variant) {
	path = file.path("results/models/un_variants", variant)
	features = readRDS(file.path(path, "chain_features.rds"))[[variant]]
	summary = read_analysis_table(file.path(path, "variant_summary.csv"))
	continuous = do.call(rbind, lapply(features, `[[`, "continuous"))
	for (parameter in c("sigma_uv", "rho_uv", "pi_snap", "ve", "va", "intercept")) {
		add_check(
			paste(variant, parameter),
			mean(continuous[, parameter]),
			summary[quantity == parameter, mean]
		)
	}
	terminal = read_analysis_table(file.path(path, "terminal_probabilities.csv"))
	for (actor_value in terminal$actor) {
		id = paste(actor_value, terminal[actor == actor_value, year], sep = "_")
		probability = mean(
			vapply(
				features,
				function(chain) mean(plogis(chain$actor_log_odds[, id])),
				0.0
			)
		)
		add_check(
			paste(variant, actor_value, "probability"),
			probability,
			terminal[actor == actor_value, probability]
		)
	}
	manifest = read_analysis_table(file.path(path, "fit_manifest.csv"))
	add_check(
		paste(variant, "completed chains"),
		as.numeric(all(manifest$complete)),
		1
	)
}

features = readRDS("results/models/un/chain_features.rds")$full
joint_probability = mean(
	vapply(
		features,
		function(chain) {
			columns = grep("^USA_", colnames(chain$actor_log_odds))
			probability = plogis(chain$actor_log_odds[, columns, drop = FALSE])
			terminal = which(colnames(probability) == "USA_2025")
			mean(
				probability[, terminal] * apply(1 - probability[, -terminal, drop = FALSE], 1, prod)
			)
		},
		0.0
	)
)
fwrite(
	data.table(
		quantity = "american_first_break_2025",
		probability = joint_probability
	),
	file.path(output, "joint_probability.csv")
)

outcomes = read_analysis_table(
	"results/proposal_comparison/reading_agreement/outcome_agreement.csv"
)[scope == "all_substantive_pairs"]
agreement = read_analysis_table(
	"results/proposal_comparison/reading_agreement/batch_agreement.csv"
)[scope == "all_substantive_pairs"]
observed = outcomes[first_outcome == independent_outcome, sum(n)] / sum(outcomes$n)
first = outcomes[, .(count = sum(n)), by = first_outcome]
second = outcomes[, .(count = sum(n)), by = independent_outcome]
chance = merge(first, second, by.x = "first_outcome", by.y = "independent_outcome")[, sum(count.x * count.y)] / sum(outcomes$n)^2
add_check("outcome agreement", observed * 100, agreement$percent_agreement)
add_check(
	"outcome kappa",
	(observed - chance) / (1 - chance),
	agreement$kappa
)

votes_environment = new.env()
load_vote_archive("data/un_votes_long.rda", envir = votes_environment)
votes = as.data.table(votes_environment$un_votes_long)
score = c(Y = 1, A = 0, N = -1)
switches = read_analysis_table("results/observed/followers/us_switches.csv")
responses = read_analysis_table("results/observed/followers/switch_votes.csv")
before_votes = votes[year == 2024L, .(
	resolution_2024 = resolution,
	iso3 = voting_country_code,
	earlier_vote = recorded_vote
)]
after_votes = votes[year == 2025L, .(
	resolution_2025 = resolution,
	iso3 = voting_country_code,
	later_vote = recorded_vote
)]
paired_votes = merge(
	responses[, .(resolution_2024, resolution_2025, iso3, us_vote_2025)],
	before_votes,
	by = c("resolution_2024", "iso3"),
	all.x = TRUE
)
paired_votes = merge(
	paired_votes,
	after_votes,
	by = c("resolution_2025", "iso3"),
	all.x = TRUE
)
earlier = unname(score[paired_votes$earlier_vote])
later = unname(score[paired_votes$later_vote])
new_us = unname(score[paired_votes$us_vote_2025])
paired = is.finite(earlier) & is.finite(later)
toward = paired & abs(later - new_us) < abs(earlier - new_us)
away = paired & abs(later - new_us) > abs(earlier - new_us)
add_check("raw paired votes", sum(paired), sum(responses$paired))
add_check(
	"raw unchanged votes",
	sum(paired & earlier == later, na.rm = TRUE),
	sum(responses$response == "unchanged")
)
add_check(
	"raw votes toward Washington",
	sum(toward, na.rm = TRUE),
	sum(responses$response == "toward new US vote")
)
add_check(
	"raw votes away from Washington",
	sum(away, na.rm = TRUE),
	sum(responses$response == "away from new US vote")
)
add_check("fourteen switches", nrow(switches), 14L)
add_check(
	"paired response count",
	sum(responses$response != "not paired"),
	2307L
)
add_check(
	"unchanged response count",
	sum(responses$response == "unchanged"),
	2241L
)
# reproduce the observed no-vote share from the voting archive
profile = read_analysis_table("results/observed/annual/us_vote_profile.csv")
for (year_value in c(2024L, 2025L)) {
	us = votes[voting_country_code == "USA" & year == year_value & recorded_vote %in% names(score)]
	add_check(
		paste("American no share", year_value),
		mean(us$recorded_vote == "N"),
		profile[year == year_value, no_share_recorded]
	)
}

result = rbindlist(check_store$rows)
fwrite(result, file.path(output, "numerical_checks.csv"))
write_json(
	list(
		script = "code/20_numerical_checks.R",
		checks = nrow(result),
		passed = all(result$passed),
		variants = variants$variant,
		feature_hashes = setNames(
			lapply(
				variants$variant,
				function(variant) {
					digest(
						file = file.path("results/models/un_variants", variant, "chain_features.rds"),
						algo = "sha256"
					)
				}
			),
			variants$variant
		),
		script_sha256 = digest(file = "code/20_numerical_checks.R", algo = "sha256"),
		R = R.version.string
	),
	file.path(output, "numerical_checks_provenance.json"),
	pretty = TRUE,
	auto_unbox = TRUE
)
stopifnot(all(result$passed))
