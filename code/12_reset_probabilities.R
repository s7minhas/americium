# calculate country reset probabilities and sensitivity checks

helper_file = if (sys.nframe() > 0L && !is.null(sys.frame(1)$ofile)) {
	sys.frame(1)$ofile
} else {
	sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
}

source(
	file.path(dirname(normalizePath(helper_file)), "00_data_helpers.R"),
	local = TRUE
)
suppressMessages(
	{
		library(data.table)
	}
)
data_root = "data"
results_root = "results"
model_root = file.path(results_root, "models", "un")
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
regular_share = .75
share_cutoffs = c(0, .5, .75, .9)
prior_rate = .001
labels = data.table(
	actor = c(
		"BGR",
		"CSK",
		"HUN",
		"POL",
		"ROU",
		"SUN",
		"UKR",
		"BLR",
		"ARG",
		"RUS",
		"BLR",
		"KAZ",
		"SRB",
		"ARG",
		"USA"
	),
	year = c(
		1990L,
		1990L,
		1990L,
		1990L,
		1990L,
		1991L,
		1991L,
		1991L,
		1991L,
		2000L,
		2000L,
		2000L,
		2008L,
		2024L,
		2025L
	),
	name = c(
		"Bulgaria",
		"Czechoslovakia",
		"Hungary",
		"Poland",
		"Romania",
		"Soviet Union",
		"Ukraine",
		"Belarus",
		"Argentina",
		"Russia",
		"Belarus",
		"Kazakhstan",
		"Serbia",
		"Argentina",
		"United States"
	)
)

vote_env = new.env()
load_vote_archive(file.path(data_root, "un_votes_long.rda"), envir = vote_env)
votes = as.data.table(vote_env$un_votes_long)[, .(actor = voting_country_code, year, recorded_vote)]
participation = votes[, .(share_cast = mean(recorded_vote %in% c("Y", "N", "A"))), by = .(actor, year)]
# conditional reset probabilities by chain, then averaged over chains
features = readRDS(file.path(model_root, "chain_features.rds"))$full
by_chain = rbindlist(
	lapply(
		seq_along(features),
		function(chain) {
			log_odds = features[[chain]]$actor_log_odds
			data.table(
				chain = chain,
				id = colnames(log_odds),
				probability = colMeans(plogis(log_odds))
			)
		}
	)
)
by_chain[, c("actor", "year") := tstrsplit(id, "_")]
by_chain[, year := as.integer(year)]
by_chain = merge(by_chain, participation, by = c("actor", "year"), all.x = TRUE)
probs = by_chain[, .(
	probability = mean(probability),
	chain_min = min(probability),
	chain_max = max(probability)
), by = .(actor, year, share_cast)]
probs[, group := fifelse(
	actor == "USA",
	"United States",
	fifelse(actor %in% allies, "Fourteen allies", "Other countries")
)]
probs[, regular := share_cast >= regular_share]
fwrite(
	probs[order(year, -probability)],
	file.path(model_root, "reset_probabilities.csv")
)
# participation in country-years with and without resets
probs[, reset := probability > .5]
repeat_actors = probs[reset == TRUE, .N, by = actor][N >= 3, actor]
participation_summary = rbind(
	probs[, .(
		country_years = .N,
		countries = uniqueN(actor),
		median_share_cast = median(share_cast),
		share_below_half = mean(share_cast < .5)
	), by = .(subset = fifelse(reset, "reset", "no reset"))],
	probs[reset == TRUE & actor %in% repeat_actors, .(
		subset = "reset, countries with three or more",
		country_years = .N,
		countries = uniqueN(actor),
		median_share_cast = median(share_cast),
		share_below_half = mean(share_cast < .5)
	)],
	probs[reset == TRUE & regular == TRUE, .(
		subset = "reset, regular voters",
		country_years = .N,
		countries = uniqueN(actor),
		median_share_cast = median(share_cast),
		share_below_half = mean(share_cast < .5)
	)]
)
fwrite(
	participation_summary,
	file.path(model_root, "reset_participation.csv")
)

counts = rbindlist(
	lapply(
		share_cutoffs,
		function(cutoff) {
			probs[actor != "USA" & share_cast >= cutoff, .(expected = sum(probability), countries = .N), by = year][, cutoff := cutoff]
		}
	)
)
counts[, rank := frank(-expected, ties.method = "min"), by = cutoff]
chain_counts = rbindlist(
	lapply(
		share_cutoffs,
		function(cutoff) {
			by_chain[actor != "USA" & share_cast >= cutoff, .(expected = sum(probability)), by = .(chain, year)][, cutoff := cutoff]
		}
	)
)
chain_spread = chain_counts[, .(
	chain_min = min(expected),
	chain_max = max(expected),
	chain_range = max(expected) - min(expected)
), by = .(cutoff, year)]
counts = merge(counts, chain_spread, by = c("cutoff", "year"))
fwrite(counts[order(cutoff, year)], file.path(model_root, "reset_counts.csv"))
# reset probabilities under other reset rates, holding other draws fixed
regular_rate = probs[regular == TRUE & actor != "USA", sum(probability) / .N]
reference_rates = c(prior_mean = prior_rate, regular_voters = regular_rate)
rate_check = rbindlist(
	lapply(
		seq_along(features),
		function(chain) {
			log_odds = features[[chain]]$actor_log_odds[, c("USA_2025", paste0(allies, "_2025"))]
			rate = features[[chain]]$continuous[, "pi_snap"]
			rbindlist(
				lapply(
					names(reference_rates),
					function(reference) {
						shifted = log_odds - qlogis(rate) + qlogis(reference_rates[[reference]])
						data.table(
							chain = chain,
							actor = sub("_2025", "", colnames(log_odds)),
							estimated_rate = mean(rate),
							reference = reference,
							reference_rate = reference_rates[[reference]],
							fitted = colMeans(plogis(log_odds)),
							at_reference_rate = colMeans(plogis(shifted))
						)
					}
				)
			)
		}
	)
)
fwrite(rate_check, file.path(model_root, "reset_rate_check.csv"))
# selected country-years with high reset probabilities
fwrite(
	merge(labels, probs, by = c("actor", "year"))[order(year)],
	file.path(model_root, "reset_selected_cases.csv")
)
# detection of planted resets in the final period of the simulated networks
control_root = file.path(results_root, "models")
truth = readRDS(file.path(control_root, "control_truth.rds"))
cases = read_analysis_table(file.path(control_root, "control_cases.csv"))
control_probs = read_analysis_table(file.path(control_root, "controls", "actors.csv"))[, .(probability = mean(probability)), by = .(actor, specification)]
final_period = rbindlist(
	lapply(
		seq_len(nrow(cases)),
		function(row) {
			planted = truth[[as.character(cases$case_id[row])]]$truth$reset_indicator[, "2012"]
			fitted = control_probs[specification == cases$specification[row]]
			if (!nrow(fitted)) {
				return(NULL)
			}
			fitted[, planted := planted[actor]]
			data.table(
				case = cases$case[row],
				event_period = cases$event_period[row],
				reset_sd = cases$reset_sd[row],
				planted_count = sum(fitted$planted),
				planted_mean = fitted[planted == 1, mean(probability)],
				planted_above_half = fitted[planted == 1, sum(probability > .5)],
				unplanted_max = fitted[planted == 0, max(probability)]
			)
		}
	)
)
detection = final_period[event_period == 12L, .(
	arrays = .N,
	planted = sum(planted_count),
	planted_mean = mean(planted_mean, na.rm = TRUE),
	planted_above_half = sum(planted_above_half),
	unplanted_max = max(unplanted_max)
), by = .(case, reset_sd)]
fwrite(
	detection,
	file.path(control_root, "controls", "final_period_detection.csv")
)
