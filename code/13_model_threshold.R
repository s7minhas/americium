# calculate displacements and reset thresholds

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
		library(jsonlite)
		library(digest)
	}
)
set.seed(6886)
script_path = file.path("code", "13_model_threshold.R")
model_root = file.path("results", "models", "un")
variant_root = file.path("results", "models", "un_variants")
threshold_root = file.path(variant_root, "threshold")
decomposition_root = file.path(variant_root, "decomposition")
kappa_root = file.path(variant_root, "kappa_sensitivity")
for (path in c(threshold_root, decomposition_root, kappa_root)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
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
focal = c("USA", allies)
kappa_fitted = 2
kappa_values = c(.5, 1, 4)
interval = c(.025, .975)
final_year = "2025"
input_paths = c(
	draws = file.path(model_root, "space", "draws.rds"),
	features = file.path(model_root, "chain_features.rds"),
	diagnostics = file.path(model_root, "diagnostics.csv"),
	positions = file.path(model_root, "space", "positions.csv"),
	probabilities = file.path(model_root, "reset_probabilities.csv")
)
space = readRDS(input_paths[["draws"]])
features = readRDS(input_paths[["features"]])$full
diagnostics = read_analysis_table(input_paths[["diagnostics"]])[specification == "full"]
positions = read_analysis_table(input_paths[["positions"]])
probabilities = read_analysis_table(input_paths[["probabilities"]])
observed = space$observed
actors = rownames(observed)
years = colnames(observed)
periods = length(years)
n = length(actors)
eligible = observed
eligible[, 1] = FALSE
eligible[, -1] = observed[, -1] & observed[, -periods]
latent_rank = dim(space$raw_positions[[1]])[2]
t_final = match(final_year, years)
stopifnot(
	latent_rank == 2L,
	all(eligible[focal, t_final]),
	all(vapply(features, function(f) f$prior$snap_kappa == kappa_fitted, TRUE))
)
# posterior means of the global parameters from the released diagnostics
global_parameters = c("sigma_uv", "pi_snap", "rho_uv", "va", "ve", "intercept")
parameter_means = diagnostics[family == "continuous" & parameter %in% global_parameters, .(parameter, mean, sd)]
sigma_mean = parameter_means[parameter == "sigma_uv", mean]
pi_mean = parameter_means[parameter == "pi_snap", mean]
rho_mean = parameter_means[parameter == "rho_uv", mean]

break_log_odds = function(adjusted_square, length_square, sigma, pi, kappa) {
	adjusted_square / (2 * sigma^2) - length_square / (2 * kappa^2) + qlogis(pi) + latent_rank * log(sigma / kappa)
}

# displacement at which the log-odds cross zero
break_threshold = function(length_square, sigma, pi, kappa) {
	sqrt(
		2 * sigma^2 * (length_square / (2 * kappa^2) - qlogis(pi) - latent_rank * log(sigma / kappa))
	)
}

shift_log_odds = function(log_odds, length_square, kappa) {
	log_odds + length_square / (2 * kappa_fitted^2) - length_square / (2 * kappa^2) +
		latent_rank * log(kappa_fitted / kappa)
}

chain_geometry = function(coordinates, rho) {
	draws = dim(coordinates)[4]
	length_square = raw_square = adjusted_square = array(0, c(n, periods, draws), dimnames = list(actors, years, NULL))
	for (t in seq_len(periods)) {
		for (k in seq_len(latent_rank)) {
			current = matrix(coordinates[, k, t, ], n, draws)
			length_square[, t, ] = length_square[, t, ] + current^2
			if (t == 1L) next
			previous = matrix(coordinates[, k, t - 1L, ], n, draws)
			raw_square[, t, ] = raw_square[, t, ] + (current - previous)^2
			adjusted_square[, t, ] = adjusted_square[, t, ] + (current - sweep(previous, 2, rho, "*"))^2
		}
	}
	list(
		length_square = length_square,
		raw_square = raw_square,
		adjusted_square = adjusted_square
	)
}

partner_change = function(coordinates, actor, t) {
	draws = dim(coordinates)[4]
	partners = which(observed[, t] & observed[, t - 1L] & actors != actor)
	self = match(actor, actors)
	inner = function(s) {
		value = matrix(0, length(partners), draws)
		for (k in seq_len(latent_rank)) {
			partner_coordinates = matrix(coordinates[partners, k, s, ], length(partners), draws)
			value = value + sweep(partner_coordinates, 2, coordinates[self, k, s, ], "*")
		}
		value
	}
	list(partners = actors[partners], delta = inner(t) - inner(t - 1L))
}

# uniform and directional components of a partner-change matrix, by draw
decompose_change = function(change) {
	delta = change$delta
	group = fifelse(
		change$partners %in% allies,
		"allies",
		fifelse(change$partners %in% rivals, "rivals", "others")
	)
	uniform = colMeans(delta)
	directional = sweep(delta, 2, uniform, "-")
	ss_uniform = nrow(delta) * uniform^2
	ss_directional = colSums(directional^2)
	out = data.table(
		partners = nrow(delta),
		uniform_mean = uniform,
		directional_sd = sqrt(ss_directional / (nrow(delta) - 1)),
		ss_share_uniform = ss_uniform / (ss_uniform + ss_directional),
		ss_share_directional = ss_directional / (ss_uniform + ss_directional)
	)
	for (g in c("allies", "rivals", "others")) {
		rows = group == g
		set(out, j = paste0("partners_", g), value = sum(rows))
		set(
			out,
			j = paste0("total_mean_", g),
			value = colMeans(delta[rows, , drop = FALSE])
		)
		set(
			out,
			j = paste0("directional_mean_", g),
			value = colMeans(directional[rows, , drop = FALSE])
		)
	}
	out
}

# rank of the final transition among a country's own transitions, by draw
final_rank = function(values, actor) {
	keep = which(eligible[actor, ])
	apply(
		values[actor, keep, ],
		2,
		function(x) rank(-x, ties.method = "min")[length(keep)]
	)
}

# draws of several chains joined along the draw dimension
bind_draws = function(pieces) {
	draws = vapply(pieces, function(x) dim(x)[3], 1L)
	out = array(
		NA_real_,
		c(dim(pieces[[1]])[1:2], sum(draws)),
		dimnames = list(dimnames(pieces[[1]])[[1]], dimnames(pieces[[1]])[[2]], NULL)
	)
	start = 0L
	for (piece in pieces) {
		out[, , start + seq_len(dim(piece)[3])] = piece
		start = start + dim(piece)[3]
	}
	out
}

summarize_draws = function(x, prefix) {
	q = quantile(x, interval, names = FALSE)
	out = data.table(mean(x), sd(x), q[1], q[2])
	setnames(out, paste(prefix, c("mean", "sd", "lower", "upper"), sep = "_"))
}

summarize_rank = function(x, prefix) {
	q = quantile(x, interval, names = FALSE)
	out = data.table(mean(x), mean(x == 1), q[1], q[2])
	setnames(
		out,
		paste(prefix, c("draw_mean", "share_first", "lower", "upper"), sep = "_")
	)
}

chain_results = lapply(
	names(space$raw_positions),
	function(seed) {
		chain = which(vapply(features, function(f) f$seed == as.integer(seed), TRUE))
		stopifnot(length(chain) == 1L)
		continuous = features[[chain]]$continuous
		coordinates = space$raw_positions[[seed]]
		draws = dim(coordinates)[4]
		stopifnot(
			nrow(continuous) == draws,
			identical(dimnames(coordinates)[[1]], actors),
			identical(dimnames(coordinates)[[3]], years)
		)
		sigma = continuous[, "sigma_uv"]
		pi = continuous[, "pi_snap"]
		rho = continuous[, "rho_uv"]
		geometry = chain_geometry(coordinates, rho)
		difference = max(
			vapply(
				2:periods,
				function(t) {
					keep = which(eligible[, t])
					id = paste(actors[keep], years[t], sep = "_")
					recomputed = break_log_odds(
						matrix(geometry$adjusted_square[keep, t, ], length(keep), draws),
						matrix(geometry$length_square[keep, t, ], length(keep), draws),
						rep(sigma, each = length(keep)),
						rep(pi, each = length(keep)),
						kappa_fitted
					)
					max(abs(t(recomputed) - features[[chain]]$actor_log_odds[, id, drop = FALSE]))
				},
				1
			)
		)
		raw = sqrt(geometry$raw_square[focal, , , drop = FALSE])
		adjusted = sqrt(geometry$adjusted_square[focal, , , drop = FALSE])
		sigma_units = sweep(adjusted, 3, sigma, "/")
		length_final = matrix(
			geometry$length_square[focal, t_final, ],
			length(focal),
			draws,
			dimnames = list(focal, NULL)
		)
		member_length = colMeans(geometry$length_square[observed[, t_final], t_final, ])
		# posterior mean of the squared length of every country-year in this chain
		length_mean = matrix(
			rowMeans(matrix(geometry$length_square, n * periods, draws)),
			n,
			periods,
			dimnames = list(actors, years)
		)
		ranks = list(
			adjusted = t(vapply(focal, function(actor) final_rank(adjusted, actor), numeric(draws))),
			raw = t(vapply(focal, function(actor) final_rank(raw, actor), numeric(draws)))
		)
		id = paste(focal, final_year, sep = "_")
		base = features[[chain]]$actor_log_odds[, id, drop = FALSE]
		recomputed_final = t(
			break_log_odds(
				geometry$adjusted_square[focal, t_final, ],
				length_final,
				rep(sigma, each = length(focal)),
				rep(pi, each = length(focal)),
				kappa_fitted
			)
		)
		colnames(recomputed_final) = id
		shifted = lapply(
			kappa_values,
			function(kappa) shift_log_odds(base, t(length_final), kappa)
		)
		names(shifted) = as.character(kappa_values)
		by_year = rbindlist(
			lapply(
				2:periods,
				function(t) {
					cbind(
						year = as.integer(years[t]),
						decompose_change(partner_change(coordinates, "USA", t))
					)
				}
			)
		)
		list(
			seed = as.integer(seed),
			chain = chain,
			draws = draws,
			difference = difference,
			raw = raw,
			adjusted = adjusted,
			sigma_units = sigma_units,
			length_final = length_final,
			member_length = member_length,
			length_mean = length_mean,
			ranks = ranks,
			sigma = sigma,
			pi = pi,
			rho = rho,
			base = base,
			recomputed_final = recomputed_final,
			shifted = shifted,
			by_year = by_year,
			threshold = break_threshold(length_final["USA", ], sigma, pi, kappa_fitted)
		)
	}
)
# pooled draws across chains
raw_all = bind_draws(lapply(chain_results, `[[`, "raw"))
adjusted_all = bind_draws(lapply(chain_results, `[[`, "adjusted"))
sigma_units_all = bind_draws(lapply(chain_results, `[[`, "sigma_units"))
rank_adjusted = do.call(cbind, lapply(chain_results, function(x) x$ranks$adjusted))
rank_raw = do.call(cbind, lapply(chain_results, function(x) x$ranks$raw))
length_final_all = do.call(cbind, lapply(chain_results, `[[`, "length_final"))
threshold_all = unlist(lapply(chain_results, `[[`, "threshold"))
base_all = do.call(rbind, lapply(chain_results, `[[`, "base"))
recomputed_all = do.call(rbind, lapply(chain_results, `[[`, "recomputed_final"))
member_length_all = unlist(lapply(chain_results, `[[`, "member_length"))
draws_by_chain = vapply(chain_results, `[[`, 1L, "draws")
weighted_lengths = Map(function(x, weight) x$length_mean * weight, chain_results, draws_by_chain)
length_mean_all = Reduce(`+`, weighted_lengths) / sum(draws_by_chain)

observed_lengths = data.table(
	actor = rep(actors, periods),
	year = rep(as.integer(years), each = n),
	is_observed = as.vector(observed),
	length_square = as.vector(length_mean_all)
)[is_observed == TRUE]
final_lengths = observed_lengths[year == as.integer(final_year)]
smallest = observed_lengths[which.min(length_square)]
largest = observed_lengths[which.max(length_square)]
smallest_final = final_lengths[which.min(length_square)]
largest_final = final_lengths[which.max(length_square)]
# posterior mean log-odds and 2025 rank from the released diagnostics
log_odds_table = diagnostics[family == "actor_log_odds", .(id = parameter, log_odds_mean = mean, log_odds_sd = sd)]
log_odds_table[, c("actor", "year") := tstrsplit(id, "_")]
log_odds_table[, year := as.integer(year)]
log_odds_table[, rank_log_odds := frank(-log_odds_mean, ties.method = "min"), by = actor]

by_year = rbindlist(
	lapply(
		focal,
		function(actor) {
			rbindlist(
				lapply(
					which(eligible[actor, ]),
					function(t) {
						cbind(
							data.table(actor = actor, year = as.integer(years[t])),
							summarize_draws(raw_all[actor, t, ], "raw"),
							summarize_draws(adjusted_all[actor, t, ], "adjusted"),
							summarize_draws(sigma_units_all[actor, t, ], "sigma_units")
						)
					}
				)
			)
		}
	)
)
by_year = merge(
	by_year,
	log_odds_table[, .(actor, year, log_odds_mean, rank_log_odds)],
	by = c("actor", "year"),
	all.x = TRUE
)
by_year = merge(
	by_year,
	probabilities[, .(actor, year, probability)],
	by = c("actor", "year"),
	all.x = TRUE
)
by_year[, rank_adjusted_posterior_mean := frank(-adjusted_mean, ties.method = "min"), by = actor]
by_year[, rank_raw_posterior_mean := frank(-raw_mean, ties.method = "min"), by = actor]
by_year[, transitions := .N, by = actor]
setcolorder(by_year, c("actor", "year", "transitions"))
fwrite(
	by_year[order(match(actor, focal), year)],
	file.path(threshold_root, "displacements_by_year.csv")
)

aligned = dcast(
	positions[actor %in% focal & year %in% c(2024L, 2025L), .(actor, year, x, y)],
	actor ~ year,
	value.var = c("x", "y")
)
aligned[, aligned_positions_displacement := sqrt((x_2025 - x_2024)^2 + (y_2025 - y_2024)^2)]
chain_means = rbindlist(
	lapply(
		chain_results,
		function(x) {
			data.table(
				chain = x$chain,
				actor = focal,
				adjusted = rowMeans(x$adjusted[, t_final, ]),
				raw = rowMeans(x$raw[, t_final, ])
			)
		}
	)
)
chain_spread = chain_means[, .(
	adjusted_chain_min = min(adjusted),
	adjusted_chain_max = max(adjusted),
	raw_chain_min = min(raw),
	raw_chain_max = max(raw)
), by = actor]
final = by_year[year == as.integer(final_year)]
rank_summary = rbindlist(
	lapply(
		focal,
		function(actor) {
			cbind(
				actor = actor,
				summarize_rank(rank_adjusted[actor, ], "rank_adjusted"),
				summarize_rank(rank_raw[actor, ], "rank_raw")
			)
		}
	)
)
final = merge(final, rank_summary, by = "actor")
final = merge(final, chain_spread, by = "actor")
final = merge(
	final,
	aligned[, .(actor, aligned_positions_displacement)],
	by = "actor"
)
final[, length_square_mean := rowMeans(length_final_all)[actor]]
final[, adjusted_over_sigma_mean := adjusted_mean / sigma_mean]
final[, threshold_at_own_length := break_threshold(length_square_mean, sigma_mean, pi_mean, kappa_fitted)]
final[, adjusted_over_threshold := adjusted_mean / threshold_at_own_length]
final[, raw_share_of_us := raw_mean / raw_mean[actor == "USA"]]
final[, adjusted_share_of_us := adjusted_mean / adjusted_mean[actor == "USA"]]
final[, aligned_share_of_us := aligned_positions_displacement / aligned_positions_displacement[actor == "USA"]]
final[, group := fifelse(actor == "USA", "United States", "Ally")]
setcolorder(final, c("actor", "group", "year", "transitions"))
fwrite(
	final[order(match(actor, focal))],
	file.path(threshold_root, "displacements.csv")
)
# one reference length for the threshold table
reference_row = function(label, length_square, actor = NA_character_, year = NA_integer_) {
	data.table(
		length_reference = label,
		length_actor = actor,
		length_year = year,
		length_square = length_square
	)
}

length_references = rbindlist(
	list(
		reference_row("zero", 0),
		reference_row("unit", 1),
		reference_row(
			"usa_2025_mean",
			mean(length_final_all["USA", ]),
			"USA",
			as.integer(final_year)
		),
		reference_row(
			"member_mean_2025",
			mean(member_length_all),
			year = as.integer(final_year)
		),
		reference_row(
			"min_observed_2025",
			smallest_final$length_square,
			smallest_final$actor,
			smallest_final$year
		),
		reference_row(
			"max_observed_2025",
			largest_final$length_square,
			largest_final$actor,
			largest_final$year
		),
		reference_row(
			"min_observed_all_years",
			smallest$length_square,
			smallest$actor,
			smallest$year
		),
		reference_row(
			"max_observed_all_years",
			largest$length_square,
			largest$actor,
			largest$year
		)
	)
)
threshold = rbindlist(
	lapply(
		c(kappa_values, kappa_fitted),
		function(kappa) {
			cbind(
				kappa = kappa,
				length_references,
				sigma = sigma_mean,
				pi = pi_mean,
				rho = rho_mean
			)
		}
	)
)
threshold = threshold[order(kappa, match(length_reference, length_references$length_reference))]
threshold[, threshold_units := break_threshold(length_square, sigma, pi, kappa)]
threshold[, threshold_sigma_units := threshold_units / sigma]
threshold[, break_density_relative_to_origin := exp(-length_square / (2 * kappa^2))]
threshold_draws = summarize_draws(threshold_all, "threshold_posterior")
posterior_row = threshold$kappa == kappa_fitted & threshold$length_reference == "usa_2025_mean"
for (column in names(threshold_draws)) {
	set(
		threshold,
		j = column,
		value = fifelse(posterior_row, threshold_draws[[column]], NA_real_)
	)
}

fwrite(threshold, file.path(threshold_root, "threshold.csv"))

usa_id = paste("USA", final_year, sep = "_")
usa_adjusted = adjusted_all["USA", t_final, ]
usa_raw = raw_all["USA", t_final, ]
usa_length = length_final_all["USA", ]
usa_aligned = aligned[actor == "USA", aligned_positions_displacement]
log_odds_quantities = c(
	"diagnostics_posterior_mean",
	"recomputed_posterior_mean",
	"stored_features_posterior_mean",
	"at_posterior_means_adjusted_displacement",
	"at_posterior_means_raw_displacement",
	"at_posterior_mean_squared_adjusted_displacement",
	"at_posterior_means_aligned_displacement"
)
other_quantities = c(
	"posterior_mean_probability",
	"max_abs_difference_all_country_years",
	"posterior_mean_adjusted_displacement",
	"posterior_mean_raw_displacement",
	"posterior_mean_squared_adjusted_displacement",
	"posterior_mean_length_square",
	"aligned_positions_displacement"
)
log_odds_check = data.table(
	quantity = c(log_odds_quantities, other_quantities),
	value = c(
		log_odds_table[id == usa_id, log_odds_mean],
		mean(recomputed_all[, usa_id]),
		mean(base_all[, usa_id]),
		break_log_odds(
			mean(usa_adjusted)^2,
			mean(usa_length),
			sigma_mean,
			pi_mean,
			kappa_fitted
		),
		break_log_odds(
			mean(usa_raw)^2,
			mean(usa_length),
			sigma_mean,
			pi_mean,
			kappa_fitted
		),
		break_log_odds(
			mean(usa_adjusted^2),
			mean(usa_length),
			sigma_mean,
			pi_mean,
			kappa_fitted
		),
		break_log_odds(
			usa_aligned^2,
			mean(usa_length),
			sigma_mean,
			pi_mean,
			kappa_fitted
		),
		mean(plogis(base_all[, usa_id])),
		max(vapply(chain_results, `[[`, 1, "difference")),
		mean(usa_adjusted),
		mean(usa_raw),
		mean(usa_adjusted^2),
		mean(usa_length),
		usa_aligned
	)
)
log_odds_check[, probability := fifelse(quantity %in% log_odds_quantities, plogis(value), NA_real_)]
fwrite(log_odds_check, file.path(threshold_root, "log_odds_check.csv"))

components_draws = rbindlist(
	lapply(
		chain_results,
		function(x) {
			cbind(chain = x$chain, x$by_year[year == as.integer(final_year)])
		}
	)
)
component_names = setdiff(names(components_draws), c("chain", "year"))
components = rbindlist(
	lapply(
		component_names,
		function(column) {
			x = components_draws[[column]]
			chain_values = components_draws[, mean(get(column)), by = chain][, V1]
			cbind(
				component = column,
				summarize_draws(x, "posterior"),
				chain_min = min(chain_values),
				chain_max = max(chain_values)
			)
		}
	)
)
fwrite(components, file.path(decomposition_root, "us_2025_components.csv"))
by_year_draws = rbindlist(
	lapply(chain_results, function(x) cbind(chain = x$chain, x$by_year))
)
components_by_year = by_year_draws[, .(
	partners = partners[1],
	uniform_mean = mean(uniform_mean),
	uniform_lower = quantile(uniform_mean, interval[1]),
	uniform_upper = quantile(uniform_mean, interval[2]),
	directional_sd = mean(directional_sd),
	ss_share_uniform = mean(ss_share_uniform),
	ss_share_uniform_lower = quantile(ss_share_uniform, interval[1]),
	ss_share_uniform_upper = quantile(ss_share_uniform, interval[2]),
	total_mean_allies = mean(total_mean_allies),
	total_mean_rivals = mean(total_mean_rivals),
	total_mean_others = mean(total_mean_others),
	directional_mean_allies = mean(directional_mean_allies),
	directional_mean_rivals = mean(directional_mean_rivals),
	directional_mean_others = mean(directional_mean_others)
), by = year]
components_by_year[, rank_ss_share_uniform := frank(-ss_share_uniform, ties.method = "min")]
components_by_year[, rank_abs_uniform_mean := frank(-abs(uniform_mean), ties.method = "min")]
fwrite(
	components_by_year[order(year)],
	file.path(decomposition_root, "us_components_by_year.csv")
)

kappa_table = rbindlist(
	lapply(
		chain_results,
		function(x) {
			out = data.table(
				chain = as.character(x$chain),
				seed = x$seed,
				actor = focal,
				probability_kappa_2 = colMeans(plogis(x$base))
			)
			for (kappa in kappa_values) {
				set(
					out,
					j = paste0("probability_kappa_", kappa),
					value = colMeans(plogis(x$shifted[[as.character(kappa)]]))
				)
			}
			out
		}
	)
)
pooled = kappa_table[, lapply(.SD, mean), by = actor, .SDcols = patterns("^probability")]
pooled[, c("chain", "seed") := .("all", NA_integer_)]
kappa_table = rbind(kappa_table, pooled, use.names = TRUE)
reported = probabilities[year == as.integer(final_year), .(actor, reported_probability = probability)]
kappa_table = merge(kappa_table, reported, by = "actor")
kappa_table[, length_square_mean := rowMeans(length_final_all)[actor]]
for (kappa in c(kappa_values, kappa_fitted)) {
	column = paste0("threshold_kappa_", kappa)
	set(
		kappa_table,
		j = column,
		value = break_threshold(kappa_table$length_square_mean, sigma_mean, pi_mean, kappa)
	)
}

for (kappa in kappa_values) {
	column = paste0("log_odds_shift_kappa_", kappa)
	set(
		kappa_table,
		j = column,
		value = shift_log_odds(0, kappa_table$length_square_mean, kappa)
	)
}

setcolorder(
	kappa_table,
	c(
		"actor",
		"chain",
		"seed",
		"reported_probability",
		"probability_kappa_0.5",
		"probability_kappa_1",
		"probability_kappa_2",
		"probability_kappa_4"
	)
)
fwrite(
	kappa_table[order(match(actor, focal), chain)],
	file.path(kappa_root, "reset_probabilities_by_kappa.csv")
)
# provenance beside each output set
alignment = lapply(
	chain_results,
	function(x) {
		list(
			seed = x$seed,
			chain_features_index = x$chain,
			draws = x$draws,
			max_abs_log_odds_difference = x$difference
		)
	}
)
inputs = lapply(
	input_paths,
	function(path) list(path = path, sha256 = digest(file = path, algo = "sha256"))
)
provenance = list(
	script = script_path,
	script_sha256 = digest(file = script_path, algo = "sha256"),
	created_utc = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
	inputs = inputs,
	posterior_means = as.list(setNames(parameter_means$mean, parameter_means$parameter)),
	kappa_fitted = kappa_fitted,
	kappa_values = kappa_values,
	latent_rank = latent_rank,
	interval = interval,
	allies = allies,
	rivals = rivals,
	alignment_check = alignment,
	R = R.version.string
)
for (path in c(threshold_root, decomposition_root, kappa_root)) {
	write_json(
		provenance,
		file.path(path, "provenance.json"),
		pretty = TRUE,
		auto_unbox = TRUE,
		digits = NA
	)
}
