# compare observed and fitted trade relationships

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
specifications = c("current_baseline", "partner_reallocation_terminal")

relationship_change = function(fit, seed, specification) {
	y = fit$Y
	u = fit$U_draws
	actors = dimnames(y)[[1]]
	years = dimnames(y)[[3]]
	us = match("USA", actors)
	rows = list()
	for (k in 2:length(years)) {
		now = years[k]
		before = years[k - 1]
		seen = is.finite(y[us, , now]) & is.finite(y[us, , before])
		seen[us] = FALSE
		ally = seen & actors %in% allies
		other = seen & !(actors %in% allies)
		stopifnot(sum(ally) >= 3, sum(other) >= 3)
		observed = (mean(y[us, ally, now]) - mean(y[us, other, now])) -
			(mean(y[us, ally, before]) - mean(y[us, other, before]))
		model = vapply(
			seq_len(dim(u)[4]),
			function(d) {
				r_now = drop(u[, , now, d] %*% u[us, , now, d])
				r_before = drop(u[, , before, d] %*% u[us, , before, d])
				(mean(r_now[ally]) - mean(r_now[other])) - (mean(r_before[ally]) - mean(r_before[other]))
			},
			numeric(1)
		)
		rows[[length(rows) + 1L]] = data.frame(
			specification = specification,
			seed = seed,
			year = as.integer(now),
			allies = sum(ally),
			other_partners = sum(other),
			observed_change = observed,
			draw = seq_along(model),
			model_change = model
		)
	}
	do.call(rbind, rows)
}

draws = list()
for (specification in specifications) {
	replay = file.path(root, "results/models/trade_space_replay", specification)
	verified = read_json(file.path(replay, "verification.json"), simplifyVector = TRUE)
	stopifnot(isTRUE(verified$passed))
	runs = read.csv(file.path(replay, "fit_manifest.csv"))
	for (i in seq_len(nrow(runs))) {
		at = which(verified$fits$seed == runs$seed[i])
		stopifnot(
			length(at) == 1L,
			digest(file = runs$fit_file[i], algo = "sha256") == verified$fits$sha256[at]
		)
		fit = readRDS(runs$fit_file[i])
		draws[[length(draws) + 1L]] = relationship_change(fit, runs$seed[i], specification)
		rm(fit)
		gc(FALSE)
	}
}

draws = do.call(rbind, draws)
summary = do.call(
	rbind,
	lapply(
		split(draws, list(draws$specification, draws$year), drop = TRUE),
		function(x) {
			stopifnot(
				length(unique(x$observed_change)) == 1L,
				length(unique(x$allies)) == 1L
			)
			data.frame(
				specification = x$specification[1],
				year = x$year[1],
				allies = x$allies[1],
				other_partners = x$other_partners[1],
				observed_change = x$observed_change[1],
				model_mean = mean(x$model_change),
				model_lower = unname(quantile(x$model_change, .025)),
				model_upper = unname(quantile(x$model_change, .975)),
				draws = nrow(x)
			)
		}
	)
)
summary = summary[order(summary$specification, summary$year), ]
output = file.path(root, "results/models/trade/relations")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(
	summary,
	file.path(output, "us_ally_relationship.csv"),
	row.names = FALSE
)
write_json(
	list(
		quantity = "Change from the preceding year in the mean logged trade of the United States with the fourteen allies minus its mean logged trade with all other partners observed in both years. The observed change uses the model input; the model change uses the inner products of the latent positions in each posterior draw.",
		specifications = specifications,
		sources = setNames(
			lapply(
				c(
					"code/18_trade_relations.R",
					"code/17_trade_replay.R",
					file.path(
						"results/models/trade_space_replay",
						specifications,
						"verification.json"
					),
					file.path(
						"results/models/trade_space_replay",
						specifications,
						"fit_manifest.csv"
					)
				),
				function(p) digest(file = file.path(root, p), algo = "sha256")
			),
			c(
				"code/18_trade_relations.R",
				"code/17_trade_replay.R",
				file.path(
					"results/models/trade_space_replay",
					specifications,
					"verification.json"
				),
				file.path(
					"results/models/trade_space_replay",
					specifications,
					"fit_manifest.csv"
				)
			)
		),
		output_sha256 = digest(file = file.path(output, "us_ally_relationship.csv"), algo = "sha256")
	),
	file.path(output, "provenance.json"),
	pretty = TRUE,
	auto_unbox = TRUE
)
