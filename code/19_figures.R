# draw the manuscript and appendix figures

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
suppressPackageStartupMessages(library(ggplot2))
source(file.path(dirname(normalizePath(script_file)), "00_plot_style.R"))
paths = pipeline_paths()
graphics_root = file.path(paths$output, "figures")
dir.create(graphics_root, recursive = TRUE, showWarnings = FALSE)

write_latent_space = function(paths, manifest_path) {
	root = dirname(dirname(normalizePath(script_file)))
	runs = read_analysis_table(manifest_path)
	if ("domain" %in% names(runs)) runs = runs[domain == "un"]
	if ("specification" %in% names(runs)) runs = runs[specification == "full"]
	plan = read_analysis_table(file.path(root, "data/model_chains.csv"))[domain == "un" & specification == "full"]
	stopifnot(
		nrow(runs) == 4L,
		all(runs$complete),
		setequal(runs$seed, plan$seed),
		!anyDuplicated(runs$seed)
	)
	runs = runs[order(seed)]
	out = file.path(paths$output, "models/un/space")
	dir.create(out, recursive = TRUE, showWarnings = FALSE)
	selected = as.character(c(1989, 1992, 2000, 2002, 2016, 2020, 2024, 2025))
	first = readRDS(runs$fit_file[1])
	n = dim(first$Y)[1]
	years = dimnames(first$Y)[[3]]
	actors = dimnames(first$Y)[[1]]
	observed = vapply(
		seq_along(years),
		function(t) rowSums(is.finite(first$Y[, , t])) > 0,
		logical(n)
	)
	dimnames(observed) = list(actors, years)
	selected_index = match(selected, years)
	stopifnot(!anyNA(selected_index))
	reference = first$U_draws[, , "2024", 1]
	reference = reference %*% eigen(crossprod(reference[observed[, "2024"], ]), symmetric = TRUE)$vectors
	if (reference[match("USA", actors), 1] > 0) reference[, 1] = -reference[, 1]
	if (reference[match("USA", actors), 2] < 0) reference[, 2] = -reference[, 2]
	rm(first)
	gc(FALSE)
	position_chains = list()
	distance_chains = list()
	checks = list()
	distance_counts = NULL
	for (j in seq_len(nrow(runs))) {
		f = readRDS(runs$fit_file[j])
		job = plan[seed == runs$seed[j]]
		cfg = f$revision_build$config
		stopifnot(
			isTRUE(f$complete),
			identical(digest::digest(f$Y, algo = "sha256"), job$input_sha256),
			cfg$seed == job$seed,
			cfg$nscan == job$sampling_sweeps,
			cfg$burn == job$warmup_sweeps,
			cfg$odens == job$thin,
			identical(
				digest::digest(cfg$initial_state, algo = "sha256"),
				job$initial_state_sha256
			),
			identical(
				digest::digest(cfg$frozen_shifts, algo = "sha256"),
				job$shifts_sha256
			),
			cfg$temporal_tau == job$temporal_tau,
			cfg$indicator_every == job$indicator_every
		)
		draws = dim(f$U_draws)[4]
		stopifnot(draws == job$sampling_sweeps / job$thin)
		distance = array(
			NA_real_,
			c(draws, length(years), 4),
			dimnames = list(
				NULL,
				years,
				c("US to allies", "Among allies", "US to rivals", "Europe to Indo-Pacific")
			)
		)
		count = list()
		for (t in seq_along(years)) {
			available = actors[observed[, t]]
			allies = intersect(c(actor_groups$allies, if (years[t] == "1989") "GER"), available)
			rivals = intersect(c(actor_groups$rivals, "SUN"), available)
			europe = intersect(
				c("FRA", "DEU", "GER", "GBR", "ITA", "NLD", "ESP", "BEL", "NOR", "DNK", "POL"),
				available
			)
			indo = intersect(c("JPN", "KOR", "AUS"), available)
			pair = list(
				if ("USA" %in% available) cbind("USA", allies) else matrix(character(), 0, 2),
				if (length(allies) > 1) t(combn(allies, 2)) else matrix(character(), 0, 2),
				if ("USA" %in% available) cbind("USA", rivals) else matrix(character(), 0, 2),
				as.matrix(expand.grid(europe, indo, stringsAsFactors = FALSE))
			)
			for (g in seq_along(pair)) {
				ij = pair[[g]]
				stopifnot(nrow(ij) > 0)
				a = match(ij[, 1], actors)
				b = match(ij[, 2], actors)
				dx = matrix(f$U_draws[a, 1, t, ] - f$U_draws[b, 1, t, ], nrow(ij), draws)
				dy = matrix(f$U_draws[a, 2, t, ] - f$U_draws[b, 2, t, ], nrow(ij), draws)
				distance[, t, g] = colMeans(sqrt(dx^2 + dy^2))
				count[[length(count) + 1L]] = data.table(
					year = as.integer(years[t]),
					relationship = dimnames(distance)[[3]][g],
					pairs = nrow(ij),
					countries = length(unique(c(ij)))
				)
			}
		}
		position_chains[[as.character(job$seed)]] = f$U_draws
		distance_chains[[as.character(job$seed)]] = distance
		cnt = rbindlist(count)
		if (is.null(distance_counts)) distance_counts = cnt else stopifnot(identical(cnt, distance_counts))
		checks[[j]] = data.table(
			seed = job$seed,
			draws = draws,
			fit_sha256 = digest::digest(file = runs$fit_file[j], algo = "sha256")
		)
		rm(f, distance)
		gc(FALSE)
	}
	saveRDS(
		list(
			raw_positions = position_chains,
			distances = distance_chains,
			observed = observed,
			display_years = selected,
			pair_counts = distance_counts,
			reference = reference
		),
		file.path(out, "draws.rds"),
		compress = FALSE
	)
	summarize_latent_space(paths)
	summarize_spatial_groups(paths)
	fwrite(rbindlist(checks), file.path(out, "replay_checks.csv"))
	jsonlite::write_json(
		list(
			created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
			seeds = runs$seed,
			draws_per_chain = draws,
			display_years = as.integer(selected),
			retained_years = as.integer(years),
			alignment = "Each year is rotated/reflected to the same fixed 2024 reference on countries observed in both years; no translation or rescaling. These display transformations preserve all within-year inner products and distances and do not enter estimation or reset classification. Reference is the lowest-seed first retained 2024 draw in principal axes, with USA negative on the first axis and positive on the second.",
			ellipses = "Empirical 95th percentile of squared Mahalanobis radius about the pooled aligned mean and covariance; elliptical summaries, not independent significance tests.",
			distances = "Posterior mean of pairwise Euclidean distances, calculated within each draw before averaging; available countries, with counts retained. Not the distance between posterior-mean coordinates.",
			inputs = lapply(
				c(file.path(root, "data/model_chains.csv"), manifest_path),
				function(p) list(path = p, sha256 = digest::digest(file = p, algo = "sha256"))
			)
		),
		file.path(out, "provenance.json"),
		pretty = TRUE,
		auto_unbox = TRUE
	)
	invisible(out)
}

align_latent_space_draws = function(archive) {
	reference = archive$reference
	observed = archive$observed
	selected = archive$display_years
	actors = dimnames(archive$raw_positions[[1]])[[1]]
	n = length(actors)
	chains = list()
	checks = list()
	align = function(raw, keep) {
		v = svd(crossprod(raw[keep, ], reference[keep, ]))
		raw %*% (v$u %*% t(v$v))
	}
	test_rotation = matrix(c(cos(.37), sin(.37), -sin(.37), cos(.37)), 2, 2) %*% diag(c(1, -1))
	for (seed in names(archive$raw_positions)) {
		u = archive$raw_positions[[seed]]
		draws = dim(u)[4]
		a = array(
			NA_real_,
			c(n, 2, length(selected), draws),
			dimnames = list(actors, c("x", "y"), selected, NULL)
		)
		gram_error = invariance_error = 0
		for (t in selected) {
			keep = observed[, t] & observed[, "2024"]
			stopifnot(sum(keep) > 2)
			for (d in seq_len(draws)) {
				raw = u[, , t, d]
				z = align(raw, keep)
				a[, , t, d] = z
				if (d %in% c(1L, draws)) {
					gram_error = max(gram_error, max(abs(tcrossprod(raw) - tcrossprod(z))))
					invariance_error = max(invariance_error, max(abs(z - align(raw %*% test_rotation, keep))))
				}
			}
		}
		stopifnot(gram_error < 1e-10, invariance_error < 1e-10)
		chains[[seed]] = a
		checks[[seed]] = data.table(
			seed = as.integer(seed),
			max_gram_error = gram_error,
			max_orientation_invariance_error = invariance_error
		)
	}
	list(positions = chains, checks = rbindlist(checks))
}

summarize_latent_space = function(paths) {
	out = file.path(paths$output, "models/un/space")
	archive = readRDS(file.path(out, "draws.rds"))
	aligned = align_latent_space_draws(archive)
	position_chains = aligned$positions
	fwrite(aligned$checks, file.path(out, "geometry_checks.csv"))
	distance_chains = archive$distances
	observed = archive$observed
	distance_counts = archive$pair_counts
	actors = dimnames(position_chains[[1]])[[1]]
	selected = archive$display_years
	years = dimnames(distance_chains[[1]])[[2]]
	diagnostic = function(x) {
		data.table(
			rhat = posterior::rhat(x),
			ess_bulk = posterior::ess_bulk(x),
			ess_tail = posterior::ess_tail(x)
		)
	}
	positions = list()
	coordinate_diagnostics = list()
	for (t in seq_along(selected)) {
		for (i in which(observed[, selected[t]])) {
			x = do.call(cbind, lapply(position_chains, function(a) a[i, 1, selected[t], ]))
			y = do.call(cbind, lapply(position_chains, function(a) a[i, 2, selected[t], ]))
			a = cbind(as.vector(x), as.vector(y))
			s = cov(a)
			radius = quantile(mahalanobis(a, colMeans(a), s), .95, names = FALSE)
			positions[[length(positions) + 1L]] = data.table(
				actor = actors[i],
				year = as.integer(selected[t]),
				x = mean(x),
				y = mean(y),
				var_x = s[1, 1],
				var_y = s[2, 2],
				cov_xy = s[1, 2],
				radius_squared_95 = radius
			)
			coordinate_diagnostics[[length(coordinate_diagnostics) + 1L]] = rbind(
				cbind(
					data.table(
						actor = actors[i],
						year = as.integer(selected[t]),
						coordinate = "x"
					),
					diagnostic(x)
				),
				cbind(
					data.table(
						actor = actors[i],
						year = as.integer(selected[t]),
						coordinate = "y"
					),
					diagnostic(y)
				)
			)
		}
	}
	distances = list()
	for (t in seq_along(years)) {
		for (g in 1:4) {
			x = do.call(cbind, lapply(distance_chains, function(a) a[, t, g]))
			ci = quantile(x, c(.025, .975), names = FALSE)
			distances[[length(distances) + 1L]] = cbind(
				data.table(
					year = as.integer(years[t]),
					relationship = dimnames(distance_chains[[1]])[[3]][g],
					mean = mean(x),
					lower = ci[1],
					upper = ci[2]
				),
				diagnostic(x)
			)
		}
	}
	fwrite(rbindlist(positions), file.path(out, "positions.csv"))
	fwrite(
		rbindlist(coordinate_diagnostics),
		file.path(out, "coordinate_diagnostics.csv")
	)
	d = merge(rbindlist(distances), distance_counts, by = c("year", "relationship"))
	d[, `:=`(
		rank_descending = frank(-mean, ties.method = "min"),
		rank_ascending = frank(mean, ties.method = "min")
	), by = relationship]
	fwrite(d, file.path(out, "distances.csv"))
	comparisons = list()
	history_years = which(as.integer(years) >= 1991 & as.integer(years) < 2025)
	for (g in 1:4) {
		for (baseline in c("2023", "2024")) {
			before = do.call(cbind, lapply(distance_chains, function(a) a[, baseline, g]))
			after = do.call(cbind, lapply(distance_chains, function(a) a[, "2025", g]))
			change = after - before
			ratio = after / before
			label = dimnames(distance_chains[[1]])[[3]][g]
			historical = do.call(rbind, lapply(distance_chains, function(a) a[, history_years, g]))
			counts = distance_counts[relationship == label & year >= 1991]
			ci = quantile(change, c(.025, .975), names = FALSE)
			ri = quantile(ratio, c(.025, .975), names = FALSE)
			comparisons[[length(comparisons) + 1L]] = cbind(
				data.table(
					relationship = label,
					baseline = as.integer(baseline),
					target = 2025,
					mean_change = mean(change),
					lower_change = ci[1],
					upper_change = ci[2],
					probability_increase = mean(change > 0),
					mean_ratio = mean(ratio),
					lower_ratio = ri[1],
					upper_ratio = ri[2],
					probability_more_than_double = mean(ratio > 2),
					probability_exceeds_all_since_1991 = mean(as.vector(after) > apply(historical, 1, max)),
					probability_below_all_since_1991 = mean(as.vector(after) < apply(historical, 1, min)),
					rank_mean_since_1991 = 1L + sum(d[relationship == label & year >= 1991 & year < 2025]$mean > mean(after)),
					historical_min_pairs = min(counts$pairs),
					historical_max_pairs = max(counts$pairs)
				),
				diagnostic(change)
			)
		}
	}
	fwrite(rbindlist(comparisons), file.path(out, "comparisons.csv"))
	invisible(out)
}

# group summaries retain the country maps' reference frame and original draws.
summarize_spatial_groups = function(paths) {
	out = file.path(paths$output, "models/un/space")
	a = readRDS(file.path(out, "draws.rds"))
	u = align_latent_space_draws(a)$positions
	groups = list(
		"United States" = "USA",
		"Selected allies" = c(actor_groups$allies, "GER"),
		"Selected rivals" = c(actor_groups$rivals, "SUN"),
		"Other selected countries" = actor_groups$hedgers
	)
	diagnose = function(x) {
		data.table(
			rhat = posterior::rhat(x),
			ess_bulk = posterior::ess_bulk(x),
			ess_tail = posterior::ess_tail(x)
		)
	}
	centroids = diagnostics = list()
	for (yr in a$display_years) {
		for (g in names(groups)) {
			members = intersect(groups[[g]], rownames(a$observed)[a$observed[, yr]])
			z = lapply(
				u,
				function(v) t(apply(v[members, , yr, , drop = FALSE], c(2, 4), mean))
			)
			xy = do.call(rbind, z)
			cv = cov(xy)
			q = quantile(mahalanobis(xy, colMeans(xy), cv), .95, names = FALSE)
			centroids[[length(centroids) + 1L]] = data.table(
				year = as.integer(yr),
				group = g,
				countries = length(members),
				members = paste(members, collapse = "|"),
				x = mean(xy[, 1]),
				y = mean(xy[, 2]),
				var_x = cv[1, 1],
				var_y = cv[2, 2],
				cov_xy = cv[1, 2],
				radius_squared_95 = q
			)
			for (k in 1:2) {
				diagnostics[[length(diagnostics) + 1L]] = cbind(
					data.table(year = as.integer(yr), group = g, coordinate = c("x", "y")[k]),
					diagnose(do.call(cbind, lapply(z, function(v) v[, k])))
				)
			}
		}
	}
	regions = list(
		"Europe" = c("FRA", "DEU", "GER", "GBR", "ITA", "NLD", "ESP", "BEL", "NOR", "DNK", "POL"),
		"Indo-Pacific" = c("JPN", "KOR", "AUS")
	)
	bridge = list()
	for (yr in dimnames(a$raw_positions[[1]])[[3]]) {
		available = rownames(a$observed)[a$observed[, yr]]
		europe = intersect(regions$Europe, available)
		indo = intersect(regions[["Indo-Pacific"]], available)
		pair = list(
			"US to Europe" = cbind("USA", europe),
			"US to Indo-Pacific" = cbind("USA", indo),
			"Europe to Indo-Pacific" = as.matrix(expand.grid(europe, indo, stringsAsFactors = FALSE))
		)
		for (country in intersect(c("JPN", "KOR", "AUS"), available)) {
			label = c(
				JPN = "Japan to Europe",
				KOR = "South Korea to Europe",
				AUS = "Australia to Europe"
			)[[country]]
			pair[[label]] = cbind(country, europe)
		}
		for (g in names(pair)) {
			ij = pair[[g]]
			stopifnot(all(a$observed[unique(as.vector(ij)), yr]))
			z = lapply(
				a$raw_positions,
				function(v) {
					nd = dim(v)[4]
					dx = matrix(v[ij[, 1], 1, yr, ] - v[ij[, 2], 1, yr, ], nrow(ij), nd)
					dy = matrix(v[ij[, 1], 2, yr, ] - v[ij[, 2], 2, yr, ], nrow(ij), nd)
					colMeans(sqrt(dx * dx + dy * dy))
				}
			)
			x = do.call(cbind, z)
			ci = quantile(x, c(.025, .975), names = FALSE)
			bridge[[length(bridge) + 1L]] = cbind(
				data.table(
					year = as.integer(yr),
					relationship = g,
					pairs = nrow(ij),
					mean = mean(x),
					lower = ci[1],
					upper = ci[2]
				),
				diagnose(x)
			)
		}
	}
	c = rbindlist(centroids)
	d = rbindlist(diagnostics)
	b = rbindlist(bridge)
	stopifnot(
		all(
			is.finite(d$rhat) & d$rhat <= 1.01 & d$ess_bulk >= 400 & d$ess_tail >= 400
		),
		all(
			is.finite(b$rhat) & b$rhat <= 1.01 & b$ess_bulk >= 400 & b$ess_tail >= 400
		)
	)
	fwrite(c, file.path(out, "centroids.csv"))
	fwrite(d, file.path(out, "centroid_diagnostics.csv"))
	fwrite(b, file.path(out, "regional_distances.csv"))
	p = read_analysis_table(file.path(out, "positions.csv"))
	check = list()
	for (i in seq_len(nrow(c))) {
		row = c[i]
		members = strsplit(row$members, "|", fixed = TRUE)[[1]]
		v = p[year == row$year & actor %in% members]
		error = max(abs(c(mean(v$x), mean(v$y)) - c(row$x, row$y)))
		stopifnot(nrow(v) == row$countries, error < 1e-12)
		check[[i]] = data.table(
			check = paste(row$year, row$group, "mean_of_country_means"),
			error = error,
			passed = TRUE
		)
	}
	d0 = read_analysis_table(file.path(out, "distances.csv"))
	r = b[relationship == "Europe to Indo-Pacific"]
	r0 = d0[relationship == "Europe to Indo-Pacific" & year %in% r$year]
	error = max(
		abs(
			as.matrix(r[, .(mean, lower, upper)]) - as.matrix(r0[, .(mean, lower, upper)])
		)
	)
	stopifnot(error < 1e-12)
	check[[length(check) + 1L]] = data.table(
		check = "Regional distances agree with existing history",
		error = error,
		passed = TRUE
	)
	fwrite(rbindlist(check), file.path(out, "group_checks.csv"))
	jsonlite::write_json(
		list(
			draws_per_chain = dim(a$raw_positions[[1]])[4],
			chains = length(u),
			source_sha256 = digest::digest(file = file.path(out, "draws.rds"), algo = "sha256"),
			centroid = "Equal-country mean of aligned positions within each draw, then posterior mean. Available members at each endpoint; group counts and identifiers retained.",
			regional_distances = "Equal-pair mean of raw Euclidean distances within each draw, then posterior mean and 95% interval. European group has ten countries in 2024-2025, with West Germany used where observed historically; Indo-Pacific group is Japan, South Korea and Australia, with only observed countries entering each year. Canada is included in the fourteen-allied comparisons but not either regional subgroup.",
			schematic = "Fixed node locations are a diagram, not latent coordinates. Edge labels report posterior mean pairwise distances and intervals; line widths do not encode tie strength."
		),
		file.path(out, "groups_provenance.json"),
		auto_unbox = TRUE,
		pretty = TRUE
	)
	invisible(out)
}

trim_position_map = function(path) {
	cropped = tempfile(fileext = ".pdf")
	on.exit(unlink(cropped), add = TRUE)
	stopifnot(nzchar(Sys.which("pdfcrop")))
	status = system2(
		"pdfcrop",
		c("--margins", "2", shQuote(path), shQuote(cropped)),
		stdout = FALSE,
		stderr = FALSE
	)
	stopifnot(
		status == 0L,
		file.exists(cropped),
		file.copy(cropped, path, overwrite = TRUE)
	)
	invisible(path)
}

draw_latent_space = function(paths, graphics_root) {
	root = file.path(paths$output, "models/un/space")
	p = read_analysis_table(file.path(root, "positions.csv"))
	colors = setNames(
		c(RColorBrewer::brewer.pal(9, "Set1")[c(1, 2, 3, 9)], "grey80"),
		c(
			"United States",
			"Selected allies",
			"Selected rivals",
			"Other selected countries",
			"Other countries"
		)
	)
	p[, category := "Other countries"]
	p[actor %in% actor_groups$hedgers, category := "Other selected countries"]
	p[actor %in% c(actor_groups$rivals, "SUN"), category := "Selected rivals"]
	p[actor %in% c(actor_groups$allies, "GER"), category := "Selected allies"]
	p[actor == "USA", category := "United States"]
	p[, category := factor(category, levels = names(colors))]
	country_names = c(
		USA = "United States",
		DEU = "Germany",
		JPN = "Japan",
		RUS = "Russia",
		CHN = "China",
		IND = "India",
		TUR = "Türkiye",
		BRA = "Brazil",
		POL = "Poland",
		GER = "West Germany",
		SUN = "Soviet Union"
	)
	focus = c("USA", "DEU", "JPN", "RUS", "CHN", "IND", "TUR", "BRA")
	p[, label := ifelse(actor %in% focus, unname(country_names[actor]), "")]
	style = paper_plot_theme(base_size = 10)
	ellipse = function(r) {
		v = matrix(c(r$var_x, r$cov_xy, r$cov_xy, r$var_y), 2, 2)
		e = eigen(v, symmetric = TRUE)
		theta = seq(0, 2 * pi, length.out = 81)
		xy = cbind(cos(theta), sin(theta)) %*% (diag(sqrt(pmax(e$values, 0) * r$radius_squared_95)) %*% t(e$vectors))
		data.table(
			actor = r$actor,
			year = r$year,
			category = r$category,
			x = xy[, 1] + r$x,
			y = xy[, 2] + r$y
		)
	}
	recent = p[year %in% c(2024, 2025) & category != "Other countries"]
	ell = rbindlist(
		lapply(
			seq_len(nrow(recent[label != ""])),
			function(i) ellipse(recent[label != ""][i])
		)
	)
	map = ggplot(recent, aes(x, y, color = category))
	map = map + geom_hline(yintercept = 0, alpha = 0.3, linewidth = 0.25)
	map = map + geom_vline(xintercept = 0, alpha = 0.3, linewidth = 0.25)
	map = map + geom_path(
		data = ell,
		aes(group = interaction(actor, year)),
		linewidth = 0.3,
		alpha = 0.55,
		show.legend = FALSE
	)
	map = map + geom_point(size = 1.4, alpha = 0.85)
	map = map + ggrepel::geom_text_repel(
		data = recent[label != ""],
		aes(label = label),
		size = 3.5,
		seed = 19421,
		max.overlaps = Inf,
		box.padding = 0.35,
		min.segment.length = 0,
		show.legend = FALSE
	)
	map = map + facet_wrap(~year, nrow = 1)
	map = map + coord_fixed(ratio = 0.75)
	map = map + scale_color_manual(values = colors, drop = TRUE, name = NULL)
	map = map + labs(x = "Latent dimension 1", y = "Latent dimension 2")
	map = map + style
	map = map + theme(text = element_text(size = 12), strip.text = element_text(size = 11))
	map = map + guides(color = guide_legend(nrow = 2, byrow = TRUE))
	ggsave(
		filename = file.path(graphics_root, "cooperation_space.pdf"),
		plot = map,
		width = 9,
		height = 4.5,
		device = cairo_pdf
	)
	trim_position_map(file.path(graphics_root, "cooperation_space.pdf"))
	events = data.table(
		episode = c(
			"Post-Cold War: 1989 to 1992",
			"2000 to 2002",
			"First Trump term: 2016 to 2020",
			"2024 to 2025"
		),
		before = c(1989, 2000, 2016, 2024),
		after = c(1992, 2002, 2020, 2025)
	)
	hp = list()
	hs = list()
	for (i in seq_len(nrow(events))) {
		ev = events[i]
		a = copy(p[year == ev$before])
		b = copy(p[year == ev$after])
		a[, stage := "Before"]
		b[, stage := "After"]
		z = rbind(a, b)
		z[, episode := ev$episode]
		hp[[i]] = z
		seg = merge(
			a[category != "Other countries", .(actor, x, y, category)],
			b[, .(actor, xend = x, yend = y)],
			by = "actor"
		)
		seg[, episode := ev$episode]
		hs[[i]] = seg
	}
	hp = rbindlist(hp)
	hs = rbindlist(hs)
	hp[, episode := factor(episode, levels = events$episode)]
	hs[, episode := factor(episode, levels = events$episode)]
	hp[, stage := factor(stage, levels = c("Before", "After"))]
	labels = copy(hp[category != "Other countries"])
	labels[, label := ""]
	labels[(stage == "After" & actor %in% focus) | (stage == "Before" & actor %in% c("SUN", "GER")) | (actor == "POL" & stage == "After" & episode == events$episode[1]), label := unname(country_names[actor])]
	history = ggplot(hp, aes(x, y))
	history = history + geom_point(
		data = hp[category == "Other countries"],
		aes(shape = stage),
		color = "grey80",
		size = 0.6,
		alpha = 0.45
	)
	history = history + geom_segment(
		data = hs,
		aes(xend = xend, yend = yend, color = category),
		arrow = grid::arrow(
			length = grid::unit(
				0.065,
				"inches"
			)
		),
		linewidth = 0.35,
		alpha = 0.7
	)
	history = history + geom_point(
		data = hp[category != "Other countries"],
		aes(color = category, shape = stage),
		size = 1.15
	)
	history = history + ggrepel::geom_text_repel(
		data = labels,
		aes(label = label, color = category),
		size = 3.2,
		seed = 19421,
		max.overlaps = Inf,
		min.segment.length = 0,
		box.padding = 0.3,
		show.legend = FALSE
	)
	history = history + facet_wrap(~episode, ncol = 2)
	history = history + coord_fixed(ratio = 0.65)
	history = history + scale_color_manual(values = colors, drop = TRUE, name = NULL)
	history = history + scale_shape_manual(values = c(Before = 1, After = 16), name = NULL)
	history = history + labs(x = "Latent dimension 1", y = "Latent dimension 2")
	history = history + style
	history = history + theme(text = element_text(size = 12), strip.text = element_text(size = 11))
	history = history + guides(
		color = guide_legend(nrow = 2, byrow = TRUE, order = 1),
		shape = guide_legend(order = 2)
	)
	ggsave(
		filename = file.path(graphics_root, "historical_space_shifts.pdf"),
		plot = history,
		width = 9,
		height = 6.5,
		device = cairo_pdf
	)
	trim_position_map(file.path(graphics_root, "historical_space_shifts.pdf"))
	invisible(file.path(root, "positions.csv"))
}

if (identical(commandArgs(trailingOnly = TRUE)[3L], "space-layout")) {
	draw_latent_space(paths, graphics_root)
	quit(save = "no", status = 0L)
}

if (identical(commandArgs(trailingOnly = TRUE)[3L], "space-summarize")) {
	summarize_latent_space(paths)
	summarize_spatial_groups(paths)
	quit(save = "no", status = 0L)
}

if (identical(commandArgs(trailingOnly = TRUE)[3L], "space-extract")) {
	write_latent_space(paths, commandArgs(trailingOnly = TRUE)[4L])
	quit(save = "no", status = 0L)
}

if (identical(commandArgs(trailingOnly = TRUE)[3L], "space")) {
	space_sources = draw_latent_space(paths, graphics_root)
	record_execution(
		paths,
		"space_figure",
		space_sources,
		list(models_refitted = FALSE)
	)
	quit(save = "no", status = 0L)
}

draw_trade_figure = function(paths, graphics_root) {
	trade_sources = file.path(paths$output, "trade/census/levels.csv")
	trade_levels = read_analysis_table(trade_sources)[variant == "none" & flow %in% c("imports", "exports")]
	stopifnot(
		nrow(trade_levels) == 68L,
		!anyDuplicated(trade_levels[, .(flow, year)]),
		identical(sort(unique(trade_levels$year)), 1992:2025),
		all(is.finite(trade_levels$allied_share)),
		all(is.finite(trade_levels$allied_value_usd)),
		all(trade_levels$allies == 14)
	)
	trade_levels[, relationship := factor(flow, c("imports", "exports"), c("US imports", "US exports"))]
	trade_theme = paper_plot_theme(base_size = 10) + theme(
		legend.title = element_blank(),
		plot.title = element_text(size = 11, face = "bold"),
		plot.subtitle = element_text(size = 9),
		legend.text = element_text(size = 8)
	)
	trade_colors = setNames(
		RColorBrewer::brewer.pal(3, "Set1")[1:2],
		c("US imports", "US exports")
	)
	trade_endpoints = merge(
		trade_levels[year == 2025],
		trade_levels[year == 2024, .(flow, previous_share = allied_share, previous_value = allied_value_usd)],
		by = "flow"
	)
	trade_endpoints[, `:=`(
		share_change_pp = 100 * (allied_share - previous_share),
		value_change_pct = 100 * (allied_value_usd / previous_value - 1)
	)]
	trade_endpoints[, share_label := sprintf("%.2f%%\n(%+.2f pp)", 100 * allied_share, share_change_pp)]
	trade_endpoints[, value_label := sprintf("$%.3ftn\n(%+.2f%%)", allied_value_usd / 1e12, value_change_pct)]
	trade_a = ggplot(
		trade_levels,
		aes(year, 100 * allied_share, color = relationship, linetype = relationship)
	)
	trade_a = trade_a + geom_line(linewidth = 0.6)
	trade_a = trade_a + geom_point(data = trade_levels[year %in% c(2024, 2025)], size = 1.6)
	trade_a = trade_a + geom_text(
		data = trade_endpoints,
		aes(x = 2026, label = share_label),
		hjust = 0,
		size = 2.8,
		lineheight = 1.1,
		show.legend = FALSE
	)
	trade_a = trade_a + scale_color_manual(values = trade_colors)
	trade_a = trade_a + scale_linetype_manual(
		values = c(`US imports` = "solid", `US exports` = "dashed")
	)
	trade_a = trade_a + scale_x_continuous(
		breaks = c(1992, 2000, 2010, 2020, 2025),
		limits = c(1992, 2033)
	)
	trade_a = trade_a + labs(
		title = "A. Shares of US merchandise trade",
		subtitle = "Labels: 2025 share; change from 2024",
		x = NULL,
		y = "Share of corresponding US total (%)"
	)
	trade_a = trade_a + scale_y_continuous(expand = expansion(mult = c(0.12, 0.06)))
	trade_a = trade_a + trade_theme
	trade_b = ggplot(
		trade_levels,
		aes(
			year,
			allied_value_usd / 1e+12,
			color = relationship,
			linetype = relationship
		)
	)
	trade_b = trade_b + geom_line(linewidth = 0.6)
	trade_b = trade_b + geom_point(data = trade_levels[year %in% c(2024, 2025)], size = 1.6)
	trade_b = trade_b + geom_text(
		data = trade_endpoints,
		aes(x = 2026, label = value_label),
		hjust = 0,
		size = 2.8,
		lineheight = 1.1,
		show.legend = FALSE
	)
	trade_b = trade_b + scale_color_manual(values = trade_colors)
	trade_b = trade_b + scale_linetype_manual(
		values = c(`US imports` = "solid", `US exports` = "dashed")
	)
	trade_b = trade_b + scale_x_continuous(
		breaks = c(1992, 2000, 2010, 2020, 2025),
		limits = c(1992, 2033)
	)
	trade_b = trade_b + labs(
		title = "B. US trade with the fourteen",
		subtitle = "Labels: 2025 dollars; change from 2024",
		x = NULL,
		y = "Nominal US dollars (trillions)"
	)
	trade_b = trade_b + trade_theme
	grDevices::pdf(file = NULL)
	on.exit(grDevices::dev.off(), add = TRUE)
	ggsave(
		filename = file.path(graphics_root, "trade_shares_flows.pdf"),
		plot = gridExtra::arrangeGrob(
			trade_a,
			trade_b,
			ncol = 1
		),
		width = 6.8,
		height = 6.8,
		device = cairo_pdf
	)
	invisible(trade_sources)
}

if (identical(commandArgs(trailingOnly = TRUE)[3L], "trade")) {
	trade_sources = draw_trade_figure(paths, graphics_root)
	record_execution(
		paths,
		"trade_figure",
		trade_sources,
		list(models_refitted = FALSE, network_access = FALSE)
	)
	quit(save = "no", status = 0L)
}

draw_agreement_history = function(paths, graphics_root) {
	input = file.path(paths$output, "observed/allies/group_history.csv")
	history = read.csv(input, stringsAsFactors = FALSE)
	history = history[history$specification == "ordinal" & history$group %in% c("us_allies", "within_allies"), ]
	stopifnot(nrow(history) == 142L)
	history$relationship = factor(
		history$group,
		levels = c("us_allies", "within_allies"),
		labels = c("US–allies", "Among the fourteen allies")
	)
	history$horizon_label = factor(
		history$horizon,
		levels = 1:2,
		labels = c("One-year change", "Two-year change")
	)
	markers = data.frame(
		year = c(1991, 2001, 2017, 2025),
		label = c("Cold War's end", "2000–02", "First Trump term", "2025"),
		relationship = factor("US–allies", levels = levels(history$relationship))
	)
	p = ggplot(history, aes(target, 100 * change))
	p = p + geom_vline(
		data = markers,
		aes(xintercept = year),
		alpha = 0.3,
		linewidth = 0.3,
		linetype = "dotted"
	)
	p = p + geom_linerange(
		aes(ymin = 100 * lower_agenda_95, ymax = 100 * upper_agenda_95),
		alpha = 0.45,
		linewidth = 0.35
	)
	p = p + geom_line(linewidth = 0.4)
	p = p + geom_point(size = 1.8)
	p = p + geom_text(
		data = markers,
		aes(x = year, y = Inf, label = label),
		inherit.aes = FALSE,
		angle = 90,
		hjust = 1.1,
		vjust = -0.35,
		size = 2.5
	)
	p = p + facet_grid(relationship ~ horizon_label)
	p = p + labs(x = "Later year", y = "Agreement change (percentage points)")
	p = p + paper_plot_theme(base_size = 11)
	ggsave(
		filename = file.path(graphics_root, "allied_history.pdf"),
		plot = p,
		width = 8.5,
		height = 5.8,
		device = cairo_pdf
	)
	invisible(input)
}

if (identical(commandArgs(trailingOnly = TRUE)[3L], "agreement-history")) {
	agreement_source = draw_agreement_history(paths, graphics_root)
	record_execution(
		paths,
		"agreement_history_figure",
		agreement_source,
		list(
			models_refitted = FALSE,
			summary_tables_recomputed = FALSE,
			network_access = FALSE
		)
	)
	quit(save = "no", status = 0L)
}

input = draw_agreement_history(paths, graphics_root)
trade_sources = draw_trade_figure(paths, graphics_root)
space_sources = draw_latent_space(paths, graphics_root)
record_execution(paths, "figures", c(input, trade_sources, space_sources))
