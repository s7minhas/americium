# compare network positions with dynamic ideal points

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
		library(posterior)
	}
)
mode = commandArgs(TRUE)[1]
if (is.na(mode)) mode = "compare"
stopifnot(mode %in% c("download", "compare"))
argument = sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
root = dirname(dirname(normalizePath(argument)))
data_root = file.path(root, "data", "ideal_points")
output_root = file.path(root, "results", "models", "un_variants", "ideal_points")
doi = "doi:10.7910/DVN/LEJUQZ"
api = "https://dataverse.harvard.edu/api"
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
years = 1989:2025
set.seed(6886)
or_na = function(x) if (is.null(x)) NA else x
timestamp = function() format(Sys.time(), tz = "UTC", usetz = TRUE)
file_record = function(path) {
	list(
		path = sub(paste0("^", root, "/"), "", path),
		sha256 = digest(file = path, algo = "sha256")
	)
}

wanted = data.table(
	file_id = c(14098429L, 13642024L, 11824621L, 11824371L, 11824371L),
	filename = c(
		"IdealpointestimatesFP_2026FP.csv",
		"Codebook_IdealPointEstimates1946-2025.txt",
		"Codebook United Nations General Assembly Ideal Points.docx",
		"IdealpointsJuly2025.dta",
		"IdealpointsJuly2025.tab"
	),
	access = c("default", "default", "default", "original", "default")
)
# one row per file of the dataset version, from the api listing
version_listing = function(version) {
	rbindlist(
		lapply(
			version$files,
			function(f) {
				data.table(
					file_id = f$dataFile$id,
					filename = f$dataFile$filename,
					description = or_na(f$description),
					original_filename = or_na(f$dataFile$originalFileName),
					content_type = or_na(f$dataFile$contentType),
					bytes = f$dataFile$filesize,
					checksum_type = f$dataFile$checksum$type,
					checksum = f$dataFile$checksum$value,
					file_created = or_na(f$dataFile$creationDate),
					file_published = or_na(f$dataFile$publicationDate)
				)
			}
		)
	)
}

download_one = function(row, listing) {
	url = paste0(
		api,
		"/access/datafile/",
		row$file_id,
		if (row$access == "original") "?format=original" else ""
	)
	path = file.path(data_root, row$filename)
	download.file(url, path, mode = "wb", quiet = TRUE)
	listed = listing[file_id == row$file_id]
	local_md5 = digest(file = path, algo = "md5")
	verifiable = row$access == "original" || is.na(listed$original_filename)
	stopifnot(!verifiable || local_md5 == listed$checksum)
	list(
		filename = row$filename,
		dataverse_file_id = row$file_id,
		url = url,
		description = listed$description,
		dataverse_md5 = listed$checksum,
		local_md5 = local_md5,
		md5_verified = if (verifiable) TRUE else NA,
		sha256 = digest(file = path, algo = "sha256"),
		bytes = file.size(path),
		retrieved_utc = timestamp(),
		file_created = listed$file_created,
		file_published = listed$file_published,
		form = if (verifiable) "deposited file" else "dataverse tabular ingest of the deposited file"
	)
}

check_tabular_ingest = function() {
	tab = read_analysis_table(file.path(data_root, "IdealpointsJuly2025.tab"))
	dta = haven::read_dta(file.path(data_root, "IdealpointsJuly2025.dta"))
	stopifnot(
		nrow(tab) == nrow(dta),
		isTRUE(all.equal(as.integer(tab$year), as.integer(dta$year))),
		isTRUE(all.equal(tab$ccode, as.integer(dta$ccode))),
		isTRUE(
			all.equal(tab$idealpointfp, as.numeric(dta$idealpointfp), tolerance = 1e-6)
		),
		isTRUE(
			all.equal(tab$idealpointall, as.numeric(dta$idealpointall), tolerance = 1e-6)
		)
	)
	tab
}

# years, sessions, rows and countries covered by the two estimate files
coverage_record = function(tab) {
	current = read_analysis_table(file.path(data_root, "IdealpointestimatesFP_2026FP.csv"))
	list(
		current_file = "IdealpointestimatesFP_2026FP.csv",
		current_rows = nrow(current),
		current_first_year = min(current$year),
		current_last_year = max(current$year),
		current_sessions = paste(range(as.integer(current$session)), collapse = " to "),
		current_countries_last_year = current[year == max(year), uniqueN(iso3c)],
		current_unit = "country and calendar year; one session per year in the file",
		july_2025_file = "IdealpointsJuly2025.tab",
		july_2025_rows = nrow(tab),
		july_2025_last_year_final_passage = max(tab$year[!is.na(tab$idealpointfp)]),
		july_2025_last_year_all_votes = max(tab$year[!is.na(tab$idealpointall)])
	)
}

download_files = function() {
	dir.create(data_root, recursive = TRUE, showWarnings = FALSE)
	dataset_path = file.path(data_root, "dataverse_dataset.json")
	versions_path = file.path(data_root, "dataverse_versions.json")
	dataset_url = paste0(api, "/datasets/:persistentId/?persistentId=", doi)
	versions_url = paste0(api, "/datasets/:persistentId/versions?persistentId=", doi)
	download.file(dataset_url, dataset_path, mode = "wb", quiet = TRUE)
	download.file(versions_url, versions_path, mode = "wb", quiet = TRUE)
	dataset = fromJSON(dataset_path, simplifyVector = FALSE)$data
	version = dataset$latestVersion
	fields = version$metadataBlocks$citation$fields
	field = function(name) fields[[which(vapply(fields, function(f) f$typeName == name, logical(1)))]]$value
	listing = version_listing(version)
	stopifnot(all(wanted$file_id %in% listing$file_id))
	downloaded = lapply(seq_len(nrow(wanted)), function(i) download_one(wanted[i], listing))
	coverage = coverage_record(check_tabular_ingest())
	record = list(
		created_utc = timestamp(),
		dataset = list(
			doi = sub("^doi:", "", doi),
			persistent_url = dataset$persistentUrl,
			title = field("title"),
			authors = vapply(field("author"), function(a) a$authorName$value, character(1)),
			dataset_id = dataset$id,
			dataset_first_published = dataset$publicationDate,
			version = paste(version$versionNumber, version$versionMinorNumber, sep = "."),
			version_released = version$releaseTime,
			license = version$license$name,
			distribution_date = field("distributionDate"),
			depositor_citation_request = field("publication")[[1]]$publicationCitation$value,
			dataset_api_url = dataset_url,
			versions_api_url = versions_url,
			api_responses = c("dataverse_dataset.json", "dataverse_versions.json")
		),
		version_file_listing = listing,
		downloaded = downloaded,
		coverage = coverage
	)
	write_json(
		record,
		file.path(data_root, "retrieval.json"),
		pretty = TRUE,
		auto_unbox = TRUE,
		na = "null"
	)
}

model_actor = function(iso3c, ccode, year) {
	actor = iso3c
	actor[ccode == 365L & year <= 1991L] = "SUN"
	actor[ccode == 365L & year >= 1992L] = "RUS"
	actor[ccode == 345L & year <= 2002L] = "YUG"
	actor[ccode == 345L & year %in% 2003:2005] = "SCG"
	actor[ccode == 345L & year >= 2006L] = "SRB"
	actor[ccode == 260L & year == 1989L] = "GER"
	actor[ccode == 678L & year == 1989L] = "YEM"
	actor[ccode == 680L & year == 1989L] = "YMD"
	actor
}

annual_change = function(d, column, by = "ccode") {
	previous = paste0(column, "_previous")
	setorderv(d, c(by, "year"))
	d[, previous_year := shift(year), by = by]
	d[, (previous) := shift(get(column)), by = by]
	d[, (paste0(column, "_change")) := fifelse(
		!is.na(previous_year) & previous_year == year - 1L,
		get(column) - get(previous),
		NA_real_
	)]
	d[, previous_year := NULL]
	d
}

rank_descending = function(x) frank(-x, ties.method = "min", na.last = "keep")

correlate = function(d, model_measure, ideal_measure) {
	a = d[[model_measure]]
	b = d[[ideal_measure]]
	ok = is.finite(a) & is.finite(b)
	data.table(
		model_measure = model_measure,
		ideal_measure = ideal_measure,
		n = sum(ok),
		pearson = if (sum(ok) > 2) cor(a[ok], b[ok]) else NA_real_,
		spearman = if (sum(ok) > 2) cor(a[ok], b[ok], method = "spearman") else NA_real_
	)
}

read_ideal_points = function(path) {
	current = read_analysis_table(path)
	current = current[, .(
		ccode,
		iso3c,
		year = as.integer(year),
		ideal_point = IdealPointFP,
		ideal_point_lower = `Q5%FP`,
		ideal_point_upper = `Q95%FP`,
		ideal_point_votes = NVotesFP
	)]
	current[, actor := model_actor(iso3c, ccode, year)]
	full_series = annual_change(copy(current), "ideal_point")
	full_series[, `:=`(
		ideal_point_change_rank_abs_country_full_series = rank_descending(abs(ideal_point_change)),
		ideal_point_change_rank_signed_country_full_series = rank_descending(ideal_point_change),
		ideal_point_change_n_country_full_series = sum(!is.na(ideal_point_change))
	), by = ccode]
	full_series[, `:=`(
		ideal_point_change_rank_abs_all_full_series = rank_descending(abs(ideal_point_change)),
		ideal_point_change_rank_signed_all_full_series = rank_descending(ideal_point_change),
		ideal_point_change_n_all_full_series = sum(!is.na(ideal_point_change))
	)]
	window = annual_change(current[year %in% years], "ideal_point")
	stopifnot(!anyDuplicated(window, by = c("actor", "year")))
	list(full_series = full_series, window = window)
}

# its rows without a country label carry no estimates and are dropped
read_july_2025 = function(path) {
	july = read_analysis_table(path)
	july = july[, .(
		ccode,
		iso3c,
		year = as.integer(year),
		ideal_point_all = idealpointall,
		ideal_point_july_2025 = idealpointfp
	)]
	unlabeled = is.na(july$iso3c) | july$iso3c == ""
	stopifnot(
		all(is.na(july$ideal_point_all[unlabeled])),
		all(is.na(july$ideal_point_july_2025[unlabeled]))
	)
	july = july[!unlabeled]
	july[, actor := model_actor(iso3c, ccode, year)]
	july = july[year %in% years]
	july = annual_change(july, "ideal_point_all")
	july = annual_change(july, "ideal_point_july_2025")
	july = july[, .(
		actor,
		year,
		ideal_point_all,
		ideal_point_all_change,
		ideal_point_july_2025,
		ideal_point_july_2025_change
	)]
	stopifnot(!anyDuplicated(july, by = c("actor", "year")))
	july
}

read_diagnostics = function(path) {
	diagnostics = read_analysis_table(path)
	rho = diagnostics[specification == "full" & parameter == "rho_uv", mean]
	stopifnot(length(rho) == 1L)
	log_odds = diagnostics[specification == "full" & family == "actor_log_odds", .(
		actor = sub("_[0-9]{4}$", "", parameter),
		year = as.integer(sub("^.*_", "", parameter)),
		log_odds = mean
	)]
	list(rho = rho, log_odds = log_odds)
}

procrustes = function(raw, reference, keep) {
	v = svd(crossprod(raw[keep, ], reference[keep, ]))
	v$u %*% t(v$v)
}

align_draws = function(path, rho) {
	archive = readRDS(path)
	observed = archive$observed
	reference = archive$reference
	actors = rownames(observed)
	stopifnot(identical(colnames(observed), as.character(years)))
	n = length(actors)
	periods = length(years)
	keep = observed[, "2024"]
	chains = names(archive$raw_positions)
	draws = dim(archive$raw_positions[[1]])[4]
	total = draws * length(chains)
	level_x = level_y = matrix(0, n, periods, dimnames = list(actors, years))
	axis_one_change = axis_two_change = displacement = displacement_rho_adjusted = array(
		NA_real_,
		c(n, periods - 1L, total),
		dimnames = list(actors, years[-1], NULL)
	)
	index = 0L
	for (chain in chains) {
		u = archive$raw_positions[[chain]]
		stopifnot(identical(dimnames(u)[[1]], actors), dim(u)[4] == draws)
		for (d in seq_len(draws)) {
			index = index + 1L
			transformation = procrustes(u[, , "2024", d], reference, keep)
			z = matrix(aperm(u[, , , d], c(1, 3, 2)), n * periods, 2) %*% transformation
			x = matrix(z[, 1], n, periods)
			y = matrix(z[, 2], n, periods)
			x[!observed] = NA_real_
			y[!observed] = NA_real_
			level_x = level_x + fifelse(is.na(x), 0, x)
			level_y = level_y + fifelse(is.na(y), 0, y)
			dx = x[, -1] - x[, -periods]
			dy = y[, -1] - y[, -periods]
			axis_one_change[, , index] = dx
			axis_two_change[, , index] = dy
			displacement[, , index] = sqrt(dx^2 + dy^2)
			displacement_rho_adjusted[, , index] = sqrt((x[, -1] - rho * x[, -periods])^2 + (y[, -1] - rho * y[, -periods])^2)
		}
		rm(u)
		gc(FALSE)
	}
	rm(archive)
	gc(FALSE)
	level_x = level_x / total
	level_y = level_y / total
	level_x[!observed] = NA_real_
	level_y[!observed] = NA_real_
	list(
		actors = actors,
		chains = chains,
		draws = draws,
		level_x = level_x,
		level_y = level_y,
		axis_one_change = axis_one_change,
		axis_two_change = axis_two_change,
		displacement = displacement,
		displacement_rho_adjusted = displacement_rho_adjusted
	)
}

summarize_cells = function(a, name, chains, draws) {
	cell = function(f) apply(a, c(1, 2), f)
	quantile_cell = function(probability) {
		cell(
			function(v) if (anyNA(v)) NA_real_ else quantile(v, probability, names = FALSE)
		)
	}
	chain_matrix = function(v) matrix(v, draws, length(chains))
	out = data.table(
		actor = rep(dimnames(a)[[1]], dim(a)[2]),
		year = rep(as.integer(dimnames(a)[[2]]), each = dim(a)[1]),
		mean = as.vector(cell(mean)),
		lower = as.vector(quantile_cell(0.025)),
		upper = as.vector(quantile_cell(0.975)),
		rhat = as.vector(cell(function(v) if (anyNA(v)) NA_real_ else rhat(chain_matrix(v)))),
		ess_bulk = as.vector(
			cell(function(v) if (anyNA(v)) NA_real_ else ess_bulk(chain_matrix(v)))
		)
	)
	setnames(
		out,
		c("mean", "lower", "upper", "rhat", "ess_bulk"),
		paste0(name, c("", "_lower", "_upper", "_rhat", "_ess_bulk"))
	)
	out
}

model_changes = function(aligned, log_odds) {
	summarize = function(a, name) summarize_cells(a, name, aligned$chains, aligned$draws)
	model = merge(
		summarize(aligned$axis_one_change, "axis_one_change"),
		summarize(aligned$axis_two_change, "axis_two_change"),
		by = c("actor", "year")
	)
	model = merge(
		model,
		summarize(aligned$displacement, "displacement"),
		by = c("actor", "year")
	)
	adjusted = data.table(
		actor = rep(aligned$actors, length(years) - 1L),
		year = rep(years[-1], each = length(aligned$actors)),
		displacement_rho_adjusted = as.vector(apply(aligned$displacement_rho_adjusted, c(1, 2), mean))
	)
	model = merge(model, adjusted, by = c("actor", "year"))
	model = model[!is.na(axis_one_change)]
	model = merge(model, log_odds, by = c("actor", "year"), all.x = TRUE)
	stopifnot(!anyNA(model$log_odds))
	model
}

model_levels = function(aligned) {
	levels = data.table(
		actor = rep(aligned$actors, length(years)),
		year = rep(years, each = length(aligned$actors)),
		axis_one = as.vector(aligned$level_x),
		axis_two = as.vector(aligned$level_y)
	)
	levels[!is.na(axis_one)]
}

position_checks = function(aligned, path) {
	positions = read_analysis_table(path)
	difference = function(year_value) {
		p = positions[year == year_value]
		column = as.character(year_value)
		max(
			abs(
				c(
					p$x - aligned$level_x[p$actor, column],
					p$y - aligned$level_y[p$actor, column]
				)
			)
		)
	}
	checks = list(
		positions_2024_max_difference = difference(2024),
		positions_2025_max_difference = difference(2025)
	)
	stopifnot(checks$positions_2024_max_difference < 1e-8)
	checks
}

axis_one_orientation = function(levels) {
	year_facts = function(year_value) {
		d = levels[year == year_value]
		facts = list(
			members = nrow(d),
			members_negative = sum(d$axis_one < 0),
			usa_axis_one = d[actor == "USA", axis_one],
			usa_rank_from_negative_end = frank(d$axis_one, ties.method = "min")[d$actor == "USA"]
		)
		setNames(facts, paste0(names(facts), "_", year_value))
	}
	c(
		list(
			country_years = nrow(levels),
			country_years_negative = sum(levels$axis_one < 0)
		),
		year_facts(2024L),
		year_facts(2025L)
	)
}

country_year_table = function(levels, model, ideal, july) {
	table = merge(levels, model, by = c("actor", "year"), all = TRUE)
	table = merge(
		table,
		ideal[, .(
			actor,
			year,
			ideal_point,
			ideal_point_lower,
			ideal_point_upper,
			ideal_point_votes,
			ideal_point_previous,
			ideal_point_change
		)],
		by = c("actor", "year"),
		all = TRUE
	)
	table = merge(table, july, by = c("actor", "year"), all.x = TRUE)
	table[, `:=`(
		axis_one_change_abs = abs(axis_one_change),
		axis_two_change_abs = abs(axis_two_change),
		ideal_point_change_abs = abs(ideal_point_change),
		ideal_point_all_change_abs = abs(ideal_point_all_change)
	)]
	table[, `:=`(
		ideal_point_change_z = (ideal_point_change - mean(ideal_point_change, na.rm = TRUE)) / sd(ideal_point_change, na.rm = TRUE),
		ideal_point_change_rank_abs_country = rank_descending(ideal_point_change_abs),
		ideal_point_change_rank_signed_country = rank_descending(ideal_point_change),
		ideal_point_change_n_country = sum(!is.na(ideal_point_change)),
		axis_one_change_rank_abs_country = rank_descending(axis_one_change_abs),
		axis_one_change_rank_signed_country = rank_descending(axis_one_change),
		axis_two_change_rank_abs_country = rank_descending(axis_two_change_abs),
		axis_two_change_rank_signed_country = rank_descending(axis_two_change),
		displacement_rank_country = rank_descending(displacement),
		displacement_rho_adjusted_rank_country = rank_descending(displacement_rho_adjusted),
		log_odds_rank_country = rank_descending(log_odds),
		model_n_country = sum(!is.na(displacement))
	), by = actor]
	table[, `:=`(
		ideal_point_change_rank_abs_all = rank_descending(ideal_point_change_abs),
		ideal_point_change_rank_signed_all = rank_descending(ideal_point_change),
		ideal_point_change_n_all = sum(!is.na(ideal_point_change)),
		axis_one_change_rank_abs_all = rank_descending(axis_one_change_abs),
		axis_one_change_rank_signed_all = rank_descending(axis_one_change),
		axis_two_change_rank_abs_all = rank_descending(axis_two_change_abs),
		axis_two_change_rank_signed_all = rank_descending(axis_two_change),
		displacement_rank_all = rank_descending(displacement),
		displacement_rho_adjusted_rank_all = rank_descending(displacement_rho_adjusted),
		log_odds_rank_all = rank_descending(log_odds),
		model_n_all = sum(!is.na(displacement))
	)]
	table[, `:=`(
		ally = actor %in% allies,
		both_measures = !is.na(ideal_point_change) & !is.na(displacement)
	)]
	setorder(table, actor, year)
	table
}

correlation_table = function(table) {
	pairs = list(
		c("axis_one", "ideal_point"),
		c("axis_two", "ideal_point"),
		c("axis_one_change", "ideal_point_change"),
		c("axis_two_change", "ideal_point_change"),
		c("axis_one_change_abs", "ideal_point_change_abs"),
		c("axis_two_change_abs", "ideal_point_change_abs"),
		c("displacement", "ideal_point_change_abs"),
		c("displacement_rho_adjusted", "ideal_point_change_abs"),
		c("log_odds", "ideal_point_change_abs"),
		c("axis_two_change", "ideal_point_all_change"),
		c("displacement", "ideal_point_all_change_abs"),
		c("axis_two_change", "ideal_point_july_2025_change")
	)
	rows = list()
	for (without_usa in c(FALSE, TRUE)) {
		d = if (without_usa) table[actor != "USA"] else table
		for (pair in pairs) {
			rows[[length(rows) + 1L]] = cbind(
				data.table(scope = "pooled", year = NA_integer_, without_usa = without_usa),
				correlate(d, pair[1], pair[2])
			)
			for (year_value in years) {
				rows[[length(rows) + 1L]] = cbind(
					data.table(scope = "year", year = year_value, without_usa = without_usa),
					correlate(d[year == year_value], pair[1], pair[2])
				)
			}
		}
	}
	comparison = rbindlist(rows)[!(scope == "year" & n == 0L)]
	year_summary = comparison[scope == "year", .(
		scope = c("year_minimum", "year_median", "year_maximum"),
		year = NA_integer_,
		n = .N,
		pearson = c(
			min(pearson, na.rm = TRUE),
			median(pearson, na.rm = TRUE),
			max(pearson, na.rm = TRUE)
		),
		spearman = c(
			min(spearman, na.rm = TRUE),
			median(spearman, na.rm = TRUE),
			max(spearman, na.rm = TRUE)
		)
	), by = .(without_usa, model_measure, ideal_measure)]
	rbind(comparison, year_summary[, names(comparison), with = FALSE])
}

american_table = function(full_series, table) {
	full = full_series[actor == "USA", .(
		year,
		ideal_point_full_series = ideal_point,
		ideal_point_change_full_series = ideal_point_change,
		ideal_point_change_rank_abs_country_full_series,
		ideal_point_change_rank_signed_country_full_series,
		ideal_point_change_n_country_full_series,
		ideal_point_change_rank_abs_all_full_series,
		ideal_point_change_rank_signed_all_full_series,
		ideal_point_change_n_all_full_series
	)]
	merge(full, table[actor == "USA"], by = "year", all.x = TRUE)
}

# years of each actor present in one series only, one string per actor
years_by_actor = function(d) {
	out = d[, .(years = paste(year, collapse = " ")), by = actor][order(actor)]
	as.list(setNames(out$years, out$actor))
}

provenance_record = function(paths, retrieval, meta, rho, table, orientation, checks) {
	alignment = paste(
		"One orthogonal transformation (a rotation or a reflection) per draw, fitted by Procrustes analysis of the 2024 configuration",
		"on the countries observed in 2024 to the fixed reference stored in draws.rds (19_figures.R) and applied to all 37 years of that draw.",
		"Within-draw displacements between years are unchanged by this transformation.",
		"Axis one is the first principal axis of the uncentered reference configuration, the common direction in which nearly every member lies,",
		"so that nearly every coordinate on it is negative (axis_one_orientation); axis two is the second principal axis,",
		"with the United States positive. Both keep the orientation of positions.csv."
	)
	axis_changes = paste(
		"Posterior mean, 2.5 and 97.5 percent quantiles over all draws of the change in each coordinate from the previous year;",
		"country-years observed in both years.",
		"The level correlations in comparison.csv show which axis corresponds to the ideal-point dimension."
	)
	displacement = paste(
		"Posterior mean and quantiles of the Euclidean distance between the positions in consecutive years,",
		"as in the log-odds expression of 10_models.R;",
		"the rho-adjusted version measures the distance from rho times the previous position at the posterior mean of rho."
	)
	code_mapping = list(
		SUN = "Correlates of War 365 (RUS) in 1989 to 1991",
		RUS = "Correlates of War 365 (RUS) from 1992",
		YUG = "Correlates of War 345 (YUG, SCG) through 2002",
		SCG = "Correlates of War 345 (SCG; YUG in the July 2025 file) in 2003 to 2005",
		SRB = "Correlates of War 345 (SRB; YUG in the July 2025 file) from 2006",
		GER = "Correlates of War 260 (DEU) in 1989",
		YEM = "Correlates of War 678 (YAR) in 1989 and 679 (YEM) from 1990",
		YMD = "Correlates of War 680 (YMD; YAR in the July 2025 file) in 1989",
		other = "the iso3c label of the ideal-point file",
		annual_changes = paste(
			"computed within each Correlates of War code and, in the model, within each actor code,",
			"so a country-year compares changes only when both series carry the previous year under the same code"
		)
	)
	list(
		script = "code/15_ideal_points.R",
		mode = "compare",
		created_utc = timestamp(),
		inputs = lapply(paths, file_record),
		ideal_point_vintage = list(
			dataset_version = retrieval$dataset$version,
			version_released = retrieval$dataset$version_released,
			current_file = retrieval$coverage$current_file,
			current_file_id = retrieval$downloaded[[1]]$dataverse_file_id,
			first_year = retrieval$coverage$current_first_year,
			last_year = retrieval$coverage$current_last_year
		),
		alignment = alignment,
		axis_changes = axis_changes,
		displacement = displacement,
		axis_one_orientation = orientation,
		rho_posterior_mean = rho,
		chains = meta$chains,
		draws_per_chain = meta$draws,
		window = range(years),
		code_mapping = code_mapping,
		unmatched_country_years = list(
			model_only = years_by_actor(table[!is.na(axis_one) & is.na(ideal_point)]),
			ideal_point_only = years_by_actor(table[is.na(axis_one) & !is.na(ideal_point)])
		),
		counts = list(
			country_years_with_both_levels = sum(!is.na(table$axis_one) & !is.na(table$ideal_point)),
			country_years_with_both_changes = sum(table$both_measures),
			country_years_with_model_change = sum(!is.na(table$displacement)),
			country_years_with_ideal_point_change = sum(!is.na(table$ideal_point_change)),
			countries_with_both_changes = uniqueN(table[both_measures == TRUE, actor])
		),
		checks = checks,
		R = R.version.string
	)
}

compare = function() {
	dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
	paths = list(
		draws = file.path(root, "results", "models", "un", "space", "draws.rds"),
		diagnostics = file.path(root, "results", "models", "un", "diagnostics.csv"),
		positions = file.path(root, "results", "models", "un", "space", "positions.csv"),
		ideal_points_current = file.path(data_root, "IdealpointestimatesFP_2026FP.csv"),
		ideal_points_july_2025 = file.path(data_root, "IdealpointsJuly2025.tab"),
		retrieval = file.path(data_root, "retrieval.json")
	)
	retrieval = fromJSON(paths$retrieval, simplifyVector = FALSE)
	ideal = read_ideal_points(paths$ideal_points_current)
	july = read_july_2025(paths$ideal_points_july_2025)
	diagnostics = read_diagnostics(paths$diagnostics)
	aligned = align_draws(paths$draws, diagnostics$rho)
	levels = model_levels(aligned)
	checks = position_checks(aligned, paths$positions)
	model = model_changes(aligned, diagnostics$log_odds)
	meta = list(chains = aligned$chains, draws = aligned$draws)
	rm(aligned)
	gc(FALSE)
	table = country_year_table(levels, model, ideal$window, july)
	comparison = correlation_table(table)
	american = american_table(ideal$full_series, table)
	allied_2025 = table[year == 2025L & actor %in% c("USA", allies)][order(-ally, -displacement)]
	checks$usa_2025_log_odds_posterior_mean = diagnostics$log_odds[actor == "USA" & year == 2025L, log_odds]
	provenance = provenance_record(
		paths,
		retrieval,
		meta,
		diagnostics$rho,
		table,
		axis_one_orientation(levels),
		checks
	)
	fwrite(
		table[year > min(years)],
		file.path(output_root, "country_year_changes.csv")
	)
	fwrite(comparison, file.path(output_root, "comparison.csv"))
	fwrite(american, file.path(output_root, "american_changes.csv"))
	fwrite(allied_2025, file.path(output_root, "allied_changes_2025.csv"))
	write_json(
		provenance,
		file.path(output_root, "provenance.json"),
		pretty = TRUE,
		auto_unbox = TRUE,
		digits = NA
	)
}

if (mode == "download") download_files() else compare()
