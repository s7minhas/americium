# check the allied comparison group

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
		library(readxl)
		library(jsonlite)
		library(digest)
	}
)
argument = sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
root = dirname(dirname(normalizePath(argument)))
merge_groups = "merge" %in% commandArgs(TRUE)
output_root = file.path(root, "results", "observed", "ally_selection")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
expected_allies = c(
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
window_start = 2001
window_end = 2025
gdp_threshold = 400
share_threshold = 0.5
table_gdp_floor = 200
inputs = c(
	votes = "data/un_votes_long.rda",
	treaties = "data/ally_selection/treaties.csv",
	nato_page = "data/ally_selection/nato/nato_countries.html",
	oas_page = "data/ally_selection/oas_rio/b-29_en.html",
	class_current = "data/ally_selection/world_bank_class/CLASS.xlsx",
	class_history = "data/ally_selection/world_bank_class/OGHIST.xlsx",
	wb_gdp = "data/ally_selection/world_bank_api/NY.GDP.MKTP.CD_2024.json",
	wb_countries = "data/ally_selection/world_bank_api/countries.json",
	imf_gdp = "data/ally_selection/imf_weo/ngdpd_2024_transcribed.csv",
	dac_page = "data/ally_selection/oecd_dac/dac_body_869.html",
	scale_text = "data/ally_selection/un_scale/scale_of_assessments_1946_2027.txt",
	assessment_text = "data/ally_selection/un_scale/A_79_318_Add.1.txt",
	retrieval = "data/ally_selection/retrieval.json",
	member_groups = "data/member_groups.csv"
)
path_of = function(name) file.path(root, inputs[[name]])

simplify_name = function(x) {
	x = gsub("\\(formerly[^)]*\\)", "", x)
	x = gsub("[^ -~]", "", x, perl = TRUE)
	gsub("[^a-z]", "", tolower(x))
}

# text lines of a stored html page without tags
html_lines = function(path) {
	text = paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
	text = gsub("(?is)<(script|style)[^>]*>.*?</\\1>", " ", text, perl = TRUE)
	text = gsub("<[^>]+>", "\n", text)
	entities = c("&nbsp;" = " ", "&amp;" = "&", "&#39;" = "'", "&quot;" = '"')
	for (entity in names(entities)) {
		text = gsub(entity, entities[[entity]], text, fixed = TRUE)
	}
	lines = trimws(gsub("\\s+", " ", strsplit(text, "\n", fixed = TRUE)[[1]]))
	lines[nzchar(lines)]
}

vote_env = new.env(parent = emptyenv())
load_vote_archive(path_of("votes"), envir = vote_env)
votes_2025 = as.data.table(vote_env$un_votes_long)[year == 2025]
members = votes_2025[, .(member_name = voting_country_name[.N]), by = .(iso3 = voting_country_code)][order(iso3)]
stopifnot(nrow(members) == 193L, all(expected_allies %in% members$iso3))
member_key = members[, .(iso3, key = simplify_name(member_name))]

treaties = read_analysis_table(path_of("treaties"), na.strings = "")
stopifnot(
	ncol(treaties) == 10L,
	all(treaties$iso3 %in% members$iso3),
	all(nzchar(treaties$source_url))
)
covers_window = function(start_year, from, to) !is.na(from) & from < start_year & (is.na(to) | to > window_end)
treaties[, covers := covers_window(window_start, in_force_year, end_year)]
treaty_summary = treaties[, .(
	treaty = paste(unique(treaty), collapse = "; "),
	treaty_group = paste(unique(treaty_group), collapse = "; "),
	treaty_in_force_year = suppressWarnings(min(in_force_year, na.rm = TRUE)),
	treaty_end_year = if (any(is.na(end_year))) NA_integer_ else max(end_year),
	treaty_in_force_2001_2025 = any(covers),
	treaty_note = paste(status_note, collapse = " | ")
), by = iso3]
treaty_summary[!is.finite(treaty_in_force_year), treaty_in_force_year := NA_integer_]
# the nato accession years in treaties.csv must match the stored nato page
nato_html = paste(
	readLines(path_of("nato_page"), warn = FALSE, encoding = "UTF-8"),
	collapse = "\n"
)
nato_pieces = strsplit(nato_html, 'class="cie-country__name"', fixed = TRUE)[[1]][-1]
nato_names = sub('(?s)^.*?data-title="([^"]+)".*$', "\\1", nato_pieces, perl = TRUE)
nato_years = as.integer(
	sub(
		'(?s)^.*?cie-country__yearEntry">\\s*([0-9]{4}).*$',
		"\\1",
		nato_pieces,
		perl = TRUE
	)
)
nato_lookup = c(
	albania = "ALB",
	belgium = "BEL",
	bulgaria = "BGR",
	canada = "CAN",
	croatia = "HRV",
	czechia = "CZE",
	denmark = "DNK",
	estonia = "EST",
	finland = "FIN",
	france = "FRA",
	germany = "DEU",
	greece = "GRC",
	hungary = "HUN",
	iceland = "ISL",
	italy = "ITA",
	latvia = "LVA",
	lithuania = "LTU",
	luxembourg = "LUX",
	montenegro = "MNE",
	northmacedonia = "MKD",
	norway = "NOR",
	poland = "POL",
	portugal = "PRT",
	romania = "ROU",
	slovakia = "SVK",
	slovenia = "SVN",
	spain = "ESP",
	sweden = "SWE",
	thenetherlands = "NLD",
	trkiye = "TUR",
	unitedkingdom = "GBR",
	unitedstates = "USA"
)
nato_page = data.table(
	iso3 = unname(nato_lookup[simplify_name(nato_names)]),
	nato_year = nato_years
)
stopifnot(
	nrow(nato_page) == 32L,
	!anyNA(nato_page$iso3),
	!anyNA(nato_page$nato_year)
)
nato_check = merge(
	treaties[treaty_group == "NATO", .(iso3, in_force_year)],
	nato_page,
	by = "iso3",
	all = TRUE
)
stopifnot(
	nrow(nato_check) == 32L,
	nato_check[iso3 != "USA", all(in_force_year == nato_year)],
	nato_check[iso3 == "USA", is.na(in_force_year)]
)

oas_text = paste(
	readLines(path_of("oas_page"), warn = FALSE, encoding = "UTF-8"),
	collapse = " "
)
oas_sentinels = c(
	"09/06/02",
	"09/06/04",
	"October 17th, 2012",
	"October 17th, 2014",
	"February 19, 2014",
	"May 14, 2013",
	"August 6, 2019",
	"01/22/90",
	"12/16/91",
	"September 27, 2019",
	"3/10/20",
	"September 20, 2012",
	"12/03/48"
)
stopifnot(
	all(
		vapply(
			oas_sentinels,
			function(s) grepl(s, oas_text, fixed = TRUE),
			logical(1)
		)
	)
)

history = read_excel(
	path_of("class_history"),
	sheet = "Country Analytical History",
	col_names = FALSE,
	.name_repair = "minimal"
)
header_row = which(history[[2]] == "Bank's fiscal year:")
fiscal_years = unlist(history[header_row, ])
code_rows = which(grepl("^[A-Z]{3}$", history[[1]]))
group_labels = c(
	H = "High income",
	UM = "Upper middle income",
	LM = "Lower middle income",
	L = "Low income"
)
label_group = function(code) {
	ifelse(
		is.na(code) | !code %in% names(group_labels),
		"Not classified",
		unname(group_labels[code])
	)
}

income = data.table(
	iso3 = history[[1]][code_rows],
	income_group_previous = label_group(history[[which(fiscal_years == "FY26")]][code_rows]),
	income_group_latest = label_group(history[[which(fiscal_years == "FY27")]][code_rows])
)
current = as.data.table(
	read_excel(path_of("class_current"), sheet = "List of economies")
)
current = current[grepl("^[A-Z]{3}$", Code) & !is.na(`Income group`), .(
	iso3 = Code,
	wb_economy = Economy,
	income_group_current_list = `Income group`
)]
countries = fromJSON(path_of("wb_countries"))[[2]]
countries = data.table(
	iso3 = countries$id,
	wb_name = countries$name,
	income_group_api = countries$incomeLevel$value,
	wb_region = countries$region$value
)
countries = countries[wb_region != "Aggregates"]
income = merge(
	merge(income, current, by = "iso3", all.x = TRUE),
	countries[, .(iso3, wb_name, income_group_api)],
	by = "iso3",
	all.x = TRUE
)
income_members = income[iso3 %in% members$iso3]
stopifnot(
	nrow(income_members) == 193L,
	income_members[!is.na(income_group_current_list), all(income_group_latest == income_group_current_list)],
	income_members[!is.na(income_group_api), all(income_group_latest == income_group_api)]
)

imf = read_analysis_table(path_of("imf_gdp"))[, .(iso3, gdp_2024_imf = ngdpd_2024)]
stopifnot(anyDuplicated(imf$iso3) == 0L)
wb_gdp = fromJSON(path_of("wb_gdp"))[[2]]
wb_gdp = data.table(
	iso3 = wb_gdp$countryiso3code,
	gdp_2024_wb = round(wb_gdp$value / 1e9, 3)
)[nzchar(iso3)]
stopifnot(anyDuplicated(wb_gdp$iso3) == 0L)

dac_lines = html_lines(path_of("dac_page"))
dac_names = dac_lines[(which(dac_lines == "Members:")[1] + 1):(which(dac_lines == "Associate:")[1] - 1)]
dac_names = dac_names[!grepl("^\\[", dac_names)]
dac_lookup = c(
	australia = "AUS",
	austria = "AUT",
	belgium = "BEL",
	canada = "CAN",
	czechia = "CZE",
	denmark = "DNK",
	estonia = "EST",
	finland = "FIN",
	france = "FRA",
	germany = "DEU",
	greece = "GRC",
	hungary = "HUN",
	iceland = "ISL",
	ireland = "IRL",
	italy = "ITA",
	japan = "JPN",
	korea = "KOR",
	latvia = "LVA",
	lithuania = "LTU",
	luxembourg = "LUX",
	netherlands = "NLD",
	newzealand = "NZL",
	norway = "NOR",
	poland = "POL",
	portugal = "PRT",
	slovakrepublic = "SVK",
	slovenia = "SVN",
	spain = "ESP",
	sweden = "SWE",
	switzerland = "CHE",
	unitedkingdom = "GBR",
	unitedstates = "USA",
	eu = "EU"
)
dac_codes = unname(dac_lookup[simplify_name(dac_names)])
stopifnot(length(dac_codes) == 33L, !anyNA(dac_codes), "EU" %in% dac_codes)
dac_members = setdiff(dac_codes, "EU")

scale_aliases = c(
	bolivia = "BOL",
	czechrepublic = "CZE",
	netherlands = "NLD",
	trkiye = "TUR",
	ctedivoire = "CIV",
	unitedkingdomofgreatbritainandnorthernireland = "GBR",
	unitedstatesofamerica = "USA"
)
match_members = function(name) {
	key = simplify_name(name)
	out = member_key$iso3[match(key, member_key$key)]
	out[is.na(out)] = unname(scale_aliases[key[is.na(out)]])
	out
}

clean_name = function(name) trimws(sub("\\s+[a-z]/?$", "", trimws(gsub("\\(formerly[^)]*\\)", "", name))))
scale_lines = readLines(path_of("scale_text"), warn = FALSE, encoding = "UTF-8")
scale_lines = scale_lines[grepl("^[A-Z]", scale_lines)]
scale_rows = data.table(
	name = clean_name(
		substr(
			scale_lines,
			1,
			regexpr("\\s{2,}(\\d|-)", scale_lines, perl = TRUE) - 1
		)
	),
	last_token = vapply(
		regmatches(
			scale_lines,
			gregexpr("\\d+\\.\\d+|(?<=\\s)-(?=\\s|$)", scale_lines, perl = TRUE)
		),
		function(tokens) tokens[length(tokens)],
		character(1)
	)
)
scale_rows = scale_rows[last_token != "-" & name != "Total"]
scale_rows[, share := as.numeric(last_token)]
scale_rows[, iso3 := match_members(name)]
stopifnot(
	!anyNA(scale_rows$iso3),
	nrow(scale_rows) == 193L,
	anyDuplicated(scale_rows$iso3) == 0L,
	abs(sum(scale_rows$share) - 100) < 0.005
)
assessment_lines = readLines(path_of("assessment_text"), warn = FALSE, encoding = "UTF-8")
assessment_pattern = "^\\s*([A-Za-z][^0-9]*?)\\s{2,}(\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s*$"
assessment_lines = assessment_lines[grepl(assessment_pattern, assessment_lines, perl = TRUE)]
assessment_rows = data.table(
	name = clean_name(sub(assessment_pattern, "\\1", assessment_lines, perl = TRUE)),
	share_2025_2027 = as.numeric(sub(assessment_pattern, "\\3", assessment_lines, perl = TRUE))
)
assessment_rows = assessment_rows[!grepl("^(Total|Grand total)", name)]
assessment_rows[, iso3 := match_members(name)]
stopifnot(!anyNA(assessment_rows$iso3), nrow(assessment_rows) == 193L)
scale_check = merge(
	scale_rows[, .(iso3, share)],
	assessment_rows[, .(iso3, share_2025_2027)],
	by = "iso3"
)
stopifnot(
	nrow(scale_check) == 193L,
	all(abs(scale_check$share - scale_check$share_2025_2027) < 1e-9)
)
# one row per member with the three criteria and the alternatives
selection = merge(members, treaty_summary, by = "iso3", all.x = TRUE)
selection = merge(
	selection,
	income_members[, .(iso3, wb_name, income_group_latest, income_group_previous)],
	by = "iso3",
	all.x = TRUE
)
selection = merge(selection, imf, by = "iso3", all.x = TRUE)
selection = merge(selection, wb_gdp, by = "iso3", all.x = TRUE)
selection = merge(
	selection,
	scale_rows[, .(iso3, un_assessment_share_2025_2027 = share)],
	by = "iso3",
	all.x = TRUE
)
selection[, dac_member := iso3 %in% dac_members]
selection[is.na(treaty_in_force_2001_2025), treaty_in_force_2001_2025 := FALSE]
selection[iso3 == "USA", `:=`(
	treaty = "United States",
	treaty_group = "United States",
	treaty_in_force_2001_2025 = FALSE
)]
selection[, criterion_treaty := treaty_in_force_2001_2025]
selection[, criterion_income := income_group_latest == "High income"]
selection[, criterion_gdp := !is.na(gdp_2024_imf) & gdp_2024_imf >= gdp_threshold]
selection[, selected_ally := criterion_treaty & criterion_income & criterion_gdp]
selection[, gdp_relative_difference := (gdp_2024_imf - gdp_2024_wb) / gdp_2024_wb]
stopifnot(setequal(selection[selected_ally == TRUE, iso3], expected_allies))
stopifnot(
	selection[!is.na(treaty) & gdp_2024_imf >= table_gdp_floor & !is.na(gdp_relative_difference), all(abs(gdp_relative_difference) < 0.05)]
)
# the same rule with the alternatives put in place of one criterion
rule_set = function(treaty_ok, income_ok, size_ok) {
	sort(
		selection$iso3[treaty_ok & income_ok & size_ok & selection$iso3 != "USA"]
	)
}

describe = function(selected) {
	data.table(
		n_selected = length(selected),
		selected = paste(selected, collapse = " "),
		added = paste(setdiff(selected, expected_allies), collapse = " "),
		dropped = paste(setdiff(expected_allies, selected), collapse = " "),
		same_fourteen = setequal(selected, expected_allies)
	)
}

variant_row = function(label, selected) cbind(variant = label, describe(selected))
treaty_ok = selection$criterion_treaty
income_ok = selection$criterion_income
size_ok = selection$criterion_gdp
variants = rbindlist(
	list(
		variant_row(
			"rule: treaty window 2001 to 2025, high income (latest), imf gdp at least 400",
			rule_set(treaty_ok, income_ok, size_ok)
		),
		variant_row(
			"world bank gdp at least 400 in place of imf gdp",
			rule_set(
				treaty_ok,
				income_ok,
				!is.na(selection$gdp_2024_wb) & selection$gdp_2024_wb >= gdp_threshold
			)
		),
		variant_row(
			"previous classification (july 2025) in place of the latest",
			rule_set(treaty_ok, selection$income_group_previous == "High income", size_ok)
		),
		variant_row(
			"dac membership in place of the income criterion",
			rule_set(treaty_ok, selection$dac_member, size_ok)
		),
		variant_row(
			"assessment share at least 0.5 percent in place of gdp",
			rule_set(
				treaty_ok,
				income_ok,
				!is.na(selection$un_assessment_share_2025_2027) & selection$un_assessment_share_2025_2027 >= share_threshold
			)
		),
		variant_row(
			"treaty window and high income without the size criterion",
			rule_set(treaty_ok, income_ok, TRUE)
		),
		variant_row(
			"treaty window and size without the income criterion",
			rule_set(treaty_ok, TRUE, size_ok)
		),
		variant_row(
			"high income and size without the treaty criterion",
			rule_set(TRUE, income_ok, size_ok)
		),
		variant_row(
			"treaty in force at the end of 2025 in place of the window",
			rule_set(
				!is.na(selection$treaty_in_force_year) & is.na(selection$treaty_end_year) & selection$iso3 != "USA",
				income_ok,
				size_ok
			)
		)
	)
)

candidates = selection[criterion_treaty & criterion_income]
interval = function(value) {
	in_group = candidates$iso3 %in% expected_allies
	included = value[in_group]
	excluded = value[!in_group & !is.na(value)]
	data.table(
		lower_exclusive = max(excluded),
		upper_inclusive = min(included),
		marginal_exclusion = candidates$iso3[!in_group][which.max(value[!in_group])],
		marginal_inclusion = candidates$iso3[in_group][which.min(value[in_group])]
	)
}

thresholds = rbindlist(
	list(
		cbind(
			quantity = "imf gdp 2024, billions of dollars",
			interval(candidates$gdp_2024_imf)
		),
		cbind(
			quantity = "world bank gdp 2024, billions of dollars",
			interval(candidates$gdp_2024_wb)
		),
		cbind(
			quantity = "regular budget assessment share 2025 to 2027, percent",
			interval(candidates$un_assessment_share_2025_2027)
		)
	)
)
start_years = rbindlist(
	lapply(
		1990:2026,
		function(start_year) {
			covers = treaties[, .(covers = any(covers_window(start_year, in_force_year, end_year))), by = iso3]
			window_ok = selection$iso3 %in% covers[covers == TRUE, iso3]
			cbind(
				window_start = start_year,
				describe(rule_set(window_ok, income_ok, size_ok))
			)
		}
	)
)
admissible_years = start_years[same_fourteen == TRUE, window_start]
stopifnot(
	window_start %in% admissible_years,
	all(diff(admissible_years) == 1L)
)
robustness = rbindlist(
	list(
		data.table(
			check = "gdp threshold interval (imf)",
			value = sprintf(
				"above %.3f and up to %.3f",
				thresholds$lower_exclusive[1],
				thresholds$upper_inclusive[1]
			),
			detail = sprintf(
				"marginal exclusion %s, marginal inclusion %s",
				thresholds$marginal_exclusion[1],
				thresholds$marginal_inclusion[1]
			)
		),
		data.table(
			check = "gdp threshold interval (world bank)",
			value = sprintf(
				"above %.3f and up to %.3f",
				thresholds$lower_exclusive[2],
				thresholds$upper_inclusive[2]
			),
			detail = sprintf(
				"marginal exclusion %s, marginal inclusion %s",
				thresholds$marginal_exclusion[2],
				thresholds$marginal_inclusion[2]
			)
		),
		data.table(
			check = "assessment share interval",
			value = sprintf(
				"above %.3f and up to %.3f",
				thresholds$lower_exclusive[3],
				thresholds$upper_inclusive[3]
			),
			detail = sprintf(
				"marginal exclusion %s, marginal inclusion %s",
				thresholds$marginal_exclusion[3],
				thresholds$marginal_inclusion[3]
			)
		),
		data.table(
			check = "admissible window start years",
			value = sprintf("%d to %d", min(admissible_years), max(admissible_years)),
			detail = sprintf(
				"%d drops %s; %d adds %s",
				min(admissible_years) - 1L,
				start_years[window_start == min(admissible_years) - 1L, dropped],
				max(admissible_years) + 1L,
				start_years[window_start == max(admissible_years) + 1L, added]
			)
		),
		variants[, .(
			check = variant,
			value = fifelse(same_fourteen, "same fourteen", sprintf("%d states", n_selected)),
			detail = trimws(
				paste(
					fifelse(nzchar(added), paste("adds", added), ""),
					fifelse(nzchar(dropped), paste("drops", dropped), "")
				)
			)
		)]
	)
)
# the criterion or criteria that exclude each treaty ally outside the group
exclusions = selection[!is.na(treaty) & iso3 != "USA" & selected_ally == FALSE]
exclusions[, excluded_by := trimws(
	paste(
		fifelse(criterion_treaty, "", "treaty window;"),
		fifelse(criterion_income, "", "income;"),
		fifelse(criterion_gdp, "", "gdp;")
	)
)]
exclusions[, excluded_by := sub(";$", "", gsub(";\\s+", "; ", excluded_by))]
exclusions[, n_failed := (!criterion_treaty) + (!criterion_income) + (!criterion_gdp)]
exclusions = exclusions[order(-gdp_2024_imf, na.last = TRUE), .(
	iso3,
	member_name,
	treaty_group,
	treaty_in_force_year,
	treaty_end_year,
	treaty_in_force_2001_2025,
	income_group_latest,
	gdp_2024_imf,
	gdp_2024_wb,
	dac_member,
	un_assessment_share_2025_2027,
	excluded_by,
	n_failed,
	treaty_note
)]
necessary = exclusions[n_failed == 1L, .(states = paste(iso3, collapse = " ")), by = .(criterion = excluded_by)]
stopifnot(setequal(necessary$criterion, c("treaty window", "income", "gdp")))

display_names = c(
	TUR = 'T\\"urkiye',
	KOR = "South Korea",
	GBR = "United Kingdom",
	USA = "United States",
	BHS = "Bahamas",
	VEN = "Venezuela",
	SVK = "Slovakia",
	CZE = "Czechia",
	NLD = "Netherlands",
	EGY = "Egypt",
	IRN = "Iran",
	RUS = "Russia",
	SYR = "Syria",
	LAO = "Laos",
	YEM = "Yemen"
)
selection[, display_name := fifelse(iso3 %in% names(display_names), display_names[iso3], wb_name)]
selection[is.na(display_name), display_name := member_name]
treaty_labels = c(
	NATO = "North Atlantic",
	ANZUS = "ANZUS",
	Japan = "Japan",
	Korea = "Korea",
	Philippines = "Philippines",
	`Manila Pact` = "Manila Pact",
	`Rio Treaty` = "Rio"
)
table_rows = selection[!is.na(treaty) & iso3 != "USA" & (selected_ally | gdp_2024_imf > table_gdp_floor)][order(-selected_ally, -gdp_2024_imf)]
table_rows[, in_force := fifelse(
	is.na(treaty_end_year),
	as.character(treaty_in_force_year),
	sprintf("%d to %d", treaty_in_force_year, treaty_end_year)
)]
table_lines = table_rows[, sprintf(
	"%s & %s & %s & %s & %s & %s\\\\",
	display_name,
	treaty_labels[treaty_group],
	in_force,
	sub(" income$", "", income_group_latest),
	formatC(round(gdp_2024_imf), format = "d", big.mark = ","),
	fifelse(selected_ally, "yes", "no")
)]
selection_table = c(
	"\\begin{table}[ht]",
	"\\centering\\small",
	sprintf(
		"\\caption{Selection of the fourteen allies. The fourteen states that meet the three criteria (a collective defense treaty with the United States in force continuously from %d through %d, high income in the World Bank classification of July 2026, and a 2024 gross domestic product of at least %d billion dollars in the October 2026 World Economic Outlook) and every other state that is or has been party to such a treaty and whose 2024 gross domestic product exceeds %d billion dollars. Treaty names the North Atlantic Treaty, the ANZUS Treaty, the bilateral mutual defense treaties with Japan, Korea and the Philippines, the Manila Pact and the Rio Treaty. GDP is in billions of current dollars. In force gives the year the obligation with the United States entered into force and, where it ended, the year it ended; for New Zealand the second year is the year the United States suspended its obligations under the treaty.}",
		window_start,
		window_end,
		gdp_threshold,
		table_gdp_floor
	),
	"\\label{tab:ally_selection}",
	"\\begin{tabularx}{\\textwidth}{>{\\raggedright\\arraybackslash}X >{\\raggedright\\arraybackslash}X l l c l}\\toprule",
	"State & Treaty & In force & Income group & GDP 2024 & Included\\\\\\midrule",
	table_lines,
	"\\bottomrule",
	"\\end{tabularx}",
	"\\end{table}"
)
# outputs
selection_out = selection[, .(
	iso3,
	member_name,
	display_name,
	treaty,
	treaty_group,
	treaty_in_force_year,
	treaty_end_year,
	treaty_in_force_2001_2025,
	income_group_latest,
	income_group_previous,
	gdp_2024_imf,
	gdp_2024_wb,
	gdp_relative_difference,
	dac_member,
	un_assessment_share_2025_2027,
	criterion_treaty,
	criterion_income,
	criterion_gdp,
	selected_ally,
	treaty_note
)][order(-selected_ally, -gdp_2024_imf, na.last = TRUE)]
fwrite(selection_out, file.path(output_root, "selection.csv"))
fwrite(robustness, file.path(output_root, "robustness.csv"))
fwrite(variants, file.path(output_root, "rule_variants.csv"))
fwrite(start_years, file.path(output_root, "window_start_years.csv"))
fwrite(thresholds, file.path(output_root, "threshold_intervals.csv"))
fwrite(exclusions, file.path(output_root, "exclusions.csv"))
fwrite(necessary, file.path(output_root, "necessary_criteria.csv"))
writeLines(selection_table, file.path(output_root, "selection_table.tex"))
# the new columns for data/member_groups.csv
member_columns = selection[, .(
	iso3,
	income_group_latest,
	income_group_previous,
	gdp_2024_imf,
	gdp_2024_wb,
	dac_member,
	un_assessment_share_2025_2027,
	treaty_in_force_2001_2025,
	selected_ally
)][order(iso3)]
fwrite(
	member_columns,
	file.path(root, "data", "ally_selection", "member_columns.csv")
)
merged = FALSE
if (merge_groups && file.exists(path_of("member_groups"))) {
	groups = read_analysis_table(path_of("member_groups"))
	stopifnot(
		setequal(groups$iso3, members$iso3),
		anyDuplicated(groups$iso3) == 0L
	)
	existing = setdiff(names(groups), setdiff(names(member_columns), "iso3"))
	groups = merge(
		groups[, existing, with = FALSE],
		member_columns,
		by = "iso3",
		all.x = TRUE,
		sort = FALSE
	)
	fwrite(groups, path_of("member_groups"))
	merged = TRUE
}

retrieval = fromJSON(path_of("retrieval"), simplifyVector = FALSE)
write_json(
	list(
		script = "code/04_ally_selection.R",
		run_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
		inputs = lapply(
			inputs[names(inputs) != "member_groups"],
			function(p) list(path = p, sha256 = digest(file = file.path(root, p), algo = "sha256"))
		),
		sources_retrieved = retrieval$created_utc,
		parameters = list(
			window_start = window_start,
			window_end = window_end,
			gdp_threshold_billions = gdp_threshold,
			share_threshold_percent = share_threshold,
			table_gdp_floor_billions = table_gdp_floor,
			income_latest = "World Bank fiscal year 2027 classification (July 2026), column FY27 of OGHIST.xlsx",
			income_previous = "World Bank fiscal year 2026 classification (July 2025), column FY26 of OGHIST.xlsx",
			gdp_imf = "IMF datamapper NGDPD, World Economic Outlook October 2026, 2024 values",
			gdp_wb = "World Bank NY.GDP.MKTP.CD, 2024, World Development Indicators last updated 2026-10-08",
			assessment = "regular budget scale 2025 to 2027, resolution 79/249"
		),
		rule = "a member is an ally when a collective defense treaty with the United States was in force continuously from 1 January of the window start year through the window end year, the latest World Bank classification is high income, and 2024 GDP in the IMF series is at least the threshold",
		members = nrow(selection),
		treaty_rows = nrow(treaties),
		treaty_states = nrow(treaty_summary),
		dac_members = length(dac_members),
		selected = sort(selection[selected_ally == TRUE, iso3]),
		expected = sort(expected_allies),
		member_groups_merged = merged,
		outputs = list.files(output_root),
		R = R.version.string
	),
	file.path(output_root, "provenance.json"),
	pretty = TRUE,
	auto_unbox = TRUE
)
