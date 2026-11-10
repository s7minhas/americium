# read retained data with descriptive column names

describe_source_columns = function(table) {
	columns = c(
		ms_code = "voting_country_code",
		ms_name = "voting_country_name",
		ms_vote = "recorded_vote",
		CANTOT = "canada_import_price_index"
	)
	present = intersect(names(columns), names(table))
	if (length(present)) {
		data.table::setnames(table, present, unname(columns[present]))
	}
	table
}

read_analysis_table = function(...) {
	describe_source_columns(data.table::fread(...))
}

load_vote_archive = function(file, envir = parent.frame()) {
	base::load(file, envir = envir)
	votes = get("un_votes_long", envir = envir, inherits = FALSE)
	votes = describe_source_columns(data.table::as.data.table(votes))
	assign("un_votes_long", votes, envir = envir)
	invisible("un_votes_long")
}

read_census_table = function(...) {
	table = data.table::as.data.table(readxl::read_excel(...))
	columns = c(
		CTY_CODE = "source_country_code",
		CTYNAME = "source_country_name",
		CTY_DESC = "source_country_name",
		END_USE = "enduse_category_code",
		COMM_DESC = "enduse_category_description",
		IYR = "annual_imports_millions",
		EYR = "annual_exports_millions"
	)
	present = intersect(names(columns), names(table))
	data.table::setnames(table, present, unname(columns[present]))
	for (month in month.abb) {
		raw_name = paste0("I", toupper(month))
		if (raw_name %in% names(table)) {
			data.table::setnames(
				table,
				raw_name,
				paste0("imports_", tolower(month), "_millions")
			)
		}
	}
	table
}
