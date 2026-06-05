# download un voting data

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))

load_pkg(c('xml2', 'curl', 'data.table'))

# un bulk-vote record
record_url = 'https://digitallibrary.un.org/record/4060887?of=xm'
browser_agent = 'Mozilla/5.0'

# read marcxml
meta_handle = curl::new_handle()
curl::handle_setheaders(meta_handle, 'User-Agent' = browser_agent)
marc = xml2::read_xml(rawToChar(curl::curl_fetch_memory(record_url, handle = meta_handle)$content))

# find the bulk csv
fields_856 = xml2::xml_find_all(marc, ".//*[local-name()='datafield' and @tag='856']")
csv_url = NA_character_
csv_size = NA_real_
for(fld in fields_856){
	this_url = xml2::xml_text(xml2::xml_find_first(fld, "./*[local-name()='subfield' and @code='u']"))
	if(!is.na(this_url) && grepl('ga_voting\\.csv$', this_url)){
		csv_url = this_url
		csv_size = suppressWarnings(as.numeric(
			xml2::xml_text(xml2::xml_find_first(fld, "./*[local-name()='subfield' and @code='s']"))
		))
	}
}

# download if needed
csv_file = paste0(path_in, basename(csv_url))
is_cached = file.exists(csv_file) && (is.na(csv_size) || file.info(csv_file)$size == csv_size)
if(!is_cached){
	dl_handle = curl::new_handle()
	curl::handle_setheaders(dl_handle, 'User-Agent' = browser_agent)
	curl::curl_download(csv_url, csv_file, handle = dl_handle)
}

# keep used fields
keep_cols = c(
	'undl_id', 'ms_code', 'ms_name', 'ms_vote', 'date', 'session',
	'resolution', 'title', 'agenda_title', 'subjects',
	'total_yes', 'total_no', 'total_abstentions', 'total_non_voting', 'total_ms'
)
header = names(data.table::fread(csv_file, nrows = 0))
un_votes_long = data.table::fread(
	csv_file,
	select = intersect(keep_cols, header),
	encoding = 'UTF-8'
)

# standardize vote codes
un_votes_long[, ms_vote := toupper(trimws(ms_vote))]
un_votes_long[, date := as.Date(date)]
un_votes_long[, year := as.integer(format(date, '%Y'))]

# save coverage
un_votes_meta = list(
	csv_url = csv_url,
	file_size = file.info(csv_file)$size,
	n_rows = nrow(un_votes_long),
	year_range = range(un_votes_long$year, na.rm = TRUE),
	n_resolutions = length(unique(un_votes_long$resolution)),
	n_states = length(unique(un_votes_long$ms_code)),
	vote_codes = sort(unique(un_votes_long$ms_vote)),
	has_es11_7 = any(grepl('ES-11/7', un_votes_long$resolution))
)

save(un_votes_long, file = paste0(path_in, 'un_votes_long.rda'))
save(un_votes_meta, file = paste0(path_out, 'un_votes_meta.rda'))
