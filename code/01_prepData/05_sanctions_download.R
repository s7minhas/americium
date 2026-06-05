# download opensanctions data

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))

load_pkg(c('jsonlite', 'curl', 'countrycode', 'data.table'))

# map program prefixes
sender_lookup = c(
	US = 'USA', EU = 'EU', UN = 'UN', GB = 'GBR', CA = 'CAN', AU = 'AUS',
	UA = 'UKR', NZ = 'NZL', JP = 'JPN', CH = 'CHE', SECO = 'CHE', FR = 'FRA',
	BE = 'BEL', RU = 'RUS', CN = 'CHN', AE = 'ARE', IR = 'IRN', PK = 'PAK'
)
map_sender = function(prefix){
	out = sender_lookup[prefix]
	fallback = suppressWarnings(countrycode(prefix, 'iso2c', 'iso3c'))
	out[is.na(out)] = fallback[is.na(out)]
	unname(out)
}

# download targets file
nested_url = 'https://data.opensanctions.org/datasets/latest/sanctions/targets.nested.json'
nested_file = paste0(path_in, 'opensanctions_targets_nested.json')
if(!file.exists(nested_file)){
	dl = curl::new_handle()
	curl::handle_setheaders(dl, 'User-Agent' = 'Mozilla/5.0')
	curl::curl_download(nested_url, nested_file, handle = dl)
}

# stream designation records
con = file(nested_file, 'r')
chunk_tables = list()
chunk_size = 2000
repeat {
	lines = readLines(con, n = chunk_size, warn = FALSE)
	if(length(lines) == 0){ break }
	lines = sub(',$', '', trimws(lines))
	lines = lines[nchar(lines) > 1]
	chunk_sender = character(0)
	chunk_target = character(0)
	chunk_date = character(0)
	for(line in lines){
		obj = tryCatch(jsonlite::fromJSON(line, simplifyVector = FALSE), error = function(e) NULL)
		if(is.null(obj)){ next }
		target_iso2 = unlist(obj$properties$country)[1]
		sanctions = obj$properties$sanctions
		if(length(target_iso2) == 0 || is.na(target_iso2)){ next }
		if(is.null(sanctions) || length(sanctions) == 0){ next }
		for(s in sanctions){
			program_id = unlist(s$properties$programId)[1]
			start_date = unlist(s$properties$startDate)[1]
			if(length(program_id) == 0 || is.na(program_id)){ next }
			chunk_sender = c(chunk_sender, sub('-.*', '', program_id))
			chunk_target = c(chunk_target, target_iso2)
			chunk_date = c(chunk_date, if(length(start_date) == 0) NA_character_ else start_date)
		}
	}
	if(length(chunk_sender) > 0){
		chunk_tables[[length(chunk_tables) + 1]] = data.table::data.table(
			sender_prefix = chunk_sender,
			target_iso2 = chunk_target,
			start_date = chunk_date
		)
	}
}
close(con)

# resolve dyad years
sanctions_long = data.table::rbindlist(chunk_tables)
sanctions_long[, sender := map_sender(sender_prefix)]
sanctions_long[, target := suppressWarnings(countrycode(toupper(target_iso2), 'iso2c', 'iso3c'))]
sanctions_long[, year := as.integer(substr(start_date, 1, 4))]
sanctions_long = sanctions_long[!is.na(sender) & !is.na(target) & sender != target]

sanctions_meta = list(
	source = nested_url,
	cached_file = nested_file,
	cached_file_size = file.info(nested_file)$size,
	cached_file_mtime = as.character(file.info(nested_file)$mtime),
	cached_file_md5 = as.character(tools::md5sum(nested_file)),
	n_records = nrow(sanctions_long),
	n_dated = sum(!is.na(sanctions_long$year)),
	year_range = range(sanctions_long$year, na.rm = TRUE),
	n_senders = length(unique(sanctions_long$sender)),
	n_targets = length(unique(sanctions_long$target))
)

save(sanctions_long, file = paste0(path_in, 'sanctions_long.rda'))
save(sanctions_meta, file = paste0(path_out, 'sanctions_meta.rda'))
