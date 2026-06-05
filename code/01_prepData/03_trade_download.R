# download imf trade data

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))

load_pkg(c('imfapi', 'countrycode', 'data.table'))

# valid reporters
country_codes = as.data.frame(
	imf_get_codelists(dimension_ids = 'COUNTRY', dataflow_id = 'IMTS')
)$code
iso_countries = country_codes[
	!is.na(suppressWarnings(countrycode(country_codes, 'iso3c', 'country.name')))
]

# pull annual exports
trade_list = vector('list', length(iso_countries))
names(trade_list) = iso_countries
for(rep_code in iso_countries){
	trade_list[[rep_code]] = tryCatch(
		imf_get(
			'IMTS',
			dimensions = list(
				COUNTRY = rep_code,
				INDICATOR = 'XG_FOB_USD',
				FREQUENCY = 'A'
			),
			start_period = 1989,
			end_period = 2026
		),
		error = function(e) NULL
	)
}

# combine directed exports
trade_long = data.table::rbindlist(trade_list, fill = TRUE)
trade_long = trade_long[
	COUNTERPART_COUNTRY %in% iso_countries,
	.(
		reporter = COUNTRY,
		partner = COUNTERPART_COUNTRY,
		year = as.integer(TIME_PERIOD),
		exports_usd = as.numeric(OBS_VALUE)
	)
]
trade_long = trade_long[!is.na(exports_usd) & reporter != partner]

trade_meta = list(
	indicator = 'XG_FOB_USD',
	reporters_requested = length(iso_countries),
	n_reporters = length(unique(trade_long$reporter)),
	n_partners = length(unique(trade_long$partner)),
	year_range = range(trade_long$year, na.rm = TRUE),
	n_obs = nrow(trade_long)
)

save(trade_long, file = paste0(path_in, 'trade_long.rda'))
save(trade_meta, file = paste0(path_out, 'trade_meta.rda'))
