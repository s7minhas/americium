# build lagged vectors by country-year

lagger = function(variable, country, year, lag_length){
	order_id = order(country, year)
	out = rep(NA, length(variable))
	
	country_ordered = as.character(country[order_id])
	variable_ordered = variable[order_id]
	
	lagged_country = c(rep(NA, lag_length), country_ordered[1:(length(country_ordered) - lag_length)])
	lagged_variable = c(rep(NA, lag_length), variable_ordered[1:(length(variable_ordered) - lag_length)])
	keep = country_ordered == lagged_country
	keep[is.na(keep)] = FALSE
	
	out_ordered = rep(NA, length(variable))
	out_ordered[keep] = lagged_variable[keep]
	out[order_id] = out_ordered
	out
}

# build lagged columns by country-year
multi_lagger = function(x, country, year, lag_length, relabel = TRUE){
	if(!is.data.frame(x)){
		stop('x must be a data frame')
	}
	
	lagged_x = x
	for(i in seq_len(ncol(x))){
		lagged_x[, i] = lagger(
			variable = x[, i],
			country = country,
			year = year,
			lag_length = lag_length
		)
	}
	
	if(relabel){
		names(lagged_x) = paste0(names(x), '.l', lag_length)
	}
	
	lagged_x
}
