# heavy-tailed transition check

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))
source(paste0(path_funcs, 'project_metadata.R'))

load_pkg('lame')
set.seed(seed)

load(paste0(path_in, 'arr_list_lame.rda'))
un = arr_list_lame$un_voting
status_file = paste0(path_out, 'fattail_strain_status.log')
write_status = function(...){
	write(
		paste(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), paste(..., collapse = ' ')),
		file = status_file,
		append = TRUE
	)
}

# fit the heavy-tailed model
write_status('fitting heavy-tailed un model')
fit = lame(
	un, R = 2, family = 'normal', symmetric = TRUE,
	dynamic_uv = TRUE, dynamic_uv_kind = 't',
	nscan = 4000, burn = 2000, odens = 20,
	seed = seed, gof = FALSE, verbose = FALSE
)

# rank latent innovations
u = fit$U
years = dimnames(u)[[3]]
actors = rownames(u[, , 1])
rho = mean(fit$rho_uv)

# standardized ar innovations
innov = matrix(NA_real_, length(actors), length(years) - 1,
	dimnames = list(actors, years[-1]))
for(t in 2:length(years)){
	resid = u[, , t] - rho * u[, , t-1]
	innov[, t-1] = sqrt(rowSums(resid^2))
}

# rank the us 2025 move
all_innov = as.vector(innov)
us_2025 = innov['USA', '2025']
pct = mean(all_innov <= us_2025, na.rm = TRUE)
focal_actors = intersect(amer_focal_actors, rownames(innov))
focal_2025 = data.frame(
	actor = focal_actors,
	category = amer_actor_type(focal_actors),
	innov = as.numeric(innov[focal_actors, '2025']),
	percentile = sapply(focal_actors, function(a){
		mean(all_innov <= innov[a, '2025'], na.rm = TRUE)
	}),
	stringsAsFactors = FALSE
)
category_2025 = do.call(rbind, lapply(amer_category_labels, function(cat){
	x = focal_2025[focal_2025$category == cat, ]
	data.frame(
		category = cat,
		n = nrow(x),
		mean_percentile = mean(x$percentile, na.rm = TRUE),
		median_percentile = median(x$percentile, na.rm = TRUE),
		max_percentile = max(x$percentile, na.rm = TRUE),
		stringsAsFactors = FALSE
	)
}))

fattail_strain = list(
	us_2025_innov = us_2025,
	us_2025_percentile = pct,
	focal_2025 = focal_2025,
	category_2025 = category_2025,
	top10_innov = sort(all_innov, decreasing = TRUE)[1:10],
	us_innov_by_year = innov['USA', ],
	build = amer_build_metadata(list(comparison_model = 't_innovation_strain'))
)
save(fattail_strain, file = paste0(path_out, 'fattail_strain.rda'))
write_status('saved heavy-tailed strain')
