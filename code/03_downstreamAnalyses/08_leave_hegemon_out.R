# leave-out rupture checks

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))
source(paste0(path_funcs, 'project_metadata.R'))

load(paste0(path_out, 'snap_fit_un_voting.rda'))
sp = fit$snap_prob

# focal actor sets
order_actors = intersect(amer_focal_actors, rownames(sp))

# mean snap by actor set
index_for = function(set, yr){ mean(sp[set, as.character(yr)], na.rm = TRUE) }

years = 2018:2025
leave_hegemon_out = data.frame(
	year = years,
	full = sapply(years, function(y) index_for(order_actors, y)),
	excl_us = sapply(years, function(y) index_for(setdiff(order_actors, 'USA'), y)),
	excl_us_rus_chn = sapply(years, function(y) {
		index_for(setdiff(order_actors, c('USA', 'RUS', 'CHN')), y)
	}),
	excl_us_rivals = sapply(years, function(y) {
		index_for(setdiff(order_actors, c(amer_actor_cats$hegemon, amer_actor_cats$rivals)), y)
	})
)

save(leave_hegemon_out, file = paste0(path_out, 'leave_hegemon_out.rda'))
