# summarize hedger snaps

rm(list = ls())
pth = paste0(here::here(), '/')
source(paste0(pth, 'setup.R'))
source(paste0(path_funcs, 'project_metadata.R'))

load(paste0(path_out, 'snap_fit_un_voting.rda'))
sp = fit$snap_prob

hedgers = intersect(amer_actor_cats$hedgers, rownames(sp))
allies = intersect(amer_actor_cats$allies, rownames(sp))

pivotal_hedgers = list(
	hedger_snap = data.frame(
		hedger = hedgers,
		snap_2024 = sp[hedgers, '2024'],
		snap_2025 = sp[hedgers, '2025']
	),
	mean_hedger_2025 = mean(sp[hedgers, '2025']),
	mean_ally_2025 = mean(sp[allies, '2025']),
	us_2025 = sp['USA', '2025']
)

save(pivotal_hedgers, file = paste0(path_out, 'pivotal_hedgers.rda'))
