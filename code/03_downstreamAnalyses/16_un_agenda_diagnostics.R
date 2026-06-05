# un agenda checks

rm(list = ls())
setup_candidates = c(
	file.path(here::here(), 'setup.R'),
	file.path(here::here(), 'code', 'setup.R'),
	file.path(getwd(), 'setup.R'),
	file.path(getwd(), 'code', 'setup.R'),
	file.path(dirname(getwd()), 'code', 'setup.R')
)
setup_candidates = setup_candidates[file.exists(setup_candidates)]
if(length(setup_candidates) == 0){ stop('Could not locate setup.R.') }
source(setup_candidates[1])
source(paste0(path_funcs, 'project_metadata.R'))

load_pkg(c('data.table'))

load(paste0(path_in, 'un_votes_long.rda'))
setDT(un_votes_long)

votes = unique(un_votes_long[
	year >= 1989 & ms_vote %in% c('Y', 'N', 'A'),
	.(
		resolution, year, title, agenda_title, subjects,
		total_yes, total_no, total_abstentions, total_non_voting, total_ms
	)
])

classify_agenda = function(title, agenda_title, subjects){
	text = paste(title, agenda_title, subjects)
	out = rep('Other/institutional', length(text))
	out[grepl('DEVELOPMENT|POVERTY|FINANC|ECONOMIC|TRADE|TECHNOLOGY|ENVIRONMENT|CLIMATE|SUSTAINABLE|ISLAND|MOUNTAIN',
		text, ignore.case = TRUE)] = 'Development/economic'
	out[grepl('NUCLEAR|DISARMAMENT|WEAPON|ARMS|TERRORISM|SECURITY',
		text, ignore.case = TRUE)] = 'Security/disarmament'
	out[grepl('HUMAN RIGHTS|RIGHTS--|RACISM|MINORITIES|SELF-DETERMINATION|WOMEN',
		text, ignore.case = TRUE)] = 'Human rights/self-determination'
	out[grepl('PALESTIN|ISRAEL|MIDDLE EAST|UNRWA|JERUSALEM|GAZA',
		text, ignore.case = TRUE)] = 'Palestine/Israel'
	out[grepl('UKRAINE|RUSSIAN|CRIMEA',
		text, ignore.case = TRUE)] = 'Ukraine/Russia'
	out
}

votes[, present := total_yes + total_no + total_abstentions]
votes[, turnout := present / total_ms]
votes[, competitiveness := fifelse(
	total_yes + total_no > 0,
	1 - abs(total_yes - total_no) / (total_yes + total_no),
	0
)]
votes[, agenda_group := classify_agenda(title, agenda_title, subjects)]

vote_pos = copy(un_votes_long[
	year %in% c(2024, 2025) & ms_vote %in% c('Y', 'N', 'A'),
	.(year, resolution, ms_code, ms_vote)
])
pos_map = c(Y = 1, A = 0, N = -1)
vote_pos[, pos := pos_map[ms_vote]]
vote_pos = merge(
	vote_pos,
	votes[, .(year, resolution, agenda_group)],
	by = c('year', 'resolution'),
	all.x = TRUE
)

agenda_summary = votes[year %in% c(2024, 2025), .(
	resolutions = .N,
	mean_turnout = mean(turnout, na.rm = TRUE),
	mean_competitiveness = mean(competitiveness, na.rm = TRUE),
	ukraine_russia = sum(agenda_group == 'Ukraine/Russia'),
	palestine_israel = sum(agenda_group == 'Palestine/Israel'),
	human_rights = sum(agenda_group == 'Human rights/self-determination'),
	development = sum(agenda_group == 'Development/economic'),
	security = sum(agenda_group == 'Security/disarmament'),
	other = sum(agenda_group == 'Other/institutional')
), by = year][order(year)]

agenda_counts = votes[, .(resolutions = .N), by = year][order(-resolutions)]

mean_agreement = function(d, set_a, set_b, within = FALSE){
	set_a = intersect(set_a, unique(d$ms_code))
	set_b = intersect(set_b, unique(d$ms_code))
	if(length(set_a) == 0 || length(set_b) == 0){ return(NA_real_) }
	if(within){
		pairs = data.table::CJ(a = set_a, b = set_a)[a < b]
	} else {
		pairs = data.table::CJ(a = set_a, b = set_b)
	}
	if(nrow(pairs) == 0){ return(NA_real_) }
	da = d[ms_code %in% unique(pairs$a), .(resolution, a = ms_code, pos_a = pos)]
	db = d[ms_code %in% unique(pairs$b), .(resolution, b = ms_code, pos_b = pos)]
	x = merge(pairs, da, by = 'a', allow.cartesian = TRUE)
	x = merge(x, db, by = c('resolution', 'b'), allow.cartesian = TRUE)
	if(nrow(x) == 0){ return(NA_real_) }
	mean(1 - abs(x$pos_a - x$pos_b) / 2, na.rm = TRUE)
}

issue_groups = c(
	'Ukraine/Russia',
	'Palestine/Israel',
	'Human rights/self-determination',
	'Development/economic',
	'Security/disarmament',
	'Other/institutional'
)

issue_agreement = rbindlist(lapply(issue_groups, function(group){
	rbindlist(lapply(c(2024, 2025), function(yr){
		d = vote_pos[year == yr & agenda_group == group]
		data.frame(
			agenda_group = group,
			year = yr,
			resolutions = uniqueN(d$resolution),
			us_allies = mean_agreement(d, amer_actor_cats$hegemon, amer_actor_cats$allies),
			us_rivals = mean_agreement(d, amer_actor_cats$hegemon, amer_actor_cats$rivals),
			us_hedgers = mean_agreement(d, amer_actor_cats$hegemon, amer_actor_cats$hedgers),
			allies_internal = mean_agreement(d, amer_actor_cats$allies, amer_actor_cats$allies, within = TRUE),
			stringsAsFactors = FALSE
		)
	}))
}))

issue_agreement_wide = dcast(
	issue_agreement,
	agenda_group ~ year,
	value.var = c('resolutions', 'us_allies', 'us_rivals', 'us_hedgers', 'allies_internal')
)
issue_agreement_wide[, us_allies_change := us_allies_2025 - us_allies_2024]
issue_agreement_wide[, us_rivals_change := us_rivals_2025 - us_rivals_2024]
issue_agreement_wide[, us_hedgers_change := us_hedgers_2025 - us_hedgers_2024]
issue_agreement_wide[, allies_internal_change := allies_internal_2025 - allies_internal_2024]
issue_agreement_wide[, agenda_order := match(agenda_group, issue_groups)]
setorder(issue_agreement_wide, agenda_order)

save(
	agenda_summary, agenda_counts, issue_agreement, issue_agreement_wide,
	file = paste0(path_out, 'un_agenda_diagnostics.rda')
)
write.csv(
	agenda_summary,
	paste0(path_out, 'un_agenda_diagnostics_summary.csv'),
	row.names = FALSE
)
write.csv(
	issue_agreement_wide,
	paste0(path_out, 'un_issue_agreement_diagnostics.csv'),
	row.names = FALSE
)

fmt = function(x, digits = 3){ sprintf(paste0('%.', digits, 'f'), x) }

table_rows = c()
for(i in seq_len(nrow(agenda_summary))){
	row = agenda_summary[i]
	table_rows = c(table_rows, paste(
		row$year,
		row$resolutions,
		fmt(row$mean_turnout),
		fmt(row$mean_competitiveness),
		row$ukraine_russia,
		row$palestine_israel,
		row$human_rights,
		row$development,
		row$security,
		row$other,
		sep = ' & '
	) |> paste0(' \\\\'))
}

tex_lines = c(
	'\\begin{table}[!htbp]',
	'\\centering',
	'\\begingroup',
	'\\scriptsize',
	'\\caption{Descriptive UN-voting agenda diagnostics, 2024--2025. Agenda groups are rule-based summaries from resolution titles, agenda titles, and UN subject tags. The categories are diagnostic rather than mutually exhaustive substantive codings; the model uses the full annual voting network.}',
	'\\label{tab:un-agenda-diagnostics}',
	'\\begin{tabularx}{\\textwidth}{@{}rrrrrrrrrr@{}}',
	'\\toprule',
	'Year & Res. & Turnout & Compet. & Ukr./Rus. & Pal./Isr. & Rights & Dev./econ. & Sec./disarm. & Other \\\\',
	'\\midrule',
	table_rows,
	'\\bottomrule',
	'\\end{tabularx}',
	'\\endgroup',
	'\\end{table}'
)

writeLines(tex_lines, paste0(path_paper, 'generated_tables/un_agenda_diagnostics.tex'))

issue_rows = c()
for(i in seq_len(nrow(issue_agreement_wide))){
	row = issue_agreement_wide[i]
	issue_rows = c(issue_rows, paste(
		row$agenda_group,
		paste0(row$resolutions_2024, '/', row$resolutions_2025),
		fmt(row$us_allies_2024),
		fmt(row$us_allies_2025),
		fmt(row$us_allies_change),
		fmt(row$allies_internal_2024),
		fmt(row$allies_internal_2025),
		fmt(row$allies_internal_change),
		sep = ' & '
	) |> paste0(' \\\\'))
}

issue_tex_lines = c(
	'\\begin{table}[!htbp]',
	'\\centering',
	'\\begingroup',
	'\\scriptsize',
	'\\caption{Issue-stratified raw UN-voting agreement diagnostics, 2024--2025. Entries are Lijphart three-category agreement averaged within broad agenda groups. The diagnostic asks whether U.S.-ally agreement falls and allied internal agreement holds inside issue strata; it is not a snap-model refit.}',
	'\\label{tab:un-issue-agreement}',
	'\\begin{tabularx}{\\textwidth}{@{}>{\\raggedright\\arraybackslash}Xrrrrrrr@{}}',
	'\\toprule',
	'Agenda group & Res. 24/25 & U.S.-ally 24 & U.S.-ally 25 & Change & Ally internal 24 & Ally internal 25 & Change \\\\',
	'\\midrule',
	issue_rows,
	'\\bottomrule',
	'\\end{tabularx}',
	'\\endgroup',
	'\\end{table}'
)

writeLines(issue_tex_lines, paste0(path_paper, 'generated_tables/un_issue_agreement.tex'))
