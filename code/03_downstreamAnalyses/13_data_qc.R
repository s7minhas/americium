# array qc

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

load(paste0(path_in, 'arr_list_lame.rda'))

summarize_matrix = function(m, relation, year){
	diag_idx = seq_len(min(nrow(m), ncol(m)))
	diag_vals = m[cbind(diag_idx, diag_idx)]
	obs = !is.na(m)
	obs[cbind(diag_idx, diag_idx)] = FALSE
	data.frame(
		relation = relation,
		year = as.integer(year),
		n_actors = nrow(m),
		observed_offdiag = sum(obs),
		fully_unobserved_actors = sum(rowSums(!is.na(m), na.rm = TRUE) +
			colSums(!is.na(m), na.rm = TRUE) == 0),
		nonmissing_diagonal = sum(!is.na(diag_vals)),
		stringsAsFactors = FALSE
	)
}

data_qc_by_year = do.call(rbind, unlist(lapply(names(arr_list_lame), function(rel){
	lapply(names(arr_list_lame[[rel]]), function(yr){
		summarize_matrix(arr_list_lame[[rel]][[yr]], rel, yr)
	})
}), recursive = FALSE))

data_qc_summary = do.call(rbind, lapply(split(data_qc_by_year, data_qc_by_year$relation), function(d){
	data.frame(
		relation = d$relation[1],
		years = paste(min(d$year), max(d$year), sep = '--'),
		n_actors = max(d$n_actors),
		min_observed_offdiag = min(d$observed_offdiag),
		median_observed_offdiag = median(d$observed_offdiag),
		max_fully_unobserved_actors = max(d$fully_unobserved_actors),
		nonmissing_diagonal_total = sum(d$nonmissing_diagonal),
		stringsAsFactors = FALSE
	)
}))
rownames(data_qc_summary) = NULL

data_qc = list(
	by_year = data_qc_by_year,
	summary = data_qc_summary,
	actor_categories = amer_actor_category_table,
	build = amer_build_metadata()
)

save(data_qc, file = paste0(path_out, 'data_qc.rda'))
write.csv(data_qc_by_year, paste0(path_out, 'data_qc_by_year.csv'), row.names = FALSE)
write.csv(data_qc_summary, paste0(path_out, 'data_qc_summary.csv'), row.names = FALSE)
