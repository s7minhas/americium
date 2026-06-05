# rebuild paper outputs

rm(list = ls())

Sys.setenv(AMERICIUM_REPO_ROOT = normalizePath(getwd(), mustWork = TRUE))
source(file.path(Sys.getenv('AMERICIUM_REPO_ROOT'), 'code', 'setup.R'))

source(file.path(path_code, '03_downstreamAnalyses', '12_paper_numbers.R'))
source(file.path(path_code, '03_downstreamAnalyses', '14_snap_stress_tests.R'))
source(file.path(path_code, '00_descriptives', '01_paper_figures.R'))
