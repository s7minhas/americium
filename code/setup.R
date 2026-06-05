# americium global setup

add_slash = function(x){
	if(!grepl('/$', x)){ x = paste0(x, '/') }
	x
}

# find the repo root
repo_candidates = c(
	Sys.getenv('AMERICIUM_REPO_ROOT', unset = NA_character_),
	here::here(),
	dirname(here::here()),
	getwd(),
	dirname(getwd()),
	'/home/s7m/Research/americium',
	'~/Research/americium'
)
repo_candidates = path.expand(repo_candidates[!is.na(repo_candidates) & nzchar(repo_candidates)])
repo_candidates = unique(normalizePath(repo_candidates, mustWork = FALSE))
repo_candidates = repo_candidates[file.exists(file.path(repo_candidates, 'code', 'setup.R'))]
if(length(repo_candidates) == 0){ stop('Could not locate americium repository root.') }

path_git = add_slash(repo_candidates[1])
path_paper = paste0(path_git, 'Paper/')
path_code = paste0(path_git, 'code/')
path_funcs = paste0(path_code, 'funcs/')

# data paths
drop_env = Sys.getenv('DROPBOX_AMERICIUM', unset = NA_character_)
if(!is.na(drop_env) && nzchar(drop_env)){
	path_drop = add_slash(path.expand(drop_env))
} else if(Sys.info()['user'] %in% c('s7m')){
	path_drop = '/mnt/d/Dropbox/americium/'
} else if(Sys.info()['user'] %in% c('sminhas')){
	path_drop = '/Users/sminhas/Library/CloudStorage/Dropbox/americium/'
} else {
	path_drop = path_git
}

path_in = paste0(path_drop, 'data/')
path_out = paste0(path_drop, 'results/')
path_graphics = paste0(path_drop, 'graphics/')
for(d in c(path_in, path_out, path_graphics)){
	dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# load packages
load_pkg = function(to_load){
	for(lib in to_load){
		if(identical(lib, 'lame')){
			lame_dev_path = path.expand(Sys.getenv(
				'LAME_DEV_PATH',
				unset = '~/Research/netify_dev/lame'
			))
			if(requireNamespace('pkgload', quietly = TRUE) && dir.exists(lame_dev_path)){
				pkgload::load_all(lame_dev_path, quiet = TRUE)
				next
			}
		}
		if(!(lib %in% installed.packages()[,1])){
			install.packages(lib, repos = 'http://cran.rstudio.com/')
		}
		suppressMessages(suppressWarnings(library(lib, character.only = TRUE)))
	}
}

to_load = c(
	'here',
	'countrycode',
	'reshape2', 'tidyr', 'dplyr', 'abind'
)
load_pkg(to_load)

# utility functions
cname = function(x){ countrycode(x, 'country.name', 'country.name') }
char = function(x){ as.character(x) }
num = function(x){ as.numeric(char(x)) }
trim = function(x){ gsub("^\\s+|\\s+$", "", x) }
substr_right = function(x, n){ substr(x, nchar(x)-n+1, nchar(x)) }
sum_na = function(x){
	if(all(is.na(x))){ NA } else { sum(x, na.rm = TRUE) }
}
extract_from_id = function(x, pos){
	out = strsplit(x, '_')
	out = lapply(out, function(z){ z[pos] })
	unlist(out)
}
paste_vec = function(x, y, sep = ''){
	combos = expand.grid(x, y)
	apply(combos, 1, paste, collapse = sep)
}
stdz = function(x){ (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) }
cbind_numb = function(x){ cbind(x, 1:length(x)) }
recode_qt = function(var, pseq){
	qts = quantile(var, probs = pseq, na.rm = TRUE)
	qts[1] = -Inf
	qts = unique(qts)
	qt_labs = paste0('cat_', 1:(length(qts)-1))
	cut(var, qts, qt_labs)
}
digs = 3
get_coef = function(x){
	round(cbind(mu = fixef(x), sd = sqrt(diag(vcov(x)))), digs)
}
get_coef_b = function(x, vars, int = FALSE){
	out = cbind(
		mu = apply(x[,vars], 2, mean),
		sd = apply(x[,vars], 2, sd)
	)
	if(int){
		out = cbind(
			out,
			qt_lo95 = apply(x[,vars], 2, quantile, .025),
			qt_hi95 = apply(x[,vars], 2, quantile, .975)
		)
	}
	round(out, digs)
}

# global parameters
seed = 6886
set.seed(seed)
