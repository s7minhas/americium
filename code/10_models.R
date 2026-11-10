# prepare and fit the reported network models

helper_file = if (sys.nframe() > 0L && !is.null(sys.frame(1)$ofile)) {
	sys.frame(1)$ofile
} else {
	sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
}

source(
	file.path(dirname(normalizePath(helper_file)), "00_data_helpers.R"),
	local = TRUE
)
model_priors = function() {
	list(
		full = list(
			beta_precision = 0.99999999999999911,
			va_shape = 5,
			va_scale = 0.043620075897986083,
			s2_shape = 0.5,
			s2_scale = 0.5,
			rho_mean = 0.90000000000000002,
			rho_sd = 0.10000000000000001,
			sigma_shape = 2,
			sigma_scale = 1,
			snap_kappa = 2,
			snap_pi_a = 1,
			snap_pi_b = 999
		),
		no_us = list(
			beta_precision = 0.99999999999999911,
			va_shape = 5,
			va_scale = 0.041242807336404994,
			s2_shape = 0.5,
			s2_scale = 0.5,
			rho_mean = 0.90000000000000002,
			rho_sd = 0.10000000000000001,
			sigma_shape = 2,
			sigma_scale = 1,
			snap_kappa = 2,
			snap_pi_a = 1,
			snap_pi_b = 999
		),
		no_us_reweighted = list(
			beta_precision = 0.99999999999999911,
			va_shape = 5,
			va_scale = 0.041242807336404994,
			s2_shape = 0.5,
			s2_scale = 0.5,
			rho_mean = 0.90000000000000002,
			rho_sd = 0.10000000000000001,
			sigma_shape = 2,
			sigma_scale = 1,
			snap_kappa = 2,
			snap_pi_a = 1,
			snap_pi_b = 999
		),
		copy_2024 = list(
			beta_precision = 0.99999999999999911,
			va_shape = 5,
			va_scale = 0.043620075897986083,
			s2_shape = 0.5,
			s2_scale = 0.5,
			rho_mean = 0.90000000000000002,
			rho_sd = 0.10000000000000001,
			sigma_shape = 2,
			sigma_scale = 1,
			snap_kappa = 2,
			snap_pi_a = 1,
			snap_pi_b = 999
		),
		us_only_change = list(
			beta_precision = 0.99999999999999911,
			va_shape = 5,
			va_scale = 0.043620075897986083,
			s2_shape = 0.5,
			s2_scale = 0.5,
			rho_mean = 0.90000000000000002,
			rho_sd = 0.10000000000000001,
			sigma_shape = 2,
			sigma_scale = 1,
			snap_kappa = 2,
			snap_pi_a = 1,
			snap_pi_b = 999
		),
		trade = list(
			beta_precision = 0.059793887403341399,
			va_shape = 5,
			va_scale = 20.174490096939678,
			s2_shape = 0.5,
			s2_scale = 0.5,
			rho_mean = 0.90000000000000002,
			rho_sd = 0.10000000000000001,
			sigma_shape = 2,
			sigma_scale = 1,
			snap_kappa = 2,
			snap_pi_a = 1,
			snap_pi_b = 999
		),
		controls = list(
			beta_precision = 1,
			va_shape = 5,
			va_scale = 0.040000000000000001,
			s2_shape = 0.5,
			s2_scale = 0.5,
			rho_mean = 0.90000000000000002,
			rho_sd = 0.10000000000000001,
			sigma_shape = 2,
			sigma_scale = 1,
			snap_kappa = 2,
			snap_pi_a = 1,
			snap_pi_b = 999
		)
	)
}

# generate paired controls without using empirical reset estimates
simulate_control = function(
	seed,
	sigma,
	event_period,
	case,
	reset_sd = .5,
	actors = 40L,
	periods = 12L,
	rank = 2L
) {
	cases = c(
		"drift",
		"actor_one_reset",
		"simultaneous_25",
		"simultaneous_50",
		"asynchronous_25",
		"bloc_rotation_25",
		"common_level",
		"global_rotation",
		"agenda_only"
	)
	stopifnot(
		case %in% cases,
		sigma %in% c(.05, .10),
		event_period %in% c(7L, 12L),
		actors == 40L,
		periods == 12L,
		rank == 2L,
		reset_sd %in% c(.5, 2),
		length(seed) == 1L,
		is.finite(seed)
	)
	RNGkind("Mersenne-Twister", "Inversion", "Rejection")
	set.seed(seed)
	actor_names = c("USA", sprintf("A%02d", 2:actors))
	year_names = as.character(2001:(2000 + periods))
	dimnames_u = list(actor_names, paste0("u", 1:rank), year_names)
	dimnames_y = list(actor_names, actor_names, year_names)
	rho = .95
	residual_sd = .05
	intercept = .5
	additive = rnorm(actors, 0, .03)
	initial = matrix(rnorm(actors * rank), actors, rank) * sigma / sqrt(1 - rho^2)
	innovations = array(rnorm(actors * rank * periods), c(actors, rank, periods))
	reset_positions = matrix(rnorm(actors * rank), actors, rank) * reset_sd
	noise = array(0, c(actors, actors, periods), dimnames = dimnames_y)
	mask = array(FALSE, dim(noise), dimnames = dimnames_y)
	upper = upper.tri(noise[, , 1])
	for (t in seq_len(periods)) {
		e = matrix(0, actors, actors)
		e[upper] = rnorm(sum(upper), 0, residual_sd)
		noise[, , t] = e + t(e)
		m = matrix(FALSE, actors, actors)
		m[upper] = runif(sum(upper)) < .9
		mask[, , t] = m | t(m)
	}
	mask[(actors - 1):actors, , 1:2] = FALSE
	mask[, (actors - 1):actors, 1:2] = FALSE
	mask[actors, , periods] = FALSE
	mask[, actors, periods] = FALSE
	baseline_u = array(0, c(actors, rank, periods), dimnames = dimnames_u)
	baseline_u[, , 1] = initial
	for (t in 2:periods) {
		baseline_u[, , t] = rho * baseline_u[, , t - 1] + sigma * innovations[, , t]
	}
	mean_network = function(u, beta) {
		out = array(NA_real_, c(actors, actors, periods), dimnames = dimnames_y)
		for (t in seq_len(periods)) {
			out[, , t] = beta[t] + outer(additive, additive, "+") + tcrossprod(u[, , t])
			diag(out[, , t]) = NA_real_
		}
		out
	}
	baseline_mean = mean_network(baseline_u, rep(intercept, periods))
	baseline_y = baseline_mean + noise
	baseline_y[!mask] = NA_real_
	u = baseline_u
	beta = rep(intercept, periods)
	delta = matrix(0L, actors, periods, dimnames = list(actor_names, year_names))
	jump = delta
	level_change = rep(FALSE, periods)
	agenda_change = rep(FALSE, periods)
	q = matrix(c(cos(.8), sin(.8), -sin(.8), cos(.8)), 2, 2)
	alternative_u = NULL
	agenda = NULL
	checks = c(
		symmetric_noise = max(abs(noise - aperm(noise, c(2, 1, 3)))) < 1e-14,
		symmetric_mask = identical(mask, aperm(mask, c(2, 1, 3)))
	)
	selected = switch(
	case,
		actor_one_reset = 1L,
		simultaneous_25 = 1:10,
		simultaneous_50 = 1:20,
		asynchronous_25 = 1:10,
		integer()
	)
	if (length(selected)) {
		schedule = if (case == "asynchronous_25") {
			rep((event_period - 3):event_period, length.out = length(selected))
		} else {
			rep(event_period, length(selected))
		}
		delta[cbind(selected, schedule)] = 1L
		jump = delta
		for (t in 2:periods) {
			u[, , t] = rho * u[, , t - 1] + sigma * innovations[, , t]
			reset = which(delta[, t] == 1L)
			if (length(reset)) {
				u[reset, , t] = reset_positions[reset, , drop = FALSE]
			}
		}
	}
	if (case == "bloc_rotation_25") {
		for (t in event_period:periods) {
			u[1:10, , t] = baseline_u[1:10, , t] %*% q
		}
		jump[1:10, event_period] = 1L
		checks["within_bloc_gram"] = all(
			vapply(
				event_period:periods,
				function(t) {
					max(abs(tcrossprod(u[1:10, , t]) - tcrossprod(baseline_u[1:10, , t]))) < 1e-12
				},
				logical(1)
			)
		)
		checks["outside_bloc_unchanged"] = identical(u[-(1:10), , ], baseline_u[-(1:10), , ])
	}
	if (case == "common_level") {
		beta[event_period:periods] = beta[event_period:periods] + .25
		level_change[event_period] = TRUE
	}
	if (case == "global_rotation") {
		alternative_u = baseline_u
		for (t in event_period:periods) {
			alternative_u[, , t] = baseline_u[, , t] %*% q
		}
		checks["global_gram_unchanged"] = all(
			vapply(
				seq_len(periods),
				function(t) {
					max(abs(tcrossprod(alternative_u[, , t]) - tcrossprod(baseline_u[, , t]))) < 1e-12
				},
				logical(1)
			)
		)
	}
	mu = mean_network(u, beta)
	y = mu + noise
	y[!mask] = NA_real_
	if (case == "global_rotation") {
		y = baseline_y
		checks["global_y_bitwise_unchanged"] = identical(y, baseline_y)
	}
	if (case == "agenda_only") {
		preferences = matrix(rnorm(actors * rank, 0, .6), actors, rank)
		loadings = matrix(rnorm(200 * rank), 200, rank)
		scores = preferences %*% t(loadings) + matrix(rep(rnorm(200, 0, .25), each = actors), actors, 200)
		contested = ifelse(scores > .25, 1, ifelse(scores < -.25, -1, 0))
		votes = cbind(contested, matrix(1, actors, 200))
		dimnames(votes) = list(actor_names, paste0("item_", 1:400))
		agreement = function(items) {
			positions = votes[, items, drop = FALSE]
			abstain = 1 * (positions == 0)
			out = .5 * (1 + (tcrossprod(positions) + tcrossprod(abstain)) / length(items))
			diag(out) = NA_real_
			out
		}
		reference = actual = vector("list", periods)
		reference_mean = array(NA_real_, dim(y), dimnames = dimnames_y)
		for (t in seq_len(periods)) {
			left = sample.int(200)
			right = 200 + sample.int(200)
			reference[[t]] = c(left[1:80], right[1:20])
			actual[[t]] = if (t < event_period) reference[[t]] else c(left[1:20], right[1:80])
			reference_mean[, , t] = agreement(reference[[t]])
			mu[, , t] = agreement(actual[[t]])
		}
		y = mu + noise
		y[!mask] = NA_real_
		reference_y = reference_mean + noise
		reference_y[!mask] = NA_real_
		agenda_change[event_period] = TRUE
		agenda = list(
			fixed_preferences = preferences,
			fixed_item_votes = votes,
			selected_items = actual,
			reference_items = reference,
			reference_mean = reference_mean,
			reference_y = reference_y,
			contested_before = 80L,
			contested_after = 20L,
			high_agreement_before = 20L,
			high_agreement_after = 80L,
			statement = "Every actor-item vote is fixed; only the sampled item mix changes. Gaussian measurement errors remain paired."
		)
		checks["agenda_before_event_unchanged"] = identical(
			y[, , seq_len(event_period - 1)],
			reference_y[, , seq_len(event_period - 1)]
		)
		checks["agenda_fixed_votes"] = all(votes %in% c(-1, 0, 1))
	}
	checks["symmetric_y"] = identical(y, aperm(y, c(2, 1, 3)))
	checks["mask_retained"] = identical(is.finite(y), mask)
	checks["reset_count"] = sum(delta) == length(selected)
	checks["first_year_no_reset"] = all(delta[, 1] == 0)
	stopifnot(all(checks))
	list(
		y = y,
		mask = mask,
		truth = list(
			u = if (case == "agenda_only") NULL else u,
			baseline_u = baseline_u,
			alternative_coordinates = alternative_u,
			standard_innovations = innovations,
			independent_reset_positions = reset_positions,
			initial_positions = initial,
			intercept = beta,
			additive = additive,
			mean = mu,
			noise = noise,
			reset_indicator = delta,
			jump_indicator = jump,
			reset_indicator_interpretable = !case %in% c("bloc_rotation_25", "agenda_only"),
			level_change = level_change,
			agenda_change = agenda_change,
			baseline_y = baseline_y,
			agenda = agenda
		),
		specification = list(
			seed = seed,
			case = case,
			actors = actors,
			periods = periods,
			rank = rank,
			event_period = event_period,
			event_year = as.integer(year_names[event_period]),
			sigma = sigma,
			residual_sd = residual_sd,
			rho = rho,
			reset_sd = reset_sd,
			initial_distribution = "stationary Gaussian drift",
			static_intercept = intercept,
			additive_sd = .03,
			rotation_angle = .8,
			common_level_jump = .25,
			observed_dyad_probability = .9,
			rng_kind = RNGkind(),
			schema = "paired_controls_v1"
		),
		checks = checks,
		generation_rng = .Random.seed
	)
}

# reconstruct all reported model inputs from retained raw observations
prepare_un = function(root) {
	weighted = build_weighted_networks(file.path(root, "data"))
	un = weighted$original
	other = setdiff(rownames(un[[1]]), "USA")
	common = is.finite(un[["2024"]]) & is.finite(un[["2025"]])
	copy = un
	copy[["2025"]] = un[["2024"]]
	copy[["2025"]][!common] = NA_real_
	us_only = copy
	us_only[["2025"]]["USA", ] = un[["2025"]]["USA", ]
	us_only[["2025"]][, "USA"] = un[["2025"]][, "USA"]
	us_only[["2025"]][!common] = NA_real_
	stopifnot(
		identical(copy[["2025"]][other, other], us_only[["2025"]][other, other])
	)
	lapply(
		list(
			full = un,
			no_us = lapply(un, function(m) m[other, other]),
			no_us_reweighted = weighted$without_us,
			copy_2024 = copy,
			us_only_change = us_only
		),
		simplify2array
	)
}

prepare_trade = function(root) {
	input = new.env(parent = emptyenv())
	load(file.path(root, "data", "trade.rda"), input)
	trade = input$trade_long
	actors = readLines(file.path(root, "data", "actors.txt"))
	exports = setNames(
		lapply(
			1989:2025,
			function(year_value) {
				m = matrix(
					NA_real_,
					length(actors),
					length(actors),
					dimnames = list(actors, actors)
				)
				x = trade[trade$year == year_value & trade$reporter %in% actors & trade$partner %in% actors, ]
				m[cbind(match(x$reporter, actors), match(x$partner, actors))] = x$exports_usd
				diag(m) = NA_real_
				m
			}
		),
		1989:2025
	)
	sum_directions = function(m) {
		observed = is.finite(m) | is.finite(t(m))
		m[!is.finite(m)] = 0
		out = m + t(m)
		out[!observed] = NA_real_
		diag(out) = NA_real_
		out
	}
	raw = lapply(exports, sum_directions)
	both = lapply(
		exports,
		function(m) {
			out = sum_directions(m)
			out[!is.finite(m) | !is.finite(t(m))] = NA_real_
			out
		}
	)
	common = raw
	mask = is.finite(exports[["2024"]]) & is.finite(exports[["2025"]])
	for (year in c("2024", "2025")) {
		m = exports[[year]]
		m[!mask] = NA_real_
		common[[year]] = sum_directions(m)
	}
	copy = raw
	copy[["2025"]] = raw[["2024"]]
	missing = copy
	m = exports[["2024"]]
	m[!is.finite(exports[["2025"]])] = NA_real_
	missing[["2025"]] = sum_directions(m)
	planted = raw
	before = after = raw[["2025"]]
	targets = c("USA", "CHN", "DEU", "JPN", "BRA", "IND", "TZA", "AGO")
	donors = setdiff(actor_groups$allies, targets)
	receivers = setdiff(c(actor_groups$rivals, actor_groups$hedgers), targets)
	design = list()
	for (actor in targets) {
		d = donors[is.finite(before[actor, donors])]
		r = receivers[is.finite(before[actor, receivers])]
		transferred = .5 * sum(before[actor, d])
		multiplier = 1 + transferred / sum(before[actor, r])
		after[actor, d] = .5 * before[actor, d]
		after[actor, r] = before[actor, r] * multiplier
		after[d, actor] = after[actor, d]
		after[r, actor] = after[actor, r]
		design[[actor]] = data.frame(
			actor = actor,
			transferred_usd = transferred,
			share_reallocated = transferred / sum(before[actor, ], na.rm = TRUE),
			total_before = sum(before[actor, ], na.rm = TRUE),
			total_after = sum(after[actor, ], na.rm = TRUE)
		)
	}
	other = setdiff(actors, targets)
	stopifnot(
		identical(is.finite(before), is.finite(after)),
		identical(after, t(after)),
		identical(before[other, other], after[other, other]),
		max(
			abs(
				rowSums(before[targets, ], na.rm = TRUE) / rowSums(after[targets, ], na.rm = TRUE) - 1
			)
		) < 1e-12
	)
	planted[["2025"]] = after
	inputs = lapply(
		list(
			current_baseline = raw,
			both_directions = both,
			common_terminal_directions = common,
			copy_2024_terminal = copy,
			missingness_terminal = missing,
			partner_reallocation_terminal = planted
		),
		function(x) simplify2array(lapply(x, log))
	)
	fresh = readRDS(file.path(root, "results", "trade", "trade_refresh_inputs.rds"))
	inputs = c(inputs, lapply(fresh, simplify2array))
	write.csv(
		do.call(rbind, design),
		file.path(root, "results", "models", "trade_reallocation.csv"),
		row.names = FALSE
	)
	inputs
}

prepare_controls = function(root) {
	cases = c(
		"drift",
		"actor_one_reset",
		"simultaneous_25",
		"simultaneous_50",
		"asynchronous_25",
		"bloc_rotation_25",
		"common_level",
		"global_rotation",
		"agenda_only"
	)
	grid = expand.grid(
		seed = c(6886L, 19421L, 36583L),
		sigma = c(.05, .10),
		event_period = c(7L, 12L),
		case = cases,
		stringsAsFactors = FALSE
	)
	grid$reset_sd = .5
	matched = expand.grid(
		seed = 6886L,
		sigma = c(.05, .10),
		event_period = c(7L, 12L),
		case = c("simultaneous_25", "simultaneous_50"),
		stringsAsFactors = FALSE
	)
	matched$reset_sd = 2
	grid = rbind(grid, matched)
	inputs = list()
	truth = list()
	rows = list()
	for (i in seq_len(nrow(grid))) {
		generated = do.call(simulate_control, as.list(grid[i, ]))
		hash = digest::digest(generated$y, algo = "sha256")
		key = paste0("input_", substr(hash, 1, 12))
		if (key %in% names(inputs)) stopifnot(identical(inputs[[key]], generated$y)) else inputs[[key]] = generated$y
		truth[[as.character(i)]] = generated[c("truth", "specification", "checks")]
		rows[[i]] = data.frame(
			case_id = i,
			grid[i, ],
			specification = key,
			data_hash = hash,
			mask_hash = digest::digest(generated$mask, algo = "sha256"),
			checks_passed = all(generated$checks)
		)
	}
	mapping = do.call(rbind, rows)
	stopifnot(
		nrow(mapping) == 116L,
		length(inputs) == 92L,
		all(mapping$checks_passed)
	)
	write.csv(
		mapping,
		file.path(root, "results", "models", "control_cases.csv"),
		row.names = FALSE
	)
	saveRDS(truth, file.path(root, "results", "models", "control_truth.rds"))
	inputs
}

prepare_models = function(root) {
	dir.create(
		file.path(root, "results", "models"),
		recursive = TRUE,
		showWarnings = FALSE
	)
	inputs = list(
		un = prepare_un(root),
		trade = prepare_trade(root),
		controls = prepare_controls(root)
	)
	saveRDS(inputs, file.path(root, "results", "models", "inputs.rds"))
	manifest = do.call(
		rbind,
		lapply(
			names(inputs),
			function(domain) {
				do.call(
					rbind,
					lapply(
						names(inputs[[domain]]),
						function(specification) {
							y = inputs[[domain]][[specification]]
							data.frame(
								domain = domain,
								specification = specification,
								actors = dim(y)[1],
								periods = dim(y)[3],
								input_sha256 = digest::digest(y, algo = "sha256")
							)
						}
					)
				)
			}
		)
	)
	write.csv(
		manifest,
		file.path(root, "results", "models", "input_manifest.csv"),
		row.names = FALSE
	)
	invisible(inputs)
}

# rotation-invariant summaries of the stored trajectories
summary_backend = function(root) {
	backend = new.env(parent = baseenv())
	Rcpp::cppFunction(
		env = backend,
		cacheDir = file.path(root, "code", ".summary_cache"),
		code = '
Rcpp::List trajectory_geometry(Rcpp::NumericVector u, Rcpp::NumericVector y) {
 Rcpp::IntegerVector dim=u.attr("dim"); int n=dim[0],r=dim[1],T=dim[2],S=dim[3];
 Rcpp::NumericMatrix energy(S,T), centered(S,T), gram(S,T), distance(S,T);
 for(int t=0;t<T;++t) {
  std::vector<int> left,right;
  for(int j=1;j<n;++j)for(int i=0;i<j;++i)if(R_finite(y[i+n*j+n*n*t])){left.push_back(i);right.push_back(j);}
  for(int s=0;s<S;++s){
   const long offset=(long)n*r*(t+T*s); double e=0,meansquare=0;
   for(int k=0;k<r;++k){double sum=0;for(int i=0;i<n;++i){double x=u[offset+i+n*k];sum+=x;e+=x*x;}meansquare+=(sum/n)*(sum/n);}
   energy(s,t)=e/(n*r);centered(s,t)=energy(s,t)-meansquare/r;
   double g=0,d=0;
   for(unsigned a=0;a<left.size();++a){double ds=0;for(int k=0;k<r;++k){double x=u[offset+left[a]+n*k],z=u[offset+right[a]+n*k];g+=x*z;ds+=(x-z)*(x-z);}d+=std::sqrt(ds);}
   gram(s,t)=left.size()?g/left.size():NA_REAL;distance(s,t)=left.size()?d/left.size():NA_REAL;
  }
 }
 return Rcpp::List::create(Rcpp::_["energy"]=energy,Rcpp::_["centered_energy"]=centered,Rcpp::_["mean_gram"]=gram,Rcpp::_["mean_distance"]=distance);
}'
	)
	Rcpp::cppFunction(
		env = backend,
		cacheDir = file.path(root, "code", ".summary_cache"),
		code = '
Rcpp::List reset_count_probabilities(Rcpp::NumericMatrix odds,Rcpp::IntegerVector thresholds) {
 const int S=odds.nrow(),n=odds.ncol(),K=thresholds.size();
 Rcpp::NumericMatrix probability(S,K),log_odds(S,K);Rcpp::NumericVector share_odds(S),mean_counts(n+1);
 auto add=[](double a,double b){if(!R_finite(a))return b;if(!R_finite(b))return a;double top=std::max(a,b);return top+std::log1p(std::exp(-std::abs(a-b)));};
 auto softplus=[](double x){return x>0?x+std::log1p(std::exp(-x)):std::log1p(std::exp(x));};
 for(int s=0;s<S;++s){
  std::vector<double> dp(n+1,R_NegInf);dp[0]=0;double sum_p=R_NegInf,sum_q=R_NegInf;
  for(int i=0;i<n;++i){
   if(!R_finite(odds(s,i)))Rcpp::stop("Nonfinite conditional reset log odds.");
   double p=-softplus(-odds(s,i)),q=-softplus(odds(s,i));sum_p=add(sum_p,p);sum_q=add(sum_q,q);
   for(int k=i+1;k>=1;--k)dp[k]=add(dp[k]+q,dp[k-1]+p);
   dp[0]+=q;
  }
  share_odds[s]=sum_p-sum_q;
  double normalizer=R_NegInf;for(int j=0;j<=n;++j)normalizer=add(normalizer,dp[j]);
  for(int j=0;j<=n;++j)mean_counts[j]+=std::exp(dp[j]-normalizer)/S;
  for(int k=0;k<K;++k){
   if(thresholds[k]<1||thresholds[k]>n)Rcpp::stop("Invalid reset count threshold.");
   double lower=R_NegInf,upper=R_NegInf;
   for(int j=0;j<=n;++j){if(j<thresholds[k])lower=add(lower,dp[j]);else upper=add(upper,dp[j]);}
   double value=upper-lower;log_odds(s,k)=value;
   probability(s,k)=value>0?1/(1+std::exp(-value)):std::exp(value)/(1+std::exp(value));
  }
 }
 return Rcpp::List::create(Rcpp::_["probability"]=probability,Rcpp::_["log_odds"]=log_odds,Rcpp::_["share_log_odds"]=share_odds,Rcpp::_["mean_count_probabilities"]=mean_counts);
}'
	)
	backend
}

diagnose_chains = function(chains) {
	stopifnot(
		length(chains) >= 2L,
		length(unique(vapply(chains, nrow, 1L))) == 1L,
		all(
			vapply(
				chains,
				function(x) identical(colnames(x), colnames(chains[[1]])),
				TRUE
			)
		)
	)
	do.call(
		rbind,
		lapply(
			colnames(chains[[1]]),
			function(parameter) {
				x = do.call(cbind, lapply(chains, function(chain) chain[, parameter]))
				constant = all(is.finite(x)) && length(unique(as.vector(x))) == 1L
				data.frame(
					parameter = parameter,
					mean = mean(x),
					sd = sd(as.vector(x)),
					rhat = if (constant) NA_real_ else posterior::rhat(x),
					ess_bulk = if (constant) NA_real_ else posterior::ess_bulk(x),
					ess_tail = if (constant) NA_real_ else posterior::ess_tail(x),
					constant = constant,
					iterations_per_chain = nrow(x),
					chains = ncol(x)
				)
			}
		)
	)
}

summarize_chain = function(fit, domain, backend) {
	y = fit$Y
	u = fit$U_draws
	n = dim(u)[1]
	rank = dim(u)[2]
	periods = dim(u)[3]
	draws = dim(u)[4]
	actors = dimnames(y)[[1]]
	years = dimnames(y)[[3]]
	observed = vapply(
		seq_len(periods),
		function(t) rowSums(is.finite(y[, , t])) > 0,
		logical(n)
	)
	eligible = observed
	eligible[, 1] = FALSE
	eligible[, -1] = observed[, -1] & observed[, -periods]
	groups = if (domain == "controls") {
		list(all = actors, other_than_actor1 = actors[-1])
	} else if (domain == "un") {
		list(
			all_non_us = setdiff(actors, "USA"),
			focal_non_us = setdiff(unique(unlist(actor_groups, use.names = FALSE)), "USA")
		)
	} else {
		list(
			all = actors,
			focal = unique(unlist(actor_groups, use.names = FALSE)),
			nonfocal = setdiff(actors, unlist(actor_groups, use.names = FALSE))
		)
	}
	groups = lapply(groups, intersect, y = actors)
	geometry = backend$trajectory_geometry(u, y)
	continuous = cbind(
		fit$trace,
		intercept = fit$BETA[, 1],
		latent_energy = rowMeans(geometry$energy),
		centered_energy = rowMeans(geometry$centered_energy)
	)
	for (name in c("mean_gram", "mean_distance")) {
		x = geometry[[name]]
		colnames(x) = paste(name, years, sep = "_")
		continuous = cbind(continuous, x)
	}
	indicators = array(as.integer(fit$snap_draws), dim(fit$snap_draws))
	fractions = conditional = log_odds = events = event_log_odds = conditional_share_log_odds = list()
	denominators = list()
	count_distributions = list()
	actor_terminal = NULL
	for (t in 2:periods) {
		current_square = difference_square = matrix(0, n, draws)
		for (k in seq_len(rank)) {
			current = matrix(u[, k, t, ], n, draws)
			previous = matrix(u[, k, t - 1L, ], n, draws)
			current_square = current_square + current^2
			difference_square = difference_square + (current - sweep(previous, 2, fit$rho_uv, "*"))^2
		}
		odds = sweep(difference_square, 2, 2 * fit$sigma_uv^2, "/") - current_square / (2 * fit$prior$snap_kappa^2)
		odds = sweep(
			odds,
			2,
			qlogis(fit$pi_snap) + rank * log(fit$sigma_uv / fit$prior$snap_kappa),
			"+"
		)
		probability = t(plogis(odds))
		odds = t(odds)
		keep = which(eligible[, t])
		selected_odds = odds[, keep, drop = FALSE]
		colnames(selected_odds) = paste(actors[keep], years[t], sep = "_")
		log_odds[[years[t]]] = selected_odds
		if (t == periods) {
			actor_terminal = data.frame(
				actor = actors[keep],
				year = as.integer(years[t]),
				probability = colMeans(matrix(indicators[, keep, t, drop = FALSE], draws, length(keep))),
				conditional_probability = colMeans(probability[, keep, drop = FALSE])
			)
		}
		for (group in names(groups)) {
			members = intersect(match(groups[[group]], actors), keep)
			if (!length(members)) next
			key = paste(group, years[t], sep = "_")
			counts = rowSums(
				matrix(indicators[, members, t, drop = FALSE], draws, length(members))
			)
			fractions[[key]] = counts / length(members)
			conditional[[key]] = rowMeans(probability[, members, drop = FALSE])
			denominators[[key]] = data.frame(
				key = key,
				group = group,
				year = as.integer(years[t]),
				actors = length(members)
			)
			tails = backend$reset_count_probabilities(
				odds[, members, drop = FALSE],
				as.integer(ceiling(c(.1, .25, .5) * length(members)))
			)
			conditional_share_log_odds[[key]] = tails$share_log_odds
			count_distributions[[key]] = data.frame(
				key = key,
				group = group,
				year = as.integer(years[t]),
				actors = length(members),
				count = 0:length(members),
				probability = tails$mean_count_probabilities
			)
			for (k in seq_along(c(.1, .25, .5))) {
				name = paste(key, c(.1, .25, .5)[k], sep = "_")
				events[[name]] = tails$probability[, k]
				event_log_odds[[name]] = tails$log_odds[, k]
			}
		}
	}
	fractions = do.call(cbind, fractions)
	conditional = do.call(cbind, conditional)
	denominators = do.call(rbind, denominators)
	colnames(conditional) = paste0("conditional_share_", colnames(conditional))
	list(
		continuous = continuous,
		conditional_shares = conditional,
		conditional_share_log_odds = do.call(cbind, conditional_share_log_odds),
		event_log_odds = do.call(cbind, event_log_odds),
		reset_shares = fractions,
		events = do.call(cbind, events),
		actor_log_odds = do.call(cbind, log_odds),
		denominators = denominators,
		actor_terminal = actor_terminal,
		count_distributions = do.call(rbind, count_distributions),
		seed = fit$revision_build$seed,
		prior = fit$prior
	)
}

assess_diagnostics = function(diagnostics, permit_constant = FALSE) {
	if (!permit_constant && any(diagnostics$constant)) {
		return(FALSE)
	}
	x = diagnostics[!diagnostics$constant, ]
	if (!nrow(x)) {
		return(NA)
	}
	all(
		is.finite(x$rhat) & x$rhat <= 1.01 & is.finite(x$ess_bulk) & x$ess_bulk >= 400 &
			is.finite(x$ess_tail) & x$ess_tail >= 400
	)
}

summarize_features = function(features, specification) {
	stopifnot(length(features) >= 2L)
	stopifnot(
		all(
			vapply(
				features,
				function(x) {
					identical(x$denominators, features[[1]]$denominators) &&
						identical(
							x$count_distributions[, setdiff(names(x$count_distributions), "probability")],
							features[[1]]$count_distributions[, setdiff(names(features[[1]]$count_distributions), "probability")]
						)
				},
				TRUE
			)
		)
	)
	diagnostics = lapply(
		c(
			"continuous",
			"conditional_share_log_odds",
			"event_log_odds",
			"actor_log_odds"
		),
		function(family) {
			d = diagnose_chains(lapply(features, `[[`, family))
			d$family = family
			d
		}
	)
	names(diagnostics) = c(
		"continuous",
		"conditional_share_log_odds",
		"event_log_odds",
		"actor_log_odds"
	)
	status = vapply(
		names(diagnostics),
		function(family) assess_diagnostics(diagnostics[[family]]),
		TRUE
	)
	pooled = do.call(rbind, lapply(features, `[[`, "reset_shares"))
	conditional = do.call(rbind, lapply(features, `[[`, "conditional_shares"))
	pooled_events = do.call(rbind, lapply(features, `[[`, "events"))
	denominators = features[[1]]$denominators
	distribution = features[[1]]$count_distributions
	distribution$probability = Reduce(`+`, lapply(features, function(x) x$count_distributions$probability)) / length(features)
	annual = do.call(
		rbind,
		lapply(
			seq_len(nrow(denominators)),
			function(i) {
				d = denominators[i, ]
				x = pooled[, d$key]
				density = distribution[distribution$key == d$key, ]
				cdf = cumsum(density$probability)
				stopifnot(abs(tail(cdf, 1) - 1) < 1e-10)
				interval = vapply(
					c(.025, .975),
					function(q) density$count[which(cdf >= q)[1]] / d$actors,
					0.
				)
				data.frame(
					group = d$group,
					year = d$year,
					actors = d$actors,
					mean_share = mean(conditional[, paste0("conditional_share_", d$key)]),
					sampled_mean_share = mean(x),
					lower = interval[1],
					upper = interval[2],
					probability_10 = mean(pooled_events[, paste0(d$key, "_0.1")]),
					probability_25 = mean(pooled_events[, paste0(d$key, "_0.25")]),
					probability_50 = mean(pooled_events[, paste0(d$key, "_0.5")])
				)
			}
		)
	)
	annual$rank_descending = ave(
		-annual$mean_share,
		annual$group,
		FUN = function(x) rank(x, ties.method = "min")
	)
	pooled_events = do.call(rbind, lapply(features, `[[`, "events"))
	events = data.frame(
		event = colnames(pooled_events),
		probability = colMeans(pooled_events)
	)
	actors = do.call(
		rbind,
		lapply(
			features,
			function(x) {
				z = x$actor_terminal
				z$seed = x$seed
				z
			}
		)
	)
	list(
		specification = specification,
		features = features,
		diagnostics = do.call(rbind, diagnostics),
		annual = annual,
		events = events,
		actors = actors,
		count_distribution = distribution,
		status = status,
		complete = all(!is.na(status) & status),
		draws = nrow(features[[1]]$continuous),
		chains = length(features)
	)
}

append_chain_features = function(segments) {
	stopifnot(length(segments) >= 1L)
	answer = segments[[1]]
	stopifnot(
		all(
			vapply(
				segments,
				function(x) {
					identical(x$seed, answer$seed) && identical(x$prior, answer$prior) &&
						identical(x$denominators, answer$denominators)
				},
				TRUE
			)
		)
	)
	matrix_fields = c(
		"continuous",
		"conditional_shares",
		"conditional_share_log_odds",
		"event_log_odds",
		"reset_shares",
		"events",
		"actor_log_odds"
	)
	for (field in matrix_fields) {
		stopifnot(
			all(
				vapply(
					segments,
					function(x) identical(colnames(x[[field]]), colnames(answer[[field]])),
					TRUE
				)
			)
		)
		answer[[field]] = do.call(rbind, lapply(segments, `[[`, field))
	}
	draws = vapply(segments, function(x) nrow(x$continuous), 1L)
	weights = draws / sum(draws)
	for (field in c("count_distributions", "actor_terminal")) {
		values = if (field == "count_distributions") "probability" else c("probability", "conditional_probability")
		keys = setdiff(names(answer[[field]]), values)
		stopifnot(
			all(
				vapply(
					segments,
					function(x) identical(x[[field]][keys], answer[[field]][keys]),
					TRUE
				)
			)
		)
		for (value in values) answer[[field]][[value]] = Reduce(`+`, Map(function(x, w) x[[field]][[value]] * w, segments, weights))
	}
	answer$segment_draws = unlist(
		lapply(
			segments,
			function(x) if (is.null(x$segment_draws)) nrow(x$continuous) else x$segment_draws
		)
	)
	answer
}

summarize_models = function(root, domain, manifest, workers = 4L, output_domain = domain) {
	backend = summary_backend(root)
	output = file.path(root, "results", "models", domain)
	dir.create(output, recursive = TRUE, showWarnings = FALSE)
	results = parallel::mclapply(
		unique(manifest$specification),
		function(specification) {
			rows = manifest[manifest$specification == specification, ]
			fits = lapply(
				rows$fit_file,
				function(path) readRDS(if (grepl("^/", path)) path else file.path(root, path))
			)
			stopifnot(
				all(rows$complete),
				all(
					vapply(
						fits,
						function(f) isTRUE(f$complete) && identical(f$Y, fits[[1]]$Y),
						TRUE
					)
				)
			)
			features = lapply(fits, summarize_chain, domain = domain, backend = backend)
			summarize_features(features, specification)
		},
		mc.cores = workers,
		mc.preschedule = FALSE
	)
	stopifnot(all(vapply(results, is.list, TRUE)))
	write_model_summaries(root, output_domain, results)
}

write_model_summaries = function(root, domain, results) {
	output = file.path(root, "results", "models", domain)
	dir.create(output, recursive = TRUE, showWarnings = FALSE)
	collect = function(field) {
		do.call(
			rbind,
			lapply(
				results,
				function(result) {
					x = result[[field]]
					x$specification = result$specification
					x
				}
			)
		)
	}
	for (field in c("diagnostics", "annual", "events", "actors", "count_distribution")) {
		write.csv(
			collect(field),
			file.path(output, paste0(field, ".csv")),
			row.names = FALSE
		)
	}
	status = do.call(
		rbind,
		lapply(
			results,
			function(x) {
				data.frame(
					specification = x$specification,
					continuous = x$status["continuous"],
					conditional_share_log_odds = x$status["conditional_share_log_odds"],
					event_log_odds = x$status["event_log_odds"],
					actor_log_odds = x$status["actor_log_odds"],
					all_diagnostics_pass = x$complete,
					draws_per_chain = x$draws,
					chains = x$chains
				)
			}
		)
	)
	write.csv(status, file.path(output, "status.csv"), row.names = FALSE)
	features = setNames(
		lapply(results, `[[`, "features"),
		vapply(results, `[[`, "", "specification")
	)
	saveRDS(features, file.path(output, "chain_features.rds"))
	invisible(status)
}

run_models = function(root, domain, workers = 4L, cold = FALSE) {
	context = model_load(file.path(root, "code"))
	inputs = readRDS(file.path(root, "results", "models", "inputs.rds"))[[domain]]
	starts = readRDS(file.path(root, "data", "model_initial_states.rds"))[[domain]]
	plan = read.csv(file.path(root, "data", "model_chains.csv"))
	plan = plan[plan$domain == domain, ]
	priors = model_priors()
	output = file.path(
		root,
		"results",
		"models",
		paste0(domain, if (cold) "_cold" else "_replay")
	)
	dir.create(output, recursive = TRUE, showWarnings = FALSE)
	manifest = parallel::mclapply(
		seq_len(nrow(plan)),
		function(i) {
			job = plan[i, ]
			y = inputs[[job$specification]]
			initial = starts[[job$specification]]$chains[[as.character(job$seed)]]
			shifts = starts[[job$specification]]$shifts
			stopifnot(
				identical(digest::digest(y, algo = "sha256"), job$input_sha256),
				identical(digest::digest(initial, algo = "sha256"), job$initial_state_sha256),
				identical(digest::digest(shifts, algo = "sha256"), job$shifts_sha256)
			)
			prior = priors[[if (domain == "un") job$specification else domain]]
			fit = fit_model(
				y,
				paste(job$specification, job$seed, sep = "_"),
				job$seed,
				prior,
				context,
				output,
				nscan = job$sampling_sweeps,
				burn = job$warmup_sweeps,
				odens = job$thin,
				initial_state = if (cold) NULL else initial,
				frozen_shifts = if (cold) NULL else shifts,
				indicator_every = job$indicator_every,
				temporal_tau = job$temporal_tau
			)
			data.frame(
				domain = domain,
				specification = job$specification,
				seed = job$seed,
				complete = isTRUE(fit$complete),
				fit_file = fit$revision_build$fit_file,
				checkpoint_file = fit$revision_build$checkpoint_file
			)
		},
		mc.cores = workers,
		mc.preschedule = FALSE
	)
	stopifnot(all(vapply(manifest, is.data.frame, TRUE)))
	manifest = do.call(rbind, manifest)
	write.csv(manifest, file.path(output, "fit_manifest.csv"), row.names = FALSE)
	invisible(manifest)
}

assess_control_cases = function(root, output_domain = "controls") {
	output = file.path(root, "results", "models", output_domain)
	expected = read.csv(file.path(root, "results/models/control_cases.csv"))
	status = read.csv(file.path(output, "status.csv"))
	stopifnot(
		nrow(expected) == 116L,
		all(expected$checks_passed),
		setequal(unique(expected$specification), status$specification),
		!anyDuplicated(status$specification)
	)
	inputs = readRDS(file.path(root, "results/models/inputs.rds"))$controls
	truth = readRDS(file.path(root, "results/models/control_truth.rds"))
	annual = read.csv(file.path(output, "annual.csv"))
	rows = list()
	for (i in seq_len(nrow(expected))) {
		case = expected[i, ]
		y = inputs[[case$specification]]
		d = annual[annual$specification == case$specification, ]
		known = truth[[as.character(case$case_id)]]
		for (j in seq_len(nrow(d))) {
			row = d[j, ]
			period = match(as.character(row$year), dimnames(y)[[3]])
			eligible = vapply(
				seq_len(dim(y)[1]),
				function(a) any(is.finite(y[a, -a, period])) && any(is.finite(y[a, -a, period - 1L])),
				TRUE
			)
			members = which(eligible)
			if (row$group == "other_than_actor1") members = setdiff(members, 1L)
			stopifnot(length(members) == row$actors)
			for (k in seq_along(c(.1, .25, .5))) {
				threshold = c(.1, .25, .5)[k]
				probability = row[[c("probability_10", "probability_25", "probability_50")[k]]]
				true_count = sum(known$truth$reset_indicator[members, period])
				rows[[length(rows) + 1L]] = data.frame(
					case,
					year = row$year,
					period,
					group = row$group,
					actors = length(members),
					threshold,
					is_event_period = period == case$event_period,
					true_reset_count = true_count,
					true_widespread_reset = true_count >= ceiling(threshold * length(members)),
					reset_indicator_interpretable = known$truth$reset_indicator_interpretable,
					posterior_mean_share = row$mean_share,
					lower = row$lower,
					upper = row$upper,
					probability_widespread = probability,
					estimated_widespread = probability >= .95,
					all_diagnostics_pass = status$all_diagnostics_pass[match(case$specification, status$specification)]
				)
			}
		}
	}
	assessment = do.call(rbind, rows)
	stopifnot(nrow(assessment) == 116L * 11L * 2L * 3L)
	write.csv(
		assessment,
		file.path(output, "control_case_assessments.csv"),
		row.names = FALSE
	)
	write.csv(
		assessment[assessment$is_event_period, ],
		file.path(output, "control_event_assessments.csv"),
		row.names = FALSE
	)
	events = assessment[assessment$is_event_period & assessment$group == "all", ]
	recovery = lapply(
		unique(expected$case),
		function(case) {
			selected = events[events$case == case, ]
			row = data.frame(
				case = case,
				distinct_arrays = length(unique(selected$specification))
			)
			for (threshold in c(.1, .25, .5)) {
				x = unique(
					selected[selected$threshold == threshold, c("specification", "estimated_widespread")]
				)
				stopifnot(!anyDuplicated(x$specification), nrow(x) == row$distinct_arrays)
				row[[paste0("detected_", threshold * 100)]] = sum(x$estimated_widespread)
			}
			row
		}
	)
	write.csv(
		do.call(rbind, recovery),
		file.path(output, "control_recovery.csv"),
		row.names = FALSE
	)
	invisible(assessment)
}

archive_model_metadata = function(root, domain, execution_manifest) {
	stopifnot(domain %in% c("un", "trade", "controls"))
	output = file.path(root, "results", "models", domain)
	features = readRDS(file.path(output, "chain_features.rds"))
	if (is.character(execution_manifest)) execution_manifest = read.csv(execution_manifest)
	stopifnot(
		all(
			c("specification", "seed", "segment", "fit_file", "checkpoint_file") %in% names(execution_manifest)
		)
	)
	executed = list()
	priors = model_priors()
	plan = read.csv(file.path(root, "data", "model_chains.csv"))
	plan = plan[plan$domain == domain, ]
	stopifnot(setequal(names(features), unique(plan$specification)))
	manifest = do.call(
		rbind,
		lapply(
			seq_len(nrow(plan)),
			function(i) {
				job = plan[i, ]
				chains = features[[job$specification]]
				at = which(vapply(chains, function(x) x$seed == job$seed, TRUE))
				stopifnot(length(at) == 1L)
				stored = nrow(chains[[at]]$continuous)
				runs = execution_manifest[execution_manifest$specification == job$specification & execution_manifest$seed == job$seed, ]
				runs = runs[order(runs$segment), ]
				stopifnot(
					nrow(runs) >= 1L,
					identical(as.integer(runs$segment), seq_len(nrow(runs)))
				)
				warmup = sampling = offset = 0L
				previous_checkpoint = NULL
				for (j in seq_len(nrow(runs))) {
					run = runs[j, ]
					fit = readRDS(run$fit_file)
					config = fit$revision_build$config
					stopifnot(
						isTRUE(fit$complete),
						config$seed == job$seed,
						config$odens == job$thin,
						identical(config$input_hash, job$input_sha256),
						identical(digest::digest(fit$Y, algo = "sha256"), job$input_sha256),
						identical(
							config$prior,
							priors[[if (domain == "un") job$specification else domain]]
						),
						identical(
							digest::digest(config$frozen_shifts, algo = "sha256"),
							job$shifts_sha256
						),
						nrow(fit$trace) * config$odens == config$nscan,
						tail(fit$monitor$iteration, 1L) == config$nscan + config$burn,
						identical(unname(tools::md5sum(names(config$code))), unname(config$code))
					)
					stopifnot(
						(if (is.null(config$indicator_every)) 0L else config$indicator_every) == job$indicator_every
					)
					stopifnot(
						(if (is.null(config$temporal_tau)) 0 else config$temporal_tau) == job$temporal_tau
					)
					if (job$temporal_tau > 0) stopifnot(config$temporal_every == 5L, config$temporal_moves == 5L)
					if (j == 1L) {
						stopifnot(
							identical(
								digest::digest(config$initial_state, algo = "sha256"),
								job$initial_state_sha256
							)
						)
					} else {
						stopifnot(
							config$burn == 0L,
							identical(config$initial_state$parent_checkpoint_sha256, previous_checkpoint)
						)
						parent = readRDS(runs$checkpoint_file[j - 1L])
						fields = c(
							"u",
							"a",
							"beta",
							"s2",
							"va",
							"rho",
							"sigma",
							"pi",
							"delta",
							"rng",
							"initial_joint_done"
						)
						stopifnot(
							identical(config$initial_state[fields], parent[fields]),
							config$initial_state$bridge_phase == (parent$bridge_phase + parent$iteration) %% 5L
						)
					}
					indices = offset + seq_len(nrow(fit$trace))
					stopifnot(
						isTRUE(
							all.equal(
								chains[[at]]$continuous[indices, colnames(fit$trace), drop = FALSE],
								fit$trace,
								check.attributes = FALSE,
								tolerance = 0
							)
						),
						identical(
							as.numeric(chains[[at]]$continuous[indices, "intercept"]),
							as.numeric(fit$BETA[, 1])
						)
					)
					previous_checkpoint = digest::digest(file = run$checkpoint_file, algo = "sha256")
					assign(
						"executed",
						append(
							executed,
							list(
								list(
									specification = job$specification,
									seed = job$seed,
									segment = j,
									warmup_sweeps = config$burn,
									sampling_sweeps = config$nscan,
									thin = config$odens,
									stored_draws = nrow(fit$trace),
									fit_sha256 = digest::digest(file = run$fit_file, algo = "sha256"),
									checkpoint_sha256 = previous_checkpoint,
									initial_state_sha256 = digest::digest(config$initial_state, algo = "sha256"),
									input_sha256 = config$input_hash,
									source_md5 = as.list(config$code)
								)
							)
						),
						envir = parent.env(environment())
					)
					warmup = warmup + config$burn
					sampling = sampling + config$nscan
					offset = offset + nrow(fit$trace)
				}
				stopifnot(
					stored == offset,
					stored * job$thin == sampling,
					sampling == job$sampling_sweeps,
					warmup == job$warmup_sweeps
				)
				data.frame(
					domain = domain,
					specification = job$specification,
					seed = job$seed,
					complete = TRUE,
					warmup_sweeps = warmup,
					sampling_sweeps = sampling,
					thin = job$thin,
					indicator_every = job$indicator_every,
					temporal_tau = job$temporal_tau,
					stored_draws = stored,
					completed_sweeps = warmup + sampling,
					input_sha256 = job$input_sha256,
					initial_state_sha256 = job$initial_state_sha256,
					shifts_sha256 = job$shifts_sha256
				)
			}
		)
	)
	stopifnot(length(executed) == nrow(execution_manifest))
	write.csv(manifest, file.path(output, "fit_manifest.csv"), row.names = FALSE)
	sources = c(
		"code/00_sampler.R",
		"code/sampler.cpp",
		"code/10_models.R",
		"code/01_observed_votes.R",
		"code/vendor/lame_1.3.5_americium_full_fix.tar.gz",
		"data/model_chains.csv",
		"data/model_initial_states.rds",
		"results/models/input_manifest.csv",
		"results/models/control_cases.csv"
	)
	artifacts = file.path(
		"results",
		"models",
		domain,
		c(
			"diagnostics.csv",
			"status.csv",
			"annual.csv",
			"events.csv",
			"actors.csv",
			"count_distribution.csv",
			"fit_manifest.csv",
			"chain_features.rds"
		)
	)
	if (domain == "controls") {
		artifacts = c(
			artifacts,
			file.path(
				"results/models/controls",
				c(
					"control_case_assessments.csv",
					"control_event_assessments.csv",
					"control_recovery.csv"
				)
			)
		)
	}
	hashes = function(paths) {
		setNames(
			lapply(
				paths,
				function(path) digest::digest(file = file.path(root, path), algo = "sha256")
			),
			paths
		)
	}
	families = c(
		"continuous",
		"conditional_share_log_odds",
		"event_log_odds",
		"actor_log_odds"
	)
	parameters = setNames(
		lapply(
			features,
			function(chains) {
				setNames(
					lapply(families, function(family) I(colnames(chains[[1]][[family]]))),
					families
				)
			}
		),
		names(features)
	)
	jsonlite::write_json(
		list(
			schema_version = 1L,
			sources = hashes(sources),
			artifacts = hashes(artifacts),
			diagnostic_parameters = parameters,
			executed_segments = executed
		),
		file.path(output, "summary_manifest.json"),
		pretty = TRUE,
		auto_unbox = TRUE
	)
	invisible(manifest)
}

recheck_models = function(root, domain) {
	output = file.path(root, "results", "models", domain)
	features = readRDS(file.path(output, "chain_features.rds"))
	results = lapply(
		names(features),
		function(specification) summarize_features(features[[specification]], specification)
	)
	keys = list(
		diagnostics = c("specification", "family", "parameter"),
		annual = c("specification", "group", "year"),
		events = c("specification", "event"),
		actors = c("specification", "seed", "actor", "year"),
		count_distribution = c("specification", "key", "count")
	)
	for (field in names(keys)) {
		actual = do.call(
			rbind,
			lapply(
				results,
				function(x) {
					z = x[[field]]
					z$specification = x$specification
					z
				}
			)
		)
		saved = read.csv(file.path(output, paste0(field, ".csv")))
		actual = actual[do.call(order, actual[keys[[field]]]), ]
		saved = saved[do.call(order, saved[keys[[field]]]), ]
		rownames(actual) = rownames(saved) = NULL
		stopifnot(
			isTRUE(all.equal(actual, saved, check.attributes = FALSE, tolerance = 1e-12))
		)
	}
	status = read.csv(file.path(output, "status.csv"))
	for (result in results) {
		row = status[status$specification == result$specification, ]
		stopifnot(
			nrow(row) == 1L,
			row$all_diagnostics_pass == result$complete,
			row$draws_per_chain == result$draws,
			row$chains == result$chains,
			identical(as.logical(row[names(result$status)]), unname(result$status))
		)
	}
	metadata_file = file.path(output, "summary_manifest.json")
	if (file.exists(metadata_file)) {
		metadata = jsonlite::read_json(metadata_file, simplifyVector = FALSE)
		for (field in c("sources", "artifacts")) {
			for (path in names(metadata[[field]])) {
				stopifnot(
					identical(
						digest::digest(file = file.path(root, path), algo = "sha256"),
						metadata[[field]][[path]]
					)
				)
			}
		}
		jsonlite::write_json(
			list(
				passed = TRUE,
				sources = metadata$sources,
				artifacts = metadata$artifacts
			),
			file.path(output, "recheck.json"),
			pretty = TRUE,
			auto_unbox = TRUE
		)
	}
	invisible(TRUE)
}

if (sys.nframe() == 0L) {
	args = commandArgs(trailingOnly = TRUE)
	argument = sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
	root = dirname(dirname(normalizePath(argument)))
	source(file.path(root, "code", "01_observed_votes.R"))
	source(file.path(root, "code", "00_sampler.R"))
	mode = if (length(args)) args[1] else "prepare"
	domain = if (length(args) >= 2L) args[2] else "un"
	workers = if (length(args) >= 3L) as.integer(args[3]) else 4L
	stopifnot(
		mode %in% c("prepare", "fit", "cold", "summarize", "recheck"),
		domain %in% c("un", "trade", "controls")
	)
	if (mode == "prepare") prepare_models(root)
	if (mode %in% c("fit", "cold")) run_models(root, domain, workers, cold = mode == "cold")
	if (mode == "summarize") {
		manifest = read.csv(
			file.path(
				root,
				"results",
				"models",
				paste0(domain, "_replay"),
				"fit_manifest.csv"
			)
		)
		summarize_models(
			root,
			domain,
			manifest,
			workers,
			output_domain = paste0(domain, "_replay")
		)
		if (domain == "controls") assess_control_cases(root, output_domain = "controls_replay")
	}
	if (mode == "recheck") recheck_models(root, domain)
}
