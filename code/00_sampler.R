# load and fit the latent network sampler

model_load = function(directory, library = NULL) {
	directory = normalizePath(directory)
	archive = file.path(directory, "vendor", "lame_1.3.5_americium_full_fix.tar.gz")
	if (!requireNamespace("digest", quietly = TRUE) || !requireNamespace("Rcpp", quietly = TRUE)) {
		stop("Install the declared R dependencies before loading the model.")
	}
	if (!file.exists(archive) || !identical(
		digest::digest(file = archive, algo = "sha256"),
		"016ea875e99d0fbe6829593ce993fb456f64a387302d02424dd9104d470f35f8"
	)) {
		stop("The validated statistical source archive is absent or has changed.")
	}
	if (is.null(library)) library = Sys.getenv("AMERICIUM_MODEL_LIB", unset = file.path(directory, ".library"))
	dir.create(library, recursive = TRUE, showWarnings = FALSE)
	library = normalizePath(library)
	if (!file.exists(file.path(library, "lame", "DESCRIPTION"))) {
		utils::install.packages(archive, repos = NULL, type = "source", lib = library)
	}
	if (!file.exists(file.path(library, "lame", "DESCRIPTION"))) {
		stop(
			"The archived model dependency did not install. Install its declared Imports and LinkingTo dependencies, then retry; another installed lame version cannot substitute for it."
		)
	}
	.libPaths(c(library, .libPaths()))
	if ("lame" %in% loadedNamespaces() &&
		!identical(
			normalizePath(getNamespaceInfo(asNamespace("lame"), "path")),
			normalizePath(file.path(library, "lame"))
		)) {
		stop("Restart R before selecting a different model library.")
	}
	if (!requireNamespace("lame", quietly = TRUE)) stop("The archived model dependency failed to load.")
	if (!identical(
		normalizePath(getNamespaceInfo(asNamespace("lame"), "path")),
		normalizePath(file.path(library, "lame"))
	)) {
		stop(
			"The loaded model dependency is not the archived installation. Restart R and retry."
		)
	}
	package_files = sort(
		list.files(file.path(library, "lame"), recursive = TRUE, full.names = TRUE)
	)
	package_hashes = tools::md5sum(package_files)
	names(package_hashes) = substring(names(package_hashes), nchar(library) + 2)
	backend = new.env(parent = baseenv())
	compiled = Rcpp::sourceCpp(
		file.path(directory, "sampler.cpp"),
		env = backend,
		cacheDir = file.path(directory, ".compiled"),
		showOutput = FALSE,
		verbose = FALSE
	)
	backend$artifact = list(
		source_sha256 = digest::digest(file = file.path(directory, "sampler.cpp"), algo = "sha256"),
		rcpp_version = as.character(utils::packageVersion("Rcpp")),
		armadillo_version = as.character(utils::packageVersion("RcppArmadillo"))
	)
	list(
		directory = directory,
		archive = archive,
		library = library,
		backend = backend,
		fingerprint = digest::digest(package_hashes, algo = "sha256")
	)
}

revision_static_level_draw = function(
	statistics,
	s2,
	va,
	beta_precision,
	return_moments = FALSE,
	sample = TRUE
) {
	n = nrow(statistics$row_sums)
	periods = ncol(statistics$row_sums)
	if (n < 2 || periods < 1 || length(statistics$dyad_sums) != periods ||
		any(
			!is.finite(
				c(statistics$row_sums, statistics$dyad_sums, s2, va, beta_precision)
			)
		) ||
		any(c(s2, va, beta_precision) <= 0)) {
		stop("Invalid exact static-level statistics or prior constants.")
	}
	dyads = n * (n - 1) / 2
	rows = rowSums(statistics$row_sums)
	q_contrast = periods * (n - 2) / s2 + 1 / va
	mean_contrast = (rows - mean(rows)) / (s2 * q_contrast)
	q_common = diag(c(beta_precision, n / va)) +
		periods * dyads / s2 * matrix(c(1, 2, 2, 4), 2, 2)
	root = chol(q_common)
	linear = sum(statistics$dyad_sums) / s2 * c(1, 2)
	mean_common = backsolve(root, forwardsolve(t(root), linear))
	if (return_moments) {
		mapping = rbind(c(1, 0), cbind(rep(0, n), rep(1, n)))
		covariance = mapping %*% chol2inv(root) %*% t(mapping)
		covariance[-1, -1] = covariance[-1, -1] +
			(diag(n) - matrix(1 / n, n, n)) / q_contrast
		return(
			list(
				mean = c(mean_common[1], mean_contrast + mean_common[2]),
				covariance = covariance
			)
		)
	}
	contrast = mean_contrast
	common = mean_common
	if (sample) {
		contrast = contrast + rnorm(n, sd = 1 / sqrt(q_contrast))
		common = common + backsolve(root, rnorm(2))
	}
	contrast = contrast - mean(contrast)
	list(beta = unname(common[1]), a = contrast + common[2])
}

revision_rho_conditional_stats = function(
	u_positions,
	v_positions,
	sigma_uv,
	symmetric,
	delta_u = NULL,
	delta_v = NULL,
	prior_mean = 0.9,
	prior_sd = 0.1
) {
	if (length(dim(u_positions)) != 3 || any(!is.finite(u_positions)) || !is.finite(sigma_uv) ||
		sigma_uv <= 0 || !is.finite(prior_sd) || prior_sd <= 0) {
		stop("Invalid rho state.")
	}
	periods = dim(u_positions)[3]
	if (is.null(delta_u)) {
		delta_u = matrix(0, dim(u_positions)[1], periods)
	}
	if (!symmetric && is.null(delta_v)) {
		delta_v = matrix(0, dim(v_positions)[1], periods)
	}
	lag_square = lag_product = initial_sum_squares = n1 = 0
	for (side in if (symmetric) "u" else c("u", "v")) {
		x = if (side == "u") u_positions else v_positions
		delta = if (side == "u") delta_u else delta_v
		if (length(dim(x)) != 3 || dim(x)[3] != periods || any(!is.finite(x)) ||
			!identical(dim(delta), c(dim(x)[1], periods)) || any(!delta %in% c(0, 1))) {
			stop("Latent positions and reset indicators do not match.")
		}
		initial_sum_squares = initial_sum_squares + sum(x[, , 1]^2)
		n1 = n1 + dim(x)[1] * dim(x)[2]
		if (periods > 1) {
			for (t in 2:periods) {
				keep = delta[, t] == 0
				lag = x[keep, , t - 1, drop = FALSE]
				current = x[keep, , t, drop = FALSE]
				lag_square = lag_square + sum(lag^2)
				lag_product = lag_product + sum(lag * current)
			}
		}
	}
	variance = 1 / (lag_square / sigma_uv^2 + 1 / prior_sd^2)
	mean = variance * (lag_product / sigma_uv^2 + prior_mean / prior_sd^2)
	list(
		mean_post = mean,
		var_post = variance,
		initial_sum_squares = initial_sum_squares,
		n1 = n1,
		sigma_uv = sigma_uv
	)
}

revision_stationary_rho_slice = function(
	rho_current,
	mean_post,
	var_post,
	initial_sum_squares,
	n1,
	sigma_uv,
	bound = 0.99
) {
	if (any(
		!is.finite(
			c(rho_current, mean_post, var_post, initial_sum_squares, n1, sigma_uv, bound)
		)
	) ||
		var_post <= 0 || sigma_uv <= 0 || initial_sum_squares < 0 || n1 < 0 || bound <= 0 || bound >= 1 || abs(rho_current) > bound) {
		stop("Invalid stationary rho conditional parameters.")
	}
	quadratic = 1 / var_post - initial_sum_squares / sigma_uv^2
	linear = mean_post / var_post
	current_log_stationary = log1p(-rho_current^2)
	log_height = log(runif(1))
	lower = -bound
	upper = bound
	for (attempt in seq_len(100000L)) {
		proposal = runif(1, lower, upper)
		log_ratio = (proposal - rho_current) * (linear - 0.5 * quadratic * (proposal + rho_current)) +
			0.5 * n1 * (log1p(-proposal^2) - current_log_stationary)
		if (!is.finite(log_ratio)) {
			stop("Nonfinite stationary rho slice density.")
		}
		if (log_ratio >= log_height) {
			return(proposal)
		}
		if (proposal < rho_current) lower = proposal else upper = proposal
	}
	stop("Stationary rho slice exhausted its numeric iteration limit.")
}

revision_rho_uv_slice = function(
	u_positions,
	v_positions,
	sigma_uv,
	rho_current,
	symmetric,
	delta_u = NULL,
	delta_v = NULL,
	prior_mean = 0.9,
	prior_sd = 0.1,
	bound = 0.99
) {
	statistics = revision_rho_conditional_stats(
		u_positions,
		v_positions,
		sigma_uv,
		symmetric,
		delta_u,
		delta_v,
		prior_mean,
		prior_sd
	)
	do.call(
		revision_stationary_rho_slice,
		c(list(rho_current = rho_current, bound = bound), statistics)
	)
}

revision_sigma_uv_draw = function(u, rho, delta, prior_shape, prior_scale) {
	if (length(dim(u)) != 3L || dim(u)[3] < 2L || !is.matrix(delta) ||
		!identical(dim(delta), as.integer(dim(u)[c(1, 3)])) || any(!is.finite(u)) ||
		anyNA(delta) || any(!delta %in% c(0, 1)) || length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1 ||
		any(!is.finite(c(prior_shape, prior_scale))) || prior_shape <= 0 || prior_scale <= 0) {
		stop("Invalid drift scale conditional inputs.")
	}
	rank = dim(u)[2]
	sum_squares = (1 - rho^2) * sum(u[, , 1]^2)
	count = dim(u)[1] * rank
	for (t in 2:dim(u)[3]) {
		keep = delta[, t] == 0
		if (any(keep)) {
			innovation = u[keep, , t, drop = FALSE] - rho * u[keep, , t - 1L, drop = FALSE]
			sum_squares = sum_squares + sum(innovation^2)
			count = count + sum(keep) * rank
		}
	}
	sqrt(
		1 / rgamma(
			1,
			shape = count / 2 + prior_shape,
			rate = sum_squares / 2 + prior_scale
		)
	)
}

revision_missing_suffix_plan = function(y) {
	if (length(dim(y)) != 3L || dim(y)[1] < 2L || dim(y)[2] != dim(y)[1] || dim(y)[3] < 1L) {
		stop("Missing suffixes require a square actor-by-actor-by-period array.")
	}
	n = dim(y)[1]
	periods = dim(y)[3]
	observed = matrix(FALSE, n, periods)
	for (t in seq_len(periods)) {
		finite = is.finite(y[, , t])
		diag(finite) = FALSE
		if (!identical(finite, t(finite))) {
			stop("The missing-suffix plan requires symmetric dyad observation masks.")
		}
		observed[, t] = rowSums(finite) > 0
	}
	last = apply(
		observed,
		1,
		function(x) {
			if (any(x)) max(which(x)) else 0L
		}
	)
	list(
		n = n,
		periods = periods,
		last_observed = as.integer(last),
		actor_observed = observed
	)
}

revision_missing_suffix_draw = function(u, delta, plan, rho, sigma, kappa, log_odds) {
	if (length(dim(u)) != 3L || dim(u)[2] < 1L || !is.matrix(delta) ||
		!identical(dim(delta), as.integer(dim(u)[c(1, 3)])) ||
		!identical(as.integer(c(plan$n, plan$periods)), as.integer(dim(u)[c(1, 3)])) ||
		length(plan$last_observed) != plan$n ||
		any(!is.finite(u)) || anyNA(delta) || any(!delta %in% c(0, 1)) || any(delta[, 1] != 0)) {
		stop(
			"Latent positions, indicators, and missing-suffix plan must have matching finite dimensions."
		)
	}
	if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1 ||
		length(sigma) != 1L || !is.finite(sigma) || sigma <= 0 ||
		length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
		stop(
			"Missing trajectories require positive scales and an AR coefficient inside (-1,1)."
		)
	}
	last = plan$last_observed
	if (anyNA(last) || any(last != as.integer(last)) || any(last < 0 | last > plan$periods) ||
		!identical(dim(plan$actor_observed), as.integer(c(plan$n, plan$periods))) ||
		anyNA(plan$actor_observed)) {
		stop("The missing-suffix plan contains invalid observation information.")
	}
	actual_last = apply(
		plan$actor_observed,
		1,
		function(x) {
			if (any(x)) max(which(x)) else 0L
		}
	)
	if (!identical(as.integer(actual_last), as.integer(last))) {
		stop("The suffix boundary does not match the observation plan.")
	}
	if (length(log_odds) == 1L) {
		log_odds = rep(log_odds, plan$periods)
	}
	if (length(log_odds) != plan$periods || anyNA(log_odds[-1])) {
		stop(
			"Supply one shared log odds or one value per period, excluding the unused initial slot."
		)
	}
	softplus = function(x) {
		pmax(x, 0) + log1p(exp(-abs(x)))
	}
	log_pi = -softplus(-log_odds)
	updated_u = u
	updated_delta = delta
	rank = dim(u)[2]
	updated_actors = 0L
	updated_transitions = 0L
	updated_initial_paths = 0L
	for (i in seq_len(plan$n)) {
		if (last[i] == plan$periods) {
			next
		}
		updated_actors = updated_actors + 1L
		if (last[i] == 0L) {
			updated_u[i, , 1] = stats::rnorm(rank, sd = sigma / sqrt(1 - rho^2))
			updated_initial_paths = updated_initial_paths + 1L
		}
		first = max(2L, last[i] + 1L)
		if (first > plan$periods) {
			next
		}
		for (t in seq.int(first, plan$periods)) {
			reset = log(stats::runif(1)) < log_pi[t]
			updated_delta[i, t] = as.integer(reset)
			updated_u[i, , t] = if (reset) {
				stats::rnorm(rank, sd = kappa)
			} else {
				rho * updated_u[i, , t - 1L] + stats::rnorm(rank, sd = sigma)
			}
			updated_transitions = updated_transitions + 1L
		}
	}
	list(
		u = updated_u,
		delta = updated_delta,
		updated_actors = updated_actors,
		updated_transitions = updated_transitions,
		updated_initial_paths = updated_initial_paths
	)
}

revision_missing_bridge_blocks = function(observed) {
	rows = list()
	for (actor in seq_len(nrow(observed))) {
		seen = which(observed[actor, ])
		if (!length(seen)) {
			next
		}
		if (seen[1] > 1L) {
			rows[[length(rows) + 1L]] = data.frame(actor = actor, left = 0L, right = seen[1], kind = "leading")
		}
		if (length(seen) > 1L) {
			for (j in 2:length(seen)) {
				if (seen[j] > seen[j - 1L] + 1L) {
					rows[[length(rows) + 1L]] = data.frame(
						actor = actor,
						left = seen[j - 1L],
						right = seen[j],
						kind = "internal"
					)
				}
			}
		}
	}
	if (!length(rows)) {
		return(
			data.frame(
				actor = integer(),
				left = integer(),
				right = integer(),
				kind = character()
			)
		)
	}
	do.call(rbind, rows)
}

revision_missing_bridge_plan = function(y) {
	if (length(dim(y)) != 3L || dim(y)[1] < 2L || dim(y)[1] != dim(y)[2] || dim(y)[3] < 1L) {
		stop("Missing bridges require a square actor-by-actor-by-period array.")
	}
	n = dim(y)[1]
	periods = dim(y)[3]
	observed = matrix(FALSE, n, periods)
	for (t in seq_len(periods)) {
		finite = is.finite(y[, , t])
		diag(finite) = FALSE
		if (!identical(finite, t(finite))) {
			stop("Bridge plans require symmetric dyad observation masks.")
		}
		observed[, t] = rowSums(finite) > 0
	}
	list(
		n = n,
		periods = periods,
		actor_observed = observed,
		blocks = revision_missing_bridge_blocks(observed)
	)
}

revision_missing_bridge_logs = function(log_odds, periods) {
	if (length(log_odds) == 1L) {
		log_odds = rep(log_odds, periods)
	}
	if (length(log_odds) < periods || anyNA(log_odds[seq.int(2L, periods)])) {
		stop(
			"Bridge rates require one shared log odds or a complete eligible period vector."
		)
	}
	softplus = function(x) {
		pmax(x, 0) + log1p(exp(-abs(x)))
	}
	list(log_pi = -softplus(-log_odds), log_drift = -softplus(log_odds))
}

revision_missing_bridge_weights = function(
	left_position,
	right_position,
	left,
	right,
	rho,
	sigma,
	kappa,
	log_odds
) {
	if (length(left) != 1L || length(right) != 1L || any(!is.finite(c(left, right))) ||
		left != as.integer(left) || right != as.integer(right) || left < 0 || right < 2 ||
		(left > 0 && right <= left + 1L) || !length(right_position) || any(!is.finite(right_position)) ||
		(left > 0 && (length(left_position) != length(right_position) || any(!is.finite(left_position)))) ||
		length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1 ||
		length(sigma) != 1L || !is.finite(sigma) || sigma <= 0 ||
		length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
		stop(
			"Bridge endpoints, dimensions and hyperparameters must define a valid missing interval."
		)
	}
	logs = revision_missing_bridge_logs(log_odds, right)
	rank = length(right_position)
	first_transition = if (left == 0L) 2L else left + 1L
	candidates = c(0L, seq.int(first_transition, right))
	means = matrix(0, length(candidates), rank)
	variances = numeric(length(candidates))
	log_priors = numeric(length(candidates))
	log_likelihoods = numeric(length(candidates))
	for (k in seq_along(candidates)) {
		last_reset = candidates[k]
		if (last_reset == 0L) {
			if (left == 0L) {
				variance = sigma^2 / (1 - rho^2)
			} else {
				variance = 0
				for (t in seq.int(left + 1L, right)) {
					variance = rho^2 * variance + sigma^2
				}
				means[k, ] = rho^(right - left) * left_position
			}
			log_prior = sum(logs$log_drift[seq.int(first_transition, right)])
		} else {
			variance = kappa^2
			log_prior = logs$log_pi[last_reset]
			if (last_reset < right) {
				for (t in seq.int(last_reset + 1L, right)) {
					variance = rho^2 * variance + sigma^2
				}
				log_prior = log_prior + sum(logs$log_drift[seq.int(last_reset + 1L, right)])
			}
		}
		if (!is.finite(variance) || variance <= 0) {
			stop("A bridge endpoint variance exceeded numerical support.")
		}
		variances[k] = variance
		log_priors[k] = log_prior
		log_likelihoods[k] = sum(stats::dnorm(right_position, means[k, ], sqrt(variance), log = TRUE))
	}
	reference = which.max(variances)
	variance_reference = variances[reference]
	h_scaled = (right_position - means[reference, ]) / sqrt(variance_reference)
	h_squared = sum(h_scaled^2)
	if (!is.finite(h_squared)) {
		stop("The scaled bridge endpoint norm exceeded numerical support.")
	}
	relative_likelihood = numeric(length(candidates))
	for (k in seq_along(candidates)) {
		log_variance_ratio = log(variances[k]) - log(variance_reference)
		d_scaled = (means[k, ] - means[reference, ]) / sqrt(variances[k])
		variance_term = .5 * h_squared * ((variances[k] - variance_reference) / variances[k])
		linear_term = sum(h_scaled * d_scaled) * exp(-.5 * log_variance_ratio)
		relative_likelihood[k] = -.5 * rank * log_variance_ratio + variance_term + linear_term - .5 * sum(d_scaled^2)
	}
	if (any(!is.finite(relative_likelihood))) {
		stop("A Gaussian bridge density ratio exceeded numerical support.")
	}
	likelihood_reference = log_likelihoods[reference]
	log_weights = log_priors + relative_likelihood
	maximum = max(log_weights)
	if (!is.finite(maximum)) {
		stop("Every bridge component has numerically zero conditional support.")
	}
	log_normalizer = maximum + log(sum(exp(log_weights - maximum)))
	list(
		last_reset = candidates,
		log_weight = log_weights,
		log_probability = log_weights - log_normalizer,
		endpoint_variance = variances,
		endpoint_mean = means,
		log_prior = log_priors,
		endpoint_log_density = log_likelihoods,
		log_likelihood_reference = likelihood_reference,
		relative_log_likelihood = relative_likelihood
	)
}

revision_missing_bridge_draw = function(u, delta, plan, rho, sigma, kappa, log_odds) {
	if (length(dim(u)) != 3L || dim(u)[2] < 1L || !is.matrix(delta) ||
		!identical(dim(delta), as.integer(dim(u)[c(1, 3)])) ||
		!identical(as.integer(c(plan$n, plan$periods)), as.integer(dim(u)[c(1, 3)])) ||
		!is.logical(plan$actor_observed) || anyNA(plan$actor_observed) ||
		!identical(dim(plan$actor_observed), as.integer(c(plan$n, plan$periods))) ||
		any(!is.finite(u)) || anyNA(delta) || any(!delta %in% c(0, 1)) || any(delta[, 1] != 0)) {
		stop(
			"Bridge states and observation plans require matching finite dimensions."
		)
	}
	if (!identical(plan$blocks, revision_missing_bridge_blocks(plan$actor_observed))) {
		stop("The bridge intervals do not match the observation plan.")
	}
	if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1 ||
		length(sigma) != 1L || !is.finite(sigma) || sigma <= 0 ||
		length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
		stop(
			"Missing bridges require positive scales and an AR coefficient inside (-1,1)."
		)
	}
	updated_u = u
	updated_delta = delta
	components = plan$blocks
	components$last_reset = integer(nrow(components))
	if (!nrow(plan$blocks)) {
		return(
			list(
				u = updated_u,
				delta = updated_delta,
				updated_blocks = 0L,
				updated_positions = 0L,
				updated_indicators = 0L,
				components = components
			)
		)
	}
	logs = revision_missing_bridge_logs(log_odds, plan$periods)
	rank = dim(u)[2]
	positions = indicators = 0L
	for (block in seq_len(nrow(plan$blocks))) {
		actor = plan$blocks$actor[block]
		left = plan$blocks$left[block]
		right = plan$blocks$right[block]
		weights = revision_missing_bridge_weights(
			if (left == 0L) NULL else updated_u[actor, , left],
			updated_u[actor, , right],
			left,
			right,
			rho,
			sigma,
			kappa,
			log_odds
		)
		# exponential races avoid probability floors and cumulative tail fallbacks
		relative_log_weight = weights$log_weight - max(weights$log_weight)
		selected = which.min(log(stats::rexp(length(weights$last_reset))) - relative_log_weight)
		last_reset = weights$last_reset[selected]
		components$last_reset[block] = last_reset
		first_transition = if (left == 0L) 2L else left + 1L
		if (last_reset > 0L) {
			if (left == 0L) {
				updated_u[actor, , 1] = stats::rnorm(rank, sd = sigma / sqrt(1 - rho^2))
			}
			if (last_reset > first_transition) {
				for (t in seq.int(first_transition, last_reset - 1L)) {
					reset = log(stats::runif(1)) < logs$log_pi[t]
					updated_delta[actor, t] = as.integer(reset)
					updated_u[actor, , t] = if (reset) {
						stats::rnorm(rank, sd = kappa)
					} else {
						rho * updated_u[actor, , t - 1L] + stats::rnorm(rank, sd = sigma)
					}
				}
			}
			updated_delta[actor, last_reset] = 1L
			if (last_reset < right) {
				updated_delta[actor, seq.int(last_reset + 1L, right)] = 0L
			}
			start = last_reset
			initial_mean = rep(0, rank)
			initial_variance = kappa^2
		} else {
			updated_delta[actor, seq.int(first_transition, right)] = 0L
			start = if (left == 0L) 1L else left
			initial_mean = if (left == 0L) rep(0, rank) else updated_u[actor, , left]
			initial_variance = if (left == 0L) sigma^2 / (1 - rho^2) else 0
		}
		if (start < right) {
			length_path = right - start + 1L
			means = matrix(0, rank, length_path)
			variances = numeric(length_path)
			means[, 1] = initial_mean
			variances[1] = initial_variance
			for (j in 2:length_path) {
				means[, j] = rho * means[, j - 1L]
				variances[j] = rho^2 * variances[j - 1L] + sigma^2
			}
			lowest = if (last_reset == 0L && left > 0L) left + 1L else start
			if (lowest < right) {
				for (t in seq.int(right - 1L, lowest)) {
					j = t - start + 1L
					gain = rho * variances[j] / variances[j + 1L]
					mean = means[, j] + gain * (updated_u[actor, , t + 1L] - means[, j + 1L])
					variance = (variances[j] / variances[j + 1L]) * sigma^2
					if (!is.finite(variance) || variance <= 0) {
						stop("A bridge conditional variance exceeded numerical support.")
					}
					updated_u[actor, , t] = mean + stats::rnorm(rank, sd = sqrt(variance))
				}
			}
		}
		positions = positions + if (left == 0L) right - 1L else right - left - 1L
		indicators = indicators + right - first_transition + 1L
	}
	if (any(!is.finite(updated_u))) {
		stop("A missing bridge draw exceeded numerical support.")
	}
	list(
		u = updated_u,
		delta = updated_delta,
		updated_blocks = nrow(plan$blocks),
		updated_positions = positions,
		updated_indicators = indicators,
		components = components
	)
}

# jointly sample the declared symmetric gaussian network model
fit_model = function(
	y,
	label,
	seed,
	prior,
	context,
	output_dir,
	nscan = 6000,
	burn = 3000,
	odens = 5,
	rank = 2,
	latent_kind = "snap",
	checkpoint_every = 500,
	max_sweeps = Inf,
	subdir = "fits",
	initial_state = NULL,
	frozen_shifts = NULL,
	indicator_every = 0L,
	temporal_tau = 0
) {
	if (rank != 2L) stop("The declared model has latent rank two.")
	latent_kind = match.arg(latent_kind, c("snap", "ar1"))
	paths = list(fingerprint = context$fingerprint)
	compiled = rotation = collapsed = suffix = context$backend
	RNGkind(
		kind = "Mersenne-Twister",
		normal.kind = "Inversion",
		sample.kind = "Rejection"
	)
	if (!grepl("^[A-Za-z0-9_]+$", label) || !grepl("^[A-Za-z0-9_]+$", subdir)) {
		stop("Invalid static fit label or directory.")
	}
	if (nscan %% odens != 0 || burn %% odens != 0 || nscan < odens || burn < 0) {
		stop(
			"Sampling and warmup lengths must be valid multiples of the storage interval."
		)
	}
	if (is.list(y)) {
		actors = rownames(y[[1]])
		years = names(y)
		if (any(
			!vapply(
				y,
				function(m) {
					identical(rownames(m), actors) && identical(colnames(m), actors) &&
						identical(dim(m), c(length(actors), length(actors)))
				},
				logical(1)
			)
		)) {
			stop("Input actor order must match in every period.")
		}
		y = array(
			unlist(y, use.names = FALSE),
			c(length(actors), length(actors), length(y)),
			dimnames = list(actors, actors, years)
		)
	}
	n = dim(y)[1]
	periods = dim(y)[3]
	actors = dimnames(y)[[1]]
	years = dimnames(y)[[3]]
	if (n < 2 || periods < 2 || dim(y)[2] != n || rank < 1) {
		stop("Invalid static network dimensions.")
	}
	mask = upper.tri(matrix(0, n, n))
	observed = array(FALSE, dim(y))
	for (t in seq_len(periods)) {
		m = y[, , t]
		if (!identical(is.finite(m), t(is.finite(m))) ||
			any(abs(m[mask] - t(m)[mask]) > 1e-10, na.rm = TRUE)) {
			stop("The static sampler requires symmetric observed values and missingness.")
		}
		observed[, , t] = mask & is.finite(m)
	}
	required = c(
		"beta_precision",
		"va_shape",
		"va_scale",
		"s2_shape",
		"s2_scale",
		"rho_mean",
		"rho_sd",
		"sigma_shape",
		"sigma_scale",
		"snap_kappa",
		"snap_pi_a",
		"snap_pi_b"
	)
	if (!all(required %in% names(prior)) || any(!is.finite(unlist(prior[required]))) ||
		any(unlist(prior[setdiff(required, "rho_mean")]) <= 0)) {
		stop("All resolved fixed prior constants are required.")
	}
	rho_bound = if (is.null(prior$rho_bound)) 0.99 else prior$rho_bound
	rho_fixed = prior$rho_fixed
	if (length(rho_bound) != 1L || !is.finite(rho_bound) || rho_bound <= 0 || rho_bound >= 1 ||
		(!is.null(rho_fixed) && (length(rho_fixed) != 1L || !is.finite(rho_fixed) || abs(rho_fixed) >= rho_bound))) {
		stop(
			"The persistence bound must lie inside (0, 1) and a fixed persistence inside the bound."
		)
	}
	exact_sigma = !is.null(rho_fixed) || rho_bound > 0.99
	code_files = c(
		file.path(context$directory, "00_sampler.R"),
		file.path(context$directory, "sampler.cpp"),
		context$archive
	)
	missing_plan = revision_missing_suffix_plan(y)
	bridge_plan = revision_missing_bridge_plan(y)
	bridge_every = 5L
	input_hash = digest::digest(y, algo = "sha256")
	if (!is.null(initial_state$parent_input_hash) && !identical(initial_state$parent_input_hash, input_hash)) {
		stop("Exact continuation requires the same observed network.")
	}
	stopifnot(
		length(indicator_every) == 1L,
		indicator_every %in% c(0L, 5L),
		indicator_every == 0L || latent_kind == "snap"
	)
	stopifnot(
		length(temporal_tau) == 1L,
		is.finite(temporal_tau),
		temporal_tau >= 0,
		temporal_tau == 0 || latent_kind == "snap"
	)
	config = list(
		temporal_tau = temporal_tau,
		temporal_every = 5L,
		temporal_moves = 5L,
		indicator_every = as.integer(indicator_every),
		seed = seed,
		nscan = nscan,
		burn = burn,
		odens = odens,
		input_hash = input_hash,
		rank = rank,
		prior = prior,
		initial_state = initial_state,
		latent_kind = latent_kind,
		backend_artifact = compiled$artifact,
		runtime = paths$fingerprint,
		code = tools::md5sum(code_files),
		rotation_kernel = "conditional_O2_plus_marginal_period_and_suffix_moves",
		rotation_artifact = rotation$artifact,
		missing_suffix_kernel = "joint_conditional_prior_after_last_augmentation_use",
		missing_suffix_plan = missing_plan,
		missing_bridge_kernel = "last_reset_mixture_then_gaussian_bridge",
		missing_bridge_plan = bridge_plan,
		bridge_every = bridge_every,
		joint_actor_kernel = "observed_gaussian_static_effect_and_entire_path_conditional_on_indicators",
		centroid_kernel = "level_marginalized_centroid_MALA_then_exact_joint_levels",
		centroid_step = .75,
		centroid_moves = 5L,
		centroid_every = 5L,
		frozen_shifts = frozen_shifts,
		jump_every = 5L,
		rho_kernel = if (is.null(rho_fixed)) "bounded_full_conditional_slice" else "fixed",
		rng_kind = RNGkind(),
		r_version = R.version.string
	)
	key = digest::digest(list(y = y, config = config), algo = "sha256")
	directory = file.path(output_dir, subdir, "fits")
	dir.create(directory, recursive = TRUE, showWarnings = FALSE)
	prefix = file.path(directory, paste0(label, "_", substr(key, 1, 16)))
	fit_file = paste0(prefix, ".rds")
	checkpoint_file = paste0(prefix, ".checkpoint.rds")
	status_file = paste0(prefix, ".status")
	if (file.exists(fit_file)) {
		fit = readRDS(fit_file)
		if (!identical(fit$revision_build$key, key)) {
			stop("Static cache key mismatch.")
		}
		return(fit)
	}
	namespace = asNamespace("lame")
	draw_snap = get("rUV_dynamic_snap_fc_cpp", envir = namespace)
	draw_ar = get("rUV_dynamic_fc_cpp", envir = namespace)
	draw_sigma = get(".lame_sample_sigma_uv_drift_only", envir = namespace)
	if (file.exists(checkpoint_file)) {
		state = readRDS(checkpoint_file)
		if (!identical(state$key, key)) {
			stop("Static checkpoint key mismatch.")
		}
		assign(".Random.seed", state$rng, envir = .GlobalEnv)
	} else {
		set.seed(seed)
		draws = nscan / odens
		u = array(
			rnorm(n * rank * periods, sd = 0.1),
			c(n, rank, periods),
			dimnames = list(actors, NULL, years)
		)
		state = list(
			key = key,
			config = config,
			iteration = 0L,
			started = Sys.time(),
			bridge_phase = 0L,
			u = u,
			a = setNames(rep(0, n), actors),
			beta = mean(y[observed]),
			s2 = max(stats::var(y[observed]), 1e-4),
			va = prior$va_scale / prior$va_shape,
			rho = if (is.null(rho_fixed)) 0.9 else rho_fixed,
			sigma = 0.1,
			pi = if (latent_kind == "snap") 0.001 else 0,
			delta = matrix(0, n, periods),
			draw_index = 0L,
			u_draws = array(
				NA_real_,
				c(n, rank, periods, draws),
				dimnames = list(actors, NULL, years, NULL)
			),
			a_draws = matrix(NA_real_, n, draws, dimnames = list(actors, NULL)),
			beta_draws = matrix(NA_real_, draws, 1, dimnames = list(NULL, "intercept")),
			snap_draws = array(as.raw(0), c(draws, n, periods), dimnames = list(NULL, actors, years)),
			trace = matrix(
				NA_real_,
				draws,
				6,
				dimnames = list(
					NULL,
					c("ve", "va", "rho_uv", "sigma_uv", "pi_snap", "observed_log_likelihood")
				)
			),
			monitor = data.frame(
				iteration = integer(),
				elapsed_seconds = numeric(),
				rho_uv = numeric(),
				sigma_uv = numeric(),
				pi_snap = numeric(),
				va = numeric()
			)
		)
		if (!is.null(initial_state)) {
			sampler_initial = unserialize(serialize(initial_state, NULL))
			for (name in c(
				"u",
				"a",
				"beta",
				"s2",
				"va",
				"rho",
				"sigma",
				"pi",
				"delta",
				"bridge_phase",
				"initial_joint_done"
			)) {
				if (!is.null(sampler_initial[[name]])) {
					state[[name]] = sampler_initial[[name]]
				}
			}
		}
		if (!identical(as.integer(dim(state$u)), as.integer(c(n, rank, periods))) ||
			length(state$a) != n || length(state$beta) != 1 ||
			!identical(as.integer(dim(state$delta)), as.integer(c(n, periods)))) {
			stop("Initial state dimensions do not match the static model.")
		}
		if (latent_kind == "ar1") {
			state$pi = 0
			state$delta[] = 0
		}
		if (!is.null(rho_fixed)) {
			state$rho = rho_fixed
		}
		if (length(state$bridge_phase) != 1L || !is.finite(state$bridge_phase) || !state$bridge_phase %in% 0:(bridge_every - 1L)) {
			stop("Invalid bridge schedule phase.")
		}
		if (any(
			!is.finite(
				c(
					state$u,
					state$a,
					state$beta,
					state$s2,
					state$va,
					state$rho,
					state$sigma,
					state$pi
				)
			)
		) || any(c(state$s2, state$va, state$sigma) <= 0) ||
			abs(state$rho) >= rho_bound || any(!state$delta %in% c(0, 1)) ||
			(latent_kind == "snap" && (state$pi <= 0 || state$pi >= 1))) {
			stop("Initial state falls outside the specified posterior support.")
		}
		if (!is.null(initial_state$rng)) {
			assign(".Random.seed", initial_state$rng, envir = .GlobalEnv)
		}
	}
	write_checkpoint = function(state) {
		state$rng = get(".Random.seed", envir = .GlobalEnv)
		temporary = paste0(checkpoint_file, ".tmp")
		saveRDS(state, temporary, compress = FALSE)
		if (!file.rename(temporary, checkpoint_file)) {
			stop("Could not finalize the exact static checkpoint.")
		}
		write.csv(state$monitor, paste0(prefix, ".monitor.csv"), row.names = FALSE)
	}
	total = burn + nscan
	stop_iteration = min(total, state$iteration + max_sweeps)
	writeLines(
		paste("running", format(Sys.time(), tz = "UTC"), "from", state$iteration),
		status_file
	)
	if (state$iteration < stop_iteration) {
		if (!isTRUE(state$initial_joint_done)) {
			state$initial_joint_counts = c(accepted = 0, proposed = 0)
			if (length(config$frozen_shifts)) {
				for (proposal in seq_len(20L)) {
					shift_index = sample.int(length(config$frozen_shifts), 1L)
					jump = compiled$jump_draw(
						y,
						state$u,
						config$frozen_shifts[[shift_index]],
						state$s2,
						state$va,
						prior$beta_precision,
						state$rho,
						state$sigma,
						prior$snap_kappa,
						state$pi
					)
					state$u = jump$u
					state$a = as.vector(jump$a)
					state$beta = jump$beta
					state$delta = jump$delta
					state$initial_joint_counts = state$initial_joint_counts + c(jump$accepted, 1L)
				}
			}
			state$initial_joint_done = TRUE
		}
		for (iteration in seq.int(state$iteration + 1L, stop_iteration)) {
			beta_path = rep(state$beta, periods)
			a_path = matrix(state$a, n, periods)
			prepared = compiled$revision_gaussian_prepare_cpp(
				y,
				state$u,
				beta_path,
				a_path,
				state$s2
			)
			state$s2 = 1 / rgamma(
				1,
				shape = prior$s2_shape + prepared$n_dyads / 2,
				rate = prior$s2_scale + prepared$ssr / 2
			)
			levels = revision_static_level_draw(prepared, state$s2, state$va, prior$beta_precision)
			state$beta = levels$beta
			state$a = levels$a
			state$va = 1 / rgamma(
				1,
				shape = prior$va_shape + n / 2,
				rate = prior$va_scale + sum(state$a^2) / 2
			)
			residual = compiled$revision_gaussian_residual_cpp(
				prepared$z,
				rep(state$beta, periods),
				matrix(state$a, n, periods)
			)
			if (latent_kind == "snap") {
				uv = draw_snap(
					state$u,
					state$u,
					residual,
					state$rho,
					state$sigma,
					state$s2,
					prior$snap_kappa,
					state$pi,
					state$delta,
					state$delta,
					TRUE,
					TRUE
				)
				state$delta = uv$delta_u
				n_snaps = sum(state$delta[, -1, drop = FALSE])
				state$pi = rbeta(
					1,
					prior$snap_pi_a + n_snaps,
					prior$snap_pi_b + n * (periods - 1) - n_snaps
				)
			} else {
				uv = draw_ar(
					state$u,
					state$u,
					residual,
					state$rho,
					state$sigma,
					state$s2,
					TRUE,
					TRUE
				)
			}
			state$u = uv$U
			if (any(!is.finite(state$u))) {
				stop("Nonfinite latent state; no silent fallback is permitted.")
			}
			state$u = rotation$revision_period_rotation_stable_cpp(
				state$u,
				state$delta,
				state$rho,
				state$sigma
			)
			if (latent_kind == "snap") {
				orbit = collapsed$revision_collapsed_sweep_cpp(
					state$u,
					state$delta,
					state$rho,
					state$sigma,
					prior$snap_kappa,
					state$pi,
					1L
				)
				state$u = orbit$u
				state$delta = orbit$delta
				if (is.null(state$orbit_counts)) state$orbit_counts = rep(0, 4)
				state$orbit_counts = state$orbit_counts + as.numeric(orbit$counts)
				tail = suffix$revision_suffix_sweep_cpp(
					state$u,
					state$delta,
					state$rho,
					state$sigma,
					prior$snap_kappa,
					rep(qlogis(state$pi), periods),
					1L
				)
				state$u = tail$u
				state$delta = tail$delta
				if (is.null(state$suffix_counts)) state$suffix_counts = rep(0, 4)
				state$suffix_counts = state$suffix_counts + as.numeric(tail$counts)
			}
			if (is.null(rho_fixed)) {
				state$rho = revision_rho_uv_slice(
					state$u,
					state$u,
					state$sigma,
					state$rho,
					TRUE,
					state$delta,
					NULL,
					prior$rho_mean,
					prior$rho_sd,
					rho_bound
				)
			}
			state$sigma = if (exact_sigma) {
				revision_sigma_uv_draw(
					state$u,
					state$rho,
					state$delta,
					prior$sigma_shape,
					prior$sigma_scale
				)
			} else {
				draw_sigma(
					state$u,
					state$u,
					state$rho,
					state$sigma,
					TRUE,
					state$delta,
					NULL,
					prior$sigma_shape,
					prior$sigma_scale
				)
			}
			# discard the old augmentation before this observed-target update
			missing = revision_missing_suffix_draw(
				state$u,
				state$delta,
				missing_plan,
				state$rho,
				state$sigma,
				prior$snap_kappa,
				qlogis(state$pi)
			)
			state$u = missing$u
			state$delta = missing$delta
			if (is.null(state$missing_suffix_counts)) state$missing_suffix_counts = c(actors = 0, transitions = 0, initial_paths = 0)
			state$missing_suffix_counts = state$missing_suffix_counts +
				c(
					missing$updated_actors,
					missing$updated_transitions,
					missing$updated_initial_paths
				)
			if ((iteration + state$bridge_phase) %% bridge_every == 0L) {
				bridge = revision_missing_bridge_draw(
					state$u,
					state$delta,
					bridge_plan,
					state$rho,
					state$sigma,
					prior$snap_kappa,
					if (latent_kind == "snap") qlogis(state$pi) else -Inf
				)
				state$u = bridge$u
				state$delta = bridge$delta
				if (is.null(state$missing_bridge_counts)) state$missing_bridge_counts = c(blocks = 0, positions = 0, indicators = 0)
				state$missing_bridge_counts = state$missing_bridge_counts + c(bridge$updated_blocks, bridge$updated_positions, bridge$updated_indicators)
			}
			if (config$indicator_every > 0L && (iteration + state$bridge_phase) %% config$indicator_every == 0L) {
				joint = compiled$indicator_pair_sweep(
					y,
					state$u,
					state$a,
					state$beta,
					state$s2,
					state$va,
					state$rho,
					state$sigma,
					prior$snap_kappa,
					state$pi,
					state$delta
				)
				state$delta = joint$delta
			} else {
				joint = compiled$actor_path_sweep(
					y,
					state$u,
					state$a,
					state$beta,
					state$s2,
					state$va,
					state$rho,
					state$sigma,
					prior$snap_kappa,
					state$delta
				)
			}
			state$u = joint$u
			state$a = as.vector(joint$a)
			if ((iteration + state$bridge_phase) %% config$centroid_every == 0L) {
				centroid = compiled$centroid_prepare(
					y,
					state$u,
					state$s2,
					state$va,
					prior$beta_precision,
					state$rho,
					state$sigma,
					prior$snap_kappa,
					state$delta
				)
				m = as.vector(apply(state$u, c(2, 3), mean))
				move = compiled$centroid_draw(
					centroid,
					m,
					config$centroid_step,
					config$centroid_moves
				)
				state$u = move$u
				state$a = as.vector(move$a)
				state$beta = move$beta
				if (is.null(state$centroid_counts)) state$centroid_counts = c(accepted = 0, proposed = 0)
				state$centroid_counts = state$centroid_counts + c(move$accepted, move$proposed)
			}
			if (length(config$frozen_shifts) && (iteration + state$bridge_phase) %% config$jump_every == 0L) {
				shift_index = sample.int(length(config$frozen_shifts), 1L)
				jump = compiled$jump_draw(
					y,
					state$u,
					config$frozen_shifts[[shift_index]],
					state$s2,
					state$va,
					prior$beta_precision,
					state$rho,
					state$sigma,
					prior$snap_kappa,
					state$pi
				)
				state$u = jump$u
				state$a = as.vector(jump$a)
				state$beta = jump$beta
				state$delta = jump$delta
				if (is.null(state$jump_counts)) state$jump_counts = c(accepted = 0, proposed = 0)
				state$jump_counts = state$jump_counts + c(jump$accepted, 1L)
			}
			if (config$temporal_tau > 0 && (iteration + state$bridge_phase) %% config$temporal_every == 0L) {
				move = compiled$temporal_draw(
					y,
					state$u,
					config$temporal_tau,
					config$temporal_moves,
					state$s2,
					state$va,
					prior$beta_precision,
					state$rho,
					state$sigma,
					prior$snap_kappa,
					state$pi
				)
				state$u = move$u
				state$a = as.vector(move$a)
				state$beta = move$beta
				state$delta = move$delta
				if (is.null(state$temporal_counts)) state$temporal_counts = c(accepted = 0, proposed = 0)
				state$temporal_counts = state$temporal_counts + c(move$accepted, move$proposed)
			}
			state$iteration = iteration
			if (iteration > burn && iteration %% odens == 0) {
				index = state$draw_index + 1L
				state$draw_index = index
				state$u_draws[, , , index] = state$u
				state$a_draws[, index] = state$a
				state$beta_draws[index, ] = state$beta
				state$snap_draws[index, , ] = as.raw(state$delta)
				log_likelihood = compiled$revision_gaussian_loglik_cpp(
					y,
					state$u,
					rep(state$beta, periods),
					matrix(state$a, n, periods),
					state$s2
				)
				state$trace[index, ] = c(state$s2, state$va, state$rho, state$sigma, state$pi, log_likelihood)
			}
			if (iteration %% checkpoint_every == 0 || iteration == stop_iteration) {
				state$monitor = rbind(
					state$monitor,
					data.frame(
						iteration = iteration,
						elapsed_seconds = as.numeric(difftime(Sys.time(), state$started, units = "secs")),
						rho_uv = state$rho,
						sigma_uv = state$sigma,
						pi_snap = state$pi,
						va = state$va
					)
				)
				write_checkpoint(state)
			}
		}
	}
	if (state$iteration < total) {
		writeLines(paste("checkpointed", state$iteration, "of", total), status_file)
		return(
			list(
				complete = FALSE,
				checkpoint_file = checkpoint_file,
				iteration = state$iteration
			)
		)
	}
	fit = list(
		complete = TRUE,
		method = "joint_actor_centroid_and_fixed_symmetric_proposals",
		jump_counts = state$jump_counts,
		initial_joint_counts = state$initial_joint_counts,
		orbit_counts = state$orbit_counts,
		suffix_counts = state$suffix_counts,
		centroid_counts = state$centroid_counts,
		missing_suffix_counts = state$missing_suffix_counts,
		missing_bridge_counts = state$missing_bridge_counts,
		family = "normal",
		symmetric = TRUE,
		dynamic_uv = TRUE,
		dynamic_uv_kind = latent_kind,
		dynamic_ab = FALSE,
		dynamic_beta = FALSE,
		Y = y,
		prior = prior,
		U_draws = state$u_draws,
		U = apply(state$u_draws, c(1, 2, 3), mean),
		revision_a_draws = state$a_draws,
		APM = rowMeans(state$a_draws),
		BETA = state$beta_draws,
		snap_draws = state$snap_draws,
		snap_prob = apply(
			array(
				as.integer(state$snap_draws),
				dim(state$snap_draws),
				dimnames = dimnames(state$snap_draws)
			),
			c(2, 3),
			mean
		),
		trace = state$trace,
		VC = state$trace[, c("ve", "va"), drop = FALSE],
		rho_uv = state$trace[, "rho_uv"],
		sigma_uv = state$trace[, "sigma_uv"],
		pi_snap = state$trace[, "pi_snap"],
		monitor = state$monitor,
		tryErrorChecks = list(UV = 0L),
		revision_exact_symmetric_target = list(
			likelihood = "one gaussian observation per unordered off-diagonal dyad",
			additive_prior = "a_i iid normal(0, va); shared sender and receiver effect",
			beta_prior_precision = matrix(prior$beta_precision, 1, 1),
			va_prior_shape = prior$va_shape,
			va_prior_scale = prior$va_scale,
			s2_prior_shape = prior$s2_shape,
			s2_prior_scale = prior$s2_scale,
			diagonal_likelihood = FALSE,
			missing_dyads_mirrored = TRUE,
			missing_trajectory_kernel = config$missing_suffix_kernel,
			missing_bridge_kernel = config$missing_bridge_kernel,
			bridge_every = config$bridge_every,
			initial_latent_prior = "normal(0, sigma_uv^2 / (1 - rho_uv^2))",
			rho_support = c(-rho_bound, rho_bound),
			rho_kernel = config$rho_kernel
		),
		start_vals = state[c("u", "a", "beta", "s2", "va", "rho", "sigma", "pi", "delta")],
		revision_build = list(
			key = key,
			label = label,
			seed = seed,
			config = config,
			args = list(
				R = rank,
				family = "normal",
				symmetric = TRUE,
				dynamic_uv = TRUE,
				dynamic_uv_kind = latent_kind,
				dynamic_ab = FALSE,
				dynamic_beta = FALSE,
				nscan = nscan,
				burn = burn,
				odens = odens,
				prior = prior
			),
			input_hash = input_hash,
			runtime_hash = paths$fingerprint,
			started = state$started,
			finished = Sys.time(),
			elapsed_seconds = as.numeric(difftime(Sys.time(), state$started, units = "secs")),
			fit_file = fit_file,
			checkpoint_file = checkpoint_file
		)
	)
	if (!is.null(rho_fixed)) {
		fit$revision_exact_symmetric_target$rho_fixed = rho_fixed
	}
	if (exact_sigma) {
		fit$revision_exact_symmetric_target$sigma_kernel = "exact_stationary_inverse_gamma"
	}
	temporary = paste0(fit_file, ".tmp")
	saveRDS(fit, temporary, compress = FALSE)
	if (!file.rename(temporary, fit_file)) {
		stop("Could not finalize the static fit cache.")
	}
	writeLines(
		paste("completed", format(Sys.time(), tz = "UTC"), state$iteration),
		status_file
	)
	fit
}

continue_model = function(
	checkpoint_file,
	y,
	label,
	context,
	output_dir,
	nscan,
	burn = 0,
	subdir = "fits",
	checkpoint_every = 500,
	max_sweeps = Inf,
	allow_partial = FALSE
) {
	state = readRDS(checkpoint_file)
	config = state$config
	partial = state$iteration != config$burn + config$nscan
	if ((partial && (!isTRUE(allow_partial) || burn <= 0)) || is.null(state$rng)) {
		stop(
			"A new segment requires a completed checkpoint, or explicit warmup after a partial checkpoint."
		)
	}
	if (!identical(unname(tools::md5sum(names(config$code))), unname(config$code)) ||
		!identical(context$fingerprint, config$runtime)) {
		stop("The checkpoint code or installed dependency has changed.")
	}
	initial = state[c(
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
	)]
	initial$bridge_phase = as.integer((state$bridge_phase + state$iteration) %% 5L)
	initial$parent_checkpoint_sha256 = digest::digest(file = checkpoint_file, algo = "sha256")
	initial$parent_input_hash = config$input_hash
	if (partial) initial$partial_parent_iteration = state$iteration
	fit_model(
		y,
		label,
		config$seed,
		config$prior,
		context,
		output_dir,
		nscan,
		burn,
		config$odens,
		config$rank,
		config$latent_kind,
		checkpoint_every,
		max_sweeps,
		subdir,
		initial,
		config$frozen_shifts,
		if (is.null(config$indicator_every)) 0L else config$indicator_every,
		if (is.null(config$temporal_tau)) 0 else config$temporal_tau
	)
}
