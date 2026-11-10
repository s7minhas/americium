# validate the sampler against independent calculations

args = commandArgs(trailingOnly = TRUE)
argument = sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])[1]
root = if (length(args)) normalizePath(args[1]) else dirname(dirname(normalizePath(argument)))
source(file.path(root, "code", "00_sampler.R"))
source(file.path(root, "code", "10_models.R"))
context = model_load(file.path(root, "code"))
invisible(list2env(as.list(context$backend), envir = environment()))
output = file.path(root, "results", "checks", "model")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
# validate actor path
dense_reference = function(y, u, a, beta, s2, va, rho, sigma, kappa, delta, actor) {
	n = dim(u)[1]
	rank = dim(u)[2]
	periods = dim(u)[3]
	p = rank * periods + 1L
	precision = matrix(0, p, p)
	linear = numeric(p)
	precision[p, p] = 1 / va
	precision[seq_len(rank), seq_len(rank)] = diag((1 - rho^2) / sigma^2, rank)
	for (t in seq_len(periods)) {
		at = (t - 1L) * rank + seq_len(rank)
		if (t > 1) {
			reset = delta[actor, t] == 1
			innovation = matrix(0, rank, p)
			innovation[, at] = diag(rank)
			if (!reset) innovation[, at - rank] = -rho * diag(rank)
			precision = precision + crossprod(innovation) / if (reset) kappa^2 else sigma^2
		}
		for (j in setdiff(seq_len(n), actor)) {
			if (!is.finite(y[actor, j, t])) next
			x = numeric(p)
			x[at] = u[j, , t]
			x[p] = 1
			precision = precision + tcrossprod(x) / s2
			linear = linear + x * (y[actor, j, t] - beta - a[j]) / s2
		}
	}
	covariance = solve(precision)
	list(mean = as.vector(covariance %*% linear), covariance = covariance)
}

set.seed(130926)
cases = expand.grid(rank = 1:3, periods = c(1L, 4L), rho = c(-.8, .7))
checks = list()
moments = list()
index = 0L
for (case in seq_len(nrow(cases))) {
	rank = cases$rank[case]
	periods = cases$periods[case]
	n = 4L
	u = array(rnorm(n * rank * periods), c(n, rank, periods))
	a = rnorm(n, sd = .3)
	beta = -.2
	s2 = .4
	va = .3
	rho = cases$rho[case]
	sigma = .2
	kappa = 1.2
	delta = matrix(rbinom(n * periods, 1, .35), n, periods)
	delta[, 1] = 0
	y = array(NA_real_, c(n, n, periods))
	for (t in seq_len(periods)) {
		m = matrix(rnorm(n * n), n, n)
		m = (m + t(m)) / 2
		diag(m) = NA_real_
		if (t %% 2 == 0) {
			m[1, 3] = m[3, 1] = NA_real_
		}
		y[, , t] = m
	}
	for (actor in 1:2) {
		index = index + 1L
		reference = dense_reference(y, u, a, beta, s2, va, rho, sigma, kappa, delta, actor)
		actual = actor_path_conditional(
			y,
			u,
			a,
			beta,
			s2,
			va,
			rho,
			sigma,
			kappa,
			delta,
			actor - 1L,
			20000L
		)
		error_mean = max(abs(reference$mean - actual$mean))
		error_covariance = max(abs(reference$covariance - actual$covariance))
		checks[[index]] = data.frame(
			case = case,
			actor = actor,
			rank = rank,
			periods = periods,
			rho = rho,
			mean_error = error_mean,
			covariance_error = error_covariance,
			passed = error_mean < 1e-10 && error_covariance < 1e-10
		)
		z = sweep(actual$draws, 2, reference$mean) %*% solve(chol(reference$covariance))
		z_mean = colMeans(z) * sqrt(nrow(z))
		z_variance = (colMeans(z^2) - 1) * sqrt(nrow(z) / 2)
		products = crossprod(z) / nrow(z)
		z_cross = products[upper.tri(products)] * sqrt(nrow(z))
		all_z = c(z_mean, z_variance, z_cross)
		moments[[index]] = data.frame(
			case = case,
			actor = actor,
			moment = seq_along(all_z),
			z = all_z,
			passed = abs(all_z) < 4.5
		)
	}
	before = serialize(list(u = u, a = a), NULL)
	set.seed(700 + case)
	first = actor_path_sweep(y, u, a, beta, s2, va, rho, sigma, kappa, delta)
	seed = .Random.seed
	set.seed(700 + case)
	second = actor_path_sweep(y, u, a, beta, s2, va, rho, sigma, kappa, delta)
	stopifnot(
		identical(first, second),
		identical(.Random.seed, seed),
		identical(before, serialize(list(u = u, a = a), NULL))
	)
}

checks = do.call(rbind, checks)
moments = do.call(rbind, moments)
write.csv(
	checks,
	file.path(output, "actor_path_analytic_checks.csv"),
	row.names = FALSE
)
write.csv(
	moments,
	file.path(output, "actor_path_draw_moments.csv"),
	row.names = FALSE
)
stopifnot(all(checks$passed), all(moments$passed))

# check centroid
centroid_reference = function(y, shape, m, s2, va, beta_precision, rho, sigma, kappa, delta) {
	n = dim(shape)[1]
	rank = dim(shape)[2]
	periods = dim(shape)[3]
	p = rank * periods
	a = j = list()
	residual = numeric()
	at = 0L
	for (t in seq_len(periods)) {
		x = matrix(shape[, , t], n, rank)
		mu = m[(t - 1L) * rank + seq_len(rank)]
		for (i in seq_len(n - 1L)) {
			for (k in seq.int(i + 1L, n)) {
				if (!is.finite(y[i, k, t])) next
				at = at + 1L
				row = numeric(n + 1)
				row[c(1, i + 1, k + 1)] = 1
				a[[at]] = row
				row = numeric(p)
				row[(t - 1L) * rank + seq_len(rank)] = x[i, ] + x[k, ] + 2 * mu
				j[[at]] = row
				residual[at] = y[i, k, t] - sum((x[i, ] + mu) * (x[k, ] + mu))
			}
		}
	}
	a = do.call(rbind, a)
	j = do.call(rbind, j)
	level_prior = diag(c(1 / beta_precision, rep(va, n)))
	covariance = diag(s2, nrow(a)) + a %*% level_prior %*% t(a)
	weight = solve(covariance)
	level_mean = as.vector(level_prior %*% crossprod(a, weight %*% residual))
	level_covariance = level_prior - level_prior %*% crossprod(a, weight %*% a) %*% level_prior
	prior_precision = matrix(0, p, p)
	prior_linear = numeric(p)
	prior_value = 0
	for (t in seq_len(periods)) {
		for (i in seq_len(n)) {
			current = (t - 1L) * rank + seq_len(rank)
			design = matrix(0, rank, p)
			design[, current] = diag(rank)
			constant = shape[i, , t]
			variance = if (t == 1) sigma^2 / (1 - rho^2) else if (delta[i, t] == 1) kappa^2 else sigma^2
			if (t > 1 && delta[i, t] == 0) {
				design[, current - rank] = -rho * diag(rank)
				constant = constant - rho * shape[i, , t - 1L]
			}
			prior_precision = prior_precision + crossprod(design) / variance
			prior_linear = prior_linear - as.vector(crossprod(design, constant)) / variance
			prior_value = prior_value + sum(dnorm(constant + design %*% m, sd = sqrt(variance), log = TRUE))
		}
	}
	log_density = -.5 * (as.numeric(determinant(covariance, logarithm = TRUE)$modulus) +
		sum(residual * (weight %*% residual)) + length(residual) * log(2 * pi)) + prior_value
	gradient = as.vector(crossprod(j, weight %*% residual)) - as.vector(prior_precision %*% m) + prior_linear
	precision = prior_precision + crossprod(j, weight %*% j)
	list(
		log_density = log_density,
		gradient = gradient,
		precision = precision,
		level_mean = level_mean,
		level_covariance = level_covariance
	)
}

proposal_reference = function(to, from, evaluation, step) {
	mean = from + .5 * step^2 * solve(evaluation$precision, evaluation$gradient)
	difference = to - mean
	.5 * as.numeric(determinant(evaluation$precision, logarithm = TRUE)$modulus) -
		length(from) * log(step) - .5 * sum(difference * (evaluation$precision %*% difference)) / step^2
}

set.seed(92471)
grid = expand.grid(
	rank = 1:3,
	periods = c(1L, 3L),
	rho = c(-.7, .8),
	small_shape = c(FALSE, TRUE)
)
checks = list()
for (case in seq_len(nrow(grid))) {
	g = grid[case, ]
	n = 4L
	rank = g$rank
	periods = g$periods
	u = array(
		rnorm(n * rank * periods, sd = if (g$small_shape) 1e-7 else .4),
		c(n, rank, periods)
	)
	u = sweep(u, c(2, 3), matrix(rnorm(rank * periods, sd = .3), rank, periods), "+")
	y = array(NA_real_, c(n, n, periods))
	for (t in seq_len(periods)) {
		x = matrix(rnorm(n * n), n)
		y[, , t] = (x + t(x)) / 2
		diag(y[, , t]) = NA_real_
	}
	if (periods > 1) y[1, 2, 2] = y[2, 1, 2] = NA_real_
	delta = matrix(rbinom(n * periods, 1, .3), n, periods)
	delta[, 1] = 0
	s2 = .3
	va = .4
	beta_precision = .2
	sigma = .4
	kappa = 1.3
	rho = g$rho
	object = centroid_prepare(y, u, s2, va, beta_precision, rho, sigma, kappa, delta)
	zero = rep(0, rank * periods)
	m = rnorm(length(zero), sd = .25)
	actual = centroid_evaluate(object, m)
	initial = centroid_evaluate(object, zero)
	ref = function(at) {
		centroid_reference(
			y,
			actual$shape,
			at,
			s2,
			va,
			beta_precision,
			rho,
			sigma,
			kappa,
			delta
		)
	}
	reference = ref(m)
	reference_zero = ref(zero)
	finite = vapply(
		seq_along(m),
		function(k) {
			plus = minus = m
			plus[k] = plus[k] + 1e-6
			minus[k] = minus[k] - 1e-6
			(centroid_evaluate(object, plus)$log_density - centroid_evaluate(object, minus)$log_density) / 2e-6
		},
		0.
	)
	proposal = rnorm(length(m), sd = .3)
	step = .7
	q_actual = centroid_proposal_log(object, proposal, m, step)
	q_ref = proposal_reference(proposal, m, reference, step)
	reverse_actual = centroid_proposal_log(object, m, proposal, step)
	next_actual = centroid_evaluate(object, proposal)
	ratio = next_actual$log_density - actual$log_density + reverse_actual - q_actual
	flux_forward = reference$log_density + q_ref + min(0, ratio)
	flux_reverse = ref(proposal)$log_density + proposal_reference(m, proposal, ref(proposal), step) + min(0, -ratio)
	errors = c(
		density_difference = abs(
			actual$log_density - initial$log_density - reference$log_density + reference_zero$log_density
		),
		gradient = max(abs(actual$gradient - reference$gradient)),
		finite_difference = max(abs(actual$gradient - finite)),
		metric = max(abs(actual$precision - reference$precision)),
		level_mean = max(abs(actual$level_mean - reference$level_mean)),
		level_covariance = max(abs(actual$level_covariance - reference$level_covariance)),
		proposal = abs(q_actual - q_ref),
		full_target_flux = abs(flux_forward - flux_reverse)
	)
	checks[[case]] = data.frame(
		case = case,
		rank = rank,
		periods = periods,
		rho = rho,
		small_shape = g$small_shape,
		check = names(errors),
		error = errors,
		passed = errors < ifelse(names(errors) == "finite_difference", 1e-6, 1e-8)
	)
	before = serialize(u, NULL)
	set.seed(case)
	first = centroid_draw(object, m, step)
	rng = .Random.seed
	set.seed(case)
	second = centroid_draw(object, m, step)
	stopifnot(
		identical(first, second),
		identical(rng, .Random.seed),
		identical(serialize(u, NULL), before)
	)
}

checks = do.call(rbind, checks)
write.csv(
	checks,
	file.path(output, "centroid_dense_checks.csv"),
	row.names = FALSE
)
stopifnot(all(checks$passed))

# check centroid extreme
set.seed(33952)
checks = list()
counter = 0L
for (offset in c(0, 2, 20)) {
	for (shape_scale in c(1e-7, .4)) {
		for (rho in c(-.8, .9899)) {
			n = 7L
			rank = 2L
			periods = 4L
			u = array(
				rnorm(n * rank * periods, sd = shape_scale) + offset,
				c(n, rank, periods)
			)
			delta = matrix(rbinom(n * periods, 1, .2), n, periods)
			delta[, 1] = 0
			y = array(NA_real_, c(n, n, periods))
			for (t in seq_len(periods)) {
				a = matrix(rnorm(n * n, sd = .1), n)
				y[, , t] = .5 + (a + t(a)) / 2
				diag(y[, , t]) = NA_real_
			}
			y[1, 2, 2] = y[2, 1, 2] = NA_real_
			s2 = .008
			va = .02
			beta_precision = .8
			sigma = .08
			kappa = 2
			object = centroid_prepare(y, u, s2, va, beta_precision, rho, sigma, kappa, delta)
			m = as.vector(apply(u, c(2, 3), mean))
			actual = centroid_evaluate(object, m)
			at = m + rnorm(length(m), sd = .003)
			ref = function(z) {
				centroid_reference(
					y,
					actual$shape,
					z,
					s2,
					va,
					beta_precision,
					rho,
					sigma,
					kappa,
					delta
				)
			}
			reference = ref(m)
			errors = c(
				density_difference = abs(
					(centroid_evaluate(object, at)$log_density - actual$log_density) - (ref(at)$log_density - reference$log_density)
				),
				gradient_relative = max(
					abs(actual$gradient - reference$gradient) / pmax(1, abs(reference$gradient))
				),
				metric_relative = max(abs(actual$precision - reference$precision)) / max(abs(reference$precision)),
				level_mean_relative = max(abs(actual$level_mean - reference$level_mean)) / max(1, max(abs(reference$level_mean))),
				positive_metric = max(
					0,
					-min(eigen(actual$precision, symmetric = TRUE, only.values = TRUE)$values)
				)
			)
			counter = counter + 1L
			checks[[counter]] = data.frame(
				offset,
				shape_scale,
				rho,
				check = names(errors),
				error = errors,
				passed = errors < ifelse(names(errors) == "density_difference", 1e-5, 1e-7)
			)
		}
	}
}

checks = do.call(rbind, checks)
write.csv(
	checks,
	file.path(output, "centroid_numeric_checks.csv"),
	row.names = FALSE
)
stopifnot(all(checks$passed))

# check centroid stationarity
set.seed(602839)
checks = list()
receipts = list()
for (case in 1:6) {
	n = 3L
	rank = if (case %in% c(3, 6)) 2L else 1L
	periods = if (case %in% c(1, 4)) 1L else 2L
	p = rank * periods
	s2 = .8
	va = .4
	beta_precision = .5
	rho = if (case %% 2) .7 else -.5
	sigma = .5
	kappa = 1
	u = array(rnorm(n * p, sd = .25), c(n, rank, periods))
	for (t in seq_len(periods)) u[, , t] = sweep(matrix(u[, , t], n, rank), 2, colMeans(matrix(u[, , t], n, rank)))
	delta = matrix(0, n, periods)
	if (periods > 1) delta[1, 2] = 1
	y = array(NA_real_, c(n, n, periods))
	if (case <= 3) {
		for (t in seq_len(periods)) {
			z = matrix(rnorm(n * n, sd = .2), n)
			y[, , t] = .1 + (z + t(z)) / 2
			diag(y[, , t]) = NA_real_
		}
	}
	if (case == 2) y[1, 2, 2] = y[2, 1, 2] = NA_real_
	object = centroid_prepare(y, u, s2, va, beta_precision, rho, sigma, kappa, delta)
	shape = centroid_evaluate(object, rep(0, p))$shape
	q = matrix(0, p, p)
	h = numeric(p)
	for (t in seq_len(periods)) {
		for (i in seq_len(n)) {
			current = (t - 1L) * rank + seq_len(rank)
			d = matrix(0, rank, p)
			d[, current] = diag(rank)
			constant = shape[i, , t]
			variance = if (t == 1) sigma^2 / (1 - rho^2) else if (delta[i, t]) kappa^2 else sigma^2
			if (t > 1 && delta[i, t] == 0) {
				d[, current - rank] = -rho * diag(rank)
				constant = constant - rho * shape[i, , t - 1L]
			}
			q = q + crossprod(d) / variance
			h = h - as.vector(crossprod(d, constant)) / variance
		}
	}
	prior_mean = as.vector(solve(q, h))
	prior_root = chol(solve(q))
	edges = which(
		is.finite(y) & array(rep(upper.tri(y[, , 1]), periods), dim(y)),
		arr.ind = TRUE
	)
	a = matrix(0, nrow(edges), n + 1L)
	if (nrow(a)) a[, 1] = 1
	for (k in seq_len(nrow(edges))) a[k, edges[k, 1:2] + 1L] = 1
	level_prior = diag(c(1 / beta_precision, rep(va, n)))
	level_cov = solve(solve(level_prior) + crossprod(a) / s2)
	level_root = chol(level_cov)
	covariance = diag(s2, nrow(a)) + a %*% level_prior %*% t(a)
	weight = if (nrow(a)) solve(covariance) else matrix(numeric(), 0, 0)
	residuals = function(m) {
		answer = matrix(0, nrow(m), nrow(edges))
		for (k in seq_len(nrow(edges))) {
			i = edges[k, 1]
			j = edges[k, 2]
			t = edges[k, 3]
			now = m[, (t - 1L) * rank + seq_len(rank), drop = FALSE]
			answer[, k] = y[i, j, t] - sum(shape[i, , t] * shape[j, , t]) - as.vector(now %*% (shape[i, , t] + shape[j, , t])) - rowSums(now^2)
		}
		answer
	}
	count = 30000L
	target = list()
	total = 0L
	proposals = 0L
	while (total < count) {
		m = sweep(matrix(rnorm(5000 * p), 5000, p) %*% prior_root, 2, prior_mean, "+")
		r = residuals(m)
		log_accept = -.5 * rowSums(r * (r %*% weight))
		keep = log(runif(nrow(m))) < log_accept
		target[[length(target) + 1L]] = m[keep, , drop = FALSE]
		total = total + sum(keep)
		proposals = proposals + nrow(m)
		stopifnot(proposals < 1e7)
	}
	target = do.call(rbind, target)[seq_len(count), , drop = FALSE]
	level_mean = residuals(target) %*% a %*% level_cov / s2
	levels = level_mean + matrix(rnorm(count * (n + 1L)), count, n + 1L) %*% level_root
	before = cbind(target, levels)
	after = matrix(0, count, ncol(before))
	acceptance = numeric(count)
	step = c(.3, .8, 1.3)[(case - 1L) %% 3 + 1L]
	for (k in seq_len(count)) {
		draw = centroid_draw(object, target[k, ], step)
		after[k, ] = c(draw$centroids, draw$beta, draw$a)
		acceptance[k] = draw$mean_acceptance
	}
	features = function(x) {
		ans = cbind(x, x^2)
		colnames(ans) = c(paste0("mean_", seq_len(ncol(x))), paste0("square_", seq_len(ncol(x))))
		for (i in seq_len(ncol(x) - 1L)) {
			for (j in seq.int(i + 1L, ncol(x))) {
				ans = cbind(ans, x[, i] * x[, j])
				colnames(ans)[ncol(ans)] = paste0("cross_", i, "_", j)
			}
		}
		ans
	}
	difference = features(after) - features(before)
	mean = colMeans(difference)
	se = apply(difference, 2, sd) / sqrt(count)
	z = mean / se
	checks[[case]] = data.frame(
		case,
		rank,
		periods,
		observed_edges = nrow(edges),
		step,
		check = colnames(difference),
		difference = mean,
		mcse = se,
		z,
		passed = abs(z) < 4.5
	)
	receipts[[case]] = data.frame(
		case,
		independent_target_draws = count,
		rejection_proposals = proposals,
		step,
		mean_acceptance = mean(acceptance),
		gaussian_target = !nrow(edges)
	)
}

checks = do.call(rbind, checks)
receipts = do.call(rbind, receipts)
write.csv(
	checks,
	file.path(output, "centroid_stationarity_checks.csv"),
	row.names = FALSE
)
write.csv(
	receipts,
	file.path(output, "centroid_stationarity_targets.csv"),
	row.names = FALSE
)
stopifnot(all(checks$passed))

# check jump
set.seed(1944002)
reference = function(y, u, s2, va, beta_precision, rho, sigma, kappa, pi) {
	n = dim(u)[1]
	r = dim(u)[2]
	periods = dim(u)[3]
	level_design = matrix(0, 0, n + 1)
	response = numeric()
	for (t in 1:periods) {
		for (j in 2:n) {
			for (i in 1:(j - 1)) {
				if (is.finite(y[i, j, t])) {
					a = numeric(n + 1)
					a[c(1, i + 1, j + 1)] = 1
					level_design = rbind(level_design, a)
					response = c(response, y[i, j, t] - sum(u[i, , t] * u[j, , t]))
				}
			}
		}
	}
	level_prior_precision = diag(c(beta_precision, rep(1 / va, n)))
	level_precision = level_prior_precision + crossprod(level_design) / s2
	mean = solve(level_precision, crossprod(level_design, response) / s2)
	reference_covariance = diag(s2, nrow(level_design)) + level_design %*% solve(level_prior_precision, t(level_design))
	loglik = -.5 * drop(crossprod(response, solve(reference_covariance, response)))
	initial = sum(dnorm(u[, , 1], 0, sigma / sqrt(1 - rho^2), log = TRUE))
	prior = initial
	prob = matrix(0, n, periods)
	for (t in 2:periods) {
		for (i in 1:n) {
			d = log1p(-pi) + sum(dnorm(u[i, , t], rho * u[i, , t - 1], sigma, log = TRUE))
			z = log(pi) + sum(dnorm(u[i, , t], 0, kappa, log = TRUE))
			top = max(d, z)
			mix = top + log(exp(d - top) + exp(z - top))
			prior = prior + mix
			prob[i, t] = exp(z - mix)
		}
	}
	list(
		log_density = loglik + prior,
		level_mean = mean,
		level_precision = level_precision,
		probabilities = prob
	)
}

checks = list()
add = function(name, error, tolerance) {
	assign(
		"checks",
		append(
			checks,
			list(
				data.frame(
					check = name,
					error = error,
					tolerance = tolerance,
					passed = is.finite(error) && error < tolerance
				)
			)
		),
		envir = parent.env(environment())
	)
}

for (r in 1:2) {
	for (rho in c(-.8, .95, .9899)) {
		for (scale in c(.1, 1, 4)) {
			n = 4
			periods = 3
			u = array(rnorm(n * r * periods, sd = scale), c(n, r, periods))
			v = u + array(rnorm(n * r * periods, sd = .1), dim(u))
			y = array(NA_real_, c(n, n, periods))
			for (t in 1:periods) {
				z = matrix(rnorm(n * n), n)
				y[, , t] = (z + t(z)) / 2
				diag(y[, , t]) = NA
			}
			y[1, 2, 2] = y[2, 1, 2] = NA
			args = list(
				y = y,
				u = u,
				s2 = .7,
				va = .3,
				beta_precision = .4,
				rho = rho,
				sigma = .2,
				kappa = 2,
				pi = .04
			)
			before = .Random.seed
			u0 = u
			y0 = y
			a = do.call(jump_evaluate, args)
			ra = do.call(reference, args)
			stopifnot(identical(before, .Random.seed), identical(u, u0), identical(y, y0))
			args$u = v
			b = do.call(jump_evaluate, args)
			rb = do.call(reference, args)
			id = paste(r, rho, scale, sep = "_")
			add(
				paste0("density_difference_", id),
				abs((b$log_density - a$log_density) - (rb$log_density - ra$log_density)),
				1e-8
			)
			add(paste0("level_mean_", id), max(abs(a$level_mean - ra$level_mean)), 1e-10)
			add(
				paste0("level_precision_", id),
				max(abs(a$level_precision - ra$level_precision)),
				1e-10
			)
			add(
				paste0("reset_probability_", id),
				max(abs(a$probabilities - ra$probabilities)),
				1e-12
			)
			forward = min(0, b$log_density - a$log_density)
			reverse = min(0, a$log_density - b$log_density)
			add(
				paste0("detailed_balance_", id),
				abs((a$log_density + forward) - (b$log_density + reverse)),
				1e-8
			)
		}
	}
}

n = 3
periods = 2
r = 1
y = array(.2, c(n, n, periods))
for (t in 1:periods) diag(y[, , t]) = NA
u = array(c(-.2, .1, .3, -.1, .2, .4), c(n, r, periods))
shift = array(c(.08, -.03, .02, .04, .02, -.05), dim(u))
args = list(
	y = y,
	u = u,
	s2 = .7,
	va = .3,
	beta_precision = .4,
	rho = .7,
	sigma = .2,
	kappa = .5,
	pi = .2
)
grid = -100:100
logp = vapply(
	grid,
	function(k) {
		args$u = u + k * shift
		do.call(reference, args)$log_density
	},
	0.
)
p = exp(logp - max(logp))
p = p / sum(p)
stopifnot(sum(p[c(1, length(p))]) < 1e-20)
simulation_draws = 25000L
initial = sample(grid, simulation_draws, TRUE, p)
out = numeric(simulation_draws)
zs = matrix(NA_real_, simulation_draws, n + 1)
delta = matrix(NA_real_, simulation_draws, n)
expected_delta = matrix(NA_real_, simulation_draws, n)
for (j in seq_len(simulation_draws)) {
	args$u = u + initial[j] * shift
	z = do.call(jump_draw, c(args, list(shift = shift)))
	out[j] = initial[j] + if (z$accepted) z$direction else 0
	args$u = z$u
	a = do.call(reference, args)
	zs[j, ] = drop(chol(a$level_precision) %*% (c(z$beta, z$a) - a$level_mean))
	delta[j, ] = z$delta[, 2] - a$probabilities[, 2]
	expected_delta[j, ] = a$probabilities[, 2] * (1 - a$probabilities[, 2])
}

for (power in 1:4) {
	expected = sum(p * grid^power)
	sd = sqrt(sum(p * (grid^power - expected)^2))
	add(
		paste0("stationary_lattice_moment_", power),
		abs(mean(out^power) - expected) / (sd / sqrt(simulation_draws)),
		5
	)
}

for (j in 1:(n + 1)) {
	add(
		paste0("conditional_level_mean_", j),
		abs(mean(zs[, j])) * sqrt(simulation_draws),
		5
	)
	add(
		paste0("conditional_level_variance_", j),
		abs(mean(zs[, j]^2) - 1) * sqrt(simulation_draws / 2),
		5
	)
}

for (j in 1:n) {
	add(
		paste0("conditional_reset_mean_", j),
		abs(mean(delta[, j])) / sqrt(mean(expected_delta[, j]) / simulation_draws),
		5
	)
}

checks = do.call(rbind, checks)
write.csv(checks, file.path(output, "jump_validation.csv"), row.names = FALSE)
stopifnot(nrow(checks) == 105L, all(checks$passed))
# independent scalar and missing-trajectory calculations
local(
	{
		set.seed(91138)
		checks = list()
		for (left in c(0L, 1L)) {
			for (rho in c(-.8, .95)) {
				right = 5L
				sigma = .3
				kappa = 1.2
				pi = .07
				start = if (left == 0L) NULL else c(.2, -.1)
				end = c(.4, .15)
				first = if (left == 0L) 2L else left + 1L
				combinations = as.matrix(expand.grid(rep(list(0:1), right - first + 1L)))
				weight = numeric(nrow(combinations))
				last = integer(nrow(combinations))
				for (j in seq_len(nrow(combinations))) {
					mean = if (left == 0L) c(0, 0) else start
					variance = if (left == 0L) sigma^2 / (1 - rho^2) else 0
					for (t in first:right) {
						if (combinations[j, t - first + 1L]) {
							mean = c(0, 0)
							variance = kappa^2
							last[j] = t
						} else {
							mean = rho * mean
							variance = rho^2 * variance + sigma^2
						}
					}
					weight[j] = exp(
						sum(dbinom(combinations[j, ], 1, pi, log = TRUE)) + sum(dnorm(end, mean, sqrt(variance), log = TRUE))
					)
				}
				actual = revision_missing_bridge_weights(
					start,
					end,
					left,
					right,
					rho,
					sigma,
					kappa,
					qlogis(pi)
				)
				expected = vapply(
					actual$last_reset,
					function(k) sum(weight[last == k]) / sum(weight),
					0.
				)
				error = max(abs(exp(actual$log_probability) - expected))
				checks[[length(checks) + 1L]] = data.frame(
					check = paste("enumerated_missing_bridge", left, rho),
					error = error,
					passed = error < 1e-12
				)
			}
		}
		for (case in 1:3) {
			u = array(rnorm(4 * 2 * 3, sd = .3), c(4, 2, 3))
			delta = matrix(c(0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1), 4, 3)
			rho = c(-.8, .5, .98)[case]
			prior_shape = 2
			prior_scale = 1
			ss = (1 - rho^2) * sum(u[, , 1]^2)
			count = 8L
			for (t in 2:3) {
				keep = delta[, t] == 0
				ss = ss + sum((u[keep, , t] - rho * u[keep, , t - 1])^2)
				count = count + 2 * sum(keep)
			}
			shape = prior_shape + count / 2
			rate = prior_scale + ss / 2
			draw = get(".lame_sample_sigma_uv_drift_only", asNamespace("lame"))
			values = replicate(
				20000,
				1 / draw(u, u, rho, .3, TRUE, delta, NULL, prior_shape, prior_scale)^2
			)
			z = abs(mean(values) - shape / rate) / sqrt(shape / rate^2 / length(values))
			checks[[length(checks) + 1L]] = data.frame(
				check = paste("drift_variance_conditional", case),
				error = z,
				passed = z < 4.5
			)
		}
		checks = do.call(rbind, checks)
		stopifnot(nrow(checks) == 7L, all(checks$passed))
		write.csv(
			checks,
			file.path(output, "elementary_checks.csv"),
			row.names = FALSE
		)
	}
)
# restarting and continuing must reproduce the identical random path
for (indicator_every in c(0L, 5L)) {
	for (temporal_tau in c(0, .01)) {
		local(
	{
				set.seed(33917)
				y = array(
					NA_real_,
					c(4, 4, 5),
					dimnames = list(LETTERS[1:4], LETTERS[1:4], as.character(1:5))
				)
				for (t in 1:5) {
					m = matrix(rnorm(16), 4)
					y[, , t] = (m + t(m)) / 2
					diag(y[, , t]) = NA_real_
				}
				y[4, , 2:3] = y[, 4, 2:3] = NA_real_
				y[3, , 4:5] = y[, 3, 4:5] = NA_real_
				shifts = list(array(rnorm(40, sd = .01), c(4, 2, 5)))
				prior = model_priors()$controls
				scratch = tempfile("model_continuation_")
				dir.create(scratch)
				args = list(
					y = y,
					seed = 921,
					prior = prior,
					context = context,
					output_dir = scratch,
					nscan = 30,
					burn = 10,
					odens = 1,
					frozen_shifts = shifts,
					indicator_every = indicator_every,
					temporal_tau = temporal_tau
				)
				full = do.call(fit_model, c(args, list(label = "full")))
				partial = do.call(fit_model, c(args, list(label = "partial", max_sweeps = 17)))
				stopifnot(!partial$complete)
				resumed = do.call(fit_model, c(args, list(label = "partial")))
				args$nscan = 13
				segment = do.call(fit_model, c(args, list(label = "segment")))
				continued = continue_model(
					segment$revision_build$checkpoint_file,
					y,
					"continued",
					context,
					scratch,
					nscan = 17,
					burn = 0
				)
				fields = c("U_draws", "BETA", "revision_a_draws", "trace", "snap_draws", "start_vals")
				checks = data.frame(
					check = paste0("restart_", fields),
					passed = vapply(
						fields,
						function(field) identical(full[[field]], resumed[[field]]),
						TRUE
					)
				)
				checks = rbind(
					checks,
					data.frame(
						check = c("continued_state", "continued_trace"),
						passed = c(
							identical(full$start_vals, continued$start_vals),
							identical(full$trace[14:30, ], continued$trace)
						)
					)
				)
				stopifnot(nrow(checks) == 8L, all(checks$passed))
				write.csv(
					checks,
					file.path(
						output,
						if (temporal_tau > 0) {
							if (indicator_every == 0L) "temporal_continuation_checks.csv" else "combined_continuation_checks.csv"
						} else {
							if (indicator_every == 0L) "continuation_checks.csv" else "indicator_continuation_checks.csv"
						}
					),
					row.names = FALSE
				)
				unlink(scratch, recursive = TRUE)
			}
)
	}
}

local(
	{
		backend = summary_backend(root)
		softplus = function(x) pmax(x, 0) + log1p(exp(-abs(x)))
		logsum = function(x) {
			m = max(x)
			m + log(sum(exp(x - m)))
		}
		set.seed(94122)
		checks = list()
		for (n in 2:8) {
			states = as.matrix(expand.grid(rep(list(0:1), n)))
			for (scale in c(.1, 4, 100)) {
				odds = matrix(rnorm(n * 6, sd = scale), 6, n)
				actual = backend$reset_count_probabilities(odds, seq_len(n))
				counts = matrix(0, 6, n + 1)
				for (s in 1:6) {
					lp = -softplus(-odds[s, ])
					lq = -softplus(odds[s, ])
					logweights = drop(states %*% lp + (1 - states) %*% lq)
					distribution = vapply(0:n, function(k) logsum(logweights[rowSums(states) == k]), 0.)
					counts[s, ] = exp(distribution - logsum(distribution))
					expected = vapply(
						1:n,
						function(k) logsum(logweights[rowSums(states) >= k]) - logsum(logweights[rowSums(states) < k]),
						0.
					)
					error = max(abs(expected - actual$log_odds[s, ]))
					checks[[length(checks) + 1]] = data.frame(
						check = paste("enumerated", n, scale, s),
						error = error,
						passed = error < 1e-8
					)
				}
				error = max(abs(colMeans(counts) - actual$mean_count_probabilities))
				checks[[length(checks) + 1]] = data.frame(
					check = paste("enumerated_count_mean", n, scale),
					error = error,
					passed = error < 1e-10
				)
			}
		}
		for (n in c(40, 200)) {
			for (odds in c(-1000, -10, -1, 0, 1, 10, 1000)) {
				thresholds = as.integer(ceiling(c(.1, .25, .5) * n))
				actual = backend$reset_count_probabilities(matrix(odds, 1, n), thresholds)
				k = 0:n
				weights = lchoose(n, k) - k * softplus(-odds) - (n - k) * softplus(odds)
				expected = vapply(
					thresholds,
					function(threshold) logsum(weights[k >= threshold]) - logsum(weights[k < threshold]),
					0.
				)
				error = max(abs(expected - actual$log_odds))
				checks[[length(checks) + 1]] = data.frame(
					check = paste("homogeneous_binomial", n, odds),
					error = error,
					passed = error < 1e-7
				)
			}
		}
		checks = do.call(rbind, checks)
		stopifnot(nrow(checks) == 161L, all(checks$passed))
		write.csv(
			checks,
			file.path(output, "count_probability_checks.csv"),
			row.names = FALSE
		)
	}
)
local(
	{
		checks = list()
		check = function(label, error, tolerance) {
			assign(
				"checks",
				append(
					checks,
					list(
						data.frame(
							check = label,
							error,
							tolerance,
							passed = is.finite(error) && error <= tolerance
						)
					)
				),
				envir = parent.env(environment())
			)
			if (!tail(checks, 1)[[1]]$passed) stop(label, ": ", error)
		}
		# independent dense gaussian integration over an actor's path and static effect.
		dense = function(y, u, a, beta, s2, va, rho, sigma, kappa, pi, delta, actor) {
			n = dim(u)[1]
			r = dim(u)[2]
			periods = dim(u)[3]
			p = r * periods + 1L
			path_precision = matrix(0, p, p)
			b = numeric(p)
			path_precision[p, p] = 1 / va
			path_precision[1:r, 1:r] = diag((1 - rho^2) / sigma^2, r)
			for (t in 2:periods) {
				now = (t - 1) * r + 1:r
				before = (t - 2) * r + 1:r
				if (delta[actor, t]) {
					path_precision[now, now] = path_precision[now, now] + diag(1 / kappa^2, r)
				} else {
					innovation_design = matrix(0, r, p)
					innovation_design[, now] = diag(r)
					innovation_design[, before] = -rho * diag(r)
					path_precision = path_precision + crossprod(innovation_design) / sigma^2
				}
			}
			for (t in 1:periods) {
				for (j in setdiff(1:n, actor)) {
					if (is.finite(y[actor, j, t])) {
						x = numeric(p)
						x[(t - 1) * r + 1:r] = u[j, , t]
						x[p] = 1
						path_precision = path_precision + tcrossprod(x) / s2
						b = b + x * (y[actor, j, t] - beta - a[j]) / s2
					}
				}
			}
			reference_covariance = solve(path_precision)
			m = drop(reference_covariance %*% b)
			snaps = sum(delta[actor, -1])
			loge = .5 * sum(b * m) - .5 * as.numeric(determinant(path_precision, logarithm = TRUE)$modulus) +
				snaps * (log(pi) - r * log(kappa)) + (periods - 1 - snaps) * (log1p(-pi) - r * log(sigma))
			list(mean = m, covariance = reference_covariance, log_evidence = loge)
		}
		set.seed(390127)
		for (rank in c(1L, 2L, 3L)) {
			for (rep in 1:12) {
				n = 5L
				periods = 4L
				u = array(rnorm(n * rank * periods, sd = .4), c(n, rank, periods))
				a = rnorm(n, sd = .1)
				y = array(NA_real_, c(n, n, periods))
				for (t in 1:periods) {
					m = matrix(rnorm(n * n), n)
					m = (m + t(m)) / 2
					diag(m) = NA
					y[, , t] = m
				}
				if (rep %% 2 == 0) y[1, 2, 2] = y[2, 1, 2] = NA
				if (rep %% 3 == 0) y[1, , 3] = y[, 1, 3] = NA
				if (rep %% 4 == 0) y[1, , ] = y[, 1, ] = NA
				delta = matrix(rbinom(n * periods, 1, .3), n, periods)
				delta[, 1] = 0
				args = list(
					y = y,
					u = u,
					a = a,
					beta = .2,
					s2 = .3,
					va = .1,
					rho = if (rep %% 2 == 0) -.65 else .82,
					sigma = .15,
					kappa = .7,
					pi = .12,
					delta = delta
				)
				actual = do.call(indicator_evaluate, c(args, list(actor = 0L)))
				expected = do.call(dense, c(args, list(actor = 1L)))
				check(
					paste("dense_evidence", rank, rep),
					abs(actual$log_evidence - expected$log_evidence),
					1e-8
				)
				check(
					paste("dense_mean", rank, rep),
					max(abs(c(actual$path_mean, actual$a_mean) - expected$mean)),
					1e-9
				)
				old = do.call(
					actor_path_conditional,
					args[setdiff(names(args), "pi")] |> c(list(actor = 0L))
				)
				check(paste("old_mean", rank, rep), max(abs(old$mean - expected$mean)), 1e-9)
			}
		}
		rank = 2L
		n = 3L
		periods = 4L
		u = array(rnorm(n * rank * periods, sd = .25), c(n, rank, periods))
		a = c(.1, -.1, .05)
		y = array(NA_real_, c(n, n, periods))
		for (t in 1:periods) {
			y[, , t] = .2 + outer(a, a, "+") + tcrossprod(u[, , t])
			diag(y[, , t]) = NA
		}
		delta = matrix(0, n, periods)
		args = list(
			y = y,
			u = u,
			a = a,
			beta = .2,
			s2 = .2,
			va = .2,
			rho = .7,
			sigma = .2,
			kappa = .65,
			pi = .25,
			delta = delta
		)
		patterns = as.matrix(expand.grid(rep(list(0:1), periods - 1)))
		exact = lapply(
			1:nrow(patterns),
			function(k) {
				x = args
				x$delta[1, -1] = patterns[k, ]
				do.call(dense, c(x, list(actor = 1L)))
			}
		)
		lw = vapply(exact, "[[", 0, "log_evidence")
		weight = exp(lw - max(lw))
		weight = weight / sum(weight)
		means = do.call(rbind, lapply(exact, "[[", "mean"))
		mean_exact = colSums(means * weight)
		second_exact = Reduce(
			"+",
			Map(function(e, w) w * (e$covariance + tcrossprod(e$mean)), exact, weight)
		)
		simulation_draws = 24000L
		draws = matrix(NA_real_, simulation_draws, rank * periods + 1L)
		states = integer(simulation_draws)
		for (k in 1:simulation_draws) {
			z = do.call(indicator_pair_sweep, c(args, list(selected_actor = 0L)))
			args$u = z$u
			args$a = as.vector(z$a)
			args$delta = z$delta
			draws[k, ] = c(z$u[1, , ], z$a[1])
			states[k] = 1L + sum(z$delta[1, -1] * 2^(0:(periods - 2)))
		}
		keep = 1001:simulation_draws
		batch = rep(1:100, length.out = length(keep))
		bmeans = rowsum(draws[keep, ], batch) / as.vector(table(batch))
		se = apply(bmeans, 2, sd) / 10
		check(
			"mixture_path_means",
			max(abs(colMeans(draws[keep, ]) - mean_exact) / pmax(se, 1e-10)),
			6
		)
		freq = tabulate(states[keep], nrow(patterns)) / length(keep)
		check("mixture_pattern_probabilities", max(abs(freq - weight)), .022)
		check(
			"mixture_second_moments",
			max(abs(crossprod(draws[keep, ]) / length(keep) - second_exact)),
			.025
		)
		set.seed(713)
		one = do.call(indicator_pair_sweep, args)
		rng = .Random.seed
		set.seed(713)
		two = do.call(indicator_pair_sweep, args)
		check(
			"rng_reproduction",
			as.numeric(!identical(one, two) || !identical(rng, .Random.seed)),
			0
		)
		before = digest::digest(args)
		set.seed(913)
		invisible(do.call(indicator_pair_sweep, args))
		check(
			"input_preservation",
			as.numeric(!identical(before, digest::digest(args))),
			0
		)
		for (field in c("s2", "va", "sigma", "kappa", "pi")) {
			bad = args
			bad[[field]] = 0
			check(
				paste("invalid_scale", field),
				as.numeric(
					!inherits(try(do.call(indicator_pair_sweep, bad), silent = TRUE), "try-error")
				),
				0
			)
		}
		check(
			"invalid_actor",
			as.numeric(
				!inherits(
					try(
						do.call(indicator_pair_sweep, c(args, list(selected_actor = -2L))),
						silent = TRUE
					),
					"try-error"
				)
			),
			0
		)
		write.csv(
			do.call(rbind, checks),
			file.path(output, "indicator_checks.csv"),
			row.names = FALSE
		)
	}
)

local(
	{
		logsum = function(x) {
			m = max(x)
			m + log(sum(exp(x - m)))
		}
		lognormal = function(x, covariance) {
			root = chol(covariance)
			z = forwardsolve(t(root), x)
			-.5 * sum(z * z) - sum(log(diag(root))) - length(x) * log(2 * pi) / 2
		}
		reference = function(y, u, s2, va, beta_precision, rho, sigma, kappa, pi) {
			n = dim(u)[1]
			rank = dim(u)[2]
			tt = dim(u)[3]
			xx = list()
			rr = numeric()
			for (t in seq_len(tt)) {
				for (j in 2:n) {
					for (i in seq_len(j - 1)) {
						if (is.finite(y[i, j, t])) {
							row = numeric(n + 1)
							row[c(1, i + 1, j + 1)] = 1
							xx[[length(xx) + 1]] = row
							rr = c(rr, y[i, j, t] - sum(u[i, , t] * u[j, , t]))
						}
					}
				}
			}
			ans = 0
			if (length(rr)) {
				design = do.call(rbind, xx)
				cov = diag(s2, length(rr)) + design %*% diag(c(1 / beta_precision, rep(va, n))) %*% t(design)
				ans = lognormal(rr, cov)
			}
			patterns = if (tt > 1) as.matrix(expand.grid(rep(list(0:1), tt - 1))) else matrix(numeric(), 1, 0)
			for (i in seq_len(n)) {
				terms = numeric(nrow(patterns))
				for (k in seq_len(nrow(patterns))) {
					cov = matrix(0, tt, tt)
					cov[1, 1] = sigma^2 / (1 - rho^2)
					if (tt > 1) {
						for (t in 2:tt) {
							if (patterns[k, t - 1] == 1) {
								cov[t, t] = kappa^2
							} else {
								cov[t, seq_len(t - 1)] = rho * cov[t - 1, seq_len(t - 1)]
								cov[seq_len(t - 1), t] = cov[t, seq_len(t - 1)]
								cov[t, t] = rho^2 * cov[t - 1, t - 1] + sigma^2
							}
						}
					}
					terms[k] = sum(patterns[k, ]) * log(pi) + (tt - 1 - sum(patterns[k, ])) * log1p(-pi)
					for (axis in seq_len(rank)) terms[k] = terms[k] + lognormal(u[i, axis, ], cov)
				}
				ans = ans + logsum(terms)
			}
			ans
		}
		checks = list()
		check = function(label, error, tolerance) {
			assign(
				"checks",
				append(
					checks,
					list(
						data.frame(
							check = label,
							error,
							tolerance,
							passed = is.finite(error) && error <= tolerance
						)
					)
				),
				envir = parent.env(environment())
			)
			if (!tail(checks, 1)[[1]]$passed) stop(label, ": ", error)
		}
		set.seed(390128)
		for (rank in 1:3) {
			for (periods in c(1L, 3L)) {
				for (rep in 1:6) {
					n = 3L
					u = array(rnorm(n * rank * periods, sd = .3), c(n, rank, periods))
					grand = apply(u, 2, mean)
					min_square = sum(grand^2)
					increment = if (rep %% 2 == 0) -.3 * min_square else .2 * min_square
					z = temporal_map(u, increment)
					inverse = temporal_map(z$u, -increment)
					check(paste("inverse", rank, periods, rep), max(abs(inverse$u - u)), 1e-12)
					check(
						paste("inverse_jacobian", rank, periods, rep),
						abs(z$log_jacobian + inverse$log_jacobian),
						1e-12
					)
					state = c(u, increment)
					mapping = function(v) {
						x = temporal_map(array(v[-length(v)], dim(u)), v[length(v)])
						c(x$u, -tail(v, 1))
					}
					step = min(1e-6, min_square / 1000)
					jacobian_matrix = sapply(
						seq_along(state),
						function(k) {
							a = b = state
							a[k] = a[k] + step
							b[k] = b[k] - step
							(mapping(a) - mapping(b)) / (2 * step)
						}
					)
					check(
						paste("dense_jacobian", rank, periods, rep),
						abs(
							as.numeric(determinant(jacobian_matrix, logarithm = TRUE)$modulus) - z$log_jacobian
						),
						2e-7
					)
					y = array(NA_real_, c(n, n, periods))
					for (t in 1:periods) {
						m = matrix(rnorm(n * n), n)
						y[, , t] = (m + t(m)) / 2
						diag(y[, , t]) = NA
					}
					if (rep %% 3 == 0) y[] = NA
					if (rep %% 3 == 1) y[1, 2, 1] = y[2, 1, 1] = NA
					p = list(
						s2 = .2,
						va = .3,
						beta_precision = .7,
						rho = if (rep %% 2 == 0) -.6 else .7,
						sigma = .4,
						kappa = 1.3,
						pi = .12
					)
					old = do.call(reference, c(list(y = y, u = u), p))
					nextvalue = do.call(reference, c(list(y = y, u = z$u), p))
					actual = do.call(temporal_evaluate, c(list(y = y, u = u, increment = increment), p))
					check(
						paste("independent_marginal_ratio", rank, periods, rep),
						abs(actual$log_ratio - (nextvalue - old + z$log_jacobian)),
						1e-8
					)
				}
			}
		}
		statistics = function(z) {
			c(
				mean = mean(z$u),
				energy = sum(z$u^2),
				fourth = sum(z$u^2)^2,
				beta = z$beta,
				beta_squared = z$beta^2,
				additive = sum(z$a^2),
				resets = sum(z$delta)
			)
		}
		simulation_draws = 20000L
		for (rank in 1:3) {
			n = 3L
			periods = 3L
			p = list(
				s2 = .2,
				va = .3,
				beta_precision = .7,
				rho = .7,
				sigma = .4,
				kappa = 1.3,
				pi = .12
			)
			y = array(NA_real_, c(n, n, periods))
			difference = matrix(0, simulation_draws, 7)
			for (k in 1:simulation_draws) {
				u = array(0, c(n, rank, periods))
				u[, , 1] = rnorm(n * rank, sd = p$sigma / sqrt(1 - p$rho^2))
				delta = matrix(0, n, periods)
				for (t in 2:periods) {
					delta[, t] = rbinom(n, 1, p$pi)
					u[, , t] = p$rho * u[, , t - 1] + matrix(rnorm(n * rank, sd = p$sigma), n, rank)
					for (i in which(delta[, t] == 1)) u[i, , t] = rnorm(rank, sd = p$kappa)
				}
				old = list(
					u = u,
					delta = delta,
					beta = rnorm(1, sd = 1 / sqrt(p$beta_precision)),
					a = rnorm(n, sd = sqrt(p$va))
				)
				z = do.call(temporal_draw, c(list(y = y, u = u, tau = .1, moves = 2L), p))
				difference[k, ] = statistics(z) - statistics(old)
			}
			z = colMeans(difference) / (apply(difference, 2, sd) / sqrt(simulation_draws))
			check(paste("independent_stationarity", rank), max(abs(z)), 5)
		}
		before = digest::digest(u)
		set.seed(1213)
		x = do.call(temporal_draw, c(list(y = y, u = u, tau = .1, moves = 2L), p))
		rng = .Random.seed
		set.seed(1213)
		z = do.call(temporal_draw, c(list(y = y, u = u, tau = .1, moves = 2L), p))
		check(
			"rng_and_input_preservation",
			as.numeric(
				!identical(x, z) || !identical(rng, .Random.seed) || !identical(before, digest::digest(u))
			),
			0
		)
		check("invalid_temporal_rejected", as.numeric(temporal_map(u, -1e6)$valid), 0)
		write.csv(
			do.call(rbind, checks),
			file.path(output, "temporal_checks.csv"),
			row.names = FALSE
		)
	}
)
# garbage collection must not change a draw or its random-number state.
local(
	{
		cpp = paste(readLines(file.path(root, "code", "sampler.cpp")), collapse = "\n")
		header = "#include <RcppArmadillo.h>"
		helper = '\nstatic bool lifetime_stress = false;\n// [[Rcpp::export]]\nvoid set_lifetime_stress(bool active){lifetime_stress=active;}\nvoid lifetime_collect(unsigned n){if(lifetime_stress){Rcpp::Function collect=Rcpp::Environment::base_env()["gc"];collect();Rcpp::NumericVector reused(n,123.0);}}\n'
		stopifnot(grepl(header, cpp, fixed = TRUE))
		cpp = sub(header, paste0(header, helper), cpp, fixed = TRUE)
		old = c(
			"for(arma::uword t = 0; t < u.n_slices; ++t) {\n    Rcpp::List result",
			"for(arma::uword t = 1; t < u.n_slices; ++t) {\n    Rcpp::List result",
			"for(int k=0;k<moves;++k){\n  Rcpp::List mapping"
		)
		new = c(
			"for(arma::uword t = 0; t < u.n_slices; ++t) {\n    lifetime_collect(u.n_elem);\n    Rcpp::List result",
			"for(arma::uword t = 1; t < u.n_slices; ++t) {\n    lifetime_collect(u.n_elem);\n    Rcpp::List result",
			"for(int k=0;k<moves;++k){\n  lifetime_collect(selected.n_elem);\n  Rcpp::List mapping"
		)
		for (i in seq_along(old)) {
			matches = gregexpr(old[i], cpp, fixed = TRUE)[[1]]
			stopifnot(length(matches) == 1L, matches[1] > 0)
			cpp = sub(old[i], new[i], cpp, fixed = TRUE)
		}
		probe = new.env(parent = baseenv())
		Rcpp::sourceCpp(
			code = cpp,
			env = probe,
			cacheDir = file.path(root, "code", ".compiled", "lifetime_validation"),
			showOutput = FALSE
		)
		set.seed(21)
		u = array(rnorm(8 * 2 * 5, sd = .1), c(8, 2, 5))
		delta = matrix(0, 8, 5)
		y = array(0, c(8, 8, 5))
		for (t in 1:5) diag(y[, , t]) = NA
		checks = list()
		for (kind in c("temporal", "collapsed", "suffix")) {
			draw = if (kind == "temporal") {
				function() probe$temporal_draw(y, u + .5, .01, 15L, 1, 1, 1, .7, .2, 1, .1)
			} else if (kind == "collapsed") {
				function() probe$revision_collapsed_sweep_cpp(u, delta, .7, .2, 1, .1, 1L)
			} else {
				function() probe$revision_suffix_sweep_cpp(u, delta, .7, .2, 1, rep(qlogis(.1), 5), 1L)
			}
			before = digest::digest(list(u, y, delta))
			probe$set_lifetime_stress(FALSE)
			set.seed(502)
			baseline = draw()
			rng = .Random.seed
			probe$set_lifetime_stress(TRUE)
			set.seed(502)
			stress = draw()
			if (kind == "temporal") stopifnot(baseline$accepted > 0)
			checks[[kind]] = data.frame(
				check = paste(kind, c("draw", "rng", "input"), sep = "_"),
				passed = c(
					identical(baseline, stress),
					identical(rng, .Random.seed),
					identical(before, digest::digest(list(u, y, delta)))
				)
			)
		}
		checks = do.call(rbind, checks)
		write.csv(
			checks,
			file.path(output, "array_lifetime_checks.csv"),
			row.names = FALSE
		)
		stopifnot(all(checks$passed))
	}
)
files = list.files(output, "\\.csv$", full.names = TRUE)
checks = lapply(files, read.csv)
checks = checks[vapply(checks, function(x) "passed" %in% names(x), TRUE)]
stopifnot(
	length(checks) >= 8L,
	all(vapply(checks, function(x) nrow(x) > 0L && all(x$passed), TRUE))
)
source_files = file.path(
	root,
	"code",
	c(
		"00_sampler.R",
		"sampler.cpp",
		"10_models.R",
		"21_validate.R",
		"vendor/lame_1.3.5_americium_full_fix.tar.gz"
	)
)
jsonlite::write_json(
	list(
		passed = TRUE,
		checks = sum(vapply(checks, nrow, 1L)),
		checked_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
		evidence = setNames(
			lapply(files, function(p) digest::digest(file = p, algo = "sha256")),
			substring(files, nchar(root) + 2L)
		),
		sources = setNames(
			lapply(source_files, function(p) digest::digest(file = p, algo = "sha256")),
			basename(source_files)
		),
		R = R.version.string,
		packages = setNames(
			lapply(
				c(
					"data.table",
					"digest",
					"ggplot2",
					"netify",
					"jsonlite",
					"posterior",
					"Rcpp",
					"RcppArmadillo",
					"rvest",
					"xml2",
					"lame"
				),
				function(package) as.character(utils::packageVersion(package))
			),
			c(
				"data.table",
				"digest",
				"ggplot2",
				"netify",
				"jsonlite",
				"posterior",
				"Rcpp",
				"RcppArmadillo",
				"rvest",
				"xml2",
				"lame"
			)
		),
		interpretation = "Implementation checks; empirical convergence is assessed separately."
	),
	file.path(root, "results", "checks", "model_validation.json"),
	pretty = TRUE,
	auto_unbox = TRUE
)
