#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

namespace gaussian {




Rcpp::List revision_gaussian_prepare_cpp(const arma::cube& y,
    const arma::cube& u, const arma::vec& beta, const arma::mat& a,
    double s2) {
  const int n = y.n_rows;
  const int periods = y.n_slices;
  arma::cube z(n, n, periods, arma::fill::zeros);
  arma::mat row_sums(n, periods, arma::fill::zeros);
  arma::vec dyad_sums(periods, arma::fill::zeros);
  long double ssr = 0;
  const double sd = std::sqrt(s2);
  for(int t = 0; t < periods; ++t) {
    for(int j = 0; j < n; ++j) {
      z(j, j, t) = beta(t) + 2 * a(j, t) + arma::dot(u.slice(t).row(j), u.slice(t).row(j));
      for(int i = 0; i < j; ++i) {
        double latent = arma::dot(u.slice(t).row(i), u.slice(t).row(j));
        double eta = beta(t) + (a(i, t) + a(j, t)) + latent;
        double value = y(i, j, t);
        if(!std::isfinite(value)) value = y(j, i, t);
        if(!std::isfinite(value)) value = R::rnorm(eta, sd);
        z(i, j, t) = z(j, i, t) = value;
        double residual = value - latent;
        row_sums(i, t) += residual;
        row_sums(j, t) += residual;
        dyad_sums(t) += residual;
        double error = value - eta;
        ssr += static_cast<long double>(error * error);
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("z") = z,
    Rcpp::Named("row_sums") = row_sums,
    Rcpp::Named("dyad_sums") = dyad_sums,
    Rcpp::Named("ssr") = static_cast<double>(ssr),
    Rcpp::Named("n_dyads") = n * (n - 1) / 2 * periods);
}


arma::cube revision_gaussian_residual_cpp(const arma::cube& z,
    const arma::vec& beta, const arma::mat& a) {
  arma::cube residual(z.n_rows, z.n_cols, z.n_slices);
  for(arma::uword t = 0; t < z.n_slices; ++t) {
    for(arma::uword j = 0; j < z.n_cols; ++j) {
      for(arma::uword i = 0; i < z.n_rows; ++i) {
        residual(i, j, t) = (z(i, j, t) - beta(t)) - (a(i, t) + a(j, t));
      }
    }
  }
  return residual;
}


double revision_gaussian_loglik_cpp(const arma::cube& y,
    const arma::cube& u, const arma::vec& beta, const arma::mat& a,
    double s2) {
  long double result = 0;
  const double sd = std::sqrt(s2);
  for(arma::uword t = 0; t < y.n_slices; ++t) {
    for(arma::uword j = 0; j < y.n_cols; ++j) {
      for(arma::uword i = 0; i < j; ++i) {
        double value = y(i, j, t);
        if(!std::isfinite(value)) value = y(j, i, t);
        if(std::isfinite(value)) {
          double eta = beta(t) + (a(i, t) + a(j, t)) +
            arma::dot(u.slice(t).row(i), u.slice(t).row(j));
          result += static_cast<long double>(R::dnorm4(value, eta, sd, 1));
        }
      }
    }
  }
  return static_cast<double>(result);
}

}
// [[Rcpp::export]]
Rcpp::List revision_gaussian_prepare_cpp(const arma::cube& y,
    const arma::cube& u, const arma::vec& beta, const arma::mat& a,
    double s2) {
  return gaussian::revision_gaussian_prepare_cpp(y, u, beta, a, s2);
}

// [[Rcpp::export]]
arma::cube revision_gaussian_residual_cpp(const arma::cube& z,
    const arma::vec& beta, const arma::mat& a) {
  return gaussian::revision_gaussian_residual_cpp(z, beta, a);
}

// [[Rcpp::export]]
double revision_gaussian_loglik_cpp(const arma::cube& y,
    const arma::cube& u, const arma::vec& beta, const arma::mat& a,
    double s2) {
  return gaussian::revision_gaussian_loglik_cpp(y, u, beta, a, s2);
}

namespace rotation {



namespace {
  const double rotation_pi = std::acos(-1.0);

  double von_mises_offset(double concentration) {
    if(!std::isfinite(concentration) || concentration < 0) Rcpp::stop("Invalid angular concentration.");
    for(int attempt = 0; attempt < 1000000; ++attempt) {
      double angle;
      double log_acceptance;
      if(concentration < 1) {
        angle = R::runif(-rotation_pi, rotation_pi);
        log_acceptance = -2 * (concentration * std::pow(std::sin(angle / 2), 2));
      } else {
        angle = R::rnorm(0, rotation_pi / (2 * std::sqrt(concentration)));
        if(std::abs(angle) >= rotation_pi) continue;
        log_acceptance = 2 * (concentration *
          (std::pow(angle / rotation_pi, 2) - std::pow(std::sin(angle / 2), 2)));
      }
      if(std::log(R::runif(0, 1)) <= log_acceptance) return angle;
    }
    Rcpp::stop("Angular rejection sampler exhausted its iteration limit.");
    return NA_REAL;
  }

  double log_bessel_scaled(double concentration) {
    if(!std::isfinite(concentration) || concentration < 0) Rcpp::stop("Invalid Bessel concentration.");
    if(concentration >= 10000) {
      const double inverse = 1 / concentration;
      const double correction = inverse * (1.0 / 8 + inverse * (9.0 / 128 +
        inverse * (225.0 / 3072 + inverse * (11025.0 / 98304))));
      return -0.5 * (std::log(2 * rotation_pi) + std::log(concentration)) + std::log1p(correction);
    }
    const double value = R::bessel_i(concentration, 0, 2);
    if(!std::isfinite(value) || value <= 0) Rcpp::stop("Bessel normalizer is outside numerical support.");
    return std::log(value);
  }

  struct OrthogonalDraw {
    arma::mat q;
    double angle;
    bool reflected;
    double probability_reflection;
  };

  OrthogonalDraw orthogonal_draw(const arma::mat& m) {
    const double a0 = m(0, 0) + m(1, 1);
    const double b0 = m(1, 0) - m(0, 1);
    const double a1 = m(0, 0) - m(1, 1);
    const double b1 = m(1, 0) + m(0, 1);
    const double k0 = std::hypot(a0, b0);
    const double k1 = std::hypot(a1, b1);
    const double log_odds = (k1 - k0) + log_bessel_scaled(k1) - log_bessel_scaled(k0);
    const double probability = log_odds >= 0 ? 1 / (1 + std::exp(-log_odds)) :
      std::exp(log_odds) / (1 + std::exp(log_odds));
    const bool reflected = R::runif(0, 1) < probability;
    const double angle = std::atan2(reflected ? b1 : b0, reflected ? a1 : a0) +
      von_mises_offset(reflected ? k1 : k0);
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    arma::mat q(2, 2);
    q(0, 0) = cosine;
    q(1, 0) = sine;
    q(0, 1) = reflected ? sine : -sine;
    q(1, 1) = reflected ? -cosine : cosine;
    return OrthogonalDraw{q, angle, reflected, probability};
  }
}


Rcpp::List revision_o2_rank_draw_cpp(const arma::mat& m) {
  if(m.n_rows != 2 || m.n_cols != 2 || !m.is_finite()) Rcpp::stop("Expected a finite two-dimensional orbit matrix.");
  OrthogonalDraw result = orthogonal_draw(m);
  return Rcpp::List::create(Rcpp::Named("q") = result.q,
    Rcpp::Named("angle") = result.angle, Rcpp::Named("reflected") = result.reflected,
    Rcpp::Named("probability_reflection") = result.probability_reflection);
}


arma::vec revision_von_mises_rank_cpp(int draws, double concentration) {
  if(draws < 1) Rcpp::stop("A positive draw count is required.");
  arma::vec result(draws);
  for(int draw = 0; draw < draws; ++draw) result(draw) = von_mises_offset(concentration);
  return result;
}


arma::cube revision_period_rotation_stable_cpp(const arma::cube& input, const arma::mat& delta,
    double rho, double sigma) {
  arma::cube u(input);
  if(u.n_cols < 1 || u.n_cols > 3 || !u.is_finite() || !delta.is_finite() || delta.n_rows != u.n_rows ||
      delta.n_cols != u.n_slices || !std::isfinite(rho) || !std::isfinite(sigma) || sigma <= 0) {
    Rcpp::stop("Invalid rank-one to rank-three period rotation state.");
  }
  for(arma::uword t = 0; t < u.n_slices; ++t) {
    arma::mat neighbors(u.n_rows, u.n_cols, arma::fill::zeros);
    for(arma::uword actor = 0; actor < u.n_rows; ++actor) {
      if(t > 0 && delta(actor, t) == 0) neighbors.row(actor) += u.slice(t - 1).row(actor);
      if(t + 1 < u.n_slices && delta(actor, t + 1) == 0) neighbors.row(actor) += u.slice(t + 1).row(actor);
    }
    if(u.n_cols == 1) {
      const double m = (rho / (sigma * sigma)) * arma::accu(u.slice(t) % neighbors);
      const double log_odds = 2 * m;
      const double probability_positive = log_odds >= 0 ? 1 / (1 + std::exp(-log_odds)) :
        std::exp(log_odds) / (1 + std::exp(log_odds));
      if(R::runif(0, 1) >= probability_positive) u.slice(t) *= -1;
    } else {
      for(arma::uword left = 0; left + 1 < u.n_cols; ++left) {
        for(arma::uword right = left + 1; right < u.n_cols; ++right) {
          arma::uvec columns = {left, right};
          const arma::mat positions = u.slice(t).cols(columns);
          const arma::mat m = (rho / (sigma * sigma)) * positions.t() * neighbors.cols(columns);
          const OrthogonalDraw result = orthogonal_draw(m);
          u.slice(t).cols(columns) = positions * result.q;
        }
      }
    }
  }
  return u;
}


double revision_log_bessel_scaled_cpp(double concentration) {
  return log_bessel_scaled(concentration);
}

}
// [[Rcpp::export]]
Rcpp::List revision_o2_rank_draw_cpp(const arma::mat& m) {
  return rotation::revision_o2_rank_draw_cpp(m);
}

// [[Rcpp::export]]
arma::vec revision_von_mises_rank_cpp(int draws, double concentration) {
  return rotation::revision_von_mises_rank_cpp(draws, concentration);
}

// [[Rcpp::export]]
arma::cube revision_period_rotation_stable_cpp(const arma::cube& input, const arma::mat& delta,
    double rho, double sigma) {
  return rotation::revision_period_rotation_stable_cpp(input, delta, rho, sigma);
}

// [[Rcpp::export]]
double revision_log_bessel_scaled_cpp(double concentration) {
  return rotation::revision_log_bessel_scaled_cpp(concentration);
}

namespace collapsed {



namespace {
  const double circle_pi = std::acos(-1.0);

  double log_sum(double a, double b) {
    const double large = std::max(a, b);
    return large + std::log1p(std::exp(std::min(a, b) - large));
  }

  double probability(double log_odds) {
    return log_odds >= 0 ? 1 / (1 + std::exp(-log_odds)) :
      std::exp(log_odds) / (1 + std::exp(log_odds));
  }

  struct Transition {
    double drift;
    double reset;
  };

  Transition transition(const arma::rowvec& now, const arma::rowvec& before,
      double rho, double sigma, double kappa, double rate) {
    const double rank = now.n_elem;
    const arma::rowvec innovation = (now - rho * before) / sigma;
    const arma::rowvec standardized = now / kappa;
    const double drift = std::log1p(-rate) - rank * std::log(sigma) -
      rank * std::log(2 * circle_pi) / 2 - arma::dot(innovation, innovation) / 2;
    const double reset = std::log(rate) - rank * std::log(kappa) -
      rank * std::log(2 * circle_pi) / 2 - arma::dot(standardized, standardized) / 2;
    if(!std::isfinite(drift) || !std::isfinite(reset)) Rcpp::stop("Nonfinite component log density.");
    return Transition{drift, reset};
  }

  void check_state(const arma::cube& u, double rho, double sigma, double kappa, double rate) {
    if(u.n_cols != 2 || u.n_rows < 1 || u.n_slices < 2 || !u.is_finite() ||
        !std::isfinite(rho) || std::abs(rho) >= 1 || !std::isfinite(sigma) || sigma <= 0 ||
        !std::isfinite(kappa) || kappa <= 0 || !std::isfinite(rate) || rate <= 0 || rate >= 1) {
      Rcpp::stop("Invalid collapsed orientation state.");
    }
  }

  double orbit_log(const arma::cube& u, arma::uword period, const arma::mat& proposed,
      double rho, double sigma, double kappa, double rate) {
    double value = 0;
    for(arma::uword i = 0; i < u.n_rows; ++i) {
      if(period > 0) {
        const Transition z = transition(proposed.row(i), u.slice(period - 1).row(i), rho, sigma, kappa, rate);
        value += log_sum(z.drift, z.reset);
      }
      if(period + 1 < u.n_slices) {
        const Transition z = transition(u.slice(period + 1).row(i), proposed.row(i), rho, sigma, kappa, rate);
        value += log_sum(z.drift, z.reset);
      }
    }
    return value;
  }

  arma::mat orthogonal(double angle, bool reflected) {
    arma::mat q(2, 2);
    q(0, 0) = std::cos(angle);
    q(1, 0) = std::sin(angle);
    q(0, 1) = reflected ? std::sin(angle) : -std::sin(angle);
    q(1, 1) = reflected ? -std::cos(angle) : std::cos(angle);
    return q;
  }

  void redraw(arma::mat& delta, const arma::cube& u, arma::uword period,
      double rho, double sigma, double kappa, double rate) {
    for(arma::uword i = 0; i < u.n_rows; ++i) {
      for(arma::uword s = std::max<arma::uword>(1, period); s <= std::min(period + 1, u.n_slices - 1); ++s) {
        const Transition z = transition(u.slice(s).row(i), u.slice(s - 1).row(i), rho, sigma, kappa, rate);
        delta(i, s) = R::runif(0, 1) < probability(z.reset - z.drift);
      }
    }
    delta.col(0).zeros();
  }
}


double revision_collapsed_orbit_log_cpp(const arma::cube& u, int period, const arma::mat& q,
    double rho, double sigma, double kappa, double rate) {
  check_state(u, rho, sigma, kappa, rate);
  if(period < 1 || period > static_cast<int>(u.n_slices) || q.n_rows != 2 || q.n_cols != 2 ||
      !q.is_finite() || arma::norm(q * q.t() - arma::eye(2, 2), "fro") > 1e-10) {
    Rcpp::stop("Expected a valid period and orthogonal matrix.");
  }
  return orbit_log(u, period - 1, u.slice(period - 1) * q, rho, sigma, kappa, rate);
}


Rcpp::List revision_collapsed_period_cpp(const arma::cube& input, const arma::mat& input_delta,
    int period, double rho, double sigma, double kappa, double rate, int haar_moves = 1) {
  check_state(input, rho, sigma, kappa, rate);
  if(period < 1 || period > static_cast<int>(input.n_slices) || haar_moves < 0 ||
      input_delta.n_rows != input.n_rows || input_delta.n_cols != input.n_slices ||
      !input_delta.is_finite() || arma::any(arma::vectorise((input_delta != 0) % (input_delta != 1)))) {
    Rcpp::stop("Invalid indicator array or period.");
  }
  arma::cube u(input);
  arma::mat delta(input_delta);
  const arma::uword t = period - 1;
  arma::vec counts(4, arma::fill::zeros);
  double current = orbit_log(u, t, u.slice(t), rho, sigma, kappa, rate);
  for(int move = 0; move < 2 + haar_moves; ++move) {
    arma::mat q;
    if(move == 0) q = -arma::eye(2, 2);
    else if(move == 1) { q = arma::eye(2, 2); q(0, 0) = -1; }
    else q = orthogonal(R::runif(-circle_pi, circle_pi), R::runif(0, 1) < .5);
    const arma::mat proposed = u.slice(t) * q;
    const double next = orbit_log(u, t, proposed, rho, sigma, kappa, rate);
    const double odds = next - current;
    const bool take = move < 2 ? R::runif(0, 1) < probability(odds) :
      std::log(R::runif(0, 1)) < std::min(0.0, odds);
    counts(move < 2 ? 0 : 2) += 1;
    if(take) {
      u.slice(t) = proposed;
      current = next;
      counts(move < 2 ? 1 : 3) += 1;
    }
  }
  redraw(delta, u, t, rho, sigma, kappa, rate);
  return Rcpp::List::create(Rcpp::Named("u") = u, Rcpp::Named("delta") = delta,
    Rcpp::Named("counts") = counts);
}


Rcpp::List revision_collapsed_sweep_cpp(const arma::cube& input, const arma::mat& input_delta,
    double rho, double sigma, double kappa, double rate, int haar_moves = 1) {
  arma::cube u(input);
  arma::mat delta(input_delta);
  arma::vec counts(4, arma::fill::zeros);
  for(arma::uword t = 0; t < u.n_slices; ++t) {
    Rcpp::List result = revision_collapsed_period_cpp(u, delta, t + 1, rho, sigma, kappa, rate, haar_moves);
    // Cube conversion borrows R memory; keep an owned copy after result leaves scope.
    const arma::cube extracted_u = Rcpp::as<arma::cube>(result["u"]);
    u = extracted_u;
    delta = Rcpp::as<arma::mat>(result["delta"]);
    counts += Rcpp::as<arma::vec>(result["counts"]);
  }
  return Rcpp::List::create(Rcpp::Named("u") = u, Rcpp::Named("delta") = delta,
    Rcpp::Named("counts") = counts);
}

}
// [[Rcpp::export]]
double revision_collapsed_orbit_log_cpp(const arma::cube& u, int period, const arma::mat& q,
    double rho, double sigma, double kappa, double rate) {
  return collapsed::revision_collapsed_orbit_log_cpp(u, period, q, rho, sigma, kappa, rate);
}

// [[Rcpp::export]]
Rcpp::List revision_collapsed_period_cpp(const arma::cube& input, const arma::mat& input_delta,
    int period, double rho, double sigma, double kappa, double rate, int haar_moves = 1) {
  return collapsed::revision_collapsed_period_cpp(input, input_delta, period, rho, sigma, kappa, rate, haar_moves);
}

// [[Rcpp::export]]
Rcpp::List revision_collapsed_sweep_cpp(const arma::cube& input, const arma::mat& input_delta,
    double rho, double sigma, double kappa, double rate, int haar_moves = 1) {
  return collapsed::revision_collapsed_sweep_cpp(input, input_delta, rho, sigma, kappa, rate, haar_moves);
}

namespace suffix {



namespace {
  const double circle_pi = std::acos(-1.0);
  double softplus(double x) { return std::max(0.0, x) + std::log1p(std::exp(-std::abs(x))); }
  double log_sum(double a, double b) { return std::max(a, b) + std::log1p(std::exp(-std::abs(a - b))); }
  double probability(double x) { return std::exp(-softplus(-x)); }
  struct Terms { double drift; double reset; };
  Terms terms(const arma::rowvec& now, const arma::rowvec& before,
      double rho, double sigma, double kappa, double odds) {
    const arma::rowvec innovation = (now - rho * before) / sigma;
    const arma::rowvec standardized = now / kappa;
    const double d = -softplus(odds) - now.n_elem * (std::log(sigma) + std::log(2 * circle_pi) / 2) - arma::dot(innovation, innovation) / 2;
    const double r = -softplus(-odds) - now.n_elem * (std::log(kappa) + std::log(2 * circle_pi) / 2) - arma::dot(standardized, standardized) / 2;
    if(!std::isfinite(d) || !std::isfinite(r)) Rcpp::stop("Nonfinite suffix component density.");
    return Terms{d, r};
  }
  void check(const arma::cube& u, int period, double rho, double sigma, double kappa, const arma::vec& odds) {
    if(u.n_cols != 2 || u.n_rows < 1 || u.n_slices < 2 || !u.is_finite() ||
        period < 2 || period > static_cast<int>(u.n_slices) || !std::isfinite(rho) || std::abs(rho) >= 1 ||
        !std::isfinite(sigma) || sigma <= 0 || !std::isfinite(kappa) || kappa <= 0 ||
        odds.n_elem != u.n_slices || !odds.is_finite()) Rcpp::stop("Invalid suffix orientation state.");
  }
  double boundary_log(const arma::cube& u, arma::uword t, const arma::mat& proposed,
      double rho, double sigma, double kappa, const arma::vec& odds) {
    double value = 0;
    for(arma::uword i = 0; i < u.n_rows; ++i) {
      Terms z = terms(proposed.row(i), u.slice(t - 1).row(i), rho, sigma, kappa, odds(t));
      value += log_sum(z.drift, z.reset);
    }
    return value;
  }
  arma::mat orthogonal(double angle, bool reflected) {
    arma::mat q(2, 2);
    q(0, 0) = std::cos(angle); q(1, 0) = std::sin(angle);
    q(0, 1) = reflected ? std::sin(angle) : -std::sin(angle);
    q(1, 1) = reflected ? -std::cos(angle) : std::cos(angle);
    return q;
  }
}


double revision_suffix_log_cpp(const arma::cube& u, int period, const arma::mat& q,
    double rho, double sigma, double kappa, const arma::vec& log_odds) {
  check(u, period, rho, sigma, kappa, log_odds);
  if(q.n_rows != 2 || q.n_cols != 2 || !q.is_finite() || arma::norm(q * q.t() - arma::eye(2, 2), "fro") > 1e-10)
    Rcpp::stop("Expected an orthogonal matrix.");
  return boundary_log(u, period - 1, u.slice(period - 1) * q, rho, sigma, kappa, log_odds);
}


Rcpp::List revision_suffix_period_cpp(const arma::cube& input, const arma::mat& input_delta,
    int period, double rho, double sigma, double kappa, const arma::vec& log_odds, int haar_moves = 1) {
  check(input, period, rho, sigma, kappa, log_odds);
  if(haar_moves < 0 || input_delta.n_rows != input.n_rows || input_delta.n_cols != input.n_slices ||
      !input_delta.is_finite() || arma::any(arma::vectorise((input_delta != 0) % (input_delta != 1))))
    Rcpp::stop("Invalid suffix indicators.");
  arma::cube u(input);
  arma::mat delta(input_delta);
  arma::vec counts(4, arma::fill::zeros);
  const arma::uword t = period - 1;
  double current = boundary_log(u, t, u.slice(t), rho, sigma, kappa, log_odds);
  for(int move = 0; move < 2 + haar_moves; ++move) {
    arma::mat q;
    if(move == 0) q = -arma::eye(2, 2);
    else if(move == 1) { q = arma::eye(2, 2); q(0, 0) = -1; }
    else q = orthogonal(R::runif(-circle_pi, circle_pi), R::runif(0, 1) < .5);
    const double next = boundary_log(u, t, u.slice(t) * q, rho, sigma, kappa, log_odds);
    const double change = next - current;
    const bool take = move < 2 ? R::runif(0, 1) < probability(change) : std::log(R::runif(0, 1)) < std::min(0.0, change);
    counts(move < 2 ? 0 : 2) += 1;
    if(take) {
      for(arma::uword s = t; s < u.n_slices; ++s) u.slice(s) = u.slice(s) * q;
      current = next;
      counts(move < 2 ? 1 : 3) += 1;
    }
  }
  for(arma::uword i = 0; i < u.n_rows; ++i) {
    Terms z = terms(u.slice(t).row(i), u.slice(t - 1).row(i), rho, sigma, kappa, log_odds(t));
    delta(i, t) = R::runif(0, 1) < probability(z.reset - z.drift);
  }
  delta.col(0).zeros();
  return Rcpp::List::create(Rcpp::Named("u") = u, Rcpp::Named("delta") = delta, Rcpp::Named("counts") = counts);
}


Rcpp::List revision_suffix_sweep_cpp(const arma::cube& input, const arma::mat& input_delta,
    double rho, double sigma, double kappa, const arma::vec& log_odds, int haar_moves = 1) {
  check(input, 2, rho, sigma, kappa, log_odds);
  arma::cube u(input);
  arma::mat delta(input_delta);
  arma::vec counts(4, arma::fill::zeros);
  for(arma::uword t = 1; t < u.n_slices; ++t) {
    Rcpp::List result = revision_suffix_period_cpp(u, delta, t + 1, rho, sigma, kappa, log_odds, haar_moves);
    // Cube conversion borrows R memory; keep an owned copy after result leaves scope.
    const arma::cube extracted_u = Rcpp::as<arma::cube>(result["u"]);
    u = extracted_u;
    delta = Rcpp::as<arma::mat>(result["delta"]);
    counts += Rcpp::as<arma::vec>(result["counts"]);
  }
  return Rcpp::List::create(Rcpp::Named("u") = u, Rcpp::Named("delta") = delta, Rcpp::Named("counts") = counts);
}

}
// [[Rcpp::export]]
double revision_suffix_log_cpp(const arma::cube& u, int period, const arma::mat& q,
    double rho, double sigma, double kappa, const arma::vec& log_odds) {
  return suffix::revision_suffix_log_cpp(u, period, q, rho, sigma, kappa, log_odds);
}

// [[Rcpp::export]]
Rcpp::List revision_suffix_period_cpp(const arma::cube& input, const arma::mat& input_delta,
    int period, double rho, double sigma, double kappa, const arma::vec& log_odds, int haar_moves = 1) {
  return suffix::revision_suffix_period_cpp(input, input_delta, period, rho, sigma, kappa, log_odds, haar_moves);
}

// [[Rcpp::export]]
Rcpp::List revision_suffix_sweep_cpp(const arma::cube& input, const arma::mat& input_delta,
    double rho, double sigma, double kappa, const arma::vec& log_odds, int haar_moves = 1) {
  return suffix::revision_suffix_sweep_cpp(input, input_delta, rho, sigma, kappa, log_odds, haar_moves);
}

namespace path {



struct actor_factor {
  int rank, periods;
  std::vector<arma::mat> inverse, lower, root;
  arma::mat mean, direction;
  double a_mean, a_precision;
};

arma::mat path_solve(const actor_factor& f, const arma::mat& rhs) {
  arma::mat forward(rhs), answer(rhs.n_rows, rhs.n_cols);
  for (int t=1; t<f.periods; ++t) forward.col(t) -= f.lower[t] * forward.col(t-1);
  for (int t=0; t<f.periods; ++t) answer.col(t) = f.inverse[t] * forward.col(t);
  for (int t=f.periods-2; t>=0; --t) answer.col(t) -= f.lower[t+1].t() * answer.col(t+1);
  return answer;
}

actor_factor factor_actor(const arma::cube& y, const arma::cube& u, const arma::vec& a,
  double beta, double s2, double va, double rho, double sigma, double kappa,
  const arma::mat& delta, int actor) {
  const int n=u.n_rows, rank=u.n_cols, periods=u.n_slices;
  if (actor<0 || actor>=n || y.n_rows!=(unsigned)n || y.n_cols!=(unsigned)n ||
      y.n_slices!=(unsigned)periods || a.n_elem!=(unsigned)n || delta.n_rows!=(unsigned)n ||
      delta.n_cols!=(unsigned)periods || periods<1 || rank<1 || !u.is_finite() || !a.is_finite() ||
      !std::isfinite(beta+s2+va+rho+sigma+kappa) || s2<=0 || va<=0 || sigma<=0 || kappa<=0 || std::abs(rho)>=1)
    Rcpp::stop("Invalid actor-path conditional inputs.");
  actor_factor f; f.rank=rank; f.periods=periods;
  f.inverse.resize(periods); f.lower.resize(periods); f.root.resize(periods);
  std::vector<arma::mat> precision(periods);
  arma::mat linear(rank,periods,arma::fill::zeros), cross(rank,periods,arma::fill::zeros);
  const arma::mat identity=arma::eye(rank,rank);
  double qa=1/va, ba=0;
  for (int t=0; t<periods; ++t) {
    precision[t]=arma::zeros(rank,rank);
    for (int j=0; j<n; ++j) {
      if (j==actor || !std::isfinite(y(actor,j,t))) continue;
      arma::vec v=u.slice(t).row(j).t();
      const double r=y(actor,j,t)-beta-a[j];
      precision[t] += (v*v.t())/s2;
      linear.col(t) += v*(r/s2); cross.col(t) += v/s2;
      qa += 1/s2; ba += r/s2;
    }
  }
  precision[0] += ((1-rho*rho)/(sigma*sigma))*identity;
  for (int t=1; t<periods; ++t) {
    if (delta(actor,t)!=0 && delta(actor,t)!=1) Rcpp::stop("Invalid reset indicator.");
    const bool reset=delta(actor,t)==1;
    precision[t] += identity/(reset?kappa*kappa:sigma*sigma);
    if (!reset) precision[t-1] += (rho*rho/(sigma*sigma))*identity;
  }
  for (int t=0; t<periods; ++t) {
    f.lower[t]=arma::zeros(rank,rank);
    arma::mat d=precision[t];
    if (t>0 && delta(actor,t)==0) {
      const double h=-rho/(sigma*sigma);
      f.lower[t]=h*f.inverse[t-1];
      d -= h*h*f.inverse[t-1];
    }
    d=(d+d.t())/2;
    if (!arma::chol(f.root[t],d) || !arma::inv_sympd(f.inverse[t],d))
      Rcpp::stop("Actor-path precision is not numerically positive definite.");
  }
  f.mean=path_solve(f,linear); f.direction=path_solve(f,cross);
  f.a_precision=qa-arma::accu(cross%f.direction);
  if (!(f.a_precision>0) || !std::isfinite(f.a_precision)) Rcpp::stop("Invalid additive Schur complement.");
  f.a_mean=(ba-arma::accu(cross%f.mean))/f.a_precision;
  return f;
}

std::pair<arma::mat,double> draw_actor(const actor_factor& f) {
  const double a=f.a_mean+R::rnorm(0,1)/std::sqrt(f.a_precision);
  arma::mat noise(f.rank,f.periods);
  for (int t=0; t<f.periods; ++t) {
    arma::vec z(f.rank); for(int k=0;k<f.rank;++k) z[k]=R::rnorm(0,1);
    noise.col(t)=arma::solve(arma::trimatu(f.root[t]),z);
  }
  for (int t=f.periods-2;t>=0;--t) noise.col(t)-=f.lower[t+1].t()*noise.col(t+1);
  return std::make_pair(f.mean-f.direction*a+noise,a);
}


Rcpp::List actor_path_conditional(const arma::cube& y, const arma::cube& u, const arma::vec& a,
  double beta, double s2, double va, double rho, double sigma, double kappa,
  const arma::mat& delta, int actor, int draws=0) {
  actor_factor f=factor_actor(y,u,a,beta,s2,va,rho,sigma,kappa,delta,actor);
  const int p=f.rank*f.periods;
  arma::vec mean(p+1); mean.head(p)=arma::vectorise(f.mean-f.direction*f.a_mean); mean[p]=f.a_mean;
  arma::mat covariance(p+1,p+1,arma::fill::zeros);
  for(int j=0;j<p;++j) {
    arma::mat unit(f.rank,f.periods,arma::fill::zeros); unit[j]=1;
    covariance.col(j).head(p)=arma::vectorise(path_solve(f,unit));
  }
  arma::vec v=arma::vectorise(f.direction);
  covariance.submat(0,0,p-1,p-1)+=v*v.t()/f.a_precision;
  covariance.col(p).head(p)=-v/f.a_precision; covariance.row(p).head(p)=(-v/f.a_precision).t();
  covariance(p,p)=1/f.a_precision;
  arma::mat sample(draws,p+1);
  for(int j=0;j<draws;++j) { auto d=draw_actor(f); sample.row(j).head(p)=arma::vectorise(d.first).t(); sample(j,p)=d.second; }
  return Rcpp::List::create(Rcpp::Named("mean")=mean,Rcpp::Named("covariance")=covariance,Rcpp::Named("draws")=sample);
}


Rcpp::List actor_path_sweep(const arma::cube& y, const arma::cube& u, const arma::vec& a,
  double beta, double s2, double va, double rho, double sigma, double kappa,
  const arma::mat& delta) {
  arma::cube next(u.n_rows,u.n_cols,u.n_slices); std::copy(u.begin(),u.end(),next.begin());
  arma::vec next_a(a.n_elem);std::copy(a.begin(),a.end(),next_a.begin());
  for(unsigned i=0;i<u.n_rows;++i) {
    actor_factor f=factor_actor(y,next,next_a,beta,s2,va,rho,sigma,kappa,delta,i);
    auto d=draw_actor(f);next_a[i]=d.second;
    for(unsigned t=0;t<u.n_slices;++t)next.slice(t).row(i)=d.first.col(t).t();
  }
  return Rcpp::List::create(Rcpp::Named("u")=next,Rcpp::Named("a")=next_a);
}

}
// [[Rcpp::export]]
Rcpp::List actor_path_conditional(const arma::cube& y, const arma::cube& u, const arma::vec& a,
  double beta, double s2, double va, double rho, double sigma, double kappa,
  const arma::mat& delta, int actor, int draws=0) {
  return path::actor_path_conditional(y, u, a, beta, s2, va, rho, sigma, kappa, delta, actor, draws);
}

// [[Rcpp::export]]
Rcpp::List actor_path_sweep(const arma::cube& y, const arma::cube& u, const arma::vec& a,
  double beta, double s2, double va, double rho, double sigma, double kappa,
  const arma::mat& delta) {
  return path::actor_path_sweep(y, u, a, beta, s2, va, rho, sigma, kappa, delta);
}

namespace centroid {

struct centroid_model {
  int n, rank, periods, dimension;
  double s2;
  arma::cube shape, row_v;
  arma::uvec edge_i, edge_j, edge_t;
  arma::vec edge_c;
  arma::mat edge_v;
  arma::vec level_prior;
  arma::mat level_precision, level_root, level_inverse, prior_precision;
  arma::vec prior_linear, count, sum_c, sum_c2;
  arma::mat sum_v, sum_vc, row_count, row_c;
  std::vector<arma::mat> sum_vv;
  arma::vec start;
};

struct centroid_value {
  double log_density, logdet;
  arma::vec gradient, level_mean;
  arma::mat precision, root;
};

centroid_model prepare_centroid(const arma::cube& y, const arma::cube& u, double s2,
  double va, double beta_precision, double rho, double sigma, double kappa, const arma::mat& delta) {
  centroid_model f; f.n=u.n_rows;f.rank=u.n_cols;f.periods=u.n_slices;f.dimension=f.rank*f.periods;f.s2=s2;
  const int n=f.n,r=f.rank,T=f.periods,p=f.dimension;
  if(n<2||r<1||T<1||y.n_rows!=(unsigned)n||y.n_cols!=(unsigned)n||y.n_slices!=(unsigned)T||
    delta.n_rows!=(unsigned)n||delta.n_cols!=(unsigned)T||!u.is_finite()||!delta.is_finite()||
    !std::isfinite(s2+va+beta_precision+rho+sigma+kappa)||s2<=0||va<=0||beta_precision<=0||sigma<=0||kappa<=0||std::abs(rho)>=1)
    Rcpp::stop("Invalid centroid conditional inputs.");
  f.shape.set_size(n,r,T);f.start.set_size(p);
  const unsigned max_edges=(unsigned)n*(n-1)*T/2;unsigned edge=0;
  f.edge_i.set_size(max_edges);f.edge_j.set_size(max_edges);f.edge_t.set_size(max_edges);
  f.edge_c.set_size(max_edges);f.edge_v.set_size(r,max_edges);
  f.row_v.zeros(n,r,T);f.row_count.zeros(n,T);f.row_c.zeros(n,T);
  f.count.zeros(T);f.sum_c.zeros(T);f.sum_c2.zeros(T);f.sum_v.zeros(r,T);f.sum_vc.zeros(r,T);f.sum_vv.resize(T);
  f.level_prior=arma::vec(n+1,arma::fill::ones)/va;f.level_prior[0]=beta_precision;
  f.level_precision.zeros(n+1,n+1);f.level_precision(0,0)=beta_precision;
  for(int i=0;i<n;++i)f.level_precision(i+1,i+1)=1/va;
  for(int t=0;t<T;++t){
    arma::rowvec m=arma::mean(u.slice(t),0);f.start.subvec(t*r,(t+1)*r-1)=m.t();
    f.shape.slice(t)=u.slice(t);f.shape.slice(t).each_row()-=m;
    f.sum_vv[t].zeros(r,r);
    for(int j=1;j<n;++j)for(int i=0;i<j;++i){
      const double response=y(i,j,t);
      if(std::isfinite(response)!=std::isfinite(y(j,i,t))||(std::isfinite(response)&&std::abs(response-y(j,i,t))>1e-10))
        Rcpp::stop("Centroid updates require symmetric observed dyads.");
      if(!std::isfinite(response))continue;
      const arma::vec v=(f.shape.slice(t).row(i)+f.shape.slice(t).row(j)).t();
      const double c=response-arma::dot(f.shape.slice(t).row(i),f.shape.slice(t).row(j));
      f.edge_i[edge]=i;f.edge_j[edge]=j;f.edge_t[edge]=t;f.edge_c[edge]=c;f.edge_v.col(edge)=v;edge++;
      f.count[t]++;f.sum_c[t]+=c;f.sum_c2[t]+=c*c;f.sum_v.col(t)+=v;f.sum_vc.col(t)+=v*c;f.sum_vv[t]+=v*v.t();
      for(int a : {i,j}){f.row_count(a,t)++;f.row_c(a,t)+=c;f.row_v.slice(t).row(a)+=v.t();}
      const int indices[3]={0,i+1,j+1};
      for(int a=0;a<3;++a)for(int b=0;b<3;++b)f.level_precision(indices[a],indices[b])+=1/s2;
    }
  }
  f.edge_i.resize(edge);f.edge_j.resize(edge);f.edge_t.resize(edge);f.edge_c.resize(edge);f.edge_v.resize(r,edge);
  if(!arma::chol(f.level_root,f.level_precision)||!arma::inv_sympd(f.level_inverse,f.level_precision))
    Rcpp::stop("Centroid level precision is not positive definite.");
  f.prior_precision.zeros(p,p);f.prior_linear.zeros(p);
  const arma::mat eye=arma::eye(r,r);
  const double initial=(1-rho*rho)/(sigma*sigma);
  f.prior_precision.submat(0,0,r-1,r-1)=n*initial*eye;
  f.prior_linear.head(r)=-initial*arma::sum(f.shape.slice(0),0).t();
  for(int t=1;t<T;++t){
    int drift=0,reset=0;arma::vec drift_sum(r,arma::fill::zeros),reset_sum(r,arma::fill::zeros);
    for(int i=0;i<n;++i){
      if(delta(i,t)==0){drift++;drift_sum+=(f.shape.slice(t).row(i)-rho*f.shape.slice(t-1).row(i)).t();}
      else if(delta(i,t)==1){reset++;reset_sum+=f.shape.slice(t).row(i).t();}
      else Rcpp::stop("Invalid centroid reset indicator.");
    }
    const int at=t*r,before=(t-1)*r;
    f.prior_precision.submat(at,at,at+r-1,at+r-1)+=(drift/(sigma*sigma)+reset/(kappa*kappa))*eye;
    f.prior_precision.submat(before,before,before+r-1,before+r-1)+=(drift*rho*rho/(sigma*sigma))*eye;
    f.prior_precision.submat(at,before,at+r-1,before+r-1)-=(drift*rho/(sigma*sigma))*eye;
    f.prior_precision.submat(before,at,before+r-1,at+r-1)-=(drift*rho/(sigma*sigma))*eye;
    f.prior_linear.subvec(at,at+r-1)-=drift_sum/(sigma*sigma)+reset_sum/(kappa*kappa);
    f.prior_linear.subvec(before,before+r-1)+=rho*drift_sum/(sigma*sigma);
  }
  return f;
}

centroid_value evaluate_centroid(const centroid_model& f,const arma::vec& m){
  const int n=f.n,r=f.rank,T=f.periods,p=f.dimension;
  if(m.n_elem!=(unsigned)p||!m.is_finite())Rcpp::stop("Invalid centroid coordinates.");
  arma::vec b(n+1,arma::fill::zeros),sum_residual(T);
  arma::mat cross(n+1,p,arma::fill::zeros),jtj(p,p,arma::fill::zeros),sum_v_residual(r,T);
  for(int t=0;t<T;++t){
    const int at=t*r;const arma::vec now=m.subvec(at,at+r-1);const double norm=arma::dot(now,now);
    const double total=f.sum_c[t]-arma::dot(f.sum_v.col(t),now)-f.count[t]*norm;
    const arma::vec rows=f.row_c.col(t)-f.row_v.slice(t)*now-f.row_count.col(t)*norm;
    b[0]+=total/f.s2;b.tail(n)+=rows/f.s2;
    cross.submat(0,at,0,at+r-1)=(f.sum_v.col(t)+2*f.count[t]*now).t()/f.s2;
    cross.submat(1,at,n,at+r-1)=(f.row_v.slice(t)+2*f.row_count.col(t)*now.t())/f.s2;
    jtj.submat(at,at,at+r-1,at+r-1)=(f.sum_vv[t]+2*f.sum_v.col(t)*now.t()+
      2*now*f.sum_v.col(t).t()+4*f.count[t]*now*now.t())/f.s2;
    sum_residual[t]=total;sum_v_residual.col(t)=f.sum_vc.col(t)-f.sum_vv[t]*now-f.sum_v.col(t)*norm;
  }
  centroid_value z;z.level_mean=f.level_inverse*b;
  z.gradient=f.prior_linear-f.prior_precision*m;
  const double beta=z.level_mean[0];const arma::vec a=z.level_mean.tail(n);
  long double likelihood_quadratic=0;
  arma::vec norms(T);for(int t=0;t<T;++t)norms[t]=arma::dot(m.subvec(t*r,(t+1)*r-1),m.subvec(t*r,(t+1)*r-1));
  for(unsigned edge=0;edge<f.edge_c.n_elem;++edge){
    const int i=f.edge_i[edge],j=f.edge_j[edge],t=f.edge_t[edge],at=t*r;
    double residual=f.edge_c[edge]-norms[t]-beta-a[i]-a[j];
    for(int d=0;d<r;++d)residual-=f.edge_v(d,edge)*m[at+d];
    likelihood_quadratic+=(long double)residual*residual/f.s2;
    for(int d=0;d<r;++d)z.gradient[at+d]+=(f.edge_v(d,edge)+2*m[at+d])*(residual/f.s2);
  }
  likelihood_quadratic+=arma::dot(f.level_prior,arma::square(z.level_mean));
  z.log_density=-.5*(double)likelihood_quadratic-
    .5*arma::as_scalar(m.t()*f.prior_precision*m)+arma::dot(f.prior_linear,m);
  z.precision=f.prior_precision+jtj-cross.t()*f.level_inverse*cross;
  z.precision=(z.precision+z.precision.t())/2;
  if(!std::isfinite(z.log_density)||!z.gradient.is_finite()||!arma::chol(z.root,z.precision))
    Rcpp::stop("Centroid target or proposal precision exceeded numerical support.");
  z.logdet=2*arma::sum(arma::log(z.root.diag()));
  return z;
}

arma::vec centroid_mean(const arma::vec& m,const centroid_value& z,double step){
  return m+.5*step*step*arma::solve(arma::trimatu(z.root),arma::solve(arma::trimatl(z.root.t()),z.gradient));
}

double centroid_log_proposal(const arma::vec& to,const arma::vec& from,const centroid_value& z,double step){
  const arma::vec difference=to-centroid_mean(from,z,step);
  return .5*z.logdet-from.n_elem*std::log(step)-.5*arma::as_scalar(difference.t()*z.precision*difference)/(step*step);
}

SEXP centroid_prepare(const arma::cube& y,const arma::cube& u,double s2,double va,double beta_precision,
  double rho,double sigma,double kappa,const arma::mat& delta){
  return Rcpp::XPtr<centroid_model>(new centroid_model(prepare_centroid(y,u,s2,va,beta_precision,rho,sigma,kappa,delta)),true);
}

Rcpp::List centroid_evaluate(SEXP object,const arma::vec& m){
  Rcpp::XPtr<centroid_model> f(object);centroid_value z=evaluate_centroid(*f,m);
  return Rcpp::List::create(Rcpp::Named("log_density")=z.log_density,Rcpp::Named("gradient")=z.gradient,
    Rcpp::Named("precision")=z.precision,Rcpp::Named("level_mean")=z.level_mean,
    Rcpp::Named("level_covariance")=f->level_inverse,Rcpp::Named("shape")=f->shape,Rcpp::Named("start")=f->start);
}

double centroid_proposal_log(SEXP object,const arma::vec& to,const arma::vec& from,double step){
  Rcpp::XPtr<centroid_model> f(object);centroid_value z=evaluate_centroid(*f,from);
  return centroid_log_proposal(to,from,z,step);
}

Rcpp::List centroid_draw(SEXP object,const arma::vec& initial,double step,int moves=1){
  if(!std::isfinite(step)||step<=0||moves<1)Rcpp::stop("Invalid centroid proposal settings.");
  Rcpp::XPtr<centroid_model> f(object);arma::vec m(initial);centroid_value current=evaluate_centroid(*f,m);
  int accepted=0;double probability_sum=0;
  for(int k=0;k<moves;++k){
    arma::vec noise(m.n_elem);for(unsigned j=0;j<m.n_elem;++j)noise[j]=R::rnorm(0,1);
    arma::vec proposed=centroid_mean(m,current,step)+step*arma::solve(arma::trimatu(current.root),noise);
    centroid_value next=evaluate_centroid(*f,proposed);
    const double log_ratio=next.log_density-current.log_density+
      centroid_log_proposal(m,proposed,next,step)-centroid_log_proposal(proposed,m,current,step);
    if(!std::isfinite(log_ratio))Rcpp::stop("Nonfinite centroid MH ratio.");
    probability_sum+=std::exp(std::min(0.0,log_ratio));
    if(std::log(R::runif(0,1))<log_ratio){m=proposed;current=next;accepted++;}
  }
  arma::vec noise(f->n+1);for(int j=0;j<f->n+1;++j)noise[j]=R::rnorm(0,1);
  arma::vec levels=current.level_mean+arma::solve(arma::trimatu(f->level_root),noise);
  arma::cube u(f->shape);for(int t=0;t<f->periods;++t)u.slice(t).each_row()+=m.subvec(t*f->rank,(t+1)*f->rank-1).t();
  return Rcpp::List::create(Rcpp::Named("u")=u,Rcpp::Named("a")=levels.tail(f->n),Rcpp::Named("beta")=levels[0],
    Rcpp::Named("centroids")=m,Rcpp::Named("accepted")=accepted,Rcpp::Named("proposed")=moves,
    Rcpp::Named("mean_acceptance")=probability_sum/moves);
}

}

// [[Rcpp::export]]
SEXP centroid_prepare(const arma::cube& y,const arma::cube& u,double s2,double va,double beta_precision,
  double rho,double sigma,double kappa,const arma::mat& delta){
  return centroid::centroid_prepare(y,u,s2,va,beta_precision,rho,sigma,kappa,delta);
}
// [[Rcpp::export]]
Rcpp::List centroid_evaluate(SEXP object,const arma::vec& m){return centroid::centroid_evaluate(object,m);}
// [[Rcpp::export]]
double centroid_proposal_log(SEXP object,const arma::vec& to,const arma::vec& from,double step){
  return centroid::centroid_proposal_log(object,to,from,step);
}
// [[Rcpp::export]]
Rcpp::List centroid_draw(SEXP object,const arma::vec& initial,double step,int moves=1){
  return centroid::centroid_draw(object,initial,step,moves);
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

struct jump_value {
  double log_density;
  arma::vec level_mean;
  arma::mat level_root, probabilities;
};

jump_value evaluate_jump(const arma::cube& y,const arma::cube& u,double s2,double va,
 double beta_precision,double rho,double sigma,double kappa,double pi){
 const int n=u.n_rows,r=u.n_cols,T=u.n_slices;
 if(n<2||r<1||T<1||y.n_rows!=(unsigned)n||y.n_cols!=(unsigned)n||y.n_slices!=(unsigned)T||
  !u.is_finite()||!std::isfinite(s2+va+beta_precision+rho+sigma+kappa+pi)||
  s2<=0||va<=0||beta_precision<=0||sigma<=0||kappa<=0||pi<=0||pi>=1||std::abs(rho)>=1)
  Rcpp::stop("Invalid joint proposal inputs.");
 arma::mat precision(n+1,n+1,arma::fill::zeros);
 arma::vec linear(n+1,arma::fill::zeros),prior(n+1,arma::fill::ones);
 prior/=va;prior[0]=beta_precision;precision.diag()=prior;
 for(int t=0;t<T;++t)for(int j=1;j<n;++j)for(int i=0;i<j;++i){
  const double response=y(i,j,t);
  if(std::isfinite(response)!=std::isfinite(y(j,i,t))||
   (std::isfinite(response)&&std::abs(response-y(j,i,t))>1e-10))Rcpp::stop("Asymmetric observation.");
  if(!std::isfinite(response))continue;
  const double residual=response-arma::dot(u.slice(t).row(i),u.slice(t).row(j));
  int idx[3]={0,i+1,j+1};
  for(int a=0;a<3;++a){linear[idx[a]]+=residual/s2;for(int b=0;b<3;++b)precision(idx[a],idx[b])+=1/s2;}
 }
 jump_value out;
 if(!arma::chol(out.level_root,precision))Rcpp::stop("Invalid joint level precision.");
 out.level_mean=arma::solve(arma::trimatu(out.level_root),arma::solve(arma::trimatl(out.level_root.t()),linear));
 long double quadratic=arma::dot(prior,arma::square(out.level_mean));
 for(int t=0;t<T;++t)for(int j=1;j<n;++j)for(int i=0;i<j;++i){
  if(!std::isfinite(y(i,j,t)))continue;
  const double residual=y(i,j,t)-arma::dot(u.slice(t).row(i),u.slice(t).row(j))-
   out.level_mean[0]-out.level_mean[i+1]-out.level_mean[j+1];
  quadratic+=(long double)residual*residual/s2;
 }
 out.log_density=-.5*(double)quadratic;
 out.probabilities.zeros(n,T);
 const double initial_variance=sigma*sigma/(1-rho*rho);
 out.log_density-=.5*arma::accu(arma::square(u.slice(0)))/initial_variance;
 for(int t=1;t<T;++t)for(int i=0;i<n;++i){
  const arma::rowvec difference=u.slice(t).row(i)-rho*u.slice(t-1).row(i);
  const double drift=std::log1p(-pi)-r*std::log(sigma)-.5*arma::dot(difference,difference)/(sigma*sigma);
  const double reset=std::log(pi)-r*std::log(kappa)-.5*arma::dot(u.slice(t).row(i),u.slice(t).row(i))/(kappa*kappa);
  const double largest=std::max(drift,reset);
  const double logmix=largest+std::log(std::exp(drift-largest)+std::exp(reset-largest));
  out.log_density+=logmix;
  out.probabilities(i,t)=std::exp(reset-logmix);
 }
 if(!std::isfinite(out.log_density))Rcpp::stop("Nonfinite joint proposal density.");
 return out;
}

// [[Rcpp::export]]
Rcpp::List jump_evaluate(const arma::cube& y,const arma::cube& u,double s2,double va,
 double beta_precision,double rho,double sigma,double kappa,double pi){
 jump_value x=evaluate_jump(y,u,s2,va,beta_precision,rho,sigma,kappa,pi);
 return Rcpp::List::create(Rcpp::_["log_density"]=x.log_density,Rcpp::_["level_mean"]=x.level_mean,
 Rcpp::_["level_precision"]=x.level_root.t()*x.level_root,Rcpp::_["probabilities"]=x.probabilities);
}

// [[Rcpp::export]]
Rcpp::List jump_draw(const arma::cube& y,const arma::cube& u,const arma::cube& shift,
 double s2,double va,double beta_precision,double rho,double sigma,double kappa,double pi){
 if(shift.n_rows!=u.n_rows||shift.n_cols!=u.n_cols||shift.n_slices!=u.n_slices||!shift.is_finite())
  Rcpp::stop("Invalid frozen joint shift.");
 const double direction=R::runif(0,1)<.5?-1:1;
 const arma::cube proposal=u+direction*shift;
 jump_value current=evaluate_jump(y,u,s2,va,beta_precision,rho,sigma,kappa,pi);
 jump_value next=evaluate_jump(y,proposal,s2,va,beta_precision,rho,sigma,kappa,pi);
 const double log_ratio=next.log_density-current.log_density;
 const bool accepted=std::log(R::runif(0,1))<log_ratio;
 const jump_value& x=accepted?next:current;
 const arma::cube selected=accepted?proposal:u;
 arma::vec noise(x.level_mean.n_elem);
 for(unsigned i=0;i<noise.n_elem;++i)noise[i]=R::rnorm(0,1);
 const arma::vec levels=x.level_mean+arma::solve(arma::trimatu(x.level_root),noise);
 arma::mat delta(u.n_rows,u.n_slices,arma::fill::zeros);
 for(unsigned t=1;t<u.n_slices;++t)for(unsigned i=0;i<u.n_rows;++i)delta(i,t)=R::runif(0,1)<x.probabilities(i,t);
 return Rcpp::List::create(Rcpp::_["u"]=selected,Rcpp::_["a"]=levels.tail(u.n_rows),
 Rcpp::_["beta"]=levels[0],Rcpp::_["delta"]=delta,Rcpp::_["accepted"]=accepted,
 Rcpp::_["direction"]=direction,Rcpp::_["log_ratio"]=log_ratio);
}

// Integrate an actor path and static effect when updating reset indicators.
using path::actor_factor;
using path::path_solve;
using path::draw_actor;
void check_indicator_inputs(const arma::cube& y,const arma::cube& u,const arma::vec& a,
 double beta,double s2,double va,double rho,double sigma,double kappa,double pi,const arma::mat& delta){
 if(u.n_rows<2 || u.n_cols<1 || u.n_slices<1 || y.n_rows!=u.n_rows || y.n_cols!=u.n_rows ||
  y.n_slices!=u.n_slices || a.n_elem!=u.n_rows || delta.n_rows!=u.n_rows || delta.n_cols!=u.n_slices ||
  !u.is_finite() || !a.is_finite() || !delta.is_finite() || !std::isfinite(beta+s2+va+rho+sigma+kappa+pi) ||
  s2<=0 || va<=0 || sigma<=0 || kappa<=0 || std::abs(rho)>=1 || pi<=0 || pi>=1)
  Rcpp::stop("Invalid collapsed indicator inputs.");
 for(double value:delta)if(value!=0 && value!=1)Rcpp::stop("Invalid reset indicator.");
 if(arma::any(delta.col(0)!=0))Rcpp::stop("Initial reset indicators must be zero.");
}
struct actor_data {
 int rank, periods;
 std::vector<arma::mat> precision;
 arma::mat linear, cross;
 double qa, ba;
};

actor_data cache_actor(const arma::cube& y,const arma::cube& u,const arma::vec& a,
 double beta,double s2,double va,int actor){
 actor_data d; d.rank=u.n_cols;d.periods=u.n_slices;
 d.precision.resize(d.periods); d.linear.zeros(d.rank,d.periods);d.cross.zeros(d.rank,d.periods);
 d.qa=1/va;d.ba=0;
 for(int t=0;t<d.periods;++t){
  d.precision[t].zeros(d.rank,d.rank);
  for(unsigned j=0;j<u.n_rows;++j){
   if((int)j==actor || !std::isfinite(y(actor,j,t)))continue;
   arma::vec v=u.slice(t).row(j).t(); double r=y(actor,j,t)-beta-a[j];
   d.precision[t]+=v*v.t()/s2;d.linear.col(t)+=v*(r/s2);d.cross.col(t)+=v/s2;
   d.qa+=1/s2;d.ba+=r/s2;
  }
 }
 return d;
}

actor_factor cached_factor(const actor_data& d,const arma::rowvec& delta,double rho,double sigma,double kappa){
 actor_factor f;f.rank=d.rank;f.periods=d.periods;
 f.inverse.resize(d.periods);f.lower.resize(d.periods);f.root.resize(d.periods);
 std::vector<arma::mat> precision=d.precision;
 const arma::mat identity=arma::eye(d.rank,d.rank);
 precision[0]+=(1-rho*rho)/(sigma*sigma)*identity;
 for(int t=1;t<d.periods;++t){
  bool reset=delta[t]==1;precision[t]+=identity/(reset?kappa*kappa:sigma*sigma);
  if(!reset)precision[t-1]+=rho*rho/(sigma*sigma)*identity;
 }
 for(int t=0;t<d.periods;++t){
  f.lower[t].zeros(d.rank,d.rank);arma::mat q=precision[t];
  if(t>0 && delta[t]==0){double h=-rho/(sigma*sigma);f.lower[t]=h*f.inverse[t-1];q-=h*h*f.inverse[t-1];}
  q=(q+q.t())/2;
  if(!arma::chol(f.root[t],q)||!arma::inv_sympd(f.inverse[t],q))Rcpp::stop("Invalid collapsed actor precision.");
 }
 f.mean=path_solve(f,d.linear);f.direction=path_solve(f,d.cross);
 f.a_precision=d.qa-arma::accu(d.cross%f.direction);
 if(!(f.a_precision>0))Rcpp::stop("Invalid collapsed additive precision.");
 f.a_mean=(d.ba-arma::accu(d.cross%f.mean))/f.a_precision;
 return f;
}

double actor_evidence(const actor_data& d,const actor_factor& f,const arma::rowvec& delta,double sigma,double kappa,double pi){
 long double quadratic=0;
 for(unsigned j=0;j<d.linear.n_elem;++j)quadratic+=(long double)d.linear[j]*f.mean[j];
 quadratic+=(long double)f.a_precision*f.a_mean*f.a_mean;
 double value=.5*(double)quadratic-.5*std::log(f.a_precision);
 for(int t=0;t<d.periods;++t)value-=arma::accu(arma::log(f.root[t].diag()));
 double resets=arma::accu(delta.tail(d.periods-1));
 value+=resets*(std::log(pi)-d.rank*std::log(kappa))+(d.periods-1-resets)*(std::log1p(-pi)-d.rank*std::log(sigma));
 return value;
}

// [[Rcpp::export]]
Rcpp::List indicator_evaluate(const arma::cube& y,const arma::cube& u,const arma::vec& a,
 double beta,double s2,double va,double rho,double sigma,double kappa,double pi,
 const arma::mat& delta,int actor){
 check_indicator_inputs(y,u,a,beta,s2,va,rho,sigma,kappa,pi,delta);
 if(actor<0 || actor>=(int)u.n_rows)Rcpp::stop("Invalid actor index.");
 actor_data d=cache_actor(y,u,a,beta,s2,va,actor);arma::rowvec indicators=delta.row(actor);
 actor_factor f=cached_factor(d,indicators,rho,sigma,kappa);
 return Rcpp::List::create(Rcpp::_["log_evidence"]=actor_evidence(d,f,indicators,sigma,kappa,pi),
 Rcpp::_["path_mean"]=f.mean-f.direction*f.a_mean,Rcpp::_["a_mean"]=f.a_mean);
}

// [[Rcpp::export]]
Rcpp::List indicator_pair_sweep(const arma::cube& y,const arma::cube& u,const arma::vec& a,
 double beta,double s2,double va,double rho,double sigma,double kappa,double pi,
 const arma::mat& delta, int selected_actor=-1, int repetitions=1){
 check_indicator_inputs(y,u,a,beta,s2,va,rho,sigma,kappa,pi,delta);
 if(repetitions<1||selected_actor < -1||selected_actor>=(int)u.n_rows)Rcpp::stop("Invalid indicator update.");
 arma::cube next(u);arma::vec next_a(a);arma::mat next_delta(delta);int changed=0;
 for(int repeat=0;repeat<repetitions;++repeat)for(unsigned i=0;i<u.n_rows;++i){
  if(selected_actor>=0 && (int)i!=selected_actor)continue;
  actor_data d=cache_actor(y,next,next_a,beta,s2,va,i);arma::rowvec indicators=next_delta.row(i);
  actor_factor f=cached_factor(d,indicators,rho,sigma,kappa);
  double value=actor_evidence(d,f,indicators,sigma,kappa,pi);
  for(unsigned t=1;t<u.n_slices;++t){
   arma::rowvec proposal(indicators);proposal[t]=1-indicators[t];
   actor_factor candidate=cached_factor(d,proposal,rho,sigma,kappa);
   double proposed=actor_evidence(d,candidate,proposal,sigma,kappa,pi),difference=proposed-value;
   double probability=difference>=0?1/(1+std::exp(-difference)):std::exp(difference)/(1+std::exp(difference));
   if(R::runif(0,1)<probability){indicators=proposal;f=std::move(candidate);value=proposed;++changed;}
  }
  // Adjacent complementary patterns are a second exact two-state Gibbs block.
  // This permits a reset to move between neighboring years without crossing an unlikely intermediate pattern.
  for(unsigned t=1;t+1<u.n_slices;++t){
   arma::rowvec proposal(indicators);proposal[t]=1-indicators[t];proposal[t+1]=1-indicators[t+1];
   actor_factor candidate=cached_factor(d,proposal,rho,sigma,kappa);
   double proposed=actor_evidence(d,candidate,proposal,sigma,kappa,pi),difference=proposed-value;
   double probability=difference>=0?1/(1+std::exp(-difference)):std::exp(difference)/(1+std::exp(difference));
   if(R::runif(0,1)<probability){indicators=proposal;f=std::move(candidate);value=proposed;changed+=2;}
  }
  auto draw=draw_actor(f);next_a[i]=draw.second;next_delta.row(i)=indicators;
  for(unsigned t=0;t<u.n_slices;++t)next.slice(t).row(i)=draw.first.col(t).t();
 }
 return Rcpp::List::create(Rcpp::_["u"]=next,Rcpp::_["a"]=next_a,Rcpp::_["delta"]=next_delta,Rcpp::_["changed"]=changed);
}

// Change the grand centroid radius and inversely scale only within-actor
// temporal deviations in its direction. Static actor deviations stay fixed.
// Jacobian exponent: rank - 2 - actors * (periods - 1).
// [[Rcpp::export]]
Rcpp::List temporal_map(const arma::cube& u,double increment){
 if(!u.is_finite() || !std::isfinite(increment) || u.n_rows<2 || u.n_cols<1 || u.n_slices<1)
  Rcpp::stop("Invalid temporal transport input.");
 arma::mat actor_mean(u.n_rows,u.n_cols,arma::fill::zeros);
 for(unsigned t=0;t<u.n_slices;++t)actor_mean+=u.slice(t)/u.n_slices;
 arma::rowvec m=arma::mean(actor_mean,0);double square=arma::dot(m,m);
 if(square<=0 || square+increment<=0)return Rcpp::List::create(Rcpp::_["valid"]=false);
 double ratio=std::sqrt((square+increment)/square);
 arma::rowvec direction=m/std::sqrt(square);arma::cube proposed(u);
 for(unsigned t=0;t<u.n_slices;++t){
  arma::mat deviation=u.slice(t)-actor_mean;
  proposed.slice(t)+=(1/ratio-1)*(deviation*direction.t())*direction;
  proposed.slice(t).each_row()+=(ratio-1)*m;
 }
 double exponent=(double)u.n_cols-2-(double)u.n_rows*(u.n_slices-1);
 return Rcpp::List::create(Rcpp::_["valid"]=true,Rcpp::_["u"]=proposed,Rcpp::_["log_jacobian"]=exponent*std::log(ratio));
}

// [[Rcpp::export]]
Rcpp::List temporal_evaluate(const arma::cube& y,const arma::cube& u,double increment,
 double s2,double va,double beta_precision,double rho,double sigma,double kappa,double pi){
 Rcpp::List mapping=temporal_map(u,increment);
 if(!Rcpp::as<bool>(mapping["valid"]))return Rcpp::List::create(Rcpp::_["valid"]=false,Rcpp::_["log_ratio"]=-arma::datum::inf);
 const arma::cube proposed=Rcpp::as<arma::cube>(mapping["u"]);
 jump_value current=evaluate_jump(y,u,s2,va,beta_precision,rho,sigma,kappa,pi);
 jump_value next=evaluate_jump(y,proposed,s2,va,beta_precision,rho,sigma,kappa,pi);
 return Rcpp::List::create(Rcpp::_["valid"]=true,Rcpp::_["u"]=proposed,
  Rcpp::_["log_ratio"]=next.log_density-current.log_density+Rcpp::as<double>(mapping["log_jacobian"]),
  Rcpp::_["beta_change"]=next.level_mean[0]-current.level_mean[0]);
}

// [[Rcpp::export]]
Rcpp::List temporal_draw(const arma::cube& y,const arma::cube& u,double tau,int moves,
 double s2,double va,double beta_precision,double rho,double sigma,double kappa,double pi){
 if(!std::isfinite(tau)||!(tau>0)||moves<1)Rcpp::stop("Invalid radius proposal scale or count.");
 arma::cube selected(u);jump_value current=evaluate_jump(y,selected,s2,va,beta_precision,rho,sigma,kappa,pi);
 int accepted=0;
 for(int k=0;k<moves;++k){
  Rcpp::List mapping=temporal_map(selected,R::rnorm(0,tau));double log_uniform=std::log(R::runif(0,1));
  if(!Rcpp::as<bool>(mapping["valid"]))continue;
  arma::cube proposal=Rcpp::as<arma::cube>(mapping["u"]);
  jump_value candidate=evaluate_jump(y,proposal,s2,va,beta_precision,rho,sigma,kappa,pi);
  double ratio=candidate.log_density-current.log_density+Rcpp::as<double>(mapping["log_jacobian"]);
  // proposal borrows mapping's R array, so selected must copy its values.
  if(log_uniform<ratio){selected=proposal;current=std::move(candidate);++accepted;}
 }
 arma::vec noise(current.level_mean.n_elem);for(unsigned j=0;j<noise.n_elem;++j)noise[j]=R::rnorm(0,1);
 arma::vec levels=current.level_mean+arma::solve(arma::trimatu(current.level_root),noise);
 arma::mat delta(u.n_rows,u.n_slices,arma::fill::zeros);
 for(unsigned t=1;t<u.n_slices;++t)for(unsigned i=0;i<u.n_rows;++i)delta(i,t)=R::runif(0,1)<current.probabilities(i,t);
 return Rcpp::List::create(Rcpp::_["u"]=selected,Rcpp::_["a"]=levels.tail(u.n_rows),
 Rcpp::_["beta"]=levels[0],Rcpp::_["delta"]=delta,Rcpp::_["accepted"]=accepted,Rcpp::_["proposed"]=moves);
}
