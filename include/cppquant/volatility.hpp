#pragma once

#include "cppquant/options.hpp"

#include <vector>

namespace cppquant {

/// Raw SVI total-variance parameters:
///
///   w(k) = a + b (rho (k - m) + sqrt((k - m)^2 + sigma^2))
///
/// where `k = log(strike / forward)` and `w` is total variance
/// (`implied_vol^2 * time_to_expiry`).
struct SviParameters {
  double a{};
  double b{};
  double rho{};
  double m{};
  double sigma{};
};

/// One implied-volatility observation for a single expiry.
struct VolQuote {
  double strike{};
  double implied_volatility{};
  /// Relative weight in the least-squares objective. Must be positive.
  double weight{1.0};
};

struct SviFitSettings {
  int max_iterations{200};
  /// Cost tolerance and step tolerance for the optimizer.
  double tolerance{1e-12};
  double lambda_initial{1e-3};
};

/// A calibrated SVI slice for one forward and expiry.
struct SviSlice {
  double forward{};
  double time_to_expiry{};
  SviParameters parameters{};

  /// Raw total variance at `k = log(strike / forward)`.
  double total_variance(double log_moneyness) const noexcept;

  /// Implied volatility at `k = log(strike / forward)`.
  double implied_volatility(double log_moneyness) const noexcept;

  double implied_volatility_for_strike(double strike) const noexcept;

  /// Checks `b >= 0`, `|rho| < 1`, `sigma > 0` and the SVI lower bound
  /// `a + b sigma sqrt(1 - rho^2) >= 0`.
  bool satisfies_basic_constraints(double tolerance = 1e-12) const noexcept;
};

struct SviFitResult {
  SviSlice slice;
  /// Weighted RMSE in implied-volatility units.
  double rmse{};
  double max_error{};
  int iterations{};
  bool converged{};
};

/// Calibrate one raw SVI slice to a set of quotes for the same expiry.
///
/// The optimizer works in an unconstrained space: `b`, `sigma` and the SVI
/// minimum variance are exponentiated, `rho` is transformed with `tanh`, and
/// several deterministic starting points are tried. The best fit is returned.
///
/// Throws std::invalid_argument for non-positive forwards/expiries, fewer than
/// five distinct strikes, invalid quotes or invalid settings.
SviFitResult calibrate_svi(double forward, double time_to_expiry,
                           const std::vector<VolQuote>& quotes,
                           const SviFitSettings& settings = {});

} // namespace cppquant
