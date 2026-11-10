#pragma once

#include "cppquant/options.hpp"
#include "cppquant/surface.hpp"

#include <cstddef>

namespace cppquant {

struct PdeSettings {
  /// Number of space intervals on [0, S_max]. At least 10.
  std::size_t space_steps{400};
  /// Number of time intervals. At least 1.
  std::size_t time_steps{400};
  /// Grid width in standard deviations of log-spot.
  double width_stddev{6.0};
  /// Floor applied to local variance when it is non-finite or non-positive.
  double local_variance_floor{1e-8};
};

struct PdeResult {
  double price{};
  double delta{};
  double gamma{};
  std::size_t space_steps{};
  std::size_t time_steps{};
  /// Number of grid/time points where local variance had to be floored.
  std::size_t clamped_local_variance_points{};
};

/// Crank-Nicolson finite-difference price for a European or American option.
///
/// The Black-Scholes PDE is solved on a uniform spot grid with Dirichlet
/// boundaries. American exercise is applied by projecting the solution onto
/// the intrinsic value after every time step. The returned delta and gamma are
/// central finite differences of the final grid.
PdeResult crank_nicolson_price(const OptionSpec& spec,
                               const PdeSettings& settings = {});

/// Local-volatility PDE using Dupire local variance from a fitted surface.
///
/// The diffusion term uses `surface.local_variance(spot, tau)` at every grid
/// node and time step. Non-finite or non-positive local variance is replaced by
/// `settings.local_variance_floor` and counted in
/// `PdeResult::clamped_local_variance_points`.
///
/// The drift still comes from the option's risk-free rate and dividend yield,
/// so the surface forward should be consistent with them for the boundary
/// conditions to be accurate. The option expiry must lie inside the surface
/// time range. Below the first slice, `local_variance` extrapolates total
/// variance linearly from zero so the integral matches the surface. The grid
/// is sized from implied volatility, not the larger local volatility in the
/// SVI wings.
PdeResult crank_nicolson_local_vol_price(const OptionSpec& spec,
                                         const VolatilitySurface& surface,
                                         const PdeSettings& settings = {});

} // namespace cppquant
