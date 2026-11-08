#pragma once

#include "cppquant/options.hpp"

#include <cstddef>
#include <cstdint>

namespace cppquant {

struct LongstaffSchwartzSettings {
  std::size_t paths{100'000};
  /// Number of exercise dates after time zero (Bermudan approximation).
  std::size_t exercise_steps{50};
  std::uint64_t seed{42};
  /// Polynomial degree for the continuation regression: 1, 2 or 3.
  int polynomial_degree{2};
  /// Pair each normal draw with its negative to reduce variance.
  bool antithetic{true};
};

struct AmericanMonteCarloResult {
  double price{};
  double standard_error{};
  double european_price{};
  double early_exercise_premium{};
  std::size_t paths{};
  std::size_t effective_paths{};
  std::size_t exercise_steps{};
};

/// Longstaff-Schwartz least-squares Monte Carlo for an American option.
///
/// Paths are simulated exactly under geometric Brownian motion on a uniform
/// exercise grid. At each step, in-the-money paths regress the discounted
/// future cashflow on `1, S/K, (S/K)^2` (or degree 1/3) and exercise when the
/// intrinsic value exceeds the fitted continuation value.
///
/// Throws std::invalid_argument for European options, non-positive expiry,
/// fewer than two paths, zero exercise steps or an unsupported degree.
AmericanMonteCarloResult
longstaff_schwartz_price(const OptionSpec& spec,
                         const LongstaffSchwartzSettings& settings = {});

} // namespace cppquant
