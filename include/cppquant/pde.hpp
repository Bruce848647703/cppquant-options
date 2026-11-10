#pragma once

#include "cppquant/options.hpp"

#include <cstddef>

namespace cppquant {

struct PdeSettings {
  /// Number of space intervals on [0, S_max]. At least 10.
  std::size_t space_steps{400};
  /// Number of time intervals. At least 1.
  std::size_t time_steps{400};
  /// Grid width in standard deviations of log-spot.
  double width_stddev{6.0};
};

struct PdeResult {
  double price{};
  double delta{};
  double gamma{};
  std::size_t space_steps{};
  std::size_t time_steps{};
};

/// Crank-Nicolson finite-difference price for a European or American option.
///
/// The Black-Scholes PDE is solved on a uniform spot grid with Dirichlet
/// boundaries. American exercise is applied by projecting the solution onto
/// the intrinsic value after every time step. The returned delta and gamma are
/// central finite differences of the final grid.
PdeResult crank_nicolson_price(const OptionSpec& spec,
                               const PdeSettings& settings = {});

} // namespace cppquant
