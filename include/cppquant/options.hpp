#pragma once

#include <cstddef>
#include <cstdint>

namespace cppquant {

enum class OptionType : std::uint8_t { Call, Put };

enum class ExerciseStyle : std::uint8_t { European, American };

/// A single European or American option under the Black-Scholes-Merton model.
///
/// Times are in years, rates and volatility are annualized decimals (0.05 =
/// 5%), and `dividend_yield` is a continuous yield. Prices and Greeks are in
/// the same units as `spot` and `strike`.
struct OptionSpec {
  OptionType type{OptionType::Call};
  ExerciseStyle style{ExerciseStyle::European};
  double spot{};
  double strike{};
  double time_to_expiry{};
  double volatility{};
  double risk_free_rate{};
  double dividend_yield{};
};

/// Analytic sensitivities of a European option.
///
/// `vega` is per 1.00 change in volatility, `rho` per 1.00 change in the risk
/// free rate, and `theta` per year. Divide vega and rho by 100 for the common
/// "per 1%" convention.
struct Greeks {
  double delta{};
  double gamma{};
  double vega{};
  double theta{};
  double rho{};
};

struct SolverSettings {
  /// Price tolerance used while inverting Black-Scholes.
  double tolerance{1e-10};
  int max_iterations{100};
  double min_volatility{1e-8};
  double max_volatility{5.0};
};

struct MonteCarloResult {
  double price{};
  double standard_error{};
  /// Number of paths requested.
  std::size_t paths{};
  /// Number of independent samples actually used (paths rounded down to an
  /// even number when antithetic variates are enabled).
  std::size_t effective_samples{};
};

const char* to_string(OptionType type) noexcept;
const char* to_string(ExerciseStyle style) noexcept;

/// Throws std::invalid_argument when a field is out of range. Rates may be
/// negative; spot and strike must be positive, time and volatility
/// non-negative.
void validate_spec(const OptionSpec& spec);

/// Undiscounted exercise value at the current spot.
double intrinsic_value(const OptionSpec& spec) noexcept;

/// European Black-Scholes-Merton price. Throws when `spec.style` is American.
double black_scholes_price(const OptionSpec& spec);

/// Analytic European Greeks. Throws when `spec.style` is American.
Greeks black_scholes_greeks(const OptionSpec& spec);

/// Invert a European option price. `spec.volatility` is ignored. Throws
/// std::domain_error when the price is outside the no-arbitrage bounds or
/// std::invalid_argument for a non-European or expired option.
double implied_volatility(const OptionSpec& spec, double market_price,
                          const SolverSettings& settings = {});

/// Cox-Ross-Rubinstein binomial price. Supports European and American
/// exercise, including early exercise with a continuous dividend yield.
/// Throws when steps < 1.
double binomial_price(const OptionSpec& spec, int steps);

/// Monte Carlo price for European options under risk-neutral geometric
/// Brownian motion. Antithetic variates and a terminal-stock control variate
/// are enabled by default. Throws for American options.
MonteCarloResult monte_carlo_price(const OptionSpec& spec, std::size_t paths,
                                   std::uint64_t seed, bool antithetic = true,
                                   bool control_variate = true);

} // namespace cppquant
