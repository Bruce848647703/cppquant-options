#pragma once

#include "cppquant/options.hpp"

namespace cppquant::detail {

inline constexpr double kPi = 3.14159265358979323846;

double norm_cdf(double x) noexcept;
double norm_pdf(double x) noexcept;
double forward_price(const OptionSpec& spec) noexcept;
double option_payoff(OptionType type, double spot, double strike) noexcept;

/// Price and Greeks with an explicit volatility. Assumes the spec has been
/// validated, the style is European and time_to_expiry > 0.
double bs_price(const OptionSpec& spec, double volatility);
Greeks bs_greeks(const OptionSpec& spec, double volatility);

/// Value of a deterministic (zero volatility) European or American option by
/// scanning the exercise dates implied by `steps`.
double deterministic_price(const OptionSpec& spec, int steps);

} // namespace cppquant::detail
