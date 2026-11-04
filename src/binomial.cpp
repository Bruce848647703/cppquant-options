#include "detail.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace cppquant::detail {

double deterministic_price(const OptionSpec& spec, int steps) {
  const double t = spec.time_to_expiry;
  if (t == 0.0) {
    return option_payoff(spec.type, spec.spot, spec.strike);
  }

  if (spec.style == ExerciseStyle::European) {
    return std::exp(-spec.risk_free_rate * t) *
           option_payoff(spec.type, forward_price(spec), spec.strike);
  }

  double best = 0.0;
  for (int index = 0; index <= steps; ++index) {
    const double tau =
        t * static_cast<double>(index) / static_cast<double>(steps);
    const double spot =
        spec.spot * std::exp((spec.risk_free_rate - spec.dividend_yield) * tau);
    const double value = std::exp(-spec.risk_free_rate * tau) *
                         option_payoff(spec.type, spot, spec.strike);
    best = std::max(best, value);
  }
  return best;
}

} // namespace cppquant::detail

namespace cppquant {

double binomial_price(const OptionSpec& spec, int steps) {
  validate_spec(spec);
  if (steps < 1) {
    throw std::invalid_argument("binomial steps must be at least 1");
  }
  if (spec.time_to_expiry == 0.0) {
    return intrinsic_value(spec);
  }
  if (spec.volatility == 0.0) {
    return detail::deterministic_price(spec, steps);
  }

  const double t = spec.time_to_expiry;
  const double dt = t / static_cast<double>(steps);
  const double sqrt_dt = std::sqrt(dt);
  const double up = std::exp(spec.volatility * sqrt_dt);
  const double down = 1.0 / up;
  const double growth =
      std::exp((spec.risk_free_rate - spec.dividend_yield) * dt);
  const double probability = (growth - down) / (up - down);

  if (!std::isfinite(probability) || probability < 0.0 || probability > 1.0) {
    throw std::domain_error("CRR risk-neutral probability is outside [0, 1]; "
                            "choose a different model or time step");
  }

  const double discount = std::exp(-spec.risk_free_rate * dt);
  const bool american = spec.style == ExerciseStyle::American;
  const double up_squared = up * up;

  std::vector<double> values(static_cast<std::size_t>(steps) + 1);

  // Terminal layer. At step `steps`, node j has j up moves and (steps - j)
  // down moves, so spot = S * up^j * down^(steps-j).
  double base = spec.spot * std::pow(down, steps);
  double node_spot = base;
  for (int j = 0; j <= steps; ++j) {
    values[static_cast<std::size_t>(j)] =
        detail::option_payoff(spec.type, node_spot, spec.strike);
    node_spot *= up_squared;
  }

  for (int step = steps - 1; step >= 0; --step) {
    base *= up;
    node_spot = base;
    for (int j = 0; j <= step; ++j) {
      const auto index = static_cast<std::size_t>(j);
      double value = discount * (probability * values[index + 1] +
                                 (1.0 - probability) * values[index]);
      if (american) {
        value = std::max(
            value, detail::option_payoff(spec.type, node_spot, spec.strike));
      }
      values[index] = value;
      node_spot *= up_squared;
    }
  }

  return values[0];
}

} // namespace cppquant
