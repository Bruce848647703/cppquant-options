#include "detail.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cppquant {

const char* to_string(OptionType type) noexcept {
  return type == OptionType::Call ? "call" : "put";
}

const char* to_string(ExerciseStyle style) noexcept {
  return style == ExerciseStyle::European ? "european" : "american";
}

void validate_spec(const OptionSpec& spec) {
  if (!std::isfinite(spec.spot) || !std::isfinite(spec.strike) ||
      !std::isfinite(spec.time_to_expiry) || !std::isfinite(spec.volatility) ||
      !std::isfinite(spec.risk_free_rate) ||
      !std::isfinite(spec.dividend_yield)) {
    throw std::invalid_argument("option inputs must be finite");
  }
  if (spec.spot <= 0.0) {
    throw std::invalid_argument("spot must be positive");
  }
  if (spec.strike <= 0.0) {
    throw std::invalid_argument("strike must be positive");
  }
  if (spec.time_to_expiry < 0.0) {
    throw std::invalid_argument("time to expiry must be non-negative");
  }
  if (spec.volatility < 0.0) {
    throw std::invalid_argument("volatility must be non-negative");
  }
}

double intrinsic_value(const OptionSpec& spec) noexcept {
  return detail::option_payoff(spec.type, spec.spot, spec.strike);
}

namespace detail {

double norm_cdf(double x) noexcept {
  return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double norm_pdf(double x) noexcept {
  return std::exp(-0.5 * x * x) / std::sqrt(2.0 * kPi);
}

double forward_price(const OptionSpec& spec) noexcept {
  return spec.spot * std::exp((spec.risk_free_rate - spec.dividend_yield) *
                              spec.time_to_expiry);
}

double option_payoff(OptionType type, double spot, double strike) noexcept {
  return type == OptionType::Call ? std::max(spot - strike, 0.0)
                                  : std::max(strike - spot, 0.0);
}

double bs_price(const OptionSpec& spec, double volatility) {
  const double forward = forward_price(spec);
  const double discount = std::exp(-spec.risk_free_rate * spec.time_to_expiry);

  if (volatility <= 0.0) {
    return discount * option_payoff(spec.type, forward, spec.strike);
  }

  const double sqrt_t = std::sqrt(spec.time_to_expiry);
  const double d1 = (std::log(spec.spot / spec.strike) +
                     (spec.risk_free_rate - spec.dividend_yield +
                      0.5 * volatility * volatility) *
                         spec.time_to_expiry) /
                    (volatility * sqrt_t);
  const double d2 = d1 - volatility * sqrt_t;

  if (spec.type == OptionType::Call) {
    return spec.spot * std::exp(-spec.dividend_yield * spec.time_to_expiry) *
               norm_cdf(d1) -
           spec.strike * discount * norm_cdf(d2);
  }
  return spec.strike * discount * norm_cdf(-d2) -
         spec.spot * std::exp(-spec.dividend_yield * spec.time_to_expiry) *
             norm_cdf(-d1);
}

Greeks bs_greeks(const OptionSpec& spec, double volatility) {
  Greeks greeks;

  const double t = spec.time_to_expiry;
  const double dfr = std::exp(-spec.risk_free_rate * t);
  const double dfq = std::exp(-spec.dividend_yield * t);

  if (volatility <= 0.0) {
    const double forward = forward_price(spec);
    const bool call = spec.type == OptionType::Call;
    const bool itm = call ? forward > spec.strike : forward < spec.strike;
    const bool atm = forward == spec.strike;

    if (atm) {
      greeks.delta = call ? 0.5 * dfq : -0.5 * dfq;
    } else if (itm) {
      greeks.delta = call ? dfq : -dfq;
    }

    const double carry = spec.risk_free_rate * spec.strike * dfr;
    if (call) {
      greeks.theta = itm ? spec.dividend_yield * spec.spot * dfq - carry : 0.0;
      greeks.rho = itm ? spec.strike * t * dfr : 0.0;
    } else {
      greeks.theta = itm ? carry - spec.dividend_yield * spec.spot * dfq : 0.0;
      greeks.rho = itm ? -spec.strike * t * dfr : 0.0;
    }
    return greeks;
  }

  const double sqrt_t = std::sqrt(t);
  const double d1 = (std::log(spec.spot / spec.strike) +
                     (spec.risk_free_rate - spec.dividend_yield +
                      0.5 * volatility * volatility) *
                         t) /
                    (volatility * sqrt_t);
  const double d2 = d1 - volatility * sqrt_t;
  const double pdf = norm_pdf(d1);

  greeks.gamma = dfq * pdf / (spec.spot * volatility * sqrt_t);
  greeks.vega = spec.spot * dfq * pdf * sqrt_t;

  if (spec.type == OptionType::Call) {
    const double nd1 = norm_cdf(d1);
    const double nd2 = norm_cdf(d2);
    greeks.delta = dfq * nd1;
    greeks.theta = -spec.spot * dfq * pdf * volatility / (2.0 * sqrt_t) -
                   spec.risk_free_rate * spec.strike * dfr * nd2 +
                   spec.dividend_yield * spec.spot * dfq * nd1;
    greeks.rho = spec.strike * t * dfr * nd2;
  } else {
    const double nminusd1 = norm_cdf(-d1);
    const double nminusd2 = norm_cdf(-d2);
    greeks.delta = -dfq * nminusd1;
    greeks.theta = -spec.spot * dfq * pdf * volatility / (2.0 * sqrt_t) +
                   spec.risk_free_rate * spec.strike * dfr * nminusd2 -
                   spec.dividend_yield * spec.spot * dfq * nminusd1;
    greeks.rho = -spec.strike * t * dfr * nminusd2;
  }

  return greeks;
}

} // namespace detail

double black_scholes_price(const OptionSpec& spec) {
  validate_spec(spec);
  if (spec.style != ExerciseStyle::European) {
    throw std::invalid_argument("Black-Scholes prices European options only");
  }
  if (spec.time_to_expiry == 0.0) {
    return intrinsic_value(spec);
  }
  return detail::bs_price(spec, spec.volatility);
}

Greeks black_scholes_greeks(const OptionSpec& spec) {
  validate_spec(spec);
  if (spec.style != ExerciseStyle::European) {
    throw std::invalid_argument(
        "Black-Scholes Greeks are defined for European options only");
  }

  if (spec.time_to_expiry == 0.0) {
    Greeks greeks;
    const bool call = spec.type == OptionType::Call;
    if (call) {
      greeks.delta = spec.spot > spec.strike
                         ? 1.0
                         : (spec.spot == spec.strike ? 0.5 : 0.0);
    } else {
      greeks.delta = spec.spot < spec.strike
                         ? -1.0
                         : (spec.spot == spec.strike ? -0.5 : 0.0);
    }
    return greeks;
  }

  return detail::bs_greeks(spec, spec.volatility);
}

} // namespace cppquant
