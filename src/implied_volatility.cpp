#include "detail.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cppquant {
namespace {

double no_arbitrage_lower(const OptionSpec& spec, double dfr, double dfq) {
  if (spec.type == OptionType::Call) {
    return std::max(spec.spot * dfq - spec.strike * dfr, 0.0);
  }
  return std::max(spec.strike * dfr - spec.spot * dfq, 0.0);
}

double no_arbitrage_upper(const OptionSpec& spec, double dfr, double dfq) {
  return spec.type == OptionType::Call ? spec.spot * dfq : spec.strike * dfr;
}

} // namespace

double implied_volatility(const OptionSpec& spec, double market_price,
                          const SolverSettings& settings) {
  validate_spec(spec);
  if (spec.style != ExerciseStyle::European) {
    throw std::invalid_argument(
        "implied_volatility supports European options only");
  }
  if (spec.time_to_expiry <= 0.0) {
    throw std::invalid_argument(
        "cannot imply volatility for an expired option");
  }
  if (!std::isfinite(market_price) || market_price < 0.0) {
    throw std::invalid_argument("market price must be finite and non-negative");
  }
  if (!std::isfinite(settings.tolerance) || settings.tolerance <= 0.0 ||
      settings.max_iterations < 1 || !std::isfinite(settings.max_volatility) ||
      settings.max_volatility <= 0.0) {
    throw std::invalid_argument("invalid implied volatility solver settings");
  }

  const double t = spec.time_to_expiry;
  const double dfq = std::exp(-spec.dividend_yield * t);
  const double dfr = std::exp(-spec.risk_free_rate * t);
  const double lower = no_arbitrage_lower(spec, dfr, dfq);
  const double upper = no_arbitrage_upper(spec, dfr, dfq);
  const double scale = std::max({1.0, lower, upper});
  const double slack = 1e-12 * scale;

  if (market_price < lower - slack || market_price > upper + slack) {
    throw std::domain_error("market price is outside the no-arbitrage bounds");
  }
  if (market_price <= lower + slack) {
    return 0.0;
  }
  if (market_price >= upper - slack) {
    throw std::domain_error("market price implies infinite volatility");
  }

  double lo = 0.0;
  double hi = settings.max_volatility;
  if (detail::bs_price(spec, hi) < market_price) {
    throw std::domain_error(
        "market price requires volatility above the solver maximum");
  }

  double volatility = std::clamp(0.2, lo, hi);
  for (int iteration = 0; iteration < settings.max_iterations; ++iteration) {
    const double difference = detail::bs_price(spec, volatility) - market_price;
    if (std::fabs(difference) <= settings.tolerance) {
      return volatility;
    }

    if (difference > 0.0) {
      hi = volatility;
    } else {
      lo = volatility;
    }
    if (hi - lo <= 1e-12) {
      break;
    }

    const double vega = detail::bs_greeks(spec, volatility).vega;
    double candidate =
        vega > 1e-12 ? volatility - difference / vega : 0.5 * (lo + hi);
    if (!std::isfinite(candidate) || candidate <= lo || candidate >= hi) {
      candidate = 0.5 * (lo + hi);
    }
    volatility = candidate;
  }

  return 0.5 * (lo + hi);
}

} // namespace cppquant
