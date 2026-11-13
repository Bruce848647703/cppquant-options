#include "cppquant/barrier.hpp"

#include "detail.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cppquant {
namespace {

double european_price(const OptionSpec& spec, OptionType type, double spot,
                      double strike, double volatility) {
  OptionSpec adjusted = spec;
  adjusted.type = type;
  adjusted.style = ExerciseStyle::European;
  adjusted.spot = spot;
  adjusted.strike = strike;
  return detail::bs_price(adjusted, volatility);
}

double raw_price(const OptionSpec& spec, double spot, double volatility,
                 double barrier) {
  const double time = spec.time_to_expiry;
  const double rate = spec.risk_free_rate;
  const double dividend = spec.dividend_yield;
  const double sigma = volatility;
  const double lambda =
      (rate - dividend + 0.5 * sigma * sigma) / (sigma * sigma);
  const double sqrt_time = std::sqrt(time);
  const double step = sigma * sqrt_time;
  const double ratio = barrier / spot;
  const double threshold = spec.strike;

  if (spec.type == OptionType::Call) {
    if (barrier <= threshold) {
      const double reflected = barrier * barrier / spot;
      return european_price(spec, OptionType::Call, spot, threshold, sigma) -
             std::pow(ratio, 2.0 * lambda - 2.0) *
                 european_price(spec, OptionType::Call, reflected, threshold,
                                sigma);
    }

    const double x1 = std::log(spot / barrier) / step + lambda * step;
    const double y1 = std::log(barrier / spot) / step + lambda * step;
    return spot * std::exp(-dividend * time) * detail::norm_cdf(x1) -
           threshold * std::exp(-rate * time) * detail::norm_cdf(x1 - step) -
           spot * std::exp(-dividend * time) * std::pow(ratio, 2.0 * lambda) *
               detail::norm_cdf(y1) +
           threshold * std::exp(-rate * time) *
               std::pow(ratio, 2.0 * lambda - 2.0) *
               detail::norm_cdf(y1 - step);
  }

  // A down-and-out put with the barrier at or above the strike is knocked out
  // before it can finish in the money.
  if (barrier >= threshold) {
    return 0.0;
  }

  const double drift = rate - dividend - 0.5 * sigma * sigma;
  const double b = std::log(barrier / spot);
  const double reflected = barrier * barrier / spot;
  const double coefficient = std::pow(ratio, 2.0 * lambda - 2.0);
  const double db = (b - drift * time) / step;
  const double dl = (-b - drift * time) / step;

  auto put_tail = [&](double level, double d) {
    return threshold * std::exp(-rate * time) * detail::norm_cdf(d) -
           level * std::exp(-dividend * time) * detail::norm_cdf(d - step);
  };

  const double first =
      european_price(spec, OptionType::Put, spot, threshold, sigma) -
      put_tail(spot, db);
  const double second =
      european_price(spec, OptionType::Put, reflected, threshold, sigma) -
      put_tail(reflected, dl);
  return first - coefficient * second;
}

double zero_volatility_price(const OptionSpec& spec, double barrier) {
  const double drift = spec.risk_free_rate - spec.dividend_yield;
  const double terminal = spec.spot * std::exp(drift * spec.time_to_expiry);
  if (drift < 0.0 && terminal <= barrier) {
    return 0.0;
  }
  return std::exp(-spec.risk_free_rate * spec.time_to_expiry) *
         detail::option_payoff(spec.type, terminal, spec.strike);
}

} // namespace

void validate_down_and_out(const OptionSpec& spec, double barrier) {
  validate_spec(spec);
  if (spec.style != ExerciseStyle::European) {
    throw std::invalid_argument(
        "down-and-out barrier options price European exercise only");
  }
  if (!std::isfinite(barrier) || barrier <= 0.0) {
    throw std::invalid_argument("barrier must be positive and finite");
  }
}

double down_and_out_price(const OptionSpec& spec, double barrier) {
  validate_down_and_out(spec, barrier);
  if (spec.spot <= barrier) {
    return 0.0;
  }
  if (spec.time_to_expiry == 0.0) {
    return intrinsic_value(spec);
  }
  if (spec.volatility == 0.0) {
    return zero_volatility_price(spec, barrier);
  }
  return std::max(0.0, raw_price(spec, spec.spot, spec.volatility, barrier));
}

BarrierGreeks down_and_out_greeks(const OptionSpec& spec, double barrier) {
  validate_down_and_out(spec, barrier);
  BarrierGreeks result;
  if (spec.spot <= barrier) {
    return result;
  }

  result.price = down_and_out_price(spec, barrier);
  if (spec.time_to_expiry == 0.0) {
    const bool call = spec.type == OptionType::Call;
    const bool in_the_money =
        call ? spec.spot > spec.strike : spec.spot < spec.strike;
    result.delta = in_the_money ? (call ? 1.0 : -1.0) : 0.0;
    return result;
  }

  const double spot_step = std::max(1e-4, 1e-4 * spec.spot);
  const double spot_up = spec.spot + spot_step;
  const double spot_down = spec.spot - spot_step;
  const double price_up =
      down_and_out_price(OptionSpec{spec.type, spec.style, spot_up, spec.strike,
                                    spec.time_to_expiry, spec.volatility,
                                    spec.risk_free_rate, spec.dividend_yield},
                         barrier);
  if (spec.spot - spot_step <= barrier) {
    result.delta = (price_up - result.price) / spot_step;
  } else {
    const double price_down = down_and_out_price(
        OptionSpec{spec.type, spec.style, spot_down, spec.strike,
                   spec.time_to_expiry, spec.volatility, spec.risk_free_rate,
                   spec.dividend_yield},
        barrier);
    result.delta = (price_up - price_down) / (2.0 * spot_step);
  }

  if (spec.volatility > 0.0) {
    const double vol_step = std::max(1e-4, 0.01 * spec.volatility);
    const double vol_up = spec.volatility + vol_step;
    const double vol_down = std::max(0.0, spec.volatility - vol_step);
    OptionSpec bumped = spec;
    bumped.volatility = vol_up;
    const double up_price = down_and_out_price(bumped, barrier);
    bumped.volatility = vol_down;
    const double down_price = down_and_out_price(bumped, barrier);
    result.vega = (up_price - down_price) / (vol_up - vol_down);
  }
  return result;
}

} // namespace cppquant