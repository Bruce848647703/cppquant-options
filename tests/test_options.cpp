#include "cppquant/options.hpp"
#include "test_framework.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

using namespace cppquant;

namespace {

OptionSpec make_spec(OptionType type, double spot, double strike, double time,
                     double volatility, double rate, double dividend = 0.0,
                     ExerciseStyle style = ExerciseStyle::European) {
  return OptionSpec{type, style,      spot, strike,
                    time, volatility, rate, dividend};
}

OptionSpec with_spot(const OptionSpec& base, double spot) {
  OptionSpec copy = base;
  copy.spot = spot;
  return copy;
}

OptionSpec with_volatility(const OptionSpec& base, double volatility) {
  OptionSpec copy = base;
  copy.volatility = volatility;
  return copy;
}

OptionSpec with_time(const OptionSpec& base, double time) {
  OptionSpec copy = base;
  copy.time_to_expiry = time;
  return copy;
}

OptionSpec with_rate(const OptionSpec& base, double rate) {
  OptionSpec copy = base;
  copy.risk_free_rate = rate;
  return copy;
}

constexpr double kPriceTolerance = 1e-9;

} // namespace

TEST(black_scholes_matches_atm_reference_values) {
  const OptionSpec call = make_spec(OptionType::Call, 100, 100, 1.0, 0.2, 0.05);
  const OptionSpec put = make_spec(OptionType::Put, 100, 100, 1.0, 0.2, 0.05);

  CHECK_NEAR(black_scholes_price(call), 10.450583572186, kPriceTolerance);
  CHECK_NEAR(black_scholes_price(put), 5.573526022257, kPriceTolerance);

  const Greeks call_greeks = black_scholes_greeks(call);
  CHECK_NEAR(call_greeks.delta, 0.636830651176, 1e-9);
  CHECK_NEAR(call_greeks.gamma, 0.018762017346, 1e-9);
  CHECK_NEAR(call_greeks.vega, 37.524034691694, 1e-9);
  CHECK_NEAR(call_greeks.theta, -6.414027546438, 1e-9);
  CHECK_NEAR(call_greeks.rho, 53.232481545376, 1e-9);

  const Greeks put_greeks = black_scholes_greeks(put);
  CHECK_NEAR(put_greeks.delta, -0.363169348824, 1e-9);
  CHECK_NEAR(put_greeks.gamma, 0.018762017346, 1e-9);
  CHECK_NEAR(put_greeks.vega, 37.524034691694, 1e-9);
  CHECK_NEAR(put_greeks.theta, -1.657880423935, 1e-9);
  CHECK_NEAR(put_greeks.rho, -41.890460904695, 1e-9);
}

TEST(black_scholes_matches_dividend_reference_values) {
  const OptionSpec call =
      make_spec(OptionType::Call, 100, 110, 0.5, 0.25, 0.03, 0.02);
  const OptionSpec put =
      make_spec(OptionType::Put, 100, 110, 0.5, 0.25, 0.03, 0.02);

  CHECK_NEAR(black_scholes_price(call), 3.553525293024, kPriceTolerance);
  CHECK_NEAR(black_scholes_price(put), 12.910855274444, kPriceTolerance);

  const Greeks call_greeks = black_scholes_greeks(call);
  CHECK_NEAR(call_greeks.delta, 0.332989587821, 1e-9);
  CHECK_NEAR(call_greeks.gamma, 0.020435395970, 1e-9);
  CHECK_NEAR(call_greeks.vega, 25.544244962323, 1e-9);
  CHECK_NEAR(call_greeks.theta, -6.612445069611, 1e-9);
  CHECK_NEAR(call_greeks.rho, 14.872716744529, 1e-9);

  const Greeks put_greeks = black_scholes_greeks(put);
  CHECK_NEAR(put_greeks.delta, -0.657060245928, 1e-9);
  CHECK_NEAR(put_greeks.theta, -5.341675336419, 1e-9);
  CHECK_NEAR(put_greeks.rho, -39.308439933639, 1e-9);
}

TEST(black_scholes_supports_negative_rates) {
  const OptionSpec call =
      make_spec(OptionType::Call, 100, 90, 2.0, 0.3, -0.01, 0.01);
  const OptionSpec put =
      make_spec(OptionType::Put, 100, 90, 2.0, 0.3, -0.01, 0.01);

  CHECK_NEAR(black_scholes_price(call), 19.232835880878, 1e-9);
  CHECK_NEAR(black_scholes_price(put), 13.031089152611, 1e-9);

  const Greeks greeks = black_scholes_greeks(call);
  CHECK_NEAR(greeks.delta, 0.630157665008, 1e-9);
  CHECK_NEAR(greeks.gamma, 0.008619256755, 1e-9);
  CHECK_NEAR(greeks.vega, 51.715540530682, 1e-9);
  CHECK_NEAR(greeks.theta, -2.810678568594, 1e-9);
  CHECK_NEAR(greeks.rho, 87.565861239844, 1e-9);
}

TEST(put_call_parity_holds_with_a_dividend) {
  const double spot = 100.0;
  const double rate = 0.03;
  const double dividend = 0.02;
  const double time = 0.5;

  for (const double strike : {80.0, 90.0, 100.0, 110.0, 120.0}) {
    const OptionSpec call =
        make_spec(OptionType::Call, spot, strike, time, 0.25, rate, dividend);
    const OptionSpec put =
        make_spec(OptionType::Put, spot, strike, time, 0.25, rate, dividend);

    const double left = black_scholes_price(call) - black_scholes_price(put);
    const double right =
        spot * std::exp(-dividend * time) - strike * std::exp(-rate * time);
    CHECK_NEAR(left, right, 1e-10);
  }
}

TEST(greeks_match_finite_differences) {
  const OptionSpec base =
      make_spec(OptionType::Call, 100, 105, 0.75, 0.28, 0.04, 0.015);
  const Greeks greeks = black_scholes_greeks(base);

  const double spot_step = 0.01;
  const double delta_numeric =
      (black_scholes_price(with_spot(base, base.spot + spot_step)) -
       black_scholes_price(with_spot(base, base.spot - spot_step))) /
      (2.0 * spot_step);
  CHECK_NEAR(greeks.delta, delta_numeric, 1e-5);

  const double gamma_step = 1.0;
  const double gamma_numeric =
      (black_scholes_price(with_spot(base, base.spot + gamma_step)) -
       2.0 * black_scholes_price(base) +
       black_scholes_price(with_spot(base, base.spot - gamma_step))) /
      (gamma_step * gamma_step);
  CHECK_NEAR(greeks.gamma, gamma_numeric, 1e-4);

  const double vol_step = 1e-4;
  const double vega_numeric =
      (black_scholes_price(with_volatility(base, base.volatility + vol_step)) -
       black_scholes_price(with_volatility(base, base.volatility - vol_step))) /
      (2.0 * vol_step);
  CHECK_NEAR(greeks.vega, vega_numeric, 1e-4);

  const double time_step = 1e-4;
  const double theta_numeric =
      -(black_scholes_price(with_time(base, base.time_to_expiry + time_step)) -
        black_scholes_price(with_time(base, base.time_to_expiry - time_step))) /
      (2.0 * time_step);
  CHECK_NEAR(greeks.theta, theta_numeric, 1e-3);

  const double rate_step = 1e-5;
  const double rho_numeric =
      (black_scholes_price(with_rate(base, base.risk_free_rate + rate_step)) -
       black_scholes_price(with_rate(base, base.risk_free_rate - rate_step))) /
      (2.0 * rate_step);
  CHECK_NEAR(greeks.rho, rho_numeric, 1e-3);
}

TEST(greeks_match_finite_differences_for_puts) {
  const OptionSpec base =
      make_spec(OptionType::Put, 100, 95, 1.25, 0.22, 0.02, 0.01);
  const Greeks greeks = black_scholes_greeks(base);

  const double time_step = 1e-4;
  const double theta_numeric =
      -(black_scholes_price(with_time(base, base.time_to_expiry + time_step)) -
        black_scholes_price(with_time(base, base.time_to_expiry - time_step))) /
      (2.0 * time_step);
  CHECK_NEAR(greeks.theta, theta_numeric, 1e-3);

  const double rate_step = 1e-5;
  const double rho_numeric =
      (black_scholes_price(with_rate(base, base.risk_free_rate + rate_step)) -
       black_scholes_price(with_rate(base, base.risk_free_rate - rate_step))) /
      (2.0 * rate_step);
  CHECK_NEAR(greeks.rho, rho_numeric, 1e-3);

  CHECK(greeks.delta < 0.0);
  CHECK(greeks.delta > -1.0);
  CHECK(greeks.gamma > 0.0);
  CHECK(greeks.vega > 0.0);
}

TEST(expired_options_return_intrinsic_value) {
  const OptionSpec call = make_spec(OptionType::Call, 110, 100, 0.0, 0.2, 0.05);
  const OptionSpec put = make_spec(OptionType::Put, 90, 100, 0.0, 0.2, 0.05);

  CHECK_NEAR(black_scholes_price(call), 10.0, 0.0);
  CHECK_NEAR(black_scholes_price(put), 10.0, 0.0);

  const Greeks call_greeks = black_scholes_greeks(call);
  CHECK_NEAR(call_greeks.delta, 1.0, 0.0);
  CHECK_NEAR(call_greeks.gamma, 0.0, 0.0);
  CHECK_NEAR(call_greeks.vega, 0.0, 0.0);

  const Greeks put_greeks = black_scholes_greeks(put);
  CHECK_NEAR(put_greeks.delta, -1.0, 0.0);
}

TEST(zero_volatility_uses_the_discounted_forward) {
  const OptionSpec call =
      make_spec(OptionType::Call, 100, 90, 1.0, 0.0, 0.05, 0.02);
  const OptionSpec put =
      make_spec(OptionType::Put, 100, 90, 1.0, 0.0, 0.05, 0.02);

  const double forward = 100.0 * std::exp((0.05 - 0.02) * 1.0);
  const double discount = std::exp(-0.05);
  CHECK_NEAR(black_scholes_price(call),
             discount * std::max(forward - 90.0, 0.0), 1e-10);
  CHECK_NEAR(black_scholes_price(put), discount * std::max(90.0 - forward, 0.0),
             1e-10);

  const Greeks greeks = black_scholes_greeks(call);
  CHECK_NEAR(greeks.delta, std::exp(-0.02), 1e-12);
  CHECK_NEAR(greeks.gamma, 0.0, 0.0);
  CHECK_NEAR(greeks.vega, 0.0, 0.0);
}

TEST(black_scholes_rejects_invalid_inputs) {
  const double nan = std::numeric_limits<double>::quiet_NaN();

  CHECK_THROW(
      black_scholes_price(make_spec(OptionType::Call, -1, 100, 1, 0.2, 0.05)),
      std::invalid_argument);
  CHECK_THROW(
      black_scholes_price(make_spec(OptionType::Call, 100, 0, 1, 0.2, 0.05)),
      std::invalid_argument);
  CHECK_THROW(
      black_scholes_price(make_spec(OptionType::Call, 100, 100, -1, 0.2, 0.05)),
      std::invalid_argument);
  CHECK_THROW(
      black_scholes_price(make_spec(OptionType::Call, 100, 100, 1, -0.2, 0.05)),
      std::invalid_argument);
  CHECK_THROW(
      black_scholes_price(make_spec(OptionType::Call, nan, 100, 1, 0.2, 0.05)),
      std::invalid_argument);
  CHECK_THROW(
      black_scholes_price(make_spec(OptionType::Call, 100, 100, 1, 0.2, 0.05,
                                    0.0, ExerciseStyle::American)),
      std::invalid_argument);
}

TEST(implied_volatility_round_trips) {
  const double spot = 100.0;
  const double rate = 0.03;
  const double dividend = 0.01;

  for (const auto type : {OptionType::Call, OptionType::Put}) {
    for (const double strike : {80.0, 100.0, 120.0}) {
      for (const double time : {0.1, 1.0, 3.0}) {
        for (const double volatility : {0.05, 0.2, 0.8}) {
          const OptionSpec spec =
              make_spec(type, spot, strike, time, volatility, rate, dividend);
          const double price = black_scholes_price(spec);
          const double dfr = std::exp(-rate * time);
          const double dfq = std::exp(-dividend * time);
          const double lower = type == OptionType::Call
                                   ? std::max(spot * dfq - strike * dfr, 0.0)
                                   : std::max(strike * dfr - spot * dfq, 0.0);
          if (price - lower < 1e-9 * std::max(1.0, spot)) {
            // Deep in-the-money prices sit at the no-arbitrage bound to
            // machine precision; volatility is not identifiable from them.
            continue;
          }
          const double implied = implied_volatility(spec, price);
          CHECK_NEAR(implied, volatility, 1e-7);
        }
      }
    }
  }
}

TEST(implied_volatility_round_trips_with_negative_rates) {
  const OptionSpec spec =
      make_spec(OptionType::Call, 100, 110, 2.0, 0.35, -0.01, 0.02);
  const double price = black_scholes_price(spec);
  CHECK_NEAR(implied_volatility(spec, price), 0.35, 1e-7);
}

TEST(implied_volatility_handles_low_volatility) {
  const OptionSpec spec =
      make_spec(OptionType::Put, 100, 100, 0.5, 0.001, 0.0, 0.0);
  const double price = black_scholes_price(spec);
  const double implied = implied_volatility(spec, price);
  CHECK_NEAR(implied, 0.001, 1e-5);
}

TEST(implied_volatility_rejects_out_of_bounds_prices) {
  const OptionSpec call =
      make_spec(OptionType::Call, 100, 80, 1.0, 0.2, 0.05, 0.0);

  CHECK_THROW(implied_volatility(call, 10.0), std::domain_error);
  CHECK_THROW(implied_volatility(call, 101.0), std::domain_error);
  CHECK_THROW(implied_volatility(call, -1.0), std::invalid_argument);

  const OptionSpec expired =
      make_spec(OptionType::Call, 100, 100, 0.0, 0.2, 0.05, 0.0);
  CHECK_THROW(implied_volatility(expired, 5.0), std::invalid_argument);

  const OptionSpec american = make_spec(OptionType::Call, 100, 100, 1.0, 0.2,
                                        0.05, 0.0, ExerciseStyle::American);
  CHECK_THROW(implied_volatility(american, 10.0), std::invalid_argument);
}

TEST(binomial_european_converges_to_black_scholes) {
  const OptionSpec call =
      make_spec(OptionType::Call, 100, 110, 0.5, 0.25, 0.03, 0.02);
  const OptionSpec put =
      make_spec(OptionType::Put, 100, 110, 0.5, 0.25, 0.03, 0.02);

  CHECK_NEAR(binomial_price(call, 1000), black_scholes_price(call), 5e-3);
  CHECK_NEAR(binomial_price(put, 1000), black_scholes_price(put), 5e-3);
}

TEST(binomial_american_put_exercises_early) {
  const OptionSpec european =
      make_spec(OptionType::Put, 100, 100, 1.0, 0.2, 0.05);
  const OptionSpec american = make_spec(OptionType::Put, 100, 100, 1.0, 0.2,
                                        0.05, 0.0, ExerciseStyle::American);

  const double european_price = black_scholes_price(european);
  const double american_price = binomial_price(american, 1000);

  CHECK(american_price > european_price + 0.4);
  CHECK(american_price < european_price + 0.7);
}

TEST(binomial_american_call_without_dividends_has_no_early_exercise) {
  const OptionSpec european =
      make_spec(OptionType::Call, 100, 100, 1.0, 0.2, 0.05);
  const OptionSpec american = make_spec(OptionType::Call, 100, 100, 1.0, 0.2,
                                        0.05, 0.0, ExerciseStyle::American);

  const double european_price = black_scholes_price(european);
  const double american_price = binomial_price(american, 1000);
  CHECK_NEAR(american_price, european_price, 1e-2);
}

TEST(binomial_american_call_with_dividends_exercises_early) {
  const OptionSpec european =
      make_spec(OptionType::Call, 100, 100, 1.0, 0.2, 0.03, 0.10);
  const OptionSpec american = make_spec(OptionType::Call, 100, 100, 1.0, 0.2,
                                        0.03, 0.10, ExerciseStyle::American);

  const double european_price = black_scholes_price(european);
  const double american_price = binomial_price(american, 1000);
  CHECK(american_price > european_price + 0.6);
}

TEST(binomial_zero_volatility_american_put_scans_exercise_dates) {
  const OptionSpec american = make_spec(OptionType::Put, 80, 100, 1.0, 0.0,
                                        0.05, 0.0, ExerciseStyle::American);

  // The put is worth exercising immediately: 100 - 80 = 20.
  CHECK_NEAR(binomial_price(american, 100), 20.0, 1e-12);
}

TEST(binomial_rejects_invalid_steps) {
  const OptionSpec spec = make_spec(OptionType::Call, 100, 100, 1.0, 0.2, 0.05);
  CHECK_THROW(binomial_price(spec, 0), std::invalid_argument);
  CHECK_THROW(binomial_price(spec, -5), std::invalid_argument);
}

TEST(monte_carlo_matches_black_scholes) {
  const OptionSpec call =
      make_spec(OptionType::Call, 100, 110, 0.5, 0.25, 0.03, 0.02);
  const OptionSpec put =
      make_spec(OptionType::Put, 100, 110, 0.5, 0.25, 0.03, 0.02);

  const auto call_result = monte_carlo_price(call, 400'000, 42);
  const auto put_result = monte_carlo_price(put, 400'000, 42);

  CHECK_NEAR(call_result.price, black_scholes_price(call),
             std::max(4.0 * call_result.standard_error, 0.01));
  CHECK_NEAR(put_result.price, black_scholes_price(put),
             std::max(4.0 * put_result.standard_error, 0.01));
  CHECK(call_result.standard_error > 0.0);
  CHECK(put_result.standard_error > 0.0);
}

TEST(monte_carlo_control_variate_reduces_standard_error) {
  const OptionSpec call = make_spec(OptionType::Call, 100, 100, 1.0, 0.2, 0.05);
  const OptionSpec put = make_spec(OptionType::Put, 100, 100, 1.0, 0.2, 0.05);

  for (const auto& spec : {call, put}) {
    const auto plain = monte_carlo_price(spec, 100'000, 7, false, false);
    const auto controlled = monte_carlo_price(spec, 100'000, 7, false, true);
    CHECK(controlled.standard_error < plain.standard_error);

    const auto antithetic_plain =
        monte_carlo_price(spec, 100'000, 7, true, false);
    const auto antithetic_controlled =
        monte_carlo_price(spec, 100'000, 7, true, true);
    CHECK(antithetic_controlled.standard_error <
          antithetic_plain.standard_error);
  }
}

TEST(monte_carlo_antithetic_variates_reduce_standard_error) {
  const OptionSpec spec = make_spec(OptionType::Call, 100, 100, 1.0, 0.2, 0.05);
  const auto plain = monte_carlo_price(spec, 100'000, 11, false, false);
  const auto antithetic = monte_carlo_price(spec, 100'000, 11, true, false);
  CHECK(antithetic.standard_error < plain.standard_error);
}

TEST(monte_carlo_is_reproducible) {
  const OptionSpec spec =
      make_spec(OptionType::Call, 100, 105, 0.75, 0.3, 0.04, 0.01);
  const auto first = monte_carlo_price(spec, 50'000, 123);
  const auto second = monte_carlo_price(spec, 50'000, 123);
  CHECK_EQ(first.price, second.price);
  CHECK_EQ(first.standard_error, second.standard_error);
}

TEST(monte_carlo_rounds_to_even_samples_for_antithetic) {
  const OptionSpec spec = make_spec(OptionType::Call, 100, 100, 1.0, 0.2, 0.05);
  const auto result = monte_carlo_price(spec, 10'001, 1);
  CHECK_EQ(result.paths, 10'001u);
  CHECK_EQ(result.effective_samples, 10'000u);
}

TEST(monte_carlo_handles_a_two_path_run) {
  const OptionSpec spec =
      make_spec(OptionType::Call, 100, 100, 1.0, 0.2, 0.05);
  const auto result = monte_carlo_price(spec, 2, 1);
  CHECK(std::isfinite(result.price));
  CHECK_EQ(result.effective_samples, 2u);
  CHECK_NEAR(result.standard_error, 0.0, 0.0);
}

TEST(monte_carlo_rejects_unsupported_inputs) {
  const OptionSpec american = make_spec(OptionType::Call, 100, 100, 1.0, 0.2,
                                        0.05, 0.0, ExerciseStyle::American);
  CHECK_THROW(monte_carlo_price(american, 1000, 1), std::invalid_argument);

  const OptionSpec spec = make_spec(OptionType::Call, 100, 100, 1.0, 0.2, 0.05);
  CHECK_THROW(monte_carlo_price(spec, 1, 1), std::invalid_argument);
}

TEST(deep_in_the_money_and_out_of_the_money_prices_stay_finite) {
  for (const double strike : {1.0, 10'000.0}) {
    const OptionSpec call =
        make_spec(OptionType::Call, 100, strike, 1.0, 0.2, 0.05);
    const OptionSpec put =
        make_spec(OptionType::Put, 100, strike, 1.0, 0.2, 0.05);

    for (const auto& spec : {call, put}) {
      const double price = black_scholes_price(spec);
      CHECK(std::isfinite(price));
      CHECK(price >= 0.0);
    }
  }

  const OptionSpec call = make_spec(OptionType::Call, 100, 50, 1.0, 0.2, 0.05);
  const OptionSpec put = make_spec(OptionType::Put, 100, 50, 1.0, 0.2, 0.05);
  const double parity = black_scholes_price(call) - black_scholes_price(put);
  CHECK_NEAR(parity, 100.0 - 50.0 * std::exp(-0.05), 1e-10);
}

TEST(to_string_covers_enums) {
  CHECK_EQ(std::string(to_string(OptionType::Call)), std::string("call"));
  CHECK_EQ(std::string(to_string(OptionType::Put)), std::string("put"));
  CHECK_EQ(std::string(to_string(ExerciseStyle::European)),
           std::string("european"));
  CHECK_EQ(std::string(to_string(ExerciseStyle::American)),
           std::string("american"));
}

TEST(intrinsic_value_is_the_undiscounted_payoff) {
  const OptionSpec call = make_spec(OptionType::Call, 110, 100, 1, 0.2, 0.05);
  const OptionSpec put = make_spec(OptionType::Put, 90, 100, 1, 0.2, 0.05);
  CHECK_NEAR(intrinsic_value(call), 10.0, 0.0);
  CHECK_NEAR(intrinsic_value(put), 10.0, 0.0);
}
