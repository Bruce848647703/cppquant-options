#include "cppquant/volatility.hpp"
#include "test_framework.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

using namespace cppquant;

namespace {

SviSlice make_slice(const SviParameters& parameters, double forward = 100.0,
                    double time = 1.0) {
  return SviSlice{forward, time, parameters};
}

std::vector<VolQuote> make_quotes(const SviSlice& slice, double k_min,
                                  double k_max, int count,
                                  double weight = 1.0) {
  std::vector<VolQuote> quotes;
  quotes.reserve(static_cast<std::size_t>(count));
  for (int index = 0; index < count; ++index) {
    const double t =
        static_cast<double>(index) / static_cast<double>(count - 1);
    const double k = k_min + (k_max - k_min) * t;
    const double strike = slice.forward * std::exp(k);
    quotes.push_back(VolQuote{strike, slice.implied_volatility(k), weight});
  }
  return quotes;
}

} // namespace

TEST(svi_total_variance_matches_the_raw_formula) {
  const SviParameters parameters{0.02, 0.12, -0.4, -0.1, 0.15};
  const SviSlice slice = make_slice(parameters, 100.0, 2.0);

  for (const double k : {-1.0, -0.5, -0.1, 0.0, 0.3, 1.0}) {
    const double d = k - parameters.m;
    const double expected =
        parameters.a +
        parameters.b * (parameters.rho * d +
                        std::sqrt(d * d + parameters.sigma * parameters.sigma));
    CHECK_NEAR(slice.total_variance(k), expected, 1e-15);
    CHECK_NEAR(slice.implied_volatility(k),
               std::sqrt(expected / slice.time_to_expiry), 1e-15);
  }
}

TEST(svi_calibration_recovers_a_symmetric_smile) {
  const SviParameters truth{0.015, 0.12, -0.3, 0.0, 0.18};
  const SviSlice truth_slice = make_slice(truth);
  const auto quotes = make_quotes(truth_slice, -1.0, 1.0, 25);

  const SviFitResult result = calibrate_svi(100.0, 1.0, quotes);

  CHECK(result.converged);
  CHECK(result.rmse < 1e-10);
  CHECK(result.max_error < 1e-9);
  CHECK(result.slice.satisfies_basic_constraints());
  CHECK(result.slice.parameters.b > 0.0);
  CHECK(std::abs(result.slice.parameters.rho) < 1.0);
  CHECK_NEAR(result.slice.parameters.b, truth.b, 1e-3);
  CHECK_NEAR(result.slice.parameters.rho, truth.rho, 1e-3);
  CHECK_NEAR(result.slice.parameters.m, truth.m, 1e-3);
  CHECK_NEAR(result.slice.parameters.sigma, truth.sigma, 1e-3);
}

TEST(svi_calibration_recovers_an_asymmetric_smile) {
  const SviParameters truth{0.01, 0.15, 0.45, -0.25, 0.2};
  const SviSlice truth_slice = make_slice(truth);
  const auto quotes = make_quotes(truth_slice, -1.0, 1.0, 31);

  const SviFitResult result = calibrate_svi(100.0, 1.0, quotes);

  CHECK(result.converged);
  CHECK(result.max_error < 1e-9);
  CHECK(result.slice.parameters.rho > 0.0);
  CHECK(result.slice.parameters.m < 0.0);
  CHECK_NEAR(result.slice.parameters.rho, truth.rho, 1e-2);
  CHECK_NEAR(result.slice.parameters.b, truth.b, 1e-2);
}

TEST(svi_calibration_round_trips_each_quote) {
  const SviParameters truth{0.02, 0.1, -0.2, 0.05, 0.12};
  const SviSlice truth_slice = make_slice(truth, 100.0, 0.5);
  const auto quotes = make_quotes(truth_slice, -0.8, 0.8, 21);

  const SviFitResult result = calibrate_svi(100.0, 0.5, quotes);
  for (const auto& quote : quotes) {
    CHECK_NEAR(result.slice.implied_volatility_for_strike(quote.strike),
               quote.implied_volatility, 1e-8);
  }
}

TEST(svi_calibration_is_deterministic) {
  const SviParameters truth{0.018, 0.11, 0.25, 0.1, 0.16};
  const auto quotes = make_quotes(make_slice(truth), -0.9, 0.9, 23);

  const SviFitResult first = calibrate_svi(100.0, 1.0, quotes);
  const SviFitResult second = calibrate_svi(100.0, 1.0, quotes);

  CHECK_EQ(first.rmse, second.rmse);
  CHECK_EQ(first.iterations, second.iterations);
  CHECK_NEAR(first.slice.parameters.a, second.slice.parameters.a, 0.0);
  CHECK_NEAR(first.slice.parameters.b, second.slice.parameters.b, 0.0);
  CHECK_NEAR(first.slice.parameters.rho, second.slice.parameters.rho, 0.0);
  CHECK_NEAR(first.slice.parameters.m, second.slice.parameters.m, 0.0);
  CHECK_NEAR(first.slice.parameters.sigma, second.slice.parameters.sigma, 0.0);
}

TEST(svi_calibration_handles_a_narrow_strike_range) {
  const SviParameters truth{0.02, 0.1, -0.35, 0.0, 0.15};
  const auto quotes = make_quotes(make_slice(truth), -0.1, 0.1, 11);

  const SviFitResult result = calibrate_svi(100.0, 1.0, quotes);
  CHECK(result.converged);
  CHECK(result.max_error < 1e-8);
  CHECK(result.slice.satisfies_basic_constraints());
}

TEST(svi_calibration_rejects_invalid_inputs) {
  const auto quotes = make_quotes(
      make_slice(SviParameters{0.02, 0.1, 0.0, 0.0, 0.15}), -1.0, 1.0, 10);

  CHECK_THROW(calibrate_svi(0.0, 1.0, quotes), std::invalid_argument);
  CHECK_THROW(calibrate_svi(100.0, 0.0, quotes), std::invalid_argument);
  CHECK_THROW(calibrate_svi(100.0, -1.0, quotes), std::invalid_argument);

  std::vector<VolQuote> too_few(quotes.begin(), quotes.begin() + 4);
  CHECK_THROW(calibrate_svi(100.0, 1.0, too_few), std::invalid_argument);

  std::vector<VolQuote> duplicate(5, VolQuote{100.0, 0.2, 1.0});
  CHECK_THROW(calibrate_svi(100.0, 1.0, duplicate), std::invalid_argument);

  auto bad_vol = quotes;
  bad_vol[0].implied_volatility = -0.1;
  CHECK_THROW(calibrate_svi(100.0, 1.0, bad_vol), std::invalid_argument);

  auto bad_weight = quotes;
  bad_weight[0].weight = 0.0;
  CHECK_THROW(calibrate_svi(100.0, 1.0, bad_weight), std::invalid_argument);

  SviFitSettings settings;
  settings.max_iterations = 0;
  CHECK_THROW(calibrate_svi(100.0, 1.0, quotes, settings),
              std::invalid_argument);
}

TEST(svi_basic_constraints_reject_invalid_parameters) {
  const SviSlice valid = make_slice(SviParameters{0.02, 0.1, 0.0, 0.0, 0.15});
  CHECK(valid.satisfies_basic_constraints());

  SviParameters bad = valid.parameters;
  bad.rho = 1.0;
  CHECK(!make_slice(bad).satisfies_basic_constraints());

  bad = valid.parameters;
  bad.b = -0.1;
  CHECK(!make_slice(bad).satisfies_basic_constraints());

  bad = valid.parameters;
  bad.sigma = 0.0;
  CHECK(!make_slice(bad).satisfies_basic_constraints());

  bad = SviParameters{-0.1, 0.1, 0.0, 0.0, 0.15};
  CHECK(!make_slice(bad).satisfies_basic_constraints());
}

TEST(svi_variance_stays_non_negative_inside_the_fitted_region) {
  const SviParameters truth{0.02, 0.12, -0.5, -0.1, 0.2};
  const auto quotes = make_quotes(make_slice(truth), -1.5, 1.5, 31);
  const SviFitResult result = calibrate_svi(100.0, 1.0, quotes);

  for (int index = -50; index <= 50; ++index) {
    const double k = 0.05 * static_cast<double>(index);
    const double variance = result.slice.total_variance(k);
    CHECK(variance >= -1e-12);
    CHECK(std::isfinite(result.slice.implied_volatility(k)));
  }
}
