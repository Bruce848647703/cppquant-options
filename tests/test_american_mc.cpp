#include "cppquant/american_mc.hpp"
#include "test_framework.hpp"

#include <cmath>
#include <stdexcept>

using namespace cppquant;

namespace {

OptionSpec american(OptionType type, double spot, double strike, double time,
                    double volatility, double rate, double dividend = 0.0) {
  return OptionSpec{
      type,    ExerciseStyle::American, spot, strike, time, volatility, rate,
      dividend};
}

OptionSpec european(const OptionSpec& spec) {
  OptionSpec copy = spec;
  copy.style = ExerciseStyle::European;
  return copy;
}

LongstaffSchwartzSettings settings(std::size_t paths = 200'000,
                                   std::size_t steps = 50,
                                   std::uint64_t seed = 42) {
  LongstaffSchwartzSettings value;
  value.paths = paths;
  value.exercise_steps = steps;
  value.seed = seed;
  value.antithetic = true;
  value.polynomial_degree = 2;
  return value;
}

} // namespace

TEST(lsm_american_put_matches_binomial) {
  const OptionSpec spec = american(OptionType::Put, 100, 100, 1.0, 0.2, 0.05);
  const double tree = binomial_price(spec, 1000);
  const auto result = longstaff_schwartz_price(spec, settings());

  const double tolerance = std::max(0.01 * tree, 4.0 * result.standard_error);
  CHECK_NEAR(result.price, tree, tolerance);
  CHECK(result.early_exercise_premium > 0.2);
  CHECK(result.standard_error > 0.0);
  CHECK_NEAR(result.european_price, black_scholes_price(european(spec)), 1e-12);
}

TEST(lsm_american_call_without_dividends_matches_european) {
  const OptionSpec spec = american(OptionType::Call, 100, 100, 1.0, 0.2, 0.05);
  const double european_price = black_scholes_price(european(spec));
  const auto result = longstaff_schwartz_price(spec, settings());

  CHECK_NEAR(result.price, european_price,
             std::max(0.01 * european_price, 4.0 * result.standard_error));
  CHECK_NEAR(result.early_exercise_premium, 0.0,
             std::max(0.01 * european_price, 4.0 * result.standard_error));
}

TEST(lsm_american_call_with_dividends_has_a_premium) {
  const OptionSpec spec =
      american(OptionType::Call, 100, 100, 1.0, 0.2, 0.03, 0.10);
  const double tree = binomial_price(spec, 1000);
  const auto result = longstaff_schwartz_price(spec, settings());

  CHECK_NEAR(result.price, tree,
             std::max(0.01 * tree, 4.0 * result.standard_error));
  CHECK(result.early_exercise_premium > 0.5);
}

TEST(lsm_is_reproducible) {
  const OptionSpec spec = american(OptionType::Put, 100, 100, 1.0, 0.2, 0.05);
  const auto first = longstaff_schwartz_price(spec, settings(50'000, 30, 7));
  const auto second = longstaff_schwartz_price(spec, settings(50'000, 30, 7));

  CHECK_EQ(first.price, second.price);
  CHECK_EQ(first.standard_error, second.standard_error);
  CHECK_EQ(first.effective_paths, second.effective_paths);
}

TEST(lsm_rounds_to_even_paths_for_antithetic) {
  const OptionSpec spec = american(OptionType::Put, 100, 100, 1.0, 0.2, 0.05);
  const auto result = longstaff_schwartz_price(spec, settings(10'001, 20, 1));

  CHECK_EQ(result.paths, 10'001u);
  CHECK_EQ(result.effective_paths, 10'000u);
}

TEST(lsm_handles_zero_volatility) {
  const OptionSpec spec = american(OptionType::Put, 80, 100, 1.0, 0.0, 0.05);
  const auto result = longstaff_schwartz_price(spec, settings(10'000, 100));

  CHECK_NEAR(result.price, 20.0, 1e-12);
  CHECK_NEAR(result.standard_error, 0.0, 0.0);
  CHECK_NEAR(result.european_price, black_scholes_price(european(spec)), 1e-12);
}

TEST(lsm_handles_immediate_exercise) {
  const OptionSpec spec = american(OptionType::Put, 50, 100, 1.0, 0.2, 0.05);
  const auto result = longstaff_schwartz_price(spec, settings(50'000, 30));

  CHECK(result.price >= 50.0 - 1e-9);
  CHECK(result.early_exercise_premium > 4.0);
}

TEST(lsm_handles_expiry) {
  const OptionSpec spec = american(OptionType::Put, 90, 100, 0.0, 0.2, 0.05);
  const auto result = longstaff_schwartz_price(spec, settings());

  CHECK_NEAR(result.price, 10.0, 1e-12);
  CHECK_NEAR(result.european_price, 10.0, 1e-12);
  CHECK_NEAR(result.early_exercise_premium, 0.0, 0.0);
}

TEST(lsm_rejects_invalid_inputs) {
  const OptionSpec put = american(OptionType::Put, 100, 100, 1.0, 0.2, 0.05);
  CHECK_THROW(longstaff_schwartz_price(european(put), settings()),
              std::invalid_argument);

  auto bad = settings();
  bad.paths = 1;
  CHECK_THROW(longstaff_schwartz_price(put, bad), std::invalid_argument);

  bad = settings();
  bad.exercise_steps = 0;
  CHECK_THROW(longstaff_schwartz_price(put, bad), std::invalid_argument);

  bad = settings();
  bad.polynomial_degree = 0;
  CHECK_THROW(longstaff_schwartz_price(put, bad), std::invalid_argument);

  bad = settings();
  bad.polynomial_degree = 4;
  CHECK_THROW(longstaff_schwartz_price(put, bad), std::invalid_argument);
}
