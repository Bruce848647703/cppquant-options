#include "cppquant/pde.hpp"
#include "test_framework.hpp"

#include <cmath>
#include <stdexcept>

using namespace cppquant;

namespace {

OptionSpec make_spec(OptionType type, ExerciseStyle style, double spot,
                     double strike, double time, double volatility, double rate,
                     double dividend = 0.0) {
  return OptionSpec{type, style,      spot, strike,
                    time, volatility, rate, dividend};
}

OptionSpec european(const OptionSpec& spec) {
  OptionSpec copy = spec;
  copy.style = ExerciseStyle::European;
  return copy;
}

} // namespace

TEST(pde_european_matches_black_scholes) {
  const OptionSpec cases[] = {
      make_spec(OptionType::Call, ExerciseStyle::European, 100, 100, 1.0, 0.2,
                0.05),
      make_spec(OptionType::Put, ExerciseStyle::European, 100, 100, 1.0, 0.2,
                0.05),
      make_spec(OptionType::Call, ExerciseStyle::European, 100, 110, 0.5, 0.25,
                0.03, 0.02),
      make_spec(OptionType::Put, ExerciseStyle::European, 90, 100, 2.0, 0.3,
                -0.01, 0.01),
  };

  PdeSettings accurate;
  accurate.space_steps = 3'000;
  accurate.time_steps = 3'000;

  for (const auto& spec : cases) {
    const auto result = crank_nicolson_price(spec, accurate);
    CHECK_NEAR(result.price, black_scholes_price(spec), 1e-4);
    CHECK(std::isfinite(result.delta));
    CHECK(std::isfinite(result.gamma));
  }
}

TEST(pde_american_put_matches_binomial) {
  const OptionSpec spec = make_spec(OptionType::Put, ExerciseStyle::American,
                                    100, 100, 1.0, 0.2, 0.05);
  const double tree = binomial_price(spec, 1000);
  const auto result = crank_nicolson_price(spec);

  CHECK_NEAR(result.price, tree, std::max(0.01 * tree, 0.01));
  CHECK(result.price > black_scholes_price(european(spec)) + 0.3);
}

TEST(pde_american_call_without_dividends_matches_european) {
  const OptionSpec spec = make_spec(OptionType::Call, ExerciseStyle::American,
                                    100, 100, 1.0, 0.2, 0.05);
  const double european_price = black_scholes_price(european(spec));
  const auto result = crank_nicolson_price(spec);

  CHECK_NEAR(result.price, european_price,
             std::max(0.01 * european_price, 0.01));
}

TEST(pde_american_call_with_dividends_has_a_premium) {
  const OptionSpec spec = make_spec(OptionType::Call, ExerciseStyle::American,
                                    100, 100, 1.0, 0.2, 0.03, 0.10);
  const double tree = binomial_price(spec, 1000);
  const auto result = crank_nicolson_price(spec);

  CHECK_NEAR(result.price, tree, std::max(0.01 * tree, 0.01));
  CHECK(result.price > black_scholes_price(european(spec)) + 0.5);
}

TEST(pde_converges_with_grid_refinement) {
  const OptionSpec spec = make_spec(OptionType::Call, ExerciseStyle::European,
                                    100, 100, 1.0, 0.2, 0.05);
  const double exact = black_scholes_price(spec);

  PdeSettings coarse;
  coarse.space_steps = 100;
  coarse.time_steps = 100;
  PdeSettings medium = coarse;
  medium.space_steps = 200;
  medium.time_steps = 200;
  PdeSettings fine = coarse;
  fine.space_steps = 400;
  fine.time_steps = 400;

  const double coarse_error =
      std::abs(crank_nicolson_price(spec, coarse).price - exact);
  const double medium_error =
      std::abs(crank_nicolson_price(spec, medium).price - exact);
  const double fine_error =
      std::abs(crank_nicolson_price(spec, fine).price - exact);

  CHECK(medium_error < coarse_error);
  CHECK(fine_error < medium_error);
  CHECK(fine_error < 2e-3);

  PdeSettings very_fine;
  very_fine.space_steps = 2'000;
  very_fine.time_steps = 2'000;
  const double very_fine_error =
      std::abs(crank_nicolson_price(spec, very_fine).price - exact);
  CHECK(very_fine_error < 1e-4);
}

TEST(pde_greeks_match_analytic_values) {
  const OptionSpec spec = make_spec(OptionType::Call, ExerciseStyle::European,
                                    100, 100, 1.0, 0.2, 0.05);
  const Greeks analytic = black_scholes_greeks(spec);
  const auto result = crank_nicolson_price(spec);

  CHECK_NEAR(result.delta, analytic.delta, 1e-3);
  CHECK_NEAR(result.gamma, analytic.gamma, 1e-3);
}

TEST(pde_handles_expiry_and_zero_volatility) {
  const OptionSpec expired = make_spec(OptionType::Put, ExerciseStyle::American,
                                       90, 100, 0.0, 0.2, 0.05);
  CHECK_NEAR(crank_nicolson_price(expired).price, 10.0, 1e-12);

  const OptionSpec deterministic = make_spec(
      OptionType::Put, ExerciseStyle::American, 80, 100, 1.0, 0.0, 0.05);
  CHECK_NEAR(crank_nicolson_price(deterministic).price,
             binomial_price(deterministic, 400), 1e-12);
}

TEST(pde_interpolates_between_grid_nodes) {
  const OptionSpec spec = make_spec(OptionType::Call, ExerciseStyle::European,
                                    100.7, 100, 1.0, 0.2, 0.05);
  const auto result = crank_nicolson_price(spec);
  CHECK_NEAR(result.price, black_scholes_price(spec), 2e-3);
}

TEST(pde_rejects_invalid_settings) {
  const OptionSpec spec = make_spec(OptionType::Call, ExerciseStyle::European,
                                    100, 100, 1.0, 0.2, 0.05);

  PdeSettings settings;
  settings.space_steps = 5;
  CHECK_THROW(crank_nicolson_price(spec, settings), std::invalid_argument);

  settings = PdeSettings{};
  settings.time_steps = 0;
  CHECK_THROW(crank_nicolson_price(spec, settings), std::invalid_argument);

  settings = PdeSettings{};
  settings.width_stddev = 0.0;
  CHECK_THROW(crank_nicolson_price(spec, settings), std::invalid_argument);
}
