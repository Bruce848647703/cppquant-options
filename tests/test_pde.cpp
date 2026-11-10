#include "cppquant/pde.hpp"
#include "test_framework.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

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

namespace {

SurfaceSlice flat_slice(double spot, double rate, double dividend, double time,
                        double volatility) {
  const double forward = spot * std::exp((rate - dividend) * time);
  const SviParameters parameters{volatility * volatility * time, 0.0, 0.0, 0.0,
                                 0.1};
  return SurfaceSlice{SviSlice{forward, time, parameters},
                      std::exp(-rate * time)};
}

VolatilitySurface flat_surface(double spot, double rate, double dividend,
                               double volatility) {
  std::vector<SurfaceSlice> slices;
  for (const double time : {0.25, 0.5, 1.0, 2.0}) {
    slices.push_back(flat_slice(spot, rate, dividend, time, volatility));
  }
  return VolatilitySurface(std::move(slices));
}

VolatilitySurface
term_structure_surface(double spot, double rate, double dividend,
                       const std::vector<std::pair<double, double>>& expiries) {
  std::vector<SurfaceSlice> slices;
  for (const auto& [time, volatility] : expiries) {
    slices.push_back(flat_slice(spot, rate, dividend, time, volatility));
  }
  return VolatilitySurface(std::move(slices));
}

} // namespace

TEST(local_vol_pde_matches_black_scholes_on_a_flat_surface) {
  const double spot = 100.0;
  const double rate = 0.05;
  const double volatility = 0.2;
  const VolatilitySurface surface = flat_surface(spot, rate, 0.0, volatility);

  PdeSettings settings;
  settings.space_steps = 800;
  settings.time_steps = 800;

  for (const auto type : {OptionType::Call, OptionType::Put}) {
    for (const double strike : {80.0, 100.0, 120.0}) {
      const OptionSpec spec = make_spec(type, ExerciseStyle::European, spot,
                                        strike, 1.0, volatility, rate);
      const auto result =
          crank_nicolson_local_vol_price(spec, surface, settings);
      CHECK_NEAR(result.price, black_scholes_price(spec), 2e-3);
      CHECK_EQ(result.clamped_local_variance_points, 0u);
      CHECK(std::isfinite(result.delta));
      CHECK(std::isfinite(result.gamma));
    }
  }
}

TEST(local_vol_pde_integrates_a_term_structure_to_total_variance) {
  const double spot = 100.0;
  const double rate = 0.03;
  const VolatilitySurface surface = term_structure_surface(
      spot, rate, 0.0, {{0.05, 0.15}, {0.5, 0.15}, {1.0, 0.3}});

  const OptionSpec spec = make_spec(OptionType::Call, ExerciseStyle::European,
                                    spot, 100.0, 1.0, 0.3, rate);
  PdeSettings settings;
  settings.space_steps = 800;
  settings.time_steps = 800;

  const auto result = crank_nicolson_local_vol_price(spec, surface, settings);
  CHECK_NEAR(result.price, black_scholes_price(spec), 3e-3);
  CHECK_EQ(result.clamped_local_variance_points, 0u);
}

TEST(local_vol_pde_prices_an_american_put) {
  const double spot = 100.0;
  const double rate = 0.05;
  const double volatility = 0.2;
  const VolatilitySurface surface = flat_surface(spot, rate, 0.0, volatility);
  const OptionSpec spec = make_spec(OptionType::Put, ExerciseStyle::American,
                                    spot, 100.0, 1.0, volatility, rate);

  PdeSettings settings;
  settings.space_steps = 600;
  settings.time_steps = 600;
  const auto result = crank_nicolson_local_vol_price(spec, surface, settings);
  const double tree = binomial_price(spec, 1000);

  CHECK_NEAR(result.price, tree, std::max(0.01 * tree, 0.01));
  CHECK(result.price > black_scholes_price(european(spec)) + 0.3);
}

TEST(local_vol_pde_floors_negative_variance) {
  const double spot = 100.0;
  const double rate = 0.03;
  const VolatilitySurface surface =
      term_structure_surface(spot, rate, 0.0, {{0.5, 0.3}, {1.0, 0.2}});
  CHECK(surface.local_variance(100.0, 0.75) < 0.0);

  const OptionSpec spec = make_spec(OptionType::Call, ExerciseStyle::European,
                                    spot, 100.0, 1.0, 0.3, rate);
  PdeSettings settings;
  settings.space_steps = 200;
  settings.time_steps = 200;
  const auto result = crank_nicolson_local_vol_price(spec, surface, settings);

  CHECK(result.clamped_local_variance_points > 0u);
  CHECK(std::isfinite(result.price));
  CHECK(std::isfinite(result.delta));
  CHECK(std::isfinite(result.gamma));
}

TEST(local_vol_pde_handles_expiry_and_validates_inputs) {
  const VolatilitySurface surface = flat_surface(100.0, 0.05, 0.0, 0.2);
  const OptionSpec expired = make_spec(OptionType::Put, ExerciseStyle::European,
                                       90.0, 100.0, 0.0, 0.2, 0.05);
  CHECK_NEAR(crank_nicolson_local_vol_price(expired, surface).price, 10.0,
             1e-12);

  const OptionSpec too_short = make_spec(
      OptionType::Call, ExerciseStyle::European, 100.0, 100.0, 0.1, 0.2, 0.05);
  CHECK_THROW(crank_nicolson_local_vol_price(too_short, surface),
              std::domain_error);

  const OptionSpec too_long = make_spec(
      OptionType::Call, ExerciseStyle::European, 100.0, 100.0, 2.5, 0.2, 0.05);
  CHECK_THROW(crank_nicolson_local_vol_price(too_long, surface),
              std::domain_error);

  const SviParameters parameters{0.04, 0.0, 0.0, 0.0, 0.1};
  std::vector<SurfaceSlice> single{
      SurfaceSlice{SviSlice{100.0, 1.0, parameters}, 0.95},
  };
  const VolatilitySurface one_slice(std::move(single));
  const OptionSpec spec = make_spec(OptionType::Call, ExerciseStyle::European,
                                    100.0, 100.0, 1.0, 0.2, 0.05);
  CHECK_THROW(crank_nicolson_local_vol_price(spec, one_slice),
              std::domain_error);

  PdeSettings settings;
  settings.local_variance_floor = 0.0;
  CHECK_THROW(crank_nicolson_local_vol_price(spec, surface, settings),
              std::invalid_argument);
}

TEST(local_vol_pde_reproduces_a_skewed_surface_price) {
  std::vector<SurfaceSlice> slices;
  for (const double time : {0.5, 1.0, 2.0}) {
    const SviParameters parameters{0.01 * time, 0.12 * time, -0.4, 0.0, 0.15};
    slices.push_back(SurfaceSlice{SviSlice{100.0, time, parameters}, 1.0});
  }
  const VolatilitySurface surface(std::move(slices));

  PdeSettings settings;
  settings.space_steps = 800;
  settings.time_steps = 800;
  for (const double strike : {90.0, 100.0, 110.0}) {
    const OptionSpec spec = make_spec(OptionType::Call, ExerciseStyle::European,
                                      100.0, strike, 1.0, 0.2, 0.0);
    const double reference = surface.price(strike, 1.0, OptionType::Call);
    const auto result = crank_nicolson_local_vol_price(spec, surface, settings);
    CHECK_NEAR(result.price, reference, std::max(2e-3, 2e-4 * reference));
  }
}

TEST(pde_deep_in_the_money_boundaries_match_black_scholes) {
  PdeSettings settings;
  settings.space_steps = 1200;
  settings.time_steps = 1200;

  const OptionSpec deep_put = make_spec(
      OptionType::Put, ExerciseStyle::European, 60.0, 100.0, 1.0, 0.2, 0.3);
  const OptionSpec deep_call =
      make_spec(OptionType::Call, ExerciseStyle::European, 140.0, 100.0, 1.0,
                0.2, 0.0, 0.3);
  CHECK_NEAR(crank_nicolson_price(deep_put, settings).price,
             black_scholes_price(deep_put), 1e-3);
  CHECK_NEAR(crank_nicolson_price(deep_call, settings).price,
             black_scholes_price(deep_call), 1e-3);
}
