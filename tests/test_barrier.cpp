#include "cppquant/barrier.hpp"
#include "test_framework.hpp"

#include <cmath>
#include <stdexcept>

using namespace cppquant;

namespace {

OptionSpec barrier_spec(OptionType type, double spot = 100.0,
                        double strike = 100.0, double time = 1.0,
                        double volatility = 0.2, double rate = 0.05,
                        double dividend = 0.02) {
  return OptionSpec{
      type,    ExerciseStyle::European, spot, strike, time, volatility, rate,
      dividend};
}

} // namespace

TEST(barrier_approaching_zero_recovers_vanilla) {
  for (const auto type : {OptionType::Call, OptionType::Put}) {
    const OptionSpec spec = barrier_spec(type);
    CHECK_NEAR(down_and_out_price(spec, 1e-8), black_scholes_price(spec), 1e-9);
  }
}

TEST(barrier_analytic_matches_independent_reflection_values) {
  struct Case {
    OptionType type;
    double strike;
    double barrier;
    double volatility;
    double expected;
  };
  const Case cases[] = {
      {OptionType::Call, 100.0, 60.0, 0.2, 9.22700491123071},
      {OptionType::Put, 100.0, 60.0, 0.2, 5.976976773767874},
      {OptionType::Call, 100.0, 80.0, 0.2, 9.133306436496742},
      {OptionType::Put, 100.0, 80.0, 0.2, 1.7326777631695307},
      {OptionType::Call, 100.0, 95.0, 0.2, 4.883524498744151},
      {OptionType::Put, 100.0, 95.0, 0.2, 0.010277121695278453},
      {OptionType::Call, 100.0, 99.0, 0.2, 1.1703614888352933},
      {OptionType::Call, 90.0, 60.0, 0.25, 16.635130216191072},
      {OptionType::Put, 90.0, 60.0, 0.25, 3.060898727852254},
      {OptionType::Call, 90.0, 95.0, 0.25, 6.462628165181819},
  };

  for (const auto& item : cases) {
    const OptionSpec spec =
        barrier_spec(item.type, 100.0, item.strike, 1.0, item.volatility);
    CHECK_NEAR(down_and_out_price(spec, item.barrier), item.expected, 1e-6);
  }
}

TEST(barrier_knockout_cases_are_zero) {
  const OptionSpec call = barrier_spec(OptionType::Call);
  const OptionSpec put = barrier_spec(OptionType::Put);

  CHECK_NEAR(down_and_out_price(call, 100.0), 0.0, 0.0);
  CHECK_NEAR(down_and_out_price(call, 105.0), 0.0, 0.0);
  CHECK_NEAR(down_and_out_price(put, 100.0), 0.0, 0.0);
  CHECK_NEAR(down_and_out_price(put, 105.0), 0.0, 0.0);

  const BarrierGreeks knocked = down_and_out_greeks(call, 105.0);
  CHECK_NEAR(knocked.price, 0.0, 0.0);
  CHECK_NEAR(knocked.delta, 0.0, 0.0);
  CHECK_NEAR(knocked.vega, 0.0, 0.0);
}

TEST(barrier_greeks_are_finite) {
  for (const auto type : {OptionType::Call, OptionType::Put}) {
    const BarrierGreeks greeks = down_and_out_greeks(barrier_spec(type), 80.0);
    CHECK(std::isfinite(greeks.price));
    CHECK(std::isfinite(greeks.delta));
    CHECK(std::isfinite(greeks.vega));
    CHECK(greeks.price > 0.0);
  }
  CHECK(down_and_out_greeks(barrier_spec(OptionType::Call), 80.0).delta > 0.0);
  CHECK(down_and_out_greeks(barrier_spec(OptionType::Put), 60.0).delta < 0.0);
}

TEST(barrier_pde_matches_the_analytic_price) {
  PdeSettings settings;
  settings.space_steps = 800;
  settings.time_steps = 800;

  for (const auto type : {OptionType::Call, OptionType::Put}) {
    for (const double barrier : {60.0, 80.0, 95.0, 99.0}) {
      const OptionSpec spec = barrier_spec(type);
      const double analytic = down_and_out_price(spec, barrier);
      const auto pde =
          crank_nicolson_down_and_out_price(spec, barrier, settings);
      CHECK_NEAR(pde.price, analytic, 2e-3);
      CHECK_EQ(pde.clamped_local_variance_points, 0u);
    }
  }
}

TEST(barrier_pde_recovers_vanilla_at_a_tiny_barrier) {
  PdeSettings settings;
  settings.space_steps = 800;
  settings.time_steps = 800;

  for (const auto type : {OptionType::Call, OptionType::Put}) {
    const OptionSpec spec = barrier_spec(type);
    const auto pde = crank_nicolson_down_and_out_price(spec, 1e-4, settings);
    CHECK_NEAR(pde.price, black_scholes_price(spec), 2e-3);
  }
}

TEST(barrier_handles_expiry_and_zero_volatility) {
  OptionSpec expired = barrier_spec(OptionType::Put, 90.0, 100.0, 0.0, 0.2);
  CHECK_NEAR(down_and_out_price(expired, 80.0), 10.0, 1e-12);
  CHECK_NEAR(down_and_out_price(expired, 95.0), 0.0, 1e-12);

  OptionSpec deterministic =
      barrier_spec(OptionType::Call, 100.0, 100.0, 1.0, 0.0, 0.05, 0.0);
  CHECK_NEAR(down_and_out_price(deterministic, 80.0),
             std::exp(-0.05) * (100.0 * std::exp(0.05) - 100.0), 1e-12);
  CHECK_NEAR(down_and_out_price(deterministic, 103.0), 0.0, 1e-12);
}

TEST(barrier_rejects_invalid_inputs) {
  const OptionSpec call = barrier_spec(OptionType::Call);
  CHECK_THROW(validate_down_and_out(call, 0.0), std::invalid_argument);
  CHECK_THROW(validate_down_and_out(call, -1.0), std::invalid_argument);
  CHECK_THROW(validate_down_and_out(call, std::nan("")), std::invalid_argument);

  OptionSpec american = call;
  american.style = ExerciseStyle::American;
  CHECK_THROW(validate_down_and_out(american, 80.0), std::invalid_argument);

  PdeSettings settings;
  settings.space_steps = 5;
  CHECK_THROW(crank_nicolson_down_and_out_price(call, 80.0, settings),
              std::invalid_argument);
}