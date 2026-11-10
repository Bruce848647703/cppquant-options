#include "cppquant/surface.hpp"
#include "test_framework.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

using namespace cppquant;

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

} // namespace

TEST(flat_surface_reproduces_black_scholes) {
  const double spot = 100.0;
  const double rate = 0.05;
  const double dividend = 0.02;
  const double volatility = 0.25;
  const VolatilitySurface surface =
      flat_surface(spot, rate, dividend, volatility);

  for (const double time : {0.25, 0.5, 1.0, 2.0}) {
    for (const double strike : {70.0, 85.0, 100.0, 115.0, 130.0}) {
      for (const auto type : {OptionType::Call, OptionType::Put}) {
        const OptionSpec option{type, ExerciseStyle::European,
                                spot, strike,
                                time, volatility,
                                rate, dividend};
        CHECK_NEAR(surface.price(strike, time, type),
                   black_scholes_price(option), 1e-10);
        CHECK_NEAR(surface.price(option), black_scholes_price(option), 1e-10);
        CHECK_NEAR(surface.implied_volatility(strike, time), volatility, 1e-12);
      }
    }
  }
}

TEST(surface_interpolates_total_variance_in_time) {
  const double forward = 100.0;
  const SviParameters low{0.04 * 0.5, 0.0, 0.0, 0.0, 0.1};
  const SviParameters high{0.09 * 1.5, 0.0, 0.0, 0.0, 0.1};
  std::vector<SurfaceSlice> slices{
      SurfaceSlice{SviSlice{forward, 0.5, low}, 0.98},
      SurfaceSlice{SviSlice{forward, 1.5, high}, 0.95},
  };
  const VolatilitySurface surface(std::move(slices));

  const double expected = 0.5 * (0.04 * 0.5 + 0.09 * 1.5);
  CHECK_NEAR(surface.total_variance(0.0, 1.0), expected, 1e-14);
  CHECK_NEAR(surface.implied_volatility(forward, 1.0),
             std::sqrt(expected / 1.0), 1e-14);
}

TEST(surface_interpolates_forward_and_discount_factor) {
  const SviParameters parameters{0.04, 0.0, 0.0, 0.0, 0.1};
  std::vector<SurfaceSlice> slices{
      SurfaceSlice{SviSlice{105.0, 1.0, parameters}, 0.95},
      SurfaceSlice{SviSlice{110.0, 2.0, parameters}, 0.90},
  };
  const VolatilitySurface surface(std::move(slices));

  CHECK_NEAR(surface.forward(1.5), std::sqrt(105.0 * 110.0), 1e-12);
  CHECK_NEAR(surface.discount_factor(1.5), std::sqrt(0.95 * 0.90), 1e-12);
  CHECK_NEAR(surface.forward(1.0), 105.0, 1e-12);
  CHECK_NEAR(surface.discount_factor(2.0), 0.90, 1e-12);
}

TEST(surface_obeys_forward_put_call_parity) {
  const VolatilitySurface surface = flat_surface(100.0, 0.04, 0.01, 0.3);

  for (const double time : {0.25, 0.5, 1.0, 2.0}) {
    for (const double strike : {80.0, 95.0, 100.0, 110.0}) {
      const double call = surface.price(strike, time, OptionType::Call);
      const double put = surface.price(strike, time, OptionType::Put);
      const double expected =
          surface.discount_factor(time) * (surface.forward(time) - strike);
      CHECK_NEAR(call - put, expected, 1e-10);
    }
  }
}

TEST(surface_detects_calendar_arbitrage) {
  const SviParameters low{0.04, 0.0, 0.0, 0.0, 0.1};
  const SviParameters decreasing{0.01, 0.0, 0.0, 0.0, 0.1};

  std::vector<SurfaceSlice> bad{
      SurfaceSlice{SviSlice{100.0, 1.0, low}, 0.95},
      SurfaceSlice{SviSlice{100.0, 2.0, decreasing}, 0.90},
  };
  const VolatilitySurface bad_surface(std::move(bad));
  CHECK(!bad_surface.is_calendar_arbitrage_free());
  CHECK(!bad_surface.is_arbitrage_free());

  std::vector<SurfaceSlice> good{
      SurfaceSlice{SviSlice{100.0, 1.0, low}, 0.95},
      SurfaceSlice{
          SviSlice{100.0, 2.0, SviParameters{0.09, 0.0, 0.0, 0.0, 0.1}}, 0.90},
  };
  const VolatilitySurface good_surface(std::move(good));
  CHECK(good_surface.is_calendar_arbitrage_free());
}

TEST(surface_detects_butterfly_arbitrage) {
  const SviParameters bad{0.01, 0.5, -0.9, 0.0, 0.01};
  std::vector<SurfaceSlice> slices{
      SurfaceSlice{SviSlice{100.0, 1.0, bad}, 0.95},
  };
  const VolatilitySurface surface(std::move(slices));

  CHECK(!surface.is_butterfly_arbitrage_free());
  CHECK(!surface.is_arbitrage_free());

  const VolatilitySurface flat = flat_surface(100.0, 0.05, 0.0, 0.2);
  CHECK(flat.is_butterfly_arbitrage_free());
  CHECK(flat.butterfly_g(0.0, 1.0) > 0.0);
}

TEST(surface_sorts_slices_and_validates_inputs) {
  const SviParameters parameters{0.04, 0.0, 0.0, 0.0, 0.1};
  std::vector<SurfaceSlice> slices{
      SurfaceSlice{SviSlice{105.0, 2.0, parameters}, 0.90},
      SurfaceSlice{SviSlice{100.0, 1.0, parameters}, 0.95},
  };
  const VolatilitySurface surface(std::move(slices));
  CHECK_NEAR(surface.minimum_time(), 1.0, 0.0);
  CHECK_NEAR(surface.maximum_time(), 2.0, 0.0);
  CHECK_NEAR(surface.forward(1.0), 100.0, 1e-12);

  CHECK_THROW(VolatilitySurface(std::vector<SurfaceSlice>{}),
              std::invalid_argument);

  std::vector<SurfaceSlice> duplicate{
      SurfaceSlice{SviSlice{100.0, 1.0, parameters}, 0.95},
      SurfaceSlice{SviSlice{100.0, 1.0 + 1e-14, parameters}, 0.95},
  };
  CHECK_THROW(VolatilitySurface(std::move(duplicate)), std::invalid_argument);

  std::vector<SurfaceSlice> bad_forward{
      SurfaceSlice{SviSlice{0.0, 1.0, parameters}, 0.95},
  };
  CHECK_THROW(VolatilitySurface(std::move(bad_forward)), std::invalid_argument);

  std::vector<SurfaceSlice> bad_discount{
      SurfaceSlice{SviSlice{100.0, 1.0, parameters}, 0.0},
  };
  CHECK_THROW(VolatilitySurface(std::move(bad_discount)),
              std::invalid_argument);

  const VolatilitySurface valid = flat_surface(100.0, 0.05, 0.0, 0.2);
  CHECK_THROW(valid.price(100.0, 3.0, OptionType::Call), std::domain_error);
  CHECK_THROW(valid.implied_volatility(100.0, 0.1), std::domain_error);
  CHECK_THROW(valid.price(0.0, 1.0, OptionType::Call), std::invalid_argument);
  CHECK_THROW(valid.is_calendar_arbitrage_free(1e-8, 1), std::invalid_argument);

  OptionSpec american;
  american.style = ExerciseStyle::American;
  american.strike = 100.0;
  american.time_to_expiry = 1.0;
  CHECK_THROW(valid.price(american), std::invalid_argument);
}

TEST(surface_supports_a_single_slice) {
  const SviParameters parameters{0.04, 0.0, 0.0, 0.0, 0.1};
  std::vector<SurfaceSlice> slices{
      SurfaceSlice{SviSlice{100.0, 1.0, parameters}, 0.95},
  };
  const VolatilitySurface surface(std::move(slices));

  CHECK_NEAR(surface.price(100.0, 1.0, OptionType::Call),
             surface.price(100.0, 1.0, OptionType::Call), 0.0);
  CHECK_THROW(surface.price(100.0, 0.5, OptionType::Call), std::domain_error);
}

TEST(surface_handles_zero_variance) {
  const SviParameters parameters{0.0, 0.0, 0.0, 0.0, 0.1};
  std::vector<SurfaceSlice> slices{
      SurfaceSlice{SviSlice{100.0, 1.0, parameters}, 1.0},
  };
  const VolatilitySurface surface(std::move(slices));

  CHECK_NEAR(surface.price(90.0, 1.0, OptionType::Call), 10.0, 1e-12);
  CHECK_NEAR(surface.price(110.0, 1.0, OptionType::Put), 10.0, 1e-12);
  CHECK_NEAR(surface.implied_volatility(100.0, 1.0), 0.0, 1e-12);
}

TEST(surface_implied_volatility_matches_each_slice) {
  const SviParameters parameters{0.02, 0.12, -0.35, -0.1, 0.18};
  const SviSlice slice{100.0, 1.0, parameters};
  std::vector<SurfaceSlice> slices{SurfaceSlice{slice, 0.95}};
  const VolatilitySurface surface(std::move(slices));

  for (int index = -10; index <= 10; ++index) {
    const double k = 0.1 * static_cast<double>(index);
    const double strike = slice.forward * std::exp(k);
    CHECK_NEAR(surface.implied_volatility(strike, 1.0),
               slice.implied_volatility(k), 1e-12);
  }
}

TEST(surface_local_variance_of_a_flat_surface_is_the_implied_variance) {
  const double volatility = 0.25;
  const VolatilitySurface surface = flat_surface(100.0, 0.05, 0.02, volatility);

  for (const double time : {0.25, 0.5, 1.0, 2.0}) {
    for (const double strike : {70.0, 85.0, 100.0, 115.0, 130.0}) {
      CHECK_NEAR(surface.local_variance(strike, time), volatility * volatility,
                 1e-10);
      CHECK_NEAR(surface.local_volatility(strike, time), volatility, 1e-10);
    }
  }
}

TEST(surface_local_variance_follows_a_term_structure) {
  const SviParameters low{0.04 * 0.5, 0.0, 0.0, 0.0, 0.1};
  const SviParameters high{0.09 * 1.0, 0.0, 0.0, 0.0, 0.1};
  std::vector<SurfaceSlice> slices{
      SurfaceSlice{SviSlice{100.0, 0.5, low}, 1.0},
      SurfaceSlice{SviSlice{100.0, 1.0, high}, 1.0},
  };
  const VolatilitySurface surface(std::move(slices));

  const double expected = (0.09 * 1.0 - 0.04 * 0.5) / (1.0 - 0.5);
  for (const double time : {0.5, 0.75, 1.0}) {
    CHECK_NEAR(surface.local_variance(100.0, time), expected, 1e-12);
    CHECK_NEAR(surface.local_volatility(100.0, time), std::sqrt(expected),
               1e-12);
  }
}

TEST(surface_local_volatility_floors_negative_variance) {
  const SviParameters high{0.09 * 0.5, 0.0, 0.0, 0.0, 0.1};
  const SviParameters low{0.04 * 1.0, 0.0, 0.0, 0.0, 0.1};
  std::vector<SurfaceSlice> slices{
      SurfaceSlice{SviSlice{100.0, 0.5, high}, 1.0},
      SurfaceSlice{SviSlice{100.0, 1.0, low}, 1.0},
  };
  const VolatilitySurface surface(std::move(slices));

  CHECK(surface.local_variance(100.0, 0.75) < 0.0);
  CHECK_NEAR(surface.local_volatility(100.0, 0.75, 1e-4), 0.01, 1e-12);
}

TEST(surface_local_variance_validates_inputs) {
  const SviParameters parameters{0.04, 0.0, 0.0, 0.0, 0.1};
  std::vector<SurfaceSlice> single{
      SurfaceSlice{SviSlice{100.0, 1.0, parameters}, 0.95},
  };
  const VolatilitySurface one_slice(std::move(single));
  CHECK_THROW(one_slice.local_variance(100.0, 1.0), std::domain_error);

  const VolatilitySurface surface = flat_surface(100.0, 0.05, 0.0, 0.2);
  CHECK_THROW(surface.local_variance(0.0, 1.0), std::invalid_argument);
  CHECK_THROW(surface.local_variance(-10.0, 1.0), std::invalid_argument);
  CHECK_NEAR(surface.local_variance(100.0, 0.1), 0.04, 1e-10);
  CHECK_NEAR(surface.local_variance(100.0, 0.0), 0.04, 1e-10);
  CHECK_THROW(surface.local_variance(100.0, -0.1), std::invalid_argument);
  CHECK_THROW(surface.local_variance(100.0, 3.0), std::domain_error);
  CHECK_THROW(surface.local_volatility(100.0, 1.0, 0.0), std::invalid_argument);
}

TEST(surface_local_variance_matches_numerical_dupire) {
  std::vector<SurfaceSlice> slices;
  for (const double time : {0.5, 1.0, 2.0}) {
    const SviParameters parameters{0.01 * time, 0.12 * time, -0.4, 0.0, 0.15};
    slices.push_back(SurfaceSlice{SviSlice{100.0, time, parameters}, 1.0});
  }
  const VolatilitySurface surface(std::move(slices));

  const double time = 1.0;
  const double dt = 1e-4;
  for (const double strike : {80.0, 90.0, 100.0, 110.0, 120.0}) {
    const double price = surface.price(strike, time, OptionType::Call);
    const double time_up = surface.price(strike, time + dt, OptionType::Call);
    const double time_down = surface.price(strike, time - dt, OptionType::Call);
    const double dc_dt = (time_up - time_down) / (2.0 * dt);

    const double step = 0.01 * strike;
    const double strike_up =
        surface.price(strike + step, time, OptionType::Call);
    const double strike_down =
        surface.price(strike - step, time, OptionType::Call);
    const double d2c_dk2 =
        (strike_up - 2.0 * price + strike_down) / (step * step);

    const double numerical = dc_dt / (0.5 * strike * strike * d2c_dk2);
    const double analytic = surface.local_variance(strike, time);
    CHECK_NEAR(analytic, numerical, 0.02 * numerical);
  }
}
