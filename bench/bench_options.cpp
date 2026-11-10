#include "cppquant/american_mc.hpp"
#include "cppquant/barrier.hpp"
#include "cppquant/options.hpp"
#include "cppquant/pde.hpp"
#include "cppquant/surface.hpp"
#include "cppquant/volatility.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace cppquant;

namespace {

double sink = 0.0;

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

std::vector<OptionSpec> make_specs(std::size_t count, std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::vector<OptionSpec> specs;
  specs.reserve(count);

  for (std::size_t index = 0; index < count; ++index) {
    OptionSpec spec;
    spec.type = (index & 1U) != 0 ? OptionType::Put : OptionType::Call;
    spec.style = ExerciseStyle::European;
    spec.spot = 50.0 + static_cast<double>(rng() % 151);
    spec.strike = 50.0 + static_cast<double>(rng() % 151);
    spec.time_to_expiry = 0.05 + static_cast<double>(rng() % 500) / 100.0;
    spec.volatility = 0.05 + static_cast<double>(rng() % 80) / 100.0;
    spec.risk_free_rate = static_cast<double>(rng() % 8) / 100.0;
    spec.dividend_yield = static_cast<double>(rng() % 5) / 100.0;
    specs.push_back(spec);
  }

  return specs;
}

void run_black_scholes(const std::vector<OptionSpec>& specs,
                       std::size_t repetitions) {
  const auto start = Clock::now();
  std::size_t operations = 0;
  for (std::size_t repeat = 0; repeat < repetitions; ++repeat) {
    for (const auto& spec : specs) {
      const double price = black_scholes_price(spec);
      const Greeks greeks = black_scholes_greeks(spec);
      sink += price + greeks.delta + greeks.gamma + greeks.vega + greeks.theta +
              greeks.rho;
      ++operations;
    }
  }
  const double seconds = seconds_since(start);
  std::printf(
      "black-scholes price+greeks  %10zu opts  %8.1f ns/opt  %6.2f M opts/s\n",
      operations, seconds * 1e9 / static_cast<double>(operations),
      static_cast<double>(operations) / seconds / 1e6);
}

void run_implied_volatility(const std::vector<OptionSpec>& specs,
                            std::size_t operations) {
  std::vector<double> prices(specs.size());
  for (std::size_t index = 0; index < specs.size(); ++index) {
    prices[index] = black_scholes_price(specs[index]);
  }

  const auto start = Clock::now();
  for (std::size_t index = 0; index < operations; ++index) {
    const std::size_t slot = index % specs.size();
    sink += implied_volatility(specs[slot], prices[slot]);
  }
  const double seconds = seconds_since(start);
  std::printf(
      "implied volatility          %10zu opts  %8.1f ns/opt  %6.2f M opts/s\n",
      operations, seconds * 1e9 / static_cast<double>(operations),
      static_cast<double>(operations) / seconds / 1e6);
}

void run_svi(std::size_t operations) {
  const SviParameters truth{0.01, 0.12, -0.4, -0.1, 0.15};
  const SviSlice slice{100.0, 1.0, truth};
  std::vector<VolQuote> quotes;
  quotes.reserve(25);
  for (int index = 0; index < 25; ++index) {
    const double k = -1.0 + 2.0 * static_cast<double>(index) / 24.0;
    const double strike = slice.forward * std::exp(k);
    quotes.push_back(VolQuote{strike, slice.implied_volatility(k), 1.0});
  }

  const auto start = Clock::now();
  for (std::size_t index = 0; index < operations; ++index) {
    const auto result = calibrate_svi(100.0, 1.0, quotes);
    sink += result.rmse + result.slice.parameters.b +
            static_cast<double>(result.iterations);
  }
  const double seconds = seconds_since(start);
  std::printf(
      "SVI calibration (25 quotes) %10zu fits  %8.1f us/fit  %6.1f k fits/s\n",
      operations, seconds * 1e6 / static_cast<double>(operations),
      static_cast<double>(operations) / seconds / 1e3);
}

void run_surface(std::size_t operations) {
  std::vector<SurfaceSlice> slices;
  for (const double time : {0.25, 0.5, 1.0, 2.0}) {
    const double forward = 100.0 * std::exp(0.03 * time);
    const SviParameters parameters{0.04 * time, 0.0, 0.0, 0.0, 0.1};
    slices.push_back(SurfaceSlice{SviSlice{forward, time, parameters},
                                  std::exp(-0.03 * time)});
  }
  const VolatilitySurface surface(std::move(slices));

  const auto start = Clock::now();
  for (std::size_t index = 0; index < operations; ++index) {
    const double moneyness =
        -0.6 + 1.2 * static_cast<double>(index % 121) / 120.0;
    const double time = 0.25 + 1.75 * static_cast<double>(index % 4) / 3.0;
    const double strike = surface.forward(time) * std::exp(moneyness);
    sink += surface.price(strike, time, OptionType::Call);
    sink += surface.implied_volatility(strike, time);
  }
  const double seconds = seconds_since(start);
  std::printf(
      "surface price+IV (4 slices) %9zu ops  %8.1f ns/op  %6.2f M ops/s\n",
      operations, seconds * 1e9 / static_cast<double>(operations),
      static_cast<double>(operations) / seconds / 1e6);
}

void run_lsm(std::size_t paths, std::size_t steps) {
  OptionSpec spec{OptionType::Put,
                  ExerciseStyle::American,
                  100.0,
                  100.0,
                  1.0,
                  0.2,
                  0.05,
                  0.0};
  LongstaffSchwartzSettings settings;
  settings.paths = paths;
  settings.exercise_steps = steps;
  settings.seed = 42;

  const auto start = Clock::now();
  const auto result = longstaff_schwartz_price(spec, settings);
  const double seconds = seconds_since(start);
  sink += result.price + result.standard_error;

  std::printf("LSM American put (%zu paths, %zu steps)  price %.4f  se %.4f  "
              "%.1f ms  %.1f M paths/s\n",
              result.effective_paths, steps, result.price,
              result.standard_error, seconds * 1e3,
              static_cast<double>(result.effective_paths) / seconds / 1e6);
}

void run_lsm_convergence() {
  OptionSpec spec{OptionType::Put,
                  ExerciseStyle::American,
                  100.0,
                  100.0,
                  1.0,
                  0.2,
                  0.05,
                  0.0};
  std::printf("\nLSM convergence (ATM American put, 50 steps)\n");
  std::printf("%10s %8s %12s %10s %12s %10s\n", "paths", "control", "price",
              "se", "ms", "beta");
  for (std::size_t paths : {10'000, 50'000, 200'000, 1'000'000}) {
    for (const bool control : {false, true}) {
      LongstaffSchwartzSettings settings;
      settings.paths = paths;
      settings.exercise_steps = 50;
      settings.seed = 42;
      settings.antithetic = true;
      settings.control_variate = control;

      const auto start = Clock::now();
      const auto result = longstaff_schwartz_price(spec, settings);
      const double seconds = seconds_since(start);
      sink += result.price + result.standard_error;
      std::printf("%10zu %8s %12.6f %10.6f %12.1f %10.4f\n",
                  result.effective_paths, control ? "yes" : "no", result.price,
                  result.standard_error, seconds * 1e3,
                  result.control_variate_beta);
    }
  }
}

void run_pde(std::size_t space_steps, std::size_t time_steps,
             std::size_t operations) {
  OptionSpec spec{OptionType::Call,
                  ExerciseStyle::European,
                  100.0,
                  100.0,
                  1.0,
                  0.2,
                  0.05,
                  0.0};
  PdeSettings settings;
  settings.space_steps = space_steps;
  settings.time_steps = time_steps;

  const auto start = Clock::now();
  for (std::size_t index = 0; index < operations; ++index) {
    const auto result = crank_nicolson_price(spec, settings);
    sink += result.price + result.delta + result.gamma;
  }
  const double seconds = seconds_since(start);
  std::printf("PDE Crank-Nicolson (%zux%zu) %7zu prices  %8.2f us/price  %6.1f "
              "k prices/s\n",
              space_steps, time_steps, operations,
              seconds * 1e6 / static_cast<double>(operations),
              static_cast<double>(operations) / seconds / 1e3);
}

void run_local_vol_pde(std::size_t space_steps, std::size_t time_steps,
                       std::size_t operations) {
  std::vector<SurfaceSlice> slices;
  for (const double time : {0.25, 0.5, 1.0, 2.0}) {
    const double forward = 100.0 * std::exp(0.05 * time);
    const SviParameters parameters{0.04 * time, 0.0, 0.0, 0.0, 0.1};
    slices.push_back(SurfaceSlice{SviSlice{forward, time, parameters},
                                  std::exp(-0.05 * time)});
  }
  const VolatilitySurface surface(std::move(slices));

  OptionSpec spec{OptionType::Call,
                  ExerciseStyle::European,
                  100.0,
                  100.0,
                  1.0,
                  0.2,
                  0.05,
                  0.0};
  PdeSettings settings;
  settings.space_steps = space_steps;
  settings.time_steps = time_steps;

  std::size_t clamped = 0;
  const auto start = Clock::now();
  for (std::size_t index = 0; index < operations; ++index) {
    const auto result = crank_nicolson_local_vol_price(spec, surface, settings);
    sink += result.price + result.delta + result.gamma;
    clamped += result.clamped_local_variance_points;
  }
  const double seconds = seconds_since(start);
  std::printf(
      "PDE Dupire local vol (%zux%zu) %6zu prices  %8.2f us/price  %6.1f "
      "prices/s  clamped=%zu\n",
      space_steps, time_steps, operations,
      seconds * 1e6 / static_cast<double>(operations),
      static_cast<double>(operations) / seconds, clamped);
}

void run_barrier(std::size_t operations) {
  const OptionSpec spec{OptionType::Call,
                        ExerciseStyle::European,
                        100.0,
                        100.0,
                        1.0,
                        0.2,
                        0.05,
                        0.02};
  const auto start = Clock::now();
  for (std::size_t index = 0; index < operations; ++index) {
    sink += down_and_out_price(spec, 80.0);
  }
  const double seconds = seconds_since(start);
  std::printf(
      "barrier down-and-out price %9zu ops  %8.1f ns/op  %6.2f M ops/s\n",
      operations, seconds * 1e9 / static_cast<double>(operations),
      static_cast<double>(operations) / seconds / 1e6);

  const auto greeks_start = Clock::now();
  for (std::size_t index = 0; index < operations; ++index) {
    const auto greeks = down_and_out_greeks(spec, 80.0);
    sink += greeks.price + greeks.delta + greeks.vega;
  }
  const double greeks_seconds = seconds_since(greeks_start);
  std::printf(
      "barrier price+delta+vega  %9zu ops  %8.1f ns/op  %6.2f M ops/s\n",
      operations, greeks_seconds * 1e9 / static_cast<double>(operations),
      static_cast<double>(operations) / greeks_seconds / 1e6);
}

void run_barrier_pde(std::size_t space_steps, std::size_t time_steps,
                     std::size_t operations) {
  const OptionSpec spec{OptionType::Call,
                        ExerciseStyle::European,
                        100.0,
                        100.0,
                        1.0,
                        0.2,
                        0.05,
                        0.02};
  PdeSettings settings;
  settings.space_steps = space_steps;
  settings.time_steps = time_steps;

  const auto start = Clock::now();
  for (std::size_t index = 0; index < operations; ++index) {
    const auto result = crank_nicolson_down_and_out_price(spec, 80.0, settings);
    sink += result.price + result.delta + result.gamma;
  }
  const double seconds = seconds_since(start);
  std::printf("PDE barrier down-and-out (%zux%zu) %5zu prices  %8.2f us/price  "
              "%6.1f prices/s\n",
              space_steps, time_steps, operations,
              seconds * 1e6 / static_cast<double>(operations),
              static_cast<double>(operations) / seconds);
}

void run_binomial(const OptionSpec& spec, int steps, std::size_t operations) {
  const auto start = Clock::now();
  for (std::size_t index = 0; index < operations; ++index) {
    sink += binomial_price(spec, steps);
  }
  const double seconds = seconds_since(start);
  std::printf("binomial american (%d steps) %8zu opts  %10.1f us/opt\n", steps,
              operations, seconds * 1e6 / static_cast<double>(operations));
}

void run_monte_carlo(const OptionSpec& spec, std::size_t paths,
                     std::size_t runs, std::uint64_t seed) {
  const auto start = Clock::now();
  std::size_t total_paths = 0;
  for (std::size_t index = 0; index < runs; ++index) {
    const auto result = monte_carlo_price(spec, paths, seed + index);
    sink += result.price + result.standard_error;
    total_paths += result.effective_samples;
  }
  const double seconds = seconds_since(start);
  std::printf(
      "monte carlo (%zu paths)     %10zu runs  %8.2f ms/run  %6.2f M paths/s\n",
      paths, runs, seconds * 1e3 / static_cast<double>(runs),
      static_cast<double>(total_paths) / seconds / 1e6);
}

} // namespace

int main(int argc, char** argv) {
  std::size_t bs_repetitions = 1000;
  std::size_t iv_operations = 20'000;
  std::size_t svi_operations = 2'000;
  std::size_t surface_operations = 1'000'000;
  std::size_t lsm_paths = 200'000;
  std::size_t lsm_steps = 50;
  std::size_t pde_space = 400;
  std::size_t pde_time = 400;
  std::size_t pde_operations = 500;
  bool convergence = false;
  int binomial_steps = 500;
  std::size_t binomial_operations = 200;
  std::size_t mc_paths = 200'000;
  std::size_t mc_runs = 20;
  std::size_t barrier_operations = 200'000;
  std::size_t barrier_pde_operations = 200;

  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--quick") {
      bs_repetitions = 50;
      iv_operations = 1'000;
      svi_operations = 200;
      surface_operations = 100'000;
      lsm_paths = 50'000;
      lsm_steps = 25;
      pde_space = 100;
      pde_time = 100;
      pde_operations = 50;
      binomial_operations = 20;
      mc_paths = 50'000;
      mc_runs = 5;
      barrier_operations = 20'000;
      barrier_pde_operations = 20;
    } else if (argument == "--convergence") {
      convergence = true;
    } else if (argument == "--help" || argument == "-h") {
      std::printf("usage: %s [--quick]\n", argv[0]);
      return 0;
    } else {
      std::fprintf(stderr, "unknown argument: %s\n", argument.c_str());
      return 2;
    }
  }

  std::printf("cppquant-options benchmark\n");
  std::printf("fixed seeds, release build recommended\n\n");

  const auto specs = make_specs(1024, 7);
  run_black_scholes(specs, bs_repetitions);
  run_implied_volatility(specs, iv_operations);
  run_svi(svi_operations);
  run_surface(surface_operations);
  run_lsm(lsm_paths, lsm_steps);
  run_pde(pde_space, pde_time, pde_operations);
  run_local_vol_pde(pde_space, pde_time, pde_operations / 5 + 1);
  run_barrier(barrier_operations);
  run_barrier_pde(pde_space, pde_time, barrier_pde_operations);
  if (convergence) {
    run_lsm_convergence();
  }

  OptionSpec american;
  american.type = OptionType::Put;
  american.style = ExerciseStyle::American;
  american.spot = 100.0;
  american.strike = 100.0;
  american.time_to_expiry = 1.0;
  american.volatility = 0.2;
  american.risk_free_rate = 0.05;
  run_binomial(american, binomial_steps, binomial_operations);

  OptionSpec european = american;
  european.style = ExerciseStyle::European;
  run_monte_carlo(european, mc_paths, mc_runs, 42);

  std::printf("\nchecksum: %.6f\n", sink);
  return 0;
}
