# cppquant-options

European and American option pricing in modern C++20, with analytic
Black-Scholes-Merton Greeks, a robust implied-volatility solver, a
Cox-Ross-Rubinstein binomial tree and a Monte Carlo engine with variance
reduction. No third-party dependencies.

[![CI](https://github.com/Bruce848647703/cppquant-options/actions/workflows/ci.yml/badge.svg)](https://github.com/Bruce848647703/cppquant-options/actions/workflows/ci.yml)

## Features

- **Black-Scholes-Merton** prices and analytic Greeks for European calls and
  puts, with continuous dividend yield and support for negative rates.
- **Implied volatility** via a Newton iteration guarded by a bisection
  bracket, with explicit no-arbitrage bounds.
- **Binomial CRR tree** for European and American options, including early
  exercise with dividends.
- **Monte Carlo** for European options with antithetic variates, a
  terminal-stock control variate, standard errors and deterministic seeds.
- **Edge cases**: expired options, zero volatility, deep in/out-of-the-money
  strikes, invalid inputs.
- Tests with reference values, put-call parity, finite-difference Greeks and
  Monte Carlo convergence; CI with GCC, Clang, Apple clang and ASan/UBSan.

## Quick start

Requirements: a C++20 compiler and CMake 3.20+. Ninja is optional.

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Sanitizer build:

```bash
cmake -S . -B build-asan -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug -DCPPQUANT_SANITIZE=ON
cmake --build build-asan
ctest --test-dir build-asan --output-on-failure
```

Install and consume:

```bash
cmake --install build --prefix /tmp/cppquant
```

```cmake
find_package(cppquant_options CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE cppquant::options)
```

## Library usage

```cpp
#include "cppquant/options.hpp"

using namespace cppquant;

OptionSpec call;
call.type = OptionType::Call;
call.spot = 100.0;
call.strike = 100.0;
call.time_to_expiry = 1.0;
call.volatility = 0.20;
call.risk_free_rate = 0.05;
call.dividend_yield = 0.0;

const double price = black_scholes_price(call);   // 10.450583572...
const Greeks greeks = black_scholes_greeks(call); // delta, gamma, vega...
const double iv = implied_volatility(call, price);

OptionSpec american = call;
american.type = OptionType::Put;
american.style = ExerciseStyle::American;
const double american_put = binomial_price(american, 1000);

const auto mc = monte_carlo_price(call, 200'000, /*seed=*/42);
// mc.price, mc.standard_error, mc.effective_samples
```

### Units

- `time_to_expiry` is in years, `volatility` and rates are annualized decimals
  (`0.20` = 20%).
- `vega` is per 1.00 change in volatility; divide by 100 for "per 1%".
- `rho` is per 1.00 change in the risk-free rate; divide by 100 for "per 1%".
- `theta` is per year of calendar time.
- Prices and Greeks use the same units as `spot` and `strike`.

## Command-line pricer

```bash
./build/option_pricer \
  --type call --spot 100 --strike 100 --time 1 \
  --vol 0.2 --rate 0.05 --market 10.450583572186
```

```text
cppquant-options
call european  spot=100.0000  strike=100.0000  T=1.0000  vol=0.2000  r=0.0500  q=0.0000

black-scholes (european)
  price       10.450584
  delta        0.636831
  gamma        0.018762
  vega        37.524035
  theta       -6.414028
  rho         53.232482

binomial (european, 500 steps)
  price       10.446585

implied volatility (market 10.450584)
  vol          0.200000

monte carlo (200000 paths, antithetic + control variate)
  price       10.451654  (se 0.006111, 200000 samples)
```

Use `--style american` for the binomial American price. `--no-mc` skips Monte
Carlo and `--help` lists every flag.

## Benchmark

```bash
./build/options_bench          # full run
./build/options_bench --quick  # CI smoke test
```

Reference numbers from an Apple M-series laptop (Apple clang 21, `-O3`):

| Workload | Volume | Latency | Throughput |
| --- | ---: | ---: | ---: |
| Black-Scholes price + Greeks | 1,024,000 options | 80 ns/option | 12.5 M options/s |
| Implied volatility | 20,000 options | 379 ns/option | 2.6 M options/s |
| Binomial American, 500 steps | 200 prices | 106 µs/price | 9.5 k prices/s |
| Monte Carlo, 200,000 paths | 20 runs | 3.0 ms/run | 66 M paths/s |

These are indicative; run the benchmark on your own hardware.

## Accuracy notes

- Black-Scholes prices agree with reference values to better than `1e-9`.
- The implied-volatility solver converges to `1e-7` in volatility across a
  grid of strikes, maturities and volatilities, including negative rates.
- The CRR tree converges to Black-Scholes as steps increase; for 1,000 steps
  the European error is below `5e-3` on the tested cases.
- Monte Carlo uses the sample standard error of the estimator after variance
  reduction. Antithetic pairs are averaged before estimating the control
  variate, so both techniques reduce the variance of the reported estimator.
- Deep in-the-money prices sit at the no-arbitrage bound to machine precision.
  Implied volatility is not numerically identifiable there and the solver
  returns zero; the tests skip those degenerate round trips.

## Design

```text
include/cppquant/options.hpp   public API and data types
src/black_scholes.cpp          validation, normal CDF/PDF, analytic prices
src/implied_volatility.cpp     bracketed Newton solver
src/binomial.cpp               CRR tree, deterministic zero-vol scan
src/monte_carlo.cpp            GBM sampling, antithetic + control variate
examples/option_pricer.cpp     CLI
```

See [docs/MATH.md](docs/MATH.md) for the formulas, the Greeks, the CRR
recursion and the variance-reduction details.

### Modelling limits

- Black-Scholes assumes constant volatility and log-normal returns.
- Monte Carlo supports European exercise only; American options use the
  binomial tree.
- The CRR tree rejects parameters where the risk-neutral probability falls
  outside `[0, 1]`, rather than silently returning a misleading price.
- There is no local/stochastic volatility, jump diffusion, discrete dividend
  schedule or volatility-surface calibration.

## Tests

The repository ships a tiny header-only test harness, so the build has no
external test dependency. The 28 tests cover:

- three sets of reference prices and Greeks (ATM, dividends, negative rates);
- put-call parity and finite-difference checks of every Greek;
- expiry, zero volatility and invalid-input behaviour;
- implied-volatility round trips across 50+ combinations and no-arbitrage
  rejection;
- binomial convergence, American early exercise and zero-vol exercise scans;
- Monte Carlo accuracy, reproducibility, variance reduction and argument
  validation.

```bash
./build/options_tests             # all tests
./build/options_tests implied     # filter by substring
```

## License

MIT. See [LICENSE](LICENSE).
