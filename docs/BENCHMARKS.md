# Benchmarks

Reference environment: Apple M-series (arm64), Apple clang 21.0.0,
`-O3 -DNDEBUG`, 2026-11-04.

| Workload | Volume | Latency | Throughput |
| --- | ---: | ---: | ---: |
| Black-Scholes price + Greeks | 1,024,000 options | 70.0 ns/option | 14.28 M options/s |
| Implied volatility | 20,000 options | 375 ns/option | 2.66 M options/s |
| SVI calibration, 25 quotes | 2,000 fits | 16.5 µs/fit | 60.5 k fits/s |
| Surface price + IV, 4 slices | 1,000,000 ops | 76.1 ns/op | 13.14 M ops/s |
| LSM American put | 200,000 paths × 50 steps | 354 ms | 0.60 M paths/s |
| PDE Crank-Nicolson | 400 × 400 grid | 0.87 ms/price | 1.2 k prices/s |
| PDE Dupire local vol | 400 × 400 grid | 6.4 ms/price | 157 prices/s |
| Barrier analytic price | 200,000 prices | 49 ns/price | 20.2 M prices/s |
| Barrier price + delta + vega | 200,000 prices | 256 ns/price | 3.9 M prices/s |
| PDE barrier down-and-out | 400 × 400 grid | 0.87 ms/price | 1.2 k prices/s |
| Binomial American, 500 steps | 200 prices | 107.1 µs/price | 9.3 k prices/s |
| Monte Carlo, 200,000 paths | 20 runs | 3.25 ms/run | 61.46 M paths/s |

Commands:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/options_bench
```

Monte Carlo uses antithetic variates and the terminal-stock control variate.
Numbers are indicative; run the benchmark on your own hardware.

## W1: SVI slice calibration

Raw SVI slices are fitted in total-variance space with a constrained
Levenberg-Marquardt optimizer. The benchmark fits 25 synthetic quotes with four
deterministic starting points.

| Configuration | Fits | Latency | Throughput |
| --- | ---: | ---: | ---: |
| 25 quotes, 5 parameters | 2,000 | 17.4 µs/fit | 57.6 k fits/s |

Accuracy on the benchmark slice is at machine precision because the quotes are
noiseless. The tests cover symmetric and asymmetric parameter recovery, narrow
strike ranges and invalid input.

## W2: volatility surface

The surface benchmark builds four flat SVI slices (0.25y to 2y), then performs
1,000,000 price + implied-volatility lookups at interpolated expiries.

| Configuration | Operations | Latency | Throughput |
| --- | ---: | ---: | ---: |
| 4 slices, price + IV | 1,000,000 | 81.3 ns/op | 12.3 M ops/s |

The arbitrage checks are separate from the hot lookup path and are intended
for validation and calibration workflows, not per-tick pricing.

## W3: Longstaff-Schwartz

The LSM benchmark prices an ATM American put (S=100, K=100, T=1, vol=20%,
r=5%) with 200,000 antithetic paths and 50 exercise steps.

| Configuration | Paths × steps | Time | Throughput | Price | SE |
| --- | ---: | ---: | ---: | ---: | ---: |
| LSM American put | 200,000 × 50 | 331 ms | 0.60 M paths/s | 6.0523 | 0.0083 |

The 1,000-step binomial reference is 6.0896. The difference is within 1% and
within four standard errors.

## W4: LSM variance reduction

The European payoff is used as a control variate with the Black-Scholes price
as its known expectation. Pair averages are formed before estimating beta when
antithetic variates are enabled.

At 200,000 paths and 50 exercise steps on the ATM American put:

| Estimator | Price | Standard error |
| --- | ---: | ---: |
| Plain, no antithetic, no control | 6.0499 | 0.0160 |
| Antithetic only | 6.0531 | 0.0098 |
| Antithetic + control | 6.0523 | 0.0083 |

The combined estimator cuts the standard error by about 48% versus plain Monte
Carlo and by about 15% on top of antithetic variates alone.

`options_bench --convergence` produces the full path-count table:

| Paths | Antithetic price / SE | Control price / SE |
| ---: | ---: | ---: |
| 10,000 | 6.0167 / 0.0436 | 6.0529 / 0.0375 |
| 50,000 | 6.0398 / 0.0196 | 6.0432 / 0.0167 |
| 200,000 | 6.0531 / 0.0098 | 6.0523 / 0.0083 |
| 1,000,000 | 6.0558 / 0.0044 | 6.0563 / 0.0037 |

The control variate reduces noise, not the discrete-exercise or in-sample
regression bias.

## W5: Crank-Nicolson PDE

The PDE benchmark prices an ATM European call with Rannacher smoothing on a
400 × 400 grid.

| Grid | Operations | Latency | Throughput |
| ---: | ---: | ---: | ---: |
| 400 × 400 | 500 prices | 0.87 ms/price | 1.2 k prices/s |

European error versus Black-Scholes:

| Grid | Error |
| ---: | ---: |
| 100 × 100 | 2.8e-2 |
| 200 × 200 | 6.9e-3 |
| 400 × 400 | 1.7e-3 |
| 800 × 800 | 4.3e-4 |
| 1,600 × 1,600 | 1.1e-4 |

The error scales as O(dS^2). The tests use a 3,000-node grid for the `1e-4`
accuracy check. American prices are within 1% of a 1,000-step binomial tree.

## W6: Dupire local volatility

The benchmark builds a flat 20% surface with four expiries and prices the ATM
European call through the local-vol PDE. Every interior node and time step
queries the surface, so the cost is dominated by Dupire variance evaluation
rather than the tridiagonal solve.

| Grid | Operations | Latency | Throughput | Floored points |
| ---: | ---: | ---: | ---: | ---: |
| 400 × 400 | 101 prices | 6.36 ms/price | 157 prices/s | 0 |

The constant-volatility PDE runs at about 0.87 ms on the same grid, so local
variance adds roughly 7x. The flat surface reproduces Black-Scholes within `2e-3` on an
800 × 800 grid; a two-slice term structure integrates to the same total
variance.

## W7: barrier options

The analytic benchmark prices the same European down-and-out call 200,000 times.
The Greeks benchmark adds central-difference delta and vega. The PDE benchmark
places the barrier on the lower grid boundary.

| Workload | Operations | Latency | Throughput |
| --- | ---: | ---: | ---: |
| Analytic price | 200,000 | 49.4 ns/op | 20.2 M ops/s |
| Price + delta + vega | 200,000 | 256 ns/op | 3.9 M ops/s |
| PDE, 400 × 400 | 200 | 0.87 ms/price | 1.15 k prices/s |

The 800 × 800 PDE agrees with the analytic price within `7e-4` across barriers
from 60 to 99. A barrier near zero recovers the vanilla Black-Scholes price.

## W8: release

The W1-W7 numbers above are the final v0.5.0 baselines. The release adds no
hot-path changes; the barrier and local-vol entries are the reference for the
later pricing paths.
