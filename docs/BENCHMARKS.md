# Benchmarks

Reference environment: Apple M-series (arm64), Apple clang 21.0.0,
`-O3 -DNDEBUG`, 2026-11-04.

| Workload | Volume | Latency | Throughput |
| --- | ---: | ---: | ---: |
| Black-Scholes price + Greeks | 1,024,000 options | 80.3 ns/option | 12.46 M options/s |
| Implied volatility | 20,000 options | 379 ns/option | 2.64 M options/s |
| SVI calibration, 25 quotes | 2,000 fits | 17.4 µs/fit | 57.6 k fits/s |
| Surface price + IV, 4 slices | 1,000,000 ops | 81.3 ns/op | 12.3 M ops/s |
| LSM American put | 200,000 paths × 50 steps | 311 ms | 0.64 M paths/s |
| Binomial American, 500 steps | 200 prices | 105.7 µs/price | 9.5 k prices/s |
| Monte Carlo, 200,000 paths | 20 runs | 3.03 ms/run | 66.0 M paths/s |

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
| LSM American put | 200,000 × 50 | 311 ms | 0.64 M paths/s | 6.0531 | 0.0161 |

The 1,000-step binomial reference is 6.0896. The difference is within 1% and
within four standard errors. The next milestone adds a European control
variate to reduce the Monte Carlo standard error.

## Planned baselines

- W4: LSM price/standard error versus path count, basis choice and control
  variate.
- W5: PDE grid size versus error and runtime.
- W6: Dupire local-vol build time.
