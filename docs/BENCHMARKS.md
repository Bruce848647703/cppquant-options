# Benchmarks

Reference environment: Apple M-series (arm64), Apple clang 21.0.0,
`-O3 -DNDEBUG`, 2026-11-04.

| Workload | Volume | Latency | Throughput |
| --- | ---: | ---: | ---: |
| Black-Scholes price + Greeks | 1,024,000 options | 80.3 ns/option | 12.46 M options/s |
| Implied volatility | 20,000 options | 379 ns/option | 2.64 M options/s |
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

## Planned baselines

- W1–W2: SVI calibration time and surface IV lookup latency.
- W3–W4: LSM price/standard error versus path count and basis choice.
- W5: PDE grid size versus error and runtime.
- W6: Dupire local-vol build time.