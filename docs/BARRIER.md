# Down-and-out barrier options

## Contract

A European down-and-out option pays the vanilla payoff at expiry unless the
spot touches the barrier at any time before expiry. The implementation uses a
continuous barrier and no rebate.

```cpp
OptionSpec spec{OptionType::Call, ExerciseStyle::European,
                100.0, 100.0, 1.0, 0.2, 0.05, 0.02};

const BarrierGreeks greeks = down_and_out_greeks(spec, 80.0);
// greeks.price, greeks.delta, greeks.vega
```

The spot must be above the barrier. If it is at or below the barrier, the
option is already knocked out and the price, delta and vega are zero. American
exercise is rejected because the analytic formulas assume European exercise.

## Notation

```text
lambda = (r - q + sigma^2 / 2) / sigma^2
S'     = H^2 / S
c      = (H / S)^(2 lambda - 2)
```

`r` is the risk-free rate, `q` the dividend yield, `sigma` the volatility, `T`
time to expiry, `H` the barrier and `S` the spot.

## Down-and-out call

When `H <= K`:

```text
C = C_BS(S, K) - c C_BS(S', K)
```

When `H > K`, the barrier lies above the strike and the reflected payoff needs
the full expression:

```text
x1 = log(S / H) / (sigma sqrt(T)) + lambda sigma sqrt(T)
y1 = log(H / S) / (sigma sqrt(T)) + lambda sigma sqrt(T)

C = S e^(-qT) N(x1) - K e^(-rT) N(x1 - sigma sqrt(T))
    - S e^(-qT) (H / S)^(2 lambda) N(y1)
    + K e^(-rT) (H / S)^(2 lambda - 2) N(y1 - sigma sqrt(T))
```

## Down-and-out put

When `H >= K`, the put is knocked out before it can finish in the money, so its
value is zero. When `H < K`, the price follows from the reflection density with
truncated normal terms:

```text
mu    = r - q - sigma^2 / 2
s     = sigma sqrt(T)
b     = log(H / S)
d_b   = ( b - mu T) / s
d_-b  = (-b - mu T) / s
T(S,d)= K e^(-rT) N(d) - S e^(-qT) N(d - s)

I1 = P_BS(S, K) - T(S, d_b)
I2 = P_BS(S', K) - T(S', d_-b)

P = I1 - c I2
```

The implementation was checked against direct numerical integration of the
reflection density, not only against the same closed form.

## Greeks

Delta and vega are central finite differences of the analytic price. Vega is
per unit of volatility. Barrier vega can be negative for options close to the
barrier, because an increase in volatility raises the knock-out probability.

## PDE cross-check

`crank_nicolson_down_and_out_price` solves the same contract on a spot grid
whose lower boundary is the barrier itself. The spot is aligned to a grid node,
the barrier node is held at zero and American projection is not used.

```cpp
PdeSettings settings;
settings.space_steps = 800;
settings.time_steps = 800;

const PdeResult pde =
    crank_nicolson_down_and_out_price(spec, 80.0, settings);
```

On an 800 x 800 grid the PDE agrees with the analytic price within `7e-4` for
the tested barriers:

| Barrier | Call analytic | Call PDE | Put analytic | Put PDE |
| ---: | ---: | ---: | ---: | ---: |
| 60 | 9.227005 | 9.226729 | 5.976977 | 5.976319 |
| 80 | 9.133306 | 9.133036 | 1.732678 | 1.732335 |
| 95 | 4.883524 | 4.883467 | 0.010277 | 0.010243 |
| 99 | 1.170361 | 1.170354 | 0.000016 | 0.000014 |

A barrier close to zero converges to the vanilla Black-Scholes price; a barrier
at or above the spot returns zero.

## Benchmark

On the 400 x 400 grid:

| Workload | Operations | Latency | Throughput |
| --- | ---: | ---: | ---: |
| Analytic down-and-out price | 200,000 | 49.4 ns/op | 20.2 M ops/s |
| Price + delta + vega | 200,000 | 256 ns/op | 3.9 M ops/s |
| PDE down-and-out | 200 | 0.87 ms/price | 1.15 k prices/s |

## Limitations

- Continuous monitoring only; no discrete barrier adjustment.
- No rebate.
- No American barrier exercise.
- The PDE grid is shifted to the barrier but not locally refined near it.
- Only down-and-out contracts are implemented; up-and-out and double barriers
  are not.