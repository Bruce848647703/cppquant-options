# Crank-Nicolson PDE

## Equation

Let `tau = T - t` be time remaining. The Black-Scholes PDE becomes

```text
V_tau = 0.5 sigma^2 S^2 V_SS + (r - q) S V_S - r V
```

The implementation solves this on a uniform spot grid with Crank-Nicolson in
time and Dirichlet boundaries.

## Grid

The spot grid is `S_i = i * dS` for `i = 0..M`, with

```text
S_max = max(spot, strike) * exp(width_stddev * sigma * sqrt(T))
```

The grid is shifted so the current spot lies exactly on a node. This removes an
O(dS) interpolation error from the price and Greeks near the strike. The
number of space and time intervals is configurable.

## Crank-Nicolson

For each interior node:

```text
a_i = 0.5 sigma^2 S_i^2 / dS^2 - 0.5 (r - q) S_i / dS
b_i = -sigma^2 S_i^2 / dS^2 - r
c_i = 0.5 sigma^2 S_i^2 / dS^2 + 0.5 (r - q) S_i / dS
```

The scheme is

```text
(I - theta dt L) V^{n+1} = (I + (1 - theta) dt L) V^n
```

with `theta = 0.5` (Crank-Nicolson). Each step solves a tridiagonal system with
the Thomas algorithm.

### Rannacher smoothing

The payoff has a kink at the strike. Pure Crank-Nicolson can oscillate around
that kink. The implementation takes two half-sized implicit-Euler steps
(`theta = 1`) before the Crank-Nicolson steps. This is the standard Rannacher
start and preserves second-order accuracy in smooth regions.

## Boundaries

- `S = 0`: call `0`; European put `K e^(-r tau)`; American put `K`.
- `S = S_max`: put `0`; call European value
  `S_max e^(-q tau) - K e^(-r tau)`, or `max(..., S_max - K)` for an American
  call.

American exercise is applied after each time step by projecting the solution
onto the intrinsic value:

```text
V_i = max(V_i, payoff(S_i))
```

This is the usual explicit projection for American finite differences.

## Spot interpolation

The spot is placed on a grid node, so the price is read directly. Delta and
gamma are central differences at the two bracketing nodes and linearly
interpolated if necessary.

## API

```cpp
#include "cppquant/pde.hpp"

OptionSpec spec{OptionType::Put, ExerciseStyle::American,
                100.0, 100.0, 1.0, 0.2, 0.05, 0.0};

PdeSettings settings;
settings.space_steps = 400;
settings.time_steps = 400;

const PdeResult result = crank_nicolson_price(spec, settings);
// result.price, result.delta, result.gamma
```

## Accuracy

European prices are second order in `dS` and `dt` for smooth regions. On the
ATM call benchmark:

| Grid | European error vs Black-Scholes |
| ---: | ---: |
| 100 × 100 | 2.8e-2 |
| 200 × 200 | 6.9e-3 |
| 400 × 400 | 1.7e-3 |
| 800 × 800 | 4.3e-4 |
| 1,600 × 1,600 | 1.1e-4 |
| 3,000 × 3,000 | ~1.5e-5 (deep ITM test case) |

The test suite uses a 3,000-node grid when checking the `1e-4` European
accuracy target. American results are compared with a 1,000-step binomial tree
within 1%.

## Degenerate cases

- `T = 0` returns intrinsic value.
- `volatility = 0` falls back to the deterministic binomial price.
- Fewer than 10 space steps, zero time steps or a non-positive width are
  rejected.

## Limitations

- Uniform grid; no local refinement around the strike or barriers.
- Explicit American projection is first-order accurate near the exercise
  boundary, even though the European scheme is second order in smooth regions.
- No local volatility, discrete dividends or barrier conditions yet.
- `S_max` is chosen from a log-normal width; very heavy tails may need a wider
  grid.

## Benchmark

`options_bench` prices the ATM European call on a 400 × 400 grid in about
0.85 ms per price (1.2k prices/s) on the reference machine.
