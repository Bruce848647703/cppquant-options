# SVI calibration

## Parameterization

For one expiry with forward `F` and time to expiry `T`, define log-moneyness

```text
k = log(K / F)
```

and total variance `w(k) = implied_vol(k)^2 T`. Raw SVI models the slice as

```text
w(k) = a + b (rho (k - m) + sqrt((k - m)^2 + sigma^2))
```

with `b >= 0`, `|rho| < 1` and `sigma > 0`. The minimum total variance is

```text
w_min = a + b sigma sqrt(1 - rho^2)
```

so a slice is free of the most basic static arbitrage when `w_min >= 0`.
Butterfly and calendar conditions across slices are handled at the surface
level in the next milestone.

## Why fit total variance

Quotes arrive as implied volatilities. Converting them to total variance makes
the objective a quadratic in the SVI parameters and avoids dividing by
`T` during optimization. The reported errors are converted back to volatility
units so they are directly comparable to quote conventions.

## Unconstrained optimization

The constraints are enforced by reparameterization:

```text
a_min = exp(u_a)                  # minimum total variance
b     = exp(u_b)
rho   = tanh(u_rho)
sigma = exp(u_sigma)
a     = a_min - b sigma sqrt(1 - rho^2)
```

`a_min` is guaranteed positive, `b` and `sigma` are positive and `rho` stays
inside `(-1, 1)`. The free parameters are `u_a`, `u_b`, `u_rho`, `m` and
`u_sigma`.

The optimizer is Levenberg-Marquardt on the weighted residuals

```text
r_i = sqrt(weight_i) * (w_model(k_i) - w_quote_i)
```

with a numerical central-difference Jacobian and a 5x5 Gaussian elimination
with partial pivoting. A step is accepted only if it lowers the cost; otherwise
the damping factor is increased and the step is retried. The fit starts from
several deterministic initial points (different `rho` and `sigma`) and keeps the
lowest-cost result, which avoids most local minima on typical smiles.

## API

```cpp
#include "cppquant/volatility.hpp"

using namespace cppquant;

std::vector<VolQuote> quotes = {
    {90.0, 0.28, 1.0},
    {95.0, 0.24, 1.0},
    {100.0, 0.22, 1.0},
    {105.0, 0.23, 1.0},
    {110.0, 0.26, 1.0},
};

const SviFitResult fit = calibrate_svi(/*forward=*/100.0,
                                       /*time=*/0.5, quotes);

const SviSlice& slice = fit.slice;
const double atm_vol = slice.implied_volatility_for_strike(100.0);
const double wing_var = slice.total_variance(std::log(120.0 / 100.0));

// fit.rmse, fit.max_error, fit.iterations, fit.converged
```

`SviFitResult::converged` reports whether the optimizer met its step or cost
tolerance. A non-converged fit is still returned so the caller can inspect it;
real quote data often has noise and may not drive the cost to zero.

## Validation rules

`calibrate_svi` rejects:

- non-positive or non-finite forward, expiry, strike, volatility or weight;
- fewer than five quotes or fewer than five distinct strikes;
- non-positive iteration counts, tolerance or initial damping.

The fit itself cannot produce a negative `w_min` because `a_min` is
exponentiated. `SviSlice::satisfies_basic_constraints()` checks the result
independently, which is useful when slices are constructed or loaded from
external data.

## Limitations

- One expiry per slice; the surface and calendar-arbitrage checks are W2.
- The objective is variance-space least squares. Vega-weighted or robust
  losses can be added later if real quote data requires them.
- No bid/ask spread constraints, no penalization of extreme parameters and no
  cross-expiry regularization yet.

## Tests

- The raw formula and implied-vol conversion are checked against a direct
  calculation.
- A 25-strike synthetic smile recovers the generating parameters to `1e-3`.
- An asymmetric smile with `rho > 0` and `m < 0` recovers the sign and shape.
- Narrow strike ranges, weighted quotes and repeated runs are deterministic.
- Invalid inputs and invalid constraints are rejected.
- Fitted variance stays non-negative well outside the quoted strike range.
