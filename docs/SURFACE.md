# Volatility surface

## Construction

`VolatilitySurface` takes one or more `SurfaceSlice` values:

```cpp
SurfaceSlice {
  SviSlice slice;          // forward, expiry, raw SVI parameters
  double discount_factor;  // df(T) for pricing
};
```

Slices are sorted by expiry. Duplicate expiries, non-positive forwards,
non-positive discount factors and slices that violate the basic SVI constraints
are rejected.

## Interpolation

For an expiry `T` between two slices `T_i <= T <= T_{i+1}` with fraction

```text
f = (T - T_i) / (T_{i+1} - T_i)
```

the surface uses:

```text
log F(T) = (1 - f) log F_i + f log F_{i+1}
log df(T) = (1 - f) log df_i + f log df_{i+1}
w(k, T) = (1 - f) w_i(k) + f w_{i+1}(k)
```

Total variance is interpolated linearly in time, which preserves calendar
monotonicity when the two endpoint slices are themselves monotone. The forward
and discount factor are interpolated in log space, equivalent to piecewise
constant continuous rates.

Queries outside `[min T, max T]` throw `std::domain_error`; the surface never
silently extrapolates.

## Pricing

For a strike `K`, expiry `T`, forward `F = F(T)`, discount factor `df` and total
variance `w = w(log(K/F), T)`:

```text
d1 = (log(F/K) + w/2) / sqrt(w)
d2 = d1 - sqrt(w)

call = df (F N(d1) - K N(d2))
put  = df (K N(-d2) - F N(-d1))
```

If `w = 0`, the option is worth its discounted forward intrinsic value. Using
the forward and discount factor directly means the surface price does not need
a separate spot or dividend yield.

## Calendar arbitrage

Calendar arbitrage exists when total variance decreases with maturity for some
log-moneyness. The check samples each adjacent pair on a grid:

```text
w_{i+1}(k) - w_i(k) >= -tolerance      for k in [-8, 8]
```

`is_calendar_arbitrage_free(tolerance, samples)` returns false on the first
violation. This is a numerical check over the specified grid, not a symbolic
proof over the entire real line.

## Butterfly arbitrage

Gatheral's density function for total variance `w(k)` is

```text
g(k) = (1 - k w' / (2w))^2
       - (w'^2 / 4) (1/w + 1/4)
       + w'' / 2
```

A negative `g(k)` means a negative risk-neutral density. The implementation
computes derivatives with central differences at a configurable step and
samples `g(k)` on `[-8, 8]`.

`is_butterfly_arbitrage_free` checks every slice expiry and the midpoint
between adjacent expiries. The midpoint check catches violations introduced by
linear interpolation even when both endpoints are arbitrage-free.

`is_arbitrage_free` combines the calendar and butterfly checks.

## Dupire local variance

Dupire's formula converts the implied total-variance surface `w(k, T)` into a
local variance:

```text
sigma_loc^2(K, T) = (dw/dT) / g(k, T)
k = log(K / F(T))
```

`g(k, T)` is Gatheral's butterfly density function above. The implementation
differentiates total variance in time across the two bracketing slices and
computes `g` with the same central differences as the butterfly check.

```cpp
const double local_var = surface.local_variance(105.0, 0.75);
const double local_vol = surface.local_volatility(105.0, 0.75);
```

`local_variance` returns `NaN` when `g` is non-positive or the denominator is
numerically zero. `local_volatility` replaces non-finite or non-positive
variance with `sqrt(floor)`. Both require at least two slices and a non-negative
time. Times above the last slice throw. Below the first slice, total variance is
extrapolated linearly from zero, `w(k, t) = (t / t_min) w(k, t_min)`, and
Dupire's formula is applied to that extrapolated surface. This is finite at
`t = 0` and makes the integral of local variance match the surface total
variance at every calibrated expiry.

## API

```cpp
#include "cppquant/surface.hpp"

std::vector<SurfaceSlice> slices;
for (double expiry : {0.25, 0.5, 1.0, 2.0}) {
  SviParameters params = calibrated_parameters_for(expiry);
  const double forward = 100.0 * std::exp(0.03 * expiry);
  const double df = std::exp(-0.03 * expiry);
  slices.push_back(SurfaceSlice{
      SviSlice{forward, expiry, params}, df});
}

VolatilitySurface surface(std::move(slices));

const double atm_vol = surface.implied_volatility(100.0, 0.75);
const double call = surface.price(105.0, 0.75, OptionType::Call);

if (!surface.is_arbitrage_free()) {
  // inspect butterfly_g or the slice fits
}
```

## Diagnostics demo

`surface_demo` builds a calendar-arbitrage-free synthetic SVI chain, adds
configurable noise, calibrates every expiry and writes summary CSVs:

```bash
./build/surface_demo outputs
```

The output covers weight schemes, quote counts, in/out-of-sample error and
butterfly/calendar violation rates. It only uses the public API.

## Validation and tests

- A flat surface with constant total variance reproduces Black-Scholes to
  better than `1e-10` for calls and puts across strikes and expiries.
- Total variance, forward and discount factors match hand-computed
  interpolations at midpoints.
- Forward put-call parity holds: `call - put = df (F - K)`.
- A constructed decreasing-variance surface is detected as calendar
  arbitrage.
- A constructed negative-density SVI slice is detected as butterfly
  arbitrage.
- Invalid expiries, forwards, discount factors, duplicate slices, out-of-range
  queries and American options are rejected.
- A flat surface has local variance equal to the implied variance to `1e-10`;
  a two-slice term structure produces the expected time derivative, and a
  decreasing-variance surface produces negative local variance which the local
  volatility floors.
- On a skewed SVI surface, local variance matches a finite-difference Dupire
  calculation within 2% across strikes, including the linear extrapolation
  below the first slice.

## Limitations

- Arbitrage checks are grid-based; extreme wings outside `[-8, 8]` are not
  sampled.
- Linear total-variance interpolation is simple but not guaranteed to preserve
  every no-arbitrage condition between arbitrary endpoints; the midpoint
  butterfly check catches common failures.
- Local variance is piecewise constant in time between calibrated slices; it
  is not smoothed across expiries.
- No SSVI global parameterization yet.
- The surface owns forwards and discount factors; it does not build them from a
  spot, rate curve or dividend schedule.
