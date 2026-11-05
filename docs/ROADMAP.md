# cppquant-options roadmap

8-week plan: **2026-11-09 – 2027-01-03**. Goal: grow from single-option
pricing into a surface, American Monte Carlo and PDE toolkit.

## Milestones

| Milestone | Date | Scope |
| --- | --- | --- |
| v0.2 | 2026-11-22 | SVI slice, volatility surface |
| v0.3 | 2026-12-06 | Longstaff-Schwartz American Monte Carlo |
| v0.4 | 2026-12-20 | Crank-Nicolson PDE, Dupire local volatility |
| v0.5 | 2027-01-03 | barrier options, validation, release |

## Weekly plan

### W1 (11-09 – 11-15): SVI slice calibration

- Raw SVI: `w(k) = a + b(rho(k-m) + sqrt((k-m)^2 + sigma^2))`.
- Dependency-free least squares (Levenberg-Marquardt or coordinate descent).
- Constraints: `b >= 0`, `|rho| < 1`, `a + b sigma sqrt(1-rho^2) >= 0`.
- Recover synthetic parameters from noiseless quotes.

Acceptance: noiseless RMSE below `1e-6`; noisy fits within a documented
tolerance; invalid parameters rejected or projected.

### W2 (11-16 – 11-22): volatility surface

- Store SVI slices by expiry; interpolate total variance in time and
  log-moneyness.
- Calendar and butterfly arbitrage checks.
- `implied_vol(strike, T)` and `price(option)` on the surface.
- Flat surface must reproduce Black-Scholes exactly.

Acceptance: flat-surface price error below `1e-10`; constructed arbitrage is
detected; IV round trips within `1e-7`.

### W3 (11-23 – 11-29): Longstaff-Schwartz

- Reuse GBM path generation.
- Regress only in-the-money paths on `1, k, k^2`.
- Backward induction with optimal exercise tracking.
- Compare American puts with a 1000-step binomial tree.

Acceptance: price within 1% or 4 standard errors of the tree; early-exercise
premium is positive; runs are reproducible.

### W4 (11-30 – 12-06): LSM variance reduction

- European Black-Scholes control variate.
- Combine antithetic variates with the control variate on pair averages.
- Study basis choice and in-sample/out-of-sample overfitting.
- Convergence run from 10k to 1M paths.

Acceptance: standard error falls by at least 50% versus plain LSM; price
converges without visible bias; results are recorded in `docs/BENCHMARKS.md`.

### W5 (12-07 – 12-13): Crank-Nicolson PDE

- Thomas tridiagonal solver.
- Crank-Nicolson discretization of the Black-Scholes PDE with boundary
  conditions.
- American exercise projection at every time step.
- Space/time convergence-order tests.

Acceptance: European error below `1e-4` on a fine grid; American within 1% of
the binomial tree; convergence close to second order.

### W6 (12-14 – 12-20): Dupire local volatility

- Differentiate total variance in time and log-moneyness.
- Smooth or floor negative local variances and record diagnostics.
- Price with the local-vol PDE; flat surface must recover Black-Scholes.

Acceptance: flat-surface error below `1e-4`; synthetic non-flat cases within a
documented tolerance; a clear fallback when the surface is unstable.

### W7 (12-21 – 12-27): barrier options

- Analytic down-and-out call/put under Black-Scholes.
- Cross-check with the PDE solver; refine the grid near the barrier.
- Add delta and vega.

Acceptance: analytic versus PDE within 1%; barrier limits recover vanilla
options; boundary handling is tested.

### W8 (12-28 – 01-03): release

- Collate accuracy, convergence and performance results.
- Document surface, LSM, PDE, local volatility and barriers.
- Tag `v0.5.0` and update `CHANGELOG.md`.

## Cut lines

1. SVI surface and Longstaff-Schwartz.
2. Crank-Nicolson PDE and barrier options.
3. Dupire local volatility (highest numerical risk).

## Working agreement

- No backdated commits or empty changes.
- Every numerical feature ships with a convergence or reference test.
- PRs merge only with green Linux CI and ASan/UBSan.
- Update `CHANGELOG.md` at every milestone.