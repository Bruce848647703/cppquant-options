# Changelog

All notable changes to this project are documented here.

## [0.5.0] - 2026-11-11

### Added

- Raw SVI slice parameterization and constrained Levenberg-Marquardt
  calibration (`calibrate_svi`).
- `SviSlice` helpers for total variance, implied volatility by log-moneyness
  or strike, and basic no-arbitrage checks.
- SVI tests for parameter recovery, quote round-trips, determinism, narrow
  strike ranges and invalid input.
- Dedicated SVI calibration benchmark and mathematical notes.
- `VolatilitySurface`: multi-expiry SVI slices with linear total-variance
  interpolation, log-forward and log-discount interpolation, and direct
  forward pricing.
- Calendar arbitrage checks across adjacent slices and Gatheral butterfly
  density checks at slice expiries and midpoints.
- Surface tests for flat Black-Scholes reproduction, interpolation, forward
  parity, arbitrage detection, zero variance and invalid inputs.
- Surface benchmark and mathematical notes.
- Longstaff-Schwartz American Monte Carlo with in-the-money polynomial
  regression, antithetic variates, immediate-exercise handling and standard
  errors.
- LSM tests against a 1,000-step binomial tree, including dividend-paying calls
  and zero-volatility/expiry edge cases.
- LSM benchmark and mathematical notes.
- European-payoff control variate for LSM, estimated on antithetic pair
  averages and applied with a known Black-Scholes expectation.
- LSM standard-error, path-count convergence and basis-degree tests.
- `--convergence` LSM benchmark table and updated accuracy notes.
- Crank-Nicolson PDE solver with Rannacher smoothing, Dirichlet boundaries,
  American exercise projection and finite-difference delta/gamma.
- Spot-aligned grid, tridiagonal Thomas solver and European/american PDE tests.
- PDE convergence study, benchmark and `docs/PDE.md`.
- Dupire local variance and local volatility from an SVI surface, with
  negative-density diagnostics and a per-point variance floor.
- Local-volatility Crank-Nicolson PDE sharing the Rannacher start, boundaries,
  American projection and tridiagonal solver.
- Local-vol PDE tests for flat-surface Black-Scholes agreement, term-structure
  total variance, skewed-surface reproduction, American puts, negative-variance
  flooring and validation.
- Local-vol benchmark, surface demo data generator and expanded surface/PDE
  documentation.
- Analytic European down-and-out call and put prices with delta and vega.
- Down-and-out Crank-Nicolson PDE with the barrier on the lower grid boundary.
- Barrier tests against independent reflection integrals, knockout limits,
  zero volatility, Greeks and PDE prices.
- Barrier benchmark and `docs/BARRIER.md`.
- Barrier Greek tests for expired options, a zero-volatility knock-out and a
  spot close to the barrier; removed an unreachable one-sided vega branch.
- Final release benchmark refresh, package version 0.5.0, and release notes
  covering SVI, surfaces, LSM, PDE, local volatility and barriers.

### Fixed

- PDE boundaries now discount over time-to-expiry rather than calendar time;
  deep in-the-money put/call tests cover both boundaries.

### Planned

- **0.2** — SVI slice calibration and a no-arbitrage volatility surface.
- **0.3** — Longstaff-Schwartz American Monte Carlo with variance reduction.
- **0.4** — Crank-Nicolson PDE and Dupire local volatility.
- **0.5** — barrier options, validation suite and release.

## [0.1.0] - 2026-11-04

### Added

- Black-Scholes-Merton European prices and analytic Greeks.
- Continuous dividend yield and negative-rate support.
- Bracketed Newton implied-volatility solver with no-arbitrage bounds.
- Cox-Ross-Rubinstein binomial tree with American exercise.
- Monte Carlo with antithetic variates and a terminal-stock control variate.
- CLI pricer, deterministic benchmark, tests and documentation.
