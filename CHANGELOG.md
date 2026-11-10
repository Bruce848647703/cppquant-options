# Changelog

All notable changes to this project are documented here.

## [Unreleased]

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
