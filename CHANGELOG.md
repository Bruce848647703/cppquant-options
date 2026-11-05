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
