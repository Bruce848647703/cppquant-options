# Longstaff-Schwartz American Monte Carlo

## Model

The underlying follows risk-neutral geometric Brownian motion:

```text
dS = (r - q) S dt + sigma S dW
```

The option can be exercised on a uniform grid of `exercise_steps` dates after
time zero. This is a Bermudan approximation to continuous American exercise;
the binomial tree remains the reference for convergence tests.

## Algorithm

1. Simulate `effective_paths` spot paths on the exercise grid. Antithetic mode
   pairs each normal draw `z` with `-z`, so the path count is rounded down to an
   even number.
2. At maturity, set each path's cashflow to its intrinsic value and its
   exercise date to `T`.
3. Walk backwards from the last exercise date to the first:
   - keep only in-the-money paths;
   - regress the discounted future cashflow on `1, S/K, (S/K)^2`;
   - if intrinsic value exceeds the fitted continuation value, exercise:
     replace the cashflow and move the exercise date to the current step.
4. Discount each path's final cashflow to time zero and take the mean.
5. If immediate exercise is worth more than the mean continuation estimate, use
   the intrinsic value and set the standard error to zero.

The regression is a least-squares fit solved through the normal equations with
Gaussian elimination and partial pivoting. Paths with too few in-the-money
observations skip the regression and continue.

## Why regress only in-the-money paths?

Out-of-the-money paths have zero intrinsic value and should not be exercised.
Including them biases the continuation estimate toward zero. The standard
Longstaff-Schwartz algorithm therefore fits the continuation function only on
paths where exercise is currently meaningful.

## Standard error

The reported standard error is the sample standard deviation of the discounted
cashflows divided by `sqrt(effective_paths)`. It measures Monte Carlo noise, not
the discretization bias from a finite exercise grid or from the regression
basis.

A 95% interval is approximately `price +/- 1.96 * standard_error`.

## API

```cpp
#include "cppquant/american_mc.hpp"

OptionSpec spec{OptionType::Put, ExerciseStyle::American,
                100.0, 100.0, 1.0, 0.2, 0.05, 0.0};

LongstaffSchwartzSettings settings;
settings.paths = 200'000;
settings.exercise_steps = 50;
settings.seed = 42;
settings.polynomial_degree = 2;
settings.antithetic = true;

const auto result = longstaff_schwartz_price(spec, settings);
// result.price, result.standard_error, result.european_price,
// result.early_exercise_premium, result.effective_paths
```

The function rejects European options. Use `black_scholes_price` for European
pricing and `binomial_price` when a deterministic tree is preferred.

## Degenerate cases

- `T = 0` returns intrinsic value with zero standard error.
- `volatility = 0` returns the deterministic binomial price with zero standard
  error.
- Immediate exercise dominates the estimated continuation: the result is the
  intrinsic value with zero standard error.
- `paths < 2`, `exercise_steps < 1` and `polynomial_degree` outside `[1, 3]`
  are rejected.

## Accuracy

The test suite compares against a 1,000-step binomial tree at 200,000 paths and
50 exercise steps:

| Case | LSM | Binomial | Tolerance |
| --- | ---: | ---: | --- |
| ATM American put | ~6.05 | 6.0896 | 1% or 4 standard errors |
| American call, no dividend | ~10.45 | Black-Scholes 10.4506 | 1% or 4 standard errors |
| American call, `q = 10%` | ~5.43 | 5.4462 | 1% or 4 standard errors |

The ATM put LSM price is slightly below the tree because the exercise grid is
discrete and the regression is fitted in-sample. Increasing `exercise_steps`
reduces the first source; the next milestone adds a European control variate to
reduce the second.

## Limitations

- Regression uses a low-degree polynomial in `S/K`; no Laguerre or spline basis.
- In-sample regression introduces a small bias.
- No control variate yet (planned for W4).
- Exercise dates are uniform; no discrete dividends or event dates.
- Memory is `O(effective_paths * (exercise_steps + 1))`.
- The reported standard error does not include discretization or regression
  bias.

## Tests

- ATM American put versus a 1,000-step binomial tree.
- American call without dividends versus Black-Scholes.
- American call with a dividend yield and a positive early-exercise premium.
- Reproducibility with the same seed and settings.
- Antithetic path rounding.
- Zero volatility, expiry and immediate exercise.
- Invalid option style, paths, steps and polynomial degree.
