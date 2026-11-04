# Mathematical notes

## Black-Scholes-Merton

Under the risk-neutral measure the spot follows geometric Brownian motion

```text
dS = (r - q) S dt + sigma S dW
```

with risk-free rate `r` and continuous dividend yield `q`. For a European
option with strike `K` and time to expiry `T`,

```text
d1 = (ln(S/K) + (r - q + sigma^2 / 2) T) / (sigma sqrt(T))
d2 = d1 - sigma sqrt(T)

call = S e^(-qT) N(d1) - K e^(-rT) N(d2)
put  = K e^(-rT) N(-d2) - S e^(-qT) N(-d1)
```

`N(x) = 0.5 erfc(-x / sqrt(2))` is used instead of a polynomial
approximation, which keeps accuracy near `1e-15` for moderate arguments.

Put-call parity follows directly:

```text
call - put = S e^(-qT) - K e^(-rT)
```

## Analytic Greeks

With `phi(x) = exp(-x^2 / 2) / sqrt(2 pi)`:

```text
delta_call = e^(-qT) N(d1)
delta_put  = -e^(-qT) N(-d1)
gamma      = e^(-qT) phi(d1) / (S sigma sqrt(T))
vega       = S e^(-qT) phi(d1) sqrt(T)

theta_call = -S e^(-qT) phi(d1) sigma / (2 sqrt(T))
             - r K e^(-rT) N(d2)
             + q S e^(-qT) N(d1)

theta_put  = -S e^(-qT) phi(d1) sigma / (2 sqrt(T))
             + r K e^(-rT) N(-d2)
             - q S e^(-qT) N(-d1)

rho_call   = K T e^(-rT) N(d2)
rho_put    = -K T e^(-rT) N(-d2)
```

The implementation reports `vega` and `rho` per unit change of volatility and
rate, and `theta` per year. This is the mathematically clean convention; UIs
usually divide vega and rho by 100.

When `T = 0` the option is worth its intrinsic value and the Greeks are the
appropriate limits. When `sigma = 0` the payoff is deterministic:

```text
F = S e^((r - q) T)
call = e^(-rT) max(F - K, 0)
put  = e^(-rT) max(K - F, 0)
```

## Implied volatility

For fixed `S, K, T, r, q`, the European price is strictly increasing in
`sigma`, with limits equal to the discounted intrinsic value as `sigma -> 0`
and to `S e^(-qT)` (call) or `K e^(-rT)` (put) as `sigma -> infinity`.
Prices outside those bounds have no finite implied volatility and are rejected.

The solver brackets the root in `[0, 5]`, then runs Newton's method:

```text
sigma <- sigma - (price(sigma) - market) / vega(sigma)
```

If the step leaves the bracket or vega is too small, it falls back to the
bisection midpoint. Every iteration shrinks the bracket, so the method is
robust for deep in/out-of-the-money options while retaining Newton's fast
convergence near the root.

A price at (or numerically indistinguishable from) the lower bound is
consistent with `sigma = 0`; the solver returns zero. A price at the upper
bound implies infinite volatility and raises `std::domain_error`.

## Cox-Ross-Rubinstein tree

Split `T` into `n` steps of length `dt = T/n`:

```text
u = exp(sigma sqrt(dt))
d = 1/u
p = (exp((r - q) dt) - d) / (u - d)
discount = exp(-r dt)
```

At step `i`, node `j` (with `j` up moves) has spot

```text
S(i, j) = S u^j d^(i-j)
```

Start with terminal payoffs and roll back:

```text
V(i, j) = discount * (p V(i+1, j+1) + (1-p) V(i+1, j))
```

For American options, replace `V(i, j)` with

```text
max(V(i, j), payoff(S(i, j)))
```

The tree is computed in a single `O(n)` vector, so memory is `O(n)` and time
is `O(n^2)`. If `p` is outside `[0, 1]` the CRR parameterization is not
arbitrage-free for those inputs; the implementation rejects them instead of
clamping.

When `sigma = 0` the tree degenerates. The code instead scans the exercise
dates directly:

```text
American value = max_i e^(-r t_i) payoff(S e^((r - q) t_i))
European value = e^(-rT) payoff(S e^((r - q) T))
```

This is exact for a deterministic path.

## Monte Carlo

European payoffs are sampled from the exact GBM solution:

```text
S_T = S exp((r - q - sigma^2/2) T + sigma sqrt(T) Z),   Z ~ N(0, 1)
```

The estimator is the discounted sample mean. Two variance reductions are
available.

### Antithetic variates

Each draw `Z` is paired with `-Z`. The pair average

```text
0.5 * (payoff(Z) + payoff(-Z))
```

usually has lower variance than two independent draws because the payoff is
monotone in `Z`.

### Control variate

The terminal stock price has a known expectation:

```text
E[S_T] = S e^((r - q) T)
```

The implementation uses the discounted control `X = e^(-rT) S_T` with known
mean `S e^(-qT)` and regresses the discounted payoff `Y` on `X`:

```text
beta = Cov(Y, X) / Var(X)
Y* = Y - beta (X - E[X])
```

The adjusted mean is the reported price and the sample standard deviation of
`Y*` gives the standard error.

When antithetic variates are enabled, pair averages are formed **before**
estimating `beta`, so the control variate reduces the variance of the
estimator that is actually reported. Estimating `beta` on individual samples
and then averaging pairs can increase the standard error, because antithetic
pairing already removes much of the linear component.

### Standard error

With `n` independent samples (or `n` pair averages), the reported standard
error is

```text
se = sd(Y*) / sqrt(n)
```

For a 95% interval, use `price +/- 1.96 se`. The tests accept Monte Carlo
prices within four standard errors of the analytic value, or within `0.01`,
whichever is larger.

## Numerical edge cases

| Case | Behaviour |
| --- | --- |
| `T = 0` | Intrinsic value; delta is a step function |
| `sigma = 0` | Deterministic discounted payoff |
| `S` or `K` <= 0, negative `T` or `sigma` | `std::invalid_argument` |
| American passed to Black-Scholes or Monte Carlo | `std::invalid_argument` |
| Implied vol price outside bounds | `std::domain_error` |
| CRR probability outside `[0, 1]` | `std::domain_error` |
| Monte Carlo paths < 2 | `std::invalid_argument` |
