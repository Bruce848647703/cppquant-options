#pragma once

#include "cppquant/volatility.hpp"

#include <vector>

namespace cppquant {

/// One expiry in a volatility surface. The surface owns the forward and
/// discount factor; `OptionSpec::spot` and rates are not used for pricing.
struct SurfaceSlice {
  SviSlice slice;
  double discount_factor{1.0};
};

/// A multi-expiry SVI surface.
///
/// Total variance and log-forward are interpolated linearly in time between
/// adjacent slices. Prices are computed directly from the forward, discount
/// factor and interpolated total variance:
///
///   call = df (F N(d1) - K N(d2))
///   put  = df (K N(-d2) - F N(-d1))
///   d1 = (log(F/K) + w/2) / sqrt(w), d2 = d1 - sqrt(w)
///
/// Queries outside the slice time range throw std::domain_error rather than
/// extrapolating.
class VolatilitySurface {
public:
  explicit VolatilitySurface(std::vector<SurfaceSlice> slices);

  double minimum_time() const noexcept;
  double maximum_time() const noexcept;

  double forward(double time_to_expiry) const;
  double discount_factor(double time_to_expiry) const;

  /// Interpolated total variance at `k = log(strike / forward)`.
  double total_variance(double log_moneyness, double time_to_expiry) const;

  /// Interpolated implied volatility for a strike and expiry.
  double implied_volatility(double strike, double time_to_expiry) const;

  /// European option price from the surface forward and discount factor.
  double price(double strike, double time_to_expiry, OptionType type) const;

  /// Convenience overload. Uses the surface forward/discount factor; the
  /// option's spot is ignored. Only European options are supported.
  double price(const OptionSpec& option) const;

  /// Gatheral's butterfly density function at one point. Negative values
  /// indicate a negative risk-neutral density.
  double butterfly_g(double log_moneyness, double time_to_expiry) const;

  /// Dupire local variance at a strike and expiry. Requires at least two
  /// slices. Times below the first slice extrapolate total variance linearly
  /// from zero. Returns NaN when the butterfly denominator is non-positive or
  /// the time derivative cannot be formed.
  double local_variance(double strike, double time_to_expiry) const;

  /// Local volatility, floored at sqrt(floor). Non-finite or non-positive
  /// local variance returns sqrt(floor).
  double local_volatility(double strike, double time_to_expiry,
                          double floor = 1e-8) const;

  /// Sample total variance across adjacent slices. `samples` is the number of
  /// log-moneyness points per interval (at least 2).
  bool is_calendar_arbitrage_free(double tolerance = 1e-8,
                                  int samples = 401) const;

  /// Check Gatheral's g function at each slice and at the midpoint between
  /// adjacent slices.
  bool is_butterfly_arbitrage_free(double tolerance = 1e-8,
                                   int samples = 401) const;

  bool is_arbitrage_free(double tolerance = 1e-8, int samples = 401) const;

  const std::vector<SurfaceSlice>& slices() const noexcept {
    return slices_;
  }

private:
  struct Bracket {
    std::size_t lower{};
    std::size_t upper{};
    double fraction{};
  };

  Bracket bracket(double time_to_expiry) const;
  double checked_time(double time_to_expiry) const;

  std::vector<SurfaceSlice> slices_;
};

} // namespace cppquant
