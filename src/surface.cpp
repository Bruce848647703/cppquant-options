#include "cppquant/surface.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cppquant {
namespace {

double norm_cdf(double value) noexcept {
  return 0.5 * std::erfc(-value / std::sqrt(2.0));
}

} // namespace

VolatilitySurface::VolatilitySurface(std::vector<SurfaceSlice> slices)
    : slices_(std::move(slices)) {
  if (slices_.empty()) {
    throw std::invalid_argument(
        "volatility surface requires at least one slice");
  }

  for (const auto& entry : slices_) {
    const SviSlice& slice = entry.slice;
    if (!std::isfinite(slice.time_to_expiry) || slice.time_to_expiry <= 0.0) {
      throw std::invalid_argument(
          "surface expiries must be positive and finite");
    }
    if (!std::isfinite(slice.forward) || slice.forward <= 0.0) {
      throw std::invalid_argument(
          "surface forwards must be positive and finite");
    }
    if (!std::isfinite(entry.discount_factor) || entry.discount_factor <= 0.0) {
      throw std::invalid_argument(
          "surface discount factors must be positive and finite");
    }
    if (!slice.satisfies_basic_constraints()) {
      throw std::invalid_argument(
          "surface slice violates the basic SVI constraints");
    }
  }

  std::sort(slices_.begin(), slices_.end(),
            [](const SurfaceSlice& left, const SurfaceSlice& right) {
              return left.slice.time_to_expiry < right.slice.time_to_expiry;
            });

  for (std::size_t index = 1; index < slices_.size(); ++index) {
    const double previous = slices_[index - 1].slice.time_to_expiry;
    const double current = slices_[index].slice.time_to_expiry;
    if (current - previous <= 1e-12 * std::max(1.0, std::abs(current))) {
      throw std::invalid_argument("volatility surface has duplicate expiries");
    }
  }
}

double VolatilitySurface::minimum_time() const noexcept {
  return slices_.front().slice.time_to_expiry;
}

double VolatilitySurface::maximum_time() const noexcept {
  return slices_.back().slice.time_to_expiry;
}

double VolatilitySurface::checked_time(double time_to_expiry) const {
  const double minimum = minimum_time();
  const double maximum = maximum_time();
  const double slack = 1e-12 * std::max(1.0, maximum);
  if (!std::isfinite(time_to_expiry) || time_to_expiry < minimum - slack ||
      time_to_expiry > maximum + slack) {
    throw std::domain_error("time to expiry is outside the surface range");
  }
  return std::clamp(time_to_expiry, minimum, maximum);
}

VolatilitySurface::Bracket
VolatilitySurface::bracket(double time_to_expiry) const {
  const double time = checked_time(time_to_expiry);
  if (slices_.size() == 1) {
    return Bracket{0, 0, 0.0};
  }
  if (time <= slices_.front().slice.time_to_expiry) {
    return Bracket{0, 0, 0.0};
  }
  if (time >= slices_.back().slice.time_to_expiry) {
    const std::size_t last = slices_.size() - 1;
    return Bracket{last - 1, last, 1.0};
  }

  for (std::size_t index = 0; index + 1 < slices_.size(); ++index) {
    const double left = slices_[index].slice.time_to_expiry;
    const double right = slices_[index + 1].slice.time_to_expiry;
    if (time <= right) {
      return Bracket{index, index + 1, (time - left) / (right - left)};
    }
  }
  return Bracket{slices_.size() - 1, slices_.size() - 1, 0.0};
}

double VolatilitySurface::forward(double time_to_expiry) const {
  const Bracket bracket_ = bracket(time_to_expiry);
  const double left = slices_[bracket_.lower].slice.forward;
  if (bracket_.lower == bracket_.upper) {
    return left;
  }
  const double right = slices_[bracket_.upper].slice.forward;
  const double log_forward =
      std::log(left) + (std::log(right) - std::log(left)) * bracket_.fraction;
  return std::exp(log_forward);
}

double VolatilitySurface::discount_factor(double time_to_expiry) const {
  const Bracket bracket_ = bracket(time_to_expiry);
  const double left = slices_[bracket_.lower].discount_factor;
  if (bracket_.lower == bracket_.upper) {
    return left;
  }
  const double right = slices_[bracket_.upper].discount_factor;
  const double log_discount =
      std::log(left) + (std::log(right) - std::log(left)) * bracket_.fraction;
  return std::exp(log_discount);
}

double VolatilitySurface::total_variance(double log_moneyness,
                                         double time_to_expiry) const {
  if (!std::isfinite(log_moneyness)) {
    throw std::invalid_argument("log-moneyness must be finite");
  }
  const Bracket bracket_ = bracket(time_to_expiry);
  const double left =
      slices_[bracket_.lower].slice.total_variance(log_moneyness);
  if (bracket_.lower == bracket_.upper) {
    return left;
  }
  const double right =
      slices_[bracket_.upper].slice.total_variance(log_moneyness);
  return left + (right - left) * bracket_.fraction;
}

double VolatilitySurface::implied_volatility(double strike,
                                             double time_to_expiry) const {
  if (!std::isfinite(strike) || strike <= 0.0) {
    throw std::invalid_argument("strike must be positive and finite");
  }
  const double time = checked_time(time_to_expiry);
  const double log_moneyness = std::log(strike / forward(time));
  const double variance = total_variance(log_moneyness, time);
  if (variance <= 0.0) {
    return 0.0;
  }
  return std::sqrt(variance / time);
}

double VolatilitySurface::price(double strike, double time_to_expiry,
                                OptionType type) const {
  if (!std::isfinite(strike) || strike <= 0.0) {
    throw std::invalid_argument("strike must be positive and finite");
  }
  const double time = checked_time(time_to_expiry);
  const double forward_price = forward(time);
  const double discount = discount_factor(time);
  const double log_moneyness = std::log(strike / forward_price);
  const double variance = total_variance(log_moneyness, time);

  if (variance <= 1e-18) {
    return discount * (type == OptionType::Call
                           ? std::max(forward_price - strike, 0.0)
                           : std::max(strike - forward_price, 0.0));
  }

  const double sqrt_variance = std::sqrt(variance);
  const double d1 =
      (std::log(forward_price / strike) + 0.5 * variance) / sqrt_variance;
  const double d2 = d1 - sqrt_variance;

  if (type == OptionType::Call) {
    return discount * (forward_price * norm_cdf(d1) - strike * norm_cdf(d2));
  }
  return discount * (strike * norm_cdf(-d2) - forward_price * norm_cdf(-d1));
}

double VolatilitySurface::price(const OptionSpec& option) const {
  if (option.style != ExerciseStyle::European) {
    throw std::invalid_argument(
        "volatility surface prices European options only");
  }
  return price(option.strike, option.time_to_expiry, option.type);
}

double VolatilitySurface::butterfly_g(double log_moneyness,
                                      double time_to_expiry) const {
  if (!std::isfinite(log_moneyness)) {
    throw std::invalid_argument("log-moneyness must be finite");
  }
  const double time = checked_time(time_to_expiry);
  const double variance = total_variance(log_moneyness, time);
  if (variance <= 0.0) {
    return -std::numeric_limits<double>::infinity();
  }

  const double step = 1e-4 * std::max(1.0, std::abs(log_moneyness));
  const double variance_minus = total_variance(log_moneyness - step, time);
  const double variance_plus = total_variance(log_moneyness + step, time);
  const double first = (variance_plus - variance_minus) / (2.0 * step);
  const double second =
      (variance_plus - 2.0 * variance + variance_minus) / (step * step);

  const double term = 1.0 - log_moneyness * first / (2.0 * variance);
  return term * term - (first * first / 4.0) * (1.0 / variance + 0.25) +
         second / 2.0;
}

bool VolatilitySurface::is_calendar_arbitrage_free(double tolerance,
                                                   int samples) const {
  if (!std::isfinite(tolerance) || tolerance < 0.0 || samples < 2) {
    throw std::invalid_argument("invalid arbitrage check settings");
  }
  if (slices_.size() < 2) {
    return true;
  }

  for (std::size_t index = 0; index + 1 < slices_.size(); ++index) {
    for (int sample = 0; sample < samples; ++sample) {
      const double fraction =
          static_cast<double>(sample) / static_cast<double>(samples - 1);
      const double k = -8.0 + 16.0 * fraction;
      const double earlier = slices_[index].slice.total_variance(k);
      const double later = slices_[index + 1].slice.total_variance(k);
      if (later < earlier - tolerance) {
        return false;
      }
    }
  }
  return true;
}

bool VolatilitySurface::is_butterfly_arbitrage_free(double tolerance,
                                                    int samples) const {
  if (!std::isfinite(tolerance) || tolerance < 0.0 || samples < 2) {
    throw std::invalid_argument("invalid arbitrage check settings");
  }

  auto check_time = [&](double time) {
    for (int sample = 0; sample < samples; ++sample) {
      const double fraction =
          static_cast<double>(sample) / static_cast<double>(samples - 1);
      const double k = -8.0 + 16.0 * fraction;
      if (butterfly_g(k, time) < -tolerance) {
        return false;
      }
    }
    return true;
  };

  for (const auto& entry : slices_) {
    if (!check_time(entry.slice.time_to_expiry)) {
      return false;
    }
  }
  for (std::size_t index = 0; index + 1 < slices_.size(); ++index) {
    const double left = slices_[index].slice.time_to_expiry;
    const double right = slices_[index + 1].slice.time_to_expiry;
    if (!check_time(0.5 * (left + right))) {
      return false;
    }
  }
  return true;
}

bool VolatilitySurface::is_arbitrage_free(double tolerance, int samples) const {
  return is_calendar_arbitrage_free(tolerance, samples) &&
         is_butterfly_arbitrage_free(tolerance, samples);
}

} // namespace cppquant
