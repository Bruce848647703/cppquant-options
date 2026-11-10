#include "cppquant/pde.hpp"

#include "cppquant/barrier.hpp"

#include "detail.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace cppquant {
namespace {

void solve_tridiagonal(const std::vector<double>& lower,
                       const std::vector<double>& diagonal,
                       const std::vector<double>& upper,
                       const std::vector<double>& rhs,
                       std::vector<double>& solution) {
  const std::size_t size = diagonal.size();
  solution.assign(size, 0.0);
  if (size == 0) {
    return;
  }

  std::vector<double> modified_upper(size, 0.0);
  std::vector<double> modified_rhs(size, 0.0);

  double pivot = diagonal[0];
  if (std::abs(pivot) < 1e-18) {
    throw std::runtime_error("singular tridiagonal system");
  }
  modified_upper[0] = upper[0] / pivot;
  modified_rhs[0] = rhs[0] / pivot;

  for (std::size_t index = 1; index < size; ++index) {
    pivot = diagonal[index] - lower[index] * modified_upper[index - 1];
    if (std::abs(pivot) < 1e-18) {
      throw std::runtime_error("singular tridiagonal system");
    }
    modified_upper[index] = upper[index] / pivot;
    modified_rhs[index] =
        (rhs[index] - lower[index] * modified_rhs[index - 1]) / pivot;
  }

  solution[size - 1] = modified_rhs[size - 1];
  for (std::size_t offset = 1; offset < size; ++offset) {
    const std::size_t index = size - 1 - offset;
    solution[index] =
        modified_rhs[index] - modified_upper[index] * solution[index + 1];
  }
}

double left_boundary(const OptionSpec& spec, double tau_remaining) {
  if (spec.type == OptionType::Call) {
    return 0.0;
  }
  if (spec.style == ExerciseStyle::American) {
    return spec.strike;
  }
  return spec.strike * std::exp(-spec.risk_free_rate * tau_remaining);
}

double right_boundary(const OptionSpec& spec, double s_max,
                      double tau_remaining) {
  if (spec.type == OptionType::Put) {
    return 0.0;
  }
  const double european =
      s_max * std::exp(-spec.dividend_yield * tau_remaining) -
      spec.strike * std::exp(-spec.risk_free_rate * tau_remaining);
  if (spec.style == ExerciseStyle::American) {
    return std::max(european, s_max - spec.strike);
  }
  return european;
}

void validate_settings(const PdeSettings& settings) {
  if (settings.space_steps < 10) {
    throw std::invalid_argument("PDE requires at least 10 space steps");
  }
  if (settings.time_steps < 1) {
    throw std::invalid_argument("PDE requires at least 1 time step");
  }
  if (!std::isfinite(settings.width_stddev) || settings.width_stddev <= 0.0) {
    throw std::invalid_argument("PDE grid width must be positive");
  }
  if (!std::isfinite(settings.local_variance_floor) ||
      settings.local_variance_floor <= 0.0) {
    throw std::invalid_argument("local variance floor must be positive");
  }
}

/// Shared Crank-Nicolson solver. `variance(spot, tau_remaining)` returns the
/// instantaneous variance at a grid node. `volatility_estimate` sizes the spot
/// grid and is not used for the coefficients.
template <typename VarianceFunction>
PdeResult solve_pde(const OptionSpec& spec, const PdeSettings& settings,
                    VarianceFunction&& variance, double volatility_estimate,
                    double spot_min = 0.0, bool zero_lower_boundary = false) {
  PdeResult result;
  result.space_steps = settings.space_steps;
  result.time_steps = settings.time_steps;

  if (spec.time_to_expiry == 0.0) {
    result.price = intrinsic_value(spec);
    const bool call = spec.type == OptionType::Call;
    result.delta = call ? (spec.spot > spec.strike ? 1.0 : 0.0)
                        : (spec.spot < spec.strike ? -1.0 : 0.0);
    return result;
  }

  const double time = spec.time_to_expiry;
  const double width =
      std::max(spec.spot, spec.strike) *
      std::exp(settings.width_stddev * std::max(volatility_estimate, 1e-6) *
               std::sqrt(time));
  double s_max = std::max(width, spec.spot * 1.01);
  s_max = std::max(s_max, spot_min + (spec.spot - spot_min) * 1.01);
  const std::size_t space_steps = settings.space_steps;
  const std::size_t time_steps = settings.time_steps;
  double ds = (s_max - spot_min) / static_cast<double>(space_steps);

  // Put the spot exactly on a grid node. This removes the O(dS) error from
  // interpolating the price and Greeks near the strike. The lower grid point is
  // either zero or an exact barrier.
  const auto spot_node =
      static_cast<std::size_t>(std::llround((spec.spot - spot_min) / ds));
  if (spot_node >= 1 && spot_node < space_steps) {
    ds = (spec.spot - spot_min) / static_cast<double>(spot_node);
    s_max = spot_min + ds * static_cast<double>(space_steps);
  }
  const double dt = time / static_cast<double>(time_steps);
  const double floor = settings.local_variance_floor;

  auto node_variance = [&](double spot, double tau_remaining) {
    double value = variance(spot, tau_remaining);
    if (!std::isfinite(value) || value <= floor) {
      ++result.clamped_local_variance_points;
      return floor;
    }
    return value;
  };

  struct TimeStep {
    double size{};
    double theta{};
  };
  std::vector<TimeStep> schedule;
  schedule.reserve(time_steps + 1);
  // Rannacher smoothing: two half-sized implicit-Euler steps damp the
  // payoff kink at the strike before switching to Crank-Nicolson.
  schedule.push_back(TimeStep{0.5 * dt, 1.0});
  schedule.push_back(TimeStep{0.5 * dt, 1.0});
  for (std::size_t index = 1; index < time_steps; ++index) {
    schedule.push_back(TimeStep{dt, 0.5});
  }

  const std::size_t interior = space_steps - 1;
  std::vector<double> lower(interior, 0.0);
  std::vector<double> diagonal(interior, 0.0);
  std::vector<double> upper(interior, 0.0);
  std::vector<double> rhs(interior, 0.0);
  std::vector<double> solution;

  std::vector<double> a(space_steps + 1, 0.0);
  std::vector<double> b(space_steps + 1, 0.0);
  std::vector<double> c(space_steps + 1, 0.0);

  std::vector<double> values(space_steps + 1, 0.0);
  for (std::size_t index = 0; index <= space_steps; ++index) {
    values[index] = detail::option_payoff(
        spec.type, spot_min + static_cast<double>(index) * ds, spec.strike);
  }
  if (zero_lower_boundary) {
    values[0] = 0.0;
  }

  double tau = 0.0;
  for (const TimeStep& step : schedule) {
    tau += step.size;
    const double tau_remaining = std::max(0.0, time - tau);
    const double new_left =
        zero_lower_boundary ? 0.0 : left_boundary(spec, tau);
    const double new_right = right_boundary(spec, s_max, tau);
    const double theta = step.theta;
    const double explicit_weight = 1.0 - theta;

    for (std::size_t index = 1; index < space_steps; ++index) {
      const double spot = spot_min + static_cast<double>(index) * ds;
      const double diffusion =
          0.5 * node_variance(spot, tau_remaining) * spot * spot / (ds * ds);
      const double convection =
          0.5 * (spec.risk_free_rate - spec.dividend_yield) * spot / ds;
      a[index] = diffusion - convection;
      b[index] = -2.0 * diffusion - spec.risk_free_rate;
      c[index] = diffusion + convection;
    }

    for (std::size_t row = 0; row < interior; ++row) {
      const std::size_t index = row + 1;
      diagonal[row] = 1.0 - theta * step.size * b[index];
      lower[row] = row == 0 ? 0.0 : -theta * step.size * a[index];
      upper[row] = row + 1 == interior ? 0.0 : -theta * step.size * c[index];

      rhs[row] =
          explicit_weight * step.size * a[index] * values[index - 1] +
          (1.0 + explicit_weight * step.size * b[index]) * values[index] +
          explicit_weight * step.size * c[index] * values[index + 1];
    }
    rhs[0] -= (-theta * step.size * a[1]) * new_left;
    rhs[interior - 1] -= (-theta * step.size * c[space_steps - 1]) * new_right;

    solve_tridiagonal(lower, diagonal, upper, rhs, solution);
    values[0] = new_left;
    values[space_steps] = new_right;
    for (std::size_t row = 0; row < interior; ++row) {
      values[row + 1] = solution[row];
    }

    if (spec.style == ExerciseStyle::American) {
      for (std::size_t index = 0; index <= space_steps; ++index) {
        values[index] =
            std::max(values[index],
                     detail::option_payoff(
                         spec.type, spot_min + static_cast<double>(index) * ds,
                         spec.strike));
      }
      if (zero_lower_boundary) {
        values[0] = 0.0;
      }
    }
  }

  const double grid_position = (spec.spot - spot_min) / ds;
  std::size_t lower_index = static_cast<std::size_t>(std::floor(grid_position));
  if (lower_index >= space_steps) {
    lower_index = space_steps - 1;
  }
  const double weight = grid_position - static_cast<double>(lower_index);
  result.price = values[lower_index] +
                 (values[lower_index + 1] - values[lower_index]) * weight;

  auto node_delta = [&](std::size_t index) {
    if (index == 0) {
      return (values[1] - values[0]) / ds;
    }
    if (index >= space_steps) {
      return (values[space_steps] - values[space_steps - 1]) / ds;
    }
    return (values[index + 1] - values[index - 1]) / (2.0 * ds);
  };
  auto node_gamma = [&](std::size_t index) {
    if (index == 0 || index >= space_steps) {
      return 0.0;
    }
    return (values[index + 1] - 2.0 * values[index] + values[index - 1]) /
           (ds * ds);
  };

  result.delta =
      node_delta(lower_index) +
      (node_delta(lower_index + 1) - node_delta(lower_index)) * weight;
  result.gamma =
      node_gamma(lower_index) +
      (node_gamma(lower_index + 1) - node_gamma(lower_index)) * weight;
  return result;
}

} // namespace

PdeResult crank_nicolson_price(const OptionSpec& spec,
                               const PdeSettings& settings) {
  validate_spec(spec);
  validate_settings(settings);

  if (spec.time_to_expiry > 0.0 && spec.volatility == 0.0) {
    PdeResult result;
    result.space_steps = settings.space_steps;
    result.time_steps = settings.time_steps;
    result.price = binomial_price(spec, static_cast<int>(settings.time_steps));
    return result;
  }

  const double variance = spec.volatility * spec.volatility;
  return solve_pde(
      spec, settings, [variance](double, double) { return variance; },
      spec.volatility);
}

PdeResult crank_nicolson_local_vol_price(const OptionSpec& spec,
                                         const VolatilitySurface& surface,
                                         const PdeSettings& settings) {
  validate_spec(spec);
  validate_settings(settings);

  if (spec.time_to_expiry == 0.0) {
    return solve_pde(spec, settings, [](double, double) { return 0.0; }, 0.0);
  }

  if (spec.time_to_expiry > 0.0) {
    if (surface.slices().size() < 2) {
      throw std::domain_error(
          "local volatility requires at least two surface expiries");
    }
    if (spec.time_to_expiry < surface.minimum_time() ||
        spec.time_to_expiry > surface.maximum_time()) {
      throw std::domain_error(
          "option expiry is outside the surface time range");
    }
  }

  auto local_variance = [&](double spot, double tau_remaining) {
    return surface.local_variance(spot, std::max(0.0, tau_remaining));
  };

  // Size the grid from implied volatility rather than the local volatility in
  // the wings. Raw SVI wings can imply local volatilities several times the ATM
  // level; using them would make S_max huge and destroy spot resolution.
  double volatility_estimate =
      surface.implied_volatility(spec.strike, spec.time_to_expiry);
  for (const auto& entry : surface.slices()) {
    const double time = entry.slice.time_to_expiry;
    if (time > spec.time_to_expiry) {
      break;
    }
    volatility_estimate =
        std::max(volatility_estimate,
                 surface.implied_volatility(entry.slice.forward, time));
  }

  return solve_pde(spec, settings, local_variance, volatility_estimate);
}

PdeResult crank_nicolson_down_and_out_price(const OptionSpec& spec,
                                            double barrier,
                                            const PdeSettings& settings) {
  validate_down_and_out(spec, barrier);
  validate_settings(settings);

  if (spec.spot <= barrier) {
    PdeResult result;
    result.space_steps = settings.space_steps;
    result.time_steps = settings.time_steps;
    return result;
  }

  const double variance = spec.volatility * spec.volatility;
  return solve_pde(
      spec, settings, [variance](double, double) { return variance; },
      spec.volatility, barrier, true);
}

} // namespace cppquant
