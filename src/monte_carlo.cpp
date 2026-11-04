#include "detail.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cppquant {

MonteCarloResult monte_carlo_price(const OptionSpec& spec, std::size_t paths,
                                   std::uint64_t seed, bool antithetic,
                                   bool control_variate) {
  validate_spec(spec);
  if (spec.style != ExerciseStyle::European) {
    throw std::invalid_argument(
        "Monte Carlo pricing supports European options only");
  }
  if (paths < 2) {
    throw std::invalid_argument("Monte Carlo requires at least 2 paths");
  }

  const std::size_t effective = antithetic ? paths - (paths % 2) : paths;
  if (effective < 2) {
    throw std::invalid_argument("Monte Carlo requires at least 2 samples");
  }

  MonteCarloResult result;
  result.paths = paths;
  result.effective_samples = effective;

  if (spec.time_to_expiry == 0.0 || spec.volatility == 0.0) {
    result.price = detail::bs_price(spec, 0.0);
    return result;
  }

  const double t = spec.time_to_expiry;
  const double drift = (spec.risk_free_rate - spec.dividend_yield -
                        0.5 * spec.volatility * spec.volatility) *
                       t;
  const double diffusion = spec.volatility * std::sqrt(t);
  const double discount = std::exp(-spec.risk_free_rate * t);
  const double expected_discounted_spot =
      spec.spot * std::exp(-spec.dividend_yield * t);

  std::vector<double> payoffs;
  std::vector<double> controls;
  payoffs.reserve(effective);
  if (control_variate) {
    controls.reserve(effective);
  }

  std::mt19937_64 rng(seed);
  std::normal_distribution<double> normal(0.0, 1.0);

  auto add_sample = [&](double z) {
    const double terminal_spot = spec.spot * std::exp(drift + diffusion * z);
    payoffs.push_back(discount * detail::option_payoff(spec.type, terminal_spot,
                                                       spec.strike));
    if (control_variate) {
      controls.push_back(discount * terminal_spot);
    }
  };

  if (antithetic) {
    for (std::size_t index = 0; index < effective / 2; ++index) {
      const double z = normal(rng);
      add_sample(z);
      add_sample(-z);
    }
  } else {
    for (std::size_t index = 0; index < effective; ++index) {
      add_sample(normal(rng));
    }
  }

  std::vector<double> sample_payoffs;
  std::vector<double> sample_controls;
  if (antithetic) {
    const std::size_t pairs = effective / 2;
    sample_payoffs.reserve(pairs);
    if (control_variate) {
      sample_controls.reserve(pairs);
    }
    for (std::size_t index = 0; index < pairs; ++index) {
      sample_payoffs.push_back(0.5 *
                               (payoffs[2 * index] + payoffs[2 * index + 1]));
      if (control_variate) {
        sample_controls.push_back(
            0.5 * (controls[2 * index] + controls[2 * index + 1]));
      }
    }
  } else {
    sample_payoffs = std::move(payoffs);
    if (control_variate) {
      sample_controls = std::move(controls);
    }
  }

  const std::size_t samples = sample_payoffs.size();
  std::vector<double> adjusted = sample_payoffs;
  if (control_variate) {
    double mean_payoff = 0.0;
    double mean_control = 0.0;
    for (std::size_t index = 0; index < samples; ++index) {
      mean_payoff += sample_payoffs[index];
      mean_control += sample_controls[index];
    }
    mean_payoff /= static_cast<double>(samples);
    mean_control /= static_cast<double>(samples);

    double covariance = 0.0;
    double control_variance = 0.0;
    for (std::size_t index = 0; index < samples; ++index) {
      const double dy = sample_payoffs[index] - mean_payoff;
      const double dx = sample_controls[index] - mean_control;
      covariance += dy * dx;
      control_variance += dx * dx;
    }
    covariance /= static_cast<double>(samples - 1);
    control_variance /= static_cast<double>(samples - 1);

    const double beta =
        control_variance > 0.0 ? covariance / control_variance : 0.0;
    for (std::size_t index = 0; index < samples; ++index) {
      adjusted[index] =
          sample_payoffs[index] -
          beta * (sample_controls[index] - expected_discounted_spot);
    }
  }

  double mean = 0.0;
  for (const double value : adjusted) {
    mean += value;
  }
  mean /= static_cast<double>(samples);

  double standard_error = 0.0;
  if (samples >= 2) {
    double variance = 0.0;
    for (const double value : adjusted) {
      const double deviation = value - mean;
      variance += deviation * deviation;
    }
    variance /= static_cast<double>(samples - 1);
    standard_error = std::sqrt(variance / static_cast<double>(samples));
  }

  result.price = mean;
  result.standard_error = standard_error;
  return result;
}

} // namespace cppquant
