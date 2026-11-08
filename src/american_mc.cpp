#include "cppquant/american_mc.hpp"

#include "detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cppquant {
namespace {

bool solve_linear(std::vector<std::vector<double>> matrix,
                  std::vector<double> rhs, std::vector<double>& solution) {
  const std::size_t size = matrix.size();
  for (std::size_t column = 0; column < size; ++column) {
    std::size_t pivot = column;
    double best = std::abs(matrix[column][column]);
    for (std::size_t row = column + 1; row < size; ++row) {
      if (std::abs(matrix[row][column]) > best) {
        best = std::abs(matrix[row][column]);
        pivot = row;
      }
    }
    if (best < 1e-12) {
      return false;
    }
    std::swap(matrix[column], matrix[pivot]);
    std::swap(rhs[column], rhs[pivot]);

    const double diagonal = matrix[column][column];
    for (std::size_t row = column + 1; row < size; ++row) {
      const double factor = matrix[row][column] / diagonal;
      if (factor == 0.0) {
        continue;
      }
      for (std::size_t entry = column; entry < size; ++entry) {
        matrix[row][entry] -= factor * matrix[column][entry];
      }
      rhs[row] -= factor * rhs[column];
    }
  }

  solution.assign(size, 0.0);
  for (std::size_t offset = 0; offset < size; ++offset) {
    const std::size_t row = size - 1 - offset;
    double value = rhs[row];
    for (std::size_t column = row + 1; column < size; ++column) {
      value -= matrix[row][column] * solution[column];
    }
    solution[row] = value / matrix[row][row];
    if (!std::isfinite(solution[row])) {
      return false;
    }
  }
  return true;
}

bool fit_continuation(const std::vector<double>& moneyness,
                      const std::vector<double>& discounted_cashflow,
                      int degree, std::vector<double>& coefficients) {
  const auto size = static_cast<std::size_t>(degree + 1);
  if (moneyness.size() <= size) {
    return false;
  }

  std::vector<std::vector<double>> normal(size, std::vector<double>(size, 0.0));
  std::vector<double> rhs(size, 0.0);

  for (std::size_t index = 0; index < moneyness.size(); ++index) {
    std::array<double, 4> powers{1.0, 0.0, 0.0, 0.0};
    for (std::size_t power = 1; power < size; ++power) {
      powers[power] = powers[power - 1] * moneyness[index];
    }
    for (std::size_t row = 0; row < size; ++row) {
      rhs[row] += powers[row] * discounted_cashflow[index];
      for (std::size_t column = 0; column < size; ++column) {
        normal[row][column] += powers[row] * powers[column];
      }
    }
  }

  return solve_linear(std::move(normal), std::move(rhs), coefficients);
}

} // namespace

AmericanMonteCarloResult
longstaff_schwartz_price(const OptionSpec& spec,
                         const LongstaffSchwartzSettings& settings) {
  validate_spec(spec);
  if (spec.style != ExerciseStyle::American) {
    throw std::invalid_argument(
        "Longstaff-Schwartz requires an American option");
  }
  if (settings.paths < 2) {
    throw std::invalid_argument("Longstaff-Schwartz requires at least 2 paths");
  }
  if (settings.exercise_steps < 1) {
    throw std::invalid_argument(
        "Longstaff-Schwartz requires at least 1 exercise step");
  }
  if (settings.polynomial_degree < 1 || settings.polynomial_degree > 3) {
    throw std::invalid_argument("polynomial degree must be between 1 and 3");
  }

  OptionSpec european = spec;
  european.style = ExerciseStyle::European;

  AmericanMonteCarloResult result;
  result.paths = settings.paths;
  result.exercise_steps = settings.exercise_steps;

  if (spec.time_to_expiry == 0.0) {
    result.price = intrinsic_value(spec);
    result.european_price = result.price;
    return result;
  }

  if (spec.volatility == 0.0) {
    result.price =
        binomial_price(spec, static_cast<int>(settings.exercise_steps));
    result.european_price = black_scholes_price(european);
    result.early_exercise_premium = result.price - result.european_price;
    return result;
  }

  const std::size_t effective = settings.antithetic
                                    ? settings.paths - settings.paths % 2
                                    : settings.paths;
  result.effective_paths = effective;

  const std::size_t steps = settings.exercise_steps;
  const std::size_t stride = steps + 1;
  const double time = spec.time_to_expiry;
  const double step_time = time / static_cast<double>(steps);
  const double drift = (spec.risk_free_rate - spec.dividend_yield -
                        0.5 * spec.volatility * spec.volatility) *
                       step_time;
  const double diffusion = spec.volatility * std::sqrt(step_time);

  std::vector<double> spots(effective * stride, 0.0);
  std::mt19937_64 generator(settings.seed);
  std::normal_distribution<double> normal(0.0, 1.0);

  if (settings.antithetic) {
    const std::size_t pairs = effective / 2;
    for (std::size_t pair = 0; pair < pairs; ++pair) {
      const std::size_t first = 2 * pair;
      const std::size_t second = first + 1;
      spots[first * stride] = spec.spot;
      spots[second * stride] = spec.spot;
      for (std::size_t step = 1; step <= steps; ++step) {
        const double draw = normal(generator);
        spots[first * stride + step] = spots[first * stride + step - 1] *
                                       std::exp(drift + diffusion * draw);
        spots[second * stride + step] = spots[second * stride + step - 1] *
                                        std::exp(drift - diffusion * draw);
      }
    }
  } else {
    for (std::size_t row = 0; row < effective; ++row) {
      spots[row * stride] = spec.spot;
      for (std::size_t step = 1; step <= steps; ++step) {
        const double draw = normal(generator);
        spots[row * stride + step] =
            spots[row * stride + step - 1] * std::exp(drift + diffusion * draw);
      }
    }
  }

  std::vector<double> cashflow(effective, 0.0);
  std::vector<std::size_t> exercise_step(effective, steps);
  for (std::size_t row = 0; row < effective; ++row) {
    cashflow[row] = detail::option_payoff(
        spec.type, spots[row * stride + steps], spec.strike);
  }

  std::vector<double> moneyness;
  std::vector<double> discounted;
  std::vector<std::size_t> itm;
  moneyness.reserve(effective);
  discounted.reserve(effective);
  itm.reserve(effective);

  for (std::size_t step = steps - 1; step >= 1; --step) {
    moneyness.clear();
    discounted.clear();
    itm.clear();

    for (std::size_t row = 0; row < effective; ++row) {
      const double spot = spots[row * stride + step];
      const double intrinsic =
          detail::option_payoff(spec.type, spot, spec.strike);
      if (intrinsic <= 0.0) {
        continue;
      }
      const double discount =
          std::exp(-spec.risk_free_rate *
                   static_cast<double>(exercise_step[row] - step) * step_time);
      moneyness.push_back(spot / spec.strike);
      discounted.push_back(cashflow[row] * discount);
      itm.push_back(row);
    }

    std::vector<double> coefficients;
    if (!fit_continuation(moneyness, discounted, settings.polynomial_degree,
                          coefficients)) {
      continue;
    }

    for (std::size_t index = 0; index < itm.size(); ++index) {
      const std::size_t row = itm[index];
      const double spot = spots[row * stride + step];
      const double intrinsic =
          detail::option_payoff(spec.type, spot, spec.strike);
      double continuation = 0.0;
      double power = 1.0;
      for (const double coefficient : coefficients) {
        continuation += coefficient * power;
        power *= spot / spec.strike;
      }
      if (intrinsic > continuation) {
        cashflow[row] = intrinsic;
        exercise_step[row] = step;
      }
    }
  }

  double sum = 0.0;
  for (std::size_t row = 0; row < effective; ++row) {
    const double discounted_cashflow =
        cashflow[row] *
        std::exp(-spec.risk_free_rate *
                 static_cast<double>(exercise_step[row]) * step_time);
    cashflow[row] = discounted_cashflow;
    sum += discounted_cashflow;
  }
  const double mean = sum / static_cast<double>(effective);

  double variance = 0.0;
  for (const double value : cashflow) {
    const double difference = value - mean;
    variance += difference * difference;
  }
  variance /= static_cast<double>(effective - 1);

  result.european_price = black_scholes_price(european);
  result.standard_error = std::sqrt(variance / static_cast<double>(effective));
  result.price = mean;

  const double immediate =
      detail::option_payoff(spec.type, spec.spot, spec.strike);
  if (immediate > result.price) {
    result.price = immediate;
    result.standard_error = 0.0;
  }
  result.early_exercise_premium = result.price - result.european_price;
  return result;
}

} // namespace cppquant
