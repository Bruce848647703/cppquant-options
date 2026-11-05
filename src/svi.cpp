#include "cppquant/volatility.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cppquant {
namespace {

struct Observation {
  double log_moneyness{};
  double total_variance{};
  double weight{};
  double volatility{};
};

double svi_variance(const SviParameters& parameters,
                    double log_moneyness) noexcept {
  const double d = log_moneyness - parameters.m;
  return parameters.a +
         parameters.b *
             (parameters.rho * d +
              std::sqrt(d * d + parameters.sigma * parameters.sigma));
}

SviParameters decode(const std::array<double, 5>& state) noexcept {
  const double a_min = std::exp(std::clamp(state[0], -20.0, 20.0));
  const double b = std::exp(std::clamp(state[1], -20.0, 20.0));
  const double rho = std::tanh(state[2]);
  const double m = state[3];
  const double sigma = std::exp(std::clamp(state[4], -20.0, 20.0));
  const double a =
      a_min - b * sigma * std::sqrt(std::max(0.0, 1.0 - rho * rho));
  return SviParameters{a, b, rho, m, sigma};
}

std::array<double, 5> encode(const SviParameters& parameters) noexcept {
  const double rho = std::clamp(parameters.rho, -0.999999, 0.999999);
  const double a_min =
      parameters.a + parameters.b * parameters.sigma *
                         std::sqrt(std::max(0.0, 1.0 - rho * rho));
  return std::array<double, 5>{std::log(std::max(a_min, 1e-12)),
                               std::log(std::max(parameters.b, 1e-12)),
                               std::atanh(rho), parameters.m,
                               std::log(std::max(parameters.sigma, 1e-8))};
}

bool valid_state(const std::array<double, 5>& state,
                 double moneyness_scale) noexcept {
  for (const double value : state) {
    if (!std::isfinite(value)) {
      return false;
    }
  }
  if (state[0] < -20.0 || state[0] > 20.0 || state[1] < -20.0 ||
      state[1] > 20.0 || state[4] < -20.0 || state[4] > 20.0) {
    return false;
  }
  if (std::abs(state[2]) > 12.0) {
    return false;
  }
  return std::abs(state[3]) <= 10.0 * moneyness_scale + 1.0;
}

double cost(const SviParameters& parameters,
            const std::vector<Observation>& observations) noexcept {
  double total = 0.0;
  for (const auto& observation : observations) {
    const double difference =
        svi_variance(parameters, observation.log_moneyness) -
        observation.total_variance;
    total += observation.weight * difference * difference;
  }
  return total;
}

bool solve5(std::array<std::array<double, 5>, 5> matrix,
            std::array<double, 5> rhs, std::array<double, 5>& solution) {
  for (int column = 0; column < 5; ++column) {
    int pivot = column;
    double best = std::abs(matrix[column][column]);
    for (int row = column + 1; row < 5; ++row) {
      if (std::abs(matrix[row][column]) > best) {
        best = std::abs(matrix[row][column]);
        pivot = row;
      }
    }
    if (best < 1e-18) {
      return false;
    }
    std::swap(matrix[column], matrix[pivot]);
    std::swap(rhs[column], rhs[pivot]);

    const double diagonal = matrix[column][column];
    for (int row = column + 1; row < 5; ++row) {
      const double factor = matrix[row][column] / diagonal;
      if (factor == 0.0) {
        continue;
      }
      for (int entry = column; entry < 5; ++entry) {
        matrix[row][entry] -= factor * matrix[column][entry];
      }
      rhs[row] -= factor * rhs[column];
    }
  }

  for (int row = 4; row >= 0; --row) {
    double value = rhs[row];
    for (int column = row + 1; column < 5; ++column) {
      value -= matrix[row][column] * solution[column];
    }
    solution[row] = value / matrix[row][row];
    if (!std::isfinite(solution[row])) {
      return false;
    }
  }
  return true;
}

struct LmResult {
  std::array<double, 5> state{};
  double cost{std::numeric_limits<double>::infinity()};
  int iterations{};
  bool converged{};
};

LmResult run_levenberg_marquardt(std::array<double, 5> state,
                                 const std::vector<Observation>& observations,
                                 const SviFitSettings& settings,
                                 double moneyness_scale) {
  LmResult result;
  result.state = state;
  result.cost = cost(decode(state), observations);

  double lambda = settings.lambda_initial;
  const auto count = static_cast<int>(observations.size());

  for (; result.iterations < settings.max_iterations; ++result.iterations) {
    const SviParameters parameters = decode(state);
    std::vector<double> residual(static_cast<std::size_t>(count));
    for (int index = 0; index < count; ++index) {
      const auto& observation = observations[static_cast<std::size_t>(index)];
      residual[static_cast<std::size_t>(index)] =
          std::sqrt(observation.weight) *
          (svi_variance(parameters, observation.log_moneyness) -
           observation.total_variance);
    }

    std::vector<std::array<double, 5>> jacobian(
        static_cast<std::size_t>(count));
    for (int parameter = 0; parameter < 5; ++parameter) {
      const double step = 1e-6 * std::max(1.0, std::abs(state[parameter]));
      auto plus = state;
      auto minus = state;
      plus[static_cast<std::size_t>(parameter)] += step;
      minus[static_cast<std::size_t>(parameter)] -= step;

      if (!valid_state(plus, moneyness_scale)) {
        plus[static_cast<std::size_t>(parameter)] =
            state[static_cast<std::size_t>(parameter)];
      }
      if (!valid_state(minus, moneyness_scale)) {
        minus[static_cast<std::size_t>(parameter)] =
            state[static_cast<std::size_t>(parameter)];
      }
      const double actual_step = plus[static_cast<std::size_t>(parameter)] -
                                 minus[static_cast<std::size_t>(parameter)];
      if (actual_step == 0.0) {
        return result;
      }

      const SviParameters plus_parameters = decode(plus);
      const SviParameters minus_parameters = decode(minus);
      for (int index = 0; index < count; ++index) {
        const double k =
            observations[static_cast<std::size_t>(index)].log_moneyness;
        jacobian[static_cast<std::size_t>(index)]
                [static_cast<std::size_t>(parameter)] =
                    (svi_variance(plus_parameters, k) -
                     svi_variance(minus_parameters, k)) /
                    actual_step;
      }
    }

    std::array<std::array<double, 5>, 5> normal{};
    std::array<double, 5> gradient{};
    for (int index = 0; index < count; ++index) {
      const auto& row = jacobian[static_cast<std::size_t>(index)];
      const double value = residual[static_cast<std::size_t>(index)];
      for (int left = 0; left < 5; ++left) {
        gradient[static_cast<std::size_t>(left)] +=
            row[static_cast<std::size_t>(left)] * value;
        for (int right = 0; right < 5; ++right) {
          normal[static_cast<std::size_t>(left)]
                [static_cast<std::size_t>(right)] +=
              row[static_cast<std::size_t>(left)] *
              row[static_cast<std::size_t>(right)];
        }
      }
    }

    bool accepted = false;
    double step_norm = 0.0;
    for (int attempt = 0; attempt < 12; ++attempt) {
      auto damped = normal;
      for (int diagonal = 0; diagonal < 5; ++diagonal) {
        damped[static_cast<std::size_t>(diagonal)]
              [static_cast<std::size_t>(diagonal)] +=
            lambda * (normal[static_cast<std::size_t>(diagonal)]
                            [static_cast<std::size_t>(diagonal)] +
                      1e-12);
      }

      std::array<double, 5> rhs{};
      for (int index = 0; index < 5; ++index) {
        rhs[static_cast<std::size_t>(index)] =
            -gradient[static_cast<std::size_t>(index)];
      }

      std::array<double, 5> step{};
      if (!solve5(damped, rhs, step)) {
        lambda *= 10.0;
        continue;
      }

      auto candidate = state;
      step_norm = 0.0;
      for (int index = 0; index < 5; ++index) {
        candidate[static_cast<std::size_t>(index)] +=
            step[static_cast<std::size_t>(index)];
        step_norm += step[static_cast<std::size_t>(index)] *
                     step[static_cast<std::size_t>(index)];
      }
      step_norm = std::sqrt(step_norm);

      if (!valid_state(candidate, moneyness_scale)) {
        lambda *= 10.0;
        continue;
      }

      const double candidate_cost = cost(decode(candidate), observations);
      if (candidate_cost < result.cost) {
        state = candidate;
        result.cost = candidate_cost;
        lambda = std::max(lambda / 5.0, 1e-12);
        accepted = true;
        break;
      }
      lambda *= 10.0;
    }

    if (!accepted) {
      break;
    }

    double state_norm = 0.0;
    for (const double value : state) {
      state_norm += value * value;
    }
    state_norm = std::sqrt(state_norm);
    if (step_norm <= settings.tolerance * (1.0 + state_norm) ||
        result.cost <= 1e-24) {
      result.converged = true;
      ++result.iterations;
      break;
    }
  }

  result.state = state;
  return result;
}

SviParameters
initial_guess(const std::vector<Observation>& observations) noexcept {
  double min_variance = std::numeric_limits<double>::infinity();
  double max_variance = -std::numeric_limits<double>::infinity();
  double min_moneyness = std::numeric_limits<double>::infinity();
  double max_moneyness = -std::numeric_limits<double>::infinity();
  double moneyness_at_min = 0.0;
  double max_abs_moneyness = 0.0;

  for (const auto& observation : observations) {
    if (observation.total_variance < min_variance) {
      min_variance = observation.total_variance;
      moneyness_at_min = observation.log_moneyness;
    }
    max_variance = std::max(max_variance, observation.total_variance);
    min_moneyness = std::min(min_moneyness, observation.log_moneyness);
    max_moneyness = std::max(max_moneyness, observation.log_moneyness);
    max_abs_moneyness =
        std::max(max_abs_moneyness, std::abs(observation.log_moneyness));
  }

  const double sigma = std::max(0.05, 0.25 * (max_moneyness - min_moneyness));
  const double slope =
      (max_variance - min_variance) / std::max(1e-3, max_abs_moneyness);
  const double b = std::clamp(slope, 1e-3, 10.0);
  const double a_min = std::max(0.5 * min_variance, 1e-8);
  const double a = a_min - b * sigma;

  return SviParameters{a, b, 0.0, moneyness_at_min, sigma};
}

} // namespace

double SviSlice::total_variance(double log_moneyness) const noexcept {
  return svi_variance(parameters, log_moneyness);
}

double SviSlice::implied_volatility(double log_moneyness) const noexcept {
  if (!(time_to_expiry > 0.0)) {
    return 0.0;
  }
  const double variance = std::max(total_variance(log_moneyness), 0.0);
  return std::sqrt(variance / time_to_expiry);
}

double SviSlice::implied_volatility_for_strike(double strike) const noexcept {
  if (!(strike > 0.0) || !(forward > 0.0)) {
    return 0.0;
  }
  return implied_volatility(std::log(strike / forward));
}

bool SviSlice::satisfies_basic_constraints(double tolerance) const noexcept {
  const auto& p = parameters;
  if (!(p.b >= 0.0) || !(std::abs(p.rho) < 1.0) || !(p.sigma > 0.0)) {
    return false;
  }
  const double a_min =
      p.a + p.b * p.sigma * std::sqrt(std::max(0.0, 1.0 - p.rho * p.rho));
  return std::isfinite(a_min) && a_min >= -tolerance;
}

SviFitResult calibrate_svi(double forward, double time_to_expiry,
                           const std::vector<VolQuote>& quotes,
                           const SviFitSettings& settings) {
  if (!std::isfinite(forward) || forward <= 0.0) {
    throw std::invalid_argument("forward must be positive and finite");
  }
  if (!std::isfinite(time_to_expiry) || time_to_expiry <= 0.0) {
    throw std::invalid_argument("time to expiry must be positive and finite");
  }
  if (settings.max_iterations < 1 || !std::isfinite(settings.tolerance) ||
      settings.tolerance <= 0.0 || !std::isfinite(settings.lambda_initial) ||
      settings.lambda_initial <= 0.0) {
    throw std::invalid_argument("invalid SVI fit settings");
  }
  if (quotes.size() < 5) {
    throw std::invalid_argument("SVI calibration requires at least 5 quotes");
  }

  std::vector<Observation> observations;
  observations.reserve(quotes.size());
  std::vector<double> moneyness_values;
  moneyness_values.reserve(quotes.size());

  for (const auto& quote : quotes) {
    if (!std::isfinite(quote.strike) || quote.strike <= 0.0 ||
        !std::isfinite(quote.implied_volatility) ||
        quote.implied_volatility <= 0.0 || !std::isfinite(quote.weight) ||
        quote.weight <= 0.0) {
      throw std::invalid_argument("invalid SVI quote");
    }
    const double k = std::log(quote.strike / forward);
    observations.push_back(Observation{
        k, quote.implied_volatility * quote.implied_volatility * time_to_expiry,
        quote.weight, quote.implied_volatility});
    moneyness_values.push_back(k);
  }

  std::sort(moneyness_values.begin(), moneyness_values.end());
  std::size_t distinct = 0;
  for (const double value : moneyness_values) {
    if (distinct == 0 || std::abs(value - moneyness_values[distinct - 1]) >
                             1e-9 * std::max(1.0, std::abs(value))) {
      moneyness_values[distinct++] = value;
    }
  }
  if (distinct < 5) {
    throw std::invalid_argument(
        "SVI calibration requires at least 5 distinct strikes");
  }

  double moneyness_scale = 1.0;
  for (const auto& observation : observations) {
    moneyness_scale =
        std::max(moneyness_scale, std::abs(observation.log_moneyness));
  }

  const SviParameters base = initial_guess(observations);
  std::vector<SviParameters> starts{base};

  auto asymmetric = base;
  asymmetric.rho = 0.5;
  asymmetric.m = 0.0;
  starts.push_back(asymmetric);

  asymmetric = base;
  asymmetric.rho = -0.5;
  asymmetric.m = 0.0;
  starts.push_back(asymmetric);

  asymmetric = base;
  asymmetric.sigma *= 2.0;
  asymmetric.m = 0.0;
  starts.push_back(asymmetric);

  LmResult best;
  for (const auto& start : starts) {
    LmResult candidate = run_levenberg_marquardt(encode(start), observations,
                                                 settings, moneyness_scale);
    if (candidate.cost < best.cost) {
      best = candidate;
    }
  }

  const SviParameters parameters = decode(best.state);
  if (!std::isfinite(parameters.a) || !std::isfinite(parameters.b) ||
      !std::isfinite(parameters.rho) || !std::isfinite(parameters.m) ||
      !std::isfinite(parameters.sigma)) {
    throw std::runtime_error("SVI calibration produced non-finite parameters");
  }

  SviFitResult result;
  result.slice = SviSlice{forward, time_to_expiry, parameters};
  result.iterations = best.iterations;
  result.converged = best.converged;

  double weighted_squared_error = 0.0;
  double weight_sum = 0.0;
  for (const auto& observation : observations) {
    const double model_volatility =
        result.slice.implied_volatility(observation.log_moneyness);
    const double error = model_volatility - observation.volatility;
    weighted_squared_error += observation.weight * error * error;
    weight_sum += observation.weight;
    result.max_error = std::max(result.max_error, std::abs(error));
  }
  result.rmse = std::sqrt(weighted_squared_error / std::max(weight_sum, 1e-18));
  return result;
}

} // namespace cppquant
