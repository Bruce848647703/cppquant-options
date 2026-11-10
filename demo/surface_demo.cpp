// Surface demo: generate a synthetic option chain from a known arbitrage-free
// SVI surface, calibrate each expiry, and write the data needed to judge the
// fit (surface grid, residuals, out-of-sample error, arbitrage violations).
//
// The output is plain CSV so that plotting stays in Python. Nothing in the
// library is monkey-patched: the demo only calls the public API.

#include "cppquant/options.hpp"
#include "cppquant/surface.hpp"
#include "cppquant/volatility.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using cppquant::ExerciseStyle;
using cppquant::OptionSpec;
using cppquant::OptionType;
using cppquant::SurfaceSlice;
using cppquant::SviFitResult;
using cppquant::SviParameters;
using cppquant::SviSlice;
using cppquant::VolatilitySurface;
using cppquant::VolQuote;

constexpr double kSpot = 100.0;
constexpr double kRate = 0.03;

/// True total variance scales linearly with time, so the true surface is
/// calendar-arbitrage free by construction.
SviParameters scaled_parameters(double time) {
  // Mild upward term structure keeps the surface inside the SVI family while
  // making it visually distinguishable from a single extruded smile.
  const double scale = time * (1.0 + 0.20 * std::log1p(time));
  SviParameters parameters;
  parameters.a = 0.010 * scale;
  parameters.b = 0.15 * scale;
  parameters.rho = -0.45;
  parameters.m = 0.02;
  parameters.sigma = 0.20;
  return parameters;
}

struct TrueSlice {
  double time{};
  double forward{};
  double discount_factor{};
  SviSlice slice{};
};

std::vector<TrueSlice> build_truth(const std::vector<double>& times) {
  std::vector<TrueSlice> truth;
  truth.reserve(times.size());
  for (double time : times) {
    TrueSlice entry;
    entry.time = time;
    entry.forward = kSpot * std::exp(kRate * time);
    entry.discount_factor = std::exp(-kRate * time);
    entry.slice.forward = entry.forward;
    entry.slice.time_to_expiry = time;
    entry.slice.parameters = scaled_parameters(time);
    truth.push_back(entry);
  }
  return truth;
}

std::vector<double> linspace(double low, double high, std::size_t count) {
  if (count < 2) {
    throw std::invalid_argument("linspace needs at least two points");
  }
  std::vector<double> values(count);
  for (std::size_t i = 0; i < count; ++i) {
    const double fraction =
        static_cast<double>(i) / static_cast<double>(count - 1);
    values[i] = low + (high - low) * fraction;
  }
  return values;
}

enum class WeightScheme { Uniform, Vega, Wing };

const char* scheme_name(WeightScheme scheme) {
  switch (scheme) {
  case WeightScheme::Uniform:
    return "uniform";
  case WeightScheme::Vega:
    return "vega";
  case WeightScheme::Wing:
    return "wing";
  }
  return "unknown";
}

double vega_weight(double strike, double time, double volatility) {
  const OptionSpec spec{OptionType::Call,
                        ExerciseStyle::European,
                        kSpot,
                        strike,
                        time,
                        volatility,
                        kRate,
                        0.0};
  return std::max(cppquant::black_scholes_greeks(spec).vega, 1e-4);
}

double quote_weight(WeightScheme scheme, double log_moneyness, double strike,
                    double time, double volatility) {
  switch (scheme) {
  case WeightScheme::Uniform:
    return 1.0;
  case WeightScheme::Vega:
    return vega_weight(strike, time, volatility);
  case WeightScheme::Wing:
    return 1.0 + 3.0 * std::abs(log_moneyness);
  }
  return 1.0;
}

/// Fraction of grid points where Gatheral's density is negative.
double butterfly_violation_rate(const VolatilitySurface& surface,
                                const std::vector<double>& times) {
  const std::vector<double> grid = linspace(-1.0, 1.0, 201);
  std::size_t violations = 0;
  std::size_t points = 0;
  for (double time : times) {
    for (double log_moneyness : grid) {
      if (surface.butterfly_g(log_moneyness, time) < -1e-7) {
        ++violations;
      }
      ++points;
    }
  }
  return points == 0
             ? 0.0
             : static_cast<double>(violations) / static_cast<double>(points);
}

/// Fraction of grid points where total variance decreases with maturity.
double calendar_violation_rate(const VolatilitySurface& surface,
                               const std::vector<double>& times) {
  const std::vector<double> grid = linspace(-1.0, 1.0, 201);
  std::size_t violations = 0;
  std::size_t points = 0;
  for (std::size_t i = 0; i + 1 < times.size(); ++i) {
    for (double log_moneyness : grid) {
      const double near_variance =
          surface.total_variance(log_moneyness, times[i]);
      const double far_variance =
          surface.total_variance(log_moneyness, times[i + 1]);
      if (far_variance < near_variance - 1e-6) {
        ++violations;
      }
      ++points;
    }
  }
  return points == 0
             ? 0.0
             : static_cast<double>(violations) / static_cast<double>(points);
}

double root_mean_square(const std::vector<double>& errors) {
  if (errors.empty()) {
    return 0.0;
  }
  double sum = 0.0;
  for (double error : errors) {
    sum += error * error;
  }
  return std::sqrt(sum / static_cast<double>(errors.size()));
}

struct ScenarioResult {
  double rmse_fit{};
  double rmse_true_in{};
  double rmse_true_oos{};
  double butterfly_rate{};
  double calendar_rate{};
  bool all_converged{true};
};

struct ScenarioRun {
  ScenarioResult result;
  std::vector<SviFitResult> fits;
  std::vector<std::vector<double>> observed;
  std::optional<VolatilitySurface> surface;
};

/// Calibrate every expiry using `selected` quotes, then score the fit.
/// Strikes outside `selected` are the out-of-sample set.
ScenarioRun run_scenario(const std::vector<double>& times,
                         const std::vector<TrueSlice>& truth,
                         const std::vector<std::vector<double>>& true_iv,
                         const std::vector<double>& moneyness,
                         const std::vector<std::size_t>& selected, double noise,
                         WeightScheme scheme, std::uint64_t seed_base) {
  std::seed_seq seed{static_cast<unsigned>(seed_base)};
  std::mt19937_64 generator(seed);
  std::normal_distribution<double> normal(0.0, 1.0);

  ScenarioRun run;
  run.observed.resize(times.size());
  std::vector<double> in_sample_errors;
  std::vector<double> out_of_sample_errors;

  for (std::size_t i = 0; i < times.size(); ++i) {
    run.observed[i].resize(moneyness.size());

    std::vector<VolQuote> quotes;
    std::vector<std::size_t> quote_indices;
    quotes.reserve(selected.size());

    for (std::size_t j = 0; j < moneyness.size(); ++j) {
      const double clean = true_iv[i][j];
      const double noisy = std::max(1e-4, clean + noise * normal(generator));
      run.observed[i][j] = noisy;
    }

    for (std::size_t j : selected) {
      const double strike = truth[i].forward * std::exp(moneyness[j]);
      quotes.push_back(VolQuote{
          strike, run.observed[i][j],
          quote_weight(scheme, moneyness[j], strike, times[i], true_iv[i][j])});
      quote_indices.push_back(j);
    }

    cppquant::SviFitSettings settings;
    settings.max_iterations = 300;
    const SviFitResult fit =
        cppquant::calibrate_svi(truth[i].forward, times[i], quotes, settings);
    if (!fit.converged) {
      run.result.all_converged = false;
    }

    for (std::size_t j : selected) {
      in_sample_errors.push_back(fit.slice.implied_volatility(moneyness[j]) -
                                 true_iv[i][j]);
    }
    for (std::size_t j = 0; j < moneyness.size(); ++j) {
      if (std::find(selected.begin(), selected.end(), j) != selected.end()) {
        continue;
      }
      out_of_sample_errors.push_back(
          fit.slice.implied_volatility(moneyness[j]) - true_iv[i][j]);
    }

    run.result.rmse_fit += fit.rmse;
    run.fits.push_back(fit);
  }

  run.result.rmse_fit /= static_cast<double>(times.size());
  run.result.rmse_true_in = root_mean_square(in_sample_errors);
  run.result.rmse_true_oos = root_mean_square(out_of_sample_errors);

  std::vector<SurfaceSlice> slices;
  slices.reserve(run.fits.size());
  for (std::size_t i = 0; i < run.fits.size(); ++i) {
    slices.push_back(SurfaceSlice{run.fits[i].slice, truth[i].discount_factor});
  }
  try {
    run.surface.emplace(slices);
  } catch (const std::exception&) {
    run.surface.reset();
  }
  if (run.surface) {
    run.result.butterfly_rate = butterfly_violation_rate(*run.surface, times);
    run.result.calendar_rate = calendar_violation_rate(*run.surface, times);
  }
  return run;
}

std::ofstream open_csv(const std::filesystem::path& path) {
  std::ofstream stream(path);
  if (!stream) {
    throw std::runtime_error("cannot open " + path.string());
  }
  stream << std::setprecision(10);
  return stream;
}

/// Draw a random subset of `pool` without replacement.
std::vector<std::size_t> random_subset(const std::vector<std::size_t>& pool,
                                       std::size_t count, std::uint64_t seed) {
  std::vector<std::size_t> shuffled = pool;
  std::mt19937_64 generator(seed);
  std::shuffle(shuffled.begin(), shuffled.end(), generator);
  shuffled.resize(count);
  std::sort(shuffled.begin(), shuffled.end());
  return shuffled;
}

void write_base_details(const std::filesystem::path& output_dir,
                        const std::vector<double>& times,
                        const std::vector<TrueSlice>& truth,
                        const std::vector<std::vector<double>>& true_iv,
                        const std::vector<double>& moneyness,
                        const ScenarioRun& run) {
  {
    std::ofstream grid_csv = open_csv(output_dir / "surface_grid.csv");
    grid_csv << "time,log_moneyness,implied_vol\n";
    if (run.surface) {
      for (double time : linspace(times.front(), times.back(), 25)) {
        const double forward = run.surface->forward(time);
        for (double log_moneyness : linspace(-0.60, 0.40, 41)) {
          const double strike = forward * std::exp(log_moneyness);
          grid_csv << time << ',' << log_moneyness << ','
                   << run.surface->implied_volatility(strike, time) << '\n';
        }
      }
    }
  }
  {
    std::ofstream quote_csv = open_csv(output_dir / "market_vs_model.csv");
    quote_csv << "time,log_moneyness,strike,iv_true,iv_market,iv_fit\n";
    for (std::size_t i = 0; i < times.size(); ++i) {
      for (std::size_t j = 0; j < moneyness.size(); ++j) {
        const double strike = truth[i].forward * std::exp(moneyness[j]);
        quote_csv << times[i] << ',' << moneyness[j] << ',' << strike << ','
                  << true_iv[i][j] << ',' << run.observed[i][j] << ','
                  << run.fits[i].slice.implied_volatility(moneyness[j]) << '\n';
      }
    }
  }
  {
    std::ofstream summary_csv =
        open_csv(output_dir / "calibration_summary.csv");
    summary_csv << "time,forward,rmse,max_error,converged,iterations\n";
    for (std::size_t i = 0; i < times.size(); ++i) {
      summary_csv << times[i] << ',' << truth[i].forward << ','
                  << run.fits[i].rmse << ',' << run.fits[i].max_error << ','
                  << (run.fits[i].converged ? 1 : 0) << ','
                  << run.fits[i].iterations << '\n';
    }
  }
}

} // namespace

int main(int argc, char** argv) {
  const std::filesystem::path output_dir = argc > 1 ? argv[1] : "outputs";
  std::filesystem::create_directories(output_dir);

  const std::vector<double> times{0.25, 0.50, 1.00, 2.00};
  const std::vector<double> moneyness = linspace(-0.60, 0.40, 21);
  const std::vector<TrueSlice> truth = build_truth(times);

  // Clean implied vols per (expiry, strike); noise is added per scenario.
  std::vector<std::vector<double>> true_iv(times.size());
  for (std::size_t i = 0; i < times.size(); ++i) {
    true_iv[i].resize(moneyness.size());
    for (std::size_t j = 0; j < moneyness.size(); ++j) {
      true_iv[i][j] = truth[i].slice.implied_volatility(moneyness[j]);
    }
  }

  // The weight sweep uses every strike except a held-out third.
  std::vector<std::size_t> dense_selection;
  for (std::size_t j = 0; j < moneyness.size(); ++j) {
    if (j % 3 != 2) {
      dense_selection.push_back(j);
    }
  }

  const std::vector<double> noise_levels{0.0, 0.002, 0.004, 0.008, 0.016};
  const std::vector<WeightScheme> schemes{
      WeightScheme::Uniform, WeightScheme::Vega, WeightScheme::Wing};

  {
    std::ofstream sweep_csv = open_csv(output_dir / "sweep_weights.csv");
    sweep_csv << "noise,scheme,rmse_fit,rmse_true_in,rmse_true_oos,"
                 "butterfly_rate,calendar_rate,all_converged\n";
    for (std::size_t ni = 0; ni < noise_levels.size(); ++ni) {
      for (std::size_t si = 0; si < schemes.size(); ++si) {
        const ScenarioRun run =
            run_scenario(times, truth, true_iv, moneyness, dense_selection,
                         noise_levels[ni], schemes[si], 7919 + 37 * ni + si);
        sweep_csv << noise_levels[ni] << ',' << scheme_name(schemes[si]) << ','
                  << run.result.rmse_fit << ',' << run.result.rmse_true_in
                  << ',' << run.result.rmse_true_oos << ','
                  << run.result.butterfly_rate << ','
                  << run.result.calendar_rate << ','
                  << (run.result.all_converged ? 1 : 0) << '\n';

        if (ni == 2 && si == 0) {
          write_base_details(output_dir, times, truth, true_iv, moneyness, run);
        }
      }
    }
  }

  // How many strikes do you need before the fitted slice stays arbitrage free?
  // A fixed hold-out set keeps the out-of-sample score comparable across
  // counts.
  std::vector<std::size_t> sparse_pool;
  for (std::size_t j = 0; j < moneyness.size(); ++j) {
    if (j % 5 != 3) {
      sparse_pool.push_back(j);
    }
  }
  const std::vector<std::size_t> quote_counts{5, 7, 9, 13, 17};
  const std::vector<double> sparse_noise{0.0, 0.004, 0.008};
  {
    std::ofstream sparse_csv = open_csv(output_dir / "sweep_sparse.csv");
    sparse_csv << "n_quotes,noise,rmse_true_in,rmse_true_oos,failure_rate,"
                  "butterfly_rate,calendar_rate,all_converged\n";
    // Repeat each setting over random strike subsets: a single subset says more
    // about which strikes were picked than about the sample size. Medians are
    // used because sparse, ill-conditioned fits can diverge catastrophically.
    constexpr int kRepeats = 40;
    constexpr double kFailureThreshold = 0.05; // 5 vol points
    for (std::size_t count : quote_counts) {
      for (std::size_t ni = 0; ni < sparse_noise.size(); ++ni) {
        std::vector<double> in_values;
        std::vector<double> oos_values;
        double butterfly = 0.0;
        double calendar = 0.0;
        int failures = 0;
        bool converged = true;
        for (int repeat = 0; repeat < kRepeats; ++repeat) {
          const std::uint64_t seed = 104729 + 53 * count + 11 * ni +
                                     7919 * static_cast<std::uint64_t>(repeat);
          const std::vector<std::size_t> selection =
              random_subset(sparse_pool, count, seed);
          const ScenarioRun run =
              run_scenario(times, truth, true_iv, moneyness, selection,
                           sparse_noise[ni], WeightScheme::Uniform, seed);
          in_values.push_back(run.result.rmse_true_in);
          oos_values.push_back(run.result.rmse_true_oos);
          if (run.result.rmse_true_oos > kFailureThreshold) {
            ++failures;
          }
          butterfly += run.result.butterfly_rate;
          calendar += run.result.calendar_rate;
          converged = converged && run.result.all_converged;
        }
        std::sort(in_values.begin(), in_values.end());
        std::sort(oos_values.begin(), oos_values.end());
        const std::size_t middle = oos_values.size() / 2;
        const double scale = 1.0 / static_cast<double>(kRepeats);
        sparse_csv << count << ',' << sparse_noise[ni] << ','
                   << in_values[middle] << ',' << oos_values[middle] << ','
                   << static_cast<double>(failures) * scale << ','
                   << butterfly * scale << ',' << calendar * scale << ','
                   << (converged ? 1 : 0) << '\n';
      }
    }
  }

  std::printf("wrote demo data to %s\n",
              std::filesystem::absolute(output_dir).c_str());
  return 0;
}
