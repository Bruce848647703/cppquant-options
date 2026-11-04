#include "cppquant/options.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>

using namespace cppquant;

namespace {

struct Config {
  OptionSpec spec{OptionType::Call,
                  ExerciseStyle::European,
                  100.0,
                  100.0,
                  1.0,
                  0.2,
                  0.05,
                  0.0};
  int steps{500};
  double market_price{std::nan("")};
  std::size_t monte_carlo_paths{200'000};
  std::uint64_t seed{42};
  bool run_monte_carlo{true};
};

double next_value(int& index, int argc, char** argv, const char* flag) {
  if (index + 1 >= argc) {
    throw std::invalid_argument(std::string("missing value for ") + flag);
  }
  return std::stod(argv[++index]);
}

void print_usage(const char* program) {
  std::printf(
      "usage: %s [options]\n"
      "\n"
      "  --type call|put          option type (default call)\n"
      "  --style european|american  exercise style (default european)\n"
      "  --spot S                 spot price (default 100)\n"
      "  --strike K               strike price (default 100)\n"
      "  --time T                 time to expiry in years (default 1)\n"
      "  --vol V                  annualized volatility (default 0.2)\n"
      "  --rate R                 risk-free rate (default 0.05)\n"
      "  --dividend Q             continuous dividend yield (default 0)\n"
      "  --steps N                binomial steps (default 500)\n"
      "  --market P               market price for implied volatility\n"
      "  --mc-paths N             Monte Carlo paths (default 200000)\n"
      "  --seed N                 Monte Carlo seed (default 42)\n"
      "  --no-mc                  skip Monte Carlo\n"
      "  --help                   show this message\n",
      program);
}

Config parse(int argc, char** argv) {
  Config config;

  for (int index = 1; index < argc; ++index) {
    const std::string flag = argv[index];

    if (flag == "--help" || flag == "-h") {
      print_usage(argv[0]);
      std::exit(0);
    } else if (flag == "--type") {
      if (index + 1 >= argc) {
        throw std::invalid_argument("missing value for --type");
      }
      const std::string value = argv[++index];
      if (value == "call") {
        config.spec.type = OptionType::Call;
      } else if (value == "put") {
        config.spec.type = OptionType::Put;
      } else {
        throw std::invalid_argument("--type must be call or put");
      }
    } else if (flag == "--style") {
      if (index + 1 >= argc) {
        throw std::invalid_argument("missing value for --style");
      }
      const std::string value = argv[++index];
      if (value == "european") {
        config.spec.style = ExerciseStyle::European;
      } else if (value == "american") {
        config.spec.style = ExerciseStyle::American;
      } else {
        throw std::invalid_argument("--style must be european or american");
      }
    } else if (flag == "--spot") {
      config.spec.spot = next_value(index, argc, argv, "--spot");
    } else if (flag == "--strike") {
      config.spec.strike = next_value(index, argc, argv, "--strike");
    } else if (flag == "--time") {
      config.spec.time_to_expiry = next_value(index, argc, argv, "--time");
    } else if (flag == "--vol") {
      config.spec.volatility = next_value(index, argc, argv, "--vol");
    } else if (flag == "--rate") {
      config.spec.risk_free_rate = next_value(index, argc, argv, "--rate");
    } else if (flag == "--dividend") {
      config.spec.dividend_yield = next_value(index, argc, argv, "--dividend");
    } else if (flag == "--steps") {
      config.steps = static_cast<int>(next_value(index, argc, argv, "--steps"));
    } else if (flag == "--market") {
      config.market_price = next_value(index, argc, argv, "--market");
    } else if (flag == "--mc-paths") {
      config.monte_carlo_paths =
          static_cast<std::size_t>(next_value(index, argc, argv, "--mc-paths"));
    } else if (flag == "--seed") {
      config.seed =
          static_cast<std::uint64_t>(next_value(index, argc, argv, "--seed"));
    } else if (flag == "--no-mc") {
      config.run_monte_carlo = false;
    } else {
      throw std::invalid_argument("unknown argument: " + flag);
    }
  }

  return config;
}

} // namespace

int main(int argc, char** argv) {
  try {
    const Config config = parse(argc, argv);
    const OptionSpec& spec = config.spec;
    validate_spec(spec);

    std::printf("cppquant-options\n");
    std::printf(
        "%s %s  spot=%.4f  strike=%.4f  T=%.4f  vol=%.4f  r=%.4f  q=%.4f\n\n",
        to_string(spec.type), to_string(spec.style), spec.spot, spec.strike,
        spec.time_to_expiry, spec.volatility, spec.risk_free_rate,
        spec.dividend_yield);

    OptionSpec european = spec;
    european.style = ExerciseStyle::European;

    if (spec.style == ExerciseStyle::European) {
      const Greeks greeks = black_scholes_greeks(spec);
      std::printf("black-scholes (european)\n");
      std::printf("  price  %14.6f\n", black_scholes_price(spec));
      std::printf("  delta  %14.6f\n", greeks.delta);
      std::printf("  gamma  %14.6f\n", greeks.gamma);
      std::printf("  vega   %14.6f\n", greeks.vega);
      std::printf("  theta  %14.6f\n", greeks.theta);
      std::printf("  rho    %14.6f\n\n", greeks.rho);
    } else {
      std::printf("black-scholes european reference\n");
      std::printf("  price  %14.6f\n\n", black_scholes_price(european));
    }

    std::printf("binomial (%s, %d steps)\n", to_string(spec.style),
                config.steps);
    std::printf("  price  %14.6f\n\n", binomial_price(spec, config.steps));

    if (std::isfinite(config.market_price)) {
      try {
        const double implied =
            implied_volatility(european, config.market_price);
        std::printf("implied volatility (market %.6f)\n", config.market_price);
        std::printf("  vol    %14.6f\n\n", implied);
      } catch (const std::exception& error) {
        std::printf("implied volatility (market %.6f)\n", config.market_price);
        std::printf("  unavailable: %s\n\n", error.what());
      }
    }

    if (spec.style == ExerciseStyle::European && config.run_monte_carlo) {
      const auto result = monte_carlo_price(spec, config.monte_carlo_paths,
                                            config.seed, true, true);
      std::printf("monte carlo (%zu paths, antithetic + control variate)\n",
                  result.paths);
      std::printf("  price  %14.6f  (se %.6f, %zu samples)\n", result.price,
                  result.standard_error, result.effective_samples);
    }
  } catch (const std::exception& error) {
    std::fprintf(stderr, "error: %s\n", error.what());
    return 1;
  }

  return 0;
}
