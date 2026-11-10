#pragma once

#include "cppquant/options.hpp"
#include "cppquant/pde.hpp"

namespace cppquant {

struct BarrierGreeks {
  double price{};
  double delta{};
  double vega{};
};

/// European down-and-out call or put with a continuous barrier.
///
/// The option is knocked out as soon as the spot touches `barrier`; there is no
/// rebate. `spec.style` must be European. Returns zero when the spot is at or
/// below the barrier. Delta and vega are central finite differences of the
/// analytic price.
BarrierGreeks down_and_out_greeks(const OptionSpec& spec, double barrier);

/// Price-only convenience overload.
double down_and_out_price(const OptionSpec& spec, double barrier);

/// Crank-Nicolson PDE for a European down-and-out option. The barrier is the
/// lower grid boundary, so it is represented exactly rather than by
/// interpolation. `barrier` must be below the spot.
PdeResult crank_nicolson_down_and_out_price(const OptionSpec& spec,
                                            double barrier,
                                            const PdeSettings& settings = {});

/// Check the common validation rules for a down-and-out barrier.
void validate_down_and_out(const OptionSpec& spec, double barrier);

} // namespace cppquant