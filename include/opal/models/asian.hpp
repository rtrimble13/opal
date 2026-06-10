// Asian (average price) option analytics.
//  - Geometric average: exact closed form (continuous averaging).
//  - Arithmetic average: Turnbull-Wakeman (1991) moment-matched approximation.
#pragma once

#include <cmath>

#include "opal/core/types.hpp"
#include "opal/models/black_scholes.hpp"

namespace opal {

/// Continuously averaged geometric Asian option (exact closed form).
/// The geometric average of a GBM is lognormal with
///   sigma_G = sigma / sqrt(3),  b_G = (b - sigma^2/6) / 2.
inline double geometric_asian_price(OptionType type, double S, double K, double T,
                                    double r, double q, double sigma) {
    double b = r - q;
    double sigma_g = sigma / std::sqrt(3.0);
    double b_g = 0.5 * (b - sigma * sigma / 6.0);
    return gbs_price(type, S, K, T, r, b_g, sigma_g);
}

/// Continuously averaged arithmetic Asian option, Turnbull-Wakeman
/// approximation (matches first two moments of the average to a lognormal).
/// Accurate to a few basis points of spot for typical equity parameters.
inline double arithmetic_asian_price(OptionType type, double S, double K, double T,
                                     double r, double q, double sigma) {
    require(T > 0.0, "asian: time to expiry must be positive");
    double b = r - q;
    double vsq = sigma * sigma;

    double M1, M2;
    if (std::fabs(b) > 1e-8) {
        M1 = (std::exp(b * T) - 1.0) / (b * T);
        M2 = 2.0 * std::exp((2.0 * b + vsq) * T) /
                 ((b + vsq) * (2.0 * b + vsq) * T * T) +
             2.0 / (b * T * T) *
                 (1.0 / (2.0 * b + vsq) - std::exp(b * T) / (b + vsq));
    } else {
        // b -> 0 limits
        M1 = 1.0;
        M2 = 2.0 * std::exp(vsq * T) / (vsq * vsq * T * T) -
             2.0 * (1.0 + vsq * T) / (vsq * vsq * T * T);
    }

    double b_a = std::log(M1) / T;
    double var_a = std::log(M2) / T - 2.0 * b_a;
    if (var_a < 1e-12) var_a = 1e-12;
    double sigma_a = std::sqrt(var_a);
    return gbs_price(type, S, K, T, r, b_a, sigma_a);
}

}  // namespace opal
