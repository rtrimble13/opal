// Options on stocks paying discrete cash dividends.
//
// Uses the escrowed-dividend model: the stochastic component is
// S* = S - PV(dividends paid before expiry), which follows GBM with the
// quoted volatility. Standard desk practice; for very large dividends a
// vol adjustment may be warranted (not applied here).
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "opal/core/dividends.hpp"
#include "opal/core/types.hpp"
#include "opal/engines/lattice.hpp"
#include "opal/models/black_scholes.hpp"

namespace opal {

/// European option with discrete cash dividends (escrowed-spot BSM).
/// A continuous yield q may be supplied in addition to the cash schedule.
inline double bsm_discrete_div_price(OptionType type, double S, double K,
                                     double T, double r,
                                     const DividendSchedule& divs, double sigma,
                                     double q = 0.0) {
    validate(divs, T);
    double S_esc = S - pv_dividends(divs, r, 0.0, T);
    require(S_esc > 0.0,
            "discrete dividends: PV of dividends exceeds the spot price");
    return bsm_price(type, S_esc, K, T, r, q, sigma);
}

/// American (or European) option with discrete cash dividends on a
/// Leisen-Reimer tree built on the escrowed spot. At each node the exercise
/// value uses the full spot: escrowed node value plus the PV of dividends
/// still to be paid before expiry.
inline double binomial_lr_discrete_div_price(OptionType type,
                                             ExerciseStyle style, double S,
                                             double K, double T, double r,
                                             const DividendSchedule& divs,
                                             double sigma, int steps = 501,
                                             double q = 0.0) {
    validate(divs, T);
    require(steps >= 3, "leisen-reimer: steps must be >= 3");
    require(S > 0.0 && K > 0.0 && sigma > 0.0 && T > 0.0,
            "leisen-reimer: S, K, sigma, T must be positive");
    double S_esc = S - pv_dividends(divs, r, 0.0, T);
    require(S_esc > 0.0,
            "discrete dividends: PV of dividends exceeds the spot price");
    if (steps % 2 == 0) ++steps;
    double b = r - q;
    double dt = T / steps;
    auto [d1, d2, sqrtT] = detail::d12(S_esc, K, T, b, sigma);
    (void)sqrtT;
    double p = detail::peizer_pratt(d2, steps);
    double pp = detail::peizer_pratt(d1, steps);
    double ebdt = std::exp(b * dt);
    double u = ebdt * pp / p;
    double d = (ebdt - p * u) / (1.0 - p);
    double disc = std::exp(-r * dt);
    double phi = type_sign(type);

    std::vector<double> v(steps + 1);
    for (int j = 0; j <= steps; ++j) {
        double ST = S_esc * std::pow(u, j) * std::pow(d, steps - j);
        v[j] = std::max(phi * (ST - K), 0.0);
    }
    for (int n = steps - 1; n >= 0; --n) {
        double t_n = n * dt;
        double pv_rem = pv_dividends(divs, r, t_n, T);
        for (int j = 0; j <= n; ++j) {
            v[j] = disc * (p * v[j + 1] + (1.0 - p) * v[j]);
            if (style == ExerciseStyle::American) {
                double Sn = S_esc * std::pow(u, j) * std::pow(d, n - j) + pv_rem;
                v[j] = std::max(v[j], phi * (Sn - K));
            }
        }
    }
    return v[0];
}

}  // namespace opal
