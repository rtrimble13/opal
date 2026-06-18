// Lattice (tree) pricing engines for European and American vanilla options:
//   - Cox-Ross-Rubinstein binomial
//   - Leisen-Reimer binomial (smooth, fast convergence; preferred)
//   - Kamrad-Ritchken style trinomial
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "opal/core/types.hpp"
#include "opal/models/black_scholes.hpp"

namespace opal {

/// Cox-Ross-Rubinstein binomial tree.
inline double binomial_crr_price(OptionType type, ExerciseStyle style, double S,
                                 double K, double T, double r, double q,
                                 double sigma, int steps = 512) {
    require(steps >= 1, "binomial: steps must be >= 1");
    require(S > 0.0 && K > 0.0 && sigma > 0.0 && T > 0.0,
            "binomial: S, K, sigma, T must be positive");
    double b = r - q;
    // The CRR tree is only arbitrage-free while e^{b*dt} stays inside (d, u),
    // i.e. while |b|*sqrt(dt) < sigma. A large carry-to-vol ratio with too few
    // steps breaks 0 < p < 1. Rather than reject otherwise-valid market inputs,
    // raise the step count to twice the no-arbitrage minimum (the bare minimum
    // leaves p on the degenerate boundary p->1; the 2x margin keeps the tree
    // well inside the admissible region) and cap it to bound work. Note CRR is
    // a slow O(1/n) cross-check engine: in extreme carry/vol regimes prefer
    // Leisen-Reimer or the analytic price, or raise --steps for convergence.
    if (b != 0.0) {
        double need = 2.0 * T * b * b / (sigma * sigma);
        int min_steps = static_cast<int>(std::ceil(need)) + 1;
        if (min_steps > steps) steps = std::min(min_steps, 2'000'000);
    }
    double dt = T / steps;
    double u = std::exp(sigma * std::sqrt(dt));
    double d = 1.0 / u;
    double p = (std::exp(b * dt) - d) / (u - d);
    require(p > 0.0 && p < 1.0,
            "binomial: arbitrage in tree even at the step cap; check inputs "
            "(extreme rate-to-volatility ratio)");
    double disc = std::exp(-r * dt);
    double phi = type_sign(type);

    std::vector<double> v(steps + 1);
    // terminal payoffs; node j has price S * u^(2j - steps)
    for (int j = 0; j <= steps; ++j) {
        double ST = S * std::pow(u, 2 * j - steps);
        v[j] = std::max(phi * (ST - K), 0.0);
    }
    for (int n = steps - 1; n >= 0; --n) {
        for (int j = 0; j <= n; ++j) {
            v[j] = disc * (p * v[j + 1] + (1.0 - p) * v[j]);
            if (style == ExerciseStyle::American) {
                double Sn = S * std::pow(u, 2 * j - n);
                v[j] = std::max(v[j], phi * (Sn - K));
            }
        }
    }
    return v[0];
}

namespace detail {
/// Peizer-Pratt method-2 inversion used by Leisen-Reimer.
inline double peizer_pratt(double z, int n) {
    double t = z / (n + 1.0 / 3.0 + 0.1 / (n + 1.0));
    double sign = (z >= 0.0) ? 1.0 : -1.0;
    return 0.5 + sign * 0.5 *
                     std::sqrt(1.0 - std::exp(-t * t * (n + 1.0 / 6.0)));
}
}  // namespace detail

/// Leisen-Reimer binomial tree (steps forced odd). Converges ~O(1/n^2) with
/// no oscillation; the institutional default for American vanillas.
inline double binomial_lr_price(OptionType type, ExerciseStyle style, double S,
                                double K, double T, double r, double q,
                                double sigma, int steps = 251) {
    require(steps >= 3, "leisen-reimer: steps must be >= 3");
    require(S > 0.0 && K > 0.0 && sigma > 0.0 && T > 0.0,
            "leisen-reimer: S, K, sigma, T must be positive");
    if (steps % 2 == 0) ++steps;
    double b = r - q;
    double dt = T / steps;
    auto [d1, d2, sqrtT] = detail::d12(S, K, T, b, sigma);
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
        double ST = S * std::pow(u, j) * std::pow(d, steps - j);
        v[j] = std::max(phi * (ST - K), 0.0);
    }
    for (int n = steps - 1; n >= 0; --n) {
        for (int j = 0; j <= n; ++j) {
            v[j] = disc * (p * v[j + 1] + (1.0 - p) * v[j]);
            if (style == ExerciseStyle::American) {
                double Sn = S * std::pow(u, j) * std::pow(d, n - j);
                v[j] = std::max(v[j], phi * (Sn - K));
            }
        }
    }
    return v[0];
}

/// Trinomial tree with log-space spacing dx = sigma * sqrt(3 dt).
inline double trinomial_price(OptionType type, ExerciseStyle style, double S,
                              double K, double T, double r, double q,
                              double sigma, int steps = 400) {
    require(steps >= 1, "trinomial: steps must be >= 1");
    require(S > 0.0 && K > 0.0 && sigma > 0.0 && T > 0.0,
            "trinomial: S, K, sigma, T must be positive");
    double b = r - q;
    double nu = b - 0.5 * sigma * sigma;
    // The middle branch probability pm = 2/3 - nu^2*dt/(3*sigma^2) goes
    // negative once dt > 2*sigma^2/nu^2. Raise the step count to twice that
    // no-arbitrage minimum (keeping pm comfortably positive rather than on the
    // boundary) instead of rejecting valid inputs, capped to bound work. As
    // with CRR this is a cross-check engine; prefer Leisen-Reimer/analytic or
    // raise --steps for convergence in extreme carry/vol regimes.
    if (nu != 0.0) {
        double need = 2.0 * T * nu * nu / (2.0 * sigma * sigma);
        int min_steps = static_cast<int>(std::ceil(need)) + 1;
        if (min_steps > steps) steps = std::min(min_steps, 2'000'000);
    }
    double dt = T / steps;
    double dx = sigma * std::sqrt(3.0 * dt);
    double edx = std::exp(dx);
    double pu = 0.5 * ((sigma * sigma * dt + nu * nu * dt * dt) / (dx * dx) +
                       nu * dt / dx);
    double pd = 0.5 * ((sigma * sigma * dt + nu * nu * dt * dt) / (dx * dx) -
                       nu * dt / dx);
    double pm = 1.0 - pu - pd;
    require(pu > 0.0 && pd > 0.0 && pm > 0.0,
            "trinomial: invalid probabilities even at the step cap; check "
            "inputs (extreme rate-to-volatility ratio)");
    double disc = std::exp(-r * dt);
    double phi = type_sign(type);

    int width = 2 * steps + 1;
    std::vector<double> v(width);
    for (int j = 0; j < width; ++j) {
        double ST = S * std::pow(edx, j - steps);
        v[j] = std::max(phi * (ST - K), 0.0);
    }
    for (int n = steps - 1; n >= 0; --n) {
        int w = 2 * n + 1;
        for (int j = 0; j < w; ++j) {
            // node j at level n maps to children j, j+1, j+2 at level n+1
            v[j] = disc * (pd * v[j] + pm * v[j + 1] + pu * v[j + 2]);
            if (style == ExerciseStyle::American) {
                double Sn = S * std::pow(edx, j - n);
                v[j] = std::max(v[j], phi * (Sn - K));
            }
        }
    }
    return v[0];
}

}  // namespace opal
