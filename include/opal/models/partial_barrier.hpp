// Partial-time-start (window) barrier options: the barrier H is monitored only
// over [0, t1]; if it is not breached in that window the contract continues as
// a plain vanilla paying at T2 > t1 (a.k.a. Heynen-Kat "type A" / forward-start
// barrier). Knock-in is the complement (vanilla = knock-in + knock-out over the
// same window).
//
// Rather than the long Heynen-Kat bivariate-normal closed form (easy to get
// subtly wrong), this prices by exact one-dimensional quadrature of the [0, t1]
// "no-hit" reflection density of Y_{t1} = ln(S_{t1}/S) against the
// Black-Scholes continuation value of the surviving vanilla:
//
//   KO = e^{-r t1} * integral over y of  psi_nohit(y) * gbs(S e^y, K, T2-t1)
//
// where psi_nohit is the driftful-Brownian-motion density restricted to paths
// that never touch the barrier (image/reflection term). This uses only the
// tested gbs_price and adaptive-Simpson primitives and is validated against
// windowed-monitoring Monte Carlo.
#pragma once

#include <algorithm>
#include <cmath>

#include "opal/core/types.hpp"
#include "opal/math/normal.hpp"
#include "opal/math/solvers.hpp"
#include "opal/models/black_scholes.hpp"

namespace opal {

enum class PartialBarrierType { DownOut, UpOut, DownIn, UpIn };

inline std::string to_string(PartialBarrierType b) {
    switch (b) {
        case PartialBarrierType::DownOut: return "down-out";
        case PartialBarrierType::UpOut: return "up-out";
        case PartialBarrierType::DownIn: return "down-in";
        case PartialBarrierType::UpIn: return "up-in";
    }
    return "?";
}

namespace detail {

inline bool partial_is_down(PartialBarrierType b) {
    return b == PartialBarrierType::DownOut || b == PartialBarrierType::DownIn;
}
inline bool partial_is_in(PartialBarrierType b) {
    return b == PartialBarrierType::DownIn || b == PartialBarrierType::UpIn;
}

// Knock-out value for a barrier monitored only on [0, t1].
inline double partial_time_start_ko(OptionType type, bool down, double S, double K,
                                    double H, double t1, double T2, double r,
                                    double b, double sigma) {
    double tau = T2 - t1;
    double m = b - 0.5 * sigma * sigma;
    double sd = sigma * std::sqrt(t1);
    double mean = m * t1;
    double level = std::log(H / S);  // barrier in Y = ln(S_t/S)

    // Already knocked out at inception => worthless.
    if (down && H >= S) return 0.0;
    if (!down && H <= S) return 0.0;

    // Reflection (image) factor and mean-shift for the no-hit density.
    double img = std::exp(2.0 * m * level / (sigma * sigma));
    double shift = -2.0 * level;  // image mean is mean + shift
    auto nohit_density = [&](double y) {
        double a1 = (y - mean) / sd;
        double a2 = (y + shift - mean) / sd;
        return (math::norm_pdf(a1) - img * math::norm_pdf(a2)) / sd;
    };
    auto integrand = [&](double y) {
        return nohit_density(y) *
               gbs_price(type, S * std::exp(y), K, tau, r, b, sigma);
    };
    double lo, hi;
    if (down) {
        lo = level;             // survival requires Y_{t1} > ln(H/S)
        hi = mean + 12.0 * sd;
    } else {
        lo = mean - 12.0 * sd;
        hi = level;             // survival requires Y_{t1} < ln(H/S)
    }
    if (lo >= hi) return 0.0;
    return std::exp(-r * t1) * math::integrate(integrand, lo, hi, 1e-12);
}

}  // namespace detail

/// Partial-time-start barrier option: barrier H active only on [0, t1], vanilla
/// payoff at T2. Knock-out survives if the barrier is never touched in the
/// window; knock-in is the complement.
inline double partial_time_start_barrier_price(OptionType type,
                                               PartialBarrierType bt, double S,
                                               double K, double H, double t1,
                                               double T2, double r, double q,
                                               double sigma) {
    require(S > 0.0 && K > 0.0 && H > 0.0,
            "partial barrier: spot, strike and barrier must be positive");
    require(sigma > 0.0, "partial barrier: volatility must be positive");
    require(t1 > 0.0 && T2 > t1,
            "partial barrier: require 0 < t1 < T2 (window ends before expiry)");
    double b = r - q;
    bool down = detail::partial_is_down(bt);
    double ko = detail::partial_time_start_ko(type, down, S, K, H, t1, T2, r, b,
                                              sigma);
    if (!detail::partial_is_in(bt)) return std::max(ko, 0.0);
    // Knock-in = vanilla - knock-out (same window/barrier).
    double vanilla = gbs_price(type, S, K, T2, r, b, sigma);
    return std::max(vanilla - ko, 0.0);
}

}  // namespace opal
