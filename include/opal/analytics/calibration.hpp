// Least-squares calibration of the SABR and Heston models to quoted markets —
// the inverse of the pricers, recovering parameters from observed vols/prices.
// Both use the derivative-free Nelder-Mead minimizer (math/optimize.hpp) with
// parameters clamped into their admissible ranges inside the objective.
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "opal/core/types.hpp"
#include "opal/math/optimize.hpp"
#include "opal/models/heston.hpp"
#include "opal/models/sabr.hpp"

namespace opal {

struct SabrCalibration {
    SabrParams params;
    double rmse;       // root-mean-square vol error across the quoted strikes
    bool converged;
};

struct HestonCalibration {
    HestonParams params;
    double rmse;       // root-mean-square price error across the quotes
    bool converged;
};

namespace detail {
inline double clampv(double x, double lo, double hi) {
    return std::min(std::max(x, lo), hi);
}
}  // namespace detail

/// Fit SABR (alpha, rho, nu; beta held fixed) to a quoted lognormal-vol smile
/// at one expiry. `forward` and `expiry` are the smile's; `strikes` and
/// `market_vols` are the quotes (same length). Returns the fitted params and
/// the RMSE of the fitted vs. quoted vols.
inline SabrCalibration calibrate_sabr(double forward, double expiry,
                                      const std::vector<double>& strikes,
                                      const std::vector<double>& market_vols,
                                      double beta = 1.0) {
    require(strikes.size() == market_vols.size() && !strikes.empty(),
            "calibrate_sabr: strikes and market_vols must be non-empty and equal");
    require(beta >= 0.0 && beta <= 1.0, "calibrate_sabr: beta must be in [0, 1]");
    double n = static_cast<double>(strikes.size());

    // ATM-implied starting alpha: ATM SABR vol ~ alpha / forward^(1-beta).
    double atm = 0.0;
    for (double v : market_vols) atm += v;
    atm /= n;
    double alpha0 = atm * std::pow(forward, 1.0 - beta);

    auto unpack = [&](const std::vector<double>& x) {
        return SabrParams{detail::clampv(x[0], 1e-8, 10.0), beta,
                          detail::clampv(x[1], -0.999, 0.999),
                          detail::clampv(x[2], 1e-8, 10.0)};
    };
    auto sse = [&](const std::vector<double>& x) {
        SabrParams p = unpack(x);
        double s = 0.0;
        for (std::size_t i = 0; i < strikes.size(); ++i) {
            double d = sabr_lognormal_vol(forward, strikes[i], expiry, p) -
                       market_vols[i];
            s += d * d;
        }
        return s;
    };

    std::vector<double> x0{alpha0, -0.2, 0.4};
    std::vector<double> step{0.25 * alpha0 + 1e-4, 0.3, 0.3};
    math::OptResult r = math::nelder_mead(sse, x0, step, 1e-14, 4000);
    SabrCalibration out;
    out.params = unpack(r.x);
    out.rmse = std::sqrt(r.fval / n);
    out.converged = r.converged;
    return out;
}

/// Fit Heston {v0, kappa, theta, xi, rho} to a set of European option quotes.
/// Each quote i is option `type` struck at `strikes[i]`, expiring at
/// `expiries[i]`, with market price `market_prices[i]` (same length arrays).
/// The Feller condition (2 kappa theta >= xi^2) is applied as a soft penalty.
/// Returns the fitted params and the RMSE of fitted vs. quoted prices.
inline HestonCalibration calibrate_heston(OptionType type, double S, double r,
                                          double q,
                                          const std::vector<double>& strikes,
                                          const std::vector<double>& expiries,
                                          const std::vector<double>& market_prices,
                                          HestonParams guess = {0.04, 1.5, 0.04, 0.4,
                                                                -0.5}) {
    require(strikes.size() == expiries.size() &&
                strikes.size() == market_prices.size() && !strikes.empty(),
            "calibrate_heston: strikes, expiries, market_prices must match and be "
            "non-empty");
    double n = static_cast<double>(strikes.size());

    auto unpack = [&](const std::vector<double>& x) {
        return HestonParams{detail::clampv(x[0], 1e-6, 4.0),    // v0
                            detail::clampv(x[1], 1e-3, 30.0),   // kappa
                            detail::clampv(x[2], 1e-6, 4.0),    // theta
                            detail::clampv(x[3], 1e-3, 10.0),   // xi
                            detail::clampv(x[4], -0.999, 0.999)};
    };
    auto objective = [&](const std::vector<double>& x) {
        HestonParams p = unpack(x);
        double s = 0.0;
        for (std::size_t i = 0; i < strikes.size(); ++i) {
            double d = heston_price(type, S, strikes[i], expiries[i], r, q, p) -
                       market_prices[i];
            s += d * d;
        }
        // Soft Feller penalty: discourage 2 kappa theta < xi^2 without forbidding
        // it (the market often sits slightly outside Feller).
        double feller = p.xi * p.xi - 2.0 * p.kappa * p.theta;
        if (feller > 0.0) s += 0.01 * feller * feller;
        return s;
    };

    std::vector<double> x0{guess.v0, guess.kappa, guess.theta, guess.xi, guess.rho};
    std::vector<double> step{0.02, 0.5, 0.02, 0.2, 0.2};
    math::OptResult res = math::nelder_mead(objective, x0, step, 1e-12, 6000);
    HestonCalibration out;
    out.params = unpack(res.x);
    // Report the pure price RMSE (excluding the penalty) as the fit quality.
    double s = 0.0;
    for (std::size_t i = 0; i < strikes.size(); ++i) {
        double d = heston_price(type, S, strikes[i], expiries[i], r, q, out.params) -
                   market_prices[i];
        s += d * d;
    }
    out.rmse = std::sqrt(s / n);
    out.converged = res.converged;
    return out;
}

}  // namespace opal
