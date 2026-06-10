// Model-free Greeks via finite-difference bump-and-revalue. Works with any
// pricing function of the form f(S, T, r, sigma).
#pragma once

#include <cmath>
#include <functional>

#include "opal/core/types.hpp"

namespace opal {

/// Pricer signature for numerical Greeks: price as a function of
/// (spot, time_to_expiry, rate, vol).
using PricerFn = std::function<double(double S, double T, double r, double vol)>;

struct BumpSizes {
    double spot_rel = 1e-4;   // relative spot bump
    double vol_abs = 1e-4;    // absolute vol bump
    double rate_abs = 1e-4;   // absolute rate bump
    double time_abs = 1.0 / 365.0;  // one calendar day
};

/// Central-difference Greeks for an arbitrary pricer.
inline Greeks numerical_greeks(const PricerFn& f, double S, double T, double r,
                               double vol, const BumpSizes& h = {}) {
    Greeks g;
    double ds = S * h.spot_rel;
    double base = f(S, T, r, vol);
    double up = f(S + ds, T, r, vol);
    double dn = f(S - ds, T, r, vol);
    g.price = base;
    g.delta = (up - dn) / (2.0 * ds);
    g.gamma = (up - 2.0 * base + dn) / (ds * ds);
    g.vega = (f(S, T, r, vol + h.vol_abs) - f(S, T, r, vol - h.vol_abs)) /
             (2.0 * h.vol_abs);
    g.rho = (f(S, T, r + h.rate_abs, vol) - f(S, T, r - h.rate_abs, vol)) /
            (2.0 * h.rate_abs);
    double dt = std::min(h.time_abs, 0.5 * T);
    g.theta = (f(S, T - dt, r, vol) - base) / dt;  // per year, decay convention
    return g;
}

}  // namespace opal
