// First-class Heston (1993) greeks: spot sensitivities plus sensitivities to
// the model's own parameters, computed by central finite differences over the
// semi-analytic `heston_price`. This is the constructive counterpart to the
// generic `numerical_greeks` path, which only understands a single lognormal
// `vol` knob and so cannot expose v0/theta/xi/rho risk.
//
// Conventions (all per 1.00 of the bumped quantity unless noted):
//   delta  = dV/dS,        gamma = d2V/dS2
//   vega   = dV/dsigma     via a *parallel shift of the variance level*: a vol
//            bump sigma -> sigma +- h maps to a variance change applied to BOTH
//            v0 and theta, anchored at the instantaneous vol sqrt(v0). This is
//            the same "vega" convention adopted for the scenario/numerical-greeks
//            tooling (issue #3): at h = 0 it is a no-op, so the base price is
//            unchanged.
//   theta  = time decay, per year, decay convention: (V(T-dt) - V(T)) / dt.
//   rho    = dV/dr.
//   dv0    = dV/dv0,  dtheta = dV/dtheta  (per unit *variance*)
//   dxi    = dV/dxi   (per unit vol-of-vol)
//   drho   = dV/drho  (per unit correlation)
#pragma once

#include <algorithm>
#include <cmath>

#include "opal/core/types.hpp"
#include "opal/models/heston.hpp"

namespace opal {

/// Full set of Heston sensitivities. Spot/rate/time greeks mirror the generic
/// `Greeks` conventions; the trailing four are sensitivities to the Heston
/// parameters themselves.
struct HestonGreeks {
    double price = 0.0;
    double delta = 0.0;   // dV/dS
    double gamma = 0.0;   // d2V/dS2
    double vega = 0.0;    // dV/dsigma via parallel variance shift (v0 & theta)
    double theta = 0.0;   // time decay, per year
    double rho = 0.0;     // dV/dr
    // Sensitivities to the Heston parameters.
    double dv0 = 0.0;     // dV/dv0     (initial variance)
    double dtheta = 0.0;  // dV/dtheta  (long-run variance)
    double dxi = 0.0;     // dV/dxi     (vol of vol)
    double drho = 0.0;    // dV/drho    (spot/vol correlation)
};

struct HestonBumpSizes {
    double spot_rel = 1e-4;         // relative spot bump
    double vol_abs = 1e-4;          // absolute vol bump for the parallel "vega"
    double rate_abs = 1e-4;         // absolute rate bump
    double time_abs = 1.0 / 365.0;  // one calendar day
    double v0_abs = 1e-4;           // absolute variance bump
    double theta_abs = 1e-4;        // absolute variance bump
    double xi_abs = 1e-4;           // absolute vol-of-vol bump
    double rho_abs = 1e-4;          // absolute correlation bump
};

/// First-class greeks for a European option under Heston.
inline HestonGreeks heston_greeks(OptionType type, double S, double K, double T,
                                  double r, double q, const HestonParams& p,
                                  const HestonBumpSizes& h = {}) {
    auto price = [&](double s, double rr, double tt, const HestonParams& hp) {
        return heston_price(type, s, K, tt, rr, q, hp);
    };

    HestonGreeks g;
    g.price = price(S, r, T, p);

    // Spot delta / gamma (central second difference).
    double ds = S * h.spot_rel;
    double up = price(S + ds, r, T, p);
    double dn = price(S - ds, r, T, p);
    g.delta = (up - dn) / (2.0 * ds);
    g.gamma = (up - 2.0 * g.price + dn) / (ds * ds);

    // Rate rho.
    g.rho = (price(S, r + h.rate_abs, T, p) - price(S, r - h.rate_abs, T, p)) /
            (2.0 * h.rate_abs);

    // Theta (per year, decay convention; clamp the step near expiry).
    double dt = std::min(h.time_abs, 0.5 * T);
    g.theta = (price(S, r, T - dt, p) - g.price) / dt;

    // Vega: parallel shift of the variance level, anchored at the instantaneous
    // vol sqrt(v0). A vol bump +-hv maps to a variance change dvar applied to
    // both v0 and theta.
    double sig = std::sqrt(std::max(p.v0, 0.0));
    double hv = h.vol_abs;
    auto shift_var = [&](double dvar) {
        HestonParams hp = p;
        hp.v0 = std::max(p.v0 + dvar, 0.0);
        hp.theta = std::max(p.theta + dvar, 0.0);
        return hp;
    };
    if (sig > hv) {
        double dvar_up = (sig + hv) * (sig + hv) - sig * sig;
        double dvar_dn = (sig - hv) * (sig - hv) - sig * sig;
        g.vega = (price(S, r, T, shift_var(dvar_up)) -
                  price(S, r, T, shift_var(dvar_dn))) /
                 (2.0 * hv);
    } else {
        // Near-zero vol: one-sided to avoid a sign flip in the squared bump.
        double dvar_up = (sig + hv) * (sig + hv) - sig * sig;
        g.vega = (price(S, r, T, shift_var(dvar_up)) - g.price) / hv;
    }

    // Heston-parameter sensitivities. Each is a central difference with the
    // bump clamped to the parameter's admissible range (v0,theta >= 0; xi > 0;
    // rho in [-1, 1]); the denominator uses the realized (up - dn) span so a
    // clamped one-sided bump still yields a consistent slope.
    auto param_deriv = [&](double base, double lo, double hi, double hb,
                           auto&& setter) {
        double pu = std::min(base + hb, hi);
        double pd = std::max(base - hb, lo);
        HestonParams hpu = p, hpd = p;
        setter(hpu, pu);
        setter(hpd, pd);
        return (price(S, r, T, hpu) - price(S, r, T, hpd)) / (pu - pd);
    };
    g.dv0 = param_deriv(p.v0, 0.0, 1e9, h.v0_abs,
                        [](HestonParams& hp, double x) { hp.v0 = x; });
    g.dtheta = param_deriv(p.theta, 0.0, 1e9, h.theta_abs,
                          [](HestonParams& hp, double x) { hp.theta = x; });
    g.dxi = param_deriv(p.xi, 1e-8, 1e9, h.xi_abs,
                       [](HestonParams& hp, double x) { hp.xi = x; });
    g.drho = param_deriv(p.rho, -0.999999, 0.999999, h.rho_abs,
                        [](HestonParams& hp, double x) { hp.rho = x; });
    return g;
}

}  // namespace opal
