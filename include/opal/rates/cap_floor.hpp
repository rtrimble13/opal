// Interest rate caps and floors priced caplet-by-caplet under Black-76
// (lognormal forward rates) or Bachelier (normal forward rates).
#pragma once

#include <cmath>
#include <vector>

#include "opal/core/types.hpp"
#include "opal/models/black_scholes.hpp"
#include "opal/rates/curve.hpp"

namespace opal {

enum class RateVolType { Lognormal, Normal };

struct CapletDetail {
    double fixing;    // T1
    double payment;   // T2
    double forward;   // simple forward rate over [T1, T2]
    double price;     // caplet/floorlet PV per unit notional
};

struct CapFloorResult {
    double price = 0.0;
    std::vector<CapletDetail> caplets;
};

/// Multi-curve cap (is_cap=true) or floor on the simple forward rate:
/// forwards are projected off `forward_curve`, cashflows discounted on
/// `discount_curve` (OIS). Strike K, flat volatility `vol` (Black or normal
/// per vol_type), unit notional. Caplets fix at t, t+tau, ... and pay at
/// fixing + tau, with the first fixing at `first_fixing` (the period
/// starting today is conventionally excluded) and the last payment no later
/// than `maturity`.
inline CapFloorResult cap_floor_price(const DiscountCurve& discount_curve,
                                      const DiscountCurve& forward_curve,
                                      double K, double vol, double first_fixing,
                                      double maturity, double tau, bool is_cap,
                                      RateVolType vol_type = RateVolType::Lognormal) {
    require(tau > 0.0, "cap: tau must be positive");
    require(maturity > first_fixing, "cap: maturity must exceed first fixing");
    OptionType type = is_cap ? OptionType::Call : OptionType::Put;

    CapFloorResult res;
    for (double t1 = first_fixing; t1 + tau <= maturity + 1e-10; t1 += tau) {
        double t2 = t1 + tau;
        double F = forward_curve.forward_rate(t1, t2);
        double df = discount_curve.discount(t2);
        double undiscounted;
        if (vol_type == RateVolType::Lognormal) {
            require(F > 0.0, "cap: lognormal vol requires positive forwards; "
                             "use normal vol");
            undiscounted = gbs_price(type, F, K, t1, 0.0, 0.0, vol);
        } else {
            undiscounted = bachelier_price(type, F, K, t1, 0.0, vol);
        }
        double pv = tau * df * undiscounted;
        res.caplets.push_back({t1, t2, F, pv});
        res.price += pv;
    }
    return res;
}

/// Single-curve convenience: projection and discounting on the same curve.
inline CapFloorResult cap_floor_price(const DiscountCurve& curve, double K,
                                      double vol, double first_fixing,
                                      double maturity, double tau, bool is_cap,
                                      RateVolType vol_type = RateVolType::Lognormal) {
    return cap_floor_price(curve, curve, K, vol, first_fixing, maturity, tau,
                           is_cap, vol_type);
}

}  // namespace opal
