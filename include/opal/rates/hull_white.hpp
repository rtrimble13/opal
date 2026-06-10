// Hull-White one-factor short rate model:
//   dr = a (theta(t) - r) dt + sigma dW
// fitted to an initial discount curve. Analytic prices for zero-coupon bond
// options and caplets/floorlets (Jamshidian decomposition of a caplet into a
// ZCB option).
#pragma once

#include <cmath>

#include "opal/core/types.hpp"
#include "opal/math/normal.hpp"
#include "opal/rates/curve.hpp"

namespace opal {

struct HullWhiteParams {
    double a;      // mean reversion speed (> 0)
    double sigma;  // short-rate volatility
};

namespace detail {
/// Std dev of the log ZCB price P(T1, T2) under the T1-forward measure.
inline double hw_sigma_p(const HullWhiteParams& p, double T1, double T2) {
    double B = (1.0 - std::exp(-p.a * (T2 - T1))) / p.a;
    return p.sigma * B * std::sqrt((1.0 - std::exp(-2.0 * p.a * T1)) / (2.0 * p.a));
}
}  // namespace detail

/// European option expiring T1 on a zero-coupon bond maturing T2 (> T1),
/// strike K on the bond price, unit face value.
inline double hw_zcb_option_price(OptionType type, const DiscountCurve& curve,
                                  const HullWhiteParams& p, double T1, double T2,
                                  double K) {
    require(T2 > T1 && T1 > 0.0, "hull-white: need T2 > T1 > 0");
    require(K > 0.0, "hull-white: strike must be positive");
    require(p.a > 0.0 && p.sigma > 0.0, "hull-white: a and sigma must be positive");
    double P1 = curve.discount(T1);
    double P2 = curve.discount(T2);
    double sp = detail::hw_sigma_p(p, T1, T2);
    double h = std::log(P2 / (K * P1)) / sp + 0.5 * sp;
    if (type == OptionType::Call)
        return P2 * math::norm_cdf(h) - K * P1 * math::norm_cdf(h - sp);
    return K * P1 * math::norm_cdf(-h + sp) - P2 * math::norm_cdf(-h);
}

/// Caplet on the simple forward rate for [T1, T2], strike K, unit notional,
/// paid at T2. Priced as (1 + K tau) ZCB puts (Jamshidian).
inline double hw_caplet_price(const DiscountCurve& curve, const HullWhiteParams& p,
                              double T1, double T2, double K) {
    double tau = T2 - T1;
    double Kzcb = 1.0 / (1.0 + K * tau);
    return (1.0 + K * tau) *
           hw_zcb_option_price(OptionType::Put, curve, p, T1, T2, Kzcb);
}

/// Floorlet on the simple forward rate for [T1, T2], strike K, unit notional.
inline double hw_floorlet_price(const DiscountCurve& curve,
                                const HullWhiteParams& p, double T1, double T2,
                                double K) {
    double tau = T2 - T1;
    double Kzcb = 1.0 / (1.0 + K * tau);
    return (1.0 + K * tau) *
           hw_zcb_option_price(OptionType::Call, curve, p, T1, T2, Kzcb);
}

/// Cap (sum of caplets) with first fixing at `first_fixing` and payments
/// every `tau` years until `maturity`.
inline double hw_cap_price(const DiscountCurve& curve, const HullWhiteParams& p,
                           double first_fixing, double maturity, double tau,
                           double K, bool is_cap = true) {
    require(tau > 0.0, "hull-white: tau must be positive");
    double price = 0.0;
    for (double t1 = first_fixing; t1 + tau <= maturity + 1e-10; t1 += tau) {
        double t2 = t1 + tau;
        price += is_cap ? hw_caplet_price(curve, p, t1, t2, K)
                        : hw_floorlet_price(curve, p, t1, t2, K);
    }
    return price;
}

}  // namespace opal
