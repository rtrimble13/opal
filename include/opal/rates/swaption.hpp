// European swaptions under Black-76 (lognormal swap rate), Bachelier
// (normal swap rate) and SABR-implied vols.
#pragma once

#include <cmath>

#include "opal/core/types.hpp"
#include "opal/models/black_scholes.hpp"
#include "opal/models/sabr.hpp"
#include "opal/rates/curve.hpp"

namespace opal {

enum class SwaptionType { Payer, Receiver };

struct SwaptionResult {
    double price = 0.0;
    double forward_swap_rate = 0.0;
    double annuity = 0.0;  // PV01 of the fixed leg per unit rate
};

namespace detail {
/// Multi-curve swap metrics: the annuity is discounted on the OIS curve and
/// the par swap rate equates the OIS-discounted floating leg (forwards from
/// the projection curve) to the fixed leg. With a single curve this reduces
/// to the textbook (P(T0) - P(Tn)) / annuity.
inline SwaptionResult swap_metrics(const DiscountCurve& discount_curve,
                                   const DiscountCurve& forward_curve,
                                   double expiry, double tenor, double pay_freq) {
    require(expiry > 0.0 && tenor > 0.0, "swaption: expiry and tenor must be > 0");
    require(pay_freq > 0.0, "swaption: payment frequency must be positive");
    double tau = 1.0 / pay_freq;
    double annuity = 0.0;
    double float_pv = 0.0;
    // Index payment dates by an explicit integer period count rather than
    // accumulating `t += tau`, which drifts for non-power-of-two tau (e.g.
    // monthly) over long tenors and makes the final-period inclusion depend on
    // a float tolerance (#10).
    int n = static_cast<int>(std::lround(tenor / tau));
    for (int i = 1; i <= n; ++i) {
        double t = expiry + i * tau;
        double df = discount_curve.discount(t);
        annuity += tau * df;
        float_pv += tau * df * forward_curve.forward_rate(t - tau, t);
    }
    SwaptionResult res;
    res.annuity = annuity;
    res.forward_swap_rate = float_pv / annuity;
    return res;
}

inline SwaptionResult swap_metrics(const DiscountCurve& curve, double expiry,
                                   double tenor, double pay_freq) {
    return swap_metrics(curve, curve, expiry, tenor, pay_freq);
}
}  // namespace detail

/// European swaption: option expiring at `expiry` (years) to enter a swap of
/// `tenor` years with fixed rate K, `pay_freq` payments per year, unit
/// notional. Volatility is Black (lognormal) or normal per vol_type.
/// Multi-curve form: discounting on `discount_curve` (OIS), forwards from
/// `forward_curve`.
inline SwaptionResult swaption_price(const DiscountCurve& discount_curve,
                                     const DiscountCurve& forward_curve,
                                     SwaptionType st, double K, double vol,
                                     double expiry, double tenor,
                                     double pay_freq = 2.0,
                                     RateVolType vol_type = RateVolType::Lognormal) {
    SwaptionResult res =
        detail::swap_metrics(discount_curve, forward_curve, expiry, tenor, pay_freq);
    OptionType type =
        (st == SwaptionType::Payer) ? OptionType::Call : OptionType::Put;
    double undiscounted;
    if (vol_type == RateVolType::Lognormal) {
        require(res.forward_swap_rate > 0.0,
                "swaption: lognormal vol requires positive forward swap rate");
        undiscounted = gbs_price(type, res.forward_swap_rate, K, expiry, 0.0, 0.0, vol);
    } else {
        undiscounted = bachelier_price(type, res.forward_swap_rate, K, expiry, 0.0, vol);
    }
    res.price = res.annuity * undiscounted;
    return res;
}

/// Single-curve convenience: projection and discounting on the same curve.
inline SwaptionResult swaption_price(const DiscountCurve& curve, SwaptionType st,
                                     double K, double vol, double expiry,
                                     double tenor, double pay_freq = 2.0,
                                     RateVolType vol_type = RateVolType::Lognormal) {
    return swaption_price(curve, curve, st, K, vol, expiry, tenor, pay_freq,
                          vol_type);
}

/// European swaption priced with a SABR vol smile on the swap rate.
inline SwaptionResult sabr_swaption_price(const DiscountCurve& curve,
                                          SwaptionType st, double K,
                                          const SabrParams& p, double expiry,
                                          double tenor, double pay_freq = 2.0) {
    SwaptionResult res = detail::swap_metrics(curve, expiry, tenor, pay_freq);
    double vol = sabr_lognormal_vol(res.forward_swap_rate, K, expiry, p);
    OptionType type =
        (st == SwaptionType::Payer) ? OptionType::Call : OptionType::Put;
    res.price = res.annuity *
                gbs_price(type, res.forward_swap_rate, K, expiry, 0.0, 0.0, vol);
    return res;
}

}  // namespace opal
