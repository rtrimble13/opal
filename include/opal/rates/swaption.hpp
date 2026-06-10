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
inline SwaptionResult swap_metrics(const DiscountCurve& curve, double expiry,
                                   double tenor, double pay_freq) {
    require(expiry > 0.0 && tenor > 0.0, "swaption: expiry and tenor must be > 0");
    require(pay_freq > 0.0, "swaption: payment frequency must be positive");
    double tau = 1.0 / pay_freq;
    double annuity = 0.0;
    for (double t = expiry + tau; t <= expiry + tenor + 1e-10; t += tau)
        annuity += tau * curve.discount(t);
    double swap_rate =
        (curve.discount(expiry) - curve.discount(expiry + tenor)) / annuity;
    SwaptionResult res;
    res.annuity = annuity;
    res.forward_swap_rate = swap_rate;
    return res;
}
}  // namespace detail

/// European swaption: option expiring at `expiry` (years) to enter a swap of
/// `tenor` years with fixed rate K, `pay_freq` payments per year, unit
/// notional. Volatility is Black (lognormal) or normal per vol_type.
inline SwaptionResult swaption_price(const DiscountCurve& curve, SwaptionType st,
                                     double K, double vol, double expiry,
                                     double tenor, double pay_freq = 2.0,
                                     RateVolType vol_type = RateVolType::Lognormal) {
    SwaptionResult res = detail::swap_metrics(curve, expiry, tenor, pay_freq);
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
