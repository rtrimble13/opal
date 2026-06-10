// Digital (binary) option analytics under generalized Black-Scholes.
#pragma once

#include <cmath>

#include "opal/core/types.hpp"
#include "opal/math/normal.hpp"
#include "opal/models/black_scholes.hpp"

namespace opal {

/// Cash-or-nothing: pays `cash` at expiry if the option finishes in the money.
inline double cash_or_nothing_price(OptionType type, double S, double K, double T,
                                    double r, double q, double sigma,
                                    double cash = 1.0) {
    require(S > 0.0 && K > 0.0, "digital: spot and strike must be positive");
    double phi = type_sign(type);
    double b = r - q;
    if (T == 0.0 || sigma == 0.0) {
        double fwd = S * std::exp(b * T);
        return std::exp(-r * T) * cash * (phi * (fwd - K) > 0.0 ? 1.0 : 0.0);
    }
    auto [d1, d2, sqrtT] = detail::d12(S, K, T, b, sigma);
    (void)d1;
    (void)sqrtT;
    return cash * std::exp(-r * T) * math::norm_cdf(phi * d2);
}

/// Asset-or-nothing: pays S_T at expiry if the option finishes in the money.
inline double asset_or_nothing_price(OptionType type, double S, double K, double T,
                                     double r, double q, double sigma) {
    require(S > 0.0 && K > 0.0, "digital: spot and strike must be positive");
    double phi = type_sign(type);
    double b = r - q;
    if (T == 0.0 || sigma == 0.0) {
        double fwd = S * std::exp(b * T);
        return std::exp(-r * T) * fwd * (phi * (fwd - K) > 0.0 ? 1.0 : 0.0);
    }
    auto [d1, d2, sqrtT] = detail::d12(S, K, T, b, sigma);
    (void)d2;
    (void)sqrtT;
    return S * std::exp((b - r) * T) * math::norm_cdf(phi * d1);
}

/// Gap option: pays (S_T - payoff_strike) if S_T crosses trigger_strike.
inline double gap_option_price(OptionType type, double S, double trigger_strike,
                               double payoff_strike, double T, double r, double q,
                               double sigma) {
    return asset_or_nothing_price(type, S, trigger_strike, T, r, q, sigma) *
               type_sign(type) -
           type_sign(type) * payoff_strike *
               cash_or_nothing_price(type, S, trigger_strike, T, r, q, sigma);
}

}  // namespace opal
