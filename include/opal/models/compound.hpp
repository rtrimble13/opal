// Compound options (option-on-option) via the Geske (1979) / Rubinstein
// closed form, built on the bivariate normal CDF (math::bivar_norm_cdf).
//
// A compound option is an `outer` option, expiring at t1 with strike K1, on an
// `inner` vanilla option, expiring at T2 >= t1 with strike K2, on an asset S.
// All four combinations are supported (call/put on call/put). Cost of carry
// b = r - q, matching gbs_price.
#pragma once

#include <algorithm>
#include <cmath>

#include "opal/core/types.hpp"
#include "opal/math/normal.hpp"
#include "opal/math/solvers.hpp"
#include "opal/models/black_scholes.hpp"

namespace opal {

/// Price of a compound option. `outer` is the right held today (expiring at
/// t1, strike K1); `inner` is the option delivered on exercise (expiring at
/// T2 > t1, strike K2). Returns the present value.
inline double compound_option_price(OptionType outer, OptionType inner, double S,
                                    double K1, double K2, double t1, double T2,
                                    double r, double q, double sig) {
    require(S > 0.0 && K1 > 0.0 && K2 > 0.0,
            "compound: spot and strikes must be positive");
    require(sig > 0.0, "compound: volatility must be positive");
    require(t1 > 0.0 && T2 > t1,
            "compound: require 0 < t1 < T2 (outer expiry before inner)");

    double b = r - q;
    double tau = T2 - t1;  // inner option's remaining life at outer expiry

    // Inner option value at t1 as a function of the asset price there.
    auto inner_value = [&](double x) {
        return gbs_price(inner, x, K2, tau, r, b, sig);
    };

    // Critical asset price I at t1 where the inner option is worth exactly K1
    // (the outer exercise boundary). For an inner call the value is increasing
    // in x (from 0 to +inf); for an inner put it is decreasing (from
    // K2 e^{-r tau} down to 0). If K1 lies outside that range the outer option
    // is never exercised and the compound is worthless.
    double i_star;
    {
        auto f = [&](double x) { return inner_value(x) - K1; };
        double lo = 1e-8 * S, hi = S;
        if (inner == OptionType::Call) {
            while (f(hi) < 0.0 && hi < 1e12 * S) hi *= 2.0;  // grow to bracket
            if (f(hi) < 0.0) return 0.0;                     // never exercised
        } else {
            // Inner put: max value at x -> 0 is K2 e^{-r tau}.
            if (K2 * std::exp(-r * tau) <= K1) return 0.0;   // never exercised
            while (f(lo) < 0.0 && lo > 1e-300) lo *= 0.5;
        }
        i_star = math::brent(f, lo, hi);
    }

    double sqrt1 = std::sqrt(t1), sqrt2 = std::sqrt(T2);
    double rho = sqrt1 / sqrt2;
    double y1 = (std::log(S / i_star) + (b + 0.5 * sig * sig) * t1) / (sig * sqrt1);
    double y2 = y1 - sig * sqrt1;
    double z1 = (std::log(S / K2) + (b + 0.5 * sig * sig) * T2) / (sig * sqrt2);
    double z2 = z1 - sig * sqrt2;

    double df = std::exp((b - r) * T2);  // S e^{(b-r)T2} = S e^{-q T2}
    double dr2 = std::exp(-r * T2), dr1 = std::exp(-r * t1);
    using math::bivar_norm_cdf;
    using math::norm_cdf;

    if (outer == OptionType::Call && inner == OptionType::Call)
        return S * df * bivar_norm_cdf(z1, y1, rho) -
               K2 * dr2 * bivar_norm_cdf(z2, y2, rho) - K1 * dr1 * norm_cdf(y2);
    if (outer == OptionType::Put && inner == OptionType::Call)
        return K2 * dr2 * bivar_norm_cdf(z2, -y2, -rho) -
               S * df * bivar_norm_cdf(z1, -y1, -rho) + K1 * dr1 * norm_cdf(-y2);
    if (outer == OptionType::Call && inner == OptionType::Put)
        return K2 * dr2 * bivar_norm_cdf(-z2, -y2, rho) -
               S * df * bivar_norm_cdf(-z1, -y1, rho) - K1 * dr1 * norm_cdf(-y2);
    // Put on put.
    return S * df * bivar_norm_cdf(-z1, y1, -rho) -
           K2 * dr2 * bivar_norm_cdf(-z2, y2, -rho) + K1 * dr1 * norm_cdf(y2);
}

}  // namespace opal
