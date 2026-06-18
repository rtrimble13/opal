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
    // (the outer exercise boundary). The inner value is monotone in x
    // (increasing 0->inf for a call, decreasing K2 e^{-r tau}->0 for a put), so
    // bracket the root generically by expanding until the endpoints straddle
    // it. The only unbracketable case is "inner value < K1 for every x" (the
    // inner is never worth the outer strike): then the outer call is never
    // exercised (worthless) while the outer put is *always* exercised -- the
    // holder sells the inner for K1 -- worth K1 e^{-r t1} minus today's inner
    // value. (Both limits keep outer put-call parity exact.)
    auto f = [&](double x) { return inner_value(x) - K1; };
    double lo = 1e-8 * S, hi = S;
    double f_lo = f(lo), f_hi = f(hi);
    while (std::signbit(f_hi) == std::signbit(f_lo) && hi < 1e12 * S) {
        hi *= 2.0;
        f_hi = f(hi);
    }
    if (std::signbit(f_hi) == std::signbit(f_lo)) {  // no boundary in range
        if (outer == OptionType::Call) return 0.0;
        return K1 * std::exp(-r * t1) - gbs_price(inner, S, K2, T2, r, b, sig);
    }
    double i_star = math::brent(f, lo, hi);

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
