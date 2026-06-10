// Lookback option analytics (continuous monitoring).
//  - Floating strike: Goldman-Sosin-Gatto (1979).
//  - Fixed strike: Conze-Viswanathan (1991).
#pragma once

#include <cmath>

#include "opal/core/types.hpp"
#include "opal/math/normal.hpp"
#include "opal/models/black_scholes.hpp"

namespace opal {

/// Floating-strike lookback. Call pays S_T - S_min, put pays S_max - S_T.
/// `extremum` is the running min (call) or max (put) observed so far;
/// defaults to S for a newly issued option.
inline double floating_lookback_price(OptionType type, double S, double T,
                                      double r, double q, double sigma,
                                      double extremum = 0.0) {
    require(S > 0.0, "lookback: spot must be positive");
    require(sigma > 0.0 && T > 0.0, "lookback: sigma and T must be positive");
    double b = r - q;
    if (std::fabs(b) < 1e-8) b = 1e-8;  // formulas have removable b=0 singularity
    double M = (extremum > 0.0) ? extremum : S;
    double vsq = sigma * sigma;
    double sqrtT = std::sqrt(T);
    double a1 = (std::log(S / M) + (b + 0.5 * vsq) * T) / (sigma * sqrtT);
    double a2 = a1 - sigma * sqrtT;
    double carry = std::exp((b - r) * T);
    double df = std::exp(-r * T);

    if (type == OptionType::Call) {
        // M = running minimum
        return S * carry * math::norm_cdf(a1) - M * df * math::norm_cdf(a2) +
               S * df * (vsq / (2.0 * b)) *
                   (std::pow(S / M, -2.0 * b / vsq) *
                        math::norm_cdf(-a1 + 2.0 * b * sqrtT / sigma) -
                    std::exp(b * T) * math::norm_cdf(-a1));
    }
    // Put: M = running maximum
    double b1 = (std::log(S / M) + (b + 0.5 * vsq) * T) / (sigma * sqrtT);
    double b2 = b1 - sigma * sqrtT;
    return M * df * math::norm_cdf(-b2) - S * carry * math::norm_cdf(-b1) +
           S * df * (vsq / (2.0 * b)) *
               (-std::pow(S / M, -2.0 * b / vsq) *
                    math::norm_cdf(b1 - 2.0 * b * sqrtT / sigma) +
                std::exp(b * T) * math::norm_cdf(b1));
}

/// Fixed-strike lookback. Call pays max(S_max - K, 0), put pays max(K - S_min, 0).
/// `extremum` is the running max (call) / min (put) observed so far;
/// defaults to S for a newly issued option.
inline double fixed_lookback_price(OptionType type, double S, double K, double T,
                                   double r, double q, double sigma,
                                   double extremum = 0.0) {
    require(S > 0.0 && K > 0.0, "lookback: spot and strike must be positive");
    require(sigma > 0.0 && T > 0.0, "lookback: sigma and T must be positive");
    double b = r - q;
    if (std::fabs(b) < 1e-8) b = 1e-8;  // formulas have removable b=0 singularity
    double M = (extremum > 0.0) ? extremum : S;
    double vsq = sigma * sigma;
    double sqrtT = std::sqrt(T);
    double carry = std::exp((b - r) * T);
    double df = std::exp(-r * T);

    if (type == OptionType::Call) {
        if (K > M) {
            double d1 =
                (std::log(S / K) + (b + 0.5 * vsq) * T) / (sigma * sqrtT);
            double d2 = d1 - sigma * sqrtT;
            return S * carry * math::norm_cdf(d1) - K * df * math::norm_cdf(d2) +
                   S * df * (vsq / (2.0 * b)) *
                       (-std::pow(S / K, -2.0 * b / vsq) *
                            math::norm_cdf(d1 - 2.0 * b * sqrtT / sigma) +
                        std::exp(b * T) * math::norm_cdf(d1));
        }
        double e1 = (std::log(S / M) + (b + 0.5 * vsq) * T) / (sigma * sqrtT);
        double e2 = e1 - sigma * sqrtT;
        return df * (M - K) + S * carry * math::norm_cdf(e1) -
               M * df * math::norm_cdf(e2) +
               S * df * (vsq / (2.0 * b)) *
                   (-std::pow(S / M, -2.0 * b / vsq) *
                        math::norm_cdf(e1 - 2.0 * b * sqrtT / sigma) +
                    std::exp(b * T) * math::norm_cdf(e1));
    }

    // Put
    if (K < M) {
        double d1 = (std::log(S / K) + (b + 0.5 * vsq) * T) / (sigma * sqrtT);
        double d2 = d1 - sigma * sqrtT;
        return K * df * math::norm_cdf(-d2) - S * carry * math::norm_cdf(-d1) +
               S * df * (vsq / (2.0 * b)) *
                   (std::pow(S / K, -2.0 * b / vsq) *
                        math::norm_cdf(-d1 + 2.0 * b * sqrtT / sigma) -
                    std::exp(b * T) * math::norm_cdf(-d1));
    }
    double f1 = (std::log(S / M) + (b + 0.5 * vsq) * T) / (sigma * sqrtT);
    double f2 = f1 - sigma * sqrtT;
    return df * (K - M) - S * carry * math::norm_cdf(-f1) +
           M * df * math::norm_cdf(-f2) +
           S * df * (vsq / (2.0 * b)) *
               (std::pow(S / M, -2.0 * b / vsq) *
                    math::norm_cdf(-f1 + 2.0 * b * sqrtT / sigma) -
                std::exp(b * T) * math::norm_cdf(-f1));
}

}  // namespace opal
