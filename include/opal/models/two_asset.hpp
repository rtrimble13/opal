// Two-asset closed-form analytics built on the bivariate normal CDF
// (math::bivar_norm_cdf):
//   * exchange option (Margrabe 1978)
//   * options on the maximum / minimum of two assets (Stulz 1982)
//   * two-asset correlation options
//
// Conventions match the rest of the library: each asset has spot S_i,
// continuous dividend yield q_i and lognormal vol sig_i; rho is the
// instantaneous correlation of the two assets' returns; r is the risk-free
// rate. Internally the cost of carry is b_i = r - q_i.
#pragma once

#include <algorithm>
#include <cmath>

#include "opal/core/types.hpp"
#include "opal/math/normal.hpp"

namespace opal {

/// Margrabe (1978) option to exchange asset 2 for asset 1; pays
/// max(S1 - S2, 0) at T. Independent of the risk-free rate (the discount
/// cancels), so no rate argument is taken.
inline double exchange_option_price(double S1, double S2, double T, double q1,
                                    double q2, double sig1, double sig2,
                                    double rho) {
    require(S1 > 0.0 && S2 > 0.0, "exchange: spots must be positive");
    require(sig1 >= 0.0 && sig2 >= 0.0, "exchange: vols must be non-negative");
    require(rho >= -1.0 && rho <= 1.0, "exchange: rho must be in [-1, 1]");
    require(T >= 0.0, "exchange: time to expiry must be non-negative");
    double f1 = S1 * std::exp(-q1 * T);  // PV of receiving asset 1 at T
    double f2 = S2 * std::exp(-q2 * T);
    double sig = std::sqrt(sig1 * sig1 + sig2 * sig2 - 2.0 * rho * sig1 * sig2);
    if (T == 0.0 || sig == 0.0) return std::max(f1 - f2, 0.0);
    double sd = sig * std::sqrt(T);
    double d1 = (std::log(f1 / f2) + 0.5 * sig * sig * T) / sd;
    double d2 = d1 - sd;
    return f1 * math::norm_cdf(d1) - f2 * math::norm_cdf(d2);
}

namespace detail {

// Shared Stulz quantities for the call formulas.
struct StulzCall {
    double cmax;
    double cmin;
};

inline StulzCall stulz_calls(double S1, double S2, double K, double T, double r,
                             double q1, double q2, double sig1, double sig2,
                             double rho) {
    double b1 = r - q1, b2 = r - q2;
    double sqrtT = std::sqrt(T);
    double sig = std::sqrt(sig1 * sig1 + sig2 * sig2 - 2.0 * rho * sig1 * sig2);
    double rho1 = (sig1 - rho * sig2) / sig;
    double rho2 = (sig2 - rho * sig1) / sig;
    double d = (std::log(S1 / S2) + (b1 - b2 + 0.5 * sig * sig) * T) / (sig * sqrtT);
    double y1 = (std::log(S1 / K) + (b1 + 0.5 * sig1 * sig1) * T) / (sig1 * sqrtT);
    double y2 = (std::log(S2 / K) + (b2 + 0.5 * sig2 * sig2) * T) / (sig2 * sqrtT);
    using math::bivar_norm_cdf;
    double df1 = std::exp((b1 - r) * T), df2 = std::exp((b2 - r) * T);
    double dfr = std::exp(-r * T);

    double cmax = S1 * df1 * bivar_norm_cdf(y1, d, rho1) +
                  S2 * df2 * bivar_norm_cdf(y2, -d + sig * sqrtT, rho2) -
                  K * dfr * (1.0 - bivar_norm_cdf(-y1 + sig1 * sqrtT,
                                                  -y2 + sig2 * sqrtT, rho));
    double cmin = S1 * df1 * bivar_norm_cdf(y1, -d, -rho1) +
                  S2 * df2 * bivar_norm_cdf(y2, d - sig * sqrtT, -rho2) -
                  K * dfr * bivar_norm_cdf(y1 - sig1 * sqrtT, y2 - sig2 * sqrtT, rho);
    return {cmax, cmin};
}

}  // namespace detail

/// Stulz (1982) option on the maximum of two assets, strike K.
inline double option_on_max_price(OptionType type, double S1, double S2, double K,
                                  double T, double r, double q1, double q2,
                                  double sig1, double sig2, double rho) {
    require(S1 > 0.0 && S2 > 0.0 && K > 0.0,
            "option_on_max: spots and strike must be positive");
    require(sig1 > 0.0 && sig2 > 0.0, "option_on_max: vols must be positive");
    require(rho >= -1.0 && rho <= 1.0, "option_on_max: rho must be in [-1, 1]");
    require(T > 0.0, "option_on_max: time to expiry must be positive");
    double cmax = detail::stulz_calls(S1, S2, K, T, r, q1, q2, sig1, sig2, rho).cmax;
    if (type == OptionType::Call) return std::max(cmax, 0.0);
    // Put via parity: Pmax = Cmax - PV[max(S1,S2)] + K e^{-rT}.
    double exch = exchange_option_price(S1, S2, T, q1, q2, sig1, sig2, rho);
    double max_fwd = S2 * std::exp(-q2 * T) + exch;  // e^{-rT} E[max(S1,S2)]
    return std::max(cmax - max_fwd + K * std::exp(-r * T), 0.0);
}

/// Stulz (1982) option on the minimum of two assets, strike K.
inline double option_on_min_price(OptionType type, double S1, double S2, double K,
                                  double T, double r, double q1, double q2,
                                  double sig1, double sig2, double rho) {
    require(S1 > 0.0 && S2 > 0.0 && K > 0.0,
            "option_on_min: spots and strike must be positive");
    require(sig1 > 0.0 && sig2 > 0.0, "option_on_min: vols must be positive");
    require(rho >= -1.0 && rho <= 1.0, "option_on_min: rho must be in [-1, 1]");
    require(T > 0.0, "option_on_min: time to expiry must be positive");
    double cmin = detail::stulz_calls(S1, S2, K, T, r, q1, q2, sig1, sig2, rho).cmin;
    if (type == OptionType::Call) return std::max(cmin, 0.0);
    // Put via parity: Pmin = Cmin - PV[min(S1,S2)] + K e^{-rT}.
    double exch = exchange_option_price(S1, S2, T, q1, q2, sig1, sig2, rho);
    double min_fwd = S1 * std::exp(-q1 * T) - exch;  // e^{-rT} E[min(S1,S2)]
    return std::max(cmin - min_fwd + K * std::exp(-r * T), 0.0);
}

/// Two-asset correlation option: a call pays max(S2 - K2, 0) only if S1 > K1
/// at expiry; a put pays max(K2 - S2, 0) only if S1 < K1.
inline double two_asset_correlation_price(OptionType type, double S1, double S2,
                                          double K1, double K2, double T, double r,
                                          double q1, double q2, double sig1,
                                          double sig2, double rho) {
    require(S1 > 0.0 && S2 > 0.0 && K1 > 0.0 && K2 > 0.0,
            "two_asset_correlation: spots and strikes must be positive");
    require(sig1 > 0.0 && sig2 > 0.0,
            "two_asset_correlation: vols must be positive");
    require(rho >= -1.0 && rho <= 1.0,
            "two_asset_correlation: rho must be in [-1, 1]");
    require(T > 0.0, "two_asset_correlation: time to expiry must be positive");
    double b1 = r - q1, b2 = r - q2;
    double sqrtT = std::sqrt(T);
    double y1 = (std::log(S1 / K1) + (b1 - 0.5 * sig1 * sig1) * T) / (sig1 * sqrtT);
    double y2 = (std::log(S2 / K2) + (b2 - 0.5 * sig2 * sig2) * T) / (sig2 * sqrtT);
    double df2 = std::exp((b2 - r) * T), dfr = std::exp(-r * T);
    using math::bivar_norm_cdf;
    if (type == OptionType::Call)
        return S2 * df2 *
                   bivar_norm_cdf(y2 + sig2 * sqrtT, y1 + rho * sig2 * sqrtT, rho) -
               K2 * dfr * bivar_norm_cdf(y2, y1, rho);
    return K2 * dfr * bivar_norm_cdf(-y2, -y1, rho) -
           S2 * df2 *
               bivar_norm_cdf(-y2 - sig2 * sqrtT, -y1 - rho * sig2 * sqrtT, rho);
}

}  // namespace opal
