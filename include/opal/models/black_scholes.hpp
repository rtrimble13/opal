// Generalized Black-Scholes-Merton analytic pricing and Greeks.
//
// The generalized model prices on (S, K, T, r, b, sigma) where b is the cost
// of carry:  b = r - q   for an equity with continuous dividend yield q,
//            b = 0       for an option on a futures/forward (Black-76),
//            b = r - rf  for an FX option (rf = foreign rate).
#pragma once

#include <cmath>

#include "opal/core/types.hpp"
#include "opal/math/normal.hpp"

namespace opal {

namespace detail {
struct D12 {
    double d1, d2, sqrtT;
};
inline D12 d12(double S, double K, double T, double b, double sigma) {
    double sqrtT = std::sqrt(T);
    double d1 = (std::log(S / K) + (b + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    return {d1, d1 - sigma * sqrtT, sqrtT};
}
}  // namespace detail

/// Generalized Black-Scholes price (cost-of-carry form).
inline double gbs_price(OptionType type, double S, double K, double T, double r,
                        double b, double sigma) {
    require(S > 0.0 && K > 0.0, "gbs: spot and strike must be positive");
    require(sigma >= 0.0, "gbs: volatility must be non-negative");
    require(T >= 0.0, "gbs: time to expiry must be non-negative");
    double phi = type_sign(type);
    if (T == 0.0 || sigma == 0.0) {
        // Deterministic forward; price is discounted intrinsic on the forward.
        double fwd = S * std::exp(b * T);
        return std::exp(-r * T) * std::max(phi * (fwd - K), 0.0);
    }
    auto [d1, d2, sqrtT] = detail::d12(S, K, T, b, sigma);
    return phi * (S * std::exp((b - r) * T) * math::norm_cdf(phi * d1) -
                  K * std::exp(-r * T) * math::norm_cdf(phi * d2));
}

/// Full analytic Greeks for the generalized Black-Scholes model.
/// `rho` is the sensitivity to r holding the carry spread (r - b) fixed,
/// i.e. the standard equity rho when b = r - q.
inline Greeks gbs_greeks(OptionType type, double S, double K, double T, double r,
                         double b, double sigma) {
    Greeks g;
    g.price = gbs_price(type, S, K, T, r, b, sigma);
    if (T <= 0.0 || sigma <= 0.0) {
        double fwd = S * std::exp(b * T);
        double phi = type_sign(type);
        g.delta = (phi * (fwd - K) > 0.0) ? phi * std::exp((b - r) * T) : 0.0;
        return g;
    }
    auto [d1, d2, sqrtT] = detail::d12(S, K, T, b, sigma);
    double phi = type_sign(type);
    double df_carry = std::exp((b - r) * T);  // e^{-qT} in equity terms
    double df = std::exp(-r * T);
    double nd1 = math::norm_pdf(d1);
    double Nd1 = math::norm_cdf(phi * d1);
    double Nd2 = math::norm_cdf(phi * d2);

    g.delta = phi * df_carry * Nd1;
    g.gamma = df_carry * nd1 / (S * sigma * sqrtT);
    g.vega = S * df_carry * nd1 * sqrtT;
    g.theta = -S * df_carry * nd1 * sigma / (2.0 * sqrtT) -
              phi * (b - r) * S * df_carry * Nd1 - phi * r * K * df * Nd2;
    g.rho = phi * K * T * df * Nd2;
    g.vanna = -df_carry * nd1 * d2 / sigma;
    g.volga = S * df_carry * nd1 * sqrtT * d1 * d2 / sigma;
    g.charm = -phi * (b - r) * df_carry * Nd1 -
              df_carry * nd1 * (2.0 * b * T - d2 * sigma * sqrtT) /
                  (2.0 * T * sigma * sqrtT);
    return g;
}

/// Black-Scholes-Merton on an equity paying continuous dividend yield q.
inline double bsm_price(OptionType type, double S, double K, double T, double r,
                        double q, double sigma) {
    return gbs_price(type, S, K, T, r, r - q, sigma);
}

inline Greeks bsm_greeks(OptionType type, double S, double K, double T, double r,
                         double q, double sigma) {
    return gbs_greeks(type, S, K, T, r, r - q, sigma);
}

/// Black-76: option on a forward/futures price F, discounted at r.
inline double black76_price(OptionType type, double F, double K, double T,
                            double r, double sigma) {
    return gbs_price(type, F, K, T, r, 0.0, sigma);
}

inline Greeks black76_greeks(OptionType type, double F, double K, double T,
                             double r, double sigma) {
    return gbs_greeks(type, F, K, T, r, 0.0, sigma);
}

/// Bachelier (normal) model: option on forward F with absolute volatility
/// sigma_n (same units as F). Standard for rates markets in low/negative
/// rate environments.
inline double bachelier_price(OptionType type, double F, double K, double T,
                              double r, double sigma_n) {
    require(sigma_n >= 0.0, "bachelier: volatility must be non-negative");
    require(T >= 0.0, "bachelier: time to expiry must be non-negative");
    double phi = type_sign(type);
    double df = std::exp(-r * T);
    if (T == 0.0 || sigma_n == 0.0) return df * std::max(phi * (F - K), 0.0);
    double sd = sigma_n * std::sqrt(T);
    double d = (F - K) / sd;
    return df * (phi * (F - K) * math::norm_cdf(phi * d) + sd * math::norm_pdf(d));
}

inline Greeks bachelier_greeks(OptionType type, double F, double K, double T,
                               double r, double sigma_n) {
    Greeks g;
    g.price = bachelier_price(type, F, K, T, r, sigma_n);
    if (T <= 0.0 || sigma_n <= 0.0) return g;
    double phi = type_sign(type);
    double df = std::exp(-r * T);
    double sqrtT = std::sqrt(T);
    double sd = sigma_n * sqrtT;
    double d = (F - K) / sd;
    g.delta = phi * df * math::norm_cdf(phi * d);
    g.gamma = df * math::norm_pdf(d) / sd;
    g.vega = df * math::norm_pdf(d) * sqrtT;
    g.theta = -df * math::norm_pdf(d) * sigma_n / (2.0 * sqrtT) + r * g.price;
    g.rho = -T * g.price;
    return g;
}

}  // namespace opal
