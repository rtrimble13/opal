// Implied volatility solvers (Newton with safeguarded Brent fallback).
#pragma once

#include <algorithm>
#include <cmath>

#include "opal/core/types.hpp"
#include "opal/math/normal.hpp"
#include "opal/math/solvers.hpp"
#include "opal/models/black_scholes.hpp"

namespace opal {

/// Implied Black-Scholes-Merton volatility from a price.
inline double implied_vol_bsm(OptionType type, double price, double S, double K,
                              double T, double r, double q,
                              double tol = 1e-10) {
    require(price > 0.0, "implied vol: price must be positive");
    require(T > 0.0, "implied vol: time to expiry must be positive");
    double b = r - q;
    double phi = type_sign(type);
    double intrinsic =
        std::exp(-r * T) * std::max(phi * (S * std::exp(b * T) - K), 0.0);
    require(price >= intrinsic - 1e-12,
            "implied vol: price below discounted intrinsic value");
    double upper_bound = (type == OptionType::Call)
                             ? S * std::exp((b - r) * T)
                             : K * std::exp(-r * T);
    require(price <= upper_bound + 1e-12,
            "implied vol: price above no-arbitrage bound");

    auto f = [&](double v) { return gbs_price(type, S, K, T, r, b, v) - price; };
    auto df = [&](double v) { return gbs_greeks(type, S, K, T, r, b, v).vega; };

    // Corrado-Miller style starting guess, clamped.
    double guess = std::sqrt(2.0 * math::PI / T) * price / S;
    if (guess < 0.05) guess = 0.05;
    if (guess > 2.0) guess = 2.0;
    return math::newton_safe(f, df, guess, 1e-9, 20.0, tol);
}

/// Implied Black-76 volatility (forward-based).
inline double implied_vol_black76(OptionType type, double price, double F,
                                  double K, double T, double r,
                                  double tol = 1e-10) {
    return implied_vol_bsm(type, price, F, K, T, r, r, tol);
}

/// Implied Bachelier (normal) volatility.
inline double implied_vol_bachelier(OptionType type, double price, double F,
                                    double K, double T, double r,
                                    double tol = 1e-10) {
    require(price > 0.0, "implied vol: price must be positive");
    require(T > 0.0, "implied vol: time to expiry must be positive");
    double phi = type_sign(type);
    double intrinsic = std::exp(-r * T) * std::max(phi * (F - K), 0.0);
    require(price >= intrinsic - 1e-12,
            "implied vol: price below discounted intrinsic value");
    auto f = [&](double v) { return bachelier_price(type, F, K, T, r, v) - price; };
    // Bracket: normal vol rarely exceeds ~10 * |F| in any sane market.
    double hi = 10.0 * (std::fabs(F) + std::fabs(K) + 1.0);
    return math::brent(f, 1e-12, hi, tol);
}

}  // namespace opal
