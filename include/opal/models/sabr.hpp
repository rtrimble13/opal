// SABR stochastic volatility model: Hagan et al. (2002) asymptotic implied
// volatility, the market standard for interest rate vol cubes.
//
//   dF = alpha F^beta dW1
//   dalpha = nu alpha dW2,   d<W1,W2> = rho dt
#pragma once

#include <cmath>

#include "opal/core/types.hpp"
#include "opal/models/black_scholes.hpp"

namespace opal {

struct SabrParams {
    double alpha;  // initial vol level
    double beta;   // CEV exponent in [0, 1]
    double rho;    // correlation in [-1, 1]
    double nu;     // vol of vol
};

/// Hagan (2002) lognormal (Black) implied volatility under SABR.
inline double sabr_lognormal_vol(double F, double K, double T,
                                 const SabrParams& p) {
    require(F > 0.0 && K > 0.0, "sabr: forward and strike must be positive");
    require(p.beta >= 0.0 && p.beta <= 1.0, "sabr: beta must be in [0, 1]");
    require(p.rho > -1.0 && p.rho < 1.0, "sabr: rho must be in (-1, 1)");
    require(p.alpha > 0.0, "sabr: alpha must be positive");

    double one_m_beta = 1.0 - p.beta;
    double logFK = std::log(F / K);
    double FK_pow = std::pow(F * K, 0.5 * one_m_beta);

    // Common correction term (order-T expansion).
    double term1 = one_m_beta * one_m_beta * p.alpha * p.alpha /
                   (24.0 * FK_pow * FK_pow);
    double term2 = 0.25 * p.rho * p.beta * p.nu * p.alpha / FK_pow;
    double term3 = (2.0 - 3.0 * p.rho * p.rho) * p.nu * p.nu / 24.0;
    double correction = 1.0 + (term1 + term2 + term3) * T;

    if (std::fabs(logFK) < 1e-12) {
        // ATM formula
        return p.alpha / std::pow(F, one_m_beta) * correction;
    }

    double denom = FK_pow * (1.0 + one_m_beta * one_m_beta * logFK * logFK / 24.0 +
                             std::pow(one_m_beta, 4) * std::pow(logFK, 4) / 1920.0);
    double z = (p.nu / p.alpha) * FK_pow * logFK;
    double zx;
    if (std::fabs(z) < 1e-8) {
        zx = 1.0 - 0.5 * p.rho * z;  // series expansion of z/x(z)
    } else {
        double x = std::log((std::sqrt(1.0 - 2.0 * p.rho * z + z * z) + z - p.rho) /
                            (1.0 - p.rho));
        zx = z / x;
    }
    return (p.alpha / denom) * zx * correction;
}

/// Hagan (2002) normal (Bachelier) implied volatility under SABR.
inline double sabr_normal_vol(double F, double K, double T, const SabrParams& p) {
    require(p.beta >= 0.0 && p.beta < 1.0 + 1e-12,
            "sabr: beta must be in [0, 1]");
    // For beta < 1 use the lognormal vol converted via the standard
    // moment-matching approximation; for beta == 0 SABR is exactly normal.
    if (p.beta == 0.0) {
        double zeta = (p.nu / p.alpha) * (F - K);
        double zx;
        if (std::fabs(zeta) < 1e-8) {
            zx = 1.0 - 0.5 * p.rho * zeta;
        } else {
            double x = std::log(
                (std::sqrt(1.0 - 2.0 * p.rho * zeta + zeta * zeta) + zeta - p.rho) /
                (1.0 - p.rho));
            zx = zeta / x;
        }
        double correction =
            1.0 + (2.0 - 3.0 * p.rho * p.rho) * p.nu * p.nu / 24.0 * T;
        return p.alpha * zx * correction;
    }
    // General beta: convert Black vol to normal vol via ATM-equivalent mapping.
    double black = sabr_lognormal_vol(F, K, T, p);
    double Fmid = std::sqrt(F * K);
    return black * Fmid * (1.0 - black * black * T / 24.0);
}

/// Convenience: price a European option on forward F using SABR-implied
/// Black vol in the Black-76 formula.
inline double sabr_price(OptionType type, double F, double K, double T, double r,
                         const SabrParams& p) {
    double vol = sabr_lognormal_vol(F, K, T, p);
    return black76_price(type, F, K, T, r, vol);
}

}  // namespace opal
