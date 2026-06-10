// Heston (1993) stochastic volatility model, semi-analytic pricing via the
// characteristic function with the Albrecher et al. "little Heston trap"
// formulation (numerically stable for long maturities).
//
//   dS = (r - q) S dt + sqrt(v) S dW1
//   dv = kappa (theta - v) dt + xi sqrt(v) dW2,   d<W1,W2> = rho dt
#pragma once

#include <cmath>
#include <complex>

#include "opal/core/types.hpp"
#include "opal/math/solvers.hpp"

namespace opal {

struct HestonParams {
    double v0;     // initial variance
    double kappa;  // mean reversion speed
    double theta;  // long-run variance
    double xi;     // vol of vol
    double rho;    // spot/vol correlation

    bool feller_satisfied() const { return 2.0 * kappa * theta >= xi * xi; }
};

namespace detail {

inline std::complex<double> heston_cf(double phi_re, int j, double S, double T,
                                      double r, double q, const HestonParams& p) {
    using cd = std::complex<double>;
    const cd i(0.0, 1.0);
    cd phi(phi_re, 0.0);

    double u = (j == 1) ? 0.5 : -0.5;
    double bj = (j == 1) ? p.kappa - p.rho * p.xi : p.kappa;
    double a = p.kappa * p.theta;
    double xi2 = p.xi * p.xi;

    cd beta = bj - p.rho * p.xi * i * phi;
    cd d = std::sqrt(beta * beta - xi2 * (2.0 * u * i * phi - phi * phi));
    cd g = (beta - d) / (beta + d);  // little-trap form (uses -d)

    cd exp_dT = std::exp(-d * T);
    cd C = (r - q) * i * phi * T +
           (a / xi2) * ((beta - d) * T -
                        2.0 * std::log((1.0 - g * exp_dT) / (1.0 - g)));
    cd D = ((beta - d) / xi2) * (1.0 - exp_dT) / (1.0 - g * exp_dT);

    return std::exp(C + D * p.v0 + i * phi * std::log(S));
}

inline double heston_prob(int j, double S, double K, double T, double r, double q,
                          const HestonParams& p) {
    double lnK = std::log(K);
    auto integrand = [&](double phi) -> double {
        if (phi < 1e-10) phi = 1e-10;
        std::complex<double> i(0.0, 1.0);
        std::complex<double> val =
            std::exp(-i * phi * lnK) * heston_cf(phi, j, S, T, r, q, p) / (i * phi);
        return val.real();
    };
    // The integrand decays like exp(-c*phi); 200 is far past machine epsilon
    // for typical parameters.
    double integral = math::integrate(integrand, 1e-8, 200.0, 1e-10);
    return 0.5 + integral / M_PI;
}

}  // namespace detail

/// European option under Heston via the semi-analytic formula.
inline double heston_price(OptionType type, double S, double K, double T,
                           double r, double q, const HestonParams& p) {
    require(S > 0.0 && K > 0.0, "heston: spot and strike must be positive");
    require(p.v0 >= 0.0 && p.theta >= 0.0 && p.kappa > 0.0 && p.xi > 0.0,
            "heston: invalid parameters");
    require(p.rho >= -1.0 && p.rho <= 1.0, "heston: rho must be in [-1, 1]");
    require(T > 0.0, "heston: time to expiry must be positive");

    double P1 = detail::heston_prob(1, S, K, T, r, q, p);
    double P2 = detail::heston_prob(2, S, K, T, r, q, p);
    double call =
        S * std::exp(-q * T) * P1 - K * std::exp(-r * T) * P2;
    if (type == OptionType::Call) return std::max(call, 0.0);
    // put-call parity
    double put = call - S * std::exp(-q * T) + K * std::exp(-r * T);
    return std::max(put, 0.0);
}

}  // namespace opal
