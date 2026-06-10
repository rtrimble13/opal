// Single-barrier option analytics: Reiner-Rubinstein (1991) closed forms
// for continuously monitored barriers, as standardized in Haug (2007).
#pragma once

#include <cmath>

#include "opal/core/types.hpp"
#include "opal/math/normal.hpp"
#include "opal/models/black_scholes.hpp"

namespace opal {

namespace detail {

struct BarrierTerms {
    double S, K, H, T, r, b, sigma, rebate;
    double phi;  // +1 call, -1 put
    double eta;  // +1 down barrier, -1 up barrier
    double mu, lambda, sqrtT, vsq;

    BarrierTerms(OptionType type, BarrierType bt, double S_, double K_, double H_,
                 double T_, double r_, double b_, double sigma_, double rebate_)
        : S(S_), K(K_), H(H_), T(T_), r(r_), b(b_), sigma(sigma_), rebate(rebate_) {
        phi = type_sign(type);
        eta = is_down_barrier(bt) ? 1.0 : -1.0;
        vsq = sigma * sigma;
        mu = (b - 0.5 * vsq) / vsq;
        lambda = std::sqrt(mu * mu + 2.0 * r / vsq);
        sqrtT = std::sqrt(T);
    }

    double N(double x) const { return math::norm_cdf(x); }

    double A() const {
        double x1 = std::log(S / K) / (sigma * sqrtT) + (1.0 + mu) * sigma * sqrtT;
        return phi * S * std::exp((b - r) * T) * N(phi * x1) -
               phi * K * std::exp(-r * T) * N(phi * (x1 - sigma * sqrtT));
    }
    double B() const {
        double x2 = std::log(S / H) / (sigma * sqrtT) + (1.0 + mu) * sigma * sqrtT;
        return phi * S * std::exp((b - r) * T) * N(phi * x2) -
               phi * K * std::exp(-r * T) * N(phi * (x2 - sigma * sqrtT));
    }
    double C() const {
        double y1 = std::log(H * H / (S * K)) / (sigma * sqrtT) +
                    (1.0 + mu) * sigma * sqrtT;
        return phi * S * std::exp((b - r) * T) * std::pow(H / S, 2.0 * (mu + 1.0)) *
                   N(eta * y1) -
               phi * K * std::exp(-r * T) * std::pow(H / S, 2.0 * mu) *
                   N(eta * (y1 - sigma * sqrtT));
    }
    double D() const {
        double y2 = std::log(H / S) / (sigma * sqrtT) + (1.0 + mu) * sigma * sqrtT;
        return phi * S * std::exp((b - r) * T) * std::pow(H / S, 2.0 * (mu + 1.0)) *
                   N(eta * y2) -
               phi * K * std::exp(-r * T) * std::pow(H / S, 2.0 * mu) *
                   N(eta * (y2 - sigma * sqrtT));
    }
    // Rebate paid at expiry if a knock-in option never knocks in.
    double E() const {
        if (rebate == 0.0) return 0.0;
        double x2 = std::log(S / H) / (sigma * sqrtT) + (1.0 + mu) * sigma * sqrtT;
        double y2 = std::log(H / S) / (sigma * sqrtT) + (1.0 + mu) * sigma * sqrtT;
        return rebate * std::exp(-r * T) *
               (N(eta * (x2 - sigma * sqrtT)) -
                std::pow(H / S, 2.0 * mu) * N(eta * (y2 - sigma * sqrtT)));
    }
    // Rebate paid at the hit time if a knock-out option knocks out.
    double F() const {
        if (rebate == 0.0) return 0.0;
        double z = std::log(H / S) / (sigma * sqrtT) + lambda * sigma * sqrtT;
        return rebate * (std::pow(H / S, mu + lambda) * N(eta * z) +
                         std::pow(H / S, mu - lambda) *
                             N(eta * (z - 2.0 * lambda * sigma * sqrtT)));
    }
};

}  // namespace detail

/// Continuously monitored single-barrier option (Reiner-Rubinstein).
/// `rebate` is paid at hit for knock-outs, at expiry for unexercised knock-ins.
/// q is a continuous dividend yield (cost of carry b = r - q).
inline double barrier_price(OptionType type, BarrierType bt, double S, double K,
                            double H, double T, double r, double q, double sigma,
                            double rebate = 0.0) {
    require(S > 0.0 && K > 0.0 && H > 0.0, "barrier: S, K, H must be positive");
    require(sigma > 0.0, "barrier: volatility must be positive");
    require(T > 0.0, "barrier: time to expiry must be positive");
    double b = r - q;

    // Already-triggered barriers collapse to vanilla / rebate.
    bool down = is_down_barrier(bt);
    bool hit = down ? (S <= H) : (S >= H);
    if (hit) {
        if (is_knock_in(bt)) return gbs_price(type, S, K, T, r, b, sigma);
        return rebate;  // knocked out immediately; rebate paid now
    }

    detail::BarrierTerms t(type, bt, S, K, H, T, r, b, sigma, rebate);
    bool call = (type == OptionType::Call);
    bool k_above = (K > H);

    double price = 0.0;
    switch (bt) {
        case BarrierType::DownIn:
            if (call)
                price = k_above ? t.C() + t.E() : t.A() - t.B() + t.D() + t.E();
            else
                price = k_above ? t.B() - t.C() + t.D() + t.E() : t.A() + t.E();
            break;
        case BarrierType::UpIn:
            if (call)
                price = k_above ? t.A() + t.E() : t.B() - t.C() + t.D() + t.E();
            else
                price = k_above ? t.A() - t.B() + t.D() + t.E() : t.C() + t.E();
            break;
        case BarrierType::DownOut:
            if (call)
                price = k_above ? t.A() - t.C() + t.F() : t.B() - t.D() + t.F();
            else
                price = k_above ? t.A() - t.B() + t.C() - t.D() + t.F() : t.F();
            break;
        case BarrierType::UpOut:
            if (call)
                price = k_above ? t.F() : t.A() - t.B() + t.C() - t.D() + t.F();
            else
                price = k_above ? t.B() - t.D() + t.F() : t.A() - t.C() + t.F();
            break;
    }
    return price;
}

}  // namespace opal
