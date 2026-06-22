// Standard normal distribution utilities.
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace opal::math {

inline constexpr double PI = 3.1415926535897932384626434;
inline constexpr double INV_SQRT_2PI = 0.3989422804014326779399461;
inline constexpr double SQRT_2 = 1.4142135623730950488016887;

/// Standard normal probability density function.
inline double norm_pdf(double x) {
    return INV_SQRT_2PI * std::exp(-0.5 * x * x);
}

/// Standard normal cumulative distribution function (double precision via erfc).
inline double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / SQRT_2);
}

/// Inverse standard normal CDF (Acklam's algorithm with one Halley refinement).
/// Accurate to ~1e-15 over (0, 1). Public utility (quantiles, normal-variate
/// generation, VaR); intentionally provided even though the bundled pricers do
/// not all call it (#12).
inline double norm_ppf(double p) {
    if (p <= 0.0 || p >= 1.0) {
        if (p == 0.0) return -std::numeric_limits<double>::infinity();
        if (p == 1.0) return std::numeric_limits<double>::infinity();
        throw std::domain_error("norm_ppf: p must be in (0, 1)");
    }

    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
                               -2.759285104469687e+02, 1.383577518672690e+02,
                               -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
                               -1.556989798598866e+02, 6.680131188771972e+01,
                               -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
                               -2.400758277161838e+00, -2.549732539343734e+00,
                               4.374664141464968e+00,  2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
                               2.445134137142996e+00, 3.754408661907416e+00};

    const double p_low = 0.02425;
    double x;
    if (p < p_low) {
        double q = std::sqrt(-2.0 * std::log(p));
        x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    } else if (p <= 1.0 - p_low) {
        double q = p - 0.5;
        double r = q * q;
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    } else {
        double q = std::sqrt(-2.0 * std::log(1.0 - p));
        x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }

    // One step of Halley's method against the high-precision CDF.
    double e = norm_cdf(x) - p;
    double u = e / norm_pdf(x);
    x -= u / (1.0 + 0.5 * x * u);
    return x;
}

/// Bivariate normal CDF: P(X < a, Y < b) with correlation rho.
/// Genz's algorithm (1989-2004), accurate to ~1e-14. Used by the two-asset
/// (Stulz/correlation) and compound-option closed forms in models/ (#12).
inline double bivar_norm_cdf(double a, double b, double rho) {
    if (rho > 1.0) rho = 1.0;
    if (rho < -1.0) rho = -1.0;

    static const double w20[] = {0.01761400713915212, 0.04060142980038694,
                                 0.06267204833410906, 0.08327674157670475,
                                 0.1019301198172404,  0.1181945319615184,
                                 0.1316886384491766,  0.1420961093183821,
                                 0.1491729864726037,  0.1527533871307259};
    static const double x20[] = {0.9931285991850949, 0.9639719272779138,
                                 0.9122344282513259, 0.8391169718222188,
                                 0.7463319064601508, 0.6360536807265150,
                                 0.5108670019508271, 0.3737060887154196,
                                 0.2277858511416451, 0.07652652113349733};

    double h = -a, k = -b;
    double bvn = 0.0;
    if (std::fabs(rho) < 0.925) {
        double hs = (h * h + k * k) / 2.0;
        double asr = std::asin(rho);
        for (int i = 0; i < 10; ++i) {
            for (int sign = -1; sign <= 1; sign += 2) {
                double sn = std::sin(asr * (sign * x20[i] + 1.0) / 2.0);
                bvn += w20[i] * std::exp((sn * h * k - hs) / (1.0 - sn * sn));
            }
        }
        bvn = bvn * asr / (4.0 * PI) + norm_cdf(-h) * norm_cdf(-k);
    } else {
        if (rho < 0.0) { k = -k; }
        if (std::fabs(rho) < 1.0) {
            double as = (1.0 - rho) * (1.0 + rho);
            double a2 = std::sqrt(as);
            double bs = (h - k) * (h - k);
            double c = (4.0 - h * k) / 8.0;
            double d = (12.0 - h * k) / 16.0;
            double asr = -(bs / as + h * k) / 2.0;
            if (asr > -100.0) {
                bvn = a2 * std::exp(asr) *
                      (1.0 - c * (bs - as) * (1.0 - d * bs / 5.0) / 3.0 +
                       c * d * as * as / 5.0);
            }
            if (-h * k < 100.0) {
                double bsq = std::sqrt(bs);
                bvn -= std::exp(-h * k / 2.0) * std::sqrt(2.0 * PI) *
                       norm_cdf(-bsq / a2) * bsq *
                       (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
            }
            a2 /= 2.0;
            for (int i = 0; i < 10; ++i) {
                for (int sign = -1; sign <= 1; sign += 2) {
                    double xs = a2 * (sign * x20[i] + 1.0);
                    xs = xs * xs;
                    double rs = std::sqrt(1.0 - xs);
                    double asr2 = -(bs / xs + h * k) / 2.0;
                    if (asr2 > -100.0) {
                        bvn += a2 * w20[i] * std::exp(asr2) *
                               (std::exp(-h * k * (1.0 - rs) / (2.0 * (1.0 + rs))) / rs -
                                (1.0 + c * xs * (1.0 + d * xs)));
                    }
                }
            }
            bvn = -bvn / (2.0 * PI);
        }
        if (rho > 0.0) {
            bvn += norm_cdf(-std::max(h, k));
        } else {
            bvn = -bvn;
            if (k > h) bvn += norm_cdf(k) - norm_cdf(h);
        }
    }
    return bvn;
}

}  // namespace opal::math
