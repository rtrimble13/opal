// Root finding and numerical integration.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace opal::math {

/// Brent's method for root finding on a bracketing interval [lo, hi].
/// f(lo) and f(hi) must have opposite signs.
inline double brent(const std::function<double(double)>& f, double lo, double hi,
                    double tol = 1e-12, int max_iter = 200) {
    double a = lo, b = hi;
    double fa = f(a), fb = f(b);
    if (fa * fb > 0.0)
        throw std::runtime_error("brent: root not bracketed");
    if (std::fabs(fa) < std::fabs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }
    double c = a, fc = fa, d = b - a, s = b, fs = fb;
    bool mflag = true;
    for (int i = 0; i < max_iter; ++i) {
        if (fb == 0.0 || std::fabs(b - a) < tol) return b;
        if (fa != fc && fb != fc) {
            // inverse quadratic interpolation
            s = a * fb * fc / ((fa - fb) * (fa - fc)) +
                b * fa * fc / ((fb - fa) * (fb - fc)) +
                c * fa * fb / ((fc - fa) * (fc - fb));
        } else {
            s = b - fb * (b - a) / (fb - fa);  // secant
        }
        double lo2 = (3.0 * a + b) / 4.0;
        bool cond = !((s > std::min(lo2, b) && s < std::max(lo2, b))) ||
                    (mflag && std::fabs(s - b) >= std::fabs(b - c) / 2.0) ||
                    (!mflag && std::fabs(s - b) >= std::fabs(c - d) / 2.0) ||
                    (mflag && std::fabs(b - c) < tol) ||
                    (!mflag && std::fabs(c - d) < tol);
        if (cond) {
            s = (a + b) / 2.0;  // bisection
            mflag = true;
        } else {
            mflag = false;
        }
        fs = f(s);
        d = c;
        c = b;
        fc = fb;
        if (fa * fs < 0.0) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        if (std::fabs(fa) < std::fabs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }
    return b;
}

/// Newton's method with Brent fallback. Returns root of f given derivative df.
inline double newton_safe(const std::function<double(double)>& f,
                          const std::function<double(double)>& df, double guess,
                          double lo, double hi, double tol = 1e-12,
                          int max_iter = 100) {
    double x = guess;
    for (int i = 0; i < max_iter; ++i) {
        double fx = f(x);
        if (std::fabs(fx) < tol) return x;
        double dfx = df(x);
        if (dfx == 0.0 || !std::isfinite(dfx)) break;
        double step = fx / dfx;
        double xn = x - step;
        if (xn <= lo || xn >= hi || !std::isfinite(xn)) break;
        if (std::fabs(step) < tol * std::max(1.0, std::fabs(x))) return xn;
        x = xn;
    }
    return brent(f, lo, hi, tol, 200);
}

/// Adaptive Simpson quadrature on [a, b].
namespace detail {
inline double simpson_step(const std::function<double(double)>& f, double a,
                           double b, double fa, double fm, double fb, double whole,
                           double tol, int depth) {
    double m = 0.5 * (a + b);
    double lm = 0.5 * (a + m), rm = 0.5 * (m + b);
    double flm = f(lm), frm = f(rm);
    double left = (m - a) / 6.0 * (fa + 4.0 * flm + fm);
    double right = (b - m) / 6.0 * (fm + 4.0 * frm + fb);
    double delta = left + right - whole;
    if (depth <= 0 || std::fabs(delta) <= 15.0 * tol)
        return left + right + delta / 15.0;
    return simpson_step(f, a, m, fa, flm, fm, left, tol / 2.0, depth - 1) +
           simpson_step(f, m, b, fm, frm, fb, right, tol / 2.0, depth - 1);
}
}  // namespace detail

inline double integrate(const std::function<double(double)>& f, double a, double b,
                        double tol = 1e-10, int max_depth = 30) {
    double m = 0.5 * (a + b);
    double fa = f(a), fm = f(m), fb = f(b);
    double whole = (b - a) / 6.0 * (fa + 4.0 * fm + fb);
    return detail::simpson_step(f, a, b, fa, fm, fb, whole, tol, max_depth);
}

/// Solve the dense linear system A x = b in place via Gaussian elimination
/// with partial pivoting. A is row-major n x n. Returns false if singular.
inline bool solve_linear(std::vector<double>& A, std::vector<double>& b,
                         std::size_t n) {
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t piv = col;
        for (std::size_t row = col + 1; row < n; ++row)
            if (std::fabs(A[row * n + col]) > std::fabs(A[piv * n + col]))
                piv = row;
        if (std::fabs(A[piv * n + col]) < 1e-14) return false;
        if (piv != col) {
            for (std::size_t k = 0; k < n; ++k)
                std::swap(A[col * n + k], A[piv * n + k]);
            std::swap(b[col], b[piv]);
        }
        for (std::size_t row = col + 1; row < n; ++row) {
            double f = A[row * n + col] / A[col * n + col];
            if (f == 0.0) continue;
            for (std::size_t k = col; k < n; ++k) A[row * n + k] -= f * A[col * n + k];
            b[row] -= f * b[col];
        }
    }
    for (std::size_t row = n; row-- > 0;) {
        double s = b[row];
        for (std::size_t k = row + 1; k < n; ++k) s -= A[row * n + k] * b[k];
        b[row] = s / A[row * n + row];
    }
    return true;
}

}  // namespace opal::math
