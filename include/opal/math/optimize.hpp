// Derivative-free multidimensional minimization (Nelder-Mead simplex), used by
// the model calibrators. Robust for the small (3-5 parameter), possibly
// non-smooth objectives that arise when fitting SABR/Heston to quotes, where a
// gradient is awkward and a simplex is more forgiving of clamped parameters.
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

namespace opal::math {

struct OptResult {
    std::vector<double> x;  // best point found
    double fval;            // objective there
    int iters;              // iterations performed
    bool converged;         // simplex shrank below tolerance before max_iter
};

/// Minimize `f` from `x0` with a per-coordinate initial step `step` (sized to
/// each parameter's scale). Standard Nelder-Mead with reflection/expansion/
/// contraction/shrink; stops when the simplex spread (in both domain and
/// objective) falls below `tol`.
inline OptResult nelder_mead(const std::function<double(const std::vector<double>&)>& f,
                             const std::vector<double>& x0,
                             const std::vector<double>& step, double tol = 1e-10,
                             int max_iter = 4000) {
    const std::size_t n = x0.size();
    std::vector<std::vector<double>> s(n + 1, x0);  // n+1 vertices
    std::vector<double> fv(n + 1);
    for (std::size_t i = 0; i < n; ++i) s[i + 1][i] += step[i];
    for (std::size_t i = 0; i <= n; ++i) fv[i] = f(s[i]);

    const double alpha = 1.0, gamma = 2.0, rho = 0.5, sigma = 0.5;
    auto order = [&]() {
        std::vector<std::size_t> idx(n + 1);
        for (std::size_t i = 0; i <= n; ++i) idx[i] = i;
        std::sort(idx.begin(), idx.end(),
                  [&](std::size_t a, std::size_t b) { return fv[a] < fv[b]; });
        std::vector<std::vector<double>> ss(n + 1);
        std::vector<double> ff(n + 1);
        for (std::size_t i = 0; i <= n; ++i) { ss[i] = s[idx[i]]; ff[i] = fv[idx[i]]; }
        s.swap(ss);
        fv.swap(ff);
    };

    int it = 0;
    for (; it < max_iter; ++it) {
        order();
        // Convergence: tight spread in both objective and domain.
        double fspread = std::fabs(fv[n] - fv[0]);
        double xspread = 0.0;
        for (std::size_t j = 0; j < n; ++j)
            xspread = std::max(xspread, std::fabs(s[n][j] - s[0][j]));
        if (fspread <= tol && xspread <= tol)
            return {s[0], fv[0], it, true};

        // Centroid of all but the worst.
        std::vector<double> c(n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) c[j] += s[i][j] / static_cast<double>(n);

        auto combine = [&](double t) {
            std::vector<double> x(n);
            for (std::size_t j = 0; j < n; ++j) x[j] = c[j] + t * (c[j] - s[n][j]);
            return x;
        };
        std::vector<double> xr = combine(alpha);
        double fr = f(xr);
        if (fr < fv[0]) {
            std::vector<double> xe = combine(alpha * gamma);
            double fe = f(xe);
            if (fe < fr) { s[n] = xe; fv[n] = fe; } else { s[n] = xr; fv[n] = fr; }
        } else if (fr < fv[n - 1]) {
            s[n] = xr;
            fv[n] = fr;
        } else {
            // Contraction towards the better of (worst, reflection).
            std::vector<double> xc(n);
            bool outside = fr < fv[n];
            for (std::size_t j = 0; j < n; ++j)
                xc[j] = c[j] + rho * ((outside ? xr[j] : s[n][j]) - c[j]);
            double fc = f(xc);
            if (fc < (outside ? fr : fv[n])) {
                s[n] = xc;
                fv[n] = fc;
            } else {
                // Shrink towards the best vertex.
                for (std::size_t i = 1; i <= n; ++i) {
                    for (std::size_t j = 0; j < n; ++j)
                        s[i][j] = s[0][j] + sigma * (s[i][j] - s[0][j]);
                    fv[i] = f(s[i]);
                }
            }
        }
    }
    order();
    return {s[0], fv[0], it, false};
}

}  // namespace opal::math
