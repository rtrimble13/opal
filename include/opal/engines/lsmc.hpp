// Longstaff-Schwartz least-squares Monte Carlo for American options under
// GBM and Heston dynamics. Continuation values are regressed on polynomial
// bases over in-the-money paths only (the standard estimator; produces a
// slightly low-biased price).
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include "opal/core/types.hpp"
#include "opal/math/solvers.hpp"
#include "opal/models/heston.hpp"

namespace opal {

struct LsmcConfig {
    std::size_t paths = 50000;
    int steps = 50;          // total exercise opportunities over the option life
    std::uint64_t seed = 42;
};

namespace detail {

/// Regression of y on basis columns over the selected rows, returning the
/// fitted values for those rows. Basis is row-major paths x nb. Uses normal
/// equations with a tiny ridge for numerical safety.
inline void regress_fit(const std::vector<double>& basis, std::size_t nb,
                        const std::vector<std::size_t>& rows,
                        const std::vector<double>& y,
                        std::vector<double>& fitted) {
    std::vector<double> XtX(nb * nb, 0.0), Xty(nb, 0.0);
    for (std::size_t r : rows) {
        const double* x = &basis[r * nb];
        for (std::size_t i = 0; i < nb; ++i) {
            Xty[i] += x[i] * y[r];
            for (std::size_t j = i; j < nb; ++j) XtX[i * nb + j] += x[i] * x[j];
        }
    }
    for (std::size_t i = 0; i < nb; ++i) {
        for (std::size_t j = 0; j < i; ++j) XtX[i * nb + j] = XtX[j * nb + i];
        XtX[i * nb + i] += 1e-10 * (XtX[i * nb + i] + 1.0);  // ridge
    }
    if (!math::solve_linear(XtX, Xty, nb)) {
        // Degenerate regression: fall back to the mean of y over rows.
        double m = 0.0;
        for (std::size_t r : rows) m += y[r];
        m = rows.empty() ? 0.0 : m / rows.size();
        for (std::size_t r : rows) fitted[r] = m;
        return;
    }
    for (std::size_t r : rows) {
        const double* x = &basis[r * nb];
        double f = 0.0;
        for (std::size_t i = 0; i < nb; ++i) f += Xty[i] * x[i];
        fitted[r] = f;
    }
}

/// Backward induction given simulated spot paths (and optional variance
/// paths for the Heston basis). spots is row-major: paths x steps.
inline McResult lsmc_backward(OptionType type, double S0, double K, double T,
                              double r, const std::vector<double>& spots,
                              const std::vector<double>* variances,
                              std::size_t n_paths, int n_steps) {
    double phi = type_sign(type);
    double dt = T / n_steps;
    auto payoff = [&](double s) { return std::max(phi * (s - K), 0.0); };

    // cashflow[p] realized at time index tau[p] (1-based step index).
    std::vector<double> cashflow(n_paths);
    std::vector<int> tau(n_paths, n_steps);
    for (std::size_t p = 0; p < n_paths; ++p)
        cashflow[p] = payoff(spots[p * n_steps + (n_steps - 1)]);

    std::size_t nb = variances ? 6 : 4;
    std::vector<double> basis(n_paths * nb), y(n_paths), fitted(n_paths);
    std::vector<std::size_t> itm;

    for (int n = n_steps - 1; n >= 1; --n) {
        itm.clear();
        for (std::size_t p = 0; p < n_paths; ++p) {
            double s = spots[p * n_steps + (n - 1)];
            if (payoff(s) <= 0.0) continue;
            itm.push_back(p);
            double x = s / K;  // normalized for conditioning
            double* row = &basis[p * nb];
            row[0] = 1.0;
            row[1] = x;
            row[2] = x * x;
            if (variances) {
                // Heston basis: {1, x, x^2, v, v^2, x v}
                double v = (*variances)[p * n_steps + (n - 1)];
                row[3] = v;
                row[4] = v * v;
                row[5] = x * v;
            } else {
                // GBM basis: {1, x, x^2, x^3}
                row[3] = x * x * x;
            }
            y[p] = cashflow[p] * std::exp(-r * dt * (tau[p] - n));
        }
        if (itm.size() < nb + 2) continue;
        regress_fit(basis, nb, itm, y, fitted);
        for (std::size_t p : itm) {
            double exercise = payoff(spots[p * n_steps + (n - 1)]);
            if (exercise > fitted[p]) {
                cashflow[p] = exercise;
                tau[p] = n;
            }
        }
    }

    double sum = 0.0, sum2 = 0.0;
    for (std::size_t p = 0; p < n_paths; ++p) {
        double v = cashflow[p] * std::exp(-r * dt * tau[p]);
        sum += v;
        sum2 += v * v;
    }
    double n = static_cast<double>(n_paths);
    double mean = sum / n;
    double var = std::max(0.0, sum2 / n - mean * mean);
    McResult res;
    res.price = std::max(mean, payoff(S0));  // immediate exercise floor
    res.std_error = std::sqrt(var / n);
    res.paths = n_paths;
    return res;
}

}  // namespace detail

/// American option under GBM via Longstaff-Schwartz.
inline McResult lsmc_american(OptionType type, double S, double K, double T,
                              double r, double q, double sigma,
                              const LsmcConfig& cfg = {}) {
    require(cfg.paths > 1 && cfg.steps > 0, "lsmc: paths and steps must be positive");
    require(S > 0.0 && K > 0.0 && sigma > 0.0 && T > 0.0,
            "lsmc: S, K, sigma, T must be positive");
    int N = cfg.steps;
    double dt = T / N;
    double drift = (r - q - 0.5 * sigma * sigma) * dt;
    double volstep = sigma * std::sqrt(dt);
    std::mt19937_64 rng(cfg.seed);
    std::normal_distribution<double> norm(0.0, 1.0);

    std::vector<double> spots(cfg.paths * N);
    for (std::size_t p = 0; p < cfg.paths; ++p) {
        double s = S;
        for (int i = 0; i < N; ++i) {
            s *= std::exp(drift + volstep * norm(rng));
            spots[p * N + i] = s;
        }
    }
    return detail::lsmc_backward(type, S, K, T, r, spots, nullptr, cfg.paths, N);
}

/// American option under Heston dynamics via Longstaff-Schwartz
/// (full-truncation Euler paths; variance enters the regression basis).
inline McResult lsmc_american_heston(OptionType type, double S, double K,
                                     double T, double r, double q,
                                     const HestonParams& hp,
                                     const LsmcConfig& cfg = {}) {
    require(cfg.paths > 1 && cfg.steps > 0, "lsmc: paths and steps must be positive");
    require(S > 0.0 && K > 0.0 && T > 0.0, "lsmc: S, K, T must be positive");
    int N = cfg.steps;
    double dt = T / N;
    double sqdt = std::sqrt(dt);
    double rho_c = std::sqrt(1.0 - hp.rho * hp.rho);
    std::mt19937_64 rng(cfg.seed);
    std::normal_distribution<double> norm(0.0, 1.0);

    std::vector<double> spots(cfg.paths * N), vars(cfg.paths * N);
    for (std::size_t p = 0; p < cfg.paths; ++p) {
        double s = S, v = hp.v0;
        for (int i = 0; i < N; ++i) {
            double z1 = norm(rng);
            double z2 = hp.rho * z1 + rho_c * norm(rng);
            double vp = std::max(v, 0.0);
            s *= std::exp((r - q - 0.5 * vp) * dt + std::sqrt(vp) * sqdt * z1);
            v += hp.kappa * (hp.theta - vp) * dt + hp.xi * std::sqrt(vp) * sqdt * z2;
            spots[p * N + i] = s;
            vars[p * N + i] = std::max(v, 0.0);
        }
    }
    return detail::lsmc_backward(type, S, K, T, r, spots, &vars, cfg.paths, N);
}

}  // namespace opal
