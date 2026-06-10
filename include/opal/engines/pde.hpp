// Finite difference (PDE) engine: Crank-Nicolson in log-spot with Rannacher
// start-up (two fully implicit half-steps to damp payoff-kink oscillation).
// Supports European and American vanillas and continuously monitored
// knock-out barriers (knock-ins via in-out parity).
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "opal/core/types.hpp"
#include "opal/models/black_scholes.hpp"

namespace opal {

struct PdeGrid {
    int space_nodes = 400;
    int time_steps = 400;
    double n_std = 6.0;  // domain half-width in standard deviations
};

namespace detail {

/// Thomas algorithm: solves tridiagonal system in-place. a=sub, b=diag, c=super.
inline void thomas_solve(std::vector<double>& a, std::vector<double>& b,
                         std::vector<double>& c, std::vector<double>& d) {
    std::size_t n = b.size();
    for (std::size_t i = 1; i < n; ++i) {
        double w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }
    d[n - 1] /= b[n - 1];
    for (std::size_t i = n - 1; i-- > 0;) d[i] = (d[i] - c[i] * d[i + 1]) / b[i];
}

struct PdeSetup {
    std::vector<double> x;   // log-spot grid
    std::vector<double> s;   // spot grid
    double dx;
    int M;
};

inline PdeSetup make_grid(double S, double K, double T, double sigma,
                          const PdeGrid& g, double lo_override,
                          double hi_override) {
    double width = g.n_std * sigma * std::sqrt(T) + 1e-8;
    double x0 = std::log(S);
    double lo = std::min(x0, std::log(K)) - width;
    double hi = std::max(x0, std::log(K)) + width;
    if (!std::isnan(lo_override)) lo = lo_override;
    if (!std::isnan(hi_override)) hi = hi_override;
    PdeSetup gs;
    gs.M = g.space_nodes;
    gs.dx = (hi - lo) / gs.M;
    gs.x.resize(gs.M + 1);
    gs.s.resize(gs.M + 1);
    for (int i = 0; i <= gs.M; ++i) {
        gs.x[i] = lo + i * gs.dx;
        gs.s[i] = std::exp(gs.x[i]);
    }
    return gs;
}

/// Core solver. boundary_lo/hi give Dirichlet values as a function of time
/// to expiry tau; payoff gives the terminal/intrinsic value at spot s.
template <class Payoff, class BoundLo, class BoundHi>
double pde_solve(double S, double T, double r, double q, double sigma,
                 ExerciseStyle style, const PdeGrid& g, const PdeSetup& gs,
                 Payoff payoff, BoundLo bound_lo, BoundHi bound_hi) {
    int M = gs.M;
    int N = g.time_steps;
    double dt = T / N;
    double b = r - q;
    double nu = b - 0.5 * sigma * sigma;
    double s2 = sigma * sigma;
    double dx = gs.dx;

    // L V = 0.5 s2 V_xx + nu V_x - r V on interior nodes.
    double alpha = 0.5 * s2 / (dx * dx);
    double beta = 0.5 * nu / dx;
    double l_sub = alpha - beta;
    double l_diag = -s2 / (dx * dx) - r;
    double l_sup = alpha + beta;

    std::vector<double> v(M + 1), intrinsic(M + 1);
    for (int i = 0; i <= M; ++i) {
        intrinsic[i] = payoff(gs.s[i]);
        v[i] = intrinsic[i];
    }

    int n_interior = M - 1;
    std::vector<double> a(n_interior), bd(n_interior), c(n_interior),
        rhs(n_interior);

    // Rannacher: first two steps as two implicit half-steps each.
    int extra = 2;  // number of initial steps replaced
    for (int step = 0; step < N; ++step) {
        bool rannacher = (step < extra);
        int substeps = rannacher ? 2 : 1;
        double sub_dt = rannacher ? 0.5 * dt : dt;
        double theta = rannacher ? 1.0 : 0.5;  // implicit vs CN

        for (int ss = 0; ss < substeps; ++ss) {
            // Marching backward from the payoff: tau is time to expiry at the
            // current grid time.
            double tau = (step + (ss + 1.0) / substeps) * dt;
            double vlo = bound_lo(tau);
            double vhi = bound_hi(tau);

            for (int i = 0; i < n_interior; ++i) {
                a[i] = -theta * sub_dt * l_sub;
                bd[i] = 1.0 - theta * sub_dt * l_diag;
                c[i] = -theta * sub_dt * l_sup;
                double expl = (1.0 - theta) * sub_dt *
                              (l_sub * v[i] + l_diag * v[i + 1] + l_sup * v[i + 2]);
                rhs[i] = v[i + 1] + expl;
            }
            rhs[0] += theta * sub_dt * l_sub * vlo;
            rhs[n_interior - 1] += theta * sub_dt * l_sup * vhi;

            thomas_solve(a, bd, c, rhs);
            for (int i = 0; i < n_interior; ++i) v[i + 1] = rhs[i];
            v[0] = vlo;
            v[M] = vhi;

            if (style == ExerciseStyle::American)
                for (int i = 0; i <= M; ++i) v[i] = std::max(v[i], intrinsic[i]);
        }
    }

    // Interpolate at S (quadratic through nearest nodes).
    double x0 = std::log(S);
    int k = static_cast<int>((x0 - gs.x[0]) / dx);
    if (k < 1) k = 1;
    if (k > M - 1) k = M - 1;
    double xm = gs.x[k - 1], xc = gs.x[k], xp = gs.x[k + 1];
    double lm = (x0 - xc) * (x0 - xp) / ((xm - xc) * (xm - xp));
    double lc = (x0 - xm) * (x0 - xp) / ((xc - xm) * (xc - xp));
    double lp = (x0 - xm) * (x0 - xc) / ((xp - xm) * (xp - xc));
    return lm * v[k - 1] + lc * v[k] + lp * v[k + 1];
}

}  // namespace detail

/// European/American vanilla via Crank-Nicolson.
inline double pde_price(OptionType type, ExerciseStyle style, double S, double K,
                        double T, double r, double q, double sigma,
                        const PdeGrid& grid = {}) {
    require(S > 0.0 && K > 0.0 && sigma > 0.0 && T > 0.0,
            "pde: S, K, sigma, T must be positive");
    auto gs = detail::make_grid(S, K, T, sigma, grid, std::nan(""), std::nan(""));
    double phi = type_sign(type);
    auto payoff = [&](double s) { return std::max(phi * (s - K), 0.0); };
    double s_lo = gs.s.front(), s_hi = gs.s.back();
    auto bound_lo = [&](double tau) {
        if (type == OptionType::Put) {
            if (style == ExerciseStyle::American) return K - s_lo;
            return K * std::exp(-r * tau) - s_lo * std::exp(-q * tau);
        }
        return 0.0;
    };
    auto bound_hi = [&](double tau) {
        if (type == OptionType::Call) {
            if (style == ExerciseStyle::American) return s_hi - K;
            return s_hi * std::exp(-q * tau) - K * std::exp(-r * tau);
        }
        return 0.0;
    };
    return detail::pde_solve(S, T, r, q, sigma, style, grid, gs, payoff, bound_lo,
                             bound_hi);
}

/// Continuously monitored knock-out barrier via PDE with the barrier as a
/// Dirichlet boundary (value = rebate). Knock-ins are priced by parity:
/// in = vanilla - out (rebate must be zero for parity pricing).
inline double pde_barrier_price(OptionType type, BarrierType bt, double S,
                                double K, double H, double T, double r, double q,
                                double sigma, double rebate = 0.0,
                                const PdeGrid& grid = {}) {
    require(H > 0.0, "pde: barrier must be positive");
    if (is_knock_in(bt)) {
        require(rebate == 0.0,
                "pde: knock-in rebates not supported (use analytic engine)");
        BarrierType out_bt = is_down_barrier(bt) ? BarrierType::DownOut
                                                 : BarrierType::UpOut;
        double vanilla =
            pde_price(type, ExerciseStyle::European, S, K, T, r, q, sigma, grid);
        return vanilla -
               pde_barrier_price(type, out_bt, S, K, H, T, r, q, sigma, 0.0, grid);
    }

    bool down = is_down_barrier(bt);
    bool hit = down ? (S <= H) : (S >= H);
    if (hit) return rebate;

    double lo_ovr = down ? std::log(H) : std::nan("");
    double hi_ovr = down ? std::nan("") : std::log(H);
    auto gs = detail::make_grid(S, K, T, sigma, grid, lo_ovr, hi_ovr);
    double phi = type_sign(type);
    auto payoff = [&](double s) { return std::max(phi * (s - K), 0.0); };
    double s_lo = gs.s.front(), s_hi = gs.s.back();

    auto bound_lo = [&](double tau) -> double {
        if (down) return rebate;  // knocked out at the barrier
        if (type == OptionType::Put)
            return K * std::exp(-r * tau) - s_lo * std::exp(-q * tau);
        return 0.0;
    };
    auto bound_hi = [&](double tau) -> double {
        if (!down) return rebate;
        if (type == OptionType::Call)
            return s_hi * std::exp(-q * tau) - K * std::exp(-r * tau);
        return 0.0;
    };
    return detail::pde_solve(S, T, r, q, sigma, ExerciseStyle::European, grid, gs,
                             payoff, bound_lo, bound_hi);
}

}  // namespace opal
