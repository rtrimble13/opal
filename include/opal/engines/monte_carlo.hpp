// Monte Carlo pricing engines.
//   - GBM path simulation (exact log-normal stepping)
//   - antithetic variates; control variates for arithmetic Asians
//   - Heston dynamics via full-truncation Euler
// Path payoffs receive the simulated path EXCLUDING S0 (path[i] = S(t_{i+1}),
// t_i = (i+1) T / steps) and return the undiscounted payoff at T.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "opal/core/types.hpp"
#include "opal/math/normal.hpp"
#include "opal/models/black_scholes.hpp"
#include "opal/models/heston.hpp"

namespace opal {

struct McConfig {
    std::size_t paths = 100000;
    int steps = 252;
    std::uint64_t seed = 42;
    bool antithetic = true;
};

using PathPayoff = std::function<double(const std::vector<double>&)>;

/// Generic GBM Monte Carlo. Returns discounted price and standard error.
inline McResult mc_gbm(const PathPayoff& payoff, double S, double T, double r,
                       double q, double sigma, const McConfig& cfg = {}) {
    require(cfg.paths > 0 && cfg.steps > 0, "mc: paths and steps must be positive");
    require(S > 0.0 && sigma >= 0.0 && T > 0.0, "mc: invalid market inputs");
    double dt = T / cfg.steps;
    double drift = (r - q - 0.5 * sigma * sigma) * dt;
    double volstep = sigma * std::sqrt(dt);
    double df = std::exp(-r * T);

    std::mt19937_64 rng(cfg.seed);
    std::normal_distribution<double> norm(0.0, 1.0);

    std::vector<double> path(cfg.steps), anti(cfg.steps), z(cfg.steps);
    double sum = 0.0, sum2 = 0.0;
    for (std::size_t p = 0; p < cfg.paths; ++p) {
        double s1 = S, s2 = S;
        for (int i = 0; i < cfg.steps; ++i) {
            z[i] = norm(rng);
            s1 *= std::exp(drift + volstep * z[i]);
            path[i] = s1;
            if (cfg.antithetic) {
                s2 *= std::exp(drift - volstep * z[i]);
                anti[i] = s2;
            }
        }
        double y = payoff(path);
        if (cfg.antithetic) y = 0.5 * (y + payoff(anti));
        sum += y;
        sum2 += y * y;
    }
    double n = static_cast<double>(cfg.paths);
    double mean = sum / n;
    double var = std::max(0.0, sum2 / n - mean * mean);
    McResult res;
    res.price = df * mean;
    res.std_error = df * std::sqrt(var / n);
    res.paths = cfg.paths;
    return res;
}

// ---------------------------------------------------------------------------
// Closed form for the DISCRETELY monitored geometric Asian (used both as a
// pricer and as the control variate mean for arithmetic Asians).
// ---------------------------------------------------------------------------
inline double discrete_geometric_asian_price(OptionType type, double S, double K,
                                             double T, double r, double q,
                                             double sigma, int n_obs) {
    require(n_obs >= 1, "asian: need at least one observation");
    double b = r - q;
    double nu = b - 0.5 * sigma * sigma;
    double n = static_cast<double>(n_obs);
    // Observations at t_i = i T / n, i = 1..n.
    double mu = std::log(S) + nu * T * (n + 1.0) / (2.0 * n);
    double var = sigma * sigma * T * (n + 1.0) * (2.0 * n + 1.0) / (6.0 * n * n);
    double sd = std::sqrt(var);
    double d1 = (mu + var - std::log(K)) / sd;
    double d2 = d1 - sd;
    double phi = type_sign(type);
    double fwd_g = std::exp(mu + 0.5 * var);
    return std::exp(-r * T) * phi *
           (fwd_g * math::norm_cdf(phi * d1) - K * math::norm_cdf(phi * d2));
}

namespace detail {
inline double avg(const std::vector<double>& v) {
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}
inline double geo_avg(const std::vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += std::log(x);
    return std::exp(s / v.size());
}
}  // namespace detail

/// Arithmetic Asian with the discrete geometric Asian as control variate.
inline McResult mc_arithmetic_asian(OptionType type, double S, double K, double T,
                                    double r, double q, double sigma,
                                    const McConfig& cfg = {},
                                    StrikeStyle strike = StrikeStyle::Fixed) {
    double phi = type_sign(type);
    bool use_cv = (strike == StrikeStyle::Fixed);
    double cv_mean =
        use_cv ? discrete_geometric_asian_price(type, S, K, T, r, q, sigma,
                                                cfg.steps)
               : 0.0;
    double df = std::exp(-r * T);

    double dt = T / cfg.steps;
    double drift = (r - q - 0.5 * sigma * sigma) * dt;
    double volstep = sigma * std::sqrt(dt);
    std::mt19937_64 rng(cfg.seed);
    std::normal_distribution<double> norm(0.0, 1.0);

    std::vector<double> ys, xs;
    ys.reserve(cfg.paths);
    if (use_cv) xs.reserve(cfg.paths);
    std::vector<double> path(cfg.steps), anti(cfg.steps), z(cfg.steps);

    auto target = [&](const std::vector<double>& pth) {
        double A = detail::avg(pth);
        if (strike == StrikeStyle::Fixed) return std::max(phi * (A - K), 0.0);
        return std::max(phi * (pth.back() - A), 0.0);  // floating strike
    };
    auto control = [&](const std::vector<double>& pth) {
        double G = detail::geo_avg(pth);
        return std::max(phi * (G - K), 0.0);
    };

    for (std::size_t p = 0; p < cfg.paths; ++p) {
        double s1 = S, s2 = S;
        for (int i = 0; i < cfg.steps; ++i) {
            z[i] = norm(rng);
            s1 *= std::exp(drift + volstep * z[i]);
            path[i] = s1;
            if (cfg.antithetic) {
                s2 *= std::exp(drift - volstep * z[i]);
                anti[i] = s2;
            }
        }
        double y = target(path);
        double x = use_cv ? control(path) : 0.0;
        if (cfg.antithetic) {
            y = 0.5 * (y + target(anti));
            if (use_cv) x = 0.5 * (x + control(anti));
        }
        ys.push_back(df * y);
        if (use_cv) xs.push_back(df * x);
    }

    double n = static_cast<double>(ys.size());
    double ymean = std::accumulate(ys.begin(), ys.end(), 0.0) / n;
    McResult res;
    res.paths = cfg.paths;
    if (!use_cv) {
        double s2sum = 0.0;
        for (double y : ys) s2sum += (y - ymean) * (y - ymean);
        res.price = ymean;
        res.std_error = std::sqrt(s2sum / (n * (n - 1.0)));
        return res;
    }
    double xmean = std::accumulate(xs.begin(), xs.end(), 0.0) / n;
    double cov = 0.0, varx = 0.0;
    for (std::size_t i = 0; i < ys.size(); ++i) {
        cov += (ys[i] - ymean) * (xs[i] - xmean);
        varx += (xs[i] - xmean) * (xs[i] - xmean);
    }
    double beta = (varx > 0.0) ? cov / varx : 0.0;
    double adj_sum = 0.0, adj_sum2 = 0.0;
    for (std::size_t i = 0; i < ys.size(); ++i) {
        double a = ys[i] - beta * (xs[i] - cv_mean);
        adj_sum += a;
        adj_sum2 += a * a;
    }
    double amean = adj_sum / n;
    double avar = std::max(0.0, adj_sum2 / n - amean * amean);
    res.price = amean;
    res.std_error = std::sqrt(avar / n);
    return res;
}

// ---------------------------------------------------------------------------
// Convenience payoff builders for the generic engine.
// ---------------------------------------------------------------------------

inline PathPayoff vanilla_payoff(OptionType type, double K) {
    double phi = type_sign(type);
    return [phi, K](const std::vector<double>& p) {
        return std::max(phi * (p.back() - K), 0.0);
    };
}

inline PathPayoff digital_payoff(OptionType type, double K, double cash = 1.0) {
    double phi = type_sign(type);
    return [phi, K, cash](const std::vector<double>& p) {
        return phi * (p.back() - K) > 0.0 ? cash : 0.0;
    };
}

inline PathPayoff geometric_asian_payoff(OptionType type, double K) {
    double phi = type_sign(type);
    return [phi, K](const std::vector<double>& p) {
        return std::max(phi * (detail::geo_avg(p) - K), 0.0);
    };
}

inline PathPayoff lookback_payoff(OptionType type, StrikeStyle strike, double K,
                                  double S0) {
    double phi = type_sign(type);
    return [phi, strike, K, S0, type](const std::vector<double>& p) {
        double mx = S0, mn = S0;
        for (double s : p) {
            if (s > mx) mx = s;
            if (s < mn) mn = s;
        }
        if (strike == StrikeStyle::Fixed) {
            double extr = (type == OptionType::Call) ? mx : mn;
            return std::max(phi * (extr - K), 0.0);
        }
        // floating: call pays S_T - min, put pays max - S_T
        return (type == OptionType::Call) ? p.back() - mn : mx - p.back();
    };
}

/// Discretely monitored barrier payoff (monitoring at every step).
/// Rebate conventions match the analytic engine: knock-out rebates are paid
/// at the hit time (pass r and T so the rebate is accrued to expiry, where
/// the engine discounts it back); knock-in rebates are paid at expiry if the
/// option never knocks in. With the default r = T = 0 the knock-out rebate
/// is treated as paid at expiry instead.
inline PathPayoff barrier_payoff(OptionType type, BarrierType bt, double K,
                                 double H, double S0, double rebate = 0.0,
                                 double r = 0.0, double T = 0.0) {
    double phi = type_sign(type);
    bool down = is_down_barrier(bt);
    bool in = is_knock_in(bt);
    return [=](const std::vector<double>& p) {
        // hit_idx: -2 = never hit, -1 = hit at inception, else step index.
        int hit_idx = (down ? (S0 <= H) : (S0 >= H)) ? -1 : -2;
        if (hit_idx == -2)
            for (std::size_t i = 0; i < p.size(); ++i) {
                if (down ? (p[i] <= H) : (p[i] >= H)) {
                    hit_idx = static_cast<int>(i);
                    break;
                }
            }
        bool hit = (hit_idx != -2);
        double vanilla = std::max(phi * (p.back() - K), 0.0);
        if (in) return hit ? vanilla : rebate;
        if (!hit) return vanilla;
        double t_hit =
            (hit_idx < 0) ? 0.0
                          : (hit_idx + 1) * T / static_cast<double>(p.size());
        return rebate * std::exp(r * (T - t_hit));
    };
}

// ---------------------------------------------------------------------------
// Path-payoff dispatch: a single source of truth for selecting a path payoff
// by name. The CLI and the Python bindings previously each duplicated this
// selection and had drifted apart (e.g. the Python barrier payoff dropped the
// rebate / r / T arguments the CLI passed); routing both through here keeps
// them in lock-step.
// ---------------------------------------------------------------------------
struct PayoffParams {
    OptionType type = OptionType::Call;
    double K = 0.0;
    double cash = 1.0;                            // digital cash payout
    BarrierType barrier_type = BarrierType::DownOut;
    double barrier = 0.0;                         // barrier level H
    double S0 = 0.0;                              // spot at inception
    double rebate = 0.0;                          // barrier rebate
    double r = 0.0;                               // knock-out rebate accrual
    double T = 0.0;                               // knock-out rebate accrual
};

/// Build a path payoff for `kind`, one of: "vanilla", "digital",
/// "asian-arith" (plain arithmetic average), "asian-geo", "lookback-fixed",
/// "lookback-float", "barrier". Throws on an unknown kind.
inline PathPayoff make_path_payoff(const std::string& kind,
                                   const PayoffParams& p) {
    if (kind == "vanilla") return vanilla_payoff(p.type, p.K);
    if (kind == "digital") return digital_payoff(p.type, p.K, p.cash);
    if (kind == "asian-geo") return geometric_asian_payoff(p.type, p.K);
    if (kind == "asian-arith") {
        double phi = type_sign(p.type);
        double K = p.K;
        return [phi, K](const std::vector<double>& path) {
            double sum = 0.0;
            for (double x : path) sum += x;
            return std::max(phi * (sum / path.size() - K), 0.0);
        };
    }
    if (kind == "lookback-fixed")
        return lookback_payoff(p.type, StrikeStyle::Fixed, p.K, p.S0);
    if (kind == "lookback-float")
        return lookback_payoff(p.type, StrikeStyle::Floating, 0.0, p.S0);
    if (kind == "barrier")
        return barrier_payoff(p.type, p.barrier_type, p.K, p.barrier, p.S0,
                              p.rebate, p.r, p.T);
    throw std::invalid_argument("unknown path payoff '" + kind + "'");
}

// ---------------------------------------------------------------------------
// Heston Monte Carlo (full-truncation Euler).
// ---------------------------------------------------------------------------
inline McResult mc_heston(const PathPayoff& payoff, double S, double T, double r,
                          double q, const HestonParams& hp,
                          const McConfig& cfg = {}) {
    require(cfg.paths > 0 && cfg.steps > 0, "mc: paths and steps must be positive");
    // Mirror mc_gbm's market-input guard and validate the Heston parameters (as
    // heston_price does), so bad inputs raise a clear error rather than
    // silently producing garbage (#11).
    require(S > 0.0 && T > 0.0, "mc_heston: spot and expiry must be positive");
    require(hp.v0 >= 0.0 && hp.theta >= 0.0 && hp.kappa > 0.0 && hp.xi > 0.0,
            "mc_heston: invalid Heston parameters (need v0,theta>=0, kappa,xi>0)");
    require(hp.rho >= -1.0 && hp.rho <= 1.0, "mc_heston: rho must be in [-1, 1]");
    double dt = T / cfg.steps;
    double sqdt = std::sqrt(dt);
    double df = std::exp(-r * T);
    double rho_c = std::sqrt(1.0 - hp.rho * hp.rho);

    std::mt19937_64 rng(cfg.seed);
    std::normal_distribution<double> norm(0.0, 1.0);

    std::vector<double> path(cfg.steps);
    double sum = 0.0, sum2 = 0.0;
    for (std::size_t p = 0; p < cfg.paths; ++p) {
        double s = S, v = hp.v0;
        for (int i = 0; i < cfg.steps; ++i) {
            double z1 = norm(rng);
            double z2 = hp.rho * z1 + rho_c * norm(rng);
            double vp = std::max(v, 0.0);
            s *= std::exp((r - q - 0.5 * vp) * dt + std::sqrt(vp) * sqdt * z1);
            v += hp.kappa * (hp.theta - vp) * dt + hp.xi * std::sqrt(vp) * sqdt * z2;
            path[i] = s;
        }
        double y = payoff(path);
        sum += y;
        sum2 += y * y;
    }
    double n = static_cast<double>(cfg.paths);
    double mean = sum / n;
    double var = std::max(0.0, sum2 / n - mean * mean);
    McResult res;
    res.price = df * mean;
    res.std_error = df * std::sqrt(var / n);
    res.paths = cfg.paths;
    return res;
}

}  // namespace opal
