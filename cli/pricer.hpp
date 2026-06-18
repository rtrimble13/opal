// Instrument specification and pricing dispatch shared by the CLI commands.
#pragma once

#include <map>
#include <stdexcept>
#include <string>

#include "args.hpp"
#include "opal/opal.hpp"

namespace opal::cli {

/// Equity-style trade: everything needed to price, re-priceable under
/// bumped market data (for scenarios and numerical greeks).
struct EquityTrade {
    std::string instrument = "european";
    std::string model = "bsm";      // bsm | black76 | bachelier | heston
    std::string method = "auto";    // auto | analytic | crr | lr | trinomial | pde | mc
    OptionType type = OptionType::Call;

    double S = 0, K = 0, T = 0, r = 0, q = 0, vol = 0;
    double barrier = 0, rebate = 0, cash = 1.0, strike2 = 0, extremum = 0;

    long steps = 0;       // 0 = engine default
    long paths = 100000;
    long seed = 42;
    long grid = 400;

    DividendSchedule dividends;
    HestonParams heston{0.04, 1.5, 0.04, 0.3, -0.5};

    static EquityTrade from_args(const Args& a) {
        EquityTrade t;
        t.instrument = a.get_str("instrument", "european");
        t.model = a.get_str("model", "bsm");
        t.method = a.get_str("method", "auto");
        std::string ty = a.get_str("type", "call");
        if (ty == "call" || ty == "c")
            t.type = OptionType::Call;
        else if (ty == "put" || ty == "p")
            t.type = OptionType::Put;
        else
            throw std::runtime_error("--type must be call or put, got '" + ty + "'");

        t.S = a.require_num("spot");
        t.K = a.get_num("strike", t.S);  // ATM default
        t.T = a.get_expiry();
        t.r = a.get_num("rate", 0.0);
        t.q = a.get_num("div", 0.0);
        t.vol = a.get_num("vol", 0.0);
        if (t.model != "heston" && t.vol <= 0.0)
            throw std::runtime_error("missing required option --vol");

        t.barrier = a.get_num("barrier", 0.0);
        t.rebate = a.get_num("rebate", 0.0);
        t.cash = a.get_num("cash", 1.0);
        t.strike2 = a.get_num("payoff-strike", t.K);
        t.extremum = a.get_num("extremum", 0.0);

        t.steps = a.get_int("steps", 0);
        t.paths = a.get_int("paths", 100000);
        t.seed = a.get_int("seed", 42);
        t.grid = a.get_int("grid", 400);

        // --dividends t:amount[,t:amount...] (discrete cash dividends)
        std::string divspec = a.get_str("dividends");
        if (!divspec.empty()) {
            std::size_t pos = 0;
            while (pos < divspec.size()) {
                std::size_t comma = divspec.find(',', pos);
                std::string item = divspec.substr(
                    pos, comma == std::string::npos ? std::string::npos
                                                    : comma - pos);
                std::size_t colon = item.find(':');
                if (colon == std::string::npos)
                    throw std::runtime_error(
                        "--dividends: expected t:amount[,t:amount...], got '" +
                        divspec + "'");
                t.dividends.push_back({std::stod(item.substr(0, colon)),
                                       std::stod(item.substr(colon + 1))});
                if (comma == std::string::npos) break;
                pos = comma + 1;
            }
            if (t.instrument != "european" && t.instrument != "american")
                throw std::runtime_error(
                    "--dividends is supported for european and american "
                    "instruments only");
            if (t.model != "bsm")
                throw std::runtime_error(
                    "--dividends requires the bsm model");
        }

        if (t.model == "heston") {
            t.heston.v0 = a.get_num("v0", t.vol > 0 ? t.vol * t.vol : 0.04);
            t.heston.kappa = a.get_num("kappa", 1.5);
            t.heston.theta = a.get_num("theta", t.heston.v0);
            t.heston.xi = a.get_num("xi", 0.3);
            t.heston.rho = a.get_num("rho", -0.5);
            // Anchor a base vol level for vol-bump based tooling (numerical
            // greeks, scenario grid). When --vol is omitted under Heston we use
            // the initial instantaneous vol sqrt(v0) so a vega/scenario bump has
            // a well-defined point to perturb around.
            if (t.vol <= 0.0) t.vol = std::sqrt(t.heston.v0);
        }
        if (t.instrument.rfind("barrier", 0) == 0 && t.barrier <= 0.0)
            throw std::runtime_error("missing required option --barrier");
        return t;
    }

    BarrierType barrier_type() const {
        static const std::map<std::string, BarrierType> m = {
            {"barrier-down-in", BarrierType::DownIn},
            {"barrier-down-out", BarrierType::DownOut},
            {"barrier-up-in", BarrierType::UpIn},
            {"barrier-up-out", BarrierType::UpOut},
        };
        auto it = m.find(instrument);
        if (it == m.end())
            throw std::runtime_error("unknown barrier instrument '" + instrument + "'");
        return it->second;
    }

    /// Effective numerical method after resolving "auto".
    std::string resolved_method() const {
        if (method != "auto") return method;
        if (model == "heston")
            return (instrument == "european") ? "analytic" : "mc";
        if (instrument == "american") return "lr";
        if (instrument == "asian-arith") return "mc";
        return "analytic";
    }

    double price() const { return price_at(S, vol, r, T); }

    /// Price under bumped market data (spot, vol, rate, expiry).
    double price_at(double s, double v, double rr, double tt) const {
        std::string m = resolved_method();
        McConfig cfg;
        cfg.paths = static_cast<std::size_t>(paths);
        cfg.steps = steps > 0 ? static_cast<int>(steps) : 252;
        cfg.seed = static_cast<std::uint64_t>(seed);
        PdeGrid pg{static_cast<int>(grid), static_cast<int>(grid), 6.0};
        int lattice_steps = steps > 0 ? static_cast<int>(steps) : 0;

        // Model remapping for forward-based models.
        double qq = q;
        if (model == "black76") qq = rr;  // b = 0
        if (model == "bachelier") {
            if (instrument != "european")
                throw std::runtime_error(
                    "bachelier model supports european options only");
            return bachelier_price(type, s, K, tt, rr, v);
        }

        if (model == "heston") {
            // The generic risk tooling (numerical greeks, scenario grid)
            // perturbs the `v` argument. Heston is not driven by a single
            // lognormal vol, so map a vol bump to a parallel shift of the
            // variance level: move v0 and theta together by the change in
            // variance relative to the trade's base vol. The resulting "vega"
            // is the sensitivity to a parallel shift of the vol level, the
            // closest analogue to BSM vega. At v == vol this is a no-op, so the
            // base price is unchanged.
            HestonParams hp = heston;
            if (v > 0.0 && vol > 0.0) {
                double dvar = v * v - vol * vol;
                hp.v0 = std::max(heston.v0 + dvar, 1e-8);
                hp.theta = std::max(heston.theta + dvar, 1e-8);
            }
            return price_heston(m, s, rr, tt, cfg, hp);
        }

        if (instrument == "european" || instrument == "american") {
            ExerciseStyle st = (instrument == "american")
                                   ? ExerciseStyle::American
                                   : ExerciseStyle::European;
            if (!dividends.empty()) {
                if (st == ExerciseStyle::European && m == "analytic")
                    return bsm_discrete_div_price(type, s, K, tt, rr, dividends,
                                                  v, qq);
                if (m == "lr" || m == "analytic")
                    return binomial_lr_discrete_div_price(
                        type, st, s, K, tt, rr, dividends, v,
                        lattice_steps ? lattice_steps : 501, qq);
                throw std::runtime_error(
                    "discrete dividends support --method analytic (european) "
                    "or lr");
            }
            if (m == "analytic") {
                if (st == ExerciseStyle::American)
                    throw std::runtime_error(
                        "no analytic price for american options; use --method "
                        "lr, crr, trinomial, pde or mc (Longstaff-Schwartz)");
                return bsm_price(type, s, K, tt, rr, qq, v);
            }
            if (m == "crr")
                return binomial_crr_price(type, st, s, K, tt, rr, qq, v,
                                          lattice_steps ? lattice_steps : 512);
            if (m == "lr")
                return binomial_lr_price(type, st, s, K, tt, rr, qq, v,
                                         lattice_steps ? lattice_steps : 501);
            if (m == "trinomial")
                return trinomial_price(type, st, s, K, tt, rr, qq, v,
                                       lattice_steps ? lattice_steps : 400);
            if (m == "pde") return pde_price(type, st, s, K, tt, rr, qq, v, pg);
            if (m == "mc") {
                if (st == ExerciseStyle::American) {
                    LsmcConfig lc;
                    lc.paths = static_cast<std::size_t>(paths);
                    lc.steps = steps > 0 ? static_cast<int>(steps) : 50;
                    lc.seed = static_cast<std::uint64_t>(seed);
                    return lsmc_american(type, s, K, tt, rr, qq, v, lc).price;
                }
                cfg.steps = 1;
                return mc_gbm(make_path_payoff(payoff_kind(),
                                               payoff_params(s, rr, tt)),
                              s, tt, rr, qq, v, cfg)
                    .price;
            }
        } else if (instrument == "digital-cash") {
            if (m == "mc") {
                cfg.steps = 1;
                return mc_gbm(make_path_payoff(payoff_kind(),
                                               payoff_params(s, rr, tt)),
                              s, tt, rr, qq, v, cfg)
                    .price;
            }
            return cash_or_nothing_price(type, s, K, tt, rr, qq, v, cash);
        } else if (instrument == "digital-asset") {
            return asset_or_nothing_price(type, s, K, tt, rr, qq, v);
        } else if (instrument == "gap") {
            return gap_option_price(type, s, K, strike2, tt, rr, qq, v);
        } else if (instrument.rfind("barrier-", 0) == 0) {
            BarrierType bt = barrier_type();
            if (m == "pde")
                return pde_barrier_price(type, bt, s, K, barrier, tt, rr, qq, v,
                                         rebate, pg);
            if (m == "mc")
                return mc_gbm(make_path_payoff(payoff_kind(),
                                               payoff_params(s, rr, tt)),
                              s, tt, rr, qq, v, cfg)
                    .price;
            return barrier_price(type, bt, s, K, barrier, tt, rr, qq, v, rebate);
        } else if (instrument == "asian-arith") {
            if (m == "analytic")
                return arithmetic_asian_price(type, s, K, tt, rr, qq, v);
            return mc_arithmetic_asian(type, s, K, tt, rr, qq, v, cfg).price;
        } else if (instrument == "asian-geo") {
            if (m == "mc")
                return mc_gbm(make_path_payoff(payoff_kind(),
                                               payoff_params(s, rr, tt)),
                              s, tt, rr, qq, v, cfg)
                    .price;
            if (steps > 0)
                return discrete_geometric_asian_price(type, s, K, tt, rr, qq, v,
                                                      static_cast<int>(steps));
            return geometric_asian_price(type, s, K, tt, rr, qq, v);
        } else if (instrument == "lookback-float") {
            if (m == "mc")
                return mc_gbm(make_path_payoff(payoff_kind(),
                                               payoff_params(s, rr, tt)),
                              s, tt, rr, qq, v, cfg)
                    .price;
            return floating_lookback_price(type, s, tt, rr, qq, v, extremum);
        } else if (instrument == "lookback-fixed") {
            if (m == "mc")
                return mc_gbm(make_path_payoff(payoff_kind(),
                                               payoff_params(s, rr, tt)),
                              s, tt, rr, qq, v, cfg)
                    .price;
            return fixed_lookback_price(type, s, K, tt, rr, qq, v, extremum);
        }
        throw std::runtime_error("unknown instrument '" + instrument +
                                 "' (see `opal help`)");
    }

    /// Monte Carlo standard error when the resolved method is mc, else NaN.
    /// Mirrors the MC dispatch in price_at for every instrument that has a
    /// Monte Carlo path, so `opal price --method mc` reports the error bar that
    /// justifies the method (previously only asian-arith did).
    double mc_std_error() const {
        if (resolved_method() != "mc") return std::nan("");
        // Instruments/models with no GBM path payoff are always priced
        // analytically regardless of --method.
        if (model == "bachelier" || instrument == "digital-asset" ||
            instrument == "gap")
            return std::nan("");

        McConfig cfg;
        cfg.paths = static_cast<std::size_t>(paths);
        cfg.steps = steps > 0 ? static_cast<int>(steps) : 252;
        cfg.seed = static_cast<std::uint64_t>(seed);
        LsmcConfig lc;
        lc.paths = static_cast<std::size_t>(paths);
        lc.steps = steps > 0 ? static_cast<int>(steps) : 50;
        lc.seed = static_cast<std::uint64_t>(seed);

        if (model == "heston") {
            if (instrument == "american")
                return lsmc_american_heston(type, S, K, T, r, q, heston, lc)
                    .std_error;
            cfg.antithetic = false;
            return mc_heston(make_path_payoff(payoff_kind(), payoff_params(S, r, T)),
                             S, T, r, q, heston, cfg)
                .std_error;
        }

        double qq = (model == "black76") ? r : q;
        if (instrument == "american")
            return lsmc_american(type, S, K, T, r, qq, vol, lc).std_error;
        if (instrument == "asian-arith")
            return mc_arithmetic_asian(type, S, K, T, r, qq, vol, cfg).std_error;
        if (instrument == "european" || instrument == "digital-cash") cfg.steps = 1;
        return mc_gbm(make_path_payoff(payoff_kind(), payoff_params(S, r, T)), S, T,
                      r, qq, vol, cfg)
            .std_error;
    }

    /// First-class Heston greeks for this trade. A European priced
    /// semi-analytically uses the exact `heston_greeks`; every other Heston
    /// instrument (American/exotic, or a European forced onto Monte Carlo)
    /// finite-differences the MC/Longstaff-Schwartz price under common random
    /// numbers (a fixed seed), so the first-order greeks and parameter
    /// sensitivities are stable. Only valid when model == "heston".
    HestonGreeks heston_trade_greeks() const {
        std::string m = resolved_method();
        if (instrument == "european" && m == "analytic")
            return opal::heston_greeks(type, S, K, T, r, q, heston);
        McConfig cfg;
        cfg.paths = static_cast<std::size_t>(paths);
        cfg.steps = steps > 0 ? static_cast<int>(steps) : 252;
        cfg.seed = static_cast<std::uint64_t>(seed);
        HestonPricerFn f = [&](double s, double tt, double rr,
                               const HestonParams& hp) {
            // Fixed seed across bumps => common random numbers.
            return price_heston(m, s, rr, tt, cfg, hp);
        };
        return opal::heston_greeks_fd(f, S, T, r, heston);
    }

    /// True when this trade is priced by a noisy Monte Carlo / LSMC engine, so
    /// its greeks carry sampling noise (used to caption the risk report).
    bool is_monte_carlo() const {
        return model == "heston" ? !(instrument == "european" &&
                                     resolved_method() == "analytic")
                                 : resolved_method() == "mc";
    }

private:
    /// Canonical path-payoff kind for this instrument (CLI vocabulary ->
    /// library vocabulary). Throws for instruments with no path payoff.
    std::string payoff_kind() const {
        if (instrument == "european" || instrument == "american") return "vanilla";
        if (instrument == "digital-cash") return "digital";
        if (instrument == "asian-arith") return "asian-arith";
        if (instrument == "asian-geo") return "asian-geo";
        if (instrument == "lookback-float") return "lookback-float";
        if (instrument == "lookback-fixed") return "lookback-fixed";
        if (instrument.rfind("barrier-", 0) == 0) return "barrier";
        throw std::runtime_error("instrument '" + instrument +
                                 "' has no Monte Carlo path payoff");
    }

    /// Path-payoff parameters under bumped market data (S0 = bumped spot;
    /// r/T carry the knock-out rebate accrual convention for barriers).
    PayoffParams payoff_params(double s, double rr, double tt) const {
        PayoffParams pp;
        pp.type = type;
        pp.K = K;
        pp.cash = cash;
        pp.S0 = s;
        pp.barrier = barrier;
        pp.rebate = rebate;
        pp.r = rr;
        pp.T = tt;
        if (instrument.rfind("barrier-", 0) == 0) pp.barrier_type = barrier_type();
        return pp;
    }

    double price_heston(const std::string& m, double s, double rr, double tt,
                        McConfig cfg, const HestonParams& hp) const {
        if (instrument == "american") {
            LsmcConfig lc;
            lc.paths = static_cast<std::size_t>(paths);
            lc.steps = steps > 0 ? static_cast<int>(steps) : 50;
            lc.seed = static_cast<std::uint64_t>(seed);
            return lsmc_american_heston(type, s, K, tt, rr, q, hp, lc).price;
        }
        if (instrument == "european" && m != "mc")
            return heston_price(type, s, K, tt, rr, q, hp);
        cfg.antithetic = false;
        // Centralized payoff dispatch (#5) with the vol-bump-adjusted Heston
        // params (#3): hp carries the parallel variance-level shift used by the
        // numerical greeks / scenario tooling.
        return mc_heston(make_path_payoff(payoff_kind(), payoff_params(s, rr, tt)),
                         s, tt, rr, q, hp, cfg)
            .price;
    }
};

inline bool is_rates_instrument(const std::string& inst) {
    return inst == "cap" || inst == "floor" || inst == "caplet" ||
           inst == "floorlet" || inst == "swaption" || inst == "zcb-option";
}

}  // namespace opal::cli
