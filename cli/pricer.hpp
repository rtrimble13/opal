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

        if (t.model == "heston") {
            t.heston.v0 = a.get_num("v0", t.vol > 0 ? t.vol * t.vol : 0.04);
            t.heston.kappa = a.get_num("kappa", 1.5);
            t.heston.theta = a.get_num("theta", t.heston.v0);
            t.heston.xi = a.get_num("xi", 0.3);
            t.heston.rho = a.get_num("rho", -0.5);
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

        if (model == "heston") return price_heston(m, s, rr, tt, cfg);

        if (instrument == "european" || instrument == "american") {
            ExerciseStyle st = (instrument == "american")
                                   ? ExerciseStyle::American
                                   : ExerciseStyle::European;
            if (m == "analytic") {
                if (st == ExerciseStyle::American)
                    throw std::runtime_error(
                        "no analytic price for american options; use --method "
                        "lr, crr, trinomial or pde");
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
                if (st == ExerciseStyle::American)
                    throw std::runtime_error("mc engine is European-style only");
                cfg.steps = 1;
                return mc_gbm(vanilla_payoff(type, K), s, tt, rr, qq, v, cfg).price;
            }
        } else if (instrument == "digital-cash") {
            if (m == "mc") {
                cfg.steps = 1;
                return mc_gbm(digital_payoff(type, K, cash), s, tt, rr, qq, v, cfg)
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
                return mc_gbm(
                           barrier_payoff(type, bt, K, barrier, s, rebate, rr, tt),
                           s, tt, rr, qq, v, cfg)
                    .price;
            return barrier_price(type, bt, s, K, barrier, tt, rr, qq, v, rebate);
        } else if (instrument == "asian-arith") {
            if (m == "analytic")
                return arithmetic_asian_price(type, s, K, tt, rr, qq, v);
            return mc_arithmetic_asian(type, s, K, tt, rr, qq, v, cfg).price;
        } else if (instrument == "asian-geo") {
            if (m == "mc")
                return mc_gbm(geometric_asian_payoff(type, K), s, tt, rr, qq, v,
                              cfg)
                    .price;
            if (steps > 0)
                return discrete_geometric_asian_price(type, s, K, tt, rr, qq, v,
                                                      static_cast<int>(steps));
            return geometric_asian_price(type, s, K, tt, rr, qq, v);
        } else if (instrument == "lookback-float") {
            if (m == "mc")
                return mc_gbm(lookback_payoff(type, StrikeStyle::Floating, 0.0, s),
                              s, tt, rr, qq, v, cfg)
                    .price;
            return floating_lookback_price(type, s, tt, rr, qq, v, extremum);
        } else if (instrument == "lookback-fixed") {
            if (m == "mc")
                return mc_gbm(lookback_payoff(type, StrikeStyle::Fixed, K, s), s,
                              tt, rr, qq, v, cfg)
                    .price;
            return fixed_lookback_price(type, s, K, tt, rr, qq, v, extremum);
        }
        throw std::runtime_error("unknown instrument '" + instrument +
                                 "' (see `opal help`)");
    }

    /// Monte Carlo standard error when the resolved method is mc, else NaN.
    double mc_std_error() const {
        if (resolved_method() != "mc") return std::nan("");
        McConfig cfg;
        cfg.paths = static_cast<std::size_t>(paths);
        cfg.steps = steps > 0 ? static_cast<int>(steps) : 252;
        cfg.seed = static_cast<std::uint64_t>(seed);
        if (instrument == "asian-arith" && model != "heston")
            return mc_arithmetic_asian(type, S, K, T, r, q, vol, cfg).std_error;
        return std::nan("");
    }

private:
    double price_heston(const std::string& m, double s, double rr, double tt,
                        McConfig cfg) const {
        if (instrument == "european" && m != "mc")
            return heston_price(type, s, K, tt, rr, q, heston);
        cfg.antithetic = false;
        PathPayoff payoff;
        if (instrument == "european")
            payoff = vanilla_payoff(type, K);
        else if (instrument == "digital-cash")
            payoff = digital_payoff(type, K, cash);
        else if (instrument.rfind("barrier-", 0) == 0)
            payoff =
                barrier_payoff(type, barrier_type(), K, barrier, s, rebate, rr, tt);
        else if (instrument == "asian-arith")
            payoff = [phi = type_sign(type), K = K](const std::vector<double>& p) {
                double sum = 0;
                for (double x : p) sum += x;
                return std::max(phi * (sum / p.size() - K), 0.0);
            };
        else if (instrument == "asian-geo")
            payoff = geometric_asian_payoff(type, K);
        else if (instrument == "lookback-float")
            payoff = lookback_payoff(type, StrikeStyle::Floating, 0.0, s);
        else if (instrument == "lookback-fixed")
            payoff = lookback_payoff(type, StrikeStyle::Fixed, K, s);
        else
            throw std::runtime_error("instrument '" + instrument +
                                     "' not supported under heston");
        return mc_heston(payoff, s, tt, rr, q, heston, cfg).price;
    }
};

inline bool is_rates_instrument(const std::string& inst) {
    return inst == "cap" || inst == "floor" || inst == "caplet" ||
           inst == "floorlet" || inst == "swaption" || inst == "zcb-option";
}

}  // namespace opal::cli
