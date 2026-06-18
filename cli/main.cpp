// opal - institutional option pricing command line tool.
#include <fstream>
#include <iostream>
#include <sstream>

#include "args.hpp"
#include "format.hpp"
#include "pricer.hpp"

namespace opal::cli {

// ---------------------------------------------------------------------------
// help
// ---------------------------------------------------------------------------
const char* HELP = R"(opal - institutional option pricing and risk analytics

USAGE
  opal <command> [options]

COMMANDS
  price       Price a single instrument
  greeks      Price with a full risk (greeks) report
  implied     Implied volatility from a market price
  chain       Price a ladder of strikes (option chain)
  scenario    Spot x vol P&L scenario grid for risk management
  portfolio   Price and aggregate risk for a CSV portfolio file
  version     Print version
  help        Show this message

EQUITY INSTRUMENTS  (--instrument / -i)
  european (default)   american             digital-cash    digital-asset
  gap                  barrier-down-in      barrier-down-out
  barrier-up-in        barrier-up-out       asian-arith     asian-geo
  lookback-fixed       lookback-float

RATES INSTRUMENTS
  cap        floor      caplet     floorlet    swaption    zcb-option

MARKET DATA
  -S, --spot       Spot (or forward for black76/bachelier)
  -K, --strike     Strike (default: ATM = spot; a rate for caps/swaptions)
  -T, --expiry     Years (0.5) or date (2026-12-18, ACT/365F from today)
  -r, --rate       Continuously compounded risk-free rate (0.05 or 5%)
  -q, --div        Continuous dividend yield
  -v, --vol        Volatility (lognormal; absolute for bachelier)
  -t, --type       call | put (payer | receiver for swaptions)

MODELS  (--model)
  bsm (default)    Black-Scholes-Merton on spot with dividend yield
  black76          Black-76 on a forward/futures price
  bachelier        Normal model on a forward (absolute vol)
  heston           Heston stochastic vol: --v0 --kappa --theta --xi --rho
  sabr             SABR (rates: swaptions): --alpha --beta --nu --rho
  hullwhite        Hull-White 1F (rates): --mean-rev --hw-sigma

METHODS  (--method, default auto)
  analytic   Closed form (default where available)
  lr         Leisen-Reimer binomial (default for american)   [--steps]
  crr        Cox-Ross-Rubinstein binomial                    [--steps]
  trinomial  Trinomial tree                                  [--steps]
  pde        Crank-Nicolson finite differences               [--grid]
  mc         Monte Carlo (default for asian-arith); for american
             instruments runs Longstaff-Schwartz LSMC        [--paths --steps --seed]

EXOTIC PARAMETERS
  --barrier / -H   Barrier level        --rebate          Barrier rebate
  --cash           Digital cash payout  --payoff-strike   Gap payoff strike
  --extremum       Running min/max for seasoned lookbacks
  --dividends      Discrete cash dividends t:amount[,t:amount...]
                   (european/american, bsm model; escrowed-spot convention)

RATES PARAMETERS
  --tenor      Swap tenor in years (swaption) / T2 of the bond (zcb-option)
  --freq       Payments per year (default 4 for caps, 2 for swaptions)
  --first-fixing   First caplet fixing time (default 1/freq)
  --notional / -n  Notional (default 100)
  --vol-type   lognormal (default) | normal
  --ois-rate   Flat OIS discount rate for multi-curve pricing
               (forwards projected off --rate, cashflows discounted on OIS)

COMMAND-SPECIFIC
  implied:    --price <market price>
  chain:      --strikes lo:hi:step
  scenario:   --spot-range lo:hi:step (%)  --vol-range lo:hi:step (vol pts)
  portfolio:  --file positions.csv

OUTPUT
  -o, --output   table (default) | json | csv (chain/scenario/portfolio)

EXAMPLES
  opal price -S 100 -K 105 -T 0.5 -r 4% -v 22% -t call
  opal price -i american -t put -S 50 -K 55 -T 1 -r 3% -q 1% -v 30%
  opal price -i barrier-up-out -S 100 -K 100 -H 120 -T 1 -r 5% -v 25%
  opal price -i asian-arith -S 100 -K 100 -T 1 -r 5% -v 30% --paths 200000
  opal price --model heston -S 100 -K 100 -T 1 -r 3% --v0 0.04 --kappa 1.5 \
             --theta 0.04 --xi 0.4 --rho -0.6
  opal price -i american -t call -S 100 -K 100 -T 1 -r 4% -v 25% \
             --dividends 0.25:1.50,0.75:1.50
  opal price -i american -t put -S 100 -K 100 -T 1 -r 5% -v 25% --method mc
  opal price -i cap -K 4% -T 5 -r 4.2% -v 30% --freq 4
  opal price -i cap -K 4% -T 5 -r 4.2% --ois-rate 3.8% -v 30% --freq 4
  opal price -i swaption -t payer -K 4% -T 1 --tenor 5 -r 4% -v 25%
  opal greeks -S 100 -K 100 -T 0.5 -r 5% -v 20%
  opal implied -S 100 -K 105 -T 0.5 -r 4% --price 3.85
  opal chain -S 100 -T 0.5 -r 4% -v 25% --strikes 80:120:5
  opal scenario -S 100 -K 100 -T 0.5 -r 4% -v 25% \
                --spot-range -20:20:5 --vol-range -10:10:5
  opal portfolio --file book.csv -o json
)";

// ---------------------------------------------------------------------------
// rates pricing
// ---------------------------------------------------------------------------
void cmd_price_rates(const Args& a, OutputFormat out) {
    std::string inst = a.get_str("instrument");
    std::string model = a.get_str("model", "black76");
    double notional = a.get_num("notional", 100.0);
    double r = a.require_num("rate");
    DiscountCurve curve(r);  // projection curve
    // OIS discounting: --ois-rate sets a separate (flat) discount curve;
    // defaults to single-curve pricing on --rate.
    double ois = a.get_num("ois-rate", r);
    DiscountCurve disc(ois);
    bool dual = a.has("ois-rate");
    Report rep;
    rep.add("instrument", inst);
    rep.add("model", model);
    if (dual) rep.add("ois_rate", ois, 6);

    if (inst == "cap" || inst == "floor") {
        bool is_cap = (inst == "cap");
        double K = a.require_num("strike");
        double T = a.get_expiry();
        double freq = a.get_num("freq", 4.0);
        double tau = 1.0 / freq;
        double first = a.get_num("first-fixing", tau);
        double price;
        if (model == "hullwhite") {
            HullWhiteParams hw{a.get_num("mean-rev", 0.05),
                               a.get_num("hw-sigma", 0.01)};
            price = hw_cap_price(curve, hw, first, T, tau, K, is_cap);
        } else {
            double vol = a.require_num("vol");
            RateVolType vt = a.get_str("vol-type", "lognormal") == "normal"
                                 ? RateVolType::Normal
                                 : RateVolType::Lognormal;
            auto res =
                cap_floor_price(disc, curve, K, vol, first, T, tau, is_cap, vt);
            price = res.price;
            rep.add("caplets", static_cast<double>(res.caplets.size()), 0);
        }
        rep.add("strike", K, 4);
        rep.add("maturity_years", T, 4);
        rep.add("price", notional * price, 6);
        rep.add("notional", notional, 2);
    } else if (inst == "caplet" || inst == "floorlet") {
        double K = a.require_num("strike");
        double T1 = a.get_expiry();
        double freq = a.get_num("freq", 4.0);
        double T2 = T1 + 1.0 / freq;
        bool is_cap = (inst == "caplet");
        double price;
        if (model == "hullwhite") {
            HullWhiteParams hw{a.get_num("mean-rev", 0.05),
                               a.get_num("hw-sigma", 0.01)};
            price = is_cap ? hw_caplet_price(curve, hw, T1, T2, K)
                           : hw_floorlet_price(curve, hw, T1, T2, K);
        } else {
            double vol = a.require_num("vol");
            RateVolType vt = a.get_str("vol-type", "lognormal") == "normal"
                                 ? RateVolType::Normal
                                 : RateVolType::Lognormal;
            auto res =
                cap_floor_price(disc, curve, K, vol, T1, T2, T2 - T1, is_cap, vt);
            price = res.price;
        }
        rep.add("fixing", T1, 4);
        rep.add("payment", T2, 4);
        rep.add("forward", curve.forward_rate(T1, T2), 6);
        rep.add("price", notional * price, 6);
        rep.add("notional", notional, 2);
    } else if (inst == "swaption") {
        std::string ty = a.get_str("type", "payer");
        SwaptionType st =
            (ty == "receiver" || ty == "r") ? SwaptionType::Receiver
                                            : SwaptionType::Payer;
        double K = a.require_num("strike");
        double T = a.get_expiry();
        double tenor = a.require_num("tenor");
        double freq = a.get_num("freq", 2.0);
        SwaptionResult res;
        if (model == "sabr") {
            SabrParams sp{a.get_num("alpha", 0.2), a.get_num("beta", 0.5),
                          a.get_num("rho", -0.3), a.get_num("nu", 0.4)};
            res = sabr_swaption_price(curve, st, K, sp, T, tenor, freq);
            rep.add("sabr_vol",
                    sabr_lognormal_vol(res.forward_swap_rate, K, T, sp), 6);
        } else {
            double vol = a.require_num("vol");
            RateVolType vt = a.get_str("vol-type", "lognormal") == "normal"
                                 ? RateVolType::Normal
                                 : RateVolType::Lognormal;
            res = swaption_price(disc, curve, st, K, vol, T, tenor, freq, vt);
        }
        rep.add("style", st == SwaptionType::Payer ? "payer" : "receiver");
        rep.add("strike", K, 6);
        rep.add("forward_swap_rate", res.forward_swap_rate, 6);
        rep.add("annuity", res.annuity, 6);
        rep.add("price", notional * res.price, 6);
        rep.add("notional", notional, 2);
    } else if (inst == "zcb-option") {
        double K = a.require_num("strike");
        double T1 = a.get_expiry();
        double T2 = a.require_num("tenor");  // bond maturity
        std::string ty = a.get_str("type", "call");
        OptionType ot = (ty == "put" || ty == "p") ? OptionType::Put
                                                   : OptionType::Call;
        HullWhiteParams hw{a.get_num("mean-rev", 0.05), a.get_num("hw-sigma", 0.01)};
        double price = hw_zcb_option_price(ot, curve, hw, T1, T2, K);
        rep.add("type", ty);
        rep.add("option_expiry", T1, 4);
        rep.add("bond_maturity", T2, 4);
        rep.add("price", notional * price, 6);
        rep.add("notional", notional, 2);
    } else {
        throw std::runtime_error("unknown rates instrument '" + inst + "'");
    }
    rep.print(out, "Opal | " + inst);
}

// ---------------------------------------------------------------------------
// price / greeks
// ---------------------------------------------------------------------------
Greeks trade_greeks(const EquityTrade& t) {
    // Analytic greeks for plain BSM europeans, numerical otherwise.
    if (t.instrument == "european" && t.model == "bsm" &&
        t.resolved_method() == "analytic")
        return bsm_greeks(t.type, t.S, t.K, t.T, t.r, t.q, t.vol);
    auto pricer = [&](double s, double tt, double rr, double v) {
        return t.price_at(s, v, rr, tt);
    };
    return numerical_greeks(pricer, t.S, t.T, t.r, t.vol);
}

void cmd_price(const Args& a, bool with_greeks) {
    OutputFormat out = parse_output(a.get_str("output", "table"));
    std::string inst = a.get_str("instrument", "european");
    if (is_rates_instrument(inst)) {
        cmd_price_rates(a, out);
        return;
    }
    EquityTrade t = EquityTrade::from_args(a);
    Report rep;
    rep.add("instrument", t.instrument);
    rep.add("type", to_string(t.type));
    rep.add("model", t.model);
    rep.add("method", t.resolved_method());
    rep.add("spot", t.S, 4);
    rep.add("strike", t.K, 4);
    rep.add("expiry_years", t.T, 6);

    if (with_greeks) {
        // First-class Heston risk for the semi-analytic European: real spot,
        // vega (parallel variance shift) and v0/theta/xi/rho sensitivities,
        // rather than the generic vol-bump path (#14).
        if (t.model == "heston" && t.instrument == "european" &&
            t.resolved_method() == "analytic") {
            HestonGreeks hg =
                heston_greeks(t.type, t.S, t.K, t.T, t.r, t.q, t.heston);
            rep.add("price", hg.price, 6);
            rep.add("delta", hg.delta, 6);
            rep.add("gamma", hg.gamma, 6);
            rep.add("vega", hg.vega / 100.0, 6);    // per vol point
            rep.add("theta", hg.theta / 365.0, 6);  // per calendar day
            rep.add("rho", hg.rho / 100.0, 6);      // per 1% rate
            rep.add("dV_dv0", hg.dv0, 6);
            rep.add("dV_dtheta", hg.dtheta, 6);
            rep.add("dV_dxi", hg.dxi, 6);
            rep.add("dV_drho", hg.drho, 6);
            rep.print(out, "Opal | risk report");
            if (out == OutputFormat::Table) {
                std::cout << "  (vega per vol pt, theta per day, rho per 1% "
                             "rate)\n";
                std::cout << "  (dV_dv0/dV_dtheta per unit variance; "
                             "dV_dxi/dV_drho per unit)\n";
            }
            return;
        }
        Greeks g = trade_greeks(t);
        rep.add("price", g.price, 6);
        rep.add("delta", g.delta, 6);
        rep.add("gamma", g.gamma, 6);
        rep.add("vega", g.vega / 100.0, 6);       // per vol point
        rep.add("theta", g.theta / 365.0, 6);     // per calendar day
        rep.add("rho", g.rho / 100.0, 6);         // per rate bp*100
        if (!std::isnan(g.vanna)) rep.add("vanna", g.vanna, 6);
        if (!std::isnan(g.volga)) rep.add("volga", g.volga, 6);
        if (!std::isnan(g.charm)) rep.add("charm", g.charm, 6);
        rep.print(out, "Opal | risk report");
        if (out == OutputFormat::Table)
            std::cout << "  (vega per vol pt, theta per day, rho per 1% rate)\n";
    } else {
        rep.add("price", t.price(), 6);
        double se = t.mc_std_error();
        if (!std::isnan(se)) rep.add("mc_std_error", se, 6);
        rep.print(out, "Opal | price");
    }
}

// ---------------------------------------------------------------------------
// implied vol
// ---------------------------------------------------------------------------
void cmd_implied(const Args& a) {
    OutputFormat out = parse_output(a.get_str("output", "table"));
    double price = a.require_num("price");
    double S = a.require_num("spot");
    double K = a.get_num("strike", S);
    double T = a.get_expiry();
    double r = a.get_num("rate", 0.0);
    double q = a.get_num("div", 0.0);
    std::string ty = a.get_str("type", "call");
    OptionType t = (ty == "put" || ty == "p") ? OptionType::Put : OptionType::Call;
    std::string model = a.get_str("model", "bsm");

    double vol;
    if (model == "bachelier")
        vol = implied_vol_bachelier(t, price, S, K, T, r);
    else if (model == "black76")
        vol = implied_vol_black76(t, price, S, K, T, r);
    else
        vol = implied_vol_bsm(t, price, S, K, T, r, q);

    Report rep;
    rep.add("model", model);
    rep.add("type", to_string(t));
    rep.add("market_price", price, 6);
    rep.add("implied_vol", vol, 6);
    if (model == "bsm") {
        Greeks g = bsm_greeks(t, S, K, T, r, q, vol);
        rep.add("vega", g.vega / 100.0, 6);
        rep.add("delta", g.delta, 6);
    }
    rep.print(out, "Opal | implied volatility");
}

// ---------------------------------------------------------------------------
// chain
// ---------------------------------------------------------------------------
void cmd_chain(const Args& a) {
    OutputFormat out = parse_output(a.get_str("output", "table"));
    EquityTrade base = EquityTrade::from_args(a);
    std::vector<double> strikes = a.get_range("strikes");
    Table tab({"strike", "call", "put", "call_delta", "put_delta", "gamma",
               "vega"});
    for (double k : strikes) {
        EquityTrade tc = base;
        tc.K = k;
        tc.type = OptionType::Call;
        EquityTrade tp = tc;
        tp.type = OptionType::Put;
        Greeks gc = trade_greeks(tc);
        Greeks gp = trade_greeks(tp);
        tab.add_row({fmt_num(k, 2), fmt_num(gc.price, 4), fmt_num(gp.price, 4),
                     fmt_num(gc.delta, 4), fmt_num(gp.delta, 4),
                     fmt_num(gc.gamma, 5), fmt_num(gc.vega / 100.0, 4)});
    }
    tab.print(out);
}

// ---------------------------------------------------------------------------
// scenario
// ---------------------------------------------------------------------------
void cmd_scenario(const Args& a) {
    OutputFormat out = parse_output(a.get_str("output", "table"));
    EquityTrade t = EquityTrade::from_args(a);
    std::vector<double> spot_pct =
        a.has("spot-range") ? a.get_range("spot-range")
                            : std::vector<double>{-20, -10, -5, 0, 5, 10, 20};
    std::vector<double> vol_pts =
        a.has("vol-range") ? a.get_range("vol-range")
                           : std::vector<double>{-5, 0, 5};
    double base = t.price();
    bool pnl = !a.flag("absolute");

    std::vector<std::string> headers{"spot\\vol"};
    for (double dv : vol_pts) headers.push_back(fmt_num(t.vol * 100 + dv, 1) + "%");
    Table tab(headers);
    for (double ds : spot_pct) {
        double s = t.S * (1.0 + ds / 100.0);
        std::vector<std::string> row{fmt_num(s, 2)};
        for (double dv : vol_pts) {
            double v = t.vol + dv / 100.0;
            if (v <= 0.0) {
                row.push_back("-");
                continue;
            }
            double p = t.price_at(s, v, t.r, t.T);
            row.push_back(fmt_num(pnl ? p - base : p, 4));
        }
        tab.add_row(row);
    }
    if (out == OutputFormat::Table) {
        std::cout << (pnl ? "P&L vs base price " : "Price grid, base ")
                  << fmt_num(base, 4) << " (rows: spot, cols: vol)\n";
    }
    tab.print(out);
}

// ---------------------------------------------------------------------------
// portfolio
// ---------------------------------------------------------------------------
std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) {
        // trim
        std::size_t b = cell.find_first_not_of(" \t\r");
        std::size_t e = cell.find_last_not_of(" \t\r");
        out.push_back(b == std::string::npos ? "" : cell.substr(b, e - b + 1));
    }
    return out;
}

void cmd_portfolio(const Args& a) {
    OutputFormat out = parse_output(a.get_str("output", "table"));
    std::string file = a.get_str("file");
    if (file.empty()) throw std::runtime_error("missing required option --file");
    std::ifstream in(file);
    if (!in) throw std::runtime_error("cannot open file '" + file + "'");

    std::string line;
    if (!std::getline(in, line))
        throw std::runtime_error("empty portfolio file");
    std::vector<std::string> header = split_csv(line);
    std::map<std::string, std::size_t> col;
    for (std::size_t i = 0; i < header.size(); ++i) col[header[i]] = i;
    for (const char* req : {"instrument", "type", "quantity", "spot", "strike",
                            "vol", "rate", "expiry"})
        if (!col.count(req))
            throw std::runtime_error(std::string("portfolio csv missing column '") +
                                     req + "'");

    Table tab({"#", "instrument", "type", "qty", "price", "value", "delta",
               "gamma", "vega", "theta"});
    double tot_value = 0, tot_delta = 0, tot_gamma = 0, tot_vega = 0,
           tot_theta = 0;
    int n = 0;
    auto cell = [&](const std::vector<std::string>& row, const std::string& c,
                    const std::string& def = "") {
        auto it = col.find(c);
        if (it == col.end() || it->second >= row.size() || row[it->second].empty())
            return def;
        return row[it->second];
    };

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> row = split_csv(line);
        ++n;
        // Numeric cells parsed with row/column context so a malformed value
        // points at the offending line instead of surfacing as a bare "stod".
        auto num = [&](const std::string& c, const std::string& def) {
            return parse_number(
                cell(row, c, def),
                "portfolio row " + std::to_string(n) + ", column '" + c + "'");
        };
        EquityTrade t;
        t.instrument = cell(row, "instrument", "european");
        std::string ty = cell(row, "type", "call");
        t.type = (ty == "put" || ty == "p") ? OptionType::Put : OptionType::Call;
        double qty = num("quantity", "0");
        t.S = num("spot", "0");
        t.K = num("strike", "0");
        t.vol = num("vol", "0");
        t.r = num("rate", "0");
        t.q = num("div", "0");
        // Expiry accepts a year fraction or a YYYY-MM-DD date, matching the
        // rest of the CLI.
        t.T = parse_expiry(
            cell(row, "expiry", "0"),
            "portfolio row " + std::to_string(n) + ", column 'expiry'");
        std::string bar = cell(row, "barrier");
        if (!bar.empty()) t.barrier = num("barrier", "0");
        std::string reb = cell(row, "rebate");
        if (!reb.empty()) t.rebate = num("rebate", "0");
        std::string meth = cell(row, "method");
        if (!meth.empty()) t.method = meth;

        Greeks g = trade_greeks(t);
        double value = qty * g.price;
        tot_value += value;
        tot_delta += qty * g.delta;
        tot_gamma += qty * g.gamma;
        tot_vega += qty * g.vega / 100.0;
        tot_theta += qty * g.theta / 365.0;
        tab.add_row({std::to_string(n), t.instrument, to_string(t.type),
                     fmt_num(qty, 0), fmt_num(g.price, 4), fmt_num(value, 2),
                     fmt_num(qty * g.delta, 2), fmt_num(qty * g.gamma, 4),
                     fmt_num(qty * g.vega / 100.0, 2),
                     fmt_num(qty * g.theta / 365.0, 2)});
    }
    tab.add_row({"", "TOTAL", "", "", "", fmt_num(tot_value, 2),
                 fmt_num(tot_delta, 2), fmt_num(tot_gamma, 4),
                 fmt_num(tot_vega, 2), fmt_num(tot_theta, 2)});
    tab.print(out);
    if (out == OutputFormat::Table)
        std::cout << "(vega per vol pt, theta per day; greeks share-equivalent)\n";
}

}  // namespace opal::cli

int main(int argc, char** argv) {
    using namespace opal::cli;
    try {
        Args args(argc, argv);
        const std::string& cmd = args.command();
        if (cmd.empty() || cmd == "help" || args.flag("help")) {
            std::cout << HELP;
            return 0;
        }
        if (cmd == "version" || args.flag("version")) {
            std::cout << "opal " << OPAL_VERSION_STRING << "\n";
            return 0;
        }
        if (cmd == "price")
            cmd_price(args, false);
        else if (cmd == "greeks")
            cmd_price(args, true);
        else if (cmd == "implied")
            cmd_implied(args);
        else if (cmd == "chain")
            cmd_chain(args);
        else if (cmd == "scenario")
            cmd_scenario(args);
        else if (cmd == "portfolio")
            cmd_portfolio(args);
        else {
            std::cerr << "opal: unknown command '" << cmd
                      << "' (try `opal help`)\n";
            return 2;
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "opal: " << e.what() << "\n";
        return 1;
    }
}
