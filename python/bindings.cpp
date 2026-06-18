// pybind11 bindings exposing the opal library to Python.
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "opal/opal.hpp"

namespace py = pybind11;
using namespace opal;

namespace {

DividendSchedule parse_dividends(
    const std::vector<std::pair<double, double>>& divs) {
    DividendSchedule out;
    out.reserve(divs.size());
    for (const auto& [t, amt] : divs) out.push_back({t, amt});
    return out;
}

OptionType parse_type(const std::string& s) {
    if (s == "call" || s == "c") return OptionType::Call;
    if (s == "put" || s == "p") return OptionType::Put;
    throw std::invalid_argument("option type must be 'call' or 'put', got '" + s + "'");
}

ExerciseStyle parse_style(const std::string& s) {
    if (s == "european") return ExerciseStyle::European;
    if (s == "american") return ExerciseStyle::American;
    throw std::invalid_argument("exercise must be 'european' or 'american'");
}

BarrierType parse_barrier(const std::string& s) {
    if (s == "down-in") return BarrierType::DownIn;
    if (s == "down-out") return BarrierType::DownOut;
    if (s == "up-in") return BarrierType::UpIn;
    if (s == "up-out") return BarrierType::UpOut;
    throw std::invalid_argument(
        "barrier type must be one of down-in, down-out, up-in, up-out");
}

}  // namespace

PYBIND11_MODULE(_opal, m) {
    m.doc() = "Opal: institutional option pricing library (C++ core)";
    m.attr("__version__") = OPAL_VERSION_STRING;

    // ----- result types ---------------------------------------------------
    py::class_<Greeks>(m, "Greeks")
        .def_readonly("price", &Greeks::price)
        .def_readonly("delta", &Greeks::delta)
        .def_readonly("gamma", &Greeks::gamma)
        .def_readonly("vega", &Greeks::vega)
        .def_readonly("theta", &Greeks::theta)
        .def_readonly("rho", &Greeks::rho)
        .def_readonly("vanna", &Greeks::vanna)
        .def_readonly("volga", &Greeks::volga)
        .def_readonly("charm", &Greeks::charm)
        .def("__repr__", [](const Greeks& g) {
            return "Greeks(price=" + std::to_string(g.price) +
                   ", delta=" + std::to_string(g.delta) +
                   ", gamma=" + std::to_string(g.gamma) +
                   ", vega=" + std::to_string(g.vega) +
                   ", theta=" + std::to_string(g.theta) +
                   ", rho=" + std::to_string(g.rho) + ")";
        });

    py::class_<HestonGreeks>(m, "HestonGreeks")
        .def_readonly("price", &HestonGreeks::price)
        .def_readonly("delta", &HestonGreeks::delta)
        .def_readonly("gamma", &HestonGreeks::gamma)
        .def_readonly("vega", &HestonGreeks::vega)
        .def_readonly("theta", &HestonGreeks::theta)
        .def_readonly("rho", &HestonGreeks::rho)
        .def_readonly("dv0", &HestonGreeks::dv0)
        .def_readonly("dtheta", &HestonGreeks::dtheta)
        .def_readonly("dxi", &HestonGreeks::dxi)
        .def_readonly("drho", &HestonGreeks::drho)
        .def("__repr__", [](const HestonGreeks& g) {
            return "HestonGreeks(price=" + std::to_string(g.price) +
                   ", delta=" + std::to_string(g.delta) +
                   ", gamma=" + std::to_string(g.gamma) +
                   ", vega=" + std::to_string(g.vega) +
                   ", theta=" + std::to_string(g.theta) +
                   ", rho=" + std::to_string(g.rho) +
                   ", dv0=" + std::to_string(g.dv0) +
                   ", dtheta=" + std::to_string(g.dtheta) +
                   ", dxi=" + std::to_string(g.dxi) +
                   ", drho=" + std::to_string(g.drho) + ")";
        });

    py::class_<McResult>(m, "McResult")
        .def_readonly("price", &McResult::price)
        .def_readonly("std_error", &McResult::std_error)
        .def_readonly("paths", &McResult::paths)
        .def("__repr__", [](const McResult& r) {
            return "McResult(price=" + std::to_string(r.price) +
                   ", std_error=" + std::to_string(r.std_error) + ")";
        });

    py::class_<HestonParams>(m, "HestonParams")
        .def(py::init([](double v0, double kappa, double theta, double xi,
                         double rho) {
                 return HestonParams{v0, kappa, theta, xi, rho};
             }),
             py::arg("v0"), py::arg("kappa"), py::arg("theta"), py::arg("xi"),
             py::arg("rho"))
        .def_readwrite("v0", &HestonParams::v0)
        .def_readwrite("kappa", &HestonParams::kappa)
        .def_readwrite("theta", &HestonParams::theta)
        .def_readwrite("xi", &HestonParams::xi)
        .def_readwrite("rho", &HestonParams::rho)
        .def("feller_satisfied", &HestonParams::feller_satisfied);

    py::class_<SabrParams>(m, "SabrParams")
        .def(py::init([](double alpha, double beta, double rho, double nu) {
                 return SabrParams{alpha, beta, rho, nu};
             }),
             py::arg("alpha"), py::arg("beta"), py::arg("rho"), py::arg("nu"))
        .def_readwrite("alpha", &SabrParams::alpha)
        .def_readwrite("beta", &SabrParams::beta)
        .def_readwrite("rho", &SabrParams::rho)
        .def_readwrite("nu", &SabrParams::nu);

    py::class_<HullWhiteParams>(m, "HullWhiteParams")
        .def(py::init([](double a, double sigma) {
                 return HullWhiteParams{a, sigma};
             }),
             py::arg("a"), py::arg("sigma"))
        .def_readwrite("a", &HullWhiteParams::a)
        .def_readwrite("sigma", &HullWhiteParams::sigma);

    py::class_<DiscountCurve>(m, "DiscountCurve")
        .def(py::init<double>(), py::arg("rate"))
        .def(py::init<std::vector<double>, std::vector<double>>(),
             py::arg("times"), py::arg("zero_rates"))
        .def("discount", &DiscountCurve::discount, py::arg("t"))
        .def("zero_rate", &DiscountCurve::zero_rate, py::arg("t"))
        .def("forward_rate", &DiscountCurve::forward_rate, py::arg("t1"),
             py::arg("t2"));

    py::class_<SwaptionResult>(m, "SwaptionResult")
        .def_readonly("price", &SwaptionResult::price)
        .def_readonly("forward_swap_rate", &SwaptionResult::forward_swap_rate)
        .def_readonly("annuity", &SwaptionResult::annuity);

    py::class_<CapletDetail>(m, "CapletDetail")
        .def_readonly("fixing", &CapletDetail::fixing)
        .def_readonly("payment", &CapletDetail::payment)
        .def_readonly("forward", &CapletDetail::forward)
        .def_readonly("price", &CapletDetail::price);

    py::class_<CapFloorResult>(m, "CapFloorResult")
        .def_readonly("price", &CapFloorResult::price)
        .def_readonly("caplets", &CapFloorResult::caplets);

    // ----- vanilla analytics ----------------------------------------------
    m.def(
        "bs_price",
        [](const std::string& t, double S, double K, double T, double r,
           double vol, double q) {
            return bsm_price(parse_type(t), S, K, T, r, q, vol);
        },
        py::arg("option_type"), py::arg("spot"), py::arg("strike"),
        py::arg("expiry"), py::arg("rate"), py::arg("vol"), py::arg("div") = 0.0,
        "Black-Scholes-Merton price with continuous dividend yield.");

    m.def(
        "bs_greeks",
        [](const std::string& t, double S, double K, double T, double r,
           double vol, double q) {
            return bsm_greeks(parse_type(t), S, K, T, r, q, vol);
        },
        py::arg("option_type"), py::arg("spot"), py::arg("strike"),
        py::arg("expiry"), py::arg("rate"), py::arg("vol"), py::arg("div") = 0.0,
        "Full analytic greeks (vega/theta/rho per unit; scale as needed).");

    m.def(
        "black76_price",
        [](const std::string& t, double F, double K, double T, double r,
           double vol) { return black76_price(parse_type(t), F, K, T, r, vol); },
        py::arg("option_type"), py::arg("forward"), py::arg("strike"),
        py::arg("expiry"), py::arg("rate"), py::arg("vol"));

    m.def(
        "bachelier_price",
        [](const std::string& t, double F, double K, double T, double r,
           double vol) { return bachelier_price(parse_type(t), F, K, T, r, vol); },
        py::arg("option_type"), py::arg("forward"), py::arg("strike"),
        py::arg("expiry"), py::arg("rate"), py::arg("normal_vol"));

    m.def(
        "implied_vol",
        [](const std::string& t, double price, double S, double K, double T,
           double r, double q) {
            return implied_vol_bsm(parse_type(t), price, S, K, T, r, q);
        },
        py::arg("option_type"), py::arg("price"), py::arg("spot"),
        py::arg("strike"), py::arg("expiry"), py::arg("rate"),
        py::arg("div") = 0.0, "Implied BSM volatility from a market price.");

    m.def(
        "implied_vol_bachelier",
        [](const std::string& t, double price, double F, double K, double T,
           double r) {
            return implied_vol_bachelier(parse_type(t), price, F, K, T, r);
        },
        py::arg("option_type"), py::arg("price"), py::arg("forward"),
        py::arg("strike"), py::arg("expiry"), py::arg("rate"));

    // ----- digitals / exotics (analytic) ------------------------------------
    m.def(
        "digital_price",
        [](const std::string& t, double S, double K, double T, double r,
           double vol, double q, double cash, bool asset) {
            return asset ? asset_or_nothing_price(parse_type(t), S, K, T, r, q, vol)
                         : cash_or_nothing_price(parse_type(t), S, K, T, r, q, vol,
                                                 cash);
        },
        py::arg("option_type"), py::arg("spot"), py::arg("strike"),
        py::arg("expiry"), py::arg("rate"), py::arg("vol"), py::arg("div") = 0.0,
        py::arg("cash") = 1.0, py::arg("asset") = false);

    m.def(
        "barrier_price",
        [](const std::string& t, const std::string& bt, double S, double K,
           double H, double T, double r, double vol, double q, double rebate) {
            return barrier_price(parse_type(t), parse_barrier(bt), S, K, H, T, r,
                                 q, vol, rebate);
        },
        py::arg("option_type"), py::arg("barrier_type"), py::arg("spot"),
        py::arg("strike"), py::arg("barrier"), py::arg("expiry"), py::arg("rate"),
        py::arg("vol"), py::arg("div") = 0.0, py::arg("rebate") = 0.0,
        "Reiner-Rubinstein continuously monitored single barrier.");

    m.def(
        "asian_price",
        [](const std::string& t, double S, double K, double T, double r,
           double vol, double q, const std::string& avg, int fixings) {
            if (avg == "geometric") {
                if (fixings > 0)
                    return discrete_geometric_asian_price(parse_type(t), S, K, T,
                                                          r, q, vol, fixings);
                return geometric_asian_price(parse_type(t), S, K, T, r, q, vol);
            }
            return arithmetic_asian_price(parse_type(t), S, K, T, r, q, vol);
        },
        py::arg("option_type"), py::arg("spot"), py::arg("strike"),
        py::arg("expiry"), py::arg("rate"), py::arg("vol"), py::arg("div") = 0.0,
        py::arg("average") = "arithmetic", py::arg("fixings") = 0,
        "Asian option: geometric exact, arithmetic via Turnbull-Wakeman.");

    m.def(
        "lookback_price",
        [](const std::string& t, double S, double T, double r, double vol,
           double K, double q, const std::string& strike_style, double extremum) {
            if (strike_style == "floating")
                return floating_lookback_price(parse_type(t), S, T, r, q, vol,
                                               extremum);
            return fixed_lookback_price(parse_type(t), S, K, T, r, q, vol,
                                        extremum);
        },
        py::arg("option_type"), py::arg("spot"), py::arg("expiry"),
        py::arg("rate"), py::arg("vol"), py::arg("strike") = 0.0,
        py::arg("div") = 0.0, py::arg("strike_style") = "floating",
        py::arg("extremum") = 0.0);

    // ----- numerical engines ------------------------------------------------
    m.def(
        "binomial_price",
        [](const std::string& t, const std::string& style, double S, double K,
           double T, double r, double vol, double q, int steps,
           const std::string& flavor) {
            if (flavor == "crr")
                return binomial_crr_price(parse_type(t), parse_style(style), S, K,
                                          T, r, q, vol, steps);
            return binomial_lr_price(parse_type(t), parse_style(style), S, K, T,
                                     r, q, vol, steps);
        },
        py::arg("option_type"), py::arg("exercise"), py::arg("spot"),
        py::arg("strike"), py::arg("expiry"), py::arg("rate"), py::arg("vol"),
        py::arg("div") = 0.0, py::arg("steps") = 501, py::arg("flavor") = "lr",
        "Binomial tree (Leisen-Reimer default, or flavor='crr').");

    m.def(
        "trinomial_price",
        [](const std::string& t, const std::string& style, double S, double K,
           double T, double r, double vol, double q, int steps) {
            return trinomial_price(parse_type(t), parse_style(style), S, K, T, r,
                                   q, vol, steps);
        },
        py::arg("option_type"), py::arg("exercise"), py::arg("spot"),
        py::arg("strike"), py::arg("expiry"), py::arg("rate"), py::arg("vol"),
        py::arg("div") = 0.0, py::arg("steps") = 400);

    m.def(
        "pde_price",
        [](const std::string& t, const std::string& style, double S, double K,
           double T, double r, double vol, double q, int grid) {
            return pde_price(parse_type(t), parse_style(style), S, K, T, r, q,
                             vol, PdeGrid{grid, grid, 6.0});
        },
        py::arg("option_type"), py::arg("exercise"), py::arg("spot"),
        py::arg("strike"), py::arg("expiry"), py::arg("rate"), py::arg("vol"),
        py::arg("div") = 0.0, py::arg("grid") = 400,
        "Crank-Nicolson finite difference price.");

    m.def(
        "mc_price",
        [](const std::string& t, double S, double K, double T, double r,
           double vol, double q, const std::string& payoff, std::size_t paths,
           int steps, std::uint64_t seed, double barrier,
           const std::string& barrier_type, double cash, double rebate) {
            McConfig cfg;
            cfg.paths = paths;
            cfg.steps = steps;
            cfg.seed = seed;
            OptionType ot = parse_type(t);
            // asian-arith uses the control-variate engine, not a plain payoff.
            if (payoff == "asian-arith")
                return mc_arithmetic_asian(ot, S, K, T, r, q, vol, cfg);
            // All other payoffs share the library's path-payoff factory so the
            // CLI and Python bindings stay in lock-step (incl. barrier rebate
            // accrual via r/T).
            if (payoff == "vanilla" || payoff == "digital") cfg.steps = 1;
            PayoffParams pp;
            pp.type = ot;
            pp.K = K;
            pp.cash = cash;
            pp.S0 = S;
            pp.barrier = barrier;
            pp.rebate = rebate;
            pp.r = r;
            pp.T = T;
            if (payoff == "barrier") pp.barrier_type = parse_barrier(barrier_type);
            return mc_gbm(make_path_payoff(payoff, pp), S, T, r, q, vol, cfg);
        },
        py::arg("option_type"), py::arg("spot"), py::arg("strike"),
        py::arg("expiry"), py::arg("rate"), py::arg("vol"), py::arg("div") = 0.0,
        py::arg("payoff") = "vanilla", py::arg("paths") = 100000,
        py::arg("steps") = 252, py::arg("seed") = 42, py::arg("barrier") = 0.0,
        py::arg("barrier_type") = "down-out", py::arg("cash") = 1.0,
        py::arg("rebate") = 0.0,
        "Monte Carlo with antithetic variates; control variate for "
        "asian-arith. Returns McResult(price, std_error).");

    m.def(
        "mc_custom",
        [](std::function<double(const std::vector<double>&)> payoff, double S,
           double T, double r, double vol, double q, std::size_t paths, int steps,
           std::uint64_t seed, bool antithetic) {
            McConfig cfg;
            cfg.paths = paths;
            cfg.steps = steps;
            cfg.seed = seed;
            cfg.antithetic = antithetic;
            return mc_gbm(payoff, S, T, r, q, vol, cfg);
        },
        py::arg("payoff"), py::arg("spot"), py::arg("expiry"), py::arg("rate"),
        py::arg("vol"), py::arg("div") = 0.0, py::arg("paths") = 50000,
        py::arg("steps") = 252, py::arg("seed") = 42, py::arg("antithetic") = true,
        "Monte Carlo on a user-supplied Python path payoff f(path)->float. "
        "The path excludes S0 and the payoff is settled (undiscounted) at "
        "expiry. Slower than the built-in payoffs due to Python callbacks.");

    // ----- discrete cash dividends -------------------------------------------
    m.def(
        "bs_discrete_div_price",
        [](const std::string& t, double S, double K, double T, double r,
           double vol, const std::vector<std::pair<double, double>>& dividends,
           double q) {
            return bsm_discrete_div_price(parse_type(t), S, K, T, r,
                                          parse_dividends(dividends), vol, q);
        },
        py::arg("option_type"), py::arg("spot"), py::arg("strike"),
        py::arg("expiry"), py::arg("rate"), py::arg("vol"),
        py::arg("dividends"), py::arg("div") = 0.0,
        "European option with discrete cash dividends [(time, amount), ...] "
        "via the escrowed-spot BSM model.");

    m.def(
        "binomial_discrete_div_price",
        [](const std::string& t, const std::string& style, double S, double K,
           double T, double r, double vol,
           const std::vector<std::pair<double, double>>& dividends, int steps,
           double q) {
            return binomial_lr_discrete_div_price(
                parse_type(t), parse_style(style), S, K, T, r,
                parse_dividends(dividends), vol, steps, q);
        },
        py::arg("option_type"), py::arg("exercise"), py::arg("spot"),
        py::arg("strike"), py::arg("expiry"), py::arg("rate"), py::arg("vol"),
        py::arg("dividends"), py::arg("steps") = 501, py::arg("div") = 0.0,
        "American/European option with discrete cash dividends on a "
        "Leisen-Reimer tree (escrowed spot, dividend add-back at exercise).");

    // ----- Longstaff-Schwartz American Monte Carlo ----------------------------
    m.def(
        "lsmc_price",
        [](const std::string& t, double S, double K, double T, double r,
           double vol, double q, std::size_t paths, int steps,
           std::uint64_t seed) {
            LsmcConfig cfg;
            cfg.paths = paths;
            cfg.steps = steps;
            cfg.seed = seed;
            return lsmc_american(parse_type(t), S, K, T, r, q, vol, cfg);
        },
        py::arg("option_type"), py::arg("spot"), py::arg("strike"),
        py::arg("expiry"), py::arg("rate"), py::arg("vol"), py::arg("div") = 0.0,
        py::arg("paths") = 50000, py::arg("steps") = 50, py::arg("seed") = 42,
        "American option under GBM via Longstaff-Schwartz least-squares "
        "Monte Carlo. Returns McResult (slightly low-biased estimator).");

    m.def(
        "lsmc_heston_price",
        [](const std::string& t, double S, double K, double T, double r,
           double q, const HestonParams& p, std::size_t paths, int steps,
           std::uint64_t seed) {
            LsmcConfig cfg;
            cfg.paths = paths;
            cfg.steps = steps;
            cfg.seed = seed;
            return lsmc_american_heston(parse_type(t), S, K, T, r, q, p, cfg);
        },
        py::arg("option_type"), py::arg("spot"), py::arg("strike"),
        py::arg("expiry"), py::arg("rate"), py::arg("div"), py::arg("params"),
        py::arg("paths") = 50000, py::arg("steps") = 50, py::arg("seed") = 42,
        "American option under Heston dynamics via Longstaff-Schwartz "
        "(variance enters the regression basis).");

    // ----- stochastic vol ---------------------------------------------------
    m.def(
        "heston_price",
        [](const std::string& t, double S, double K, double T, double r, double q,
           const HestonParams& p) {
            return heston_price(parse_type(t), S, K, T, r, q, p);
        },
        py::arg("option_type"), py::arg("spot"), py::arg("strike"),
        py::arg("expiry"), py::arg("rate"), py::arg("div"), py::arg("params"),
        "Semi-analytic Heston price (characteristic function).");

    m.def(
        "heston_greeks",
        [](const std::string& t, double S, double K, double T, double r, double q,
           const HestonParams& p) {
            return heston_greeks(parse_type(t), S, K, T, r, q, p);
        },
        py::arg("option_type"), py::arg("spot"), py::arg("strike"),
        py::arg("expiry"), py::arg("rate"), py::arg("div"), py::arg("params"),
        "First-class Heston greeks: spot delta/gamma, a parallel-variance-shift "
        "vega, theta, rate rho, plus v0/theta/xi/rho parameter sensitivities. "
        "Returns HestonGreeks.");

    m.def(
        "heston_mc",
        [](std::function<double(const std::vector<double>&)> payoff, double S,
           double T, double r, double q, const HestonParams& p,
           std::size_t paths, int steps, std::uint64_t seed) {
            McConfig cfg;
            cfg.paths = paths;
            cfg.steps = steps;
            cfg.seed = seed;
            cfg.antithetic = false;
            return mc_heston(payoff, S, T, r, q, p, cfg);
        },
        py::arg("payoff"), py::arg("spot"), py::arg("expiry"), py::arg("rate"),
        py::arg("div"), py::arg("params"), py::arg("paths") = 50000,
        py::arg("steps") = 252, py::arg("seed") = 42);

    m.def("sabr_vol", &sabr_lognormal_vol, py::arg("forward"), py::arg("strike"),
          py::arg("expiry"), py::arg("params"),
          "Hagan 2002 lognormal SABR implied vol.");
    m.def("sabr_normal_vol", &sabr_normal_vol, py::arg("forward"),
          py::arg("strike"), py::arg("expiry"), py::arg("params"));
    m.def(
        "sabr_price",
        [](const std::string& t, double F, double K, double T, double r,
           const SabrParams& p) {
            return sabr_price(parse_type(t), F, K, T, r, p);
        },
        py::arg("option_type"), py::arg("forward"), py::arg("strike"),
        py::arg("expiry"), py::arg("rate"), py::arg("params"));

    // ----- rates ------------------------------------------------------------
    m.def(
        "cap_floor_price",
        [](const DiscountCurve& curve, double K, double vol, double first_fixing,
           double maturity, double tau, bool is_cap, const std::string& vol_type) {
            return cap_floor_price(curve, K, vol, first_fixing, maturity, tau,
                                   is_cap,
                                   vol_type == "normal" ? RateVolType::Normal
                                                        : RateVolType::Lognormal);
        },
        py::arg("curve"), py::arg("strike"), py::arg("vol"),
        py::arg("first_fixing"), py::arg("maturity"), py::arg("tau") = 0.25,
        py::arg("is_cap") = true, py::arg("vol_type") = "lognormal");

    m.def(
        "cap_floor_price_ois",
        [](const DiscountCurve& discount_curve, const DiscountCurve& forward_curve,
           double K, double vol, double first_fixing, double maturity, double tau,
           bool is_cap, const std::string& vol_type) {
            return cap_floor_price(discount_curve, forward_curve, K, vol,
                                   first_fixing, maturity, tau, is_cap,
                                   vol_type == "normal" ? RateVolType::Normal
                                                        : RateVolType::Lognormal);
        },
        py::arg("discount_curve"), py::arg("forward_curve"), py::arg("strike"),
        py::arg("vol"), py::arg("first_fixing"), py::arg("maturity"),
        py::arg("tau") = 0.25, py::arg("is_cap") = true,
        py::arg("vol_type") = "lognormal",
        "Multi-curve cap/floor: forwards off forward_curve, cashflows "
        "discounted on discount_curve (OIS).");

    m.def(
        "swaption_price",
        [](const DiscountCurve& curve, const std::string& style, double K,
           double vol, double expiry, double tenor, double pay_freq,
           const std::string& vol_type) {
            SwaptionType st = (style == "receiver") ? SwaptionType::Receiver
                                                     : SwaptionType::Payer;
            return swaption_price(curve, st, K, vol, expiry, tenor, pay_freq,
                                  vol_type == "normal" ? RateVolType::Normal
                                                       : RateVolType::Lognormal);
        },
        py::arg("curve"), py::arg("style"), py::arg("strike"), py::arg("vol"),
        py::arg("expiry"), py::arg("tenor"), py::arg("pay_freq") = 2.0,
        py::arg("vol_type") = "lognormal");

    m.def(
        "swaption_price_ois",
        [](const DiscountCurve& discount_curve, const DiscountCurve& forward_curve,
           const std::string& style, double K, double vol, double expiry,
           double tenor, double pay_freq, const std::string& vol_type) {
            SwaptionType st = (style == "receiver") ? SwaptionType::Receiver
                                                     : SwaptionType::Payer;
            return swaption_price(discount_curve, forward_curve, st, K, vol,
                                  expiry, tenor, pay_freq,
                                  vol_type == "normal" ? RateVolType::Normal
                                                       : RateVolType::Lognormal);
        },
        py::arg("discount_curve"), py::arg("forward_curve"), py::arg("style"),
        py::arg("strike"), py::arg("vol"), py::arg("expiry"), py::arg("tenor"),
        py::arg("pay_freq") = 2.0, py::arg("vol_type") = "lognormal",
        "Multi-curve European swaption: OIS-discounted annuity, par swap "
        "rate from the projection curve.");

    m.def(
        "sabr_swaption_price",
        [](const DiscountCurve& curve, const std::string& style, double K,
           const SabrParams& p, double expiry, double tenor, double pay_freq) {
            SwaptionType st = (style == "receiver") ? SwaptionType::Receiver
                                                     : SwaptionType::Payer;
            return sabr_swaption_price(curve, st, K, p, expiry, tenor, pay_freq);
        },
        py::arg("curve"), py::arg("style"), py::arg("strike"), py::arg("params"),
        py::arg("expiry"), py::arg("tenor"), py::arg("pay_freq") = 2.0);

    m.def(
        "hw_zcb_option",
        [](const std::string& t, const DiscountCurve& curve,
           const HullWhiteParams& p, double T1, double T2, double K) {
            return hw_zcb_option_price(parse_type(t), curve, p, T1, T2, K);
        },
        py::arg("option_type"), py::arg("curve"), py::arg("params"),
        py::arg("option_expiry"), py::arg("bond_maturity"), py::arg("strike"));

    m.def("hw_caplet", &hw_caplet_price, py::arg("curve"), py::arg("params"),
          py::arg("fixing"), py::arg("payment"), py::arg("strike"));
    m.def("hw_floorlet", &hw_floorlet_price, py::arg("curve"), py::arg("params"),
          py::arg("fixing"), py::arg("payment"), py::arg("strike"));
    m.def("hw_cap", &hw_cap_price, py::arg("curve"), py::arg("params"),
          py::arg("first_fixing"), py::arg("maturity"), py::arg("tau"),
          py::arg("strike"), py::arg("is_cap") = true);
}
