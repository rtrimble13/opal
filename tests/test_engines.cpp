// Numerical engine validation: convergence to analytic values and
// cross-engine consistency.
#include "opal/opal.hpp"
#include "opal_test.hpp"

using namespace opal;

TEST_CASE(lattice_european_convergence) {
    double S = 100, K = 95, T = 0.75, r = 0.06, q = 0.02, sig = 0.3;
    double bs = bsm_price(OptionType::Call, S, K, T, r, q, sig);
    CHECK_CLOSE(binomial_crr_price(OptionType::Call, ExerciseStyle::European, S, K,
                                   T, r, q, sig, 2000),
                bs, 2e-3);
    CHECK_CLOSE(binomial_lr_price(OptionType::Call, ExerciseStyle::European, S, K,
                                  T, r, q, sig, 251),
                bs, 1e-4);
    CHECK_CLOSE(trinomial_price(OptionType::Call, ExerciseStyle::European, S, K, T,
                                r, q, sig, 800),
                bs, 2e-3);
}

TEST_CASE(lattice_american) {
    // American call on non-dividend stock equals European call.
    double S = 100, K = 100, T = 1.0, r = 0.06, sig = 0.25;
    double eu = bsm_price(OptionType::Call, S, K, T, r, 0.0, sig);
    double am = binomial_lr_price(OptionType::Call, ExerciseStyle::American, S, K,
                                  T, r, 0.0, sig, 501);
    CHECK_CLOSE(am, eu, 1e-4);
    // American put commands an early exercise premium.
    double eu_put = bsm_price(OptionType::Put, S, K, T, r, 0.0, sig);
    double am_put = binomial_lr_price(OptionType::Put, ExerciseStyle::American, S,
                                      K, T, r, 0.0, sig, 501);
    CHECK_TRUE(am_put > eu_put + 1e-4);
    // Cross-engine agreement on the American put.
    double am_crr = binomial_crr_price(OptionType::Put, ExerciseStyle::American, S,
                                       K, T, r, 0.0, sig, 2000);
    double am_tri = trinomial_price(OptionType::Put, ExerciseStyle::American, S, K,
                                    T, r, 0.0, sig, 1000);
    CHECK_CLOSE(am_crr, am_put, 5e-3);
    CHECK_CLOSE(am_tri, am_put, 5e-3);
    // Deep ITM American put exercises immediately: value = intrinsic.
    double deep = binomial_lr_price(OptionType::Put, ExerciseStyle::American, 40,
                                    100, 1.0, 0.08, 0.0, 0.2, 501);
    CHECK_CLOSE(deep, 60.0, 1e-9);
}

TEST_CASE(pde_european_and_american) {
    double S = 100, K = 100, T = 1.0, r = 0.05, q = 0.02, sig = 0.25;
    double bs = bsm_price(OptionType::Put, S, K, T, r, q, sig);
    PdeGrid g{600, 600, 6.0};
    CHECK_CLOSE(pde_price(OptionType::Put, ExerciseStyle::European, S, K, T, r, q,
                          sig, g),
                bs, 2e-3);
    CHECK_CLOSE(pde_price(OptionType::Call, ExerciseStyle::European, S, K, T, r, q,
                          sig, g),
                bsm_price(OptionType::Call, S, K, T, r, q, sig), 2e-3);
    // American put: PDE vs Leisen-Reimer.
    double am_lr = binomial_lr_price(OptionType::Put, ExerciseStyle::American, S,
                                     K, T, r, q, sig, 1001);
    double am_pde =
        pde_price(OptionType::Put, ExerciseStyle::American, S, K, T, r, q, sig, g);
    CHECK_CLOSE(am_pde, am_lr, 5e-3);
}

TEST_CASE(pde_barrier_vs_analytic) {
    double S = 100, K = 100, T = 0.5, r = 0.05, q = 0.0, sig = 0.25;
    PdeGrid g{800, 800, 6.0};
    double ana = barrier_price(OptionType::Call, BarrierType::DownOut, S, K, 90, T,
                               r, q, sig);
    double pde = pde_barrier_price(OptionType::Call, BarrierType::DownOut, S, K,
                                   90, T, r, q, sig, 0.0, g);
    CHECK_CLOSE(pde, ana, 5e-3);
    double ana_uo = barrier_price(OptionType::Put, BarrierType::UpOut, S, K, 115,
                                  T, r, q, sig);
    double pde_uo = pde_barrier_price(OptionType::Put, BarrierType::UpOut, S, K,
                                      115, T, r, q, sig, 0.0, g);
    CHECK_CLOSE(pde_uo, ana_uo, 5e-3);
    // Knock-in by parity.
    double ana_di = barrier_price(OptionType::Call, BarrierType::DownIn, S, K, 90,
                                  T, r, q, sig);
    double pde_di = pde_barrier_price(OptionType::Call, BarrierType::DownIn, S, K,
                                      90, T, r, q, sig, 0.0, g);
    CHECK_CLOSE(pde_di, ana_di, 5e-3);
}

TEST_CASE(monte_carlo_vanilla) {
    double S = 100, K = 105, T = 1.0, r = 0.05, q = 0.02, sig = 0.3;
    double bs = bsm_price(OptionType::Call, S, K, T, r, q, sig);
    McConfig cfg;
    cfg.paths = 60000;
    cfg.steps = 1;  // terminal value only; exact GBM step
    McResult res = mc_gbm(vanilla_payoff(OptionType::Call, K), S, T, r, q, sig, cfg);
    CHECK_CLOSE(res.price, bs, 4.0 * res.std_error + 1e-9);
    CHECK_TRUE(res.std_error < 0.1);
}

TEST_CASE(monte_carlo_asian_control_variate) {
    double S = 100, K = 100, T = 1.0, r = 0.05, q = 0.0, sig = 0.3;
    McConfig cfg;
    cfg.paths = 20000;
    cfg.steps = 64;
    // Geometric Asian MC must match its discrete closed form.
    double geo_cf = discrete_geometric_asian_price(OptionType::Call, S, K, T, r, q,
                                                   sig, cfg.steps);
    McResult geo_mc =
        mc_gbm(geometric_asian_payoff(OptionType::Call, K), S, T, r, q, sig, cfg);
    CHECK_CLOSE(geo_mc.price, geo_cf, 4.0 * geo_mc.std_error + 1e-9);
    // Arithmetic Asian with control variate: tight standard error, agrees
    // with Turnbull-Wakeman up to the discrete-vs-continuous averaging gap
    // and the TW approximation error.
    McResult arith = mc_arithmetic_asian(OptionType::Call, S, K, T, r, q, sig, cfg);
    CHECK_TRUE(arith.std_error < 0.01);
    double tw = arithmetic_asian_price(OptionType::Call, S, K, T, r, q, sig);
    CHECK_CLOSE(arith.price, tw, 0.10);
}

TEST_CASE(monte_carlo_lookback_and_barrier) {
    double S = 100, T = 0.5, r = 0.05, q = 0.0, sig = 0.3;
    McConfig cfg;
    cfg.paths = 40000;
    cfg.steps = 500;
    // Discretely monitored floating lookback converges (from below) to the
    // continuous closed form; allow a discretization gap.
    double cont = floating_lookback_price(OptionType::Call, S, T, r, q, sig);
    McResult mc = mc_gbm(lookback_payoff(OptionType::Call, StrikeStyle::Floating,
                                         0.0, S),
                         S, T, r, q, sig, cfg);
    CHECK_TRUE(mc.price < cont);
    CHECK_CLOSE(mc.price, cont, 0.05 * cont);
    // Discrete barrier knock-out is worth more than continuous (fewer chances
    // to knock out), and approaches it with dense monitoring.
    double cont_b = barrier_price(OptionType::Call, BarrierType::DownOut, S, 100,
                                  92, T, r, q, sig);
    McResult mcb = mc_gbm(
        barrier_payoff(OptionType::Call, BarrierType::DownOut, 100, 92, S), S, T,
        r, q, sig, cfg);
    CHECK_TRUE(mcb.price > cont_b);
    CHECK_CLOSE(mcb.price, cont_b, 0.08 * cont_b + 4.0 * mcb.std_error);
}

TEST_CASE(monte_carlo_barrier_rebate_convention) {
    // Knock-out rebates are paid at hit: passing r and T accrues the rebate
    // to expiry so the discounted MC value matches the analytic convention.
    double S = 100, K = 100, H = 92, T = 0.5, r = 0.08, q = 0.0, sig = 0.3,
           R = 5.0;
    McConfig cfg;
    cfg.paths = 40000;
    cfg.steps = 500;
    auto hit_paid = mc_gbm(barrier_payoff(OptionType::Call, BarrierType::DownOut,
                                          K, H, S, R, r, T),
                           S, T, r, q, sig, cfg);
    auto expiry_paid = mc_gbm(
        barrier_payoff(OptionType::Call, BarrierType::DownOut, K, H, S, R), S, T,
        r, q, sig, cfg);
    // Receiving the rebate earlier is strictly more valuable when r > 0.
    CHECK_TRUE(hit_paid.price > expiry_paid.price);
    // And the hit-time convention tracks the analytic price (up to the
    // discrete-monitoring gap).
    double ana = barrier_price(OptionType::Call, BarrierType::DownOut, S, K, H, T,
                               r, q, sig, R);
    CHECK_CLOSE(hit_paid.price, ana, 0.05 * ana + 4.0 * hit_paid.std_error);
}

TEST_CASE(monte_carlo_heston) {
    double S = 100, K = 100, T = 1.0, r = 0.03, q = 0.0;
    HestonParams hp{0.04, 1.5, 0.04, 0.3, -0.6};
    double semi = heston_price(OptionType::Call, S, K, T, r, q, hp);
    McConfig cfg;
    cfg.paths = 50000;
    cfg.steps = 200;
    cfg.antithetic = false;  // generic engine; antithetic invalid for 2-factor
    McResult mc = mc_heston(vanilla_payoff(OptionType::Call, K), S, T, r, q, hp, cfg);
    CHECK_CLOSE(mc.price, semi, 5.0 * mc.std_error + 0.05);
}

TEST_CASE(heston_greeks_mc_crn) {
    // Finite-difference Heston greeks over a Monte Carlo vanilla pricer under
    // common random numbers (a fixed seed) reproduce the analytic European
    // greeks: the CRN-stable first-order risks in magnitude, the parameter
    // sensitivities at least in sign. This exercises the engine-agnostic
    // heston_greeks_fd used for American/exotic Heston risk.
    double S = 100, K = 100, T = 1.0, r = 0.03, q = 0.0;
    HestonParams p{0.04, 1.5, 0.04, 0.3, -0.5};
    McConfig cfg;
    cfg.paths = 60000;
    cfg.steps = 64;
    cfg.seed = 7;
    cfg.antithetic = false;
    HestonPricerFn f = [&](double s, double tt, double rr,
                           const HestonParams& hp) {
        McConfig c = cfg;  // fixed seed across bumps -> common random numbers
        return mc_heston(vanilla_payoff(OptionType::Call, K), s, tt, rr, q, hp, c)
            .price;
    };
    // Mirror EquityTrade::mc_heston_bumps() so the shipped CLI configuration is
    // what gets validated: a wider spot bump for a usable gamma, defaults
    // elsewhere (common random numbers keep the first-order risks stable).
    HestonBumpSizes h;
    h.spot_rel = 1e-2;
    HestonGreeks g = heston_greeks_fd(f, S, T, r, p, h);
    HestonGreeks a = heston_greeks(OptionType::Call, S, K, T, r, q, p);
    CHECK_CLOSE(g.price, a.price, 0.3);
    CHECK_CLOSE(g.delta, a.delta, 0.05);
    CHECK_CLOSE(g.gamma, a.gamma, 5e-3);  // CRN gamma at the shipped bump
    CHECK_TRUE(g.vega > 0.0);
    CHECK_TRUE(g.dv0 > 0.0);
    CHECK_TRUE(g.dtheta > 0.0);
}

TEST_CASE(discrete_dividends) {
    double S = 100, K = 100, T = 1.0, r = 0.05, sig = 0.25;
    // Empty schedule reduces to plain BSM / plain LR tree.
    DividendSchedule none;
    CHECK_CLOSE(bsm_discrete_div_price(OptionType::Call, S, K, T, r, none, sig),
                bsm_price(OptionType::Call, S, K, T, r, 0.0, sig), 1e-12);
    CHECK_CLOSE(binomial_lr_discrete_div_price(OptionType::Put,
                                               ExerciseStyle::American, S, K, T,
                                               r, none, sig, 501),
                binomial_lr_price(OptionType::Put, ExerciseStyle::American, S, K,
                                  T, r, 0.0, sig, 501),
                1e-10);
    // European escrowed model: spot reduced by PV of dividends.
    DividendSchedule divs{{0.3, 2.0}, {0.8, 2.0}};
    double S_esc = S - 2.0 * std::exp(-r * 0.3) - 2.0 * std::exp(-r * 0.8);
    CHECK_CLOSE(bsm_discrete_div_price(OptionType::Call, S, K, T, r, divs, sig),
                bsm_price(OptionType::Call, S_esc, K, T, r, 0.0, sig), 1e-12);
    // Dividends cheapen calls and enrich puts.
    CHECK_TRUE(bsm_discrete_div_price(OptionType::Call, S, K, T, r, divs, sig) <
               bsm_price(OptionType::Call, S, K, T, r, 0.0, sig));
    CHECK_TRUE(bsm_discrete_div_price(OptionType::Put, S, K, T, r, divs, sig) >
               bsm_price(OptionType::Put, S, K, T, r, 0.0, sig));
    // A large dividend makes early exercise of the American call valuable.
    DividendSchedule big{{0.5, 8.0}};
    double am_call = binomial_lr_discrete_div_price(
        OptionType::Call, ExerciseStyle::American, S, K, T, r, big, sig, 501);
    double eu_call =
        bsm_discrete_div_price(OptionType::Call, S, K, T, r, big, sig);
    CHECK_TRUE(am_call > eu_call + 0.05);
    // American put premium survives dividends.
    double am_put = binomial_lr_discrete_div_price(
        OptionType::Put, ExerciseStyle::American, S, K, T, r, divs, sig, 501);
    double eu_put = bsm_discrete_div_price(OptionType::Put, S, K, T, r, divs, sig);
    CHECK_TRUE(am_put > eu_put);
    // PV of dividends exceeding spot must throw.
    DividendSchedule absurd{{0.5, 150.0}};
    CHECK_THROWS(bsm_discrete_div_price(OptionType::Call, S, K, T, r, absurd, sig));
}

TEST_CASE(lsmc_american_gbm) {
    double S = 100, K = 100, T = 1.0, r = 0.06, q = 0.0, sig = 0.25;
    LsmcConfig cfg;
    cfg.paths = 60000;
    cfg.steps = 50;
    // American put: LSMC tracks the Leisen-Reimer benchmark (LSMC is a
    // slightly low-biased estimator).
    double lr = binomial_lr_price(OptionType::Put, ExerciseStyle::American, S, K,
                                  T, r, q, sig, 1001);
    McResult lsmc = lsmc_american(OptionType::Put, S, K, T, r, q, sig, cfg);
    CHECK_TRUE(lsmc.price < lr + 4.0 * lsmc.std_error);
    CHECK_CLOSE(lsmc.price, lr, 0.10);
    // American call without dividends: no early exercise, equals European.
    double eu_call = bsm_price(OptionType::Call, S, K, T, r, q, sig);
    McResult lsmc_call = lsmc_american(OptionType::Call, S, K, T, r, q, sig, cfg);
    CHECK_CLOSE(lsmc_call.price, eu_call, 4.0 * lsmc_call.std_error + 0.10);
    // Deep ITM exercise floor.
    McResult deep = lsmc_american(OptionType::Put, 40, 100, T, 0.08, 0.0, 0.2, cfg);
    CHECK_TRUE(deep.price >= 60.0 - 1e-9);
}

TEST_CASE(lsmc_american_heston) {
    double S = 100, K = 100, T = 1.0, r = 0.05, q = 0.0;
    LsmcConfig cfg;
    cfg.paths = 40000;
    cfg.steps = 50;
    // Degenerate Heston (xi ~ 0, v0 = theta) matches GBM LSMC at sigma=sqrt(v0).
    HestonParams flat{0.0625, 2.0, 0.0625, 1e-4, 0.0};
    McResult hes = lsmc_american_heston(OptionType::Put, S, K, T, r, q, flat, cfg);
    double lr = binomial_lr_price(OptionType::Put, ExerciseStyle::American, S, K,
                                  T, r, q, 0.25, 1001);
    CHECK_CLOSE(hes.price, lr, 0.15);
    // Full Heston: American put carries an early exercise premium over the
    // semi-analytic European price.
    HestonParams hp{0.04, 1.5, 0.05, 0.5, -0.7};
    double eu = heston_price(OptionType::Put, S, K, T, r, q, hp);
    McResult am = lsmc_american_heston(OptionType::Put, S, K, T, r, q, hp, cfg);
    CHECK_TRUE(am.price > eu);
}

TEST_CASE(numerical_greeks_match_analytic) {
    double S = 100, K = 100, T = 1.0, r = 0.05, q = 0.02, sig = 0.25;
    Greeks ana = bsm_greeks(OptionType::Call, S, K, T, r, q, sig);
    auto pricer = [&](double s, double t, double rr, double v) {
        return bsm_price(OptionType::Call, s, K, t, rr, q, v);
    };
    Greeks num = numerical_greeks(pricer, S, T, r, sig);
    CHECK_CLOSE(num.delta, ana.delta, 1e-6);
    CHECK_CLOSE(num.gamma, ana.gamma, 1e-4);
    CHECK_CLOSE(num.vega, ana.vega, 1e-4);
    CHECK_CLOSE(num.rho, ana.rho, 1e-4);
    // One-day theta approximates the instantaneous theta.
    CHECK_CLOSE(num.theta, ana.theta, 0.05);
}
