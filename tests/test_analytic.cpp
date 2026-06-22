// Validation of analytic models against textbook reference values and
// no-arbitrage identities.
#include <random>

#include "opal/opal.hpp"
#include "opal_test.hpp"

using namespace opal;

TEST_CASE(normal_distribution) {
    CHECK_CLOSE(math::norm_cdf(0.0), 0.5, 1e-15);
    CHECK_CLOSE(math::norm_cdf(1.96), 0.9750021048517795, 1e-12);
    CHECK_CLOSE(math::norm_ppf(0.975), 1.959963984540054, 1e-9);
    CHECK_CLOSE(math::norm_ppf(math::norm_cdf(1.2345)), 1.2345, 1e-10);
    CHECK_CLOSE(math::norm_pdf(0.0), 0.3989422804014327, 1e-15);
}

TEST_CASE(black_scholes_hull_reference) {
    // Hull, "Options, Futures and Other Derivatives": S=42, K=40, r=10%,
    // sigma=20%, T=0.5 -> call 4.76, put 0.81.
    double c = bsm_price(OptionType::Call, 42, 40, 0.5, 0.10, 0.0, 0.20);
    double p = bsm_price(OptionType::Put, 42, 40, 0.5, 0.10, 0.0, 0.20);
    CHECK_CLOSE(c, 4.759422392871532, 1e-9);
    CHECK_CLOSE(p, 0.8085993729000922, 1e-9);
}

TEST_CASE(put_call_parity) {
    double S = 100, K = 95, T = 1.25, r = 0.04, q = 0.015, sig = 0.3;
    double c = bsm_price(OptionType::Call, S, K, T, r, q, sig);
    double p = bsm_price(OptionType::Put, S, K, T, r, q, sig);
    double parity = S * std::exp(-q * T) - K * std::exp(-r * T);
    CHECK_CLOSE(c - p, parity, 1e-12);
}

TEST_CASE(black_scholes_greeks) {
    double S = 100, K = 100, T = 1.0, r = 0.05, q = 0.02, sig = 0.25;
    Greeks g = bsm_greeks(OptionType::Call, S, K, T, r, q, sig);
    // Compare with central finite differences of the analytic price.
    double h = 1e-5;
    auto price = [&](double s, double t, double rr, double v) {
        return bsm_price(OptionType::Call, s, K, t, rr, q, v);
    };
    CHECK_CLOSE(g.delta, (price(S + h, T, r, sig) - price(S - h, T, r, sig)) / (2 * h), 1e-6);
    CHECK_CLOSE(g.gamma,
                (price(S + h, T, r, sig) - 2 * g.price + price(S - h, T, r, sig)) / (h * h),
                1e-4);
    CHECK_CLOSE(g.vega, (price(S, T, r, sig + h) - price(S, T, r, sig - h)) / (2 * h), 1e-5);
    CHECK_CLOSE(g.rho, (price(S, T, r + h, sig) - price(S, T, r - h, sig)) / (2 * h), 1e-5);
    CHECK_CLOSE(g.theta, -(price(S, T + h, r, sig) - price(S, T - h, r, sig)) / (2 * h), 1e-5);
    // Vanna and volga vs finite differences.
    double vanna_fd = (bsm_greeks(OptionType::Call, S, K, T, r, q, sig + h).delta -
                       bsm_greeks(OptionType::Call, S, K, T, r, q, sig - h).delta) /
                      (2 * h);
    CHECK_CLOSE(g.vanna, vanna_fd, 1e-5);
    double volga_fd = (bsm_greeks(OptionType::Call, S, K, T, r, q, sig + h).vega -
                       bsm_greeks(OptionType::Call, S, K, T, r, q, sig - h).vega) /
                      (2 * h);
    CHECK_CLOSE(g.volga, volga_fd, 1e-4);
    double charm_fd = -(bsm_greeks(OptionType::Call, S, K, T + h, r, q, sig).delta -
                        bsm_greeks(OptionType::Call, S, K, T - h, r, q, sig).delta) /
                      (2 * h);
    CHECK_CLOSE(g.charm, charm_fd, 1e-5);
}

TEST_CASE(black76_and_bachelier) {
    // Black-76 equals BSM with q = r (b = 0).
    double F = 105, K = 100, T = 0.75, r = 0.03;
    CHECK_CLOSE(black76_price(OptionType::Call, F, K, T, r, 0.2),
                bsm_price(OptionType::Call, F, K, T, r, r, 0.2), 1e-12);
    // Bachelier ATM call = df * sigma * sqrt(T) / sqrt(2 pi).
    double atm = bachelier_price(OptionType::Call, 100, 100, 1.0, 0.0, 10.0);
    CHECK_CLOSE(atm, 10.0 * 0.3989422804014327, 1e-10);
    // Bachelier put-call parity: c - p = df (F - K).
    double c = bachelier_price(OptionType::Call, 102, 100, 0.5, 0.04, 8.0);
    double p = bachelier_price(OptionType::Put, 102, 100, 0.5, 0.04, 8.0);
    CHECK_CLOSE(c - p, std::exp(-0.02) * 2.0, 1e-12);
}

TEST_CASE(digitals) {
    double S = 100, K = 100, T = 1, r = 0.05, q = 0.02, sig = 0.25;
    // cash call + cash put = discount factor
    double cc = cash_or_nothing_price(OptionType::Call, S, K, T, r, q, sig);
    double cp = cash_or_nothing_price(OptionType::Put, S, K, T, r, q, sig);
    CHECK_CLOSE(cc + cp, std::exp(-r * T), 1e-12);
    // asset call + asset put = forward PV
    double ac = asset_or_nothing_price(OptionType::Call, S, K, T, r, q, sig);
    double ap = asset_or_nothing_price(OptionType::Put, S, K, T, r, q, sig);
    CHECK_CLOSE(ac + ap, S * std::exp(-q * T), 1e-12);
    // vanilla = asset-or-nothing - K * cash-or-nothing
    double vanilla = bsm_price(OptionType::Call, S, K, T, r, q, sig);
    CHECK_CLOSE(vanilla, ac - K * cc, 1e-12);
    // gap option with equal strikes = vanilla
    CHECK_CLOSE(gap_option_price(OptionType::Call, S, K, K, T, r, q, sig), vanilla,
                1e-12);
}

TEST_CASE(barrier_haug_reference) {
    // Haug, "Complete Guide to Option Pricing Formulas", standard barrier
    // table: S=100, T=0.5, r=0.08, b=0.04 (q=0.04), sigma=0.25, rebate=3.
    double S = 100, T = 0.5, r = 0.08, q = 0.04, sig = 0.25, R = 3;
    CHECK_CLOSE(barrier_price(OptionType::Call, BarrierType::DownOut, S, 90, 95, T,
                              r, q, sig, R),
                9.0246, 5e-4);
    CHECK_CLOSE(barrier_price(OptionType::Call, BarrierType::DownOut, S, 100, 95,
                              T, r, q, sig, R),
                6.7924, 5e-4);
    CHECK_CLOSE(barrier_price(OptionType::Call, BarrierType::UpOut, S, 90, 105, T,
                              r, q, sig, R),
                2.6789, 5e-4);
    CHECK_CLOSE(barrier_price(OptionType::Call, BarrierType::DownIn, S, 90, 95, T,
                              r, q, sig, R),
                7.7627, 5e-4);
    // Verified against dense-monitoring Monte Carlo (400k paths, 20k steps:
    // 1.4536; continuous-limit analytic is slightly above the discrete MC).
    CHECK_CLOSE(barrier_price(OptionType::Put, BarrierType::UpIn, S, 90, 105, T, r,
                              q, sig, R),
                1.46531, 5e-4);
}

TEST_CASE(barrier_in_out_parity) {
    // With zero rebate: knock-in + knock-out = vanilla, for all 8 combos.
    double S = 100, K = 102, T = 0.8, r = 0.05, q = 0.01, sig = 0.3;
    for (OptionType t : {OptionType::Call, OptionType::Put}) {
        double vanilla = bsm_price(t, S, K, T, r, q, sig);
        for (double H : {90.0, 110.0}) {
            BarrierType bin = H < S ? BarrierType::DownIn : BarrierType::UpIn;
            BarrierType bout = H < S ? BarrierType::DownOut : BarrierType::UpOut;
            double vin = barrier_price(t, bin, S, K, H, T, r, q, sig);
            double vout = barrier_price(t, bout, S, K, H, T, r, q, sig);
            CHECK_CLOSE(vin + vout, vanilla, 1e-10);
        }
    }
}

TEST_CASE(asian_analytic) {
    double S = 100, K = 100, T = 1, r = 0.05, q = 0.0, sig = 0.3;
    double geo = geometric_asian_price(OptionType::Call, S, K, T, r, q, sig);
    double arith = arithmetic_asian_price(OptionType::Call, S, K, T, r, q, sig);
    double vanilla = bsm_price(OptionType::Call, S, K, T, r, q, sig);
    // Averaging reduces volatility: asian < vanilla; arithmetic > geometric.
    CHECK_TRUE(geo < vanilla);
    CHECK_TRUE(arith > geo);
    CHECK_TRUE(arith < vanilla);
    // Discrete geometric converges to continuous geometric.
    double disc = discrete_geometric_asian_price(OptionType::Call, S, K, T, r, q,
                                                 sig, 10000);
    CHECK_CLOSE(disc, geo, 2e-3);
}

TEST_CASE(lookback_analytic) {
    double S = 100, T = 0.5, r = 0.1, q = 0.0, sig = 0.3;
    // Floating lookback call >= vanilla ATM call (extra optionality).
    double flb = floating_lookback_price(OptionType::Call, S, T, r, q, sig);
    double vanilla = bsm_price(OptionType::Call, S, S, T, r, q, sig);
    CHECK_TRUE(flb > vanilla);
    // Fixed lookback call with K = S0 >= vanilla.
    double fxd = fixed_lookback_price(OptionType::Call, S, S, T, r, q, sig);
    CHECK_TRUE(fxd > vanilla);
    // Seasoned floating-strike lookback call: S=120, running min=100, T=0.5,
    // r=0.1, q=0, sigma=0.3. Verified against dense-monitoring Monte Carlo
    // (300k paths, 20k steps: 28.175; continuous analytic slightly above).
    CHECK_CLOSE(floating_lookback_price(OptionType::Call, 120, 0.5, 0.1, 0.0, 0.3,
                                        100.0),
                28.2133, 5e-4);
}

TEST_CASE(implied_vol_round_trip) {
    double S = 100, K = 110, T = 0.7, r = 0.03, q = 0.01;
    for (double sig : {0.08, 0.2, 0.55, 1.2}) {
        double c = bsm_price(OptionType::Call, S, K, T, r, q, sig);
        CHECK_CLOSE(implied_vol_bsm(OptionType::Call, c, S, K, T, r, q), sig, 1e-8);
        double p = bsm_price(OptionType::Put, S, K, T, r, q, sig);
        CHECK_CLOSE(implied_vol_bsm(OptionType::Put, p, S, K, T, r, q), sig, 1e-8);
    }
    // Bachelier round trip.
    double bp = bachelier_price(OptionType::Call, 0.03, 0.035, 1.0, 0.02, 0.008);
    CHECK_CLOSE(implied_vol_bachelier(OptionType::Call, bp, 0.03, 0.035, 1.0, 0.02),
                0.008, 1e-9);
    // Arbitrage-violating price must throw.
    CHECK_THROWS(implied_vol_bsm(OptionType::Call, 200.0, S, K, T, r, q));
}

TEST_CASE(heston_bs_limit_and_parity) {
    // With xi tiny and v0 = theta, Heston degenerates to BS at sigma=sqrt(v0).
    HestonParams p{0.04, 2.0, 0.04, 1e-4, 0.0};
    double S = 100, K = 105, T = 1.0, r = 0.04, q = 0.01;
    double h = heston_price(OptionType::Call, S, K, T, r, q, p);
    double bs = bsm_price(OptionType::Call, S, K, T, r, q, 0.2);
    CHECK_CLOSE(h, bs, 1e-4);
    // Put-call parity under full Heston.
    HestonParams p2{0.05, 1.5, 0.06, 0.4, -0.7};
    double c = heston_price(OptionType::Call, S, K, T, r, q, p2);
    double pp = heston_price(OptionType::Put, S, K, T, r, q, p2);
    CHECK_CLOSE(c - pp, S * std::exp(-q * T) - K * std::exp(-r * T), 1e-8);
}

TEST_CASE(heston_short_dated_high_vol_of_vol) {
    // Regression for issue #8: the inverse-Fourier integral must extend its
    // upper limit for short maturities / high vol-of-vol, where the integrand
    // (decaying like exp(-v0*T*phi^2/2)) is still O(1e-2) past phi = 200. The
    // former hard cap at 200 over-counted the tail and mispriced OTM calls.
    // Reference values confirmed against an extended-truncation quadrature and
    // 4M-path Monte Carlo.
    double S = 100, r = 0.02, q = 0.0;
    HestonParams p{0.04, 1.5, 0.04, 3.0, -0.7};  // xi = 3.0 (strong Feller break)

    // 36-day OTM call: true ~= 0.00841 (MC 0.00859 +/- 0.00015); the bug gave
    // 0.0242, ~3x too high.
    double c1 = heston_price(OptionType::Call, S, 120.0, 0.10, r, q, p);
    CHECK_CLOSE(c1, 0.008408, 2e-4);
    // 3.65-day deep-OTM call: true ~= 0 (MC 0.0); the bug gave 0.0256.
    double c2 = heston_price(OptionType::Call, S, 120.0, 0.01, r, q, p);
    CHECK_CLOSE(c2, 0.0, 5e-4);
    // 36-day OTM at K=110: true ~= 0.06147 (MC 0.06388 +/- 0.0004); bug: 0.0749.
    double c3 = heston_price(OptionType::Call, S, 110.0, 0.10, r, q, p);
    CHECK_CLOSE(c3, 0.061465, 1e-3);
    // Put-call parity must still hold exactly in this regime.
    double p3 = heston_price(OptionType::Put, S, 110.0, 0.10, r, q, p);
    CHECK_CLOSE(c3 - p3, S * std::exp(-q * 0.10) - 110.0 * std::exp(-r * 0.10),
                1e-9);
}

TEST_CASE(heston_greeks_bs_limit_and_signs) {
    double S = 100, K = 105, T = 1.0, r = 0.04, q = 0.01;
    // BS limit: xi tiny, v0 = theta = 0.04 (sigma = 0.2). First-class Heston
    // spot/vega greeks must collapse onto the analytic BSM greeks.
    HestonParams flat{0.04, 2.0, 0.04, 1e-4, 0.0};
    HestonGreeks hg = heston_greeks(OptionType::Call, S, K, T, r, q, flat);
    Greeks bs = bsm_greeks(OptionType::Call, S, K, T, r, q, 0.2);
    CHECK_CLOSE(hg.price, bs.price, 1e-3);
    CHECK_CLOSE(hg.delta, bs.delta, 1e-4);
    CHECK_CLOSE(hg.gamma, bs.gamma, 1e-3);
    CHECK_CLOSE(hg.vega, bs.vega, 2e-2);   // per 1.00 vol
    CHECK_CLOSE(hg.rho, bs.rho, 1e-2);     // per 1.00 rate
    // Decay convention: a call without offsetting carry loses value over time.
    CHECK_TRUE(hg.theta < 0.0);

    // Full Heston: parameter sensitivities have the expected signs. More
    // variance (v0/theta) and a vega bump all enrich a vanilla; the call's
    // rate rho is positive.
    HestonParams p{0.05, 1.5, 0.06, 0.4, -0.7};
    HestonGreeks fg = heston_greeks(OptionType::Call, S, K, T, r, q, p);
    CHECK_TRUE(fg.dv0 > 0.0);
    CHECK_TRUE(fg.dtheta > 0.0);
    CHECK_TRUE(fg.vega > 0.0);
    CHECK_TRUE(fg.rho > 0.0);
    CHECK_TRUE(fg.delta > 0.0 && fg.delta < 1.0);
    CHECK_TRUE(fg.gamma > 0.0);

    // Spot-delta put-call relation: dC/dS - dP/dS = exp(-qT).
    HestonGreeks pg = heston_greeks(OptionType::Put, S, K, T, r, q, p);
    CHECK_CLOSE(fg.delta - pg.delta, std::exp(-q * T), 1e-5);
    // Parameter sensitivities are spot-direction independent (same dv0/dxi/drho
    // for a call and put struck identically, by put-call parity).
    CHECK_CLOSE(fg.dv0, pg.dv0, 1e-4);
    CHECK_CLOSE(fg.dxi, pg.dxi, 1e-4);
    CHECK_CLOSE(fg.drho, pg.drho, 1e-4);
}

TEST_CASE(sabr_hagan) {
    // beta=1, nu->0 collapses to constant lognormal vol alpha.
    SabrParams flat{0.25, 1.0, 0.0, 1e-8};
    CHECK_CLOSE(sabr_lognormal_vol(100, 80, 2.0, flat), 0.25, 1e-6);
    CHECK_CLOSE(sabr_lognormal_vol(100, 100, 2.0, flat), 0.25, 1e-6);
    // Negative rho tilts the skew: lower strikes have higher vol.
    SabrParams skew{0.2, 0.5, -0.3, 0.4};
    double v_low = sabr_lognormal_vol(0.03, 0.02, 1.0, skew);
    double v_high = sabr_lognormal_vol(0.03, 0.04, 1.0, skew);
    CHECK_TRUE(v_low > v_high);
    // ATM formula continuity: K -> F limit matches ATM branch.
    double atm = sabr_lognormal_vol(0.03, 0.03, 1.0, skew);
    double near = sabr_lognormal_vol(0.03, 0.0300001, 1.0, skew);
    CHECK_CLOSE(near, atm, 1e-5);
}

TEST_CASE(two_asset_rainbow_and_exchange) {
    double S1 = 100, S2 = 95, K = 100, T = 1.0, r = 0.05, q1 = 0.02, q2 = 0.04,
           sig1 = 0.25, sig2 = 0.30, rho = 0.4;

    // Exact identity: for every state {max,min} = {S1,S2} as a set, so
    // (max-K)+ + (min-K)+ = (S1-K)+ + (S2-K)+. Hence call-on-max + call-on-min
    // equals the two vanillas; likewise for puts.
    double cmax = option_on_max_price(OptionType::Call, S1, S2, K, T, r, q1, q2,
                                      sig1, sig2, rho);
    double cmin = option_on_min_price(OptionType::Call, S1, S2, K, T, r, q1, q2,
                                      sig1, sig2, rho);
    double v1 = bsm_price(OptionType::Call, S1, K, T, r, q1, sig1);
    double v2 = bsm_price(OptionType::Call, S2, K, T, r, q2, sig2);
    CHECK_CLOSE(cmax + cmin, v1 + v2, 1e-9);
    double pmax = option_on_max_price(OptionType::Put, S1, S2, K, T, r, q1, q2,
                                      sig1, sig2, rho);
    double pmin = option_on_min_price(OptionType::Put, S1, S2, K, T, r, q1, q2,
                                      sig1, sig2, rho);
    double pv1 = bsm_price(OptionType::Put, S1, K, T, r, q1, sig1);
    double pv2 = bsm_price(OptionType::Put, S2, K, T, r, q2, sig2);
    CHECK_CLOSE(pmax + pmin, pv1 + pv2, 1e-9);

    // Put-call parity on the maximum: Cmax - Pmax = e^{-rT}E[max] - K e^{-rT}.
    double exch = exchange_option_price(S1, S2, T, q1, q2, sig1, sig2, rho);
    double max_fwd = S2 * std::exp(-q2 * T) + exch;
    CHECK_CLOSE(cmax - pmax, max_fwd - K * std::exp(-r * T), 1e-9);

    // Exchange degenerates to a vanilla call: swap a deterministic asset 2
    // (sig2->0) that grows at r (q2=r) for asset 1 == call on S1 struck S2.
    double exch_lim = exchange_option_price(100, 100, T, 0.0, r, sig1, 1e-9, 0.0);
    CHECK_CLOSE(exch_lim, bsm_price(OptionType::Call, 100, 100, T, r, 0.0, sig1),
                1e-9);

    // Correlation option: with the condition strike K1 -> 0 the gate S1 > K1 is
    // always met, so it reduces to a vanilla call on asset 2.
    double corr_lim = two_asset_correlation_price(OptionType::Call, S1, S2, 1e-8, K,
                                                  T, r, q1, q2, sig1, sig2, rho);
    CHECK_CLOSE(corr_lim, bsm_price(OptionType::Call, S2, K, T, r, q2, sig2), 1e-6);
    // Correlation put with gate K1 -> 0: the condition S1 < K1 never holds, so
    // it is worthless.
    double pcorr_lim = two_asset_correlation_price(OptionType::Put, S1, S2, 1e-8, K,
                                                   T, r, q1, q2, sig1, sig2, rho);
    CHECK_CLOSE(pcorr_lim, 0.0, 1e-7);

    // Monte Carlo cross-check (fixed seed -> deterministic) of the absolute
    // levels for max/min/exchange and the correlation option.
    std::mt19937_64 rng(20260618);
    std::normal_distribution<double> nd(0.0, 1.0);
    const long N = 400000;
    double K1c = 90.0;  // correlation gate on S1
    double d1 = (r - q1 - 0.5 * sig1 * sig1) * T, d2 = (r - q2 - 0.5 * sig2 * sig2) * T;
    double w1 = sig1 * std::sqrt(T), w2 = sig2 * std::sqrt(T);
    double sMax = 0, sMin = 0, sExch = 0, sCorr = 0, sCorrP = 0;
    for (long n = 0; n < N; ++n) {
        double z1 = nd(rng);
        double z2 = rho * z1 + std::sqrt(1.0 - rho * rho) * nd(rng);
        double s1 = S1 * std::exp(d1 + w1 * z1);
        double s2 = S2 * std::exp(d2 + w2 * z2);
        sMax += std::max(std::max(s1, s2) - K, 0.0);
        sMin += std::max(std::min(s1, s2) - K, 0.0);
        sExch += std::max(s1 - s2, 0.0);
        if (s1 > K1c) sCorr += std::max(s2 - K, 0.0);
        if (s1 < K1c) sCorrP += std::max(K - s2, 0.0);
    }
    double disc = std::exp(-r * T);
    CHECK_CLOSE(disc * sMax / N, cmax, 0.08);
    CHECK_CLOSE(disc * sMin / N, cmin, 0.08);
    CHECK_CLOSE(disc * sExch / N, exch, 0.08);
    double corr = two_asset_correlation_price(OptionType::Call, S1, S2, K1c, K, T, r,
                                              q1, q2, sig1, sig2, rho);
    CHECK_CLOSE(disc * sCorr / N, corr, 0.08);
    double corr_put = two_asset_correlation_price(OptionType::Put, S1, S2, K1c, K, T,
                                                  r, q1, q2, sig1, sig2, rho);
    CHECK_CLOSE(disc * sCorrP / N, corr_put, 0.08);
}

TEST_CASE(compound_geske) {
    double S = 100, K1 = 6.0, K2 = 100, t1 = 0.5, T2 = 1.0, r = 0.05, q = 0.0,
           sig = 0.25;
    double b = r - q;

    // Compound put-call parity on the outer option (same inner call): holding
    // call-on-call and shorting put-on-call locks in the inner option for K1 at
    // t1, worth c - K1 e^{-r t1} today.
    double coc = compound_option_price(OptionType::Call, OptionType::Call, S, K1, K2,
                                       t1, T2, r, q, sig);
    double poc = compound_option_price(OptionType::Put, OptionType::Call, S, K1, K2,
                                       t1, T2, r, q, sig);
    double c = gbs_price(OptionType::Call, S, K2, T2, r, b, sig);
    CHECK_CLOSE(coc - poc, c - K1 * std::exp(-r * t1), 1e-9);
    // Same parity for inner puts.
    double cop = compound_option_price(OptionType::Call, OptionType::Put, S, K1, K2,
                                       t1, T2, r, q, sig);
    double pop = compound_option_price(OptionType::Put, OptionType::Put, S, K1, K2,
                                       t1, T2, r, q, sig);
    double p = gbs_price(OptionType::Put, S, K2, T2, r, b, sig);
    CHECK_CLOSE(cop - pop, p - K1 * std::exp(-r * t1), 1e-9);

    // K1 -> 0: a call-on-call is always exercised, so it equals the inner call.
    double coc0 = compound_option_price(OptionType::Call, OptionType::Call, S, 1e-8,
                                        K2, t1, T2, r, q, sig);
    CHECK_CLOSE(coc0, c, 1e-5);

    // Deeper-ITM inner put whose value at the spot exceeds K1: exercises the
    // root-bracket expansion (the naive [lo, spot] bracket would have failed).
    // Outer put-call parity must still hold exactly.
    double Sl = 80.0;  // inner put (K2=100) worth ~20 > K1 = 6
    double cop2 = compound_option_price(OptionType::Call, OptionType::Put, Sl, K1, K2,
                                        t1, T2, r, q, sig);
    double pop2 = compound_option_price(OptionType::Put, OptionType::Put, Sl, K1, K2,
                                        t1, T2, r, q, sig);
    double p2 = gbs_price(OptionType::Put, Sl, K2, T2, r, b, sig);
    CHECK_CLOSE(cop2 - pop2, p2 - K1 * std::exp(-r * t1), 1e-9);

    // Never-bracketed regime: K1 exceeds the inner put's maximum value
    // (K2 e^{-r tau}). The call-on-put is worthless; the put-on-put is always
    // exercised, worth K1 e^{-r t1} minus today's inner put.
    double bigK1 = 200.0;
    double cop3 = compound_option_price(OptionType::Call, OptionType::Put, S, bigK1,
                                        K2, t1, T2, r, q, sig);
    double pop3 = compound_option_price(OptionType::Put, OptionType::Put, S, bigK1,
                                        K2, t1, T2, r, q, sig);
    CHECK_CLOSE(cop3, 0.0, 1e-12);
    CHECK_CLOSE(pop3,
                bigK1 * std::exp(-r * t1) - gbs_price(OptionType::Put, S, K2, T2, r,
                                                      b, sig),
                1e-9);

    // Monte Carlo cross-check: value the inner option at t1, pay max(.-K1,0).
    std::mt19937_64 rng(99);
    std::normal_distribution<double> nd(0.0, 1.0);
    const long N = 200000;
    double drift = (r - q - 0.5 * sig * sig) * t1, w = sig * std::sqrt(t1);
    double sum = 0;
    for (long n = 0; n < N; ++n) {
        double s = S * std::exp(drift + w * nd(rng));
        double iv = gbs_price(OptionType::Call, s, K2, T2 - t1, r, b, sig);
        sum += std::max(iv - K1, 0.0);
    }
    CHECK_CLOSE(std::exp(-r * t1) * sum / N, coc, 0.05);
}

TEST_CASE(partial_time_start_barrier) {
    double S = 100, K = 100, T2 = 1.0, r = 0.05, q = 0.0, sig = 0.25, t1 = 0.5;
    double Hd = 90.0, Hu = 115.0;

    // Knock-in + knock-out over the same window/barrier == vanilla (exact).
    double van = bsm_price(OptionType::Call, S, K, T2, r, q, sig);
    double ko = partial_time_start_barrier_price(OptionType::Call,
                                                 PartialBarrierType::DownOut, S, K,
                                                 Hd, t1, T2, r, q, sig);
    double ki = partial_time_start_barrier_price(OptionType::Call,
                                                 PartialBarrierType::DownIn, S, K,
                                                 Hd, t1, T2, r, q, sig);
    CHECK_CLOSE(ko + ki, van, 1e-9);

    // A far barrier almost never knocks: KO ~ vanilla.
    double ko_far = partial_time_start_barrier_price(
        OptionType::Call, PartialBarrierType::DownOut, S, K, 50.0, t1, T2, r, q, sig);
    CHECK_CLOSE(ko_far, van, 1e-3);

    // Barrier at/through the spot at inception: KO already dead, KI == vanilla.
    CHECK_CLOSE(partial_time_start_barrier_price(OptionType::Call,
                                                 PartialBarrierType::DownOut, S, K,
                                                 100.0, t1, T2, r, q, sig),
                0.0, 1e-12);

    // Window-only monitoring knocks out less often than full-life continuous
    // monitoring, so the partial-time KO is worth more than the full barrier.
    double full_ko = barrier_price(OptionType::Call, BarrierType::DownOut, S, K, Hd,
                                   T2, r, q, sig);
    CHECK_TRUE(ko > full_ko);

    // Monte Carlo cross-check, barrier monitored only on [0, t1]. The model is
    // continuous-monitoring, so the MC uses a Brownian-bridge continuity
    // correction (per-step no-crossing probability) to converge to the
    // continuous price with coarse stepping instead of paying a discretization
    // premium.
    std::mt19937_64 rng(2024);
    std::normal_distribution<double> nd(0.0, 1.0);
    const long N = 300000;
    const int n1 = 100;
    double dt = t1 / n1, dr = (r - q - 0.5 * sig * sig) * dt, vo = sig * std::sqrt(dt);
    double tau = T2 - t1, drT = (r - q - 0.5 * sig * sig) * tau, voT = sig * std::sqrt(tau);
    double v2dt = sig * sig * dt, lhd = std::log(Hd), lhu = std::log(Hu);
    double sumDO = 0, sumUO = 0;
    for (long p = 0; p < N; ++p) {
        double x = std::log(S), wD = 1.0, wU = 1.0;
        for (int i = 0; i < n1; ++i) {
            double xp = x;
            x += dr + vo * nd(rng);
            if (xp <= lhd || x <= lhd) wD = 0.0;
            else wD *= 1.0 - std::exp(-2.0 * (xp - lhd) * (x - lhd) / v2dt);
            if (xp >= lhu || x >= lhu) wU = 0.0;
            else wU *= 1.0 - std::exp(-2.0 * (lhu - xp) * (lhu - x) / v2dt);
        }
        double pay = std::max(std::exp(x + drT + voT * nd(rng)) - K, 0.0);
        sumDO += wD * pay;
        sumUO += wU * pay;
    }
    double disc = std::exp(-r * T2);
    CHECK_CLOSE(disc * sumDO / N, ko, 0.05);
    double ko_up = partial_time_start_barrier_price(
        OptionType::Call, PartialBarrierType::UpOut, S, K, Hu, t1, T2, r, q, sig);
    CHECK_CLOSE(disc * sumUO / N, ko_up, 0.05);
}
