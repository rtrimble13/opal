// Interest rate model and instrument validation.
#include "opal/opal.hpp"
#include "opal_test.hpp"

using namespace opal;

TEST_CASE(discount_curve) {
    DiscountCurve flat(0.05);
    CHECK_CLOSE(flat.discount(2.0), std::exp(-0.1), 1e-15);
    CHECK_CLOSE(flat.zero_rate(7.3), 0.05, 1e-15);
    // Simple forward over [1,2] on a flat 5% curve.
    CHECK_CLOSE(flat.forward_rate(1.0, 2.0), std::exp(0.05) - 1.0, 1e-12);

    DiscountCurve curve({1.0, 2.0, 5.0}, {0.03, 0.035, 0.04});
    CHECK_CLOSE(curve.zero_rate(1.0), 0.03, 1e-15);
    CHECK_CLOSE(curve.zero_rate(5.0), 0.04, 1e-15);
    CHECK_CLOSE(curve.zero_rate(10.0), 0.04, 1e-15);  // flat extrapolation
    // Log-linear DF interpolation: r(1.5)*1.5 = avg of pillar r*t.
    CHECK_CLOSE(curve.zero_rate(1.5) * 1.5, 0.5 * (0.03 + 0.07), 1e-12);
    // Discount factors decrease.
    CHECK_TRUE(curve.discount(2.0) < curve.discount(1.0));
}

TEST_CASE(cap_floor_black) {
    DiscountCurve curve(0.04);
    double K = 0.04, vol = 0.25;
    auto cap = cap_floor_price(curve, K, vol, 0.25, 3.0, 0.25, true);
    auto floor = cap_floor_price(curve, K, vol, 0.25, 3.0, 0.25, false);
    CHECK_TRUE(cap.price > 0.0);
    CHECK_TRUE(floor.price > 0.0);
    CHECK_TRUE(cap.caplets.size() == 11);
    // Cap - floor = swap of forwards vs K (put-call parity per caplet).
    double swap_pv = 0.0;
    for (auto& c : cap.caplets)
        swap_pv += (c.payment - c.fixing) * curve.discount(c.payment) *
                   (c.forward - K);
    CHECK_CLOSE(cap.price - floor.price, swap_pv, 1e-10);
    // Vol monotonicity.
    auto cap_hi = cap_floor_price(curve, K, vol + 0.1, 0.25, 3.0, 0.25, true);
    CHECK_TRUE(cap_hi.price > cap.price);
    // Normal vol variant prices sensibly (same order of magnitude).
    auto cap_n = cap_floor_price(curve, K, 0.01, 0.25, 3.0, 0.25, true,
                                 RateVolType::Normal);
    CHECK_TRUE(cap_n.price > 0.0);
}

TEST_CASE(swaption_black) {
    DiscountCurve curve(0.04);
    auto atm = detail::swap_metrics(curve, 1.0, 5.0, 2.0);
    double K = atm.forward_swap_rate;
    auto payer = swaption_price(curve, SwaptionType::Payer, K, 0.3, 1.0, 5.0, 2.0);
    auto recv =
        swaption_price(curve, SwaptionType::Receiver, K, 0.3, 1.0, 5.0, 2.0);
    // ATM payer = ATM receiver (parity: payer - receiver = annuity (F - K)).
    CHECK_CLOSE(payer.price, recv.price, 1e-12);
    CHECK_TRUE(payer.price > 0.0);
    CHECK_TRUE(payer.annuity > 0.0);
    // OTM payer < ATM payer.
    auto otm = swaption_price(curve, SwaptionType::Payer, K + 0.01, 0.3, 1.0, 5.0,
                              2.0);
    CHECK_TRUE(otm.price < payer.price);
    // SABR swaption with flat params (beta=1, nu~0) matches Black.
    SabrParams flat{0.3, 1.0, 0.0, 1e-8};
    auto sabr = sabr_swaption_price(curve, SwaptionType::Payer, K, flat, 1.0, 5.0,
                                    2.0);
    CHECK_CLOSE(sabr.price, payer.price, 1e-6);
}

TEST_CASE(multi_curve_ois) {
    DiscountCurve ois(0.035);
    DiscountCurve proj(0.040);  // projection curve above OIS (positive basis)
    double K = 0.04, vol = 0.25;

    // Dual-curve with identical curves collapses to single-curve exactly.
    auto single = cap_floor_price(proj, K, vol, 0.25, 3.0, 0.25, true);
    auto dual_same = cap_floor_price(proj, proj, K, vol, 0.25, 3.0, 0.25, true);
    CHECK_CLOSE(dual_same.price, single.price, 1e-15);
    auto sw_single = swaption_price(proj, SwaptionType::Payer, K, 0.3, 1.0, 5.0);
    auto sw_dual_same =
        swaption_price(proj, proj, SwaptionType::Payer, K, 0.3, 1.0, 5.0);
    CHECK_CLOSE(sw_dual_same.price, sw_single.price, 1e-15);
    CHECK_CLOSE(sw_dual_same.forward_swap_rate, sw_single.forward_swap_rate, 1e-12);

    // OIS discounting with higher projected forwards: same forwards as the
    // projection-only run, discounted on the (higher-DF) OIS curve -> caplets
    // are each worth more.
    auto dual = cap_floor_price(ois, proj, K, vol, 0.25, 3.0, 0.25, true);
    CHECK_TRUE(dual.price > single.price);
    for (std::size_t i = 0; i < dual.caplets.size(); ++i)
        CHECK_CLOSE(dual.caplets[i].forward, single.caplets[i].forward, 1e-12);

    // Multi-curve swaption: par rate reflects the projection curve, annuity
    // the OIS curve; payer/receiver parity holds at the multi-curve par rate.
    auto m = detail::swap_metrics(ois, proj, 1.0, 5.0, 2.0);
    auto m_proj = detail::swap_metrics(proj, 1.0, 5.0, 2.0);
    CHECK_CLOSE(m.forward_swap_rate, m_proj.forward_swap_rate, 5e-4);
    CHECK_TRUE(m.annuity > m_proj.annuity);  // OIS discounts less
    auto payer = swaption_price(ois, proj, SwaptionType::Payer,
                                m.forward_swap_rate, 0.3, 1.0, 5.0);
    auto recv = swaption_price(ois, proj, SwaptionType::Receiver,
                               m.forward_swap_rate, 0.3, 1.0, 5.0);
    CHECK_CLOSE(payer.price, recv.price, 1e-12);
}

TEST_CASE(hull_white) {
    DiscountCurve curve(0.05);
    HullWhiteParams p{0.1, 0.01};
    // ZCB option put-call parity: c - p = P(0,T2) - K P(0,T1).
    double K = 0.9, T1 = 1.0, T2 = 3.0;
    double c = hw_zcb_option_price(OptionType::Call, curve, p, T1, T2, K);
    double pp = hw_zcb_option_price(OptionType::Put, curve, p, T1, T2, K);
    CHECK_CLOSE(c - pp, curve.discount(T2) - K * curve.discount(T1), 1e-12);
    // Caplet/floorlet parity: caplet - floorlet = tau * df * (F - K).
    double Kr = 0.05, tau = 0.5;
    double cap = hw_caplet_price(curve, p, 1.0, 1.5, Kr);
    double flo = hw_floorlet_price(curve, p, 1.0, 1.5, Kr);
    double F = curve.forward_rate(1.0, 1.5);
    CHECK_CLOSE(cap - flo, tau * curve.discount(1.5) * (F - Kr), 1e-12);
    // Vol monotonicity.
    HullWhiteParams hi{0.1, 0.02};
    CHECK_TRUE(hw_caplet_price(curve, hi, 1.0, 1.5, Kr) > cap);
    // Full cap is the sum of its caplets and is positive.
    double capfull = hw_cap_price(curve, p, 0.5, 3.0, 0.5, Kr);
    CHECK_TRUE(capfull > 0.0);
}

TEST_CASE(cashflow_schedule_period_count) {
    DiscountCurve curve(0.03);
    // Monthly cap, tau = 1/12 over ~2y: exactly 23 caplets, with no
    // floating-point drift in the schedule (#10).
    auto monthly = cap_floor_price(curve, 0.03, 0.2, 1.0 / 12.0, 2.0, 1.0 / 12.0,
                                   true);
    CHECK_TRUE(monthly.caplets.size() == 23u);
    // Power-of-two tau is unchanged: quarterly from 0.25 to 3.0 is 11 caplets.
    auto quarterly = cap_floor_price(curve, 0.04, 0.25, 0.25, 3.0, 0.25, true);
    CHECK_TRUE(quarterly.caplets.size() == 11u);
}
