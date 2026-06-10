"""Python API smoke and consistency tests (runnable with pytest or directly)."""
import glob
import math
import os
import sys

# Prefer the in-tree package only when the extension has been built in place
# (python setup.py build_ext --inplace); otherwise use the installed package.
_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
if glob.glob(os.path.join(_src, "opal", "_opal*.so")) or glob.glob(
        os.path.join(_src, "opal", "_opal*.pyd")):
    sys.path.insert(0, _src)

import opal


def approx(a, b, tol=1e-9):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def test_bs_price_and_parity():
    c = opal.bs_price("call", spot=42, strike=40, expiry=0.5, rate=0.10, vol=0.20)
    approx(c, 4.759422392871532, 1e-9)
    p = opal.bs_price("put", spot=42, strike=40, expiry=0.5, rate=0.10, vol=0.20)
    approx(c - p, 42 - 40 * math.exp(-0.05), 1e-10)


def test_greeks():
    g = opal.bs_greeks("call", spot=100, strike=100, expiry=1.0, rate=0.05,
                       div=0.02, vol=0.25)
    assert 0.5 < g.delta < 0.7
    assert g.gamma > 0 and g.vega > 0 and g.theta < 0


def test_implied_vol_roundtrip():
    c = opal.bs_price("call", spot=100, strike=110, expiry=0.7, rate=0.03,
                      div=0.01, vol=0.31)
    iv = opal.implied_vol("call", price=c, spot=100, strike=110, expiry=0.7,
                          rate=0.03, div=0.01)
    approx(iv, 0.31, 1e-7)


def test_exotics():
    # barrier in + out = vanilla
    common = dict(spot=100, strike=100, expiry=1.0, rate=0.05, vol=0.25)
    vin = opal.barrier_price("call", "down-in", barrier=90, **common)
    vout = opal.barrier_price("call", "down-out", barrier=90, **common)
    vanilla = opal.bs_price("call", **common)
    approx(vin + vout, vanilla, 1e-9)
    # asian < vanilla
    asian = opal.asian_price("call", **common)
    assert asian < vanilla
    # lookback > vanilla
    lb = opal.lookback_price("call", spot=100, expiry=1.0, rate=0.05, vol=0.25,
                             strike_style="floating")
    assert lb > vanilla
    # digital parity
    dc = opal.digital_price("call", **common)
    dp = opal.digital_price("put", **common)
    approx(dc + dp, math.exp(-0.05), 1e-10)


def test_engines_agree():
    common = dict(spot=100, strike=105, expiry=0.75, rate=0.04, div=0.01, vol=0.3)
    bs = opal.bs_price("call", **common)
    approx(opal.binomial_price("call", "european", steps=501, **common), bs, 1e-3)
    approx(opal.trinomial_price("call", "european", steps=600, **common), bs, 3e-3)
    approx(opal.pde_price("call", "european", grid=500, **common), bs, 3e-3)
    mc = opal.mc_price("call", paths=50000, steps=1, **common)
    assert abs(mc.price - bs) < 4 * mc.std_error + 1e-9
    # american put premium
    am = opal.binomial_price("put", "american", **common)
    eu = opal.bs_price("put", **common)
    assert am > eu


def test_mc_custom_payoff():
    res = opal.mc_custom(lambda path: max(path[-1] - 100.0, 0.0), spot=100,
                         expiry=1.0, rate=0.05, vol=0.2, paths=20000, steps=4)
    bs = opal.bs_price("call", spot=100, strike=100, expiry=1.0, rate=0.05,
                       vol=0.2)
    assert abs(res.price - bs) < 5 * res.std_error + 0.01


def test_heston():
    p = opal.HestonParams(v0=0.04, kappa=2.0, theta=0.04, xi=1e-4, rho=0.0)
    h = opal.heston_price("call", spot=100, strike=105, expiry=1.0, rate=0.04,
                          div=0.01, params=p)
    bs = opal.bs_price("call", spot=100, strike=105, expiry=1.0, rate=0.04,
                       div=0.01, vol=0.2)
    approx(h, bs, 1e-4)


def test_sabr():
    p = opal.SabrParams(alpha=0.25, beta=1.0, rho=0.0, nu=1e-8)
    approx(opal.sabr_vol(100, 90, 2.0, p), 0.25, 1e-6)


def test_rates():
    curve = opal.DiscountCurve(0.04)
    approx(curve.discount(2.0), math.exp(-0.08), 1e-12)
    cap = opal.cap_floor_price(curve, strike=0.04, vol=0.25, first_fixing=0.25,
                               maturity=3.0, tau=0.25, is_cap=True)
    assert cap.price > 0 and len(cap.caplets) == 11
    sw = opal.swaption_price(curve, "payer", strike=0.04, vol=0.3, expiry=1.0,
                             tenor=5.0)
    assert sw.price > 0 and sw.annuity > 0
    hw = opal.HullWhiteParams(a=0.1, sigma=0.01)
    c = opal.hw_zcb_option("call", curve, hw, 1.0, 3.0, 0.9)
    p = opal.hw_zcb_option("put", curve, hw, 1.0, 3.0, 0.9)
    approx(c - p, curve.discount(3.0) - 0.9 * curve.discount(1.0), 1e-12)
    # zero curve
    zc = opal.DiscountCurve([1.0, 2.0, 5.0], [0.03, 0.035, 0.04])
    assert zc.discount(2.0) < zc.discount(1.0)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} python API tests passed")
