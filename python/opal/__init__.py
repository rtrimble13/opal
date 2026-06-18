"""Opal: institutional option pricing and risk analytics.

Python API over the Opal C++ core. Designed for scripting, pre-trade
analysis and Jupyter notebook workflows.

Quick start
-----------
>>> import opal
>>> opal.bs_price("call", spot=100, strike=105, expiry=0.5, rate=0.04, vol=0.22)
4.937...
>>> g = opal.bs_greeks("put", spot=100, strike=100, expiry=1.0, rate=0.05, vol=0.2)
>>> g.delta, g.vega, g.theta
(-0.378..., 38.13..., -1.69...)
>>> opal.implied_vol("call", price=3.85, spot=100, strike=105, expiry=0.5, rate=0.04)
0.1811...

Monte Carlo with a custom payoff:
>>> res = opal.mc_custom(lambda path: max(path[-1] - 100, 0.0),
...                      spot=100, expiry=1.0, rate=0.05, vol=0.2)
>>> res.price, res.std_error
(10.4..., 0.06...)

Conventions
-----------
- Rates and dividend yields are continuously compounded decimals (0.05 = 5%).
- Expiries are year fractions.
- Greeks from the core are per unit: divide ``vega`` by 100 for a per-vol-point
  number, ``theta`` by 365 for per calendar day, ``rho`` by 100 for per 1%.
"""

from opal._opal import (  # noqa: F401
    # result types
    Greeks,
    HestonGreeks,
    McResult,
    HestonParams,
    SabrParams,
    HullWhiteParams,
    DiscountCurve,
    SwaptionResult,
    CapletDetail,
    CapFloorResult,
    # vanilla analytics
    bs_price,
    bs_greeks,
    black76_price,
    bachelier_price,
    implied_vol,
    implied_vol_bachelier,
    # exotics
    digital_price,
    barrier_price,
    asian_price,
    lookback_price,
    # two-asset / compound analytics
    exchange_option_price,
    rainbow_option_price,
    two_asset_correlation_price,
    compound_option_price,
    # discrete cash dividends
    bs_discrete_div_price,
    binomial_discrete_div_price,
    # numerical engines
    binomial_price,
    trinomial_price,
    pde_price,
    mc_price,
    mc_custom,
    lsmc_price,
    lsmc_heston_price,
    # stochastic vol
    heston_price,
    heston_greeks,
    heston_mc,
    sabr_vol,
    sabr_normal_vol,
    sabr_price,
    # rates
    cap_floor_price,
    cap_floor_price_ois,
    swaption_price,
    swaption_price_ois,
    sabr_swaption_price,
    hw_zcb_option,
    hw_caplet,
    hw_floorlet,
    hw_cap,
    __version__,
)

__all__ = [
    "Greeks", "HestonGreeks", "McResult", "HestonParams", "SabrParams",
    "HullWhiteParams",
    "DiscountCurve", "SwaptionResult", "CapletDetail", "CapFloorResult",
    "bs_price", "bs_greeks", "black76_price", "bachelier_price",
    "implied_vol", "implied_vol_bachelier",
    "digital_price", "barrier_price", "asian_price", "lookback_price",
    "exchange_option_price", "rainbow_option_price",
    "two_asset_correlation_price", "compound_option_price",
    "bs_discrete_div_price", "binomial_discrete_div_price",
    "binomial_price", "trinomial_price", "pde_price", "mc_price", "mc_custom",
    "lsmc_price", "lsmc_heston_price",
    "heston_price", "heston_greeks", "heston_mc", "sabr_vol",
    "sabr_normal_vol", "sabr_price",
    "cap_floor_price", "cap_floor_price_ois", "swaption_price",
    "swaption_price_ois", "sabr_swaption_price",
    "hw_zcb_option", "hw_caplet", "hw_floorlet", "hw_cap",
    "__version__",
]
