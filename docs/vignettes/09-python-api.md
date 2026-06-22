# 9. The Python API

**Use case:** you want Opal in a notebook or script — for ad-hoc pricing,
calibration loops, custom payoffs, or an end-to-end pre-trade analysis.

```sh
pip install ./python        # builds the pybind11 extension once
```

```python
import opal
opal.__version__            # '0.2.0'
```

Every CLI capability has a Python entry point. Option type, exercise style and
barrier type are passed as strings (`"call"`, `"american"`, `"down-out"`).

## Pricers and greeks

```python
opal.bs_price("call", spot=42, strike=40, expiry=0.5, rate=0.10, vol=0.20)
# 4.759422...  (Hull's textbook value)

g = opal.bs_greeks("call", spot=100, strike=100, expiry=1.0, rate=0.05,
                   div=0.02, vol=0.25)
g.price, g.delta, g.gamma, g.vega, g.theta, g.rho
g.vanna, g.volga, g.charm           # higher-order, analytic BSM
```

Greeks are **per unit**: divide `g.vega` by 100 for per-vol-point, `g.theta` by
365 for per-day, `g.rho` by 100 for per-1%. `black76_price`, `bachelier_price`
and their `*_greeks`/implied-vol cousins follow the same shape.

## Implied vol

```python
c = opal.bs_price("call", spot=100, strike=110, expiry=0.7, rate=0.03,
                  div=0.01, vol=0.31)
opal.implied_vol("call", price=c, spot=100, strike=110, expiry=0.7,
                 rate=0.03, div=0.01)         # -> 0.31, round-trips
```

## Numerical engines

```python
# Trees and PDE (exercise = "european" | "american")
opal.binomial_price("put", "american", spot=100, strike=100, expiry=1.0,
                    rate=0.06, vol=0.25, steps=501, flavor="lr")
opal.trinomial_price("call", "european", spot=100, strike=105, expiry=0.75,
                     rate=0.04, div=0.01, vol=0.3, steps=600)
opal.pde_price("put", "american", spot=100, strike=100, expiry=1.0,
               rate=0.05, vol=0.25, grid=400)

# Monte Carlo for the built-in payoffs; returns McResult(price, std_error)
res = opal.mc_price("call", payoff="asian-arith", spot=100, strike=100,
                    expiry=1.0, rate=0.05, vol=0.3, paths=200000)
res.price, res.std_error
```

`payoff=` selects `vanilla`/`digital`/`barrier`/`asian-arith`/`asian-geo`/
`lookback-*`; barrier payoffs take `barrier=`/`barrier_type=`/`rebate=`.

## Custom path payoffs with `mc_custom`

Price *any* path-dependent structure by passing a Python function of the
simulated path. `path` is the list of simulated spot values **excluding** the
initial spot, with `path[-1]` the terminal price; the return value is settled
(undiscounted) at expiry and Opal discounts it. As a sanity check, a terminal
call payoff reproduces Black–Scholes:

```python
res = opal.mc_custom(lambda path: max(path[-1] - 100.0, 0.0),
                     spot=100, expiry=1.0, rate=0.05, vol=0.2,
                     paths=50000, steps=4)
res.price, res.std_error          # ~ opal.bs_price("call", 100, 100, 1, 0.05, 0.2)

# An up-and-in-at-110 terminal call, monitored at the `steps` path points:
res = opal.mc_custom(
    lambda path: max(path[-1] - 100.0, 0.0) if max(path) >= 110 else 0.0,
    spot=100, expiry=1.0, rate=0.05, vol=0.25, paths=100000, steps=52)
```

This is the escape hatch for payoffs Opal doesn't ship analytically. It is
slower (Python is called per path) but completely general; always read
`std_error` to size `paths`.

## American Monte Carlo (Longstaff–Schwartz)

```python
opal.lsmc_price("put", spot=100, strike=100, expiry=1.0, rate=0.06,
                vol=0.25, paths=50000, steps=50)        # McResult
```

## Stochastic vol and the smile

```python
hp = opal.HestonParams(v0=0.04, kappa=1.5, theta=0.05, xi=0.6, rho=-0.7)
opal.heston_price("call", spot=100, strike=100, expiry=1.0, rate=0.03,
                  div=0.0, params=hp)
hg = opal.heston_greeks("call", spot=100, strike=100, expiry=1.0, rate=0.03,
                        div=0.0, params=hp)
hg.delta, hg.vega, hg.dv0, hg.dtheta, hg.dxi, hg.drho    # first-class Heston risk
opal.lsmc_heston_price("put", spot=100, strike=100, expiry=1.0, rate=0.05,
                       div=0.0, params=hp)               # American under Heston

sp = opal.SabrParams(alpha=0.25, beta=0.5, rho=-0.3, nu=0.4)
opal.sabr_vol(0.03, 0.025, 1.0, sp)                      # implied vol at a strike
```

## Two-asset, compound and partial-time analytics

```python
opal.exchange_option_price(spot1=100, spot2=95, expiry=1.0, vol1=0.25,
                           vol2=0.30, rho=0.4)            # Margrabe (no rate)
opal.rainbow_option_price("call", "max", spot1=100, spot2=95, strike=100,
                          expiry=1.0, rate=0.05, vol1=0.25, vol2=0.30, rho=0.4)
opal.compound_option_price("call", "call", spot=100, strike1=6, strike2=100,
                           expiry1=0.5, expiry2=1.0, rate=0.05, vol=0.25)
opal.partial_time_barrier_price("call", "down-out", spot=100, strike=100,
                                barrier=90, window=0.5, expiry=1.0, rate=0.05,
                                vol=0.25)
```

## Rates and curves

```python
curve = opal.DiscountCurve(0.04)                         # flat 4%
curve.discount(2.0)                                      # exp(-0.08)

# A non-flat zero curve (linear in zero rates between pillars):
zc = opal.DiscountCurve([1, 2, 5, 10], [0.042, 0.040, 0.039, 0.041])

cap = opal.cap_floor_price(curve, strike=0.04, vol=0.25, first_fixing=0.25,
                           maturity=3.0, tau=0.25, is_cap=True)
cap.price, len(cap.caplets)
sw = opal.swaption_price(curve, "payer", strike=0.04, vol=0.3, expiry=1.0,
                         tenor=5.0)
sw.price, sw.forward_swap_rate, sw.annuity

hw = opal.HullWhiteParams(a=0.1, sigma=0.01)
opal.hw_zcb_option("call", curve, hw, 1.0, 3.0, 0.9)     # option on a ZCB
```

Multi-curve OIS variants are `cap_floor_price_ois(ois, proj, ...)` and
`swaption_price_ois(ois, proj, ...)`.

## End-to-end pre-trade workflow

Putting it together — value a candidate trade, see its smile-consistent price,
and stress it:

```python
import opal

S, K, T, r, q = 100, 105, 0.5, 0.04, 0.01

# 1. Flat-vol price and risk.
g = opal.bs_greeks("call", spot=S, strike=K, expiry=T, rate=r, div=q, vol=0.22)
print(f"price {g.price:.4f}  delta {g.delta:.4f}  vega/pt {g.vega/100:.4f}")

# 2. Same trade under a Heston smile (downside skew via rho < 0).
hp = opal.HestonParams(v0=0.22**2, kappa=1.5, theta=0.22**2, xi=0.5, rho=-0.6)
print("heston price", opal.heston_price("call", spot=S, strike=K, expiry=T,
                                        rate=r, div=q, params=hp))

# 3. A manual spot x vol stress grid (cf. the `scenario` CLI command).
base = opal.bs_price("call", spot=S, strike=K, expiry=T, rate=r, div=q, vol=0.22)
for ds in (-0.10, 0.0, 0.10):
    row = []
    for dv in (-0.05, 0.0, 0.05):
        p = opal.bs_price("call", spot=S*(1+ds), strike=K, expiry=T, rate=r,
                          div=q, vol=0.22+dv)
        row.append(f"{p-base:+7.4f}")
    print(f"spot {S*(1+ds):6.1f}: " + "  ".join(row))
```

See `examples/pretrade_analysis.py` and `examples/opal_walkthrough.ipynb` for
fuller, runnable versions (portfolio construction, smile cubes, scenario
grids). For the CLI equivalents of each step here, see vignettes
[5](05-greeks-and-implied-vol.md) and [6](06-scenario-analysis.md).
