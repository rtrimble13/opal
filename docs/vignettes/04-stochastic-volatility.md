# 4. Stochastic volatility

**Use case:** flat Black–Scholes vol can't fit a smile. Opal offers Heston
(stochastic variance) for equity-style options and SABR for the rates smile.

## Heston: pricing

Heston gives variance its own mean-reverting, correlated process. Pass the five
parameters:

- `--v0` initial variance, `--theta` long-run variance, `--kappa` mean-reversion
  speed, `--xi` vol-of-vol, `--rho` spot/vol correlation.

```sh
opal price --model heston -S 100 -K 100 -T 1 -r 3% --v0 0.04 --kappa 1.5 --theta 0.04 --xi 0.4 --rho -0.6
```

European Heston uses the semi-analytic characteristic-function formula (fast and
exact); American and exotic Heston fall back to Monte Carlo / Longstaff–Schwartz
(`--method mc`). A negative `--rho` produces the equity skew (downside puts
richer). Python: `opal.heston_price(...)`, with a `HestonParams(...)` object.

## Heston: first-class greeks

Generic bump-and-revalue greeks don't make sense under Heston (there is no
single lognormal vol to bump). `opal greeks --model heston` instead reports
**first-class Heston risk** — spot greeks plus sensitivities to the model's own
parameters:

```sh
opal greeks --model heston -S 100 -K 105 -T 1 -r 4% --v0 0.05 --kappa 1.5 --theta 0.06 --xi 0.4 --rho -0.7
```

```text
Opal | risk report
------------------
  model         heston
  method        analytic
  price         8.233829
  delta         0.600756
  gamma         0.018718
  vega          0.387873
  theta         -0.017752
  rho           0.518418
  dV_dv0        43.234356
  dV_dtheta     43.496602
  dV_dxi        -2.594558
  dV_drho       0.546695
  (vega per vol pt, theta per day, rho per 1% rate)
  (dV_dv0/dV_dtheta per unit variance; dV_dxi/dV_drho per unit)
```

- `delta`/`gamma` — spot sensitivities, as usual.
- `vega` — defined here as a **parallel shift of the variance level** (v0 and
  theta moved together, anchored at √v0), the closest analogue to BSM vega:
  +1 vol point adds **0.388**.
- `theta`/`rho` — time decay per day and rate sensitivity per 1%.
- `dV_dv0` / `dV_dtheta` — sensitivity to the initial and long-run **variance**
  (per unit variance, so large numbers — a 0.01 variance move is ~0.43).
- `dV_dxi` — vol-of-vol sensitivity. **Negative** here: with `--rho -0.7` more
  vol-of-vol fattens the left tail this OTM call doesn't benefit from.
- `dV_drho` — correlation sensitivity.

For American/exotic Heston the same report is produced by finite-differencing
the Monte Carlo price under common random numbers; the footer then notes that
gamma carries sampling noise. Python: `opal.heston_greeks(...)` returns a
`HestonGreeks` object with these fields.

## SABR: the smile

SABR is the market-standard parameterization of the implied-vol smile for
swaptions and rate options. `opal.sabr_vol(F, K, T, params)` returns the
lognormal implied vol; `opal.sabr_price(...)` prices through it; and SABR can
drive a swaption directly (vignette 8). Parameters: `alpha` (level), `beta`
(backbone), `rho` (skew), `nu` (smile/vol-of-vol).

```python
import opal
p = opal.SabrParams(alpha=0.25, beta=0.5, rho=-0.3, nu=0.4)
opal.sabr_vol(0.03, 0.025, 1.0, p)   # implied vol at a 2.5% strike
```

A negative `rho` tilts the skew (lower strikes carry higher vol); `nu` controls
the curvature. Next: [greeks and implied vol](05-greeks-and-implied-vol.md).
