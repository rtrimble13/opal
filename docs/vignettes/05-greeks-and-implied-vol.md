# 5. Greeks and implied vol

**Use case:** you have a position (or a market price) and need the full risk
picture — first- and second-order greeks, the implied vol behind a quote, and a
strike ladder.

## The full risk report

`opal greeks` adds the sensitivities to `price`. For analytic BSM you get the
higher-order greeks too:

```sh
opal greeks -S 100 -K 100 -T 0.5 -r 5% -v 20%
```

```text
  price         6.888729
  delta         0.597734
  gamma         0.027359
  vega          0.273587
  theta         -0.022236
  rho           0.264424
  vanna         -0.205190
  volga         3.590824
  charm         -0.095755
  (vega per vol pt, theta per day, rho per 1% rate)
```

First-order (covered in vignette 1): `delta`, `gamma`, `vega`, `theta`, `rho`.
The second-order greeks tell you how the first-order ones **move**:

- `vanna` = ∂delta/∂vol = ∂vega/∂spot. Negative here: a vol spike lowers this
  call's delta. Matters for risk-reversals and skew hedging.
- `volga` (vomma) = ∂vega/∂vol. Positive: vega grows as vol rises, so a long
  option is long "vol convexity" — the engine of volatility trading.
- `charm` = ∂delta/∂time. Tells you how your delta hedge drifts overnight even
  if spot doesn't move.

These appear only when priced analytically (BSM europeans). For numerically
priced instruments the report shows the first-order greeks; for Heston it shows
the model-parameter sensitivities instead (vignette 4); for two-asset/compound
it shows per-asset greeks (vignette 3).

The CLI prints **trader conventions** — vega per vol point, theta per calendar
day, rho per 1% rate — as the footer reminds you. The library values are per
unit (per 1.00 vol, per year, per 1.00 rate); in Python you scale yourself:
`g.vega / 100`, `g.theta / 365`, `g.rho / 100`.

```python
g = opal.bs_greeks("call", spot=100, strike=100, expiry=0.5, rate=0.05, vol=0.20)
g.delta, g.vega / 100, g.theta / 365, g.rho / 100
```

## Implied volatility from a market price

Given a traded price, `implied` inverts BSM for the vol — and throws back a
couple of greeks at that vol so you can immediately gauge sensitivity:

```sh
opal implied -S 100 -K 105 -T 0.5 -r 4% --price 3.85
```

```text
Opal | implied volatility
-------------------------
  model         bsm
  type          call
  market_price  3.850000
  implied_vol   0.181116
  vega          0.278473
  delta         0.436138
```

- `implied_vol` — the BSM vol that reproduces 3.85: **18.1%**. The solver is a
  safeguarded Newton/Brent hybrid, robust across the surface.
- `vega` — at that vol, the price moves **0.278 per vol point**, so a 0.1-point
  vol quote error is ~0.03 of premium.
- `delta` — the hedge ratio at the implied vol.

`--model black76`/`bachelier` invert those models instead. Python:
`opal.implied_vol(...)`, `opal.implied_vol_bachelier(...)`.

## Option chains

`chain` prices a ladder of strikes in one shot — the quickest way to see the
skew and the greeks across moneyness:

```sh
opal chain -S 100 -T 0.5 -r 4% -v 25% --strikes 90:110:5
```

```text
strike     call      put  call_delta  put_delta    gamma    vega
----------------------------------------------------------------
 90.00  14.1127   2.3306      0.7874    -0.2126  0.01642  0.2052
 95.00  10.7854   3.9043      0.6885    -0.3115  0.02000  0.2500
100.00   8.0080   6.0279      0.5799    -0.4201  0.02211  0.2764
105.00   5.7799   8.7008      0.4703    -0.5297  0.02251  0.2813
110.00   4.0599  11.8818      0.3678    -0.6322  0.02132  0.2665
```

`--strikes lo:hi:step` sets the ladder. Each row gives the call and put prices,
their deltas (note `call_delta − put_delta ≈ 1`, put-call parity on delta), and
the shared `gamma`/`vega`. Gamma and vega peak near the money (strike ≈ 100–105
here) — that's where the position is most sensitive to spot and vol. Add `-o
csv` to import the ladder into a spreadsheet.

Next: [scenario analysis](06-scenario-analysis.md).
