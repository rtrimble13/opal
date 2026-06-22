# 3. Exotic payoffs

**Use case:** you trade beyond vanillas — digitals, barriers, Asians,
lookbacks, and multi-asset / compound structures — and want the closed-form
price plus an understanding of how each behaves.

Every instrument below is selected with `-i`/`--instrument`. Most have an
analytic price; all path-dependent ones can also be Monte-Carlo'd with
`--method mc` (which then reports a standard error — see vignette 2).

## Digitals and gap options

A cash-or-nothing digital pays a fixed `--cash` if it finishes in the money:

```sh
opal price -i digital-cash -t call -S 100 -K 100 -T 1 -r 5% -v 25% --cash 10
```

```text
  instrument    digital-cash
  price         5.040495
```

The price is `cash × e^{-rT} × N(d2)` ≈ `10 × 0.504` — i.e. the digital pays 10
with risk-neutral probability ~0.50 (just above ½ here because of the drift),
discounted. `digital-asset` pays one unit of the asset instead of cash; `gap`
uses one strike to decide exercise and another (`--payoff-strike`) to size the
payoff. Python: `opal.digital_price(...)`.

## Barriers (all eight, with rebates)

Barrier options knock in or out when spot touches `-H`. The four out/in ×
up/down families are `barrier-{down,up}-{out,in}`; `--rebate` pays a fixed
amount if a knock-out triggers:

```sh
opal price -i barrier-down-out -S 100 -K 100 -H 90 -T 1 -r 5% -v 25% --rebate 2
```

```text
  instrument    barrier-down-out
  method        analytic
  price         10.397176
```

This down-and-out call is cheaper than the ~12.3 vanilla: it dies if spot ever
falls to 90, and you are compensated only by the small `--rebate`. A useful
identity (zero rebate): **knock-in + knock-out = vanilla** — buying both is the
same as the unprotected option.

The closed form assumes **continuous** monitoring. To value the
discrete-monitoring premium (a real contract is often monitored daily), price it
with `--method mc`, which monitors at each step and is worth slightly more for a
knock-out. Python: `opal.barrier_price(...)`.

## Asian options

Asian payoffs use the **average** price, which damps volatility, so an Asian is
cheaper than the comparable vanilla. Geometric Asians have an exact closed form;
arithmetic Asians use the Turnbull–Wakeman approximation analytically, or Monte
Carlo with a geometric **control variate** for a very tight error:

```sh
opal price -i asian-arith -S 100 -K 100 -T 1 -r 5% -v 30% --paths 200000
```

```text
  instrument    asian-arith
  method        mc
  price         7.970725
  mc_std_error  0.000795
```

The control variate is why the standard error is ~0.0008 on a ~8.0 price with
only 200k paths — three orders of magnitude tighter than a naive estimator.
`asian-geo` defaults to the closed form. Python: `opal.asian_price(...)`.

## Lookbacks

Lookbacks pay off against the path extremum, so they are worth **more** than a
vanilla. `lookback-float` strikes at the realized min/max; `lookback-fixed`
uses `-K`. Seed a seasoned position's running extremum with `--extremum`:

```sh
opal price -i lookback-float -t call -S 100 -T 0.5 -r 5% -v 30%
```

Python: `opal.lookback_price(...)`.

## Partial-time (window) barriers

A partial-time-start barrier is monitored **only** over `[0, --window]`, then
becomes a plain vanilla to `-T`. It is worth more than a full-life barrier
(fewer chances to knock out):

```sh
opal price -i partial-down-out -t call -S 100 -K 100 -H 90 --window 0.5 -T 1 -r 5% -v 25%
```

```text
  instrument    partial-down-out
  barrier_type  down-out
  window_end    0.500000
  expiry_years  1.000000
  price         9.327327
```

`partial-{down,up}-{out,in}` are available; in + out = vanilla over the same
window. Python: `opal.partial_time_barrier_price(...)`.

## Two-asset options

These price a payoff on **two** correlated assets — pass the second via
`--spot2 --vol2 --div2` and the correlation `--rho`:

- `exchange` — Margrabe: option to swap asset 2 for asset 1, `max(S1-S2,0)`
  (rate-independent).
- `rainbow-max` / `rainbow-min` — Stulz: an option on the better/worse of two
  assets, struck at `-K`.
- `correlation` — pays on asset 2 vs `-K` only if asset 1 crosses `--trigger`.

These support the full greeks report (vignette 5 shows the convention):

```sh
opal greeks -i rainbow-max -t call -S 100 --spot2 95 -K 100 -T 1 -r 5% -v 25% --vol2 30% --rho 0.4
```

```text
  instrument    rainbow-max
  price         18.653208
  delta1        0.457093
  delta2        0.387874
  gamma1        0.016030
  gamma2        0.014143
  vega1         0.313689
  vega2         0.310374
  corr_sens     -5.441684
  theta         -0.032268
  rho_rate      0.639040
```

`delta1`/`delta2` are the hedge ratios in each underlying; `corr_sens` is
∂price/∂ρ — **negative** here because more correlation shrinks the dispersion of
the maximum, lowering the option. Python: `opal.rainbow_option_price(...)`,
`opal.exchange_option_price(...)`, `opal.two_asset_correlation_price(...)`.

## Compound options (option on an option)

A compound option is an outer option (expiry `-T`, strike `--strike1`) on an
inner vanilla (expiry `--expiry2`, strike `--strike2`):

```sh
opal price -i compound --outer call --inner call -S 100 --strike1 6 --strike2 100 -T 0.5 --expiry2 1 -r 5% -v 25%
```

```text
  instrument    compound
  outer         call
  inner         call
  outer_expiry  0.500000
  inner_expiry  1.000000
  price         7.758556
```

All four `--outer/--inner` call/put combinations are supported. Compound
put-call parity holds on the outer option. Python:
`opal.compound_option_price(...)`.

Next: [stochastic volatility](04-stochastic-volatility.md).
