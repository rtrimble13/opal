# 8. Interest-rate options

**Use case:** you price caps/floors and swaptions, and need multi-curve OIS
discounting and a choice of lognormal or normal vols.

The rates instruments are `cap`, `floor`, `caplet`, `floorlet`, `swaption` and
`zcb-option`. They take `-r` as the (flat) curve rate, `--notional`/`-n`
(default 100), and produce an undiscounted-to-PV price.

## A cap with OIS discounting

A cap is a strip of caplets on a forward rate. `--freq` sets payments per year;
`--ois-rate` switches on **multi-curve** pricing — forwards are projected off
`-r` while cashflows are discounted on a separate OIS curve:

```sh
opal price -i cap -K 4% -T 5 -r 4.2% --ois-rate 3.8% -v 30% --freq 4
```

```text
Opal | cap
----------
  instrument      cap
  model           black76
  ois_rate        0.038000
  caplets         19
  strike          0.0400
  maturity_years  5.0000
  price           3.589361
  notional        100.00
```

- `model black76` — the default for caps (lognormal forward rates); see vol
  types below.
- `ois_rate` — present because we passed `--ois-rate`; its appearance confirms
  dual-curve pricing is active. Drop the flag for single-curve.
- `caplets` — the strip length (19 quarterly caplets over 5y, less the first
  fixing).
- `price` — the cap premium on the `notional` of 100.

Lowering the discount curve (OIS 3.8% < projection 4.2%) **raises** the PV
versus single-curve, because the same forward cashflows are discounted less
heavily — the practical reason multi-curve matters post-2008.

## A swaption

A swaption is an option to enter a swap. `-t payer`/`receiver`, `-T` is the
option expiry, `--tenor` the underlying swap length:

```sh
opal price -i swaption -t payer -K 4% -T 1 --tenor 5 -r 4% -v 25%
```

```text
Opal | swaption
---------------
  instrument         swaption
  model              black76
  style              payer
  strike             0.040000
  forward_swap_rate  0.040403
  annuity            4.310644
  price              1.812031
```

- `forward_swap_rate` — the par forward rate of the underlying 5y swap starting
  in 1y (4.04%); the option is essentially at the money since `-K 4%`.
- `annuity` — the PV of 1bp of the swap's fixed leg (the "DV01 carrier"); the
  swaption price scales with it.
- `price` — the premium on the notional.

`swaption --model sabr` prices through a SABR smile (pass `--alpha --beta --nu
--rho`); `--ois-rate` adds OIS discounting just like the cap.

## Vol types: lognormal vs normal

`--vol-type lognormal` (default, Black-76) treats `-v` as a relative vol;
`--vol-type normal` (Bachelier) treats it as an **absolute** rate vol — the
right choice near zero or negative rates, where lognormal breaks down. Match the
convention your vol quote was calibrated in.

## Hull–White and ZCB options

`zcb-option` prices an option on a zero-coupon bond under a one-factor
Hull–White short-rate model (`--mean-rev`, `--hw-sigma`), with `--tenor` the
bond maturity. Hull–White also drives caplets/floorlets analytically.

## In Python

```python
import opal
curve = opal.DiscountCurve(0.042)                      # flat 4.2%
ois   = opal.DiscountCurve(0.038)
cap = opal.cap_floor_price_ois(ois, curve, strike=0.04, vol=0.30,
                               first_fixing=0.25, maturity=5.0, tau=0.25)
cap.price, len(cap.caplets)                            # total + per-caplet detail
sw = opal.swaption_price(curve, "payer", strike=0.04, vol=0.25, expiry=1.0, tenor=5.0)
sw.price, sw.forward_swap_rate, sw.annuity
```

`DiscountCurve` also accepts a zero-rate term structure —
`opal.DiscountCurve([1,2,5,10], [0.042,0.040,0.039,0.041])` — for a non-flat
curve. Next: [the Python API](09-python-api.md).
