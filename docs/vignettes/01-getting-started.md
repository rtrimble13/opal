# 1. Getting started

**Use case:** you have the repo and want to price your first option and read a
risk report.

## Build the CLI

Opal's core is header-only; the CLI is a single binary built with CMake:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/opal help      # usage
./build/opal version   # -> opal 0.2.0
```

`opal help` prints the full command, instrument, model and method reference.
`opal version` prints the single-source version (see
[docs/releasing.md](../releasing.md)).

## Your first price

```sh
opal price -S 100 -K 105 -T 0.5 -r 4% -v 22% -t call
```

```text
Opal | price
------------
  instrument    european
  type          call
  model         bsm
  method        analytic
  spot          100.0000
  strike        105.0000
  expiry_years  0.500000
  price         4.937141
```

Field by field:

- `instrument` / `type` / `model` / `method` — what was priced and how. With no
  `-i`, `--model` or `--method`, Opal defaults to a European option under
  Black–Scholes–Merton priced by its closed form. `method` shows the method
  **after** `auto` resolution, so you can confirm what actually ran.
- `spot` / `strike` / `expiry_years` — the market inputs, echoed back.
  `expiry_years` is the year fraction Opal used (here `-T 0.5`).
- `price` — the option premium, per unit of the underlying, in the same
  currency as spot/strike. A 6-month 105-strike call on a 100 spot at 22% vol
  is worth **4.94**.

The equivalent Python call:

```python
import opal
opal.bs_price("call", spot=100, strike=105, expiry=0.5, rate=0.04, vol=0.22)
# 4.937141...
```

## Your first risk report

`greeks` prices **and** reports sensitivities:

```sh
opal greeks -S 100 -K 100 -T 0.5 -r 5% -v 20%
```

```text
Opal | risk report
------------------
  instrument    european
  type          call
  model         bsm
  method        analytic
  spot          100.0000
  strike        100.0000
  expiry_years  0.500000
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

- `delta` — change in price per +1.00 change in spot. **0.598** ≈ "the option
  behaves like 0.6 shares." Hedge by shorting 0.598 units of the underlying.
- `gamma` — change in delta per +1.00 spot. **0.0274** means a +1 spot move
  lifts delta by ~0.027; how often you must re-hedge.
- `vega` — price change per **+1 vol point** (here the footer tells you it is
  already scaled per point): +1% vol adds **0.274** to the premium.
- `theta` — price change per **calendar day**: the option loses **0.022/day**
  to time decay.
- `rho` — price change per **+1% rate move**: **+0.264**.
- `vanna`, `volga`, `charm` — second-order greeks (∂delta/∂vol, ∂vega/∂vol,
  ∂delta/∂time). They are populated for analytic BSM and help you anticipate how
  delta and vega drift as vol or time move. They are blank for instruments
  priced numerically.

The footer restates the scaling so the numbers are unambiguous.

## Inputs: percentages and dates

Rates, dividend yields and vols accept decimals **or** percentages:

```sh
opal price -S 100 -K 100 -T 1 -r 0.04 -v 0.2     # decimals
opal price -S 100 -K 100 -T 1 -r 4%   -v 20%     # identical
```

Expiry accepts a year fraction or a calendar date (ACT/365F from today):

```sh
opal price -S 100 -K 100 -T 0.5          -r 4% -v 20%
opal price -S 100 -K 100 -T 2026-12-18   -r 4% -v 20%
```

## Output formats: table, JSON, CSV

Every command takes `-o`/`--output`. `table` (default) is human-readable; `json`
pipes into other tools; `csv` is offered for the tabular commands
(`chain`, `scenario`, `portfolio`).

```sh
opal price -S 100 -K 105 -T 0.5 -r 4% -v 22% -t call -o json
```

```text
{"instrument": "european", "type": "call", "model": "bsm", "method": "analytic", "spot": 100, "strike": 105, "expiry_years": 0.5, "price": 4.93714}
```

JSON is the easiest way to feed Opal into a spreadsheet, notebook or another
service. Next: [pricing and numerical methods](02-pricing-and-methods.md).
