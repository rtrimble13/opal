# 2. Pricing and numerical methods

**Use case:** you need to price European and American options, and decide which
numerical method to trust for accuracy vs. speed.

## Models: bsm, black76, bachelier

`--model` selects the dynamics:

- `bsm` (default) — Black–Scholes–Merton on a spot with continuous dividend
  yield `-q`.
- `black76` — option on a forward/futures; pass the forward as `-S`. Equivalent
  to BSM with zero cost of carry.
- `bachelier` — normal (absolute-vol) model for low/negative-rate regimes; `-v`
  is an absolute vol, not a percentage of spot.

```sh
opal price --model black76 -S 105 -K 100 -T 0.75 -r 3% -v 20%   # forward = 105
opal price --model bachelier -S 100 -K 100 -T 1 -r 2% -v 15      # vol in price units
```

Python equivalents: `opal.black76_price(...)`, `opal.bachelier_price(...)`.

## European vs American, and method selection

Opal resolves `--method auto` (the default) to the best available method per
instrument, but you can force any method. For a vanilla European the analytic
BSM formula is exact and instant. For **American** exercise there is no closed
form, so Opal defaults to a Leisen–Reimer binomial tree:

```sh
opal price -i american -t put -S 50 -K 55 -T 1 -r 3% -q 1% -v 30%
```

```text
Opal | price
------------
  instrument    american
  type          put
  model         bsm
  method        lr
  spot          50.0000
  strike        55.0000
  expiry_years  1.000000
  price         8.446463
```

`method lr` confirms the Leisen–Reimer tree ran. The 8.45 price exceeds the
European put (which cannot be exercised early) — the difference is the
early-exercise premium.

### The methods, and when to pick which

| `--method` | What it is | Use it when |
|------------|------------|-------------|
| `analytic` | Closed form | A formula exists (European BSM/Black-76/Bachelier, digitals, barriers, lookbacks, geometric Asians) — always preferred. |
| `lr` | Leisen–Reimer binomial | American exercise; smooth, fast O(1/n²) convergence (the default for American). |
| `crr` | Cox–Ross–Rubinstein binomial | Cross-checking / teaching; needs more steps than LR. |
| `trinomial` | Trinomial tree | Cross-checking; stable for barriers. |
| `pde` | Crank–Nicolson finite differences | American exercise and knock-out barriers as boundary conditions; smooth greeks. |
| `mc` | Monte Carlo (LSMC for American) | Path-dependent payoffs, or when you want an error bar. Slowest; the only general option for arithmetic Asians under American exercise. |

Tune trees/PDE with `--steps`/`--grid`, and Monte Carlo with
`--paths`/`--steps`/`--seed`.

## The Monte Carlo standard error

When the resolved method is `mc`, Opal reports the **standard error** of the
estimate so you can judge whether the price is precise enough:

```sh
opal price -i american -t put -S 50 -K 55 -T 1 -r 3% -q 1% -v 30% --method mc --paths 50000
```

```text
  method        mc
  ...
  price         8.427722
  mc_std_error  0.032315
```

Read this as **8.43 ± 0.03 (1 s.e.)**. The true price lies within roughly
±2 s.e. (≈ ±0.06) about 95% of the time. Compare to the LR price (8.45 above):
the LSMC estimate agrees within ~0.6 s.e., and is a slightly low-biased
estimator as expected for least-squares Monte Carlo. Quadruple `--paths` to
halve the standard error.

## Discrete cash dividends

Real equities pay discrete cash dividends, not a smooth yield. Pass them as
`t:amount` pairs; Opal uses the escrowed-spot model (European, analytic) or a
dividend-aware Leisen–Reimer tree (American, captures early exercise around the
dividend):

```sh
opal price -i american -t call -S 100 -K 100 -T 1 -r 4% -v 25% \
           --dividends 0.25:1.50,0.75:1.50
```

`--dividends` requires the `bsm` model and a european/american instrument.

Python: `opal.bs_discrete_div_price(...)` and
`opal.binomial_discrete_div_price(...)`.

Next: [exotic payoffs](03-exotic-payoffs.md).
