# Opal

Institutional-grade option pricing and risk analytics: a header-only C++17
library, a command line tool for pre-trade analysis, position monitoring and
risk management, and a Python API for scripting and Jupyter workflows.

## Capabilities

**Models**

| Model | Use |
|---|---|
| Black–Scholes–Merton | Equity options with continuous dividend yield, full analytic greeks (incl. vanna, volga, charm) |
| Discrete cash dividends | Escrowed-spot model for Europeans; Leisen–Reimer tree with dividend add-back for American early exercise |
| Black-76 | Options on forwards/futures, caps/floors, swaptions |
| Bachelier (normal) | Rates options in low/negative rate regimes |
| Heston | Stochastic volatility, semi-analytic via characteristic function ("little trap" formulation) + Monte Carlo for exotics, Longstaff–Schwartz for American exercise |
| SABR (Hagan 2002) | Vol smile for swaptions and rate options, lognormal and normal vols |
| Hull–White 1F | Bond options, caplets/floorlets, caps/floors on a fitted curve |

Caps/floors and swaptions support **multi-curve OIS discounting**: forwards
projected off one curve, cashflows discounted on another.

**Payoffs**

European and American vanillas; cash- and asset-or-nothing digitals; gap
options; all eight single-barrier variants (up/down × in/out × call/put) with
rebates; arithmetic and geometric Asians (fixed and floating strike,
continuous and discrete monitoring); fixed- and floating-strike lookbacks
(including seasoned positions via running extrema); two-asset analytics
(Margrabe exchange, Stulz options on the max/min of two assets, two-asset
correlation options), compound options (Geske option-on-option) and
partial-time-start (window) barriers; caps, floors, caplets, floorlets;
European payer/receiver swaptions; zero-coupon bond options. Arbitrary path-dependent payoffs via the Monte Carlo engines
(including from Python).

**Engines**

| Engine | Notes |
|---|---|
| Closed form | Reiner–Rubinstein barriers, Goldman–Sosin–Gatto / Conze–Viswanathan lookbacks, Turnbull–Wakeman Asians, exact geometric Asians |
| Leisen–Reimer binomial | Smooth O(1/n²) convergence; default for American exercise |
| CRR binomial, trinomial | Cross-checking and pedagogy |
| Crank–Nicolson PDE | Rannacher start-up, American exercise, knock-out barriers as boundary conditions |
| Monte Carlo | Antithetic variates, geometric-Asian control variate, standard errors reported, Heston full-truncation Euler |
| Longstaff–Schwartz | American exercise by least-squares Monte Carlo under GBM and Heston (variance in the regression basis) |

**Analytics**: implied vol (Newton + Brent safeguarded) for BSM/Black-76/
Bachelier, analytic and bump-and-revalue greeks, scenario grids, option
chains, portfolio aggregation.

## Build

Requires CMake ≥ 3.16 and a C++17 compiler. The library itself is
header-only (`#include "opal/opal.hpp"`).

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/opal_tests          # unit tests
./build/opal help           # the CLI
```

## Command line

```text
opal <command> [options]
  price | greeks | implied | chain | scenario | portfolio
```

```sh
# Vanilla pricing and a full risk report
opal price  -S 100 -K 105 -T 0.5 -r 4% -v 22% -t call
opal greeks -S 100 -K 100 -T 0.5 -r 5% -v 20%

# American, barrier, Asian, Heston
opal price -i american -t put -S 50 -K 55 -T 1 -r 3% -q 1% -v 30%
opal price -i barrier-up-out -S 100 -K 100 -H 120 -T 1 -r 5% -v 25%
opal price -i asian-arith -S 100 -K 100 -T 1 -r 5% -v 30% --paths 200000
opal price --model heston -S 100 -K 100 -T 1 -r 3% --v0 0.04 --kappa 1.5 --theta 0.04 --xi 0.4 --rho -0.6

# Two-asset (Margrabe / Stulz / correlation) and compound (Geske) analytics
opal price -i rainbow-max -t call -S 100 --spot2 95 -K 100 -T 1 -r 5% -v 25% --vol2 30% --rho 0.4
opal price -i exchange -S 100 --spot2 95 -T 1 -r 5% -v 25% --vol2 30% --rho 0.4
opal price -i compound --outer call --inner call -S 100 --strike1 6 --strike2 100 -T 0.5 --expiry2 1 -r 5% -v 25%
opal price -i partial-down-out -t call -S 100 -K 100 -H 90 --window 0.5 -T 1 -r 5% -v 25%
# First-class Heston risk: delta/gamma, vega (parallel variance shift) and
# dV/dv0, dV/dtheta, dV/dxi, dV/drho sensitivities
opal greeks --model heston -S 100 -K 105 -T 1 -r 4% --v0 0.05 --kappa 1.5 --theta 0.06 --xi 0.4 --rho -0.7

# Discrete cash dividends (escrowed model; American captures early exercise)
opal price -i american -t call -S 100 -K 100 -T 1 -r 4% -v 25% --dividends 0.25:1.50,0.75:1.50

# American by Longstaff-Schwartz Monte Carlo (also under --model heston)
opal price -i american -t put -S 100 -K 100 -T 1 -r 5% -v 25% --method mc
opal price --model heston -i american -t put -S 100 -K 100 -T 1 -r 5% --v0 0.04 --xi 0.5 --rho -0.7

# Rates: caps, swaptions (Black-76, Bachelier, SABR, Hull-White); --ois-rate
# switches to multi-curve pricing (forwards off --rate, OIS discounting)
opal price -i cap -K 4% -T 5 -r 4.2% -v 30% --freq 4
opal price -i cap -K 4% -T 5 -r 4.2% --ois-rate 3.8% -v 30% --freq 4
opal price -i swaption -t payer -K 4% -T 1 --tenor 5 -r 4% -v 25%
opal price -i swaption --model sabr -t payer -K 4% -T 1 --tenor 5 -r 4% --alpha 0.04 --beta 0.5 --nu 0.45 --rho -0.25
opal price -i zcb-option -t call -K 0.9 -T 1 --tenor 3 -r 5% --mean-rev 0.1 --hw-sigma 0.01

# Market-implied vol, chains, scenarios, portfolio risk
opal implied -S 100 -K 105 -T 0.5 -r 4% --price 3.85
opal chain -S 100 -T 0.5 -r 4% -v 25% --strikes 80:120:5
opal chain -S 100 -T 1 -r 4% --sabr 0.2:1.0:-0.3:0.4 --strikes 80:120:5   # ladder off a SABR smile
opal scenario -S 100 -K 100 -T 0.5 -r 4% -v 25% --spot-range -20:20:5 --vol-range -10:10:5
opal portfolio --file examples/book.csv
```

Expiries accept year fractions (`-T 0.5`) or dates (`-T 2026-12-18`,
ACT/365F). Rates accept decimals or percentages (`-r 0.04` ≡ `-r 4%`).
Every command takes `-o json` (and `-o csv` for tabular output) for piping
into other systems. Numerical method selection is automatic but can be forced
with `--method analytic|lr|crr|trinomial|pde|mc`.

Example scenario output (P&L of an ATM call, spot × vol grid):

```text
P&L vs base price 8.0080 (rows: spot, cols: vol)
spot\vol    20.0%    25.0%    30.0%
-----------------------------------
   90.00  -5.7836  -4.6419  -3.4454
  100.00  -1.3809   0.0000   1.3824
  110.00   5.6875   6.8072   8.0283
```

## Python API

```sh
pip install ./python        # builds the pybind11 extension
```

```python
import opal

# Pricing and greeks
opal.bs_price("call", spot=100, strike=105, expiry=0.5, rate=0.04, vol=0.22)
g = opal.bs_greeks("put", spot=100, strike=100, expiry=1.0, rate=0.05, vol=0.2)
g.delta, g.vega / 100, g.theta / 365

# Implied vol, exotics, engines
opal.implied_vol("call", price=3.85, spot=100, strike=105, expiry=0.5, rate=0.04)
opal.barrier_price("call", "up-out", spot=100, strike=100, barrier=120,
                   expiry=1.0, rate=0.05, vol=0.25)
opal.binomial_price("put", "american", spot=50, strike=55, expiry=1.0,
                    rate=0.03, vol=0.3)

# Two-asset and compound analytics (bivariate-normal closed forms)
opal.rainbow_option_price("call", "max", spot1=100, spot2=95, strike=100,
                          expiry=1.0, rate=0.05, vol1=0.25, vol2=0.30, rho=0.4)
opal.compound_option_price("call", "call", spot=100, strike1=6, strike2=100,
                           expiry1=0.5, expiry2=1.0, rate=0.05, vol=0.25)
opal.partial_time_barrier_price("call", "down-out", spot=100, strike=100,
                                barrier=90, window=0.5, expiry=1.0, rate=0.05,
                                vol=0.25)

# Monte Carlo with standard errors — including custom Python payoffs
res = opal.mc_price("call", payoff="asian-arith", spot=100, strike=100,
                    expiry=1.0, rate=0.05, vol=0.3, paths=200000)
res.price, res.std_error

cliquet = opal.mc_custom(lambda path: max(path[-1] / path[0] - 1.0, 0.0) * 100,
                         spot=100, expiry=1.0, rate=0.04, vol=0.2)

# Discrete dividends and American exercise (incl. Longstaff-Schwartz MC)
opal.bs_discrete_div_price("call", spot=100, strike=100, expiry=1.0, rate=0.04,
                           vol=0.25, dividends=[(0.25, 1.5), (0.75, 1.5)])
opal.binomial_discrete_div_price("call", "american", spot=100, strike=100,
                                 expiry=1.0, rate=0.04, vol=0.25,
                                 dividends=[(0.5, 3.0)])
opal.lsmc_price("put", spot=100, strike=100, expiry=1.0, rate=0.06, vol=0.25)

# Stochastic vol and rates
hp = opal.HestonParams(v0=0.04, kappa=1.5, theta=0.05, xi=0.6, rho=-0.7)
opal.heston_price("call", spot=100, strike=100, expiry=1.0, rate=0.03, div=0.0, params=hp)
opal.lsmc_heston_price("put", spot=100, strike=100, expiry=1.0, rate=0.05, div=0.0, params=hp)
# First-class Heston greeks: spot delta/gamma, parallel-variance-shift vega,
# theta, rate rho, plus dV/dv0, dV/dtheta, dV/dxi and dV/drho.
hg = opal.heston_greeks("call", spot=100, strike=100, expiry=1.0, rate=0.03, div=0.0, params=hp)
hg.delta, hg.vega / 100, hg.dv0, hg.drho

# SABR smile / surface container: one calibrated set serves a vol per strike/expiry
sp = opal.SabrParams(alpha=0.2, beta=1.0, rho=-0.3, nu=0.4)
smile = opal.SabrSmile(forward=100, expiry=1.0, params=sp)
smile.vol(90), smile.vol(110)                       # skew: lower strike richer
surf = opal.VolSurface([opal.SabrSmile(100, 0.5, sp), opal.SabrSmile(100, 2.0, sp)])
surf.vol(strike=95, expiry=1.0)                     # interpolated in total variance

curve = opal.DiscountCurve([1, 2, 5, 10], [0.042, 0.040, 0.039, 0.041])
opal.swaption_price(curve, "payer", strike=0.04, vol=0.25, expiry=1.0, tenor=5.0)

# Multi-curve OIS discounting
ois = opal.DiscountCurve(0.035)
opal.swaption_price_ois(ois, curve, "payer", strike=0.04, vol=0.25, expiry=1.0, tenor=5.0)
opal.cap_floor_price_ois(ois, curve, strike=0.04, vol=0.3, first_fixing=0.25, maturity=5.0)
```

Release wheels for Linux/macOS/Windows and a self-contained sdist are built
by the `Wheels` workflow (cibuildwheel) on version tags.

See `examples/opal_walkthrough.ipynb` for a full Jupyter tour (smiles under
Heston, SABR cubes, scenario grids) and `examples/pretrade_analysis.py` for a
scripted pre-trade workflow.

## Vignettes

For narrative, worked walk-throughs — every CLI command and Python function,
with real output explained field by field, and end-to-end use cases — see the
[vignettes](docs/vignettes/README.md):

1. [Getting started](docs/vignettes/01-getting-started.md)
2. [Pricing and numerical methods](docs/vignettes/02-pricing-and-methods.md)
3. [Exotic payoffs](docs/vignettes/03-exotic-payoffs.md)
4. [Stochastic volatility](docs/vignettes/04-stochastic-volatility.md)
5. [Greeks and implied vol](docs/vignettes/05-greeks-and-implied-vol.md)
6. [Scenario analysis](docs/vignettes/06-scenario-analysis.md)
7. [Portfolio risk](docs/vignettes/07-portfolio-risk.md)
8. [Interest-rate options](docs/vignettes/08-interest-rate-options.md)
9. [The Python API](docs/vignettes/09-python-api.md)

## Portfolio files

`opal portfolio --file book.csv` consumes a CSV with columns
`instrument,type,quantity,spot,strike,vol,rate,div,expiry[,barrier,rebate,method]`
and reports per-position and aggregated NPV, delta, gamma, vega and theta —
see `examples/book.csv`.

Cell formats match the command-line flags: numeric cells accept a trailing
`%` (e.g. `rate=4%`, `vol=22%`), and `expiry` takes either a year fraction
(`0.5`) or a `YYYY-MM-DD` date (ACT/365F from today). A malformed cell is
reported with its row and column, e.g.
`portfolio row 3, column 'vol': expected a number, got 'oops'`.

## Conventions

- Rates and dividend yields are continuously compounded decimals.
- Library-level greeks are per unit of the underlying variable; the CLI
  reports trader conventions (vega per vol point, theta per calendar day,
  rho per 1% rate move).
- Heston greeks (`heston_greeks` / `opal greeks --model heston`) report a
  "vega" defined as a parallel shift of the variance level (v0 and theta moved
  together, anchored at sqrt(v0)) alongside the per-parameter sensitivities
  dV/dv0, dV/dtheta, dV/dxi and dV/drho (v0/theta per unit variance). European
  options use the semi-analytic price; American and exotic Heston instruments
  finite-difference the Monte Carlo / Longstaff–Schwartz price under common
  random numbers, so first-order greeks are stable while gamma carries the
  engine's sampling noise.
- Barrier closed forms assume continuous monitoring; the Monte Carlo engine
  monitors discretely at each step (worth more for knock-outs) — use it to
  quantify discrete-monitoring premia.

## Validation

`opal_tests` (204 checks) validates against Hull and Haug reference values,
no-arbitrage identities (put-call parity, digital parity, barrier in/out
parity), cross-engine agreement (analytic vs trees vs PDE vs Monte Carlo),
model degeneracies (Heston → BS, SABR → flat lognormal), short-dated /
high-vol-of-vol Heston pricing against Monte Carlo, and dense-monitoring
Monte Carlo for the exotic closed forms. The Python suite
(`python/tests/test_api.py`) re-checks the bindings end to end.

## Repository layout

```
include/opal/      header-only C++ library (models, engines, rates, analytics)
cli/               the opal command line tool
python/            pybind11 bindings + setuptools package
tests/             C++ unit tests (custom micro-framework, no dependencies)
examples/          sample portfolio, pre-trade script, Jupyter notebook
scripts/           release helper (scripts/release.sh)
docs/              releasing & versioning, project review
```

## Releasing

Opal uses tag-driven releases with the git tag as the single source of truth:
`scripts/release.sh X.Y.Z` (or pushing a `vX.Y.Z` tag) cuts a GitHub Release
with auto-generated notes and attached CLI binaries, wheels and sdist. The
version lives only in `include/opal/version.hpp`; CMake and the Python package
derive from it, and CI fails any release whose tag disagrees. See
[docs/releasing.md](docs/releasing.md).

## License

MIT — see [LICENSE](LICENSE).
