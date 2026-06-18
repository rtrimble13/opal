# Opal — Project Review (2026-06-18)

Reviewer: senior-staff engineering health-check (whole-project, five lenses).
Scope: entire repository at branch `claude/nice-darwin-o22o1h`.

## Verdict & summary

Opal is a genuinely high-quality, header-only C++17 quantitative options
library with a CLI and pybind11 Python API. The numerical core is careful and
correct where it counts: generalized Black–Scholes (cost-of-carry form) reused
across every analytic model, the "little Heston trap" characteristic function,
Rannacher start-up for Crank–Nicolson, Leisen–Reimer trees, antithetic +
control-variate Monte Carlo, ridge-stabilized LSMC regression, and multi-curve
OIS discounting. The test suite is a real strength — cross-engine convergence,
put–call and in–out parity, textbook reference values (Hull, Haug), and Monte
Carlo bands at k·SE.

The standing issues are concentrated in the **glue layer** (the CLI/`pricer.hpp`
dispatch and the Python bindings), not in the math. The single most important
finding is that the CLI's risk tooling silently drops the volatility argument
under the Heston model, so `greeks` and `scenario` report a zero/meaningless
vega — wrong risk numbers from a risk tool, presented without warning. The rest
are mid-tier robustness, refactoring, and enhancement items.

Overall: solid foundation, ship-worthy core, with one correctness gap worth
fixing before anyone trusts Heston greeks in anger.

## How this review was scoped

**Read in full:** the entire numerical core and public API —
`math/{normal,solvers}.hpp`, `models/{black_scholes,heston,sabr,barrier,asian,
lookback,digital,discrete_div}.hpp`, `engines/{lattice,pde,lsmc,monte_carlo}.hpp`,
`analytics/{implied_vol,greeks}.hpp`, `rates/{curve,cap_floor,hull_white,
swaption}.hpp`, `core/{types,dividends}.hpp`; the CLI (`cli/{main.cpp,args.hpp,
pricer.hpp,format.hpp}`); the Python layer (`python/{bindings.cpp,opal/__init__.py,
setup.py,pyproject.toml,MANIFEST.in}`); the C++ and Python tests; CMake and both
GitHub workflows.

**Sampled:** the example scripts / notebook were listed but not deep-read; git
history was used only for churn.

**Confidence:** high on the analytic and engine code (it is also well tested);
high on the CLI dispatch findings (traced by hand through `price_at`); medium on
the Heston-integration-truncation finding (plausible but not reproduced against a
reference). Each finding below carries its own confidence.

The data flow for the headline finding was traced end to end:
`numerical_greeks` (greeks.hpp:34) bumps `vol` → `EquityTrade::price_at`
(pricer.hpp:131) → for Heston it calls `price_heston` (pricer.hpp:150) which
ignores the `v` argument entirely.

## Findings (ordered by impact × effort priority)

Each actionable finding links to its GitHub issue.

### P0

1. **Heston risk is silently wrong: `price_at` ignores its `vol` argument, so
   CLI `greeks`/`scenario` report vega ≈ 0 / a flat vol axis.**
   Lens: hidden bug · Severity: High · Impact: High · Effort: Medium ·
   Confidence: high.
   `cli/pricer.hpp:150` dispatches Heston via `price_heston(m, s, rr, tt, cfg)`,
   dropping `v`; `price_heston` (pricer.hpp:261) prices off the fixed `heston`
   struct. `trade_greeks` (cli/main.cpp:239) routes all non-`bsm` instruments to
   `numerical_greeks`, which computes vega/vanna/volga by bumping `vol`
   (greeks.hpp:34) — a no-op under Heston, so vega is reported as exactly 0.
   `cmd_scenario` (main.cpp:371) has the same defect: every vol column is
   identical. → [#3](https://github.com/rtrimble13/opal/issues/3)

### P1

2. **CRR and trinomial trees throw on legitimate inputs (high rate / coarse
   steps / low vol) instead of degrading.**
   Lens: robustness/bug · Severity: Medium · Impact: Medium · Effort: Low ·
   Confidence: high.
   `binomial_crr_price` requires `0 < p < 1` (lattice.hpp:28); when `e^{b·dt}`
   leaves `[d,u]` (large `b·dt`) it throws "arbitrage in tree". `trinomial_price`
   requires `pu,pd,pm > 0` (lattice.hpp:116). Both are user-selectable via
   `--method crr|trinomial` and reachable with ordinary inputs. → [#4](https://github.com/rtrimble13/opal/issues/4)

3. **Payoff dispatch is triplicated and already diverging.**
   Lens: refactoring · Impact: Medium · Effort: Medium · Confidence: high.
   The same vanilla/digital/asian/barrier/lookback selection is written three
   times: `EquityTrade::price_at` (pricer.hpp:152-243), `price_heston`
   (pricer.hpp:273-296), and the Python `mc_price` binding (bindings.cpp:295-317).
   They already disagree — e.g. the `mc_price` barrier payoff omits the
   `rebate`/`r`/`T` arguments (bindings.cpp:314) that the CLI passes. → [#5](https://github.com/rtrimble13/opal/issues/5)

### P2

4. **Monte Carlo standard error is surfaced only for `asian-arith`.**
   Lens: enhancement · Impact: Medium · Effort: Low · Confidence: high.
   `EquityTrade::mc_std_error` (pricer.hpp:249-258) returns NaN for every other
   MC-priced instrument, so `opal price --method mc` on a european/barrier/
   lookback/digital hides the very error bar that justifies Monte Carlo. The
   engines already return it. → [#6](https://github.com/rtrimble13/opal/issues/6)

5. **`portfolio` CSV parsing gives no row context and silently diverges from the
   CLI's date-aware expiry handling.**
   Lens: robustness/usability · Severity: Low · Impact: Medium · Effort: Low ·
   Confidence: high.
   `cmd_portfolio` (main.cpp:439-445) calls `std::stod` on each cell; a malformed
   value surfaces as bare `opal: stod` with no line/column. `expiry` is parsed as
   a plain number (main.cpp:445), unlike `Args::get_expiry` which also accepts
   `YYYY-MM-DD` everywhere else. → [#7](https://github.com/rtrimble13/opal/issues/7)

6. **Heston semi-analytic integration uses a fixed truncation `[1e-8, 200]` and a
   single absolute tolerance.**
   Lens: robustness · Severity: Low · Impact: Medium · Effort: Medium ·
   Confidence: medium.
   `detail::heston_prob` (heston.hpp:67) integrates to a hard-coded 200 with
   `tol = 1e-10`. For short maturities or high vol-of-vol the integrand may not
   have decayed by 200, and an absolute tol on a magnitude-varying integrand can
   under-resolve. Not reproduced against a reference; flagged for verification.
   → [#8](https://github.com/rtrimble13/opal/issues/8)

7. **Project version disagrees across sources of truth.**
   Lens: bug/consistency · Severity: Low · Impact: Low · Effort: Low ·
   Confidence: high.
   `CMakeLists.txt:2` declares `project(opal VERSION 0.1.0)` while
   `include/opal/version.hpp` and `python/setup.py` both say `0.2.0`. CMake-driven
   install/package metadata will report the stale value. → [#9](https://github.com/rtrimble13/opal/issues/9)

### P3

8. **Cashflow schedules are built by floating-point accumulation.**
   Lens: robustness · Severity: Low · Impact: Low · Effort: Low ·
   Confidence: medium.
   `swaption.hpp:35`, `cap_floor.hpp:45`, and `hull_white.hpp:73` all loop
   `for (double t = ...; t += tau)` with a `1e-10` end guard. Exact for
   power-of-two `tau` (quarterly/semiannual) but accumulates error for e.g.
   monthly `1/12` over long tenors, where the final period can be inconsistently
   included/dropped. Prefer an integer period count and `i·tau`. → [#10](https://github.com/rtrimble13/opal/issues/10)

9. **`mc_heston` lacks the input validation `mc_gbm` performs.**
   Lens: robustness · Severity: Low · Impact: Low · Effort: Low ·
   Confidence: high.
   `mc_gbm` guards `S>0, sigma>=0, T>0` (monte_carlo.hpp:37); `mc_heston`
   (monte_carlo.hpp:273) guards only paths/steps, so a non-positive `S`/`T`
   silently produces garbage rather than a clear error. → [#11](https://github.com/rtrimble13/opal/issues/11)

10. **Dead/unwired math capability: `bivar_norm_cdf` (and `norm_ppf`) are used
    only by tests.**
    Lens: refactoring · Impact: Low · Effort: Low · Confidence: high.
    `math::bivar_norm_cdf` (normal.hpp:72) is a careful 70-line Genz bivariate
    normal CDF that no pricing model calls; `norm_ppf` is likewise unused outside
    tests. Either wire them into new analytics (see feature ideas) or document
    them as intentional public utilities. → [#12](https://github.com/rtrimble13/opal/issues/12)

## New feature ideas (evidence-backed)

These are product speculation rather than defects, so they are tracked
separately from the engineering backlog above. Each points at something already
in the tree, and each has now been filed as a GitHub issue for backlog tracking.

- **Two-asset / compound / partial-time-barrier analytics** — the unused
  `bivar_norm_cdf` (normal.hpp:72) is exactly the machinery needed for
  Reiner–Rubinstein compound options, two-asset rainbows, and partial-barrier
  closed forms. The capability is already built and tested; only the pricing
  wrappers are missing. (P2) → [#13](https://github.com/rtrimble13/opal/issues/13)
- **First-class Heston greeks** — the constructive counterpart to finding #1:
  expose analytic or finite-difference Heston greeks (delta/gamma plus
  sensitivities w.r.t. `v0`, `theta`, `xi`, `rho`) so the CLI/Python can report
  real Heston risk rather than a BSM-style vol bump. (P2)
  → [#14](https://github.com/rtrimble13/opal/issues/14)
- **Vol-surface / smile object for SABR** — `sabr_lognormal_vol` and the swaption
  pricer already consume per-strike vols; a thin smile/cube container would let
  `chain`/`scenario` price a whole strike ladder off one calibrated SABR set. (P3)
  → [#15](https://github.com/rtrimble13/opal/issues/15)
- **Calibration routines** — Heston/SABR parameters are inputs everywhere
  (`HestonParams`, `SabrParams`); a least-squares calibrator to a quoted vol
  smile is the obvious next layer for an "institutional" library. (P3)
  → [#16](https://github.com/rtrimble13/opal/issues/16)

## What's done well (preserve these)

- **One generalized BSM core, reused everywhere.** `gbs_price`/`gbs_greeks`
  (black_scholes.hpp) parameterized on cost-of-carry `b` cleanly powers BSM,
  Black-76, FX, Asians (via moment-matched `b_a`,`sigma_a`), and swaptions. This
  is the right abstraction and keeps the analytic surface tiny.
- **Numerically literate engines.** Little-Heston-trap CF (heston.hpp), Rannacher
  start-up in CN (pde.hpp:97), Leisen–Reimer with Peizer–Pratt inversion
  (lattice.hpp:53), antithetic + geometric-control-variate Asian MC
  (monte_carlo.hpp), and a ridge term in the LSMC normal equations (lsmc.hpp:44).
- **Tests that assert the right invariants.** Put–call and in–out parity,
  cross-engine convergence, textbook reference values (Hull, Haug), MC within
  k·standard-errors, and degenerate-limit checks (Heston→BS, SABR flat). This is
  how a pricing library should be tested.
- **Honest, well-written documentation in the code** — the rebate-convention and
  escrowed-dividend comments, and the README's model/engine matrices, are
  accurate to the implementation.
- **Clean packaging story** — header-only INTERFACE target, self-contained sdist
  that bundles `../include` (setup.py), and a wheels workflow that verifies the
  sdist installs standalone.
