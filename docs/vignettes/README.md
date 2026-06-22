# Opal vignettes

Worked, narrative walk-throughs that take you from zero to confidently using
Opal for pre-trade analysis, pricing and risk management. Each vignette states
a use case, shows the exact command (and the equivalent Python where one
exists), shows the **actual output**, and explains every field and how to act
on it.

Read them in order if you are new; jump straight to a topic if you are not.

| # | Vignette | What you'll learn |
|---|----------|-------------------|
| 1 | [Getting started](01-getting-started.md) | Build, `help`/`version`, your first `price` and `greeks`, reading the output, `%`/date inputs, `-o table/json/csv` |
| 2 | [Pricing and numerical methods](02-pricing-and-methods.md) | European/American, the `bsm`/`black76`/`bachelier` models, choosing a method (`analytic`/`lr`/`crr`/`trinomial`/`pde`/`mc`), the Monte Carlo standard error, discrete dividends |
| 3 | [Exotic payoffs](03-exotic-payoffs.md) | Digitals, gap, barriers (+rebates), Asians, lookbacks, partial-time barriers, two-asset (exchange/rainbow/correlation) and compound options |
| 4 | [Stochastic volatility](04-stochastic-volatility.md) | Heston pricing and first-class Heston greeks; the SABR smile |
| 5 | [Greeks and implied vol](05-greeks-and-implied-vol.md) | The full risk report (incl. vanna/volga/charm), implied volatility, option chains |
| 6 | [Scenario analysis](06-scenario-analysis.md) | Reading a spot × vol P&L grid |
| 7 | [Portfolio risk](07-portfolio-risk.md) | Pricing and risk-aggregating a book from CSV |
| 8 | [Interest-rate options](08-interest-rate-options.md) | Caps/floors/caplets/floorlets, swaptions, ZCB options; multi-curve OIS discounting; lognormal vs normal vols |
| 9 | [The Python API](09-python-api.md) | Pricers, greeks, implied vol, engines, `mc_custom`, LSMC, Heston/SABR/Hull-White, `DiscountCurve`, and an end-to-end pre-trade workflow |

## Conventions used throughout

- **Rates and dividend yields** are continuously compounded decimals; the CLI
  also accepts `%` (`-r 4%` ≡ `-r 0.04`).
- **Expiries** are year fractions (`-T 0.5`) or `YYYY-MM-DD` dates
  (ACT/365F from today).
- **Library greeks** are per unit of the underlying variable; the CLI prints
  trader conventions — vega per vol point, theta per calendar day, rho per 1%
  rate. In Python you scale yourself (`g.vega / 100`, `g.theta / 365`).
- Numbers below are from the current build (`opal 0.2.0`). Monte Carlo figures
  carry a reported standard error and will vary at the last digits with
  `--paths`/`--seed`.
