# 7. Portfolio risk

**Use case:** you hold a book of positions and need one consolidated view —
value and net greeks — to see your aggregate exposure and what to hedge.

## A book in CSV

`portfolio` reads a CSV where each row is a position. The header names the
columns; quantity can be negative (a short). Here is the bundled
`examples/book.csv`:

```text
instrument,type,quantity,spot,strike,vol,rate,div,expiry,barrier,rebate,method
european,call,100,100,105,0.22,0.04,0.01,0.5,,,
european,put,-50,100,95,0.25,0.04,0.01,0.5,,,
american,put,200,50,55,0.30,0.03,0.00,1.0,,,
barrier-up-out,call,150,100,100,0.25,0.05,0.00,1.0,120,0,
asian-arith,call,75,100,100,0.30,0.05,0.00,1.0,,,analytic
```

Each row carries its **own** market data (spot, vol, rate, …), so a book can mix
underlyings, models and instruments. Empty cells fall back to sensible defaults;
the optional `method` column forces a numerical method per row. Expiries accept
year fractions or `YYYY-MM-DD` dates, and `vol`/`rate`/`div` accept `%`.

## Pricing and aggregating

```sh
opal portfolio --file examples/book.csv
```

```text
#      instrument  type   qty   price    value    delta    gamma    vega  theta
-------------------------------------------------------------------------------
1        european  call   100  4.7122   471.22    44.23   2.5270   27.80  -1.99
2        european   put   -50  4.0621  -203.11    16.00  -1.0084  -12.61   0.71
3        american   put   200  8.2509  1650.18  -110.58   5.4880   38.88  -1.21
4  barrier-up-out  call   150  0.6913   103.70    -3.12  -0.3444  -10.36   0.35
5     asian-arith  call    75  7.9925   599.44    43.26   1.6244   16.69  -0.91
6    digital-cash  call  1000  0.1814   181.41    26.27   2.1098   10.55  -1.43
            TOTAL                      2802.85    16.05  10.3963   70.95  -4.48
```

Reading the table:

- `price` is per unit; `value` = `qty × price` (negative for the short put,
  row 2).
- The greek columns are **share-equivalent and quantity-weighted** — i.e.
  already multiplied by `qty` — so they sum down to a portfolio total. `vega` is
  per vol point and `theta` per calendar day (the footer states the convention).
- The **TOTAL** row is the book's net risk. Here the desk is:
  - net long **value** 2802.85;
  - slightly net long **delta** (+16) — close to delta-neutral, dominated by the
    long American puts (−110) offsetting the long calls;
  - net long **gamma** (+10.4) and **vega** (+71) — long convexity and long vol;
  - paying **theta** (−4.48/day) — the cost of carrying that optionality.

That one line tells you the hedges: sell ~16 units of underlying to flatten
delta, and you are structurally long gamma/vega financed by time decay.

## Error reporting

Malformed cells are reported with their row and column so a bad book is easy to
fix, e.g.:

```text
opal: portfolio row 3, column 'vol': expected a number, got 'oops'
```

Add `-o csv` (or `-o json`) to feed the per-position breakdown into a risk
system. Next: [interest-rate options](08-interest-rate-options.md).
