# 6. Scenario analysis

**Use case:** greeks are a local, first-order view. Before putting on a trade
you want the full **non-linear** P&L across a range of spot and vol moves — a
stress grid.

## A spot × vol P&L grid

`scenario` reprices the trade across a grid of spot shifts (rows, in %) and vol
shifts (columns, in vol points) and shows the P&L versus the base price:

```sh
opal scenario -S 100 -K 100 -T 0.5 -r 4% -v 25% --spot-range -10:10:5 --vol-range -5:5:5
```

```text
P&L vs base price 8.0080 (rows: spot, cols: vol)
spot\vol    20.0%    25.0%    30.0%
-----------------------------------
   90.00  -5.7836  -4.6419  -3.4454
   95.00  -3.9451  -2.6146  -1.2775
  100.00  -1.3809   0.0000   1.3824
  105.00   1.8669   3.1645   4.5017
  110.00   5.6875   6.8072   8.0283
```

How to read it:

- The **base price** (8.0080) is the trade at today's spot (100) and vol (25%).
  Every cell is the **change** in value from that base — a P&L, not a price.
  (Use `--absolute` to print prices instead of P&L.)
- **Rows** are spot scenarios. `--spot-range -10:10:5` means −10% to +10% in 5%
  steps, so spots 90, 95, 100, 105, 110.
- **Columns** are vol scenarios. `--vol-range -5:5:5` means −5 to +5 vol points
  in steps of 5, so the 25% base vol becomes 20%, 25%, 30%.
- The **center cell is 0.0000** by construction (base spot, base vol).

What the grid tells you that greeks alone don't:

- **Convexity (gamma):** moving up the spot column the gains accelerate
  (+1.87 at +5%, +5.69 at +10% in the base-vol column) and the losses
  decelerate — the long-gamma signature. A pure delta estimate would be
  symmetric; the grid shows the curvature.
- **Vega and cross-effects:** moving across a row, more vol always helps this
  long option, and the spot×vol corner cells (e.g. +10% spot **and** +5 vol →
  +8.03) show the combined stress, capturing the vanna interaction a single
  greek misses.
- **Worst case:** the bottom-left region (spot down, vol down) is the most
  painful — here −5.78 at −10%/−5 — exactly the scenario to size against.

`scenario` works for any equity instrument, including under `--model heston`
(the vol axis shifts the variance level, per the Heston vega convention). Add
`-o csv` to drop the grid into a risk report.

Next: [portfolio risk](07-portfolio-risk.md).
