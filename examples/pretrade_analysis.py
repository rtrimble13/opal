"""Pre-trade analysis example using the Opal Python API.

Prices a proposed 6-month 105% call overwrite, checks the vol surface
implied by market quotes, and sizes the position by scenario risk.

Run from the repo root after building the Python package:
    pip install ./python
    python examples/pretrade_analysis.py
"""
import opal

SPOT = 100.0
RATE = 0.04
DIV = 0.015
EXPIRY = 0.5

print("=== 1. Implied vols from market quotes ===")
quotes = {95: 8.10, 100: 5.32, 105: 3.22, 110: 1.81}
smile = {}
for strike, px in quotes.items():
    iv = opal.implied_vol("call", price=px, spot=SPOT, strike=strike,
                          expiry=EXPIRY, rate=RATE, div=DIV)
    smile[strike] = iv
    print(f"  K={strike:>3}  price={px:>5.2f}  implied vol={iv:.2%}")

print("\n=== 2. Risk report for the 105 call ===")
vol = smile[105]
g = opal.bs_greeks("call", spot=SPOT, strike=105, expiry=EXPIRY, rate=RATE,
                   div=DIV, vol=vol)
print(f"  price  {g.price:8.4f}")
print(f"  delta  {g.delta:8.4f}")
print(f"  gamma  {g.gamma:8.4f}")
print(f"  vega   {g.vega / 100:8.4f}  (per vol pt)")
print(f"  theta  {g.theta / 365:8.4f}  (per day)")

print("\n=== 3. Compare against an American-style listed contract ===")
amer = opal.binomial_price("call", "american", spot=SPOT, strike=105,
                           expiry=EXPIRY, rate=RATE, div=DIV, vol=vol)
print(f"  european {g.price:.4f} vs american {amer:.4f} "
      f"(early-exercise premium {amer - g.price:.4f})")

print("\n=== 4. Scenario P&L for 100 contracts (x100 multiplier), short ===")
qty = -100 * 100
base = g.price
print(f"{'spot':>8} | " + " | ".join(f"vol {v:+d}pt" for v in (-5, 0, 5)))
for ds in (-10, -5, 0, 5, 10):
    s = SPOT * (1 + ds / 100)
    row = []
    for dv in (-5, 0, 5):
        p = opal.bs_price("call", spot=s, strike=105, expiry=EXPIRY, rate=RATE,
                          div=DIV, vol=vol + dv / 100)
        row.append(f"{qty * (p - base):>10,.0f}")
    print(f"{s:8.2f} | " + " | ".join(row))

print("\n=== 5. Cheapen with a knock-out: 105 call, 120 KO barrier ===")
ko = opal.barrier_price("call", "up-out", spot=SPOT, strike=105, barrier=120,
                        expiry=EXPIRY, rate=RATE, div=DIV, vol=vol)
print(f"  vanilla {g.price:.4f} vs up-and-out {ko:.4f} "
      f"({(1 - ko / g.price):.0%} cheaper)")
