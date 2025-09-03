# Kronig–Penney Model (Clojure)

This project provides a simple Kronig–Penney model implementation using the standard dispersion relation for a periodic square potential:

- Total period `L = 2a + 2b` where x=0 is at the center of the barrier
- Barrier width `2a` (region with `V=V0`), centered at x=0
- Well width `2b` (region with `V=0`), split as width `b` on each side of barrier
- Barrier height `V0`
- `mu = ħ^2/(2m)` controls units/scaling (default 1)

Dispersion relation D(E) is used via `cos(k L) = D(E)`. Allowed bands satisfy `|D(E)| ≤ 1`.

## Requirements

- Clojure CLI (tools.deps) or Leiningen

## Run

Generate a dispersion sample CSV and band-edge CSV.

Clojure CLI:
```bash
clojure -M:run -a 0.5 -b 0.25 -V 12.0 -m 1.0 --Emax 40 --steps 10000 -o dispersion.csv
```

Leiningen:
```bash
lein run -a 0.5 -b 0.25 -V 12.0 -m 1.0 --Emax 40 --steps 10000 -o dispersion.csv
```

Outputs:
- `dispersion.csv` with columns: `E, D, allowed, k_minus, k_plus`
- `dispersion-bands.csv` with band intervals: `E_lo, E_hi`

## Notes

- Units are whatever you choose consistently; with `mu=1`, energies are in units of `ħ^2/(2m L^2)` where `L = 2a + 2b`.
- Numerical tolerances and simple bisection are used; increase `--steps` if you need tighter band edges.
