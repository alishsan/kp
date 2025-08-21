# Kronig–Penney Model (Clojure)

This project provides a simple Kronig–Penney model implementation using the standard dispersion relation for a periodic square potential:

- Lattice period `a`
- Well width `b` (region with `V=0`), barrier width `c=a-b` (region with `V=V0`)
- Barrier height `V0`
- `mu = ħ^2/(2m)` controls units/scaling (default 1)

Dispersion relation D(E) is used via `cos(k a) = D(E)`. Allowed bands satisfy `|D(E)| ≤ 1`.

## Requirements

- Clojure CLI (tools.deps) or Leiningen

## Run

Generate a dispersion sample CSV and band-edge CSV.

Clojure CLI:
```bash
clojure -M:run -a 1.0 -b 0.4 -V 12.0 -m 1.0 --Emax 40 --steps 10000 -o dispersion.csv
```

Leiningen:
```bash
lein run -a 1.0 -b 0.4 -V 12.0 -m 1.0 --Emax 40 --steps 10000 -o dispersion.csv
```

Outputs:
- `dispersion.csv` with columns: `E, D, allowed, k_minus, k_plus`
- `dispersion-bands.csv` with band intervals: `E_lo, E_hi`

## Notes

- Units are whatever you choose consistently; with `mu=1`, energies are in units of `ħ^2/(2m a^2)` if you choose `a=1`.
- Numerical tolerances and simple bisection are used; increase `--steps` if you need tighter band edges.
