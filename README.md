# Kronig–Penney Model (Clojure)

A comprehensive implementation of the Kronig–Penney model in Clojure, supporting both 1D and 2D periodic potentials.

## Features

- **1D Models**:
  - Simple Kronig–Penney model
  - Extended KP model (U1-well-U2-well structure)
  - Multilayer transfer-matrix model
- **2D Models**:
  - 2D separable Kronig–Penney model
  - 2D band structure calculations
  - 2D Brillouin zone sampling
- **Comprehensive CLI**: Single interface for all model variants
- **Multiple output formats**: CSV data and band edge information

## Project Structure

```
kp/
├── src/kp/
│   ├── core.clj           # Main CLI and utilities
│   ├── common.clj         # Shared mathematical functions
│   ├── model_1d.clj       # 1D Kronig-Penney implementations
│   └── model_2d.clj       # 2D Kronig-Penney implementations
├── test/kp/
│   ├── model_1d_test.clj  # 1D model tests
│   └── model_2d_test.clj  # 2D model tests
├── examples/
│   ├── 1d_examples/       # 1D usage examples
│   └── 2d_examples/       # 2D usage examples
└── README.md
```

## Quick Start

### 1D Simple Kronig–Penney Model

Basic periodic square potential:
- Total period `L = 2a + 2b` where x=0 is at the center of the barrier
- Barrier width `2a` (region with `V=V0`), centered at x=0
- Well width `2b` (region with `V=0`), split as width `b` on each side of barrier
- Barrier height `V0`
- `mu = ħ^2/(2m)` controls units/scaling (default 1)

Dispersion relation D(E) is used via `cos(k L) = D(E)`. Allowed bands satisfy `|D(E)| ≤ 1`.

## Extended Kronig-Penney Model

The project also implements an extended Kronig-Penney model that accounts for crystal inhomogeneity through alternating barriers of different heights:

### Structure
- **Period**: `L = 4a + 4b`
- **Structure**: U₁ barrier (width 2a) - well (width 2b) - U₂ barrier (width 2a) - well (width 2b)
- **Coordinate system**: x = 0 at center of U₁ barrier
- **Parameters**:
  - `a`: Half barrier width (each barrier width = 2a)
  - `b`: Well width on each side (total well width = 2b)
  - `U₁`: Height of first barrier
  - `U₂`: Height of second barrier (U₂ ≥ U₁)
  - `mu = ħ^2/(2m)`: Mass parameter (default 1)

### Physical Interpretation
- **U₁ barriers**: Describe main periodic component of crystal potential
- **U₂ barriers**: Model local perturbations from impurities or point defects
- **Special case**: When U₁ = U₂, the extended model reproduces the simple KP model exactly

### Usage
```bash
# Extended model with different barrier heights
clojure -M:run -a 0.3 -b 0.2 --U1 8.0 --U2 12.0 --extended --Emax 40 -o extended.csv

# Extended model reproducing simple KP (U1 = U2)
clojure -M:run -a 0.3 -b 0.2 --U1 10.0 --U2 10.0 --extended --Emax 40 -o extended-equal.csv
```

## Requirements

- Clojure CLI (tools.deps) or Leiningen

## Run

Generate a dispersion sample CSV and band-edge CSV.

### Simple KP Model
Clojure CLI:
```bash
clojure -M:run -a 0.5 -b 0.25 -V 12.0 -m 1.0 --Emax 40 --steps 10000 -o dispersion.csv
```

Leiningen:
```bash
lein run -a 0.5 -b 0.25 -V 12.0 -m 1.0 --Emax 40 --steps 10000 -o dispersion.csv
```

### Extended KP Model
```bash
# Extended model with different barrier heights
clojure -M:run -a 0.3 -b 0.2 --U1 8.0 --U2 12.0 --extended --Emax 40 -o extended.csv

# Extended model reproducing simple KP (U1 = U2)
clojure -M:run -a 0.3 -b 0.2 --U1 10.0 --U2 10.0 --extended --Emax 40 -o extended-equal.csv
```

### Multilayer Model
```bash
clojure -M:run --layers "b:0.4,U:12:0.2,b:0.4" --Emax 40 --steps 10000 -o multi.csv
```

### 2D Kronig–Penney Model
```bash
# Basic 2D calculation with 20x20 k-point grid
clojure -M:run --2d -a 0.5 -b 0.25 -V 12.0 --Lx 1.0 --Ly 1.0 --nx 20 --ny 20 -o 2d-kp.csv

# Higher resolution k-grid
clojure -M:run --2d -a 0.5 -b 0.25 -V 12.0 --Lx 1.0 --Ly 1.0 --nx 50 --ny 50 -o 2d-kp-hires.csv
```

## Output Formats

### 1D Models
- `dispersion.csv` with columns: `E, D, allowed, k_minus, k_plus`
- `dispersion-bands.csv` with band intervals: `E_lo, E_hi`

### 2D Models
- `dispersion.csv` with columns: `E, kx, ky, D, allowed, k_magnitude`
- `dispersion-surface.csv` with dispersion surface: `kx, ky, D`

## Notes

- Units are whatever you choose consistently; with `mu=1`, energies are in units of `ħ^2/(2m L^2)` where `L = 2a + 2b` for simple KP and `L = 4a + 4b` for extended KP.
- Numerical tolerances and simple bisection are used; increase `--steps` if you need tighter band edges.
- The extended model uses transfer matrix formalism for U₁ ≠ U₂ and falls back to simple KP formula when U₁ = U₂.
- All models support the same output format: dispersion data and band edge intervals.
