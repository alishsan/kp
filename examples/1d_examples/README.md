# 1D Kronig-Penney Model Examples

This directory contains examples for the 1D Kronig-Penney model variants.

## Simple KP Model

Basic Kronig-Penney model with periodic square potential:

```bash
# Generate dispersion data
lein run -a 0.5 -b 0.25 -V 12.0 --Emax 40 --steps 10000 -o simple-kp.csv

# Plot the results (requires external plotting tool)
```

## Extended KP Model

Extended model with different barrier heights (U1 ≠ U2):

```bash
# Different barrier heights
lein run -a 0.3 -b 0.2 --U1 8.0 --U2 12.0 --extended --Emax 40 -o extended-kp.csv

# Equal barrier heights (reproduces simple KP)
lein run -a 0.3 -b 0.2 --U1 10.0 --U2 10.0 --extended --Emax 40 -o extended-equal.csv
```

## Multilayer Model

Arbitrary multilayer structure using transfer matrices:

```bash
# Three-layer system: well-barrier-well
lein run --layers "b:0.4,U:12:0.2,b:0.4" --Emax 40 --steps 10000 -o multilayer.csv

# Complex multilayer: barrier-well-barrier-well
lein run --layers "U:12:0.2,b:0.4,U:8:0.2,b:0.4" --Emax 40 --steps 10000 -o complex-multilayer.csv
```

## Output Files

Each example generates:
- `*.csv`: Dispersion data with columns [E, D, allowed, k_minus, k_plus]
- `*-bands.csv`: Band edge intervals [E_lo, E_hi]

## Visualization

The CSV files can be imported into plotting software like:
- Python (matplotlib, pandas)
- R (ggplot2)
- MATLAB
- GNU Octave
- Any spreadsheet application

