# 2D Kronig-Penney Model Examples

This directory contains examples for the 2D Kronig-Penney model.

## Basic 2D KP Model

2D separable Kronig-Penney model with square lattice:

```bash
# Basic 2D calculation with 20x20 k-point grid
lein run --2d -a 0.5 -b 0.25 -V 12.0 --Lx 1.0 --Ly 1.0 --nx 20 --ny 20 -o 2d-kp.csv

# Higher resolution k-grid
lein run --2d -a 0.5 -b 0.25 -V 12.0 --Lx 1.0 --Ly 1.0 --nx 50 --ny 50 -o 2d-kp-hires.csv
```

## Custom k-space Ranges

Specify custom k-space sampling ranges:

```bash
# Sample only first quadrant of Brillouin zone
lein run --2d -a 0.5 -b 0.25 -V 12.0 --Lx 1.0 --Ly 1.0 --nx 30 --ny 30 \
  --kx-min 0.0 --kx-max 1.57 --ky-min 0.0 --ky-max 1.57 -o 2d-first-quadrant.csv

# Sample along diagonal of Brillouin zone
lein run --2d -a 0.5 -b 0.25 -V 12.0 --Lx 1.0 --Ly 1.0 --nx 50 --ny 1 \
  --kx-min 0.0 --kx-max 3.14 --ky-min 0.0 --ky-max 0.0 -o 2d-diagonal.csv
```

## Different Lattice Constants

Explore effects of different lattice constants:

```bash
# Rectangular lattice (Lx ≠ Ly)
lein run --2d -a 0.5 -b 0.25 -V 12.0 --Lx 2.0 --Ly 1.0 --nx 30 --ny 30 -o 2d-rectangular.csv

# Square lattice with different period
lein run --2d -a 0.3 -b 0.2 -V 8.0 --Lx 0.8 --Ly 0.8 --nx 40 --ny 40 -o 2d-small-period.csv
```

## Output Files

Each 2D example generates:
- `*.csv`: 2D band structure data with columns [E, kx, ky, D, allowed, k_magnitude]
- `*-surface.csv`: Dispersion surface at fixed energy with columns [kx, ky, D]

## Visualization

The 2D data can be visualized using:

### Python (matplotlib)
```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data
df = pd.read_csv('2d-kp.csv')

# Create 3D surface plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(df['kx'], df['ky'], df['D'], cmap='viridis')
plt.show()
```

### R (ggplot2)
```r
library(ggplot2)
library(plotly)

# Load data
data <- read.csv('2d-kp.csv')

# Create 3D surface
plot_ly(data, x = ~kx, y = ~ky, z = ~D, type = 'scatter3d', mode = 'markers')
```

### MATLAB/Octave
```matlab
% Load data
data = csvread('2d-kp.csv', 1, 0);  % Skip header
kx = data(:, 2);
ky = data(:, 3);
D = data(:, 4);

% Create surface plot
surf(reshape(kx, [20, 20]), reshape(ky, [20, 20]), reshape(D, [20, 20]));
```

## Physical Interpretation

- **E**: Energy (in units of ħ²/(2mL²))
- **kx, ky**: Wave vector components (in units of π/L)
- **D**: Dispersion relation value (|D| ≤ 1 for allowed states)
- **allowed**: Boolean indicating if state is in allowed band
- **k_magnitude**: Magnitude of 2D wave vector

