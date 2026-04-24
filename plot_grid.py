"""
plot_grid.py - Visualize the computational grid.

Usage:
    python plot_grid.py                    # uses grid_output.dat
    python plot_grid.py grid_output.dat    # specify input file
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from plot_style import set_latex_style

set_latex_style()

# ── Input file ──────────────────────────────────────────────────────
infile = sys.argv[1] if len(sys.argv) > 1 else "grid_output.dat"

# ── Read grid ───────────────────────────────────────────────────────
with open(infile) as f:
    nx, ny = map(int, f.readline().split())
    data = np.loadtxt(f)

x = data[:, 0].reshape(nx, ny)
y = data[:, 1].reshape(nx, ny)

# ── Plot ────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 6))

# xi-lines (varying i, fixed j)
for j in range(ny):
    ax.plot(x[:, j], y[:, j], "k-", linewidth=0.4, alpha=0.6)

# eta-lines (varying j, fixed i)
for i in range(nx):
    ax.plot(x[i, :], y[i, :], "k-", linewidth=0.4, alpha=0.6)

ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_title(f"Computational Grid (${nx} \\times {ny}$)")
ax.set_aspect("equal")
ax.grid(False)

plt.savefig("grid.png")
print(f"Saved grid.png ({nx} x {ny})")
