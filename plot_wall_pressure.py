"""
plot_wall_pressure.py - Compare numerical pressure along the bottom wall vs analytical.

Usage:
    python plot_wall_pressure.py
    python plot_wall_pressure.py scheme0_solution.dat scheme1_solution.dat scheme2_solution.dat
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from plot_style import set_latex_style, COLORS, LABELS
from analytical import get_analytical_wall

set_latex_style()

files = sys.argv[1:4] if len(sys.argv) >= 4 else [
    "scheme0_solution.dat",
    "scheme1_solution.dat",
    "scheme2_solution.dat",
]
scheme_keys = ["1st_order", "2nd_order", "tvd"]


def read_wall(fname):
    with open(fname) as f:
        nx, ny = map(int, f.readline().split())
        data = np.loadtxt(f)
    x = data[:, 0].reshape(nx, ny)[:, 0]
    p = data[:, 5].reshape(nx, ny)[:, 0]
    return x, p


x_an, _, p_an = get_analytical_wall()

fig, ax = plt.subplots(figsize=(12, 6))
ax.plot(x_an, p_an, "k--", linewidth=2, label="Analytical")

for fname, key in zip(files, scheme_keys):
    try:
        x, p = read_wall(fname)
        if not np.any(np.isnan(p)):
            ax.plot(x, p, linewidth=1.6, color=COLORS[key], label=LABELS[key])
    except Exception as e:
        print(f"Warning: could not read {fname}: {e}")

for xv in [0.0, 0.5, 1.0]:
    ax.axvline(xv, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)

ax.set_xlabel("$x$")
ax.set_ylabel("Pressure $P$")
ax.set_title("Pressure Along Bottom Wall: Numerical vs Analytical")
ax.legend(loc="best")
ax.grid(True, alpha=0.3)

plt.savefig("wall_pressure_comparison.png")
print("Saved wall_pressure_comparison.png")
