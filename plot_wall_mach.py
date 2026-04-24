"""
plot_wall_mach.py - Compare numerical Mach along the bottom wall vs analytical.

Usage:
    python plot_wall_mach.py
    python plot_wall_mach.py scheme0_solution.dat scheme1_solution.dat scheme2_solution.dat
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from plot_style import set_latex_style, COLORS, LABELS
from analytical import get_analytical_wall

set_latex_style()

# ── Input files ─────────────────────────────────────────────────────
files = sys.argv[1:4] if len(sys.argv) >= 4 else [
    "scheme0_solution.dat",
    "scheme1_solution.dat",
    "scheme2_solution.dat",
]
scheme_keys = ["1st_order", "2nd_order", "tvd"]


def read_wall(fname):
    """Read solution file and extract x, M along bottom wall (j=0)."""
    with open(fname) as f:
        nx, ny = map(int, f.readline().split())
        data = np.loadtxt(f)
    x = data[:, 0].reshape(nx, ny)[:, 0]
    M = data[:, 6].reshape(nx, ny)[:, 0]
    return x, M


# ── Analytical ──────────────────────────────────────────────────────
x_an, M_an, _ = get_analytical_wall()

# ── Plot ────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 6))

ax.plot(x_an, M_an, "k--", linewidth=2, label="Analytical")

for fname, key in zip(files, scheme_keys):
    try:
        x, M = read_wall(fname)
        if not np.any(np.isnan(M)):
            ax.plot(x, M, linewidth=1.6, color=COLORS[key], label=LABELS[key])
    except Exception as e:
        print(f"Warning: could not read {fname}: {e}")

# Vertical guidelines at airfoil features
for xv in [0.0, 0.5, 1.0]:
    ax.axvline(xv, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)

ax.set_xlabel("$x$")
ax.set_ylabel("Mach Number")
ax.set_title("Mach Number Along Bottom Wall: Numerical vs Analytical")
ax.legend(loc="best")
ax.grid(True, alpha=0.3)

plt.savefig("wall_mach_comparison.png")
print("Saved wall_mach_comparison.png")
