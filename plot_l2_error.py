"""
plot_l2_error.py - Compute and plot L2 error of each scheme vs analytical.

Usage:
    python plot_l2_error.py
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from plot_style import set_latex_style, COLORS, LABELS
from analytical import get_analytical_wall

set_latex_style()

files = [
    ("scheme0_solution.dat", "1st_order"),
    ("scheme1_solution.dat", "2nd_order"),
    ("scheme2_solution.dat", "tvd"),
]


def read_wall(fname):
    with open(fname) as f:
        nx, ny = map(int, f.readline().split())
        data = np.loadtxt(f)
    x = data[:, 0].reshape(nx, ny)[:, 0]
    p = data[:, 5].reshape(nx, ny)[:, 0]
    M = data[:, 6].reshape(nx, ny)[:, 0]
    return x, p, M


def analytical_at(x_query, x_an, f_an):
    """Interpolate the analytical step function at query points."""
    return np.interp(x_query, x_an, f_an)


# ── Compute L2 errors ──────────────────────────────────────────────
x_an, M_an, p_an = get_analytical_wall()

L2_p = {}
L2_M = {}
for fname, key in files:
    try:
        x, p, M = read_wall(fname)
        if np.any(np.isnan(p)):
            L2_p[key] = np.nan
            L2_M[key] = np.nan
            continue
        M_exact = analytical_at(x, x_an, M_an)
        p_exact = analytical_at(x, x_an, p_an)
        L2_p[key] = np.sqrt(np.mean((p - p_exact) ** 2))
        L2_M[key] = np.sqrt(np.mean((M - M_exact) ** 2))
    except Exception as e:
        print(f"Warning: could not read {fname}: {e}")
        L2_p[key] = np.nan
        L2_M[key] = np.nan

# ── Print table ────────────────────────────────────────────────────
print("=" * 55)
print(f"{'Scheme':20} {'L2(pressure)':>15} {'L2(Mach)':>15}")
print("-" * 55)
for key in ["1st_order", "2nd_order", "tvd"]:
    print(f"{LABELS[key]:20} {L2_p[key]:>15.6f} {L2_M[key]:>15.6f}")
print("=" * 55)

# ── Bar chart ──────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

keys = ["1st_order", "2nd_order", "tvd"]
labels_short = [LABELS[k] for k in keys]
cols   = [COLORS[k] for k in keys]
vals_p = [L2_p[k] for k in keys]
vals_M = [L2_M[k] for k in keys]

ax1.bar(labels_short, vals_p, color=cols, edgecolor="black")
ax1.set_ylabel("$L_2$ error")
ax1.set_title("$L_2$ Error --- Pressure")
ax1.tick_params(axis="x", labelrotation=15)
for i, v in enumerate(vals_p):
    if not np.isnan(v):
        ax1.text(i, v, f"{v:.4f}", ha="center", va="bottom", fontsize=10)

ax2.bar(labels_short, vals_M, color=cols, edgecolor="black")
ax2.set_ylabel("$L_2$ error")
ax2.set_title("$L_2$ Error --- Mach Number")
ax2.tick_params(axis="x", labelrotation=15)
for i, v in enumerate(vals_M):
    if not np.isnan(v):
        ax2.text(i, v, f"{v:.4f}", ha="center", va="bottom", fontsize=10)

plt.tight_layout()
plt.savefig("l2_error_comparison.png")
print("\nSaved l2_error_comparison.png")
