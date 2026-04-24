"""
plot_grid_convergence.py - Compare coarse vs fine grid results for 1st-order scheme.

Requires two solution files (coarse and fine grid runs of the same scheme).

Usage:
    python plot_grid_convergence.py scheme0_coarse.dat scheme0_fine.dat
    python plot_grid_convergence.py                 # uses defaults
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from plot_style import set_latex_style, COLORS
from analytical import get_analytical_wall

set_latex_style()

# ── Inputs ─────────────────────────────────────────────────────────
coarse_file = sys.argv[1] if len(sys.argv) > 1 else "scheme0_coarse.dat"
fine_file   = sys.argv[2] if len(sys.argv) > 2 else "scheme0_fine.dat"


def read_wall(fname):
    with open(fname) as f:
        nx, ny = map(int, f.readline().split())
        data = np.loadtxt(f)
    x = data[:, 0].reshape(nx, ny)[:, 0]
    p = data[:, 5].reshape(nx, ny)[:, 0]
    M = data[:, 6].reshape(nx, ny)[:, 0]
    return nx, ny, x, p, M


# ── Read data ──────────────────────────────────────────────────────
nxc, nyc, xc, pc, Mc = read_wall(coarse_file)
nxf, nyf, xf, pf, Mf = read_wall(fine_file)
x_an, M_an, p_an = get_analytical_wall()

# ── Plot ──────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Mach plot
ax1.plot(x_an, M_an, "k--", linewidth=2, label="Analytical")
ax1.plot(xc, Mc, "o-", linewidth=1.5, markersize=4,
         color=COLORS["coarse"], label=f"Coarse grid (${nxc}\\times{nyc}$)")
ax1.plot(xf, Mf, "s-", linewidth=1.5, markersize=3,
         color=COLORS["fine"], label=f"Fine grid (${nxf}\\times{nyf}$)")
for xv in [0.0, 0.5, 1.0]:
    ax1.axvline(xv, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)
ax1.set_xlabel("$x$")
ax1.set_ylabel("Mach Number")
ax1.set_title("Grid Convergence: Mach Number Along Bottom Wall")
ax1.legend(loc="best")
ax1.grid(True, alpha=0.3)

# Pressure plot
ax2.plot(x_an, p_an, "k--", linewidth=2, label="Analytical")
ax2.plot(xc, pc, "o-", linewidth=1.5, markersize=4,
         color=COLORS["coarse"], label=f"Coarse grid (${nxc}\\times{nyc}$)")
ax2.plot(xf, pf, "s-", linewidth=1.5, markersize=3,
         color=COLORS["fine"], label=f"Fine grid (${nxf}\\times{nyf}$)")
for xv in [0.0, 0.5, 1.0]:
    ax2.axvline(xv, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)
ax2.set_xlabel("$x$")
ax2.set_ylabel("Pressure $P$")
ax2.set_title("Grid Convergence: Pressure Along Bottom Wall")
ax2.legend(loc="best")
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("grid_convergence.png")

# Print summary
print(f"\nPost-expansion Mach (x = 0.7 to 0.9):")
mask_c = (xc > 0.7) & (xc < 0.9)
mask_f = (xf > 0.7) & (xf < 0.9)
print(f"  Coarse ({nxc}x{nyc}): {Mc[mask_c].mean():.4f}")
print(f"  Fine   ({nxf}x{nyf}): {Mf[mask_f].mean():.4f}")
print(f"  Analytical:           2.3717")
print("\nSaved grid_convergence.png")
