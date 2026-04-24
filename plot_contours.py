"""
plot_contours.py - Generate 4 contour plots (p, u, v, Mach) for one scheme.

Usage:
    python plot_contours.py scheme0_solution.dat "1st Order"
    python plot_contours.py scheme1_solution.dat "2nd Order"
    python plot_contours.py scheme2_solution.dat "TVD"
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from plot_style import set_latex_style

set_latex_style()

# ── Inputs ──────────────────────────────────────────────────────────
if len(sys.argv) < 3:
    print("Usage: python plot_contours.py <solution.dat> <label>")
    sys.exit(1)

infile = sys.argv[1]
label  = sys.argv[2]   # e.g., "1st Order", "2nd Order", "TVD"

# Make a filename-safe prefix from the label
prefix = label.replace(" ", "_")

# ── Read solution ───────────────────────────────────────────────────
with open(infile) as f:
    nx, ny = map(int, f.readline().split())
    data = np.loadtxt(f)

x   = data[:, 0].reshape(nx, ny)
y   = data[:, 1].reshape(nx, ny)
rho = data[:, 2].reshape(nx, ny)
u   = data[:, 3].reshape(nx, ny)
v   = data[:, 4].reshape(nx, ny)
p   = data[:, 5].reshape(nx, ny)
M   = data[:, 6].reshape(nx, ny)

# ── Plot each variable ──────────────────────────────────────────────
variables = [
    (p, "Pressure",    "$P$",            "RdBu_r"),
    (u, "u-velocity",  "$u$",            "RdBu_r"),
    (v, "v-velocity",  "$v$",            "RdBu_r"),
    (M, "Mach Number", "Mach Number",    "jet"),
]

for field, vname, colorbar_label, cmap in variables:
    fig, ax = plt.subplots(figsize=(12, 5))
    cf = ax.contourf(x, y, field, levels=50, cmap=cmap)
    cbar = plt.colorbar(cf, ax=ax)
    cbar.set_label(colorbar_label)

    # Fill the diamond airfoil with white
    airfoil_x = [0.0, 0.5, 1.0, 0.0]
    airfoil_y = [0.0, 0.5 * np.tan(np.radians(10.0)), 0.0, 0.0]
    ax.fill(airfoil_x, airfoil_y, color="white", edgecolor="black", linewidth=0.8)

    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(f"{vname} Contour --- {label}")
    ax.set_aspect("equal")
    ax.grid(False)

    # File-safe variable name
    safe_vname = vname.replace(" ", "_").replace("-", "_")
    outfile = f"{prefix}_{safe_vname}_contour.png"
    plt.savefig(outfile)
    plt.close()
    print(f"Saved {outfile}")

print(f"\nAll contours for '{label}' generated.")
