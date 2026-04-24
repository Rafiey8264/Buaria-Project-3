"""
plot_convergence.py - Plot convergence history comparing all three schemes.

Usage:
    python plot_convergence.py
    python plot_convergence.py scheme0_convergence.dat scheme1_convergence.dat scheme2_convergence.dat

If no arguments, uses the default files scheme{0,1,2}_convergence.dat.
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from plot_style import set_latex_style, COLORS, LABELS

set_latex_style()

# ── Inputs ──────────────────────────────────────────────────────────
files = sys.argv[1:4] if len(sys.argv) >= 4 else [
    "scheme0_convergence.dat",
    "scheme1_convergence.dat",
    "scheme2_convergence.dat",
]
scheme_keys = ["1st_order", "2nd_order", "tvd"]

# ── Plot ────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(11, 6))

for fname, key in zip(files, scheme_keys):
    try:
        data = np.loadtxt(fname)
        iters = data[:, 0]
        R     = data[:, 1]
        ax.loglog(iters, R,
                  color=COLORS[key],
                  label=LABELS[key],
                  linewidth=1.8)
    except Exception as e:
        print(f"Warning: could not read {fname}: {e}")

# Reference line at professor's expected convergence level
ax.axhline(1e-3, color="gray", linestyle="--", linewidth=1,
           alpha=0.6, label="$R = 10^{-3}$ (reference)")
ax.axhline(1e-6, color="gray", linestyle=":", linewidth=1,
           alpha=0.6, label="$R = 10^{-6}$ (tolerance)")

ax.set_xlabel("Iteration")
ax.set_ylabel("Relative Error $R$")
ax.set_title("Convergence History Comparison")
ax.legend(loc="lower left")
ax.grid(True, which="both", alpha=0.3)

plt.savefig("convergence_comparison.png")
print("Saved convergence_comparison.png")
