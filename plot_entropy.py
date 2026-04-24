"""
plot_entropy.py - Compare entropy along the bottom wall.

Entropy s/s_inf = (p/p_inf) / (rho/rho_inf)^gamma

For ideal gas:
  - Constant across isentropic processes (including expansion fans)
  - Increases across shocks

Usage:
    python plot_entropy.py
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from plot_style import set_latex_style, COLORS, LABELS
from analytical import (get_analytical_wall, post_shock, post_expansion,
                         M_INF, P_INF, RHO_INF, GAMMA, THETA)

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
    x   = data[:, 0].reshape(nx, ny)[:, 0]
    rho = data[:, 2].reshape(nx, ny)[:, 0]
    p   = data[:, 5].reshape(nx, ny)[:, 0]
    return x, rho, p


def entropy(p, rho, p_ref=P_INF, rho_ref=RHO_INF, gamma=GAMMA):
    """Normalized entropy: (p/p_inf) / (rho/rho_inf)^gamma"""
    return (p / p_ref) / (rho / rho_ref) ** gamma


# ── Analytical entropy profile ─────────────────────────────────────
# Region 1 (freestream): s/s_inf = 1
# Region 2 (post leading shock): entropy increases
# Region 3 (post expansion): same as region 2 (isentropic)
# Region 4 (post trailing shock): further increase

# Compute analytical densities for entropy calculation
M2, p2_p1, _ = post_shock(M_INF, THETA)
# rho2/rho1 from oblique shock
Mn1 = M_INF * np.sin(np.radians(39.3137))  # approx beta
rho2_rho1 = (GAMMA + 1) * Mn1**2 / ((GAMMA - 1) * Mn1**2 + 2)
rho2 = RHO_INF * rho2_rho1
p2 = P_INF * p2_p1
s2 = entropy(p2, rho2)

# Region 3: isentropic — same entropy
s3 = s2

# Region 4: another shock
M3, p3_p2 = post_expansion(M2, 2 * THETA)
p3 = p2 * p3_p2
# rho3 from isentropic: rho/rho0 = (p/p0)^(1/gamma)
rho3 = rho2 * (p3_p2) ** (1.0 / GAMMA)

M4, p4_p3, _ = post_shock(M3, THETA)
p4 = p3 * p4_p3
Mn3 = M3 * np.sin(np.radians(35))  # approx
rho4_rho3 = (GAMMA + 1) * Mn3**2 / ((GAMMA - 1) * Mn3**2 + 2)
rho4 = rho3 * rho4_rho3
s4 = entropy(p4, rho4)

# Analytical step function
x_an = np.array([-0.5, -0.001, 0.0, 0.499, 0.5, 0.999, 1.0, 1.5])
s_an = np.array([1.0,    1.0,  s2,  s2,   s3,   s3,   s4,  s4])

# ── Plot ──────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 6))
ax.plot(x_an, s_an, "k--", linewidth=2, label="Analytical")

for fname, key in files:
    try:
        x, rho, p = read_wall(fname)
        if np.any(np.isnan(p)) or np.any(np.isnan(rho)):
            continue
        s = entropy(p, rho)
        ax.plot(x, s, linewidth=1.6, color=COLORS[key], label=LABELS[key])
    except Exception as e:
        print(f"Warning: could not read {fname}: {e}")

for xv in [0.0, 0.5, 1.0]:
    ax.axvline(xv, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)

ax.set_xlabel("$x$")
ax.set_ylabel("Normalized Entropy $s/s_\\infty$")
ax.set_title("Entropy Along Bottom Wall")
ax.legend(loc="best")
ax.grid(True, alpha=0.3)

plt.savefig("entropy_comparison.png")
print(f"Saved entropy_comparison.png")
print(f"Analytical entropy: post-shock={s2:.4f}, post-trailing-shock={s4:.4f}")
