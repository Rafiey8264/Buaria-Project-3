"""
analytical.py - Compute analytical solution for Mach 2 flow over diamond airfoil.

Returns step-function arrays for comparison with numerical results along
the bottom wall.

Usage:
    from analytical import get_analytical_wall
    x_an, M_an, p_an = get_analytical_wall()
"""
import numpy as np
from scipy.optimize import brentq

# ── Constants ───────────────────────────────────────────────────────
GAMMA   = 1.4
M_INF   = 2.0
RHO_INF = 4.0
A_INF   = 0.5                           # Mach=2 → V=1, so a=V/M=0.5
P_INF   = RHO_INF * A_INF**2 / GAMMA    # from a^2 = gamma*p/rho
THETA   = np.radians(10.0)              # half-angle of diamond


# ── Oblique shock relations ────────────────────────────────────────
def oblique_beta(M, theta, gamma=GAMMA):
    """Find weak-shock oblique-shock angle beta from theta-beta-M relation."""
    mu = np.arcsin(1.0 / M)  # Mach angle (lower bound on beta)
    def eq(beta):
        num = M**2 * np.sin(beta)**2 - 1
        den = M**2 * (gamma + np.cos(2 * beta)) + 2
        return np.tan(theta) - 2 * (1.0 / np.tan(beta)) * num / den
    for deg in range(int(np.degrees(mu)) + 1, 90):
        if eq(np.radians(deg)) < 0:
            return brentq(eq, np.radians(deg - 1), np.radians(deg))
    raise RuntimeError("No weak-shock solution found")


def post_shock(M1, theta, gamma=GAMMA):
    """Return (M2, p2/p1) across an oblique shock with deflection theta."""
    beta = oblique_beta(M1, theta, gamma)
    Mn1 = M1 * np.sin(beta)
    Mn2 = np.sqrt((1 + 0.5 * (gamma - 1) * Mn1**2) /
                  (gamma * Mn1**2 - 0.5 * (gamma - 1)))
    M2 = Mn2 / np.sin(beta - theta)
    p_ratio = 1 + 2 * gamma / (gamma + 1) * (Mn1**2 - 1)
    return M2, p_ratio, np.degrees(beta)


def prandtl_meyer(M, gamma=GAMMA):
    """Prandtl–Meyer angle (in radians)."""
    sq = np.sqrt((gamma - 1) / (gamma + 1) * (M**2 - 1))
    return (np.sqrt((gamma + 1) / (gamma - 1)) * np.arctan(sq)
            - np.arctan(np.sqrt(M**2 - 1)))


def post_expansion(M1, delta_theta, gamma=GAMMA):
    """Return (M2, p2/p1) across an expansion fan of deflection delta_theta."""
    nu1 = prandtl_meyer(M1, gamma)
    nu2 = nu1 + delta_theta
    M2 = brentq(lambda M: prandtl_meyer(M, gamma) - nu2, 1.001, 50.0)
    # Isentropic pressure ratio
    p_ratio = ((1 + 0.5 * (gamma - 1) * M1**2) /
               (1 + 0.5 * (gamma - 1) * M2**2)) ** (gamma / (gamma - 1))
    return M2, p_ratio


# ── Full wall solution ──────────────────────────────────────────────
def get_analytical_wall():
    """
    Return (x_an, M_an, p_an) as step-function arrays along the bottom wall.

    Regions:
      freestream (x < 0)         : M_inf
      front face (0 < x < 0.5)    : after leading oblique shock
      rear face  (0.5 < x < 1.0)  : after Prandtl-Meyer expansion
      downstream (x > 1.0)        : after trailing oblique shock
    """
    # Region 2: after leading-edge shock (deflection +theta)
    M2, p2_p1, beta1 = post_shock(M_INF, THETA)
    p2 = P_INF * p2_p1

    # Region 3: after expansion fan (total turn = 2*theta)
    M3, p3_p2 = post_expansion(M2, 2 * THETA)
    p3 = p2 * p3_p2

    # Region 4: after trailing-edge shock (deflection +theta relative to M3)
    M4, p4_p3, beta2 = post_shock(M3, THETA)
    p4 = p3 * p4_p3

    # Build step-function arrays
    x_an = np.array([-0.5, -0.001, 0.0, 0.499, 0.5, 0.999, 1.0, 1.5])
    M_an = np.array([M_INF, M_INF, M2,   M2,    M3,  M3,    M4,  M4])
    p_an = np.array([P_INF, P_INF, p2,   p2,    p3,  p3,    p4,  p4])

    return x_an, M_an, p_an


def print_analytical_summary():
    """Print a summary table of analytical values."""
    M2, p2_p1, beta1 = post_shock(M_INF, THETA)
    M3, p3_p2 = post_expansion(M2, 2 * THETA)
    M4, p4_p3, beta2 = post_shock(M3, THETA)

    print("=" * 60)
    print("Analytical Solution Summary (M_inf = 2, theta = 10 deg)")
    print("=" * 60)
    print(f"Freestream:           M = {M_INF:.4f}, p/p_inf = 1.0000")
    print(f"After leading shock:  M = {M2:.4f}, p/p_inf = {p2_p1:.4f}, beta = {beta1:.2f} deg")
    print(f"After expansion fan:  M = {M3:.4f}, p/p_inf = {p2_p1*p3_p2:.4f}")
    print(f"After trailing shock: M = {M4:.4f}, p/p_inf = {p2_p1*p3_p2*p4_p3:.4f}, beta = {beta2:.2f} deg")
    print("=" * 60)


if __name__ == "__main__":
    print_analytical_summary()
