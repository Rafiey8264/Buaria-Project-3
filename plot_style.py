"""
plot_style.py - Common matplotlib styling for all plots.
Import this at the top of any plot script to get consistent LaTeX-style fonts.
"""
import matplotlib.pyplot as plt
import matplotlib as mpl

def set_latex_style():
    """Apply LaTeX-style fonts to matplotlib globally."""
    # Use mathtext for LaTeX-like rendering (no external LaTeX install needed)
    mpl.rcParams.update({
        "text.usetex": False,              # False = use mathtext (faster, no LaTeX install needed)
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",          # Computer Modern math fonts
        "axes.labelsize": 13,
        "axes.titlesize": 14,
        "legend.fontsize": 11,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "figure.titlesize": 15,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "lines.linewidth": 1.8,
        "figure.dpi": 100,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
    })

# Apply on import
set_latex_style()

# Consistent colors for all three schemes
COLORS = {
    "1st_order": "tab:blue",
    "2nd_order": "tab:orange",
    "tvd":       "tab:green",
    "analytical": "black",
    "coarse":    "tab:red",
    "fine":      "tab:blue",
}

# Consistent labels
LABELS = {
    "1st_order": "1st Order ($\\phi = 0$)",
    "2nd_order": "2nd Order ($\\phi = 1$)",
    "tvd":       "TVD Limiter",
    "analytical": "Analytical",
}
