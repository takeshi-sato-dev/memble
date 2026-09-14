"""One style for every figure of the MEMBLE article.

The categorical colours were validated with the dataviz palette validator
(light surface, 4 slots): lightness band, chroma floor, CVD separation,
normal-vision floor and contrast against the surface all pass. Grey is not a
categorical slot; it is the recessive role for series that carry no identity.

Every series also carries a marker or a dash pattern, so the figures survive
greyscale printing without relying on hue.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

BLUE    = "#2b6cb0"   # cholesterol, and the fitted curve
ORANGE  = "#c05621"   # the routes of Table 4
GREEN   = "#2f855a"   # the large box
MAGENTA = "#97266d"   # reserved
GREY    = "#718096"   # species that lie on zero, reference lines
INK     = "#1a202c"
MUTED   = "#4a5568"
SURFACE = "#ffffff"

MARKERS = ["o", "s", "^", "D"]
DASHES  = [(None, None), (5, 2), (1.5, 1.5), (7, 2, 1.5, 2)]


def use():
    mpl.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.facecolor": SURFACE,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.6,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "lines.linewidth": 1.4,
        "lines.markersize": 4,
        "grid.color": "#e2e8f0",
        "grid.linewidth": 0.5,
    })


def panel_label(ax, text, dx=-0.16, dy=1.04):
    ax.text(dx, dy, text, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="left", color=INK)
