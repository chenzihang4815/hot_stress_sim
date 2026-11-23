# -*- coding: utf-8 -*-
"""Draw Case 1 floorplan for the COMSOL modeling sheet."""

from __future__ import print_function

import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from atsim.geometry import load_case


def main():
    case = load_case(os.path.join(ROOT, "cases", "case1_2d_ic.json"))
    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    colors = {
        "Cache": "#4C78A8",
        "Uncore": "#72B7B2",
        "CoreA": "#E45756",
        "CoreB": "#F58518",
    }
    for b in case.blocks:
        ax.add_patch(
            Rectangle(
                (b.x0 * 1e3, b.y0 * 1e3),
                b.w * 1e3,
                b.h * 1e3,
                facecolor=colors.get(b.name, "#CCCCCC"),
                edgecolor="black",
                linewidth=1.4,
            )
        )
        dens = b.power / (b.area * 1e6)
        ax.text(
            (b.x0 + 0.5 * b.w) * 1e3,
            (b.y0 + 0.5 * b.h) * 1e3,
            "{0}\n{1:.1f} W\n{2:.2f} W/mm$^2$".format(b.name, b.power, dens),
            ha="center",
            va="center",
            fontsize=11,
            color="white",
            fontweight="bold",
        )
    ax.set_xlim(0, case.Lx * 1e3)
    ax.set_ylim(0, case.Ly * 1e3)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)  origin = bottom-left")
    ax.set_ylabel("y (mm)")
    ax.set_title(
        "Case 1 floorplan   8 x 8 x 0.15 mm silicon   total %.1f W\n"
        "top convection h=5000 W/m$^2$K, T$_\\infty$=45 C; other faces adiabatic"
        % case.total_power
    )
    ax.set_xticks([0, 2, 4, 6, 8])
    ax.set_yticks([0, 2, 4, 6, 8])
    ax.grid(True, linestyle=":", alpha=0.4)
    fig.tight_layout()
    out = os.path.join(ROOT, "results", "python")
    if not os.path.isdir(out):
        os.makedirs(out)
    path = os.path.join(out, "floorplan.png")
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print("wrote", path)


if __name__ == "__main__":
    main()
