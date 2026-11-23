# -*- coding: utf-8 -*-
"""Heatmaps, probe tables, mid-plane CSV dumps."""

from __future__ import print_function

import csv
import os

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False


def _ensure_dir(path):
    if path and not os.path.isdir(path):
        os.makedirs(path)


def sample_probes(result):
    T = result["T"]
    grid = result["grid"]
    case = result["case"]
    xc, yc, zc = grid["xc"], grid["yc"], grid["zc"]
    rows = []
    for p in case.probes:
        i = int(np.clip(np.searchsorted(xc, p.x), 0, len(xc) - 1))
        j = int(np.clip(np.searchsorted(yc, p.y), 0, len(yc) - 1))
        k = int(np.clip(np.searchsorted(zc, p.z), 0, len(zc) - 1))
        # nearest cell center
        i = int(np.argmin(np.abs(xc - p.x)))
        j = int(np.argmin(np.abs(yc - p.y)))
        k = int(np.argmin(np.abs(zc - p.z)))
        Tk = float(T[i, j, k])
        rows.append(
            {
                "name": p.name,
                "x_m": p.x,
                "y_m": p.y,
                "z_m": p.z,
                "T_K": Tk,
                "T_C": Tk - 273.15,
                "i": i,
                "j": j,
                "k": k,
            }
        )
    return rows


def save_probes_csv(result, path):
    rows = sample_probes(result)
    _ensure_dir(os.path.dirname(path))
    with open(path, "w") as f:
        w = csv.DictWriter(
            f, fieldnames=["name", "x_m", "y_m", "z_m", "T_K", "T_C", "i", "j", "k"]
        )
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return rows


def die_mid_k(result):
    """z-index of the die mid-plane (package-aware)."""
    case = result["case"]
    zc = result["grid"]["zc"]
    if getattr(case, "has_package", False):
        from .package_grid import die_layer
        die = die_layer(case)
        zmid = 0.5 * (die.z0 + die.z1)
        return int(np.argmin(np.abs(zc - zmid)))
    return result["T"].shape[2] // 2


def save_midplane_csv(result, path, k=None):
    """Write x,y,T_C of a z-slice (default: die mid-plane)."""
    T = result["T"]
    grid = result["grid"]
    if k is None:
        k = die_mid_k(result)
    _ensure_dir(os.path.dirname(path))
    with open(path, "w") as f:
        f.write("x_m,y_m,z_m,T_K,T_C\n")
        z = grid["zc"][k]
        solid = grid.get("solid")
        for i, x in enumerate(grid["xc"]):
            for j, y in enumerate(grid["yc"]):
                if solid is not None and not solid[i, j, k]:
                    continue
                Tk = float(T[i, j, k])
                f.write("{0:.8e},{1:.8e},{2:.8e},{3:.6f},{4:.6f}\n".format(x, y, z, Tk, Tk - 273.15))
    return k


def save_numpy(result, path):
    _ensure_dir(os.path.dirname(path))
    np.savez_compressed(
        path,
        T=result["T"],
        xc=result["grid"]["xc"],
        yc=result["grid"]["yc"],
        zc=result["grid"]["zc"],
        Pvol=result["Pvol"],
    )


def plot_midplane(result, path, title=None):
    if not HAS_MPL:
        print("[warn] matplotlib unavailable, skip plot", path)
        return
    T = result["T"]
    grid = result["grid"]
    case = result["case"]
    k = die_mid_k(result)
    Tc = T[:, :, k] - 273.15
    solid = grid.get("solid")
    if solid is not None:
        Tc = np.ma.array(Tc, mask=~solid[:, :, k])
    _ensure_dir(os.path.dirname(path))
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    xe = grid.get("x_edges")
    ye = grid.get("y_edges")
    if xe is not None and ye is not None:
        XE, YE = np.meshgrid(xe * 1e3, ye * 1e3, indexing="ij")
        im = ax.pcolormesh(XE, YE, Tc, cmap="inferno", shading="flat")
        ax.set_aspect("equal")
        ax.set_xlim(0.0, case.Lx * 1e3)
        ax.set_ylim(0.0, case.Ly * 1e3)
    else:
        extent = [0.0, case.Lx * 1e3, 0.0, case.Ly * 1e3]
        im = ax.imshow(
            Tc.T,
            origin="lower",
            extent=extent,
            cmap="inferno",
            aspect="equal",
        )
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Temperature (C)")
    ox = getattr(case, "chip_ox", 0.0) * 1e3
    oy = getattr(case, "chip_oy", 0.0) * 1e3
    for b in case.blocks:
        rect = plt.Rectangle(
            (ox + b.x0 * 1e3, oy + b.y0 * 1e3),
            b.w * 1e3,
            b.h * 1e3,
            fill=False,
            edgecolor="white",
            linewidth=0.8,
            linestyle="--",
        )
        ax.add_patch(rect)
        ax.text(
            ox + (b.x0 + 0.5 * b.w) * 1e3,
            oy + (b.y0 + 0.5 * b.h) * 1e3,
            "{0}\n{1:.1f} W".format(b.name, b.power),
            color="white",
            ha="center",
            va="center",
            fontsize=8,
        )
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    if title is None:
        title = "%s  mid-plane z=%.0f um  Tmax=%.2f C" % (
            case.name,
            grid["zc"][k] * 1e6,
            T.max() - 273.15,
        )
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
