# -*- coding: utf-8 -*-
"""Package-level warpage w(r) = uz on the interposer mid-plane.

ATMPlace compact warpage is a package-level metric: the out-of-plane field
on the silicon interposer, not die-interior von Mises. This module does not
solve a new biharmonic plate; it samples the existing 3D thermoelastic uz
on that plane (same Cartesian voxels as FVM heat / hex elasticity).

Peak-to-valley is the layout quantity. substrate_bottom_fixed already holds
rigid motion, so PV is bending, not a free-body offset.
"""

from __future__ import print_function

import csv
import os

import numpy as np


def interposer_layer(case):
    for L in getattr(case, "stack", []) or []:
        if getattr(L, "name", None) == "interposer":
            return L
    return None


def interposer_k(elast):
    case = elast["case"]
    zc = elast["grid"]["zc"]
    layer = interposer_layer(case)
    if layer is not None:
        zmid = 0.5 * (layer.z0 + layer.z1)
        return int(np.argmin(np.abs(zc - zmid))), layer
    from .visualize import die_mid_k
    return die_mid_k({"grid": elast["grid"], "case": case, "T": elast["thermal"]["T"]}), None


def package_warpage(elast):
    """ATMPlace-style package warpage from 3D uz on the interposer mid-plane.

    w(r) = uz. Metrics: peak-to-valley and corner-to-corner on that plane.
    """
    grid = elast["grid"]
    case = elast["case"]
    k, layer = interposer_k(elast)
    uz = elast["ucz"][:, :, k]
    solid = grid.get("solid")
    if solid is not None:
        m = solid[:, :, k]
    else:
        m = np.ones(uz.shape, dtype=bool)
    if layer is not None:
        xc, yc = grid["xc"], grid["yc"]
        m = m & (xc[:, None] >= layer.x0 - 1e-16) & (xc[:, None] <= layer.x0 + layer.w + 1e-16)
        m = m & (yc[None, :] >= layer.y0 - 1e-16) & (yc[None, :] <= layer.y0 + layer.h + 1e-16)
    vals = uz[m]
    if vals.size == 0:
        raise ValueError("no solid cells on interposer mid-plane for warpage")
    w_um = vals * 1e6
    ii, jj = np.where(m)
    i_hi = int(np.argmax(vals))
    i_lo = int(np.argmin(vals))
    corners = []
    if layer is not None:
        pts = [
            ("SW", layer.x0, layer.y0),
            ("SE", layer.x0 + layer.w, layer.y0),
            ("NW", layer.x0, layer.y0 + layer.h),
            ("NE", layer.x0 + layer.w, layer.y0 + layer.h),
        ]
        for name, x, y in pts:
            i = int(np.argmin(np.abs(grid["xc"] - x)))
            j = int(np.argmin(np.abs(grid["yc"] - y)))
            if not m[i, j]:
                d2 = (ii - i) ** 2 + (jj - j) ** 2
                t = int(np.argmin(d2))
                i, j = int(ii[t]), int(jj[t])
            corners.append({
                "name": name,
                "x_mm": float(grid["xc"][i] * 1e3),
                "y_mm": float(grid["yc"][j] * 1e3),
                "uz_um": float(uz[i, j] * 1e6),
            })
    corner_uz = [c["uz_um"] for c in corners] if corners else []
    T = elast["thermal"]["T"][:, :, k]
    return {
        "plane": "interposer_mid" if layer is not None else "die_mid",
        "z_mm": float(grid["zc"][k] * 1e3),
        "n_cells": int(m.sum()),
        "w_max_um": float(w_um.max()),
        "w_min_um": float(w_um.min()),
        "w_mean_um": float(w_um.mean()),
        "w_pv_um": float(w_um.max() - w_um.min()),
        "w_abs_max_um": float(np.abs(w_um).max()),
        "xy_max_mm": [float(grid["xc"][ii[i_hi]] * 1e3), float(grid["yc"][jj[i_hi]] * 1e3)],
        "xy_min_mm": [float(grid["xc"][ii[i_lo]] * 1e3), float(grid["yc"][jj[i_lo]] * 1e3)],
        "corners": corners,
        "corner_pv_um": (float(max(corner_uz) - min(corner_uz)) if corner_uz else None),
        "T_mean_C": float(T[m].mean() - 273.15),
        "T_max_C": float(T[m].max() - 273.15),
        "note": (
            "Package-level w(r)=uz on interposer mid-plane. "
            "Peak-to-valley is the layout metric; not die-interior von Mises."
        ),
        "k": k,
        "mask": m,
    }


def warpage_for_json(warp):
    return {key: warp[key] for key in warp if key not in ("mask", "k")}


def save_warpage_csv(elast, warp, path):
    grid = elast["grid"]
    k = warp["k"]
    m = warp["mask"]
    uz = elast["ucz"][:, :, k]
    T = elast["thermal"]["T"][:, :, k]
    ii, jj = np.where(m)
    parent = os.path.dirname(path)
    if parent and not os.path.isdir(parent):
        os.makedirs(parent)
    with open(path, "w") as f:
        w = csv.writer(f)
        w.writerow(["i", "j", "x_mm", "y_mm", "w_um", "T_C"])
        for t in range(ii.size):
            i, j = int(ii[t]), int(jj[t])
            w.writerow([
                i, j,
                grid["xc"][i] * 1e3,
                grid["yc"][j] * 1e3,
                uz[i, j] * 1e6,
                T[i, j] - 273.15,
            ])
