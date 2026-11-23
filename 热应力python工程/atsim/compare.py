# -*- coding: utf-8 -*-
"""Compare Python FVM field vs COMSOL exported CSV.

Expected COMSOL export (see 00_COMSOL建模手册.md):
    x_m, y_m, z_m, T_K
or  x, y, z, T   (meters, Kelvin)
Header row required. Extra columns are ignored.
"""

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


def _open_text(path):
    try:
        return open(path, "r", encoding="utf-8")
    except TypeError:
        return open(path, "r")
    except UnicodeDecodeError:
        return open(path, "r", encoding="gbk")


def load_comsol_csv(path):
    """Return dict with arrays x,y,z,T (Kelvin)."""
    with _open_text(path) as f:
        # skip COMSOL '%' comment lines
        lines = []
        for line in f:
            s = line.strip()
            if not s or s.startswith("%") or s.startswith("#"):
                continue
            lines.append(s)
    if not lines:
        raise ValueError("empty COMSOL csv: " + path)
    header = [h.strip().lower() for h in lines[0].replace(";", ",").split(",")]
    # detect if first data row is numeric (no header)
    try:
        [float(v) for v in lines[0].replace(";", ",").split(",")[:4]]
        has_header = False
    except ValueError:
        has_header = True

    rows = lines[1:] if has_header else lines
    xs, ys, zs, Ts = [], [], [], []
    for row in rows:
        parts = [p.strip() for p in row.replace(";", ",").split(",")]
        if len(parts) < 4:
            continue
        xs.append(float(parts[0]))
        ys.append(float(parts[1]))
        zs.append(float(parts[2]))
        Ts.append(float(parts[3]))
    x = np.array(xs)
    y = np.array(ys)
    z = np.array(zs)
    T = np.array(Ts)
    # auto unit: if values look like mm, convert
    if np.nanmax(np.abs(x)) > 0.1:  # probably mm
        x *= 1e-3
        y *= 1e-3
        z *= 1e-3
    if np.nanmax(T) < 200:  # probably Celsius
        T = T + 273.15
    return {"x": x, "y": y, "z": z, "T": T, "path": path}


def sample_python_at_xy(py_result, x_m, y_m, z_m=None):
    """Bilinear sample of the Python mid-plane (or nearest z-slice) at (x,y)."""
    T = py_result["T"]
    grid = py_result["grid"]
    xc, yc, zc = grid["xc"], grid["yc"], grid["zc"]
    if z_m is None:
        k = T.shape[2] // 2
    else:
        k = int(np.argmin(np.abs(zc - np.asarray(z_m).reshape(-1)[0])))
    field = T[:, :, k]
    x_m = np.clip(np.asarray(x_m, dtype=float), xc[0], xc[-1])
    y_m = np.clip(np.asarray(y_m, dtype=float), yc[0], yc[-1])
    # indices of left/bottom cell center
    i0 = np.clip(np.searchsorted(xc, x_m) - 1, 0, len(xc) - 2)
    j0 = np.clip(np.searchsorted(yc, y_m) - 1, 0, len(yc) - 2)
    x0, x1 = xc[i0], xc[i0 + 1]
    y0, y1 = yc[j0], yc[j0 + 1]
    tx = np.where(x1 > x0, (x_m - x0) / (x1 - x0), 0.0)
    ty = np.where(y1 > y0, (y_m - y0) / (y1 - y0), 0.0)
    v00 = field[i0, j0]
    v10 = field[i0 + 1, j0]
    v01 = field[i0, j0 + 1]
    v11 = field[i0 + 1, j0 + 1]
    return (1 - tx) * (1 - ty) * v00 + tx * (1 - ty) * v10 + (1 - tx) * ty * v01 + tx * ty * v11


def compare_midplanes(py_result, comsol, out_dir, z_target=None):
    """Sample Python FVM onto COMSOL export nodes and report errors.

    COMSOL Case 1 export is a coarse mid-plane cloud (~21x21). Interpolating
    that cloud onto the 64x64 Python grid inflates error; do the reverse.
    """
    Tpy = py_result["T"]
    grid = py_result["grid"]
    case = py_result["case"]
    if z_target is None:
        z_target = float(np.median(comsol["z"]))

    Tpy_at = sample_python_at_xy(py_result, comsol["x"], comsol["y"], z_target)
    Tco = comsol["T"]
    err = Tpy_at - Tco
    abs_err = np.abs(err)
    rise = np.maximum(Tco - case.T_amb, 1.0)
    rel = abs_err / rise * 100.0

    # probes: bilinear Python vs nearest COMSOL node
    probe_rows = []
    for p in case.probes:
        Tpy_p = float(sample_python_at_xy(py_result, p.x, p.y, p.z))
        dist = (comsol["x"] - p.x) ** 2 + (comsol["y"] - p.y) ** 2 + (
            (comsol["z"] - p.z) * 50.0
        ) ** 2
        i = int(np.argmin(dist))
        Tco_p = float(Tco[i])
        probe_rows.append(
            {
                "name": p.name,
                "x_mm": p.x * 1e3,
                "y_mm": p.y * 1e3,
                "python_C": Tpy_p - 273.15,
                "comsol_C": Tco_p - 273.15,
                "dT_K": Tpy_p - Tco_p,
                "comsol_node_mm": (
                    comsol["x"][i] * 1e3,
                    comsol["y"][i] * 1e3,
                    comsol["z"][i] * 1e3,
                ),
            }
        )

    stats = {
        "z_m": z_target,
        "n_comsol_nodes": int(Tco.size),
        "mae_K": float(abs_err.mean()),
        "rmse_K": float(np.sqrt((err ** 2).mean())),
        "max_abs_K": float(abs_err.max()),
        "mean_rel_pct": float(rel.mean()),
        "max_rel_pct": float(rel.max()),
        "Tmax_python_C": float(Tpy.max() - 273.15),
        "Tmax_comsol_C": float(Tco.max() - 273.15),
        "Tmin_python_C": float(Tpy.min() - 273.15),
        "Tmin_comsol_C": float(Tco.min() - 273.15),
        "dTmax_K": float((Tpy.max() - Tco.max())),
        "probes": probe_rows,
    }

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    with open(os.path.join(out_dir, "compare_stats.txt"), "w") as f:
        for kname, v in stats.items():
            if kname == "probes":
                continue
            f.write("{0}: {1}\n".format(kname, v))
        f.write("\nprobes:\n")
        f.write("name,python_C,comsol_C,dT_K\n")
        for r in probe_rows:
            f.write("{0},{1:.4f},{2:.4f},{3:.4f}\n".format(
                r["name"], r["python_C"], r["comsol_C"], r["dT_K"]))

    with open(os.path.join(out_dir, "node_errors.csv"), "w") as f:
        f.write("x_mm,y_mm,z_mm,python_C,comsol_C,dT_K,rel_pct\n")
        for i in range(Tco.size):
            f.write("{0:.6f},{1:.6f},{2:.6f},{3:.6f},{4:.6f},{5:.6f},{6:.4f}\n".format(
                comsol["x"][i] * 1e3,
                comsol["y"][i] * 1e3,
                comsol["z"][i] * 1e3,
                Tpy_at[i] - 273.15,
                Tco[i] - 273.15,
                err[i],
                rel[i],
            ))

    if HAS_MPL:
        # scatter reconstruction of both fields on the COMSOL xy cloud
        fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
        vmin = min(Tpy_at.min(), Tco.min()) - 273.15
        vmax = max(Tpy_at.max(), Tco.max()) - 273.15
        a = max(abs(err.min()), abs(err.max()), 0.05)
        sets = [
            (axes[0], Tpy_at - 273.15, "Python FVM @ COMSOL nodes (C)", "inferno", vmin, vmax),
            (axes[1], Tco - 273.15, "COMSOL (C)", "inferno", vmin, vmax),
            (axes[2], err, "Python - COMSOL (K)", "coolwarm", -a, a),
        ]
        for ax, field, title, cmap, vmin_i, vmax_i in sets:
            sc = ax.scatter(
                comsol["x"] * 1e3,
                comsol["y"] * 1e3,
                c=field,
                cmap=cmap,
                s=28,
                vmin=vmin_i,
                vmax=vmax_i,
                edgecolors="none",
            )
            ax.set_aspect("equal")
            ax.set_xlim(0, case.Lx * 1e3)
            ax.set_ylim(0, case.Ly * 1e3)
            ax.set_title(title)
            ax.set_xlabel("x (mm)")
            ax.set_ylabel("y (mm)")
            fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "compare_midplane.png"), dpi=160)
        plt.close(fig)

        # also the dense Python mid-plane for visual reference
        k = Tpy.shape[2] // 2
        fig, ax = plt.subplots(figsize=(6.0, 5.2))
        extent = [0, case.Lx * 1e3, 0, case.Ly * 1e3]
        im = ax.imshow(
            (Tpy[:, :, k] - 273.15).T,
            origin="lower",
            extent=extent,
            cmap="inferno",
            aspect="equal",
            vmin=vmin,
            vmax=vmax,
        )
        ax.scatter(comsol["x"] * 1e3, comsol["y"] * 1e3, c="white", s=6, alpha=0.35)
        ax.set_title("Python mid-plane + COMSOL nodes")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "python_with_comsol_nodes.png"), dpi=160)
        plt.close(fig)

    return stats
