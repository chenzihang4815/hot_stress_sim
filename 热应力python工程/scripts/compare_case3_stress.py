# -*- coding: utf-8 -*-
"""Compare Case 3 Python heat+stress vs a COMSOL (x,y,z,T,vm) CSV.

Default CSV:
    results/python_case3_chips/case3_die_points_stress.csv

Usage:
    python scripts/compare_case3_stress.py
    python scripts/compare_case3_stress.py path/to/comsol.csv
    python scripts/compare_case3_stress.py path.csv --case cases/case3_chips.json
"""

from __future__ import print_function

import csv
import json
import os
import sys

import numpy as np
from scipy.interpolate import griddata

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from atsim.geometry import load_case
from atsim.fvm3d import solve_fvm3d
from atsim.elasticity import solve_elasticity
from atsim.visualize import die_mid_k

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False


def _ensure(p):
    if p and not os.path.isdir(p):
        os.makedirs(p)


def _open_text(path):
    try:
        return open(path, "r", encoding="utf-8")
    except TypeError:
        return open(path, "r")
    except UnicodeDecodeError:
        return open(path, "r", encoding="gbk")


def _mpa(a):
    a = np.asarray(a, dtype=float)
    if a.size and np.nanmax(np.abs(a)) > 1e4:
        a = a * 1e-6
    return a


def load_stress_csv(path):
    xs, ys, zs, Ts, vms = [], [], [], [], []
    sxs, sys_, szs = [], [], []
    with _open_text(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("%") or s.startswith("#"):
                continue
            parts = [p.strip() for p in s.replace(";", ",").split(",")]
            try:
                float(parts[0])
            except ValueError:
                continue
            if len(parts) < 5:
                continue
            xs.append(float(parts[0]))
            ys.append(float(parts[1]))
            zs.append(float(parts[2]))
            Ts.append(float(parts[3]))
            vms.append(float(parts[4]))
            if len(parts) >= 8:
                sxs.append(float(parts[5]))
                sys_.append(float(parts[6]))
                szs.append(float(parts[7]))
    x, y, z, T, vm = map(np.array, (xs, ys, zs, Ts, vms))
    if np.nanmax(np.abs(x)) > 0.1:
        x = x * 1e-3
        y = y * 1e-3
        z = z * 1e-3
    if np.nanmax(T) < 200:
        T = T + 273.15
    out = {
        "x": x, "y": y, "z": z, "T": T, "vm_MPa": _mpa(vm), "path": path,
    }
    if sxs:
        out["sx_MPa"] = _mpa(sxs)
        out["sy_MPa"] = _mpa(sys_)
        out["sz_MPa"] = _mpa(szs)
    return out


def sample_xy(field, grid, x, y, z):
    xc, yc, zc = grid["xc"], grid["yc"], grid["zc"]
    i = np.clip(np.searchsorted(xc, x) - 1, 0, len(xc) - 2)
    j = np.clip(np.searchsorted(yc, y) - 1, 0, len(yc) - 2)
    x0, x1 = xc[i], xc[i + 1]
    y0, y1 = yc[j], yc[j + 1]
    tx = np.where(x1 > x0, (x - x0) / (x1 - x0), 0.0)
    ty = np.where(y1 > y0, (y - y0) / (y1 - y0), 0.0)
    kk = int(np.argmin(np.abs(zc - np.median(z))))
    v00 = field[i, j, kk]
    v10 = field[i + 1, j, kk]
    v01 = field[i, j + 1, kk]
    v11 = field[i + 1, j + 1, kk]
    return (1 - tx) * (1 - ty) * v00 + tx * (1 - ty) * v10 + (1 - tx) * ty * v01 + tx * ty * v11


def chip_mask_m(case, x, y):
    ox = getattr(case, "chip_ox", 0.0)
    oy = getattr(case, "chip_oy", 0.0)
    m = np.zeros(x.shape, dtype=bool)
    for b in case.blocks:
        m |= (
            (x >= ox + b.x0 - 1e-9)
            & (x <= ox + b.x0 + b.w + 1e-9)
            & (y >= oy + b.y0 - 1e-9)
            & (y <= oy + b.y0 + b.h + 1e-9)
        )
    return m


def chip_window_mm(case, pad_mm=3.0):
    ox = getattr(case, "chip_ox", 0.0) * 1e3
    oy = getattr(case, "chip_oy", 0.0) * 1e3
    x0 = min(ox + b.x0 * 1e3 for b in case.blocks) - pad_mm
    x1 = max(ox + (b.x0 + b.w) * 1e3 for b in case.blocks) + pad_mm
    y0 = min(oy + b.y0 * 1e3 for b in case.blocks) - pad_mm
    y1 = max(oy + (b.y0 + b.h) * 1e3 for b in case.blocks) + pad_mm
    return x0, x1, y0, y1


def draw_chips(ax, case, edgecolor="white"):
    ox = getattr(case, "chip_ox", 0.0) * 1e3
    oy = getattr(case, "chip_oy", 0.0) * 1e3
    for b in case.blocks:
        ax.add_patch(plt.Rectangle(
            (ox + b.x0 * 1e3, oy + b.y0 * 1e3), b.w * 1e3, b.h * 1e3,
            fill=False, edgecolor=edgecolor, linewidth=0.9, linestyle="--"))
        ax.text(
            ox + (b.x0 + 0.5 * b.w) * 1e3,
            oy + (b.y0 + 0.5 * b.h) * 1e3,
            b.name, color=edgecolor, ha="center", va="center", fontsize=8,
        )


def err_stats(a, b):
    d = a - b
    return {
        "n": int(a.size),
        "mae": float(np.abs(d).mean()),
        "rmse": float(np.sqrt((d ** 2).mean())),
        "max_abs": float(np.abs(d).max()),
        "mean_py": float(a.mean()),
        "mean_co": float(b.mean()),
    }


def plot_triplet(X, Y, A, B, dA, path, titles, unit, cmap, vlim, elims, case):
    if not HAS_MPL:
        return
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.2))
    fields = (A, B, dA)
    cmaps = (cmap, cmap, "RdBu_r")
    lims = (vlim, vlim, elims)
    x0, x1, y0, y1 = chip_window_mm(case)
    for ax, F, cm, lim, title in zip(axes, fields, cmaps, lims, titles):
        im = ax.pcolormesh(
            X, Y, F, cmap=cm, shading="auto",
            vmin=lim[0], vmax=lim[1],
        )
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(unit)
        draw_chips(ax, case, edgecolor="white" if cm != "RdBu_r" else "k")
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_scatter(py, co, path, xlabel, ylabel, title):
    if not HAS_MPL:
        return
    lo = min(float(py.min()), float(co.min()))
    hi = max(float(py.max()), float(co.max()))
    pad = 0.04 * (hi - lo + 1e-9)
    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    ax.scatter(co, py, s=6, alpha=0.35, c="#2c3e50", linewidths=0)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "r--", linewidth=1.0)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_probe_compare(rows, path):
    if not HAS_MPL or not rows:
        return
    names = [r["name"].replace("_center", "") for r in rows]
    x = np.arange(len(names))
    w = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))
    axes[0].bar(x - w / 2, [r["T_py"] for r in rows], w, label="Python", color="#d35400")
    axes[0].bar(x + w / 2, [r["T_co"] for r in rows], w, label="COMSOL", color="#f5cba7")
    axes[0].set_title("T (C)")
    axes[1].bar(x - w / 2, [r["vm_py"] for r in rows], w, label="Python", color="#8e44ad")
    axes[1].bar(x + w / 2, [r["vm_co"] for r in rows], w, label="COMSOL", color="#d2b4de")
    axes[1].set_title("von Mises (MPa)")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=20, ha="right")
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    csv_path = os.path.join(
        ROOT, "results", "python_case3_chips_stress", "case3_die_points_stress.csv"
    )
    case_path = os.path.join(ROOT, "cases", "case3_chips.json")
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--case" and i + 1 < len(args):
            case_path = args[i + 1]
            if not os.path.isabs(case_path):
                case_path = os.path.join(ROOT, case_path)
            i += 2
            continue
        csv_path = args[i]
        if not os.path.isabs(csv_path):
            csv_path = os.path.join(ROOT, csv_path)
        i += 1
    if not os.path.isfile(csv_path):
        print("[error] missing COMSOL csv:", csv_path)
        sys.exit(1)

    name = os.path.splitext(os.path.basename(case_path))[0]
    out = os.path.join(ROOT, "results", "compare_" + name + "_stress")
    _ensure(out)

    co = load_stress_csv(csv_path)
    print("COMSOL n=%d  T=%.2f..%.2f C  vm=%.1f..%.1f MPa" % (
        co["T"].size, co["T"].min() - 273.15, co["T"].max() - 273.15,
        co["vm_MPa"].min(), co["vm_MPa"].max()))

    case = load_case(case_path)
    print("Case:", case.name)
    print("\n-- heat --")
    th = solve_fvm3d(case, verbose=True)
    print("\n-- elasticity --")
    el = solve_elasticity(th, T_ref=case.T_amb, verbose=True)
    grid = th["grid"]

    pyT = sample_xy(th["T"], grid, co["x"], co["y"], co["z"])
    pyVm = sample_xy(el["vm"], grid, co["x"], co["y"], co["z"]) * 1e-6
    pySx = sample_xy(el["sxx"], grid, co["x"], co["y"], co["z"]) * 1e-6
    pySy = sample_xy(el["syy"], grid, co["x"], co["y"], co["z"]) * 1e-6
    pySz = sample_xy(el["szz"], grid, co["x"], co["y"], co["z"]) * 1e-6
    has_tensor = "sx_MPa" in co
    inside = chip_mask_m(case, co["x"], co["y"])
    solid = grid.get("solid")
    kmid = die_mid_k(th)

    def valid(mask):
        ok = mask.copy()
        if solid is not None:
            i = np.clip(np.searchsorted(grid["xc"], co["x"]) - 1, 0, grid["nx"] - 1)
            j = np.clip(np.searchsorted(grid["yc"], co["y"]) - 1, 0, grid["ny"] - 1)
            ok &= solid[i, j, kmid]
        return ok

    all_ok = valid(np.ones(co["T"].size, dtype=bool))
    chip_ok = valid(inside)
    st_T_all = err_stats(pyT[all_ok] - 273.15, co["T"][all_ok] - 273.15)
    st_T = err_stats(pyT[chip_ok] - 273.15, co["T"][chip_ok] - 273.15)
    st_vm_all = err_stats(pyVm[all_ok], co["vm_MPa"][all_ok])
    st_vm = err_stats(pyVm[chip_ok], co["vm_MPa"][chip_ok])

    print("\nOn COMSOL nodes, chips only (n=%d):" % st_T["n"])
    print("  T    MAE=%.3f K  RMSE=%.3f K  max|d|=%.3f K" % (
        st_T["mae"], st_T["rmse"], st_T["max_abs"]))
    print("  vm   MAE=%.2f MPa  RMSE=%.2f MPa  max|d|=%.2f MPa" % (
        st_vm["mae"], st_vm["rmse"], st_vm["max_abs"]))
    st_sx = st_sy = st_sz = None
    if has_tensor:
        st_sx = err_stats(pySx[chip_ok], co["sx_MPa"][chip_ok])
        st_sy = err_stats(pySy[chip_ok], co["sy_MPa"][chip_ok])
        st_sz = err_stats(pySz[chip_ok], co["sz_MPa"][chip_ok])
        print("  sx   MAE=%.2f MPa  RMSE=%.2f MPa  max|d|=%.2f MPa" % (
            st_sx["mae"], st_sx["rmse"], st_sx["max_abs"]))
        print("  sy   MAE=%.2f MPa  RMSE=%.2f MPa  max|d|=%.2f MPa" % (
            st_sy["mae"], st_sy["rmse"], st_sy["max_abs"]))
        print("  sz   MAE=%.2f MPa  RMSE=%.2f MPa  max|d|=%.2f MPa" % (
            st_sz["mae"], st_sz["rmse"], st_sz["max_abs"]))
    print("All exported nodes (n=%d): T MAE=%.3f K   vm MAE=%.2f MPa" % (
        st_T_all["n"], st_T_all["mae"], st_vm_all["mae"]))

    probes = []
    for p in case.probes:
        if not p.name.endswith("_center"):
            continue
        if p.name in ("Sink_top_center", "Spreader_center", "TIM_mid",
                      "Interposer_center", "Substrate_center"):
            continue
        dist = (co["x"] - p.x) ** 2 + (co["y"] - p.y) ** 2
        ic = int(dist.argmin())
        Tpy = float(sample_xy(th["T"], grid, np.array([p.x]), np.array([p.y]), co["z"])[0])
        vpy = float(sample_xy(el["vm"], grid, np.array([p.x]), np.array([p.y]), co["z"])[0] * 1e-6)
        row = {
            "name": p.name,
            "x_mm": p.x * 1e3,
            "y_mm": p.y * 1e3,
            "T_py": Tpy - 273.15,
            "T_co": float(co["T"][ic] - 273.15),
            "dT_K": Tpy - float(co["T"][ic]),
            "vm_py": vpy,
            "vm_co": float(co["vm_MPa"][ic]),
            "dvm_MPa": vpy - float(co["vm_MPa"][ic]),
            "comsol_xy_mm": (float(co["x"][ic] * 1e3), float(co["y"][ic] * 1e3)),
        }
        if has_tensor:
            sxpy = float(sample_xy(el["sxx"], grid, np.array([p.x]), np.array([p.y]), co["z"])[0] * 1e-6)
            sypy = float(sample_xy(el["syy"], grid, np.array([p.x]), np.array([p.y]), co["z"])[0] * 1e-6)
            szpy = float(sample_xy(el["szz"], grid, np.array([p.x]), np.array([p.y]), co["z"])[0] * 1e-6)
            row.update({
                "sx_py": sxpy, "sx_co": float(co["sx_MPa"][ic]),
                "dsx_MPa": sxpy - float(co["sx_MPa"][ic]),
                "sy_py": sypy, "sy_co": float(co["sy_MPa"][ic]),
                "dsy_MPa": sypy - float(co["sy_MPa"][ic]),
                "sz_py": szpy, "sz_co": float(co["sz_MPa"][ic]),
                "dsz_MPa": szpy - float(co["sz_MPa"][ic]),
            })
        probes.append(row)
    print("\n%-14s %8s %8s %8s %8s %8s %8s" % (
        "pt", "T_py", "T_co", "dT", "vm_py", "vm_co", "dvm"))
    for r in probes:
        print(" %-14s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" % (
            r["name"], r["T_py"], r["T_co"], r["dT_K"],
            r["vm_py"], r["vm_co"], r["dvm_MPa"]))
    if has_tensor:
        print("\n%-14s %8s %8s %8s %8s %8s %8s %8s %8s %8s" % (
            "pt", "sx_py", "sx_co", "dsx", "sy_py", "sy_co", "dsy", "sz_py", "sz_co", "dsz"))
        for r in probes:
            print(" %-14s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" % (
                r["name"], r["sx_py"], r["sx_co"], r["dsx_MPa"],
                r["sy_py"], r["sy_co"], r["dsy_MPa"],
                r["sz_py"], r["sz_co"], r["dsz_MPa"]))

    x0, x1, y0, y1 = chip_window_mm(case, pad_mm=2.0)
    nxg, nyg = 220, 200
    gx = np.linspace(x0, x1, nxg)
    gy = np.linspace(y0, y1, nyg)
    GX, GY = np.meshgrid(gx, gy, indexing="ij")
    pts = np.column_stack((co["x"][chip_ok] * 1e3, co["y"][chip_ok] * 1e3))
    Tco_g = griddata(pts, co["T"][chip_ok] - 273.15, (GX, GY), method="linear")
    Vco_g = griddata(pts, co["vm_MPa"][chip_ok], (GX, GY), method="linear")
    Tpy_g = sample_xy(
        th["T"], grid, GX.ravel() * 1e-3, GY.ravel() * 1e-3, co["z"]
    ).reshape(GX.shape) - 273.15
    Vpy_g = sample_xy(
        el["vm"], grid, GX.ravel() * 1e-3, GY.ravel() * 1e-3, co["z"]
    ).reshape(GX.shape) * 1e-6
    onchip = chip_mask_m(case, GX * 1e-3, GY * 1e-3)
    for arr in (Tco_g, Vco_g, Tpy_g, Vpy_g):
        arr[~onchip] = np.nan
        arr[np.isnan(Tco_g)] = np.nan

    dT = Tpy_g - Tco_g
    dV = Vpy_g - Vco_g
    tmin = np.nanmin([np.nanmin(Tpy_g), np.nanmin(Tco_g)])
    tmax = np.nanmax([np.nanmax(Tpy_g), np.nanmax(Tco_g)])
    vmax = np.nanpercentile(np.concatenate([
        Vpy_g[np.isfinite(Vpy_g)], Vco_g[np.isfinite(Vco_g)]
    ]), 98)
    eT = np.nanmax(np.abs(dT))
    eV = np.nanpercentile(np.abs(dV[np.isfinite(dV)]), 98)
    eV = max(float(eV), 1.0)

    plot_triplet(
        GX, GY, Tpy_g, Tco_g, dT,
        os.path.join(out, "compare_T.png"),
        ("Python T", "COMSOL T", "Python - COMSOL"),
        "C", "inferno", (tmin, tmax), (-eT, eT), case,
    )
    plot_triplet(
        GX, GY, Vpy_g, Vco_g, dV,
        os.path.join(out, "compare_mises.png"),
        ("Python vm", "COMSOL vm", "Python - COMSOL"),
        "MPa", "inferno", (0.0, float(vmax)), (-eV, eV), case,
    )
    plot_scatter(
        pyT[chip_ok] - 273.15, co["T"][chip_ok] - 273.15,
        os.path.join(out, "scatter_T.png"),
        "COMSOL T (C)", "Python T (C)",
        "T  MAE=%.3f K  RMSE=%.3f K" % (st_T["mae"], st_T["rmse"]),
    )
    plot_scatter(
        pyVm[chip_ok], co["vm_MPa"][chip_ok],
        os.path.join(out, "scatter_mises.png"),
        "COMSOL vm (MPa)", "Python vm (MPa)",
        "vm  MAE=%.2f MPa  RMSE=%.2f MPa" % (st_vm["mae"], st_vm["rmse"]),
    )
    if has_tensor:
        def _grid_comp(py_field, co_arr):
            cg = griddata(pts, co_arr[chip_ok], (GX, GY), method="linear")
            pg = sample_xy(
                py_field, grid, GX.ravel() * 1e-3, GY.ravel() * 1e-3, co["z"]
            ).reshape(GX.shape) * 1e-6
            pg[~onchip] = np.nan
            cg[~onchip] = np.nan
            pg[np.isnan(cg)] = np.nan
            return pg, cg

        for key, py_field, co_key, label in (
            ("sx", el["sxx"], "sx_MPa", "sx"),
            ("sy", el["syy"], "sy_MPa", "sy"),
            ("sz", el["szz"], "sz_MPa", "sz"),
        ):
            pg, cg = _grid_comp(py_field, co[co_key])
            d = pg - cg
            finite = np.concatenate([pg[np.isfinite(pg)], cg[np.isfinite(cg)]])
            lim = float(np.nanpercentile(np.abs(finite), 98))
            lim = max(lim, 1.0)
            ed = float(np.nanpercentile(np.abs(d[np.isfinite(d)]), 98))
            ed = max(ed, 1.0)
            st = {"sx": st_sx, "sy": st_sy, "sz": st_sz}[key]
            plot_triplet(
                GX, GY, pg, cg, d,
                os.path.join(out, "compare_%s.png" % key),
                ("Python %s" % label, "COMSOL %s" % label, "Python - COMSOL"),
                "MPa", "coolwarm", (-lim, lim), (-ed, ed), case,
            )
            py_arr = {"sx": pySx, "sy": pySy, "sz": pySz}[key]
            plot_scatter(
                py_arr[chip_ok], co[co_key][chip_ok],
                os.path.join(out, "scatter_%s.png" % key),
                "COMSOL %s (MPa)" % label, "Python %s (MPa)" % label,
                "%s  MAE=%.2f MPa  RMSE=%.2f MPa" % (label, st["mae"], st["rmse"]),
            )
    plot_probe_compare(probes, os.path.join(out, "compare_probes.png"))

    stats = {
        "case": case.name,
        "comsol_csv": csv_path,
        "T_chips": st_T,
        "vm_chips": st_vm,
        "sx_chips": st_sx,
        "sy_chips": st_sy,
        "sz_chips": st_sz,
        "T_all_nodes": st_T_all,
        "vm_all_nodes": st_vm_all,
        "probes": probes,
        "note": "MAE on COMSOL nodes inside chip xy. Edge von Mises is a singularity. sx/sy/sz are COMSOL solid.sx/sy/sz vs Python sxx/syy/szz.",
    }
    with open(os.path.join(out, "compare_stats.json"), "w") as f:
        json.dump(stats, f, indent=2)
    with open(os.path.join(out, "compare_probes.csv"), "w") as f:
        fields = [
            "name", "x_mm", "y_mm", "T_py", "T_co", "dT_K",
            "vm_py", "vm_co", "dvm_MPa",
        ]
        if has_tensor:
            fields.extend([
                "sx_py", "sx_co", "dsx_MPa",
                "sy_py", "sy_co", "dsy_MPa",
                "sz_py", "sz_co", "dsz_MPa",
            ])
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in probes:
            row = dict(r)
            row.pop("comsol_xy_mm", None)
            w.writerow(row)

    print("\nWrote", out)
    print("  compare_T.png  compare_mises.png  scatter_T.png  scatter_mises.png  compare_probes.png")
    if has_tensor:
        print("  compare_sx.png  compare_sy.png  compare_sz.png  scatter_sx/sy/sz.png")


if __name__ == "__main__":
    main()
