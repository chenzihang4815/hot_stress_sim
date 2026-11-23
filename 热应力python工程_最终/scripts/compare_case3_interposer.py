# -*- coding: utf-8 -*-
"""Optional Python vs COMSOL compare on the interposer mid-plane.

Not part of the Case 3 solve path. Run only when a COMSOL CSV exists:

    python scripts/compare_case3_interposer.py
    python scripts/compare_case3_interposer.py path/to/case3_interposer_uz.csv
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
from atsim.warpage import interposer_k, package_warpage

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


def _interior(x, y, layer, margin_m=1.0e-3):
    return (
        (x >= layer.x0 + margin_m) & (x <= layer.x0 + layer.w - margin_m)
        & (y >= layer.y0 + margin_m) & (y <= layer.y0 + layer.h - margin_m)
    )


def _rel_pct(mae, mean_abs_co):
    if mean_abs_co is None or mean_abs_co < 1e-12:
        return None
    return 100.0 * mae / mean_abs_co


def _draw_chips(ax, case, edgecolor="white"):
    ox = getattr(case, "chip_ox", 0.0) * 1e3
    oy = getattr(case, "chip_oy", 0.0) * 1e3
    for b in case.blocks:
        ax.add_patch(plt.Rectangle(
            (ox + b.x0 * 1e3, oy + b.y0 * 1e3), b.w * 1e3, b.h * 1e3,
            fill=False, edgecolor=edgecolor, linewidth=0.9, linestyle="--"))
        ax.text(
            ox + (b.x0 + 0.5 * b.w) * 1e3,
            oy + (b.y0 + 0.5 * b.h) * 1e3,
            b.name, color=edgecolor, ha="center", va="center", fontsize=8)


def plot_triplet(X, Y, A, B, dA, path, titles, unit, cmap, vlim, elims, case):
    if not HAS_MPL:
        return
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.2))
    fields = (A, B, dA)
    cmaps = (cmap, cmap, "RdBu_r")
    lims = (vlim, vlim, elims)
    for ax, F, cm, lim, title in zip(axes, fields, cmaps, lims, titles):
        im = ax.pcolormesh(X, Y, F, cmap=cm, shading="auto", vmin=lim[0], vmax=lim[1])
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(unit)
        _draw_chips(ax, case, edgecolor="white" if cm != "RdBu_r" else "k")
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


def _probe_at(name, x, y, z, co, py_fields, grid):
    dist = (co["x"] - x) ** 2 + (co["y"] - y) ** 2 + (co["z"] - z) ** 2
    ic = int(dist.argmin())
    xa = np.array([x])
    ya = np.array([y])
    za = np.array([z])
    row = {
        "name": name,
        "x_mm": x * 1e3,
        "y_mm": y * 1e3,
        "z_mm": z * 1e3,
        "T_py": float(sample_xy(py_fields["T"], grid, xa, ya, za)[0] - 273.15),
        "T_co": float(co["T"][ic] - 273.15),
        "vm_py": float(sample_xy(py_fields["vm"], grid, xa, ya, za)[0] * 1e-6),
        "vm_co": float(co["vm_MPa"][ic]),
    }
    row["dT_K"] = row["T_py"] - row["T_co"]
    row["dvm_MPa"] = row["vm_py"] - row["vm_co"]
    if row["vm_co"] != 0:
        row["dvm_pct"] = 100.0 * row["dvm_MPa"] / abs(row["vm_co"])
    else:
        row["dvm_pct"] = None
    for key, field, cok in (
        ("sx", py_fields["sxx"], "sx_MPa"),
        ("sy", py_fields["syy"], "sy_MPa"),
        ("sz", py_fields["szz"], "sz_MPa"),
    ):
        if cok not in co:
            continue
        pv = float(sample_xy(field, grid, xa, ya, za)[0] * 1e-6)
        cv = float(co[cok][ic])
        row[key + "_py"] = pv
        row[key + "_co"] = cv
        row["d" + key + "_MPa"] = pv - cv
        if abs(cv) > 1e-9:
            row["d" + key + "_pct"] = 100.0 * (pv - cv) / abs(cv)
        else:
            row["d" + key + "_pct"] = None
    return row


def pack(py, cv, mask):
    st = err_stats(py[mask], cv[mask])
    st["mean_abs_co"] = float(np.mean(np.abs(cv[mask])))
    st["rel_mae_pct"] = _rel_pct(st["mae"], st["mean_abs_co"])
    return st


def main():
    csv_path = os.path.join(
        ROOT, "results", "python_case3_chips_stress", "case3_interposer_uz.csv"
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

    case_tag = os.path.splitext(os.path.basename(case_path))[0]
    out = os.path.join(ROOT, "results", "compare_" + case_tag + "_interposer")
    _ensure(out)

    co = load_stress_csv(csv_path)
    print("COMSOL n=%d  z=%.3f..%.3f mm  T=%.2f..%.2f C  vm=%.1f..%.1f MPa" % (
        co["T"].size, co["z"].min() * 1e3, co["z"].max() * 1e3,
        co["T"].min() - 273.15, co["T"].max() - 273.15,
        co["vm_MPa"].min(), co["vm_MPa"].max()))
    print("CSV columns: T, vm, sx, sy, sz  (no uz in this file)")

    case = load_case(case_path)
    print("Case:", case.name)
    print("\n-- heat --")
    th = solve_fvm3d(case, verbose=True)
    print("\n-- elasticity --")
    el = solve_elasticity(th, T_ref=case.T_amb, verbose=True)
    grid = th["grid"]
    k, layer = interposer_k(el)
    print("Python interposer mid k=%d  z=%.3f mm" % (k, grid["zc"][k] * 1e3))

    pyT = sample_xy(th["T"], grid, co["x"], co["y"], co["z"])
    pyVm = sample_xy(el["vm"], grid, co["x"], co["y"], co["z"]) * 1e-6
    pySx = sample_xy(el["sxx"], grid, co["x"], co["y"], co["z"]) * 1e-6
    pySy = sample_xy(el["syy"], grid, co["x"], co["y"], co["z"]) * 1e-6
    pySz = sample_xy(el["szz"], grid, co["x"], co["y"], co["z"]) * 1e-6

    solid = grid.get("solid")
    i_idx = np.clip(np.searchsorted(grid["xc"], co["x"]) - 1, 0, len(grid["xc"]) - 1)
    j_idx = np.clip(np.searchsorted(grid["yc"], co["y"]) - 1, 0, len(grid["yc"]) - 1)
    ok = np.ones(co["T"].size, dtype=bool)
    if solid is not None:
        ok = ok & solid[i_idx, j_idx, k]
    if layer is not None:
        ok = ok & (
            (co["x"] >= layer.x0 - 2e-4) & (co["x"] <= layer.x0 + layer.w + 2e-4)
            & (co["y"] >= layer.y0 - 2e-4) & (co["y"] <= layer.y0 + layer.h + 2e-4)
        )
    interior = ok & _interior(co["x"], co["y"], layer, 1.0e-3)

    blocks = [
        ("T_C", pyT - 273.15, co["T"] - 273.15, "K"),
        ("vm", pyVm, co["vm_MPa"], "MPa"),
        ("sx", pySx, co["sx_MPa"], "MPa"),
        ("sy", pySy, co["sy_MPa"], "MPa"),
        ("sz", pySz, co["sz_MPa"], "MPa"),
    ]
    stats = {
        "case": case.name,
        "comsol_csv": csv_path,
        "plane": "interposer_mid",
        "z_mm": float(grid["zc"][k] * 1e3),
        "n_ok": int(ok.sum()),
        "n_interior": int(interior.sum()),
    }

    print("\nInterposer mid-plane  all solid n=%d   interior (>1 mm from edge) n=%d" % (
        int(ok.sum()), int(interior.sum())))
    print("%-6s  %-10s %8s %8s %8s %10s %12s" % (
        "field", "region", "MAE", "RMSE", "max|d|", "mean|CO|", "MAE/|CO|"))
    for name, py, cv, unit in blocks:
        for region, mask in (("all", ok), ("interior", interior)):
            st = pack(py, cv, mask)
            stats["%s_%s" % (name, region)] = st
            rel = st["rel_mae_pct"]
            if rel is None:
                rel_s = "n/a"
            else:
                rel_s = "%.1f%%" % rel
            print(" %-5s  %-10s %8.3f %8.3f %8.3f %10.3f %12s   %s" % (
                name, region, st["mae"], st["rmse"], st["max_abs"],
                st["mean_abs_co"], rel_s, unit))

    py_fields = {
        "T": th["T"], "vm": el["vm"],
        "sxx": el["sxx"], "syy": el["syy"], "szz": el["szz"],
    }
    probes = []
    pts = [
        ("Interposer_center", 0.039, 0.037, 0.00105),
        ("CPU_foot", 0.031, 0.031, 0.00105),
        ("GPU_foot", 0.046, 0.031, 0.00105),
        ("HBM_foot", 0.030, 0.044, 0.00105),
        ("IP_SW", 0.026, 0.026, 0.00105),
        ("IP_SE", 0.052, 0.026, 0.00105),
        ("IP_NW", 0.026, 0.048, 0.00105),
        ("IP_NE", 0.052, 0.048, 0.00105),
    ]
    print("\n%-18s %7s %7s %7s %8s %8s %8s %8s" % (
        "pt", "T_py", "T_co", "dT", "vm_py", "vm_co", "dvm", "dvm%"))
    for name, x, y, z in pts:
        r = _probe_at(name, x, y, z, co, py_fields, grid)
        probes.append(r)
        pct = r.get("dvm_pct")
        pct_s = "n/a" if pct is None else ("%.1f" % pct)
        print(" %-18s %7.2f %7.2f %7.2f %8.2f %8.2f %8.2f %8s" % (
            r["name"], r["T_py"], r["T_co"], r["dT_K"],
            r["vm_py"], r["vm_co"], r["dvm_MPa"], pct_s))
    stats["probes"] = probes
    warp = package_warpage(el)
    stats["python_warpage_pv_um"] = warp["w_pv_um"]
    stats["note"] = (
        "CSV has no uz; warpage PV cannot be compared. "
        "Relative error = MAE / mean(|COMSOL|). Interior = >1 mm from interposer edge. "
        "Edge max is a free-edge singularity."
    )

    if layer is not None and HAS_MPL:
        x0 = layer.x0 * 1e3
        x1 = (layer.x0 + layer.w) * 1e3
        y0 = layer.y0 * 1e3
        y1 = (layer.y0 + layer.h) * 1e3
        gx = np.linspace(x0, x1, 200)
        gy = np.linspace(y0, y1, 170)
        GX, GY = np.meshgrid(gx, gy, indexing="ij")
        pts_xy = np.column_stack((co["x"][ok] * 1e3, co["y"][ok] * 1e3))
        zz = co["z"]

        def grid_pair(py_field, co_arr, scale, offset):
            cg = griddata(pts_xy, co_arr[ok], (GX, GY), method="linear")
            pg = sample_xy(
                py_field, grid, GX.ravel() * 1e-3, GY.ravel() * 1e-3, zz
            ).reshape(GX.shape) * scale + offset
            pg[np.isnan(cg)] = np.nan
            return pg, cg

        Tpy_g, Tco_g = grid_pair(th["T"], co["T"], 1.0, -273.15)
        Vpy_g, Vco_g = grid_pair(el["vm"], co["vm_MPa"], 1e-6, 0.0)
        dT = Tpy_g - Tco_g
        dV = Vpy_g - Vco_g
        tmin = np.nanmin([np.nanmin(Tpy_g), np.nanmin(Tco_g)])
        tmax = np.nanmax([np.nanmax(Tpy_g), np.nanmax(Tco_g)])
        vmax = float(np.nanpercentile(np.abs(np.concatenate([
            Vpy_g[np.isfinite(Vpy_g)], Vco_g[np.isfinite(Vco_g)]
        ])), 98))
        vmax = max(vmax, 1.0)
        eT = float(np.nanmax(np.abs(dT)))
        eV = max(float(np.nanpercentile(np.abs(dV[np.isfinite(dV)]), 98)), 1.0)
        plot_triplet(
            GX, GY, Tpy_g, Tco_g, dT,
            os.path.join(out, "compare_T.png"),
            ("Python T", "COMSOL T", "Python - COMSOL"),
            "C", "inferno", (tmin, tmax), (-eT, eT), case)
        plot_triplet(
            GX, GY, Vpy_g, Vco_g, dV,
            os.path.join(out, "compare_mises.png"),
            ("Python vm", "COMSOL vm", "Python - COMSOL"),
            "MPa", "inferno", (0.0, vmax), (-eV, eV), case)
        rel_vm = stats["vm_all"]["rel_mae_pct"] or 0.0
        plot_scatter(
            pyT[ok] - 273.15, co["T"][ok] - 273.15,
            os.path.join(out, "scatter_T.png"),
            "COMSOL T (C)", "Python T (C)",
            "T  MAE=%.3f K" % stats["T_C_all"]["mae"])
        plot_scatter(
            pyVm[ok], co["vm_MPa"][ok],
            os.path.join(out, "scatter_mises.png"),
            "COMSOL vm (MPa)", "Python vm (MPa)",
            "vm MAE=%.2f MPa (%.1f%% of mean|CO|)" % (stats["vm_all"]["mae"], rel_vm))
        arrs = {"sx": pySx, "sy": pySy, "sz": pySz}
        fields = {"sx": el["sxx"], "sy": el["syy"], "sz": el["szz"]}
        coks = {"sx": "sx_MPa", "sy": "sy_MPa", "sz": "sz_MPa"}
        for key in ("sx", "sy", "sz"):
            pg, cg = grid_pair(fields[key], co[coks[key]], 1e-6, 0.0)
            d = pg - cg
            finite = np.concatenate([pg[np.isfinite(pg)], cg[np.isfinite(cg)]])
            lim = max(float(np.nanpercentile(np.abs(finite), 98)), 1.0)
            ed = max(float(np.nanpercentile(np.abs(d[np.isfinite(d)]), 98)), 1.0)
            plot_triplet(
                GX, GY, pg, cg, d,
                os.path.join(out, "compare_%s.png" % key),
                ("Python %s" % key, "COMSOL %s" % key, "Python - COMSOL"),
                "MPa", "coolwarm", (-lim, lim), (-ed, ed), case)
            rel = stats[key + "_all"]["rel_mae_pct"] or 0.0
            plot_scatter(
                arrs[key][ok], co[coks[key]][ok],
                os.path.join(out, "scatter_%s.png" % key),
                "COMSOL %s (MPa)" % key, "Python %s (MPa)" % key,
                "%s MAE=%.2f MPa (%.1f%%)" % (key, stats[key + "_all"]["mae"], rel))

    with open(os.path.join(out, "compare_stats.json"), "w") as f:
        json.dump(stats, f, indent=2)
    with open(os.path.join(out, "compare_probes.csv"), "w") as f:
        fields = [
            "name", "x_mm", "y_mm", "z_mm", "T_py", "T_co", "dT_K",
            "vm_py", "vm_co", "dvm_MPa", "dvm_pct",
            "sx_py", "sx_co", "dsx_MPa", "dsx_pct",
            "sy_py", "sy_co", "dsy_MPa", "dsy_pct",
            "sz_py", "sz_co", "dsz_MPa", "dsz_pct",
        ]
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in probes:
            w.writerow(r)
    print("\nWrote", out)


if __name__ == "__main__":
    main()
