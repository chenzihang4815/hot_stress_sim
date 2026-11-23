# -*- coding: utf-8 -*-
"""Pick representative COMSOL nodes with modest error vs Python.

Optional. Not on the Case 3 solve path.

    python scripts/pick_compare_points.py
"""

from __future__ import print_function

import csv
import json
import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from atsim.geometry import load_case
from atsim.fvm3d import solve_fvm3d
from atsim.elasticity import solve_elasticity
from atsim.warpage import interposer_k

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
            if (not s) or s.startswith("%") or s.startswith("#"):
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
    out = {"x": x, "y": y, "z": z, "T": T, "vm_MPa": _mpa(vm), "path": path}
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


def chip_rects_m(case):
    ox = getattr(case, "chip_ox", 0.0)
    oy = getattr(case, "chip_oy", 0.0)
    recs = []
    for b in case.blocks:
        recs.append({
            "name": b.name,
            "x0": ox + b.x0,
            "y0": oy + b.y0,
            "x1": ox + b.x0 + b.w,
            "y1": oy + b.y0 + b.h,
            "cx": ox + b.x0 + 0.5 * b.w,
            "cy": oy + b.y0 + 0.5 * b.h,
        })
    return recs


def in_rect(x, y, r, margin):
    return (
        (x >= r["x0"] + margin) & (x <= r["x1"] - margin)
        & (y >= r["y0"] + margin) & (y <= r["y1"] - margin)
    )


def eval_nodes(co, th, el, grid):
    pyT = sample_xy(th["T"], grid, co["x"], co["y"], co["z"]) - 273.15
    coT = co["T"] - 273.15
    pyVm = sample_xy(el["vm"], grid, co["x"], co["y"], co["z"]) * 1e-6
    pySx = sample_xy(el["sxx"], grid, co["x"], co["y"], co["z"]) * 1e-6
    pySy = sample_xy(el["syy"], grid, co["x"], co["y"], co["z"]) * 1e-6
    pySz = sample_xy(el["szz"], grid, co["x"], co["y"], co["z"]) * 1e-6
    return {
        "T_py": pyT, "T_co": coT, "dT": pyT - coT,
        "vm_py": pyVm, "vm_co": co["vm_MPa"], "dvm": pyVm - co["vm_MPa"],
        "sx_py": pySx, "sx_co": co["sx_MPa"], "dsx": pySx - co["sx_MPa"],
        "sy_py": pySy, "sy_co": co["sy_MPa"], "dsy": pySy - co["sy_MPa"],
        "sz_py": pySz, "sz_co": co["sz_MPa"], "dsz": pySz - co["sz_MPa"],
    }


def nearest_idx(co, x, y):
    d2 = (co["x"] - x) ** 2 + (co["y"] - y) ** 2
    return int(d2.argmin())


def row_from(i, name, why, co, ev):
    dvm = float(ev["dvm"][i])
    vm_co = float(ev["vm_co"][i])
    pct = None if abs(vm_co) < 1e-9 else 100.0 * dvm / abs(vm_co)
    return {
        "name": name,
        "why": why,
        "i": i,
        "x_mm": float(co["x"][i] * 1e3),
        "y_mm": float(co["y"][i] * 1e3),
        "z_mm": float(co["z"][i] * 1e3),
        "T_py": float(ev["T_py"][i]),
        "T_co": float(ev["T_co"][i]),
        "dT_K": float(ev["dT"][i]),
        "vm_py": float(ev["vm_py"][i]),
        "vm_co": vm_co,
        "dvm_MPa": dvm,
        "dvm_pct": pct,
        "sx_py": float(ev["sx_py"][i]),
        "sx_co": float(ev["sx_co"][i]),
        "dsx_MPa": float(ev["dsx"][i]),
        "sy_py": float(ev["sy_py"][i]),
        "sy_co": float(ev["sy_co"][i]),
        "dsy_MPa": float(ev["dsy"][i]),
        "sz_py": float(ev["sz_py"][i]),
        "sz_co": float(ev["sz_co"][i]),
        "dsz_MPa": float(ev["dsz"][i]),
    }


def pick_best(mask, ev, co, name, why, t_max=0.4, vm_floor=5.0):
    """Small |dvm| among interior nodes with |dT|<t_max and vm_co>vm_floor."""
    ok = mask & (np.abs(ev["dT"]) < t_max) & (ev["vm_co"] > vm_floor)
    if not np.any(ok):
        ok = mask & (ev["vm_co"] > vm_floor)
    if not np.any(ok):
        ok = mask
    if not np.any(ok):
        return None
    score = np.abs(ev["dvm"])
    score = np.where(ok, score, np.inf)
    i = int(np.argmin(score))
    return row_from(i, name, why, co, ev)


def pick_median(mask, ev, co, name, why):
    """Node whose |dvm| is closest to the interior median |dvm| (typical, not lucky)."""
    if not np.any(mask):
        return None
    ad = np.abs(ev["dvm"][mask])
    target = float(np.median(ad))
    idx = np.where(mask)[0]
    i = int(idx[np.argmin(np.abs(ad - target))])
    return row_from(i, name, why, co, ev)


def summarize_mask(ev, mask):
    if not np.any(mask):
        return None
    dvm = ev["dvm"][mask]
    dT = ev["dT"][mask]
    return {
        "n": int(mask.sum()),
        "dT_mae": float(np.abs(dT).mean()),
        "dvm_mae": float(np.abs(dvm).mean()),
        "dvm_med": float(np.median(dvm)),
        "dvm_p10": float(np.percentile(np.abs(dvm), 10)),
        "dvm_p50": float(np.percentile(np.abs(dvm), 50)),
        "rel_mae_pct": float(100.0 * np.abs(dvm).mean() / max(np.mean(np.abs(ev["vm_co"][mask])), 1e-9)),
        "mean_vm_py": float(ev["vm_py"][mask].mean()),
        "mean_vm_co": float(ev["vm_co"][mask].mean()),
    }


def print_row(r):
    pct = r["dvm_pct"]
    pct_s = "n/a" if pct is None else ("%+.1f" % pct)
    print("  %-22s  (%7.3f, %7.3f, %6.3f) mm" % (r["name"], r["x_mm"], r["y_mm"], r["z_mm"]))
    print("      T   Py %6.2f  CO %6.2f  dT %+6.3f K" % (r["T_py"], r["T_co"], r["dT_K"]))
    print("      vm  Py %6.2f  CO %6.2f  dvm %+6.2f MPa  (%s%%)" % (
        r["vm_py"], r["vm_co"], r["dvm_MPa"], pct_s))
    print("      sx  Py %6.2f  CO %6.2f  dsx %+6.2f   sy %+6.2f vs %+6.2f" % (
        r["sx_py"], r["sx_co"], r["dsx_MPa"], r["sy_py"], r["sy_co"]))
    print("      %s" % r["why"])


def plot_picks(co, ev, recs, picks, path, title, z_mm):
    if not HAS_MPL or not picks:
        return
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    sc = ax.scatter(
        co["x"] * 1e3, co["y"] * 1e3,
        c=ev["dvm"], s=6, cmap="RdBu_r", vmin=-8, vmax=8, linewidths=0)
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("vm Python - COMSOL (MPa)")
    for r in recs:
        ax.add_patch(plt.Rectangle(
            (r["x0"] * 1e3, r["y0"] * 1e3),
            (r["x1"] - r["x0"]) * 1e3, (r["y1"] - r["y0"]) * 1e3,
            fill=False, edgecolor="k", linewidth=0.9, linestyle="--"))
        ax.text(r["cx"] * 1e3, r["cy"] * 1e3, r["name"], ha="center", va="center", fontsize=8)
    for p in picks:
        ax.plot(p["x_mm"], p["y_mm"], "ko", markersize=5, mfc="yellow")
        ax.text(p["x_mm"] + 0.25, p["y_mm"] + 0.25, p["name"], fontsize=7, color="k")
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("%s  z=%.3f mm" % (title, z_mm))
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    case_path = os.path.join(ROOT, "cases", "case3_chips.json")
    die_csv = os.path.join(ROOT, "results", "python_case3_chips_stress", "case3_die_points_stress.csv")
    ip_csv = os.path.join(ROOT, "results", "python_case3_chips_stress", "case3_interposer_uz.csv")
    out = os.path.join(ROOT, "results", "compare_case3_picked_points")
    if not os.path.isdir(out):
        os.makedirs(out)

    case = load_case(case_path)
    recs = chip_rects_m(case)
    print("Case:", case.name)
    print("-- heat --")
    th = solve_fvm3d(case, verbose=True)
    print("-- elasticity --")
    el = solve_elasticity(th, T_ref=case.T_amb, verbose=True)
    grid = th["grid"]
    k_ip, layer = interposer_k(el)
    print("interposer k=%d z=%.3f mm" % (k_ip, grid["zc"][k_ip] * 1e3))

    die = load_stress_csv(die_csv)
    ip = load_stress_csv(ip_csv)
    ev_die = eval_nodes(die, th, el, grid)
    ev_ip = eval_nodes(ip, th, el, grid)

    picks = []
    stats = {}

    print("\n========== DIE mid-plane  z~1.175 mm  file=%s ==========" % os.path.basename(die_csv))
    print("n=%d  z=%.3f..%.3f mm" % (die["T"].size, die["z"].min() * 1e3, die["z"].max() * 1e3))
    for r in recs:
        core = in_rect(die["x"], die["y"], r, 1.5e-3)
        full = in_rect(die["x"], die["y"], r, 0.3e-3)
        st = summarize_mask(ev_die, core)
        stats["die_core_" + r["name"]] = st
        print("\n%s interior (margin 1.5 mm) n=%d  dT MAE=%.3f K  |dvm| med=%.2f MAE=%.2f MPa  rel=%.1f%%  mean vm Py/CO=%.1f/%.1f" % (
            r["name"], st["n"], st["dT_mae"], st["dvm_p50"], st["dvm_mae"],
            st["rel_mae_pct"], st["mean_vm_py"], st["mean_vm_co"]))
        ic = nearest_idx(die, r["cx"], r["cy"])
        picks.append(row_from(
            ic, r["name"] + "_center",
            "chip geometric center, die mid z=1.175 mm", die, ev_die))
        b = pick_best(
            core, ev_die, die, r["name"] + "_best_interior",
            "smallest |dvm| among nodes >=1.5 mm from this chip edge, |dT|<0.4 K, vm_co>5 MPa")
        if b is not None:
            picks.append(b)
        m = pick_median(
            core, ev_die, die, r["name"] + "_typical_interior",
            "interior node whose |dvm| is the median of this chip core (typical, not the luckiest)")
        if m is not None:
            picks.append(m)

    print("\n========== INTERPOSER mid-plane  z~1.050 mm  file=%s ==========" % os.path.basename(ip_csv))
    print("n=%d  z=%.3f..%.3f mm" % (ip["T"].size, ip["z"].min() * 1e3, ip["z"].max() * 1e3))
    for r in recs:
        core = in_rect(ip["x"], ip["y"], r, 1.5e-3)
        st = summarize_mask(ev_ip, core)
        stats["ip_core_" + r["name"]] = st
        print("\n%s footprint interior n=%d  dT MAE=%.3f K  |dvm| med=%.2f MAE=%.2f MPa  rel=%.1f%%  mean vm Py/CO=%.1f/%.1f" % (
            r["name"], st["n"], st["dT_mae"], st["dvm_p50"], st["dvm_mae"],
            st["rel_mae_pct"], st["mean_vm_py"], st["mean_vm_co"]))
        ic = nearest_idx(ip, r["cx"], r["cy"])
        picks.append(row_from(
            ic, r["name"] + "_foot",
            "chip center projected onto interposer mid z=1.050 mm", ip, ev_ip))
        b = pick_best(
            core, ev_ip, ip, r["name"] + "_foot_best",
            "smallest |dvm| under this chip, >=1.5 mm from chip edge, interposer mid")
        if b is not None:
            picks.append(b)

    print("\n========== SELECTED POINTS ==========")
    die_picks = [p for p in picks if abs(p["z_mm"] - 1.175) < 0.02]
    ip_picks = [p for p in picks if abs(p["z_mm"] - 1.05) < 0.02]
    print("\n-- die mid (z=1.175 mm) --")
    for p in die_picks:
        print_row(p)
        print("")
    print("-- interposer mid (z=1.050 mm) --")
    for p in ip_picks:
        print_row(p)
        print("")

    plot_picks(die, ev_die, recs, die_picks, os.path.join(out, "picked_die.png"),
               "die mid  picked points", 1.175)
    plot_picks(ip, ev_ip, recs, ip_picks, os.path.join(out, "picked_interposer.png"),
               "interposer mid  picked points", 1.050)

    fields = [
        "name", "why", "x_mm", "y_mm", "z_mm",
        "T_py", "T_co", "dT_K",
        "vm_py", "vm_co", "dvm_MPa", "dvm_pct",
        "sx_py", "sx_co", "dsx_MPa",
        "sy_py", "sy_co", "dsy_MPa",
        "sz_py", "sz_co", "dsz_MPa",
    ]
    with open(os.path.join(out, "picked_points.csv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for p in picks:
            w.writerow(p)
    with open(os.path.join(out, "picked_stats.json"), "w") as f:
        json.dump({"stats": stats, "picks": picks}, f, indent=2)
    print("Wrote", out)


if __name__ == "__main__":
    main()
