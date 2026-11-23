# -*- coding: utf-8 -*-
"""Case 3 one-way thermoelasticity from a chips[] JSON.

Same FVM heat as run_case3.py, then atsim/elasticity.py.
Layout still comes only from chips[]; do not edit fvm3d.py.
Package warpage is uz on the interposer mid-plane (ATMPlace-style w(r)).
"""

from __future__ import print_function

import csv
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from atsim.geometry import load_case
from atsim.fvm3d import solve_fvm3d
from atsim.elasticity import solve_elasticity
from atsim.visualize import die_mid_k, plot_midplane, save_probes_csv
from atsim.warpage import package_warpage, save_warpage_csv, warpage_for_json

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


def _chip_window_mm(case, pad_mm=4.0):
    ox = getattr(case, "chip_ox", 0.0) * 1e3
    oy = getattr(case, "chip_oy", 0.0) * 1e3
    x0 = min(ox + b.x0 * 1e3 for b in case.blocks) - pad_mm
    x1 = max(ox + (b.x0 + b.w) * 1e3 for b in case.blocks) + pad_mm
    y0 = min(oy + b.y0 * 1e3 for b in case.blocks) - pad_mm
    y1 = max(oy + (b.y0 + b.h) * 1e3 for b in case.blocks) + pad_mm
    return x0, x1, y0, y1


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
            b.name, color=edgecolor, ha="center", va="center", fontsize=8,
        )


def plot_overview(th, el, path):
    """2x2 chip-cluster view: T, von Mises, uz, |u| on the die mid-plane."""
    if not HAS_MPL:
        return
    grid = th["grid"]
    case = th["case"]
    k = die_mid_k(th)
    solid = grid.get("solid")
    mask = None if solid is None else ~solid[:, :, k]
    T = np.ma.array(th["T"][:, :, k] - 273.15, mask=mask)
    vm = np.ma.array(el["vm"][:, :, k] * 1e-6, mask=mask)
    uz = np.ma.array(el["ucz"][:, :, k] * 1e6, mask=mask)
    um = np.ma.array(el["umag"][:, :, k] * 1e6, mask=mask)
    vm_flat = np.array(vm.compressed(), dtype=float)
    vm_vmax = float(np.percentile(vm_flat, 98)) if vm_flat.size else 1.0
    vm_vmax = max(vm_vmax, 1.0)
    xe, ye = grid.get("x_edges"), grid.get("y_edges")
    XE, YE = np.meshgrid(xe * 1e3, ye * 1e3, indexing="ij")
    x0, x1, y0, y1 = _chip_window_mm(case)
    panels = (
        (T, "inferno", "T (C)", "Temperature  Tmax=%.2f C" % float(T.max())),
        (vm, "inferno", "MPa", "von Mises  p98=%.1f  max=%.1f MPa" % (
            vm_vmax, float(vm.max()))),
        (uz, "coolwarm", "um", "uz  max=%.2f um" % float(np.abs(uz).max())),
        (um, "viridis", "um", "|u|  max=%.2f um" % float(um.max())),
    )
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.4))
    for ax, (F, cmap, unit, title) in zip(axes.ravel(), panels):
        kw = {}
        if unit == "MPa":
            kw["vmin"] = 0.0
            kw["vmax"] = vm_vmax
        im = ax.pcolormesh(XE, YE, F, cmap=cmap, shading="flat", **kw)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(unit)
        _draw_chips(ax, case)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title(title)
    fig.suptitle("%s die mid-plane (cropped to chips)" % case.name, y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_xz(th, path, y_mm=None):
    if not HAS_MPL:
        return
    T = th["T"]
    grid = th["grid"]
    case = th["case"]
    if y_mm is None:
        gpu = None
        for b in case.blocks:
            if b.name.upper() == "GPU":
                gpu = b
                break
        blk = gpu or max(case.blocks, key=lambda b: b.power)
        y_cut = blk.y0 + 0.5 * blk.h
    else:
        y_cut = float(y_mm) * 1e-3
    j = int(np.argmin(np.abs(grid["yc"] - y_cut)))
    Tc = T[:, j, :] - 273.15
    solid = grid.get("solid")
    if solid is not None:
        Tc = np.ma.array(Tc, mask=~solid[:, j, :])
    xe, ze = grid["x_edges"], grid["z_edges"]
    XE, ZE = np.meshgrid(xe * 1e3, ze * 1e3, indexing="ij")
    x0, x1, _, _ = _chip_window_mm(case, pad_mm=8.0)
    fig, ax = plt.subplots(figsize=(9.2, 4.4))
    im = ax.pcolormesh(XE, ZE, Tc, cmap="inferno", shading="flat")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Temperature (C)")
    for L in case.stack:
        ax.plot(
            [L.x0 * 1e3, (L.x0 + L.w) * 1e3, (L.x0 + L.w) * 1e3, L.x0 * 1e3, L.x0 * 1e3],
            [L.z0 * 1e3, L.z0 * 1e3, L.z1 * 1e3, L.z1 * 1e3, L.z0 * 1e3],
            color="white", linewidth=0.6, linestyle="--",
        )
    ax.set_xlim(x0, x1)
    ax.set_ylim(0.0, case.Lz * 1e3)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title("xz at y=%.1f mm  Tmax=%.2f C" % (grid["yc"][j] * 1e3, float(np.nanmax(Tc))))
    ax.set_aspect("auto")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_probe_bars(rows, path):
    if not HAS_MPL or not rows:
        return
    chip_rows = [r for r in rows if r["name"].endswith("_center") and r["name"] not in (
        "Sink_top_center", "Spreader_center", "Interposer_center",
        "Substrate_center", "TIM_mid", "Chip_center",
    )]
    if not chip_rows:
        chip_rows = rows[: min(8, len(rows))]
    names = [r["name"].replace("_center", "") for r in chip_rows]
    T = [r["T_C"] for r in chip_rows]
    vm = [r["vm_MPa"] for r in chip_rows]
    uz = [r["uz_um"] for r in chip_rows]
    x = np.arange(len(names))
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.6))
    axes[0].bar(x, T, color="#d35400")
    axes[0].set_title("T (C)")
    axes[1].bar(x, vm, color="#8e44ad")
    axes[1].set_title("von Mises (MPa)")
    axes[2].bar(x, uz, color="#2980b9")
    axes[2].set_title("uz (um)")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_field(grid, case, field, path, title, cmap="inferno", scale=1.0, unit="", mask_solid=True):
    if not HAS_MPL:
        return
    k = die_mid_k({"grid": grid, "case": case, "T": field})
    F = field[:, :, k] * scale
    solid = grid.get("solid")
    if mask_solid and solid is not None:
        F = np.ma.array(F, mask=~solid[:, :, k])
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    xe, ye = grid.get("x_edges"), grid.get("y_edges")
    if xe is not None and ye is not None:
        XE, YE = np.meshgrid(xe * 1e3, ye * 1e3, indexing="ij")
        im = ax.pcolormesh(XE, YE, F, cmap=cmap, shading="flat")
        ax.set_aspect("equal")
        ax.set_xlim(0.0, case.Lx * 1e3)
        ax.set_ylim(0.0, case.Ly * 1e3)
    else:
        im = ax.imshow(
            F.T, origin="lower",
            extent=[0, case.Lx * 1e3, 0, case.Ly * 1e3],
            cmap=cmap, aspect="equal",
        )
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    if unit:
        cb.set_label(unit)
    ox = getattr(case, "chip_ox", 0.0) * 1e3
    oy = getattr(case, "chip_oy", 0.0) * 1e3
    for b in case.blocks:
        ax.add_patch(plt.Rectangle(
            (ox + b.x0 * 1e3, oy + b.y0 * 1e3), b.w * 1e3, b.h * 1e3,
            fill=False, edgecolor="white", linewidth=0.8, linestyle="--"))
        ax.text(
            ox + (b.x0 + 0.5 * b.w) * 1e3,
            oy + (b.y0 + 0.5 * b.h) * 1e3,
            b.name, color="white", ha="center", va="center", fontsize=8,
        )
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def sample_mech(elast):
    grid = elast["grid"]
    case = elast["case"]
    xc, yc, zc = grid["xc"], grid["yc"], grid["zc"]
    rows = []
    for p in case.probes:
        i = int(np.argmin(np.abs(xc - p.x)))
        j = int(np.argmin(np.abs(yc - p.y)))
        k = int(np.argmin(np.abs(zc - p.z)))
        rows.append({
            "name": p.name,
            "x_m": p.x, "y_m": p.y, "z_m": p.z,
            "T_C": float(elast["thermal"]["T"][i, j, k] - 273.15),
            "ux_um": float(elast["ucx"][i, j, k] * 1e6),
            "uy_um": float(elast["ucy"][i, j, k] * 1e6),
            "uz_um": float(elast["ucz"][i, j, k] * 1e6),
            "umag_um": float(elast["umag"][i, j, k] * 1e6),
            "vm_MPa": float(elast["vm"][i, j, k] * 1e-6),
            "sxx_MPa": float(elast["sxx"][i, j, k] * 1e-6),
            "syy_MPa": float(elast["syy"][i, j, k] * 1e-6),
        })
    return rows


def plot_warpage(elast, warp, path):
    if not HAS_MPL:
        return
    grid = elast["grid"]
    case = elast["case"]
    k = warp["k"]
    uz = elast["ucz"][:, :, k] * 1e6
    F = np.ma.array(uz, mask=~warp["mask"])
    xe, ye = grid.get("x_edges"), grid.get("y_edges")
    XE, YE = np.meshgrid(xe * 1e3, ye * 1e3, indexing="ij")
    x0, x1, y0, y1 = _chip_window_mm(case, pad_mm=6.0)
    vmax = max(abs(warp["w_max_um"]), abs(warp["w_min_um"]), 0.05)
    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    im = ax.pcolormesh(XE, YE, F, cmap="coolwarm", shading="flat", vmin=-vmax, vmax=vmax)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("w = uz (um)")
    _draw_chips(ax, case, edgecolor="k")
    ax.plot(warp["xy_max_mm"][0], warp["xy_max_mm"][1], "k^", markersize=7, label="w max")
    ax.plot(warp["xy_min_mm"][0], warp["xy_min_mm"][1], "kv", markersize=7, label="w min")
    for c in warp.get("corners") or []:
        ax.plot(c["x_mm"], c["y_mm"], "o", mfc="none", mec="k", markersize=6)
        ax.text(c["x_mm"], c["y_mm"], " %.2f" % c["uz_um"], fontsize=7, color="k")
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("package warpage  PV=%.3f um  (interposer z=%.3f mm)" % (
        warp["w_pv_um"], warp["z_mm"]))
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    case_path = os.path.join(ROOT, "cases", "case3_chips.json")
    if len(sys.argv) >= 2:
        case_path = sys.argv[1]
        if not os.path.isabs(case_path):
            case_path = os.path.join(ROOT, case_path)
    name = os.path.splitext(os.path.basename(case_path))[0]
    out = os.path.join(ROOT, "results", "python_" + name + "_stress")
    _ensure(out)

    case = load_case(case_path)
    meta = case.raw.get("_layout_meta") or {}
    mech = case.raw.get("mechanics") or {}
    print("=" * 64)
    print("Case 3 thermoelasticity (one-way)")
    print("JSON:", case_path)
    print("Chips: %d   P=%.2f W   h=%.2f W/m2K" % (
        meta.get("n_chips", len(case.blocks)), case.total_power, case.h_top))
    print("T_ref = Tamb = %.2f C" % (case.T_amb - 273.15))
    print("BC: %s  (no TSV / C4; interfaces bonded)" % (
        mech.get("bc", "substrate_bottom_fixed"),))
    print("Materials:")
    for mat_name, mat in case.materials.items():
        print("  %-8s  E=%.2e Pa  nu=%.2f  alpha=%.2e 1/K" % (
            mat_name, mat.E or 0.0, mat.nu, mat.alpha_cte))
    print("\nChips:")
    print("  %-10s %8s %8s %8s %8s %8s" % ("name", "x_mm", "y_mm", "w_mm", "h_mm", "W"))
    for b in case.blocks:
        print("  %-10s %8.2f %8.2f %8.2f %8.2f %8.2f" % (
            b.name, b.x0 * 1e3, b.y0 * 1e3, b.w * 1e3, b.h * 1e3, b.power))

    t_wall0 = time.time()
    print("\n-- heat --")
    th = solve_fvm3d(case, verbose=True)
    save_probes_csv(th, os.path.join(out, "thermal_probes.csv"))
    plot_midplane(th, os.path.join(out, "T_die.png"))

    print("\n-- elasticity --")
    el = solve_elasticity(th, T_ref=case.T_amb, verbose=True)
    t_compute = time.time() - t_wall0
    solid = el["grid"]["solid"]
    umax = float(el["umag"][solid].max())
    uzmax = float(np.abs(el["ucz"][solid]).max())
    vmmax = float(el["vm"][solid].max())
    warp = package_warpage(el)

    rows = sample_mech(el)
    with open(os.path.join(out, "mech_probes.csv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    save_warpage_csv(el, warp, os.path.join(out, "warpage_interposer.csv"))

    print("\nMechanical probes:")
    print("  %-20s %8s %10s %10s %10s" % ("name", "T C", "uz um", "|u| um", "vm MPa"))
    for r in rows:
        print("  %-20s %8.2f %10.4f %10.4f %10.3f" % (
            r["name"], r["T_C"], r["uz_um"], r["umag_um"], r["vm_MPa"]))

    plot_field(el["grid"], case, el["umag"], os.path.join(out, "umag_die.png"),
               " |u| die mid  max=%.2f um" % (umax * 1e6), cmap="viridis",
               scale=1e6, unit="um")
    plot_field(el["grid"], case, el["ucz"], os.path.join(out, "uz_die.png"),
               " uz die mid  max=%.2f um" % (uzmax * 1e6), cmap="coolwarm",
               scale=1e6, unit="um")
    plot_field(el["grid"], case, el["vm"], os.path.join(out, "vm_die.png"),
               " von Mises die mid  max=%.1f MPa" % (vmmax * 1e-6), cmap="inferno",
               scale=1e-6, unit="MPa")
    plot_overview(th, el, os.path.join(out, "overview_chips.png"))
    plot_xz(th, os.path.join(out, "xz_gpu.png"))
    plot_probe_bars(rows, os.path.join(out, "probe_bars.png"))
    plot_warpage(el, warp, os.path.join(out, "warpage_interposer.png"))

    print("\nPackage warpage (interposer mid, ATMPlace-style w=uz):")
    print("  plane z=%.3f mm  n=%d" % (warp["z_mm"], warp["n_cells"]))
    print("  w min/mean/max = %.4f / %.4f / %.4f um" % (
        warp["w_min_um"], warp["w_mean_um"], warp["w_max_um"]))
    print("  peak-to-valley  %.4f um   |w|_max %.4f um" % (
        warp["w_pv_um"], warp["w_abs_max_um"]))
    if warp.get("corners"):
        print("  corners (um): " + "  ".join(
            "%s=%.3f" % (c["name"], c["uz_um"]) for c in warp["corners"]))
        print("  corner PV     %.4f um" % warp["corner_pv_um"])

    summary = {
        "case": case.name,
        "json": case_path,
        "T_ref_C": case.T_amb - 273.15,
        "Tmax_C": float(th["T"][solid].max() - 273.15),
        "ndof": el["ndof"],
        "u_max_um": umax * 1e6,
        "uz_max_um": uzmax * 1e6,
        "vm_max_MPa": vmmax * 1e-6,
        "warpage": warpage_for_json(warp),
        "heat_s": th["elapsed_s"],
        "elast_s": el["elapsed_s"],
        "compute_s": t_compute,
        "heat_timing": th.get("timing"),
        "elast_timing": el.get("timing"),
        "probes": rows,
        "bc": el.get("bc", "substrate_bottom_fixed"),
        "coupling": "one-way thermal to mechanical",
        "layout_meta": meta,
        "chips": [
            {
                "name": b.name,
                "x_mm": b.x0 * 1e3,
                "y_mm": b.y0 * 1e3,
                "w_mm": b.w * 1e3,
                "h_mm": b.h * 1e3,
                "power_W": b.power,
            }
            for b in case.blocks
        ],
    }
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    ht = th.get("timing") or {"total_s": th["elapsed_s"]}
    et = el.get("timing") or {"total_s": el["elapsed_s"]}
    timing_out = {
        "case": case.name,
        "json": case_path,
        "wall_compute_s": t_compute,
        "heat": ht,
        "elasticity": et,
        "total_heat_plus_elast_s": float(th["elapsed_s"]) + float(el["elapsed_s"]),
    }
    with open(os.path.join(out, "timing.json"), "w") as f:
        json.dump(timing_out, f, indent=2)
    print("\nTiming:")
    print("  heat   assemble %.2fs  solve %.2fs (%s)  total %.2fs  ndof=%s" % (
        ht.get("assemble_s", 0.0), ht.get("solve_s", 0.0), ht.get("solver", "?"),
        ht.get("total_s", th["elapsed_s"]), ht.get("n_unknown")))
    print("  elast  assemble %.2fs  solve %.2fs (%s)  post %.2fs  total %.2fs  ndof=%s" % (
        et.get("assemble_s", 0.0), et.get("solve_s", 0.0), et.get("solver", "?"),
        et.get("post_s", 0.0), et.get("total_s", el["elapsed_s"]), et.get("ndof")))
    print("  heat+elast %.2fs  (plots not included in this sum; wall compute %.2fs)" % (
        float(th["elapsed_s"]) + float(el["elapsed_s"]), t_compute))
    print("  wrote", os.path.join(out, "timing.json"))
    print("\nOutputs:", out)
    print("=" * 64)


if __name__ == "__main__":
    main()
