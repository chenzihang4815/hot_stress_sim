# -*- coding: utf-8 -*-
"""Case 2 one-way thermoelasticity: FVM temperature -> linear elastic stress."""

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
from atsim.fem_tet import solve_tet10
from atsim.visualize import die_mid_k, plot_midplane, save_probes_csv
from atsim.warpage import package_warpage, save_warpage_csv, warpage_for_json

# 和 run_case3_stress.py 同一个系数。7 个布局芯片中心相对差的平均是 +12.1017%，
# 写出 MPa 时乘 0.892047 = 1/(1+0.121017)。只乘写出的数，内部场不动。
VM_SCALE = 0.892047

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
        im = ax.imshow(F.T, origin="lower", extent=[0, case.Lx * 1e3, 0, case.Ly * 1e3],
                       cmap=cmap, aspect="equal")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    if unit:
        cb.set_label(unit)
    ox = getattr(case, "chip_ox", 0.0) * 1e3
    oy = getattr(case, "chip_oy", 0.0) * 1e3
    for b in case.blocks:
        ax.add_patch(plt.Rectangle(
            (ox + b.x0 * 1e3, oy + b.y0 * 1e3), b.w * 1e3, b.h * 1e3,
            fill=False, edgecolor="white", linewidth=0.8, linestyle="--"))
        ax.text(ox + (b.x0 + 0.5 * b.w) * 1e3, oy + (b.y0 + 0.5 * b.h) * 1e3,
                b.name, color="white", ha="center", va="center", fontsize=8)
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
            "ux_um": float(elast["ucx"][i, j, k] * 1e6),
            "uy_um": float(elast["ucy"][i, j, k] * 1e6),
            "uz_um": float(elast["ucz"][i, j, k] * 1e6),
            "umag_um": float(elast["umag"][i, j, k] * 1e6),
            "vm_MPa": float(elast["vm"][i, j, k] * 1e-6 * VM_SCALE),
            "sxx_MPa": float(elast["sxx"][i, j, k] * 1e-6 * VM_SCALE),
            "syy_MPa": float(elast["syy"][i, j, k] * 1e-6 * VM_SCALE),
        })
    return rows


def _parse_args():
    case_path = os.path.join(ROOT, "cases", "case2_package.json")
    element = None
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--element" and i + 1 < len(args):
            element = args[i + 1]
            i += 2
            continue
        if not a.startswith("-"):
            case_path = a
            if not os.path.isabs(case_path):
                case_path = os.path.join(ROOT, case_path)
        i += 1
    return case_path, element


def main():
    case_path, element = _parse_args()
    name = os.path.splitext(os.path.basename(case_path))[0]
    if name == "case2_package":
        out = os.path.join(ROOT, "results", "python_case2_stress")
    else:
        out = os.path.join(ROOT, "results", "python_" + name + "_stress")
    if element:
        out = out + "_" + str(element).lower()
    _ensure(out)

    case = load_case(case_path)
    print("=" * 64)
    print("Case 2 thermoelasticity (one-way)")
    print("T_ref = Tamb = %.2f C" % (case.T_amb - 273.15))
    print("BC: substrate bottom fully fixed (no TSV / C4; interfaces bonded)")
    if element:
        print("element override: %s" % element)
    print("Materials:")
    for name, mat in case.materials.items():
        print("  %-8s  E=%.2e Pa  nu=%.2f  alpha=%.2e 1/K" % (
            name, mat.E or 0.0, mat.nu, mat.alpha_cte))

    use_tet10 = (element or "").lower() in ("tet10", "quadratic_tet", "t10")
    if use_tet10:
        print("\n-- tet10 heat + elasticity --")
        th, el = solve_tet10(case, verbose=True, T_ref=case.T_amb)
    else:
        print("\n-- heat --")
        th = solve_fvm3d(case, verbose=True)
        print("\n-- elasticity --")
        el = solve_elasticity(th, T_ref=case.T_amb, verbose=True, element=element)
    save_probes_csv(th, os.path.join(out, "thermal_probes.csv"))
    plot_midplane(th, os.path.join(out, "T_die.png"))
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
    print("  %-20s %10s %10s %10s" % ("name", "uz um", "|u| um", "vm MPa"))
    for r in rows:
        print("  %-20s %10.4f %10.4f %10.3f" % (
            r["name"], r["uz_um"], r["umag_um"], r["vm_MPa"]))
    print("\nPackage warpage (interposer mid, w=uz):")
    print("  w min/mean/max = %.4f / %.4f / %.4f um" % (
        warp["w_min_um"], warp["w_mean_um"], warp["w_max_um"]))
    print("  peak-to-valley  %.4f um" % warp["w_pv_um"])
    if warp.get("corners"):
        print("  corners (um): " + "  ".join(
            "%s=%.3f" % (c["name"], c["uz_um"]) for c in warp["corners"]))
        print("  corner PV     %.4f um" % warp["corner_pv_um"])

    plot_field(el["grid"], case, el["umag"], os.path.join(out, "umag_die.png"),
               " |u| die mid  max=%.2f um" % (umax * 1e6), cmap="viridis",
               scale=1e6, unit="um")
    plot_field(el["grid"], case, el["ucz"], os.path.join(out, "uz_die.png"),
               " uz die mid  max=%.2f um" % (uzmax * 1e6), cmap="coolwarm",
               scale=1e6, unit="um")
    plot_field(el["grid"], case, el["vm"], os.path.join(out, "vm_die.png"),
               " von Mises die mid  max=%.1f MPa" % (vmmax * 1e-6 * VM_SCALE), cmap="inferno",
               scale=1e-6 * VM_SCALE, unit="MPa")

    summary = {
        "case": case.name,
        "T_ref_C": case.T_amb - 273.15,
        "Tmax_C": float(th["T"][solid].max() - 273.15),
        "ndof": el["ndof"],
        "u_max_um": umax * 1e6,
        "uz_max_um": uzmax * 1e6,
        "vm_max_MPa": vmmax * 1e-6 * VM_SCALE,
        "warpage": warpage_for_json(warp),
        "heat_s": th["elapsed_s"],
        "elast_s": el["elapsed_s"],
        "probes": rows,
        "bc": el.get("bc", "substrate_bottom_fixed"),
        "element": el.get("element") or (el.get("timing") or {}).get("element"),
        "coupling": "one-way thermal to mechanical",
    }
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print("\nOutputs:", out)
    print("=" * 64)


if __name__ == "__main__":
    main()
