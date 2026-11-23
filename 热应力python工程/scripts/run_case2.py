# -*- coding: utf-8 -*-
"""Run Case 2: 2.5D substrate + interposer + die + TIM + spreader + sink."""

from __future__ import print_function

import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from atsim.geometry import load_case
from atsim.fvm3d import solve_fvm3d
from atsim.hotspot_io import write_hotspot_case
from atsim.package_grid import die_layer
from atsim.visualize import (
    die_mid_k,
    plot_midplane,
    save_midplane_csv,
    save_numpy,
    save_probes_csv,
)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False


def plot_xz(result, path):
    """xz slice through CoreA center (y = chip_oy + 6 mm)."""
    if not HAS_MPL:
        return
    T = result["T"]
    grid = result["grid"]
    case = result["case"]
    y_cut = case.chip_oy + 0.006
    j = int((abs(grid["yc"] - y_cut)).argmin())
    Tc = T[:, j, :] - 273.15
    solid = grid.get("solid")
    if solid is not None:
        Tc = __import__("numpy").ma.array(Tc, mask=~solid[:, j, :])
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    # non-uniform axes: pcolormesh with edges
    xe = grid.get("x_edges")
    ze = grid.get("z_edges")
    if xe is None:
        xe = __import__("numpy").linspace(0, case.Lx, T.shape[0] + 1)
        ze = __import__("numpy").linspace(0, case.Lz, T.shape[2] + 1)
    import numpy as np
    XE, ZE = np.meshgrid(xe * 1e3, ze * 1e3, indexing="ij")
    im = ax.pcolormesh(XE, ZE, Tc, cmap="inferno", shading="flat")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Temperature (C)")
    for L in case.stack:
        ax.plot(
            [L.x0 * 1e3, (L.x0 + L.w) * 1e3, (L.x0 + L.w) * 1e3, L.x0 * 1e3, L.x0 * 1e3],
            [L.z0 * 1e3, L.z0 * 1e3, L.z1 * 1e3, L.z1 * 1e3, L.z0 * 1e3],
            color="white",
            linewidth=0.7,
            linestyle="--",
        )
        ax.text(
            (L.x0 + 0.5 * L.w) * 1e3,
            0.5 * (L.z0 + L.z1) * 1e3,
            L.name,
            color="white",
            ha="center",
            va="center",
            fontsize=8,
        )
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title("Case 2 xz slice at y=%.1f mm (CoreA row)  Tmax=%.2f C" % (
        y_cut * 1e3, float(T[grid["solid"]].max() - 273.15) if solid is not None else T.max() - 273.15))
    ax.set_aspect("auto")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    case_path = os.path.join(ROOT, "cases", "case2_package.json")
    out_py = os.path.join(ROOT, "results", "python_case2")
    out_hs = os.path.join(ROOT, "hotspot_case2")
    if not os.path.isdir(out_py):
        os.makedirs(out_py)

    case = load_case(case_path)
    die = die_layer(case)
    print("=" * 64)
    print("Case:", case.name)
    print("Domain: %.0f x %.0f x %.2f mm" % (case.Lx * 1e3, case.Ly * 1e3, case.Lz * 1e3))
    print("Die:    %.1f x %.1f x %.0f um  origin=(%.1f, %.1f) mm" % (
        die.w * 1e3, die.h * 1e3, die.thickness * 1e6,
        case.chip_ox * 1e3, case.chip_oy * 1e3))
    print("h_sink_top = %.2f W/m2K   Tamb=%.2f C   P=%.2f W" % (
        case.h_top, case.T_amb - 273.15, case.total_power))
    print("R_conv = 1/(h A_sink) = %.4f K/W  (HotSpot r_convec=0.1)" % (
        1.0 / (case.h_top * case.Lx * case.Ly)))
    for L in case.stack:
        print("  %-8s  z=%.0f-%.0f um  %s  %dx%dx  k=%s" % (
            L.name, L.z0 * 1e6, L.z1 * 1e6, L.material,
            int(round(L.w * 1e3)), int(round(L.h * 1e3)),
            case.material(L.material).kappa0))

    print("\n-- FVM 3D package --")
    res = solve_fvm3d(case, verbose=True)
    solid = res["grid"]["solid"]
    Tsol = res["T"][solid]
    kmid = die_mid_k(res)

    probes = save_probes_csv(res, os.path.join(out_py, "probes.csv"))
    save_midplane_csv(res, os.path.join(out_py, "die_midplane.csv"), k=kmid)
    save_numpy(res, os.path.join(out_py, "T_field.npz"))
    plot_midplane(res, os.path.join(out_py, "die_midplane.png"),
                  title="%s die mid-plane  Tmax_solid=%.2f C" % (case.name, Tsol.max() - 273.15))
    plot_xz(res, os.path.join(out_py, "xz_corea.png"))

    print("\nProbes:")
    for r in probes:
        print("  %-20s  %.3f C" % (r["name"], r["T_C"]))

    summary = {
        "case": case.name,
        "solver": "fvm3d_package",
        "grid": [case.nx, case.ny, case.nz],
        "n_solid": int(solid.sum()),
        "total_power_W": case.total_power,
        "h_W_m2K": case.h_top,
        "r_convec_K_per_W": 1.0 / (case.h_top * case.Lx * case.Ly),
        "T_amb_C": case.T_amb - 273.15,
        "Tmax_solid_C": float(Tsol.max() - 273.15),
        "Tmin_solid_C": float(Tsol.min() - 273.15),
        "Tmean_solid_C": float(Tsol.mean() - 273.15),
        "elapsed_s": res["elapsed_s"],
        "probes": probes,
        "stack": [
            {"name": L.name, "z0_mm": L.z0 * 1e3, "z1_mm": L.z1 * 1e3,
             "xy_mm": [L.w * 1e3, L.h * 1e3], "material": L.material}
            for L in case.stack
        ],
    }
    with open(os.path.join(out_py, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    hs_files = write_hotspot_case(case, out_hs)
    print("\nHotSpot files:")
    for k, v in hs_files.items():
        print("  %s: %s" % (k, v))
    print("Compile HotSpot then run hotspot_case2/run_hotspot.bat")
    print("Outputs in:", out_py)
    print("=" * 64)


if __name__ == "__main__":
    main()
