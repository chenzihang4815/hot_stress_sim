# -*- coding: utf-8 -*-
"""Run Case 3 from a chips[] layout JSON. Solver is the same FVM as Case 2."""

from __future__ import print_function

import json
import os
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from atsim.geometry import load_case
from atsim.fvm3d import solve_fvm3d
from atsim.hotspot_io import write_hotspot_case
from atsim.visualize import (
    die_mid_k,
    plot_midplane,
    save_midplane_csv,
    save_numpy,
    save_probes_csv,
)


def main():
    case_path = os.path.join(ROOT, "cases", "case3_chips.json")
    if len(sys.argv) >= 2:
        case_path = sys.argv[1]
        if not os.path.isabs(case_path):
            case_path = os.path.join(ROOT, case_path)
    name = os.path.splitext(os.path.basename(case_path))[0]
    out_py = os.path.join(ROOT, "results", "python_" + name)
    out_hs = os.path.join(ROOT, "hotspot_" + name)
    if not os.path.isdir(out_py):
        os.makedirs(out_py)

    case = load_case(case_path)
    meta = case.raw.get("_layout_meta") or {}
    print("=" * 64)
    print("Case:", case.name)
    print("JSON:", case_path)
    print("Domain: %.1f x %.1f x %.2f mm" % (case.Lx * 1e3, case.Ly * 1e3, case.Lz * 1e3))
    print("Chips: %d   P=%.2f W   h=%.2f W/m2K   Tamb=%.2f C" % (
        meta.get("n_chips", len(case.blocks)),
        case.total_power, case.h_top, case.T_amb - 273.15))
    sp = meta.get("spreader_mm", 0)
    sk = meta.get("sink_mm", 0)
    if isinstance(sp, (list, tuple)):
        print("BBox %.1f x %.1f mm   Spreader %.1f x %.1f mm   Sink %.1f x %.1f mm" % (
            meta.get("bbox_mm", [0, 0])[0], meta.get("bbox_mm", [0, 0])[1],
            sp[0], sp[1], sk[0], sk[1]))
        print("Centered, layout shift (%.2f, %.2f) mm" % tuple(meta.get("shift_mm", [0, 0])))
        print("Substrate %.1f mm  Interposer %.0f um  %.1f x %.1f mm  fill_die=%s  tim_cover=%s" % (
            meta.get("t_substrate_mm", 0), meta.get("t_interposer_um", 0),
            (meta.get("interposer_mm") or [0, 0])[0],
            (meta.get("interposer_mm") or [0, 0])[1],
            meta.get("fill_die"), meta.get("tim_cover")))
    else:
        print("Spreader %.1f mm   Sink %.1f mm" % (sp, sk))
    print("\nChips (edit cases/case3_chips.json):")
    print("  %-10s %8s %8s %8s %8s %8s" % ("name", "x_mm", "y_mm", "w_mm", "h_mm", "W"))
    for b in case.blocks:
        print("  %-10s %8.2f %8.2f %8.2f %8.2f %8.2f" % (
            b.name, b.x0 * 1e3, b.y0 * 1e3, b.w * 1e3, b.h * 1e3, b.power))

    print("\n-- FVM 3D package --")
    res = solve_fvm3d(case, verbose=True)
    solid = res["grid"]["solid"]
    Tsol = res["T"][solid]
    kmid = die_mid_k(res)

    probes = save_probes_csv(res, os.path.join(out_py, "probes.csv"))
    save_midplane_csv(res, os.path.join(out_py, "die_midplane.csv"), k=kmid)
    save_numpy(res, os.path.join(out_py, "T_field.npz"))
    plot_midplane(
        res,
        os.path.join(out_py, "die_midplane.png"),
        title="%s  Tmax=%.2f C" % (case.name, Tsol.max() - 273.15),
    )

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
        "T_amb_C": case.T_amb - 273.15,
        "Tmax_solid_C": float(Tsol.max() - 273.15),
        "Tmin_solid_C": float(Tsol.min() - 273.15),
        "elapsed_s": res["elapsed_s"],
        "timing": res.get("timing"),
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
        "probes": probes,
        "layout_meta": meta,
    }
    with open(os.path.join(out_py, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    timing = res.get("timing") or {"total_s": res["elapsed_s"]}
    with open(os.path.join(out_py, "timing.json"), "w") as f:
        json.dump({
            "case": case.name,
            "json": case_path,
            "wall_s": res["elapsed_s"],
            "heat": timing,
        }, f, indent=2)
    print("\nTiming (heat only):")
    print("  grid     %.2f s" % timing.get("grid_s", 0.0))
    print("  assemble %.2f s" % timing.get("assemble_s", 0.0))
    print("  solve    %.2f s  (%s)" % (timing.get("solve_s", 0.0), timing.get("solver", "?")))
    print("  total    %.2f s   ndof=%s  solid=%s" % (
        timing.get("total_s", res["elapsed_s"]),
        timing.get("n_unknown"), timing.get("n_solid")))

    hs_files = write_hotspot_case(case, out_hs)
    print("\nHotSpot files:")
    for k, v in hs_files.items():
        print("  %s: %s" % (k, v))
    print("Outputs:", out_py)
    print("=" * 64)


if __name__ == "__main__":
    main()
