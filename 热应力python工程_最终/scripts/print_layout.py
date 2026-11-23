# -*- coding: utf-8 -*-
"""Print Case 3 chips[] expansion: absolute mm blocks for COMSOL.

    python scripts/print_layout.py
    python scripts/print_layout.py cases/case3_4chip.json
"""

from __future__ import print_function

import csv
import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from atsim.geometry import load_case, build_grid


def main():
    case_path = os.path.join(ROOT, "cases", "case3_chips.json")
    if len(sys.argv) >= 2:
        case_path = sys.argv[1]
        if not os.path.isabs(case_path):
            case_path = os.path.join(ROOT, case_path)
    case = load_case(case_path)
    grid = build_grid(case)
    meta = case.raw.get("_layout_meta") or {}
    out_dir = os.path.join(ROOT, "results", "layout_" + os.path.splitext(os.path.basename(case_path))[0])
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    print("=" * 64)
    print("JSON:", case_path)
    print("name:", case.name)
    print("Domain %.1f x %.1f x %.2f mm   P=%.2f W   h=%.2f W/m2K" % (
        case.Lx * 1e3, case.Ly * 1e3, case.Lz * 1e3, case.total_power, case.h_top))
    print("BBox %.1f x %.1f mm  Interposer %.1f x %.1f  Spreader %.1f x %.1f  Sink %.1f x %.1f" % (
        meta.get("bbox_mm", [0, 0])[0], meta.get("bbox_mm", [0, 0])[1],
        (meta.get("interposer_mm") or [0, 0])[0], (meta.get("interposer_mm") or [0, 0])[1],
        (meta.get("spreader_mm") or [0, 0])[0], (meta.get("spreader_mm") or [0, 0])[1],
        (meta.get("sink_mm") or [0, 0])[0], (meta.get("sink_mm") or [0, 0])[1]))
    n_solid = int(grid["solid"].sum()) if grid.get("solid") is not None else case.nx * case.ny * case.nz
    print("grid %d x %d x %d  solid=%d  shift (%.2f, %.2f) mm   tim_cover=%s" % (
        case.nx, case.ny, case.nz, n_solid,
        (meta.get("shift_mm") or [0, 0])[0], (meta.get("shift_mm") or [0, 0])[1],
        meta.get("tim_cover")))
    print("\nStack (COMSOL blocks, mm, origin = sink bottom-left, z from substrate bottom):")
    print("  %-16s %8s %8s %8s %8s %8s %8s %8s" % (
        "name", "x0", "y0", "z0", "W", "D", "H", "mat"))
    rows = []
    for L in case.stack:
        row = {
            "name": L.name,
            "material": L.material,
            "x0_mm": L.x0 * 1e3,
            "y0_mm": L.y0 * 1e3,
            "z0_mm": L.z0 * 1e3,
            "w_mm": L.w * 1e3,
            "d_mm": L.h * 1e3,
            "h_mm": (L.z1 - L.z0) * 1e3,
            "power_W": 0.0,
        }
        rows.append(row)
        print("  %-16s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8s" % (
            L.name, row["x0_mm"], row["y0_mm"], row["z0_mm"],
            row["w_mm"], row["d_mm"], row["h_mm"], L.material))
    print("\nChips (power):")
    print("  %-10s %8s %8s %8s %8s %8s" % ("name", "x0", "y0", "w", "h", "W"))
    for b in case.blocks:
        print("  %-10s %8.2f %8.2f %8.2f %8.2f %8.2f" % (
            b.name, b.x0 * 1e3, b.y0 * 1e3, b.w * 1e3, b.h * 1e3, b.power))
        for r in rows:
            if r["name"] == "die_" + b.name:
                r["power_W"] = b.power
    print("\nProbes (mm):")
    for p in case.probes:
        print("  %-20s  %7.3f  %7.3f  %7.3f" % (p.name, p.x * 1e3, p.y * 1e3, p.z * 1e3))

    with open(os.path.join(out_dir, "comsol_blocks.csv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(out_dir, "layout_meta.json"), "w") as f:
        json.dump({
            "json": case_path,
            "name": case.name,
            "grid": [case.nx, case.ny, case.nz],
            "layout_meta": meta,
            "probes_mm": [
                {"name": p.name, "x_mm": p.x * 1e3, "y_mm": p.y * 1e3, "z_mm": p.z * 1e3}
                for p in case.probes
            ],
        }, f, indent=2)
    print("\nWrote", out_dir)
    print("=" * 64)


if __name__ == "__main__":
    main()
