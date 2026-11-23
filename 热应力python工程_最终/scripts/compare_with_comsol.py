# -*- coding: utf-8 -*-
"""Compare Python mid-plane vs a COMSOL CSV export.

Usage:
    python scripts/compare_with_comsol.py results/comsol/case1_points.csv
"""

from __future__ import print_function

import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from atsim.compare import load_comsol_csv, compare_midplanes
from atsim.geometry import load_case
from atsim.fvm3d import solve_fvm3d


def main():
    csv_path = None
    case_path = os.path.join(ROOT, "cases", "case1_2d_ic.json")
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--case" and i + 1 < len(args):
            case_path = args[i + 1]
            if not os.path.isabs(case_path):
                case_path = os.path.join(ROOT, case_path)
            i += 2
            continue
        if csv_path is None:
            csv_path = args[i]
        i += 1
    if csv_path is None:
        csv_path = os.path.join(ROOT, "results", "comsol", "case1_points.csv")
        print("Usage: python scripts/compare_with_comsol.py <comsol.csv> [--case cases/xxx.json]")
        print("Default csv:", csv_path)
    if not os.path.isfile(csv_path):
        print("[error] COMSOL csv not found:", csv_path)
        sys.exit(1)

    case = load_case(case_path)
    print("Case:", case.name, " from", case_path)
    print("Re-solving Python FVM ...")
    res = solve_fvm3d(case, verbose=True)
    comsol = load_comsol_csv(csv_path)
    print("COMSOL points:", comsol["T"].size, " Tmax=%.3f C" % (comsol["T"].max() - 273.15))
    out = os.path.join(ROOT, "results", "compare_" + case.name)
    stats = compare_midplanes(res, comsol, out)
    print("\nComparison (Python sampled at COMSOL nodes):")
    for k, v in stats.items():
        if k == "probes":
            continue
        print("  {0:22s} {1}".format(k, v))
    if "probes" in stats:
        print("\nProbes:")
        print("  {0:16s} {1:>10s} {2:>10s} {3:>8s}".format("name", "Python C", "COMSOL C", "dT K"))
        for r in stats["probes"]:
            print("  {0:16s} {1:10.3f} {2:10.3f} {3:8.3f}".format(
                r["name"], r["python_C"], r["comsol_C"], r["dT_K"]))
    print("Wrote", os.path.join(out, "compare_stats.txt"))
    print("Plot  ", os.path.join(out, "compare_midplane.png"))


if __name__ == "__main__":
    main()
