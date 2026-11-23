# -*- coding: utf-8 -*-
"""Run Case 1: Python FVM + ATSim3D + HotSpot file dump."""

from __future__ import print_function

import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from atsim.geometry import load_case
from atsim.fvm3d import solve_fvm3d
from atsim.global_local import solve_atsim3d, solve_global_rc
from atsim.hotspot_io import write_hotspot_case
from atsim.visualize import (
    plot_midplane,
    save_midplane_csv,
    save_numpy,
    save_probes_csv,
    sample_probes,
)


def analytical_1d(case):
    """Uniform-power 1D slab: convection on top, adiabatic bottom.

    T(z) = T_amb + P/(h A) + P (Lz - z) / (k A)   wait: generation throughout.

    With uniform volumetric generation q = P / (A Lz),
      d/dz (k dT/dz) + q = 0
      z=0: dT/dz=0
      z=Lz: -k dT/dz = h (T - T_amb)
    => T(z) = T_amb + q Lz / h + q (Lz^2 - z^2) / (2k)
       Tmax at z=0: T_amb + P/(hA) + P Lz / (2 k A)
    """
    P = case.total_power
    A = case.Lx * case.Ly
    k = case.material().kappa0
    h = case.h_top
    Tmax = case.T_amb + P / (h * A) + P * case.Lz / (2.0 * k * A)
    Ttop = case.T_amb + P / (h * A)
    return Tmax, Ttop


def main():
    case_path = os.path.join(ROOT, "cases", "case1_2d_ic.json")
    out_py = os.path.join(ROOT, "results", "python")
    out_hs = os.path.join(ROOT, "hotspot_case1")
    if not os.path.isdir(out_py):
        os.makedirs(out_py)

    case = load_case(case_path)
    print("=" * 64)
    print("Case:", case.name)
    print("Die: %.1f mm x %.1f mm x %.0f um" % (case.Lx * 1e3, case.Ly * 1e3, case.Lz * 1e6))
    print("Grid: %d x %d x %d" % (case.nx, case.ny, case.nz))
    print("Total power: %.2f W   h=%.0f W/m2K   Tamb=%.2f C" % (
        case.total_power, case.h_top, case.T_amb - 273.15))
    for b in case.blocks:
        print("  %-8s  %.1f W  density=%.2f W/mm2" % (
            b.name, b.power, b.power / (b.area * 1e6)))

    Tmax_1d, Ttop_1d = analytical_1d(case)
    print("1D uniform analytic: Tmax=%.2f C  Ttop=%.2f C" % (
        Tmax_1d - 273.15, Ttop_1d - 273.15))

    print("\n-- global RC --")
    Tg, _ = solve_global_rc(case)
    for b, t in zip(case.blocks, Tg):
        print("  %-8s  %.2f C" % (b.name, t - 273.15))

    print("\n-- FVM 3D (linear) --")
    res = solve_fvm3d(case, verbose=True)

    print("\n-- ATSim3D (global-local, linear) --")
    res2 = solve_atsim3d(case, verbose=True)
    dmax = float((res["T"] - res2["T"]).max())
    dmin = float((res["T"] - res2["T"]).min())
    print("FVM vs ATSim3D  max|dT|=%.4f K" % max(abs(dmax), abs(dmin)))

    probes = save_probes_csv(res, os.path.join(out_py, "probes.csv"))
    save_midplane_csv(res, os.path.join(out_py, "midplane.csv"))
    save_numpy(res, os.path.join(out_py, "T_field.npz"))
    plot_midplane(res, os.path.join(out_py, "midplane.png"))

    print("\nProbes:")
    for r in probes:
        print("  %-16s  %.3f C" % (r["name"], r["T_C"]))

    summary = {
        "case": case.name,
        "solver": "fvm3d",
        "grid": [case.nx, case.ny, case.nz],
        "total_power_W": case.total_power,
        "h_W_m2K": case.h_top,
        "T_amb_C": case.T_amb - 273.15,
        "Tmax_C": float(res["T"].max() - 273.15),
        "Tmin_C": float(res["T"].min() - 273.15),
        "Tmean_C": float(res["T"].mean() - 273.15),
        "Tmax_1d_uniform_C": Tmax_1d - 273.15,
        "Ttop_1d_uniform_C": Ttop_1d - 273.15,
        "elapsed_s": res["elapsed_s"],
        "probes": probes,
        "block_global_RC_C": {
            b.name: float(t - 273.15) for b, t in zip(case.blocks, Tg)
        },
    }
    with open(os.path.join(out_py, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    hs_files = write_hotspot_case(case, out_hs)
    print("\nHotSpot files:")
    for k, v in hs_files.items():
        print("  %s: %s" % (k, v))
    print("Run:  hotspot_case1/run_hotspot.bat   (after compiling HotSpot)")
    print("Outputs in:", out_py)
    print("=" * 64)


if __name__ == "__main__":
    main()
