# -*- coding: utf-8 -*-
"""Compare Case 2 Python thermoelasticity vs a COMSOL (x,y,z,T,vm) CSV."""

from __future__ import print_function

import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from atsim.geometry import load_case
from atsim.fvm3d import solve_fvm3d
from atsim.elasticity import solve_elasticity


def _open_text(path):
    try:
        return open(path, "r", encoding="utf-8")
    except TypeError:
        return open(path, "r")
    except UnicodeDecodeError:
        return open(path, "r", encoding="gbk")


def load_stress_csv(path):
    xs, ys, zs, Ts, vms = [], [], [], [], []
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
    x, y, z, T, vm = map(np.array, (xs, ys, zs, Ts, vms))
    if np.nanmax(np.abs(x)) > 0.1:
        x = x * 1e-3
        y = y * 1e-3
        z = z * 1e-3
    if np.nanmax(T) < 200:
        T = T + 273.15
    if np.nanmax(vm) > 1e4:
        vm = vm * 1e-6
    return {"x": x, "y": y, "z": z, "T": T, "vm_MPa": vm}


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


def main():
    csv_path = os.path.join(ROOT, "results", "python_case2_stress", "comsol_case2_stress.csv")
    if len(sys.argv) >= 2:
        csv_path = sys.argv[1]
        if not os.path.isabs(csv_path):
            csv_path = os.path.join(ROOT, csv_path)
    if not os.path.isfile(csv_path):
        print("[error] missing", csv_path)
        sys.exit(1)

    co = load_stress_csv(csv_path)
    print("COMSOL n=%d  T=%.2f..%.2f C  vm=%.1f..%.1f MPa" % (
        co["T"].size, co["T"].min() - 273.15, co["T"].max() - 273.15,
        co["vm_MPa"].min(), co["vm_MPa"].max()))

    case = load_case(os.path.join(ROOT, "cases", "case2_package.json"))
    th = solve_fvm3d(case, verbose=True)
    el45 = solve_elasticity(th, T_ref=case.T_amb, verbose=True)
    el20 = solve_elasticity(th, T_ref=293.15, verbose=False)

    grid = th["grid"]
    pyT = sample_xy(th["T"], grid, co["x"], co["y"], co["z"])
    py45 = sample_xy(el45["vm"], grid, co["x"], co["y"], co["z"]) * 1e-6
    py20 = sample_xy(el20["vm"], grid, co["x"], co["y"], co["z"]) * 1e-6
    print("T MAE=%.3f K" % np.abs(pyT - co["T"]).mean())
    d45 = py45 - co["vm_MPa"]
    d20 = py20 - co["vm_MPa"]
    print("vm Tref=45C  mae=%.2f  rmse=%.2f  max|d|=%.2f  mean_py=%.2f mean_co=%.2f" % (
        np.abs(d45).mean(), np.sqrt((d45 ** 2).mean()), np.abs(d45).max(),
        py45.mean(), co["vm_MPa"].mean()))
    print("vm Tref=20C  mae=%.2f  (should be worse if COMSOL Tref is 45C)" % np.abs(d20).mean())

    probes = [
        ("Chip_center", 0.030, 0.030),
        ("Cache_center", 0.028, 0.028),
        ("Uncore_center", 0.032, 0.028),
        ("CoreA_center", 0.028, 0.032),
        ("CoreB_center", 0.032, 0.032),
        ("die_SW", 0.026, 0.026),
        ("die_NW", 0.026, 0.034),
        ("die_SE", 0.034, 0.026),
        ("die_NE", 0.034, 0.034),
    ]
    print("\n%-14s %8s %8s %8s %8s" % ("pt", "Tco", "vm_co", "py45", "dTref45"))
    for name, x, y in probes:
        dist = (co["x"] - x) ** 2 + (co["y"] - y) ** 2
        i = int(dist.argmin())
        print(" %-14s %8.2f %8.2f %8.2f %8.2f" % (
            name, co["T"][i] - 273.15, co["vm_MPa"][i], py45[i], py45[i] - co["vm_MPa"][i]))


if __name__ == "__main__":
    main()
