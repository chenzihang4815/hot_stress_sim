# -*- coding: utf-8 -*-
"""Sanity checks for Kirchhoff transform invertibility and 1D FVM."""

from __future__ import print_function

import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from atsim.materials import Material
from atsim.geometry import load_case
from atsim.fvm3d import solve_fvm3d


def test_roundtrip():
    si = Material("silicon", kappa0=130.0, T_ref=300.0, alpha=1.5)
    cu = Material("copper", kappa0=400.0, T_ref=300.0, alpha=0.0)
    T_g = 350.0
    T = np.array([300.0, 320.0, 350.0, 400.0, 500.0])
    for mat in (si, cu):
        psi = mat.kirchhoff_forward(T, T_g)
        T2 = mat.kirchhoff_inverse(psi, T_g)
        err = np.max(np.abs(T2 - T))
        print("%s roundtrip max|dT|=%.3e K  psi=%s" % (mat.name, err, np.round(psi, 3)))
        assert err < 1e-6, err


def test_uniform_1d():
    """If power is uniform, 3D FVM must recover the 1D analytic slab."""
    case = load_case(os.path.join(ROOT, "cases", "case1_2d_ic.json"))
    # overwrite blocks to uniform 17.5 W total still, but equal density:
    total = case.total_power
    for b in case.blocks:
        b.power = total * b.area / (case.Lx * case.Ly)
    case.nx, case.ny, case.nz = 16, 16, 8
    case.nonlinear = False
    res = solve_fvm3d(case, verbose=False)
    P = total
    A = case.Lx * case.Ly
    k = case.material().kappa0
    h = case.h_top
    Tmax_an = case.T_amb + P / (h * A) + P * case.Lz / (2.0 * k * A)
    Ttop_an = case.T_amb + P / (h * A)
    Tmax_num = res["T"].max()
    Ttop_num = res["T"][:, :, -1].mean()
    print("1D analytic Tmax=%.4f C  Ttop=%.4f C" % (Tmax_an - 273.15, Ttop_an - 273.15))
    print("FVM        Tmax=%.4f C  Ttop_mean=%.4f C" % (Tmax_num - 273.15, Ttop_num - 273.15))
    print("error Tmax=%.4f K  Ttop=%.4f K" % (Tmax_num - Tmax_an, Ttop_num - Ttop_an))
    assert abs(Tmax_num - Tmax_an) < 0.15
    assert abs(Ttop_num - Ttop_an) < 0.15


if __name__ == "__main__":
    test_roundtrip()
    test_uniform_1d()
    print("ALL TESTS PASSED")
