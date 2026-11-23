# -*- coding: utf-8 -*-
"""3D finite-volume steady heat solver on a Cartesian grid.

Discretizes

    div( kappa(r,T) grad T ) + P(r,T) = 0

with:
  - top face Robin convection: -kappa dT/dz = h (T - T_amb)
  - other faces adiabatic
  - harmonic-mean face conductivity
  - Picard iteration for kappa(T) and optional leakage

For linear (constant-k) problems this is a single sparse solve.
"""

from __future__ import print_function

import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve, cg, LinearOperator, spilu

from .geometry import build_grid, power_density_field, material_field
from .materials import harmonic_mean


def _idx(i, j, k, ny, nz):
    return (i * ny + j) * nz + k


def _as1d(v, n):
    a = np.asarray(v, dtype=float)
    if a.ndim == 0:
        return np.full(n, float(a))
    return a


def _assemble(case, grid, kappa, Pvol, h_top, T_amb):
    """Assemble G T = b for cell-centered FVM.

    Supports uniform or non-uniform spacing. Void cells (solid==False)
    are pinned to T_amb and do not exchange heat. A solid/void face is
    treated as adiabatic. Robin convection is applied on the domain top
    face of every solid cell (Case 1: die top; Case 2: sink top).
    """
    nx, ny, nz = int(grid["nx"]), int(grid["ny"]), int(grid["nz"])
    dx = _as1d(grid["dx"], nx)
    dy = _as1d(grid["dy"], ny)
    dz = _as1d(grid["dz"], nz)
    solid = grid.get("solid")
    if solid is None:
        solid = np.ones((nx, ny, nz), dtype=bool)
    n = nx * ny * nz

    rows = []
    cols = []
    data = []
    b = np.zeros(n, dtype=float)

    def add(r, c, v):
        rows.append(r)
        cols.append(c)
        data.append(v)

    for i in range(nx):
        dxi = dx[i]
        for j in range(ny):
            dyj = dy[j]
            for k in range(nz):
                p = _idx(i, j, k, ny, nz)
                if not solid[i, j, k]:
                    add(p, p, 1.0)
                    b[p] = T_amb
                    continue
                dzk = dz[k]
                vol = dxi * dyj * dzk
                diag = 0.0
                rhs = Pvol[i, j, k] * vol
                kp = kappa[i, j, k]
                Ax = dyj * dzk
                Ay = dxi * dzk
                Az = dxi * dyj

                def g_face(k_n, d_p, d_n, area):
                    # center-to-center: R = d_p/(2 k_p A) + d_n/(2 k_n A)
                    r = 0.0
                    if kp > 0:
                        r += d_p / (2.0 * kp)
                    if k_n > 0:
                        r += d_n / (2.0 * k_n)
                    if r <= 0:
                        return 0.0
                    return area / r

                if i > 0 and solid[i - 1, j, k]:
                    g = g_face(kappa[i - 1, j, k], dxi, dx[i - 1], Ax)
                    add(p, _idx(i - 1, j, k, ny, nz), -g)
                    diag += g
                if i < nx - 1 and solid[i + 1, j, k]:
                    g = g_face(kappa[i + 1, j, k], dxi, dx[i + 1], Ax)
                    add(p, _idx(i + 1, j, k, ny, nz), -g)
                    diag += g
                if j > 0 and solid[i, j - 1, k]:
                    g = g_face(kappa[i, j - 1, k], dyj, dy[j - 1], Ay)
                    add(p, _idx(i, j - 1, k, ny, nz), -g)
                    diag += g
                if j < ny - 1 and solid[i, j + 1, k]:
                    g = g_face(kappa[i, j + 1, k], dyj, dy[j + 1], Ay)
                    add(p, _idx(i, j + 1, k, ny, nz), -g)
                    diag += g
                if k > 0 and solid[i, j, k - 1]:
                    g = g_face(kappa[i, j, k - 1], dzk, dz[k - 1], Az)
                    add(p, _idx(i, j, k - 1, ny, nz), -g)
                    diag += g
                if k < nz - 1 and solid[i, j, k + 1]:
                    g = g_face(kappa[i, j, k + 1], dzk, dz[k + 1], Az)
                    add(p, _idx(i, j, k + 1, ny, nz), -g)
                    diag += g
                elif k == nz - 1:
                    r_cond = dzk / (2.0 * kp) if kp > 0 else 0.0
                    r_conv = 1.0 / h_top if h_top > 0 else 1e30
                    g_top = Az / (r_cond + r_conv)
                    diag += g_top
                    rhs += g_top * T_amb

                if diag == 0.0:
                    diag = 1.0
                    rhs = T_amb
                add(p, p, diag)
                b[p] = rhs

    G = coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    return G, b


def _solve_linear(G, b, x0=None):
    n = b.size
    used = "spsolve"
    info = 0
    # Package FVM is strongly anisotropic (dz << dx). ILU-CG hits
    # maxiter then we spsolve anyway, so try the direct factor first.
    try:
        x = spsolve(G, b)
        return x, used, info
    except Exception as e:
        print("[warn] spsolve failed (%s), CG" % e)
    ilu = spilu(G.tocsc(), drop_tol=1e-4, fill_factor=10)
    M = LinearOperator((n, n), matvec=ilu.solve)
    x, info = cg(G, b, x0=x0, M=M, tol=1e-8, maxiter=2000)
    used = "cg_after_spsolve"
    if info != 0:
        print("[warn] CG info=%s after spsolve fail" % info)
    return x, used, info


def _kappa_from_T(case, mat_names, T):
    kappa = np.zeros(T.shape, dtype=float)
    unique = set(mat_names.ravel().tolist())
    for name in unique:
        if name is None:
            continue
        mat = case.material(name)
        mask = mat_names == name
        if case.nonlinear:
            kappa[mask] = mat.kappa(T[mask])
        else:
            kappa[mask] = mat.kappa0
    return kappa


def solve_fvm3d(case, verbose=True):
    """Return a result dict with T (nx,ny,nz) in Kelvin."""
    t0 = time.time()
    grid = build_grid(case)
    t_grid = time.time()
    Pvol = power_density_field(case, grid)
    mat_names = material_field(case, grid)
    T = np.full((case.nx, case.ny, case.nz), case.T_amb, dtype=float)
    niter = 1 if not case.nonlinear else case.max_iter
    history = []
    T_flat_prev = None
    assemble_s = 0.0
    solve_s = 0.0
    heat_solver = None
    G = None

    for it in range(niter):
        kappa = _kappa_from_T(case, mat_names, T)
        ta = time.time()
        G, b = _assemble(case, grid, kappa, Pvol, case.h_top, case.T_amb)
        assemble_s += time.time() - ta
        ts = time.time()
        T_flat, heat_solver, _info = _solve_linear(G, b, x0=T_flat_prev)
        solve_s += time.time() - ts
        T_new = T_flat.reshape((case.nx, case.ny, case.nz))
        dT = float(np.max(np.abs(T_new - T)))
        history.append(dT)
        T = T_new
        T_flat_prev = T_flat
        solid = grid.get("solid")
        if solid is not None:
            Ts = T[solid]
        else:
            Ts = T
        if verbose:
            print(
                "  Picard %d  dTmax=%.4f K  Tmin=%.2f C  Tmax=%.2f C"
                % (it + 1, dT, Ts.min() - 273.15, Ts.max() - 273.15)
            )
        if (not case.nonlinear) or dT < case.tol:
            break

    elapsed = time.time() - t0
    solid = grid.get("solid")
    Ts = T[solid] if solid is not None else T
    n_unknown = int(case.nx * case.ny * case.nz)
    timing = {
        "grid_s": t_grid - t0,
        "assemble_s": assemble_s,
        "solve_s": solve_s,
        "total_s": elapsed,
        "n_iter": len(history),
        "n_cells": n_unknown,
        "n_solid": int(solid.sum()) if solid is not None else n_unknown,
        "n_unknown": n_unknown,
        "nnz": int(G.nnz) if G is not None else 0,
        "solver": heat_solver or "spsolve",
        "grid": [int(case.nx), int(case.ny), int(case.nz)],
    }
    if verbose:
        print(
            "[fvm3d] %s  cells=%dx%dx%d  solid=%d  Tmax=%.3f C  Tmin=%.3f C  %.2fs"
            % (
                case.name,
                case.nx,
                case.ny,
                case.nz,
                int(Ts.size),
                Ts.max() - 273.15,
                Ts.min() - 273.15,
                elapsed,
            )
        )
        print(
            "[fvm3d] timing  grid=%.2fs  assemble=%.2fs  solve=%.2fs (%s)  ndof=%d  nnz=%d"
            % (
                timing["grid_s"],
                timing["assemble_s"],
                timing["solve_s"],
                timing["solver"],
                timing["n_unknown"],
                timing["nnz"],
            )
        )
    return {
        "T": T,
        "grid": grid,
        "Pvol": Pvol,
        "kappa": _kappa_from_T(case, mat_names, T),
        "mat_names": mat_names,
        "elapsed_s": elapsed,
        "history": history,
        "case": case,
        "solver": "fvm3d",
        "timing": timing,
    }
