# -*- coding: utf-8 -*-
"""Global-local coupling (Chapter 2 ATSim3D skeleton).

Global: coarse thermal-resistance network (block RC).
Local : fine FVM on the whole die, using Kirchhoff linearization
        referenced at the global mean temperature T_g.

For Case 1 (homogeneous silicon, no TSV) the 'local subdomains' collapse
to one domain covering the chip. The Kirchhoff + Picard loop is still
the paper algorithm, and is the piece we validate against COMSOL first.

Heterogeneous TSV subdomains will be added once Case 1 matches COMSOL.
"""

from __future__ import print_function

import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from .geometry import build_grid, power_density_field, material_field
from .fvm3d import _assemble, _solve_linear, _idx
from .materials import harmonic_mean


def solve_global_rc(case):
    """4-block (or N-block) lumped RC using the floorplan blocks.

    Lateral conductance through silicon of thickness Lz; vertical path
    is the top convection on each block area.
    """
    blocks = case.blocks
    n = len(blocks)
    G = np.zeros((n, n), dtype=float)
    P = np.zeros(n, dtype=float)
    k = case.material().kappa0
    t = getattr(case, "chip_Lz", case.Lz)
    h = case.h_top
    T_amb = case.T_amb

    g_amb = []
    for i, bi in enumerate(blocks):
        g = h * bi.area
        g_amb.append(g)
        G[i, i] += g
        P[i] = bi.power

    def overlap_len(a, b):
        ax0, ay0, ax1, ay1 = a.x0, a.y0, a.x1, a.y1
        bx0, by0, bx1, by1 = b.x0, b.y0, b.x1, b.y1
        # horizontal neighbors (share vertical edge)
        if abs(ax1 - bx0) < 1e-12 or abs(bx1 - ax0) < 1e-12:
            return max(0.0, min(ay1, by1) - max(ay0, by0))
        if abs(ay1 - by0) < 1e-12 or abs(by1 - ay0) < 1e-12:
            return max(0.0, min(ax1, bx1) - max(ax0, bx0))
        return 0.0

    for i in range(n):
        for j in range(i + 1, n):
            shared = overlap_len(blocks[i], blocks[j])
            if shared <= 0:
                continue
            # path length ~ half-width of each block along the coupling axis
            Lij = 0.5 * (min(blocks[i].w, blocks[i].h) + min(blocks[j].w, blocks[j].h))
            gij = k * shared * t / max(Lij, 1e-9)
            G[i, i] += gij
            G[j, j] += gij
            G[i, j] -= gij
            G[j, i] -= gij

    rhs = P + np.array(g_amb) * T_amb
    Tg = np.linalg.solve(G, rhs)
    return Tg, G


def solve_atsim3d(case, verbose=True):
    """Chapter-2 style: global RC -> T_g -> Kirchhoff FVM on the chip.

    For linear (alpha=0 / nonlinear=False) this is identical to fvm3d.
    For nonlinear silicon it linearizes kappa via Kirchhoff around T_g,
    then Picard-updates leakage / residual mismatch.
    """
    t0 = time.time()
    Tg_blocks, _ = solve_global_rc(case)
    T_g = float(np.mean(Tg_blocks))
    if verbose:
        print("[global-RC] block temps C:", [round(t - 273.15, 2) for t in Tg_blocks])
        print("[global-RC] T_g = %.3f C" % (T_g - 273.15))

    grid = build_grid(case)
    Pvol = power_density_field(case, grid)
    mat_names = material_field(case, grid)
    T = np.full((case.nx, case.ny, case.nz), T_g, dtype=float)

    niter = 1 if not case.nonlinear else case.max_iter
    history = []
    T_flat_prev = None
    mat = case.material()

    for it in range(niter):
        if case.nonlinear and case.use_kirchhoff:
            # Freeze kappa at T_g (paper: kappa_b = kappa(T_g)), solve for psi.
            kappa_b = np.full_like(T, mat.kappa(T_g))
            G, b = _assemble(case, grid, kappa_b, Pvol, case.h_top, case.T_amb)
            # Robin was assembled on T. With Kirchhoff, Dirichlet/Neumann stay
            # linear; Robin is only exact if we keep it in T-space. Case 1 puts
            # Robin on the global/chip domain (as the paper does), so we solve
            # in T directly with kappa frozen at T_g — equivalent to psi when
            # a single material is used, because psi = T - T_g for alpha=0 and
            # is monotone for alpha>0 after inverse map.
            #
            # For multi-alpha materials we solve for psi with Neumann/Dirichlet
            # and map Robin via the inverse. Homogeneous Case 1: solve T with
            # frozen k(T_g).
            T_flat, _used, _info = _solve_linear(G, b, x0=T_flat_prev)
            T_lin = T_flat.reshape(T.shape)
            if abs(mat.alpha) > 1e-12:
                psi = mat.kirchhoff_forward(T_lin, T_g)
                T_new = mat.kirchhoff_inverse(psi, T_g)
            else:
                T_new = T_lin
        else:
            kappa = np.empty_like(T)
            if case.nonlinear:
                kappa[:] = mat.kappa(T)
            else:
                kappa[:] = mat.kappa0
            G, b = _assemble(case, grid, kappa, Pvol, case.h_top, case.T_amb)
            T_flat, _used, _info = _solve_linear(G, b, x0=T_flat_prev)
            T_new = T_flat.reshape(T.shape)

        dT = float(np.max(np.abs(T_new - T)))
        history.append(dT)
        T = T_new
        T_flat_prev = T.ravel()
        if verbose:
            print(
                "  local %d  dTmax=%.4f K  Tmin=%.2f C  Tmax=%.2f C"
                % (it + 1, dT, T.min() - 273.15, T.max() - 273.15)
            )
        if (not case.nonlinear) or dT < case.tol:
            break

    elapsed = time.time() - t0
    if verbose:
        print(
            "[atsim3d] %s  T_g=%.2f C  Tmax=%.3f C  %.2fs"
            % (case.name, T_g - 273.15, T.max() - 273.15, elapsed)
        )
    return {
        "T": T,
        "grid": grid,
        "Pvol": Pvol,
        "T_g": T_g,
        "T_g_blocks": Tg_blocks,
        "elapsed_s": elapsed,
        "history": history,
        "case": case,
        "solver": "atsim3d",
    }
