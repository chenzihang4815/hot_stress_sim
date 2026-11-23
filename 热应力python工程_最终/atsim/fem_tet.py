# -*- coding: utf-8 -*-
"""Quadratic tet10 heat + elasticity on the same Cartesian voxels.

Each solid FVM cell is split into 6 face-conforming tet10s (Freudenthal).
Heat and mechanics share the same nodes. T is interpolated at elasticity
Gauss points (not cell-constant). One-way: heat first, then stress.

Does not replace fvm3d.py. Call solve_tet10() when element=tet10.
"""

from __future__ import print_function

import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import cg, LinearOperator, spsolve

from .geometry import build_grid, power_density_field, material_field
from .elasticity import (
    _iso_C,
    _fill_B,
    _von_mises,
    _fix_components,
    _solve_cg,
)

_LOC = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
        (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
_TETS6 = (
    (0, 1, 2, 6),
    (0, 2, 3, 6),
    (0, 3, 7, 6),
    (0, 7, 4, 6),
    (0, 4, 5, 6),
    (0, 5, 1, 6),
)
_TET10_EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
_TET_GP = 0.5854101966249685
_TET_GN = 0.1381966011250105
_TET_GAUSS = (
    (_TET_GP, _TET_GN, _TET_GN, _TET_GN),
    (_TET_GN, _TET_GP, _TET_GN, _TET_GN),
    (_TET_GN, _TET_GN, _TET_GP, _TET_GN),
    (_TET_GN, _TET_GN, _TET_GN, _TET_GP),
)
_TET_W = 1.0 / 24.0
# tet local faces on the hex top (z=1): two triangles covering the quad
_TOP_TETS = (3, 4)
_TOP_FACE = ((1, 2, 3), (1, 2, 3))
_TRI_G = ((2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0),
          (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0),
          (1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0))


def _axis_coord(idx, edges):
    if idx % 2 == 0:
        return float(edges[idx // 2])
    return 0.5 * (float(edges[idx // 2]) + float(edges[idx // 2 + 1]))


def _tet10_nodes(i, j, k, tet_local):
    verts = []
    for a in tet_local:
        di, dj, dk = _LOC[a]
        verts.append((2 * i + 2 * di, 2 * j + 2 * dj, 2 * k + 2 * dk))
    nodes = list(verts)
    for a, b in _TET10_EDGES:
        va, vb = verts[a], verts[b]
        nodes.append((
            (va[0] + vb[0]) // 2,
            (va[1] + vb[1]) // 2,
            (va[2] + vb[2]) // 2,
        ))
    return nodes


def _tet_xyz(nodes, xe, ye, ze):
    xyz = np.empty((10, 3), dtype=float)
    for p, (I, J, K) in enumerate(nodes):
        xyz[p, 0] = _axis_coord(I, xe)
        xyz[p, 1] = _axis_coord(J, ye)
        xyz[p, 2] = _axis_coord(K, ze)
    return xyz


def _tet_N(L0, L1, L2, L3):
    return np.array([
        L0 * (2.0 * L0 - 1.0),
        L1 * (2.0 * L1 - 1.0),
        L2 * (2.0 * L2 - 1.0),
        L3 * (2.0 * L3 - 1.0),
        4.0 * L0 * L1,
        4.0 * L1 * L2,
        4.0 * L2 * L0,
        4.0 * L0 * L3,
        4.0 * L1 * L3,
        4.0 * L2 * L3,
    ], dtype=float)


def _tet_dN_ref(L0, L1, L2, L3):
    dL = np.array([
        [-1.0, -1.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ], dtype=float)
    dN = np.zeros((10, 3), dtype=float)
    Ls = (L0, L1, L2, L3)
    for a in range(4):
        dN[a, :] = (4.0 * Ls[a] - 1.0) * dL[a]
    pairs = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
    for p, (a, b) in enumerate(pairs):
        dN[4 + p, :] = 4.0 * (dL[a] * Ls[b] + Ls[a] * dL[b])
    return dN


def _tri6_N(L0, L1, L2):
    return np.array([
        L0 * (2.0 * L0 - 1.0),
        L1 * (2.0 * L1 - 1.0),
        L2 * (2.0 * L2 - 1.0),
        4.0 * L0 * L1,
        4.0 * L1 * L2,
        4.0 * L2 * L0,
    ], dtype=float)


def _ke_heat_unit(xyz, kappa=1.0):
    Ke = np.zeros((10, 10), dtype=float)
    feN = np.zeros(10, dtype=float)
    vol = 0.0
    for L0, L1, L2, L3 in _TET_GAUSS:
        dNref = _tet_dN_ref(L0, L1, L2, L3)
        J = xyz.T.dot(dNref)
        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)
        dNxyz = dNref.dot(invJ)
        wdet = abs(detJ) * _TET_W
        vol += wdet
        Ke += kappa * dNxyz.dot(dNxyz.T) * wdet
        N = _tet_N(L0, L1, L2, L3)
        feN += N * wdet
    return Ke, feN, vol


def _ke_elast(xyz, C):
    Ke = np.zeros((30, 30), dtype=float)
    B = np.zeros((6, 30), dtype=float)
    gps = []
    for L0, L1, L2, L3 in _TET_GAUSS:
        dNref = _tet_dN_ref(L0, L1, L2, L3)
        J = xyz.T.dot(dNref)
        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)
        dNxyz = dNref.dot(invJ)
        wdet = abs(detJ) * _TET_W
        B[:, :] = 0.0
        for a in range(10):
            _fill_B(B, 3 * a, dNxyz[a, 0], dNxyz[a, 1], dNxyz[a, 2])
        Ke += B.T.dot(C.dot(B)) * wdet
        N = _tet_N(L0, L1, L2, L3)
        gps.append((N.copy(), B.copy(), wdet))
    return Ke, gps


def _fine_shape(grid):
    nx, ny, nz = int(grid["nx"]), int(grid["ny"]), int(grid["nz"])
    return 2 * nx + 1, 2 * ny + 1, 2 * nz + 1


def _mark_nodes(solid):
    nx, ny, nz = solid.shape
    nnx, nny, nnz = 2 * nx + 1, 2 * ny + 1, 2 * nz + 1
    used = np.zeros((nnx, nny, nnz), dtype=bool)
    si, sj, sk = np.where(solid)
    for t in range(6):
        tet = _TETS6[t]
        for p in range(len(si)):
            nodes = _tet10_nodes(int(si[p]), int(sj[p]), int(sk[p]), tet)
            for I, J, K in nodes:
                used[I, J, K] = True
    return used


def _heat_dof(used):
    nid = -np.ones(used.shape, dtype=int)
    nd = 0
    idx = np.argwhere(used)
    for p in range(idx.shape[0]):
        i, j, k = int(idx[p, 0]), int(idx[p, 1]), int(idx[p, 2])
        nid[i, j, k] = nd
        nd += 1
    return nid, nd, idx


def _solve_heat_cg(G, b, verbose=True):
    n = b.size
    try:
        diag = G.diagonal()
        Minv = 1.0 / np.where(np.abs(diag) > 0, diag, 1.0)
        Aop = LinearOperator((n, n), matvec=lambda x: G.dot(x))
        Mop = LinearOperator((n, n), matvec=lambda x: Minv * x)
        x, info = cg(Aop, b, M=Mop, tol=1e-8, maxiter=8000)
        used = "cg"
        if info != 0:
            if verbose:
                print("[tet10 heat] CG info=%s, spsolve" % info)
            x = spsolve(G, b)
            used = "spsolve_after_cg"
        return x, used
    except Exception as e:
        if verbose:
            print("[tet10 heat] CG failed (%s), spsolve" % e)
        return spsolve(G, b), "spsolve"


def solve_tet10(case, verbose=True, T_ref=None, bc=None):
    """Return (thermal, elast) dicts sampled on Cartesian cell centers."""
    t0 = time.time()
    grid = build_grid(case)
    Pvol = power_density_field(case, grid)
    mat_names = material_field(case, grid)
    if T_ref is None:
        T_ref = case.T_amb
    mech = case.raw.get("mechanics") or {}
    if bc is None:
        bc = mech.get("bc", "substrate_bottom_fixed")
    nx, ny, nz = int(grid["nx"]), int(grid["ny"]), int(grid["nz"])
    xe, ye, ze = grid["x_edges"], grid["y_edges"], grid["z_edges"]
    solid = grid["solid"]
    si, sj, sk = np.where(solid)
    n_hex = int(si.size)
    h_top = float(case.h_top)
    T_amb = float(case.T_amb)

    used = _mark_nodes(solid)
    hid, nh, used_idx = _heat_dof(used)
    if verbose:
        print("[tet10] nodes=%d  solid hex=%d  tets=%d" % (nh, n_hex, n_hex * 6))

    ke_h_cache = {}
    n_tets = n_hex * 6
    h_nnz_est = n_tets * 100 + n_hex * 2 * 36
    hr = np.empty(h_nnz_est, dtype=np.int64)
    hc = np.empty(h_nnz_est, dtype=np.int64)
    hd = np.empty(h_nnz_est, dtype=float)
    hpos = 0
    fb = np.zeros(nh, dtype=float)

    def _add_block(gdof, Ke_s, pos, rr, cc, dd):
        free = gdof >= 0
        gv = gdof[free]
        if gv.size == 0:
            return pos
        M = Ke_s[np.ix_(free, free)]
        n = gv.size
        block = n * n
        rr[pos:pos + block] = np.repeat(gv, n)
        cc[pos:pos + block] = np.tile(gv, n)
        dd[pos:pos + block] = M.ravel()
        return pos + block

    for t in range(n_hex):
        i, j, k = int(si[t]), int(sj[t]), int(sk[t])
        name = mat_names[i, j, k]
        mat = case.material(name) if name else case.material()
        kap = float(mat.kappa0)
        dxi, dyj, dzk = float(grid["dx"][i]), float(grid["dy"][j]), float(grid["dz"][k])
        key = (round(dxi, 16), round(dyj, 16), round(dzk, 16), t % 6)
        for ti, tet in enumerate(_TETS6):
            nodes = _tet10_nodes(i, j, k, tet)
            xyz = _tet_xyz(nodes, xe, ye, ze)
            ck = key[:-1] + (ti,)
            hit = ke_h_cache.get(ck)
            if hit is None:
                Ke1, feN, vol = _ke_heat_unit(xyz, 1.0)
                ke_h_cache[ck] = (Ke1, feN, vol)
            else:
                Ke1, feN, vol = hit
            gdof = np.empty(10, dtype=np.int64)
            for a, (I, J, K) in enumerate(nodes):
                gdof[a] = hid[I, J, K]
            hpos = _add_block(gdof, Ke1 * kap, hpos, hr, hc, hd)
            P = float(Pvol[i, j, k])
            for a in range(10):
                ga = gdof[a]
                if ga >= 0:
                    fb[ga] += P * feN[a]
            if k == nz - 1 and ti in _TOP_TETS and h_top > 0:
                # quadratic triangle on hex top
                face = _TOP_FACE[0]
                c0, c1, c2 = face
                # face node order: 3 corners then 3 mids 0-1, 1-2, 2-0
                loc = [c0, c1, c2]
                # mids among tet10 nodes 4-9 corresponding to edges
                # edge (c0,c1), (c1,c2), (c2,c0)
                edge_of = {}
                for e, (ea, eb) in enumerate(_TET10_EDGES):
                    edge_of[(ea, eb)] = 4 + e
                    edge_of[(eb, ea)] = 4 + e
                loc.extend([
                    edge_of[(c0, c1)],
                    edge_of[(c1, c2)],
                    edge_of[(c2, c0)],
                ])
                x0 = xyz[c0]
                x1 = xyz[c1]
                x2 = xyz[c2]
                cr = np.cross(x1 - x0, x2 - x0)
                area = 0.5 * np.linalg.norm(cr)
                Kg = np.zeros((6, 6), dtype=float)
                fg = np.zeros(6, dtype=float)
                for gL in _TRI_G:
                    N = _tri6_N(gL[0], gL[1], gL[2])
                    w = area / 3.0
                    Kg += h_top * np.outer(N, N) * w
                    fg += h_top * T_amb * N * w
                g6 = np.array([gdof[p] for p in loc], dtype=np.int64)
                hpos = _add_block(g6, Kg, hpos, hr, hc, hd)
                for a in range(6):
                    ga = g6[a]
                    if ga >= 0:
                        fb[ga] += fg[a]

    G = coo_matrix((hd[:hpos], (hr[:hpos], hc[:hpos])), shape=(nh, nh)).tocsr()
    t_hasm = time.time()
    Tn, hsol = _solve_heat_cg(G, fb, verbose=verbose)
    t_hsol = time.time()
    if verbose:
        print("[tet10 heat] assemble=%.2fs  solve=%.2fs (%s)  ndof=%d  nnz=%d" % (
            t_hasm - t0, t_hsol - t_hasm, hsol, nh, G.nnz))

    Tnodal = np.zeros(used.shape, dtype=float)
    Tnodal[used] = T_amb
    for p in range(used_idx.shape[0]):
        i, j, k = int(used_idx[p, 0]), int(used_idx[p, 1]), int(used_idx[p, 2])
        Tnodal[i, j, k] = Tn[hid[i, j, k]]

    # --- elasticity ---
    nnx, nny, nnz = used.shape
    fix, k_bc, triple = _fix_components(bc, used)
    dof_id = -np.ones((nnx, nny, nnz, 3), dtype=int)
    ndof = 0
    n_nodes = int(used_idx.shape[0])
    for p in range(n_nodes):
        i, j, k = int(used_idx[p, 0]), int(used_idx[p, 1]), int(used_idx[p, 2])
        for c in range(3):
            if not fix[i, j, k, c]:
                dof_id[i, j, k, c] = ndof
                ndof += 1
    if verbose:
        print("[tet10 elast] bc=%s  nodes=%d  ndof=%d" % (bc, n_nodes, ndof))
        if triple is not None:
            print("[tet10 elast] 3-2-1 %s k=%d" % (triple, k_bc))
    if ndof == 0:
        raise ValueError("no free mechanical DOFs")

    ke_e_cache = {}
    chunk = 800
    cap = chunk * 30 * 30
    er = np.empty(cap, dtype=np.int32)
    ec = np.empty(cap, dtype=np.int32)
    ed = np.empty(cap, dtype=float)
    epos = 0
    K = None
    ff = np.zeros(ndof, dtype=float)
    missing = set()
    n_built = 0
    hex_tets = [None] * n_hex

    def _flush():
        nonlocal epos, K
        if epos <= 0:
            return
        Kc = coo_matrix((ed[:epos], (er[:epos], ec[:epos])), shape=(ndof, ndof)).tocsr()
        K = Kc if K is None else K + Kc
        epos = 0

    for t in range(n_hex):
        i, j, k = int(si[t]), int(sj[t]), int(sk[t])
        name = mat_names[i, j, k]
        mat = case.material(name) if name else case.material()
        if mat.E is None:
            missing.add(name or "?")
            continue
        dxi = float(grid["dx"][i])
        dyj = float(grid["dy"][j])
        dzk = float(grid["dz"][k])
        pack = []
        for ti, tet in enumerate(_TETS6):
            nodes = _tet10_nodes(i, j, k, tet)
            xyz = _tet_xyz(nodes, xe, ye, ze)
            ck = (round(dxi, 16), round(dyj, 16), round(dzk, 16), ti, mat.E, mat.nu)
            hit = ke_e_cache.get(ck)
            if hit is None:
                C = _iso_C(mat.E, mat.nu)
                Ke, gps = _ke_elast(xyz, C)
                ke_e_cache[ck] = (Ke, gps, C)
            else:
                Ke, gps, C = hit
            gdof = np.empty(30, dtype=np.int64)
            T10 = np.empty(10, dtype=float)
            for a, nd in enumerate(nodes):
                I, J, Kk = nd
                T10[a] = Tnodal[I, J, Kk]
                gdof[3 * a] = dof_id[I, J, Kk, 0]
                gdof[3 * a + 1] = dof_id[I, J, Kk, 1]
                gdof[3 * a + 2] = dof_id[I, J, Kk, 2]
            if epos + 900 > cap:
                _flush()
            epos = _add_block(gdof, Ke, epos, er, ec, ed)
            eth = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
            for N, B, wdet in gps:
                Tgp = float(N.dot(T10))
                aDT = mat.alpha_cte * (Tgp - T_ref)
                f30 = B.T.dot(C.dot(eth) * aDT) * wdet
                free = gdof >= 0
                gv = gdof[free]
                if gv.size:
                    ff[gv] += f30[free]
            pack.append((nodes, T10, ck, mat.alpha_cte))
            n_built += 1
        hex_tets[t] = pack
    if missing:
        raise ValueError("materials missing E_Pa: %s" % sorted(missing))
    _flush()
    if K is None:
        raise ValueError("empty elasticity matrix")
    t_easm = time.time()
    if verbose:
        print("[tet10 elast] assembled %d tets  nnz=%d  %.1fs" % (
            n_built, K.nnz, t_easm - t_hsol))
    u, esol = _solve_cg(K, ff, verbose=verbose)
    t_esol = time.time()

    ux = np.zeros((nnx, nny, nnz))
    uy = np.zeros((nnx, nny, nnz))
    uz = np.zeros((nnx, nny, nnz))
    for p in range(n_nodes):
        i, j, k = int(used_idx[p, 0]), int(used_idx[p, 1]), int(used_idx[p, 2])
        dxid, dyid, dzid = dof_id[i, j, k, 0], dof_id[i, j, k, 1], dof_id[i, j, k, 2]
        if dxid >= 0:
            ux[i, j, k] = u[dxid]
        if dyid >= 0:
            uy[i, j, k] = u[dyid]
        if dzid >= 0:
            uz[i, j, k] = u[dzid]

    Tcell = np.full((nx, ny, nz), T_amb)
    ucx = np.zeros((nx, ny, nz))
    ucy = np.zeros((nx, ny, nz))
    ucz = np.zeros((nx, ny, nz))
    vm = np.zeros((nx, ny, nz))
    sxx = np.zeros((nx, ny, nz))
    syy = np.zeros((nx, ny, nz))
    szz = np.zeros((nx, ny, nz))
    sxy = np.zeros((nx, ny, nz))
    Lc = 0.25
    Nc = _tet_N(Lc, Lc, Lc, Lc)
    for t in range(n_hex):
        i, j, k = int(si[t]), int(sj[t]), int(sk[t])
        pack = hex_tets[t]
        if not pack:
            continue
        vol = 0.0
        tacc = uacc = vacc = wacc = 0.0
        sacc = np.zeros(6, dtype=float)
        vmacc = 0.0
        for nodes, T10, ck, alpha in pack:
            Ke, gps, C = ke_e_cache[ck]
            u30 = np.empty(30, dtype=float)
            for a, nd in enumerate(nodes):
                I, J, Kk = nd
                u30[3 * a] = ux[I, J, Kk]
                u30[3 * a + 1] = uy[I, J, Kk]
                u30[3 * a + 2] = uz[I, J, Kk]
            v_e = 0.0
            s_e = np.zeros(6, dtype=float)
            for N, B, wdet in gps:
                Tgp = float(N.dot(T10))
                aDT = alpha * (Tgp - T_ref)
                eth = aDT * np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
                s_e += C.dot(B.dot(u30) - eth) * wdet
                v_e += wdet
            if v_e <= 0:
                continue
            s_e /= v_e
            Tcent = float(Nc.dot(T10))
            ucent = np.array([
                Nc.dot(u30[0::3]),
                Nc.dot(u30[1::3]),
                Nc.dot(u30[2::3]),
            ])
            vol += v_e
            tacc += Tcent * v_e
            uacc += ucent[0] * v_e
            vacc += ucent[1] * v_e
            wacc += ucent[2] * v_e
            sacc += s_e * v_e
            vmacc += _von_mises(s_e) * v_e
        if vol <= 0:
            continue
        Tcell[i, j, k] = tacc / vol
        ucx[i, j, k] = uacc / vol
        ucy[i, j, k] = vacc / vol
        ucz[i, j, k] = wacc / vol
        sv = sacc / vol
        sxx[i, j, k], syy[i, j, k], szz[i, j, k] = sv[0], sv[1], sv[2]
        sxy[i, j, k] = sv[3]
        vm[i, j, k] = vmacc / vol

    umag = np.sqrt(ucx ** 2 + ucy ** 2 + ucz ** 2)
    elapsed = time.time() - t0
    Ts = Tcell[solid]
    heat_timing = {
        "grid_s": 0.0,
        "assemble_s": t_hasm - t0,
        "solve_s": t_hsol - t_hasm,
        "total_s": t_hsol - t0,
        "n_iter": 1,
        "n_cells": nx * ny * nz,
        "n_solid": n_hex,
        "n_unknown": nh,
        "nnz": int(G.nnz),
        "solver": hsol,
        "grid": [nx, ny, nz],
        "element": "tet10",
    }
    elast_timing = {
        "assemble_s": t_easm - t_hsol,
        "solve_s": t_esol - t_easm,
        "post_s": elapsed - (t_esol - t0),
        "total_s": elapsed - (t_hsol - t0),
        "ndof": int(ndof),
        "n_nodes": int(n_nodes),
        "n_hex": n_hex,
        "n_tet": n_built,
        "n_solid": n_hex,
        "nnz": int(K.nnz),
        "solver": esol,
        "bc": bc,
        "element": "tet10",
    }
    if verbose:
        print("[tet10] Tmax=%.3f C  Tmin=%.3f C  |u|_max=%.3e  vm_max=%.3e  %.1fs" % (
            float(Ts.max() - 273.15), float(Ts.min() - 273.15),
            float(umag[solid].max()), float(vm[solid].max()), elapsed))
        print("[tet10 elast] assemble=%.2fs  solve=%.2fs (%s)  post=%.2fs  ndof=%d" % (
            elast_timing["assemble_s"], elast_timing["solve_s"], esol,
            elast_timing["post_s"], ndof))

    thermal = {
        "T": Tcell,
        "grid": grid,
        "Pvol": Pvol,
        "kappa": None,
        "mat_names": mat_names,
        "elapsed_s": t_hsol - t0,
        "history": [],
        "case": case,
        "solver": "tet10",
        "timing": heat_timing,
    }
    elast = {
        "ux": ux, "uy": uy, "uz": uz,
        "ucx": ucx, "ucy": ucy, "ucz": ucz, "umag": umag,
        "vm": vm, "sxx": sxx, "syy": syy, "szz": szz, "sxy": sxy,
        "T_ref": T_ref, "elapsed_s": elapsed - (t_hsol - t0), "ndof": ndof,
        "thermal": thermal, "case": case, "grid": grid, "bc": bc,
        "timing": elast_timing, "element": "tet10",
    }
    return thermal, elast
