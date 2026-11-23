# -*- coding: utf-8 -*-
"""Linear thermoelasticity on the same Cartesian voxels as FVM heat.

    div sigma = 0
    sigma = C : (eps - alpha (T - T_ref) I)
    eps = symgrad u

Default element is Wilson/QM6 8-node hex (one per solid FVM cell).
Optional hex8 (plain linear) and hex20 (serendipity quadratic on the
same voxels, edge mid-nodes). Heat stays Cartesian FVM. One-way:
T from fvm3d is a load.

Default BC is substrate_bottom_fixed. free_321 and die_bottom_fixed remain.
"""

from __future__ import print_function

import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve, cg, LinearOperator

from .geometry import material_field

_XI = np.array([-1, 1, 1, -1, -1, 1, 1, -1], dtype=float)
_ET = np.array([-1, -1, 1, 1, -1, -1, 1, 1], dtype=float)
_ZT = np.array([-1, -1, -1, -1, 1, 1, 1, 1], dtype=float)
_G = (-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0))
_LOC = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
        (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))

# hex20 serendipity: 8 corners then 12 edge mids. Offsets on a 2-step
# fine index so corners land on even nodes.
_HEX20_OFF = (
    (0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0),
    (0, 0, 2), (2, 0, 2), (2, 2, 2), (0, 2, 2),
    (1, 0, 0), (2, 1, 0), (1, 2, 0), (0, 1, 0),
    (1, 0, 2), (2, 1, 2), (1, 2, 2), (0, 1, 2),
    (0, 0, 1), (2, 0, 1), (2, 2, 1), (0, 2, 1),
)
_HEX20_XI = np.array([-1, 1, 1, -1, -1, 1, 1, -1,
                      0, 1, 0, -1, 0, 1, 0, -1, -1, 1, 1, -1], dtype=float)
_HEX20_ET = np.array([-1, -1, 1, 1, -1, -1, 1, 1,
                      -1, 0, 1, 0, -1, 0, 1, 0, -1, -1, 1, 1], dtype=float)
_HEX20_ZT = np.array([-1, -1, -1, -1, 1, 1, 1, 1,
                      -1, -1, -1, -1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=float)
_G3 = (-np.sqrt(0.6), 0.0, np.sqrt(0.6))
_W3 = (5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0)


def _iso_C(E, nu):
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    C = np.zeros((6, 6), dtype=float)
    C[0, 0] = C[1, 1] = C[2, 2] = lam + 2.0 * mu
    C[0, 1] = C[0, 2] = C[1, 0] = C[1, 2] = C[2, 0] = C[2, 1] = lam
    C[3, 3] = C[4, 4] = C[5, 5] = mu
    return C


def _fill_B(B, col, dNx, dNy, dNz):
    B[0, col] = dNx
    B[1, col + 1] = dNy
    B[2, col + 2] = dNz
    B[3, col] = dNy
    B[3, col + 1] = dNx
    B[4, col + 1] = dNz
    B[4, col + 2] = dNy
    B[5, col] = dNz
    B[5, col + 2] = dNx


def _hex_Ke_fe(dx, dy, dz, C, element="wilson"):
    """8-node hex stiffness. wilson = QM6 incompatible modes (bending)."""
    Ke_uu = np.zeros((24, 24), dtype=float)
    fe_u = np.zeros(24, dtype=float)
    detJ = dx * dy * dz / 8.0
    invJ = np.array([2.0 / dx, 2.0 / dy, 2.0 / dz])
    eth = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    Ceth = C.dot(eth)
    Bu = np.zeros((6, 24), dtype=float)
    use_inc = (element != "hex8")
    if use_inc:
        Ke_ua = np.zeros((24, 9), dtype=float)
        Ke_aa = np.zeros((9, 9), dtype=float)
        fe_a = np.zeros(9, dtype=float)
        Ba = np.zeros((6, 9), dtype=float)
    for gx in _G:
        for gy in _G:
            for gz in _G:
                Bu[:, :] = 0.0
                for a in range(8):
                    xa, ya, za = _XI[a], _ET[a], _ZT[a]
                    dN_xi = 0.125 * xa * (1.0 + ya * gy) * (1.0 + za * gz)
                    dN_et = 0.125 * ya * (1.0 + xa * gx) * (1.0 + za * gz)
                    dN_zt = 0.125 * za * (1.0 + xa * gx) * (1.0 + ya * gy)
                    _fill_B(Bu, 3 * a, dN_xi * invJ[0], dN_et * invJ[1], dN_zt * invJ[2])
                Ke_uu += Bu.T.dot(C.dot(Bu)) * detJ
                fe_u += Bu.T.dot(Ceth) * detJ
                if use_inc:
                    Ba[:, :] = 0.0
                    _fill_B(Ba, 0, -2.0 * gx * invJ[0], 0.0, 0.0)
                    _fill_B(Ba, 3, 0.0, -2.0 * gy * invJ[1], 0.0)
                    _fill_B(Ba, 6, 0.0, 0.0, -2.0 * gz * invJ[2])
                    Ke_ua += Bu.T.dot(C.dot(Ba)) * detJ
                    Ke_aa += Ba.T.dot(C.dot(Ba)) * detJ
                    fe_a += Ba.T.dot(Ceth) * detJ
    if not use_inc:
        return Ke_uu, None
    try:
        Kaa_inv = np.linalg.inv(Ke_aa)
    except np.linalg.LinAlgError:
        return Ke_uu, None
    Ke = Ke_uu - Ke_ua.dot(Kaa_inv.dot(Ke_ua.T))
    return Ke, (Kaa_inv, Ke_ua)


def _hex8_N(gx, gy, gz):
    N = np.empty(8, dtype=float)
    for a in range(8):
        N[a] = 0.125 * (1.0 + _XI[a] * gx) * (1.0 + _ET[a] * gy) * (1.0 + _ZT[a] * gz)
    return N


def _nodal_T_hex8(T, solid, T_amb):
    """Average cell-centered FVM T onto hex8 corners (COMSOL-style nodal T)."""
    nx, ny, nz = T.shape
    Tn = np.zeros((nx + 1, ny + 1, nz + 1), dtype=float)
    cnt = np.zeros((nx + 1, ny + 1, nz + 1), dtype=float)
    si, sj, sk = np.where(solid)
    for t in range(si.size):
        i, j, k = int(si[t]), int(sj[t]), int(sk[t])
        Tn[i:i + 2, j:j + 2, k:k + 2] += float(T[i, j, k])
        cnt[i:i + 2, j:j + 2, k:k + 2] += 1.0
    m = cnt > 0
    Tn[m] /= cnt[m]
    Tn[~m] = T_amb
    return Tn


def _hex_fe_from_T8(dx, dy, dz, C, T8, T_ref, alpha, element="wilson", wpack=None):
    """Thermal load with T interpolated at 2x2x2 Gauss points."""
    detJ = dx * dy * dz / 8.0
    invJ = np.array([2.0 / dx, 2.0 / dy, 2.0 / dz])
    fe_u = np.zeros(24, dtype=float)
    use_inc = (element != "hex8") and (wpack is not None)
    if use_inc:
        fe_a = np.zeros(9, dtype=float)
        Ba = np.zeros((6, 9), dtype=float)
    Bu = np.zeros((6, 24), dtype=float)
    for gx in _G:
        for gy in _G:
            for gz in _G:
                N = _hex8_N(gx, gy, gz)
                aDT = alpha * (float(N.dot(T8)) - T_ref)
                Ceth = C.dot(np.array([aDT, aDT, aDT, 0.0, 0.0, 0.0]))
                Bu[:, :] = 0.0
                for a in range(8):
                    xa, ya, za = _XI[a], _ET[a], _ZT[a]
                    dN_xi = 0.125 * xa * (1.0 + ya * gy) * (1.0 + za * gz)
                    dN_et = 0.125 * ya * (1.0 + xa * gx) * (1.0 + za * gz)
                    dN_zt = 0.125 * za * (1.0 + xa * gx) * (1.0 + ya * gy)
                    _fill_B(Bu, 3 * a, dN_xi * invJ[0], dN_et * invJ[1], dN_zt * invJ[2])
                fe_u += Bu.T.dot(Ceth) * detJ
                if use_inc:
                    Ba[:, :] = 0.0
                    _fill_B(Ba, 0, -2.0 * gx * invJ[0], 0.0, 0.0)
                    _fill_B(Ba, 3, 0.0, -2.0 * gy * invJ[1], 0.0)
                    _fill_B(Ba, 6, 0.0, 0.0, -2.0 * gz * invJ[2])
                    fe_a += Ba.T.dot(Ceth) * detJ
    if not use_inc:
        return fe_u
    Kaa_inv, Ke_ua = wpack
    return fe_u - Ke_ua.dot(Kaa_inv.dot(fe_a))


def _hex_stress_gauss(u24, dx, dy, dz, C, T8, T_ref, alpha):
    """Compatible-B stress averaged on 8 Gauss points (COMSOL linear-hex style)."""
    invJ = np.array([2.0 / dx, 2.0 / dy, 2.0 / dz])
    B = np.zeros((6, 24), dtype=float)
    acc = np.zeros(6, dtype=float)
    n = 0
    for gx in _G:
        for gy in _G:
            for gz in _G:
                N = _hex8_N(gx, gy, gz)
                aDT = alpha * (float(N.dot(T8)) - T_ref)
                B[:, :] = 0.0
                for a in range(8):
                    xa, ya, za = _XI[a], _ET[a], _ZT[a]
                    dNx = 0.125 * xa * (1.0 + ya * gy) * (1.0 + za * gz) * invJ[0]
                    dNy = 0.125 * ya * (1.0 + xa * gx) * (1.0 + za * gz) * invJ[1]
                    dNz = 0.125 * za * (1.0 + xa * gx) * (1.0 + ya * gy) * invJ[2]
                    _fill_B(B, 3 * a, dNx, dNy, dNz)
                eth = aDT * np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
                acc += C.dot(B.dot(u24) - eth)
                n += 1
    return acc / float(n)


def _normalize_element(element):
    s = str(element or "wilson").lower()
    if s in ("hex8", "linear", "q8"):
        return "hex8"
    if s in ("hex20", "quadratic", "serendipity", "q20"):
        return "hex20"
    if s in ("tet10", "quadratic_tet", "t10"):
        return "tet10"
    if s in ("wilson", "qm6", "hex8w"):
        return "wilson"
    raise ValueError(
        "unknown mechanics.element %r (use wilson, hex8, hex20, tet10)" % element
    )


def _hex20_dN(xi, et, zt):
    """dN/d(xi,eta,zeta) for serendipity hex20, shape (20, 3)."""
    dN = np.zeros((20, 3), dtype=float)
    for a in range(8):
        xa, ya, za = _HEX20_XI[a], _HEX20_ET[a], _HEX20_ZT[a]
        A = 1.0 + xa * xi
        B = 1.0 + ya * et
        C = 1.0 + za * zt
        dN[a, 0] = 0.125 * xa * B * C * (2.0 * xa * xi + ya * et + za * zt - 1.0)
        dN[a, 1] = 0.125 * ya * A * C * (xa * xi + 2.0 * ya * et + za * zt - 1.0)
        dN[a, 2] = 0.125 * za * A * B * (xa * xi + ya * et + 2.0 * za * zt - 1.0)
    for a in range(8, 20):
        xa, ya, za = _HEX20_XI[a], _HEX20_ET[a], _HEX20_ZT[a]
        if xa == 0.0:
            dN[a, 0] = -0.5 * xi * (1.0 + ya * et) * (1.0 + za * zt)
            dN[a, 1] = 0.25 * (1.0 - xi * xi) * ya * (1.0 + za * zt)
            dN[a, 2] = 0.25 * (1.0 - xi * xi) * (1.0 + ya * et) * za
        elif ya == 0.0:
            dN[a, 0] = 0.25 * (1.0 - et * et) * xa * (1.0 + za * zt)
            dN[a, 1] = -0.5 * et * (1.0 + xa * xi) * (1.0 + za * zt)
            dN[a, 2] = 0.25 * (1.0 - et * et) * (1.0 + xa * xi) * za
        else:
            dN[a, 0] = 0.25 * (1.0 - zt * zt) * xa * (1.0 + ya * et)
            dN[a, 1] = 0.25 * (1.0 - zt * zt) * (1.0 + xa * xi) * ya
            dN[a, 2] = -0.5 * zt * (1.0 + xa * xi) * (1.0 + ya * et)
    return dN


def _hex20_N(xi, et, zt):
    N = np.empty(20, dtype=float)
    for a in range(8):
        xa, ya, za = _HEX20_XI[a], _HEX20_ET[a], _HEX20_ZT[a]
        N[a] = 0.125 * (1.0 + xa * xi) * (1.0 + ya * et) * (1.0 + za * zt) * (
            xa * xi + ya * et + za * zt - 2.0)
    for a in range(8, 20):
        xa, ya, za = _HEX20_XI[a], _HEX20_ET[a], _HEX20_ZT[a]
        if xa == 0.0:
            N[a] = 0.25 * (1.0 - xi * xi) * (1.0 + ya * et) * (1.0 + za * zt)
        elif ya == 0.0:
            N[a] = 0.25 * (1.0 - et * et) * (1.0 + xa * xi) * (1.0 + za * zt)
        else:
            N[a] = 0.25 * (1.0 - zt * zt) * (1.0 + xa * xi) * (1.0 + ya * et)
    return N


def _hex20_Ke_fe(dx, dy, dz, C):
    Ke = np.zeros((60, 60), dtype=float)
    detJ = dx * dy * dz / 8.0
    invJ = np.array([2.0 / dx, 2.0 / dy, 2.0 / dz])
    B = np.zeros((6, 60), dtype=float)
    for i, gx in enumerate(_G3):
        for j, gy in enumerate(_G3):
            for k, gz in enumerate(_G3):
                w = _W3[i] * _W3[j] * _W3[k] * detJ
                dN = _hex20_dN(gx, gy, gz)
                B[:, :] = 0.0
                for a in range(20):
                    _fill_B(B, 3 * a, dN[a, 0] * invJ[0], dN[a, 1] * invJ[1], dN[a, 2] * invJ[2])
                Ke += B.T.dot(C.dot(B)) * w
    return Ke


def _hex20_fe_from_T20(dx, dy, dz, C, T20, T_ref, alpha):
    fe = np.zeros(60, dtype=float)
    detJ = dx * dy * dz / 8.0
    invJ = np.array([2.0 / dx, 2.0 / dy, 2.0 / dz])
    B = np.zeros((6, 60), dtype=float)
    for i, gx in enumerate(_G3):
        for j, gy in enumerate(_G3):
            for k, gz in enumerate(_G3):
                w = _W3[i] * _W3[j] * _W3[k] * detJ
                aDT = alpha * (float(_hex20_N(gx, gy, gz).dot(T20)) - T_ref)
                Ceth = C.dot(np.array([aDT, aDT, aDT, 0.0, 0.0, 0.0]))
                dN = _hex20_dN(gx, gy, gz)
                B[:, :] = 0.0
                for a in range(20):
                    _fill_B(B, 3 * a, dN[a, 0] * invJ[0], dN[a, 1] * invJ[1], dN[a, 2] * invJ[2])
                fe += B.T.dot(Ceth) * w
    return fe


def _hex20_stress(u60, dx, dy, dz, C, T20, T_ref, alpha):
    invJ = np.array([2.0 / dx, 2.0 / dy, 2.0 / dz])
    B = np.zeros((6, 60), dtype=float)
    acc = np.zeros(6, dtype=float)
    n = 0
    for gx in _G:
        for gy in _G:
            for gz in _G:
                aDT = alpha * (float(_hex20_N(gx, gy, gz).dot(T20)) - T_ref)
                eth = aDT * np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
                dN = _hex20_dN(gx, gy, gz)
                B[:, :] = 0.0
                for a in range(20):
                    _fill_B(B, 3 * a, dN[a, 0] * invJ[0], dN[a, 1] * invJ[1], dN[a, 2] * invJ[2])
                acc += C.dot(B.dot(u60) - eth)
                n += 1
    return acc / float(n)


def _nodal_T_hex20(T, solid, T_amb):
    nx, ny, nz = T.shape
    nnx, nny, nnz = 2 * nx + 1, 2 * ny + 1, 2 * nz + 1
    Tn = np.zeros((nnx, nny, nnz), dtype=float)
    cnt = np.zeros((nnx, nny, nnz), dtype=float)
    si, sj, sk = np.where(solid)
    for t in range(si.size):
        i, j, k = int(si[t]), int(sj[t]), int(sk[t])
        tv = float(T[i, j, k])
        for di, dj, dk in _HEX20_OFF:
            Tn[2 * i + di, 2 * j + dj, 2 * k + dk] += tv
            cnt[2 * i + di, 2 * j + dj, 2 * k + dk] += 1.0
    m = cnt > 0
    Tn[m] /= cnt[m]
    Tn[~m] = T_amb
    return Tn


def _von_mises(s):
    sxx, syy, szz, sxy, syz, sxz = s
    return np.sqrt(
        0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
        + 3.0 * (sxy ** 2 + syz ** 2 + sxz ** 2)
    )


def _lowest_used_k(node_used):
    for k in range(node_used.shape[2]):
        if node_used[:, :, k].any():
            return k
    raise ValueError("no used mechanical nodes")


def _corner_triple(ij_list):
    if len(ij_list) < 3:
        raise ValueError("need >=3 nodes for 3-2-1 BC, got %d" % len(ij_list))
    i_min = min(p[0] for p in ij_list)
    i_max = max(p[0] for p in ij_list)
    j_min = min(p[1] for p in ij_list)
    j_max = max(p[1] for p in ij_list)

    def nearest(i0, j0):
        best = ij_list[0]
        bd = (best[0] - i0) ** 2 + (best[1] - j0) ** 2
        for p in ij_list:
            d = (p[0] - i0) ** 2 + (p[1] - j0) ** 2
            if d < bd:
                best = p
                bd = d
        return best

    a = nearest(i_min, j_min)
    b = nearest(i_max, j_min)
    c = nearest(i_min, j_max)
    if a == b or a == c or b == c:
        uniq = []
        for p in sorted(ij_list):
            if p not in uniq:
                uniq.append(p)
            if len(uniq) == 3:
                break
        if len(uniq) < 3:
            raise ValueError("cannot pick 3 distinct 3-2-1 nodes")
        a, b, c = uniq[0], uniq[1], uniq[2]
    return a, b, c


def _fix_components(bc, node_used):
    """Per-node, per-component essential BC. True = constrained to 0."""
    nnx, nny, nnz = node_used.shape
    fix = np.zeros((nnx, nny, nnz, 3), dtype=bool)
    k0 = _lowest_used_k(node_used)
    if bc in ("die_bottom_fixed", "substrate_bottom_fixed", "bottom_fixed"):
        for i in range(nnx):
            for j in range(nny):
                if node_used[i, j, k0]:
                    fix[i, j, k0, 0] = True
                    fix[i, j, k0, 1] = True
                    fix[i, j, k0, 2] = True
        return fix, k0, None
    if bc in ("free_321", "321"):
        ii, jj = np.where(node_used[:, :, k0])
        pts = list(zip(ii.tolist(), jj.tolist()))
        a, b, c = _corner_triple(pts)
        fix[a[0], a[1], k0, 0] = True
        fix[a[0], a[1], k0, 1] = True
        fix[a[0], a[1], k0, 2] = True
        fix[b[0], b[1], k0, 1] = True
        fix[b[0], b[1], k0, 2] = True
        fix[c[0], c[1], k0, 2] = True
        return fix, k0, (a, b, c)
    raise ValueError(
        "unknown mechanical bc %r (use substrate_bottom_fixed, free_321, or die_bottom_fixed)" % bc
    )


def _solve_cg(K, f, verbose=True):
    ndof = f.size
    nnz = int(K.nnz)
    # Direct factorization of Wilson hex systems (~1e5 dof, ~1e7 nnz)
    # peaks well above 12 GB. Keep CG even if it hits maxiter.
    skip_spsolve = (ndof > 80000) or (nnz > 3000000)
    u = None
    mech_solver = "spsolve"
    maxiter = 12000 if skip_spsolve else (8000 if ndof > 250000 else 4000)
    try:
        diag = K.diagonal()
        Minv = 1.0 / np.where(np.abs(diag) > 0, diag, 1.0)
        Aop = LinearOperator((ndof, ndof), matvec=lambda x: K.dot(x))
        Mop = LinearOperator((ndof, ndof), matvec=lambda x: Minv * x)
        if verbose:
            print("[elast] CG ndof=%d nnz=%d maxiter=%d skip_spsolve=%s" % (
                ndof, nnz, maxiter, skip_spsolve), flush=True)
        u, info = cg(Aop, f, M=Mop, rtol=1e-8, atol=0.0, maxiter=maxiter)
        mech_solver = "cg"
        if info != 0:
            if skip_spsolve:
                if verbose:
                    print("[warn] CG info=%s  ndof=%d nnz=%d  keep iterate (no spsolve)" % (
                        info, ndof, nnz), flush=True)
                mech_solver = "cg_unconverged"
            else:
                if verbose:
                    print("[warn] CG info=%s, spsolve" % info, flush=True)
                u = None
                mech_solver = "spsolve_after_cg"
    except Exception as e:
        if verbose:
            print("[warn] CG failed (%s)%s" % (
                e, "" if skip_spsolve else ", spsolve"), flush=True)
        u = None
        mech_solver = "spsolve_after_cg"
    if u is None:
        if skip_spsolve:
            raise MemoryError(
                "elasticity CG failed on ndof=%d nnz=%d; spsolve skipped" % (
                    ndof, nnz))
        u = spsolve(K, f)
        if mech_solver == "cg":
            mech_solver = "spsolve"
    return u, mech_solver


def _solve_elasticity_hex20(thermal, T_ref, verbose, bc, t0):
    """Quadratic serendipity hex on the same voxels as FVM heat."""
    case = thermal["case"]
    grid = thermal["grid"]
    T = thermal["T"]
    nx, ny, nz = int(grid["nx"]), int(grid["ny"]), int(grid["nz"])
    dx = np.asarray(grid["dx"], dtype=float)
    dy = np.asarray(grid["dy"], dtype=float)
    dz = np.asarray(grid["dz"], dtype=float)
    if dx.ndim == 0:
        dx = np.full(nx, float(dx))
    if dy.ndim == 0:
        dy = np.full(ny, float(dy))
    if dz.ndim == 0:
        dz = np.full(nz, float(dz))
    solid = grid.get("solid")
    if solid is None:
        solid = np.ones((nx, ny, nz), dtype=bool)
    mat_names = thermal.get("mat_names")
    if mat_names is None:
        mat_names = material_field(case, grid)

    nnx, nny, nnz = 2 * nx + 1, 2 * ny + 1, 2 * nz + 1
    node_used = np.zeros((nnx, nny, nnz), dtype=bool)
    si, sj, sk = np.where(solid)
    for di, dj, dk in _HEX20_OFF:
        node_used[2 * si + di, 2 * sj + dj, 2 * sk + dk] = True

    fix, k_bc, triple = _fix_components(bc, node_used)
    dof_id = -np.ones((nnx, nny, nnz, 3), dtype=int)
    ndof = 0
    n_nodes = 0
    used_idx = np.argwhere(node_used)
    n_nodes = int(used_idx.shape[0])
    for p in range(n_nodes):
        i, j, k = int(used_idx[p, 0]), int(used_idx[p, 1]), int(used_idx[p, 2])
        for c in range(3):
            if not fix[i, j, k, c]:
                dof_id[i, j, k, c] = ndof
                ndof += 1
    if verbose:
        print("[elast] bc=%s  element=hex20  nodes=%d  ndof=%d  solid cells=%d" % (
            bc, n_nodes, ndof, int(solid.sum())))
        if triple is not None:
            print("[elast] 3-2-1 A=%s B=%s C=%s  k=%d" % (
                triple[0], triple[1], triple[2], k_bc))
    if ndof == 0:
        raise ValueError("no free mechanical DOFs")

    ke_cache = {}

    def get_ke(dxi, dyj, dzk, E, nu):
        key = (round(dxi, 16), round(dyj, 16), round(dzk, 16), E, nu)
        hit = ke_cache.get(key)
        if hit is None:
            C = _iso_C(E, nu)
            Ke = _hex20_Ke_fe(dxi, dyj, dzk, C)
            hit = (Ke, C)
            ke_cache[key] = hit
        return hit

    Tn = _nodal_T_hex20(T, solid, T_ref)
    n_el = int(si.size)
    nnz_est = n_el * 60 * 60
    rows = np.empty(nnz_est, dtype=np.int64)
    cols = np.empty(nnz_est, dtype=np.int64)
    data = np.empty(nnz_est, dtype=float)
    pos = 0
    f = np.zeros(ndof, dtype=float)
    missing = set()
    n_built = 0
    for t in range(n_el):
        i, j, k = int(si[t]), int(sj[t]), int(sk[t])
        name = mat_names[i, j, k]
        mat = case.material(name) if name else case.material()
        if mat.E is None:
            missing.add(name or "?")
            continue
        Ke, C = get_ke(dx[i], dy[j], dz[k], mat.E, mat.nu)
        T20 = np.empty(20, dtype=float)
        gdof = np.empty(60, dtype=np.int64)
        for a, (di, dj, dk) in enumerate(_HEX20_OFF):
            ii, jj, kk = 2 * i + di, 2 * j + dj, 2 * k + dk
            T20[a] = Tn[ii, jj, kk]
            gdof[3 * a] = dof_id[ii, jj, kk, 0]
            gdof[3 * a + 1] = dof_id[ii, jj, kk, 1]
            gdof[3 * a + 2] = dof_id[ii, jj, kk, 2]
        fe = _hex20_fe_from_T20(dx[i], dy[j], dz[k], C, T20, T_ref, mat.alpha_cte)
        free = gdof >= 0
        gv = gdof[free]
        if gv.size:
            f[gv] += fe[free]
            Ke_s = Ke[np.ix_(free, free)]
            n = gv.size
            block = n * n
            rows[pos:pos + block] = np.repeat(gv, n)
            cols[pos:pos + block] = np.tile(gv, n)
            data[pos:pos + block] = Ke_s.ravel()
            pos += block
        n_built += 1
    if missing:
        raise ValueError("materials missing E_Pa: %s" % sorted(missing))
    K = coo_matrix((data[:pos], (rows[:pos], cols[:pos])), shape=(ndof, ndof)).tocsr()
    t_asm = time.time()
    if verbose:
        print("[elast] assembled %d hex20  nnz=%d  %.1fs" % (
            n_built, K.nnz, t_asm - t0))

    u, mech_solver = _solve_cg(K, f, verbose=verbose)
    t_solve = time.time()

    ux = np.zeros((nnx, nny, nnz))
    uy = np.zeros((nnx, nny, nnz))
    uz = np.zeros((nnx, nny, nnz))
    for p in range(n_nodes):
        i, j, k = int(used_idx[p, 0]), int(used_idx[p, 1]), int(used_idx[p, 2])
        dxid = dof_id[i, j, k, 0]
        dyid = dof_id[i, j, k, 1]
        dzid = dof_id[i, j, k, 2]
        if dxid >= 0:
            ux[i, j, k] = u[dxid]
        if dyid >= 0:
            uy[i, j, k] = u[dyid]
        if dzid >= 0:
            uz[i, j, k] = u[dzid]

    # hex20 N at cell center (0,0,0): corners -1/4, mids +1/4
    ucx = np.zeros((nx, ny, nz))
    ucy = np.zeros((nx, ny, nz))
    ucz = np.zeros((nx, ny, nz))
    vm = np.zeros((nx, ny, nz))
    sxx = np.zeros((nx, ny, nz))
    syy = np.zeros((nx, ny, nz))
    szz = np.zeros((nx, ny, nz))
    sxy = np.zeros((nx, ny, nz))
    wN = np.empty(20, dtype=float)
    wN[:8] = -0.25
    wN[8:] = 0.25
    for t in range(n_el):
        i, j, k = int(si[t]), int(sj[t]), int(sk[t])
        u60 = np.empty(60, dtype=float)
        px = py_ = pz = 0.0
        for a, (di, dj, dk) in enumerate(_HEX20_OFF):
            ii, jj, kk = 2 * i + di, 2 * j + dj, 2 * k + dk
            uxa, uya, uza = ux[ii, jj, kk], uy[ii, jj, kk], uz[ii, jj, kk]
            u60[3 * a] = uxa
            u60[3 * a + 1] = uya
            u60[3 * a + 2] = uza
            px += wN[a] * uxa
            py_ += wN[a] * uya
            pz += wN[a] * uza
        ucx[i, j, k] = px
        ucy[i, j, k] = py_
        ucz[i, j, k] = pz
        name = mat_names[i, j, k]
        mat = case.material(name) if name else case.material()
        C = _iso_C(mat.E, mat.nu)
        T20 = np.empty(20, dtype=float)
        for a, (di, dj, dk) in enumerate(_HEX20_OFF):
            T20[a] = Tn[2 * i + di, 2 * j + dj, 2 * k + dk]
        sv = _hex20_stress(u60, dx[i], dy[j], dz[k], C, T20, T_ref, mat.alpha_cte)
        sxx[i, j, k], syy[i, j, k], szz[i, j, k] = sv[0], sv[1], sv[2]
        sxy[i, j, k] = sv[3]
        vm[i, j, k] = _von_mises(sv)

    umag = np.sqrt(ucx ** 2 + ucy ** 2 + ucz ** 2)
    elapsed = time.time() - t0
    timing = {
        "assemble_s": t_asm - t0,
        "solve_s": t_solve - t_asm,
        "post_s": elapsed - (t_solve - t0),
        "total_s": elapsed,
        "ndof": int(ndof),
        "n_nodes": int(n_nodes),
        "n_hex": int(n_built),
        "n_solid": int(solid.sum()),
        "nnz": int(K.nnz),
        "solver": mech_solver,
        "bc": bc,
        "element": "hex20",
    }
    if verbose:
        m = solid
        print("[elast] |u|_max=%.3e m  uz_max=%.3e m  vm_max=%.3e Pa  %.1fs" % (
            umag[m].max(), np.abs(ucz[m]).max(), vm[m].max(), elapsed))
        print("[elast] timing  assemble=%.2fs  solve=%.2fs (%s)  post=%.2fs  ndof=%d  nnz=%d" % (
            timing["assemble_s"], timing["solve_s"], mech_solver,
            timing["post_s"], ndof, K.nnz))
    return {
        "ux": ux, "uy": uy, "uz": uz,
        "ucx": ucx, "ucy": ucy, "ucz": ucz, "umag": umag,
        "vm": vm, "sxx": sxx, "syy": syy, "szz": szz, "sxy": sxy,
        "T_ref": T_ref, "elapsed_s": elapsed, "ndof": ndof,
        "thermal": thermal, "case": case, "grid": grid, "bc": bc,
        "timing": timing, "element": "hex20",
    }


def solve_elasticity(thermal, T_ref=None, verbose=True, bc=None, element=None):
    t0 = time.time()
    case = thermal["case"]
    grid = thermal["grid"]
    T = thermal["T"]
    if T_ref is None:
        T_ref = case.T_amb
    mech = case.raw.get("mechanics") or {}
    if bc is None:
        bc = mech.get("bc", "substrate_bottom_fixed")
    if element is None:
        element = mech.get("element", "wilson")
    element = _normalize_element(element)
    if element == "tet10":
        raise ValueError("tet10 is solve_tet10() in atsim.fem_tet, not solve_elasticity")
    if element == "hex20":
        return _solve_elasticity_hex20(thermal, T_ref, verbose, bc, t0)
    nx, ny, nz = int(grid["nx"]), int(grid["ny"]), int(grid["nz"])
    dx = np.asarray(grid["dx"], dtype=float)
    dy = np.asarray(grid["dy"], dtype=float)
    dz = np.asarray(grid["dz"], dtype=float)
    if dx.ndim == 0:
        dx = np.full(nx, float(dx))
    if dy.ndim == 0:
        dy = np.full(ny, float(dy))
    if dz.ndim == 0:
        dz = np.full(nz, float(dz))
    solid = grid.get("solid")
    if solid is None:
        solid = np.ones((nx, ny, nz), dtype=bool)
    mat_names = thermal.get("mat_names")
    if mat_names is None:
        mat_names = material_field(case, grid)

    nnx, nny, nnz = nx + 1, ny + 1, nz + 1
    node_used = np.zeros((nnx, nny, nnz), dtype=bool)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if not solid[i, j, k]:
                    continue
                node_used[i:i + 2, j:j + 2, k:k + 2] = True

    fix, k_bc, triple = _fix_components(bc, node_used)
    dof_id = -np.ones((nnx, nny, nnz, 3), dtype=int)
    ndof = 0
    n_nodes = 0
    for i in range(nnx):
        for j in range(nny):
            for k in range(nnz):
                if not node_used[i, j, k]:
                    continue
                n_nodes += 1
                for c in range(3):
                    if not fix[i, j, k, c]:
                        dof_id[i, j, k, c] = ndof
                        ndof += 1
    if verbose:
        print("[elast] bc=%s  element=%s  nodes=%d  ndof=%d  solid cells=%d" % (
            bc, element, n_nodes, ndof, int(solid.sum())))
        if triple is not None:
            print("[elast] 3-2-1 A=%s B=%s C=%s  k=%d" % (
                triple[0], triple[1], triple[2], k_bc))
    if ndof == 0:
        raise ValueError("no free mechanical DOFs")

    ke_cache = {}

    def get_ke(dxi, dyj, dzk, E, nu):
        key = (round(dxi, 16), round(dyj, 16), round(dzk, 16), E, nu, element)
        hit = ke_cache.get(key)
        if hit is None:
            C = _iso_C(E, nu)
            Ke, wpack = _hex_Ke_fe(dxi, dyj, dzk, C, element=element)
            hit = (Ke, wpack, C)
            ke_cache[key] = hit
        return hit

    Tn = _nodal_T_hex8(T, solid, T_ref)
    rows, cols, data = [], [], []
    f = np.zeros(ndof, dtype=float)
    missing = set()
    n_el = 0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if not solid[i, j, k]:
                    continue
                name = mat_names[i, j, k]
                mat = case.material(name) if name else case.material()
                if mat.E is None:
                    missing.add(name or "?")
                    continue
                Ke, wpack, C = get_ke(dx[i], dy[j], dz[k], mat.E, mat.nu)
                T8 = np.empty(8, dtype=float)
                gdof = np.empty(24, dtype=int)
                for a, (di, dj, dk) in enumerate(_LOC):
                    T8[a] = Tn[i + di, j + dj, k + dk]
                    gdof[3 * a] = dof_id[i + di, j + dj, k + dk, 0]
                    gdof[3 * a + 1] = dof_id[i + di, j + dj, k + dk, 1]
                    gdof[3 * a + 2] = dof_id[i + di, j + dj, k + dk, 2]
                fe = _hex_fe_from_T8(
                    dx[i], dy[j], dz[k], C, T8, T_ref, mat.alpha_cte,
                    element=element, wpack=wpack)
                for a in range(24):
                    ga = gdof[a]
                    if ga < 0:
                        continue
                    f[ga] += fe[a]
                    for b in range(24):
                        gb = gdof[b]
                        if gb < 0:
                            continue
                        rows.append(ga)
                        cols.append(gb)
                        data.append(Ke[a, b])
                n_el += 1
    if missing:
        raise ValueError("materials missing E_Pa: %s" % sorted(missing))

    K = coo_matrix((data, (rows, cols)), shape=(ndof, ndof)).tocsr()
    t_asm = time.time()
    if verbose:
        print("[elast] assembled %d hexes  nnz=%d  %.1fs" % (
            n_el, K.nnz, t_asm - t0))

    u, mech_solver = _solve_cg(K, f, verbose=verbose)
    t_solve = time.time()

    ux = np.zeros((nnx, nny, nnz))
    uy = np.zeros((nnx, nny, nnz))
    uz = np.zeros((nnx, nny, nnz))
    for i in range(nnx):
        for j in range(nny):
            for k in range(nnz):
                if not node_used[i, j, k]:
                    continue
                dxid = dof_id[i, j, k, 0]
                dyid = dof_id[i, j, k, 1]
                dzid = dof_id[i, j, k, 2]
                if dxid >= 0:
                    ux[i, j, k] = u[dxid]
                if dyid >= 0:
                    uy[i, j, k] = u[dyid]
                if dzid >= 0:
                    uz[i, j, k] = u[dzid]

    ucx = np.zeros((nx, ny, nz))
    ucy = np.zeros((nx, ny, nz))
    ucz = np.zeros((nx, ny, nz))
    vm = np.zeros((nx, ny, nz))
    sxx = np.zeros((nx, ny, nz))
    syy = np.zeros((nx, ny, nz))
    szz = np.zeros((nx, ny, nz))
    sxy = np.zeros((nx, ny, nz))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if not solid[i, j, k]:
                    continue
                corners = [(i + di, j + dj, k + dk) for di, dj, dk in _LOC]
                ucx[i, j, k] = sum(ux[p] for p in corners) / 8.0
                ucy[i, j, k] = sum(uy[p] for p in corners) / 8.0
                ucz[i, j, k] = sum(uz[p] for p in corners) / 8.0
                name = mat_names[i, j, k]
                mat = case.material(name) if name else case.material()
                C = _iso_C(mat.E, mat.nu)
                T8 = np.empty(8, dtype=float)
                u24 = np.empty(24)
                for a, (di, dj, dk) in enumerate(_LOC):
                    T8[a] = Tn[i + di, j + dj, k + dk]
                    u24[3 * a] = ux[i + di, j + dj, k + dk]
                    u24[3 * a + 1] = uy[i + di, j + dj, k + dk]
                    u24[3 * a + 2] = uz[i + di, j + dj, k + dk]
                sv = _hex_stress_gauss(
                    u24, dx[i], dy[j], dz[k], C, T8, T_ref, mat.alpha_cte)
                sxx[i, j, k], syy[i, j, k], szz[i, j, k] = sv[0], sv[1], sv[2]
                sxy[i, j, k] = sv[3]
                vm[i, j, k] = _von_mises(sv)

    umag = np.sqrt(ucx ** 2 + ucy ** 2 + ucz ** 2)
    elapsed = time.time() - t0
    timing = {
        "assemble_s": t_asm - t0,
        "solve_s": t_solve - t_asm,
        "post_s": elapsed - (t_solve - t0),
        "total_s": elapsed,
        "ndof": int(ndof),
        "n_nodes": int(n_nodes),
        "n_hex": int(n_el),
        "n_solid": int(solid.sum()),
        "nnz": int(K.nnz),
        "solver": mech_solver,
        "bc": bc,
        "element": element,
    }
    if verbose:
        m = solid
        print("[elast] |u|_max=%.3e m  uz_max=%.3e m  vm_max=%.3e Pa  %.1fs" % (
            umag[m].max(), np.abs(ucz[m]).max(), vm[m].max(), elapsed))
        print("[elast] timing  assemble=%.2fs  solve=%.2fs (%s)  post=%.2fs  ndof=%d  nnz=%d" % (
            timing["assemble_s"], timing["solve_s"], mech_solver,
            timing["post_s"], ndof, K.nnz))
    return {
        "ux": ux, "uy": uy, "uz": uz,
        "ucx": ucx, "ucy": ucy, "ucz": ucz, "umag": umag,
        "vm": vm, "sxx": sxx, "syy": syy, "szz": szz, "sxy": sxy,
        "T_ref": T_ref, "elapsed_s": elapsed, "ndof": ndof,
        "thermal": thermal, "case": case, "grid": grid, "bc": bc,
        "timing": timing, "element": element,
    }
