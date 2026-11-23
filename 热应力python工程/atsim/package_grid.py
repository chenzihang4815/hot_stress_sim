# -*- coding: utf-8 -*-
"""Package-stack grid (Case 2): die + TIM + spreader + sink with voids."""

from __future__ import print_function

import numpy as np


class Layer(object):
    def __init__(self, d):
        self.name = d["name"]
        self.material = d["material"]
        self.x0 = float(d["x0_m"])
        self.y0 = float(d["y0_m"])
        self.w = float(d["w_m"])
        self.h = float(d["h_m"])
        self.z0 = float(d["z0_m"])
        self.z1 = float(d["z1_m"])
        self.nz = int(d.get("nz", 2))
        self.has_power = bool(d.get("has_power", False))

    @property
    def thickness(self):
        return self.z1 - self.z0

    @property
    def x1(self):
        return self.x0 + self.w

    @property
    def y1(self):
        return self.y0 + self.h


def attach_package(case):
    """Populate CaseConfig with stack / domain fields if present in JSON."""
    d = case.raw
    origin = d.get("chip_origin_m", [0.0, 0.0])
    case.chip_ox = float(origin[0])
    case.chip_oy = float(origin[1])
    case.conv_face = d.get("convection", {}).get("face", "top")

    if "stack" not in d:
        case.has_package = False
        case.chip_Lx = case.Lx
        case.chip_Ly = case.Ly
        case.chip_Lz = case.Lz
        case.stack = [
            Layer(
                {
                    "name": "die",
                    "material": case.default_mat,
                    "x0_m": 0.0,
                    "y0_m": 0.0,
                    "w_m": case.Lx,
                    "h_m": case.Ly,
                    "z0_m": 0.0,
                    "z1_m": case.Lz,
                    "nz": case.nz,
                    "has_power": True,
                }
            )
        ]
        case.x_regions = None
        case.y_regions = None
        return case

    case.has_package = True
    case.stack = [Layer(s) for s in d["stack"]]
    dom = d["domain"]
    case.Lx = float(dom["width_m"])
    case.Ly = float(dom["height_m"])
    case.Lz = max(L.z1 for L in case.stack)
    die = None
    for L in case.stack:
        if L.has_power or L.name == "die":
            die = L
            break
    if die is None:
        die = case.stack[0]
    dies = [L for L in case.stack if L.has_power]
    if not dies:
        dies = [die]
    case.chip_Lx = max(L.x0 + L.w for L in dies) - min(L.x0 for L in dies)
    case.chip_Ly = max(L.y0 + L.h for L in dies) - min(L.y0 for L in dies)
    case.chip_Lz = max(L.thickness for L in dies)
    grid = d.get("grid", {})
    case.x_regions = grid.get("x_regions")
    case.y_regions = grid.get("y_regions")
    # nz is set when the grid is built (parallel dies share z, do not sum)
    return case


def die_layer(case):
    for L in case.stack:
        if L.has_power or L.name == "die":
            return L
    return case.stack[0]


def _edges_from_regions(regions, k0, k1, kn):
    edges = [float(regions[0][k0])]
    for r in regions:
        n = int(r[kn])
        span = float(r[k1]) - float(r[k0])
        dx = span / float(n)
        for _ in range(n):
            edges.append(edges[-1] + dx)
    return np.array(edges, dtype=float)


def _uniform_edges(length, n):
    return np.linspace(0.0, length, n + 1)


def _unique_sorted(vals, tol=1e-15):
    s = sorted(vals)
    out = [s[0]]
    for v in s[1:]:
        if abs(v - out[-1]) > tol:
            out.append(v)
    return out


def _z_edges_from_stack(stack):
    """z mesh from unique interfaces. Side-by-side dies share a z-band."""
    cuts = []
    for L in stack:
        cuts.extend([L.z0, L.z1])
    cuts = _unique_sorted(cuts)
    edges = [cuts[0]]
    for a, b in zip(cuts[:-1], cuts[1:]):
        span = b - a
        if span <= 1e-16:
            continue
        nz = 1
        for L in stack:
            if L.z0 < b - 1e-16 and L.z1 > a + 1e-16 and L.thickness > 0:
                n_here = max(1, int(round(L.nz * span / L.thickness)))
                if n_here > nz:
                    nz = n_here
        nz = max(1, min(int(nz), 8))
        dz = span / float(nz)
        for _ in range(nz):
            edges.append(edges[-1] + dz)
    return np.array(edges, dtype=float)


def build_package_grid(case):
    """Cell-centered Cartesian grid, possibly non-uniform, with a solid mask."""
    if getattr(case, "x_regions", None):
        x_edges = _edges_from_regions(case.x_regions, "x0_m", "x1_m", "n")
    else:
        x_edges = _uniform_edges(case.Lx, case.nx)
    if getattr(case, "y_regions", None):
        y_edges = _edges_from_regions(case.y_regions, "y0_m", "y1_m", "n")
    else:
        y_edges = _uniform_edges(case.Ly, case.ny)

    z_edges = _z_edges_from_stack(case.stack)

    dx = np.diff(x_edges)
    dy = np.diff(y_edges)
    dz = np.diff(z_edges)
    xc = 0.5 * (x_edges[:-1] + x_edges[1:])
    yc = 0.5 * (y_edges[:-1] + y_edges[1:])
    zc = 0.5 * (z_edges[:-1] + z_edges[1:])
    nx, ny, nz = len(xc), len(yc), len(zc)
    case.nx, case.ny, case.nz = nx, ny, nz

    X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")
    solid = np.zeros((nx, ny, nz), dtype=bool)
    mat_names = np.empty((nx, ny, nz), dtype=object)
    mat_names[:] = None
    layer_id = np.full((nx, ny, nz), -1, dtype=int)

    nstack = len(case.stack)
    for li, L in enumerate(case.stack):
        mx = (X >= L.x0 - 1e-16) & (X <= L.x0 + L.w + 1e-16)
        my = (Y >= L.y0 - 1e-16) & (Y <= L.y0 + L.h + 1e-16)
        if li == nstack - 1:
            mz = (zc[None, None, :] >= L.z0 - 1e-16) & (zc[None, None, :] <= L.z1 + 1e-16)
        else:
            mz = (zc[None, None, :] >= L.z0 - 1e-16) & (zc[None, None, :] < L.z1 - 1e-16)
        mask = mx & my & mz
        solid[mask] = True
        mat_names[mask] = L.material
        layer_id[mask] = li

    return {
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "dx": dx,
        "dy": dy,
        "dz": dz,
        "xc": xc,
        "yc": yc,
        "zc": zc,
        "x_edges": x_edges,
        "y_edges": y_edges,
        "z_edges": z_edges,
        "X": X,
        "Y": Y,
        "Z": Z,
        "solid": solid,
        "mat_names": mat_names,
        "layer_id": layer_id,
    }


def package_power_density(case, grid):
    nx, ny, nz = grid["nx"], grid["ny"], grid["nz"]
    P = np.zeros((nx, ny, nz), dtype=float)
    X = grid["X"][:, :, 0]
    Y = grid["Y"][:, :, 0]
    zc = grid["zc"]
    ox, oy = case.chip_ox, case.chip_oy
    solid = grid["solid"]
    for blk in case.blocks:
        mask_xy = blk.contains_xy(X - ox, Y - oy)
        if not np.any(mask_xy):
            raise ValueError("block %s covers no cells" % blk.name)
        z0 = blk.z0 if blk.z0 is not None else 0.0
        z1 = blk.z1 if blk.z1 is not None else die_layer(case).z1
        thick = max(z1 - z0, 1e-16)
        in_z = (zc >= z0 - 1e-16) & (zc < z1 - 1e-16)
        if not np.any(in_z):
            in_z = (zc >= z0 - 1e-16) & (zc <= z1 + 1e-16)
        q = blk.power / (blk.area * thick)
        for k in np.where(in_z)[0]:
            sl = P[:, :, k]
            sl[mask_xy & solid[:, :, k]] = q
            P[:, :, k] = sl
    P[~solid] = 0.0
    return P
