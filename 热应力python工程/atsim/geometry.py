# -*- coding: utf-8 -*-
"""Case JSON loader and structured grid helpers."""

from __future__ import print_function

import json
import os

import numpy as np

from .materials import material_from_dict


class Block(object):
    def __init__(self, d):
        self.name = d["name"]
        self.x0 = float(d["x0_m"])
        self.y0 = float(d["y0_m"])
        self.w = float(d["w_m"])
        self.h = float(d["h_m"])
        self.power = float(d["power_W"])
        self.material = d.get("material", "silicon")
        self.z0 = float(d["z0_m"]) if d.get("z0_m") is not None else None
        self.z1 = float(d["z1_m"]) if d.get("z1_m") is not None else None

    @property
    def x1(self):
        return self.x0 + self.w

    @property
    def y1(self):
        return self.y0 + self.h

    @property
    def area(self):
        return self.w * self.h

    def contains_xy(self, x, y):
        return (x >= self.x0) & (x < self.x1) & (y >= self.y0) & (y < self.y1)


class Probe(object):
    def __init__(self, d):
        self.name = d["name"]
        self.x = float(d["x_m"])
        self.y = float(d["y_m"])
        self.z = float(d["z_m"])


class CaseConfig(object):
    def __init__(self, d, path=None):
        self.raw = d
        self.path = path
        self.name = d.get("name", "unnamed")
        self.T_amb = float(d["ambient_K"])
        conv = d["convection"]
        self.h_top = float(conv["h_W_m2K"])
        self.conv_face = conv.get("face", "top")
        grid = d.get("grid", {})
        self.nx = int(grid.get("nx", 64))
        self.ny = int(grid.get("ny", 64))
        self.nz = int(grid.get("nz", 6))
        solver = d.get("solver", {})
        self.nonlinear = bool(solver.get("nonlinear", False))
        self.leakage = bool(solver.get("leakage", False))
        self.use_kirchhoff = bool(solver.get("kirchhoff", False))
        self.tol = float(solver.get("tolerance_K", 0.05))
        self.max_iter = int(solver.get("max_iter", 20))
        self.default_mat = d.get("material_default", "silicon")
        self.materials = {
            name: material_from_dict(name, spec)
            for name, spec in d.get("materials", {}).items()
        }
        self.blocks = [Block(b) for b in d["blocks"]]
        self.probes = [Probe(p) for p in d.get("probes", [])]
        self.hotspot = d.get("hotspot", {})
        if "geometry" in d:
            g = d["geometry"]
            self.Lx = float(g["width_m"])
            self.Ly = float(g["height_m"])
            self.Lz = float(g["thickness_m"])
        else:
            # package stack: domain size filled in by attach_package
            self.Lx = self.Ly = self.Lz = 0.0
        from .package_grid import attach_package
        attach_package(self)

    @property
    def total_power(self):
        return sum(b.power for b in self.blocks)

    def material(self, name=None):
        name = name or self.default_mat
        return self.materials[name]


def load_case(path):
    try:
        with open(path, "r", encoding="utf-8") as f:
            d = json.load(f)
    except TypeError:
        # Python 2 fallback
        with open(path, "r") as f:
            d = json.loads(f.read().decode("utf-8"))
    if d.get("chips"):
        from .layout import expand_chips_layout
        d = expand_chips_layout(d)
    return CaseConfig(d, path=os.path.abspath(path))


def build_grid(case):
    """Cell-centered Cartesian grid. Package cases use a solid mask + voids."""
    if getattr(case, "has_package", False):
        from .package_grid import build_package_grid
        return build_package_grid(case)
    nx, ny, nz = case.nx, case.ny, case.nz
    dx = case.Lx / nx
    dy = case.Ly / ny
    dz = case.Lz / nz
    xc = (np.arange(nx) + 0.5) * dx
    yc = (np.arange(ny) + 0.5) * dy
    zc = (np.arange(nz) + 0.5) * dz
    X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")
    solid = np.ones((nx, ny, nz), dtype=bool)
    mat_names = np.empty((nx, ny, nz), dtype=object)
    mat_names[:] = case.default_mat
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
        "X": X,
        "Y": Y,
        "Z": Z,
        "solid": solid,
        "mat_names": mat_names,
    }


def power_density_field(case, grid):
    """Volumetric power density [W/m^3] on the cell-centered grid."""
    if getattr(case, "has_package", False):
        from .package_grid import package_power_density
        return package_power_density(case, grid)
    P = np.zeros((case.nx, case.ny, case.nz), dtype=float)
    X = grid["X"][:, :, 0]
    Y = grid["Y"][:, :, 0]
    assigned = np.zeros((case.nx, case.ny), dtype=bool)
    for blk in case.blocks:
        mask = blk.contains_xy(X, Y)
        if not np.any(mask):
            raise ValueError("block %s covers no cells" % blk.name)
        q = blk.power / (blk.area * case.Lz)
        P[mask, :] = q
        assigned |= mask
    if not np.all(assigned):
        n_empty = int(np.size(assigned) - np.count_nonzero(assigned))
        print("[warn] %d cells have no power block (treated as 0 W)" % n_empty)
    return P


def material_field(case, grid):
    if "mat_names" in grid and grid["mat_names"] is not None:
        names = grid["mat_names"]
        if getattr(case, "has_package", False):
            return names
        X = grid["X"][:, :, 0]
        Y = grid["Y"][:, :, 0]
        for blk in case.blocks:
            mask = blk.contains_xy(X, Y)
            names[mask, :] = blk.material
        return names
    names = np.empty((case.nx, case.ny, case.nz), dtype=object)
    names[:] = case.default_mat
    X = grid["X"][:, :, 0]
    Y = grid["Y"][:, :, 0]
    for blk in case.blocks:
        mask = blk.contains_xy(X, Y)
        names[mask, :] = blk.material
    return names
