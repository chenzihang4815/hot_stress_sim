#!/usr/bin/env python3
"""Cheap ATSim heat+stress proxy for SA (silicon fill, no plots).

Returns chip-center von Mises mean/max and interposer warpage PV.
Default fill_die=True so the mechanical layer matches HotSpot's continuous
silicon chip layer. Do not use package-max von Mises (edge singularity).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
ATSIM = ROOT / "热应力python工程_最终"
sys.path.insert(0, str(ATSIM))

from atsim.elasticity import solve_elasticity  # noqa: E402
from atsim.fvm3d import solve_fvm3d  # noqa: E402
from atsim.geometry import load_case  # noqa: E402
from atsim.warpage import package_warpage  # noqa: E402


def _package_dict(fill_die: bool, design) -> dict:
    """热边界和 hotspot_io 不传 package 时的那一套对齐。

    芯片框 72 mm，spreader 和热沉都是 96 mm，r_convec 按 60 mm 上的 0.1
    缩到 96 mm。以前这里写死 0.1，热沉又跟着芯片外框加边，对流面积和
    HotSpot 不是同一个。基板和中介层仍留着，应力要用；HotSpot 的网格模型
    没有这两层。DEFAULT_PACKAGE 不动。
    """
    pkg = {
        "tim_um": 20, "t_spreader_mm": 1.0, "t_sink_mm": 6.9,
        "s_spreader_mm": 96.0, "s_sink_mm": 96.0, "frame_mm": 72.0,
        "tim_cover": "frame",
        "t_substrate_mm": 1.0, "t_interposer_um": 100.0,
        "r_convec_K_per_W": 0.0390625, "ambient_C": 45,
        "pitch_die_mm": 2.0, "pitch_coarse_mm": 8.0,
        "nz_die": 3, "nz_spreader": 3, "nz_sink": 4, "nz_tim": 1,
        "nz_substrate": 3, "nz_interposer": 2,
        "fill_die": bool(fill_die), "fill_xy": "chips", "fill_material": "silicon",
    }
    if design is not None:
        p, m = design.package, design.materials
        pkg.update(tim_um=p.tim_um, t_spreader_mm=p.t_spreader_mm, t_sink_mm=p.t_sink_mm,
                   s_spreader_mm=p.spreader_mm, s_sink_mm=p.sink_mm, frame_mm=p.die_mm,
                   t_substrate_mm=m.t_substrate_mm, t_interposer_um=m.t_interposer_um,
                   r_convec_K_per_W=p.r_convec_K_per_W, ambient_C=p.ambient_C,
                   fill_material=m.fill)
    return pkg


def write_proxy_json(chips, path: Path, title: str, fill_die: bool = True, package_overrides=None, design=None) -> Path:
    payload = {
        "name": path.stem,
        "description": title,
        "chips": [
            {
                "name": c.name,
                "x_mm": c.x * 1e3,
                "y_mm": c.y * 1e3,
                "w_mm": c.w * 1e3,
                "h_mm": c.h * 1e3,
                "t_um": 150.0,
                "power_W": c.power_W,
            }
            for c in chips
        ],
        "package": _package_dict(fill_die, design),
    }
    if package_overrides:
        payload["package"].update(package_overrides)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


MIN_NORMAL_CIM_AREA_MM2 = 1.0


def _is_normal_chip(name: str, area_mm2: float) -> bool:
    """CPU/GPU always; CIM only if large enough to avoid free-edge spikes."""
    n = name.upper()
    if n.startswith("CPU") or n.startswith("GPU"):
        return True
    if n.startswith("CIM"):
        return area_mm2 >= MIN_NORMAL_CIM_AREA_MM2
    return area_mm2 >= MIN_NORMAL_CIM_AREA_MM2


def _chip_interior_vm(el) -> list[dict]:
    """Per-chiplet von Mises at die mid-plane, cells whose centers lie in the chip."""
    from atsim.visualize import die_mid_k

    case = el["case"]
    grid = el["grid"]
    k = die_mid_k({"grid": grid, "case": case, "T": el["thermal"]["T"]})
    xc, yc = grid["xc"], grid["yc"]
    vm = el["vm"][:, :, k]
    solid = grid["solid"][:, :, k]
    ox = getattr(case, "chip_ox", 0.0)
    oy = getattr(case, "chip_oy", 0.0)
    rows = []
    for b in case.blocks:
        area_mm2 = float(b.w * b.h * 1e6)
        x0, y0, w, h = ox + b.x0, oy + b.y0, b.w, b.h
        ix = np.where((xc >= x0) & (xc <= x0 + w))[0]
        jy = np.where((yc >= y0) & (yc <= y0 + h))[0]
        vals = np.array([])
        if ix.size and jy.size:
            sub = vm[np.ix_(ix, jy)]
            msk = solid[np.ix_(ix, jy)]
            vals = sub[msk] * 1e-6
        rec = {
            "name": b.name,
            "area_mm2": area_mm2,
            "normal": _is_normal_chip(b.name, area_mm2),
            "n_cells": int(vals.size),
            "vm_max_MPa": float(vals.max()) if vals.size else None,
            "vm_p95_MPa": float(np.percentile(vals, 95)) if vals.size else None,
            "vm_mean_MPa": float(vals.mean()) if vals.size else None,
        }
        rows.append(rec)
    return rows


def _chip_center_rows(el) -> list[dict]:
    skip = ("Sink", "Spreader", "Interposer", "Substrate", "TIM", "Chip_center")
    rows = []
    grid = el["grid"]
    xc, yc, zc = grid["xc"], grid["yc"], grid["zc"]
    case = el["case"]
    for p in case.probes:
        if not p.name.endswith("_center"):
            continue
        if p.name.startswith(skip) or p.name.replace("_center", "") in skip:
            continue
        i = int(np.argmin(np.abs(xc - p.x)))
        j = int(np.argmin(np.abs(yc - p.y)))
        k = int(np.argmin(np.abs(zc - p.z)))
        rows.append(
            {
                "name": p.name.replace("_center", ""),
                "T_C": float(el["thermal"]["T"][i, j, k] - 273.15),
                "vm_MPa": float(el["vm"][i, j, k] * 1e-6),
                "uz_um": float(el["ucz"][i, j, k] * 1e6),
            }
        )
    return rows


def eval_atsim_proxy(json_path: Path, verbose: bool = False) -> dict:
    t0 = time.time()
    case = load_case(str(json_path))
    th = solve_fvm3d(case, verbose=verbose)
    el = solve_elasticity(th, T_ref=case.T_amb, verbose=verbose)
    solid = el["grid"]["solid"]
    warp = package_warpage(el)
    centers = _chip_center_rows(el)
    interiors = _chip_interior_vm(el)
    vms = [c["vm_MPa"] for c in centers]
    normal = [r for r in interiors if r.get("normal") and r.get("vm_max_MPa") is not None]
    worst = max(normal, key=lambda r: r["vm_max_MPa"]) if normal else None
    out = {
        "json": str(json_path),
        "fill_die": bool((case.raw.get("package") or {}).get("fill_die", False)),
        "Tmax_C": float(th["T"][solid].max() - 273.15),
        "vm_center_mean_MPa": float(np.mean(vms)) if vms else None,
        "vm_center_max_MPa": float(np.max(vms)) if vms else None,
        "vm_normal_chip_max_MPa": None if worst is None else float(worst["vm_max_MPa"]),
        "vm_normal_chip_name": None if worst is None else worst["name"],
        "vm_normal_chip_p95_MPa": None if worst is None else worst.get("vm_p95_MPa"),
        "chip_interiors": interiors,
        "w_pv_um": float(warp["w_pv_um"]),
        "u_max_um": float(el["umag"][solid].max() * 1e6),
        "vm_package_max_MPa": float(el["vm"][solid].max() * 1e-6),
        "chip_centers": centers,
        "heat_s": float(th["elapsed_s"]),
        "elast_s": float(el["elapsed_s"]),
        "elapsed_s": time.time() - t0,
        "ndof": el["ndof"],
        "element": el.get("element") or (el.get("timing") or {}).get("element"),
        "solver": (el.get("timing") or {}).get("solver"),
    }
    return out
