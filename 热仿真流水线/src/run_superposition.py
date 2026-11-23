#!/usr/bin/env python3
"""Can one chip be edited by subtracting its field instead of re-solving?

Steady heat conduction with fixed conductivity is linear, so the package field
should equal the sum of the fields of each chip powered alone. This checks that
on one real package, then times the two ways to edit a single chip:

- full: re-solve the whole package.
- delta: solve once with only the edited chip powered, subtract that field from
  the cached full field, add the field of the chip in its new state.

The grid, the spreader and the sink stay exactly as the full solve uses them,
so the two fields live on the same cells and can be added pointwise.
"""

from __future__ import annotations

import shutil
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

HS = Path("/home/czh/任务仿真产生功耗文件工具/HotSpot-master/hotspot")
OUT = ROOT / "热仿真流水线" / "runs" / "superposition"
ROWS = 64


def grid(work: Path) -> np.ndarray:
    vals = []
    for line in (work / "grid.steady").read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("Layer"):
            continue
        vals.append(float(line.split()[-1]))
        if len(vals) == ROWS * ROWS:
            break
    return np.array(vals).reshape(ROWS, ROWS)


def solve(chips, work: Path) -> np.ndarray:
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_hotspot_inputs(sa.to_rects(chips), work, rows=ROWS, cols=ROWS)
    run_hotspot(HS, work)
    return grid(work)


def load():
    import importlib.util
    spec = importlib.util.spec_from_file_location("p", SRC / "run_sa_paper_cases.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    name, n_cpu, n_gpu, path = m.CASES[3]  # C4, 7 chips, mid size
    chips, w, h = m.load(path, n_cpu, n_gpu)
    return name, chips


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    name, chips = load()
    print(f"{name}: {len(chips)} 芯片", flush=True)

    t0 = time.time()
    full = solve(chips, OUT / "full")
    t_full = time.time() - t0

    # One field per chip, every other chip powered at zero. Same geometry, so
    # the spreader and sink match the full solve.
    parts = []
    t_parts = 0.0
    for i, c in enumerate(chips):
        solo = [sa.Chip(x.name, x.kind, x.x, x.y, x.w, x.h,
                        c.power_W if x is c else 0.0, x.area_mm2) for x in chips]
        t0 = time.time()
        parts.append(solve(solo, OUT / f"solo_{i}"))
        t_parts += time.time() - t0
    summed = np.sum(parts, axis=0) - (len(chips) - 1) * 318.15

    err = summed - full
    print(f"全场峰值        {full.max() - 273.15:.3f} °C")
    print(f"叠加峰值        {summed.max() - 273.15:.3f} °C")
    print(f"逐格最大误差    {np.abs(err).max():.4f} K")
    print(f"逐格平均误差    {np.abs(err).mean():.4f} K")
    print(f"峰值位置误差    {abs(summed.max() - full.max()):.4f} K")

    # Edit one chip: move it by 3 mm and raise its power 20%.
    i = next(k for k, c in enumerate(chips) if c.kind == "GPU")
    edited = [sa.Chip(c.name, c.kind, c.x, c.y, c.w, c.h, c.power_W, c.area_mm2)
              for c in chips]
    edited[i].x += 0.003
    edited[i].power_W *= 1.2

    t0 = time.time()
    ref = solve(edited, OUT / "edited_full")
    t_ref = time.time() - t0

    # Two solves, both on the new geometry so the cells line up. The first
    # carries every chip except the edited one, the second carries only it.
    t0 = time.time()
    others = [sa.Chip(c.name, c.kind, c.x, c.y, c.w, c.h,
                      0.0 if c is edited[i] else chips[k].power_W, c.area_mm2)
              for k, c in enumerate(edited)]
    solo_new = [sa.Chip(c.name, c.kind, c.x, c.y, c.w, c.h,
                        edited[i].power_W if c is edited[i] else 0.0, c.area_mm2)
                for c in edited]
    delta = solve(others, OUT / "edit_others") + solve(solo_new, OUT / "edit_new") - 318.15
    t_delta = time.time() - t0

    derr = delta - ref
    print()
    print(f"改动 {chips[i].name}: 平移 3 mm, 功耗 +20%")
    print(f"完整重解峰值    {ref.max() - 273.15:.3f} °C   用时 {t_ref:.2f} s")
    print(f"场替换峰值      {delta.max() - 273.15:.3f} °C   用时 {t_delta:.2f} s")
    print(f"逐格最大误差    {np.abs(derr).max():.4f} K")
    print(f"峰值误差        {abs(delta.max() - ref.max()):.4f} K")
    print()
    print(f"单次全场求解    {t_full:.2f} s")
    print(f"{len(chips)} 颗芯片各解一次 {t_parts:.2f} s")
    print(f"一次改动: 完整重解 {t_ref:.2f} s, 场替换 {t_delta:.2f} s, "
          f"加速 {t_ref / t_delta:.2f} 倍")


if __name__ == "__main__":
    main()
