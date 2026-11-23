#!/usr/bin/env python3
"""Render the three layouts behind the coarse/fine comparison.

Left: the spread starting layout. Middle: the current SA, 64 grid all the way,
uniform move mix. Right: coarse-then-fine, 32 grid for the first 75% with
illegal proposals re-drawn and a jump-heavy mix, 64 grid for the last 25%.

Seed 22 for both searches. Each panel is that layout's own 64-grid HotSpot
field, all three on one colour scale.
"""

from __future__ import annotations

import importlib.util
import json
import random
import shutil
import sys
from copy import deepcopy
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

font_manager.fontManager.addfont(
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

spec = importlib.util.spec_from_file_location("v3", SRC / "run_sa_mix_v3.py")
v3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v3)

HS = Path("/home/czh/任务仿真产生功耗文件工具/HotSpot-master/hotspot")
OUT = ROOT / "热仿真流水线" / "runs" / "sa_compare_fig"
SEED = 22


def field(chips, work: Path):
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_hotspot_inputs(sa.to_rects(chips), work, rows=64, cols=64)
    t = float(run_hotspot(HS, work)["Tmax_C"])
    vals = []
    for line in (work / "grid.steady").read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("Layer"):
            continue
        vals.append(float(line.split()[-1]))
        if len(vals) == 4096:
            break
    g = np.array(vals).reshape(64, 64)[::-1] - 273.15
    xs = [c.x for c in chips] + [c.x + c.w for c in chips]
    ys = [c.y for c in chips] + [c.y + c.h for c in chips]
    return g, t, (min(xs) * 1e3, max(xs) * 1e3, min(ys) * 1e3, max(ys) * 1e3)


def uniform_search(base, w, h, work: Path):
    """The current SA: 64 grid throughout, no redraw, the 40/25/20/15 mix."""
    rng = random.Random(SEED)
    lay, *_ = v3.anneal(base, v3.STEPS, rng, 64, v3.UNIFORM, sa.T0, w, h, work, False)
    return lay


def coarse_fine_search(base, w, h, work: Path):
    """32 grid for 75% with redraw and the jump-heavy mix, then 64 grid."""
    rng = random.Random(SEED)
    n_coarse = round(v3.STEPS * v3.FRAC)
    lay, *_ = v3.anneal(base, n_coarse, rng, 32, v3.JUMP_COARSE, sa.T0, w, h, work, True)
    lay, *_ = v3.anneal(lay, v3.STEPS - n_coarse, rng, 64, v3.STAGED_FINE,
                        sa.T0 * sa.ALPHA ** n_coarse, w, h, work, False)
    return lay


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    base, w, h = v3.starting_layout()
    panels = [
        ("搜索前", base),
        ("传统 SA\n全程 64 网格", uniform_search(base, w, h, OUT / "work_uniform")),
        ("粗细网格 SA\n32 网格前 75%，64 网格后 25%",
         coarse_fine_search(base, w, h, OUT / "work_coarse")),
    ]
    rendered = [(title, *field(lay, OUT / f"field_{i}")) for i, (title, lay) in enumerate(panels)]

    lo = min(g.min() for _, g, _, _ in rendered)
    hi = max(g.max() for _, g, _, _ in rendered)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.6), constrained_layout=True)
    im = None
    for ax, (title, g, tmax, ext) in zip(axes, rendered):
        im = ax.imshow(g, origin="lower", extent=ext, vmin=lo, vmax=hi, cmap="inferno")
        ax.set_title(f"{title}\n峰值 {tmax:.2f} °C", fontsize=12)
        ax.set_xlabel("x / mm")
        ax.set_ylabel("y / mm")
    fig.colorbar(im, ax=axes, fraction=0.046, pad=0.02).set_label("温度 / °C")
    fig.suptitle("同一起点、同一随机种子下的布局温度场", fontsize=14)
    png = OUT / "sa_before_after.png"
    fig.savefig(png, dpi=140)
    print("WROTE", png)
    for title, g, tmax, ext in rendered:
        print(f"{title.splitlines()[0]:10s} 峰值 {tmax:.3f} °C  场范围 {g.min():.2f}~{g.max():.2f}")


if __name__ == "__main__":
    main()
