#!/usr/bin/env python3
"""Six packages, traditional SA against the coarse/fine SA, for the paper.

Traditional: 64 grid for all 120 steps, illegal proposals dropped, the
40/25/20/15 move mix throughout.
New: 32 grid for the first 75% with illegal proposals re-drawn until legal and
a 10/60/20/10 mix, then 64 grid for the last 25% with an 80/5/5/10 mix.

Both start from the same L-BFGS layout and the same seed, and both final
layouts are re-scored on the 64 grid. One row per run goes to results.jsonl,
and each layout's temperature field is kept as a figure.
"""

from __future__ import annotations

import importlib.util
import json
import random
import shutil
import sys
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

font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

spec = importlib.util.spec_from_file_location("v3", SRC / "run_sa_mix_v3.py")
v3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v3)

RUNS = ROOT / "热仿真流水线" / "runs"
OUT = RUNS / "sa_paper_cases"
SEED = 22

# One HBM3 stack per package. A stack pays the background once per channel and
# then the traffic of every CPU and GPU driving it, measured back to back on
# HBM3_16Gb_8hi at 6400 Mbps: 0.088 W background, 1.035 W for the CPU task,
# 0.844 W for the GPU task, 16 channels.
CHAN_BG_W, CHAN_CPU_W, CHAN_GPU_W, N_CHAN = 0.088, 1.035, 0.844, 16
HBM_AREA = 11.0 * 11.0

# Each case keeps the chips of an existing L-BFGS layout, drops its CIM dies
# down to the same count as the GPUs, and adds one HBM3 stack. All three die
# classes land in the same area order: T2 and the CIM at 274 mm^2, Xeon at
# 422, the GPUs at 392 to 471, HBM at 121.
CASES = [
    ("C1_1T2_1GPU", 1, 1, RUNS / "big_balanced_2/lbfgs_init.json"),
    ("C2_2T2_2GPU", 2, 2, RUNS / "full_chain/pkg_2t2_2gpu/lbfgs_init.json"),
    ("C3_XEON_T2_2GPU", 2, 2, RUNS / "full_chain/pkg_1xeon_1t2_2gpu/lbfgs_init.json"),
    ("C4_2XEON_2GPU", 2, 2, RUNS / "full_chain/pkg_2xeon_2gpu/lbfgs_init.json"),
    ("C5_2XEON_3GPU", 2, 3, RUNS / "full_chain/pkg_2xeon_3gpu/lbfgs_init.json"),
    ("C6_3XEON_3GPU", 3, 3, RUNS / "full_chain/pkg_3xeon_3gpu/lbfgs_init.json"),
]


def load(path: Path, n_cpu: int, n_gpu: int):
    raw = [sa.Chip(c["name"], c["kind"], c["x_m"], c["y_m"], c["w_m"], c["h_m"],
                   c["power_W"], c["area_mm2"]) for c in json.loads(path.read_text())]
    cpus = [c for c in raw if c.kind == "CPU"][:n_cpu]
    gpus = [c for c in raw if c.kind == "GPU"][:n_gpu]
    assert len(cpus) == n_cpu and len(gpus) == n_gpu, path
    # jia at 64 macros: 274.0 mm^2 and 20.5 mW, the same area order as the
    # CPUs and GPUs. The CIM dies in the saved layouts are the 1-macro part.
    side_cim = (274.0 ** 0.5) * 1e-3
    cims = [sa.Chip(f"CIM_jia64_{i}", "CIM", 0.0, 0.0, side_cim, side_cim,
                    0.020516, 274.0) for i in range(n_gpu)]
    side = (HBM_AREA ** 0.5) * 1e-3
    # Background once, then each CPU's and each GPU's own traffic.
    hbm_w = (CHAN_BG_W + n_cpu * (CHAN_CPU_W - CHAN_BG_W)
             + n_gpu * (CHAN_GPU_W - CHAN_BG_W)) * N_CHAN
    hbm = sa.Chip("HBM3_stack", "HBM", 0.0, 0.0, side, side, hbm_w, HBM_AREA)
    chips = cpus + gpus + cims + [hbm]
    from geometry import packing_metrics
    # The saved layouts were packed around the 1-macro CIM dies. Swapping in
    # the 16.6 mm jia dies makes chips overlap, which would make every proposal
    # illegal. Push them apart until the spacing rules hold.
    sa.shift_to_origin(chips)
    for _ in range(400):
        moved = False
        for i in range(len(chips)):
            for j in range(i + 1, len(chips)):
                a, b = chips[i], chips[j]
                dx = max(0.0, max(a.x, b.x) - min(a.x + a.w, b.x + b.w))
                dy = max(0.0, max(a.y, b.y) - min(a.y + a.h, b.y + b.h))
                need = sa.GAP_GPU_M if a.kind == "GPU" and b.kind == "GPU" else sa.GAP_M
                if (dx * dx + dy * dy) ** 0.5 >= need:
                    continue
                b.x += a.w + need
                moved = True
        if not moved:
            break
    sa.shift_to_origin(chips)
    m = packing_metrics(sa.to_rects(chips))
    n_gpu = sum(1 for c in chips if c.kind == "GPU")
    w = m["width_m"] + sa.GAP_GPU_M * max(n_gpu, 1) + 0.008
    h = m["height_m"] + sa.GAP_GPU_M * max(n_gpu, 1) + 0.008
    return sa.spread_to_gap(chips, w, h), w, h


def search(base, w, h, work: Path, new: bool):
    rng = random.Random(SEED)
    if not new:
        lay, evals, acc, _, _ = v3.anneal(
            base, v3.STEPS, rng, 64, v3.UNIFORM, sa.T0, w, h, work, False)
        return lay, evals, acc
    n_coarse = round(v3.STEPS * v3.FRAC)
    lay, e1, a1, _, T = v3.anneal(
        base, n_coarse, rng, 32, v3.JUMP_COARSE, sa.T0, w, h, work, True)
    lay, e2, a2, _, _ = v3.anneal(
        lay, v3.STEPS - n_coarse, rng, 64, v3.STAGED_FINE, T, w, h, work, False)
    return lay, e1 + e2, a1 + a2


def field(chips, work: Path):
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_hotspot_inputs(sa.to_rects(chips), work, rows=64, cols=64)
    t = float(run_hotspot(sa_bin(), work)["Tmax_C"])
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


def sa_bin():
    return Path("/home/czh/任务仿真产生功耗文件工具/HotSpot-master/hotspot")


def draw(case, panels, png: Path) -> None:
    lo = min(g.min() for _, g, _, _ in panels)
    hi = max(g.max() for _, g, _, _ in panels)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4), constrained_layout=True)
    im = None
    for ax, (title, g, tmax, ext) in zip(axes, panels):
        im = ax.imshow(g, origin="lower", extent=ext, vmin=lo, vmax=hi, cmap="inferno")
        ax.set_title(f"{title}\n峰值 {tmax:.2f} °C", fontsize=12)
        ax.set_xlabel("x / mm")
        ax.set_ylabel("y / mm")
    fig.colorbar(im, ax=axes, fraction=0.046, pad=0.02).set_label("温度 / °C")
    fig.suptitle(case, fontsize=14)
    fig.savefig(png, dpi=140)
    plt.close(fig)


def main() -> None:  # cases: name, n_cpu, n_gpu, path

    import time
    OUT.mkdir(parents=True, exist_ok=True)
    log = OUT / "results.jsonl"
    for i, (name, n_cpu, n_gpu, path) in enumerate(CASES, 1):
        base, w, h = load(path, n_cpu, n_gpu)
        g0, t0, e0 = field(base, OUT / name / "field_init")
        t_a = time.time()
        lay_old, ev_old, ac_old = search(base, w, h, OUT / name / "work_old", False)
        sec_old = time.time() - t_a
        t_b = time.time()
        lay_new, ev_new, ac_new = search(base, w, h, OUT / name / "work_new", True)
        sec_new = time.time() - t_b
        g1, t_old, e1 = field(lay_old, OUT / name / "field_old")
        g2, t_new, e2 = field(lay_new, OUT / name / "field_new")
        draw(name, [("搜索前", g0, t0, e0), ("传统 SA", g1, t_old, e1),
                    ("粗细网格 SA", g2, t_new, e2)], OUT / name / "compare.png")
        row = {
            "case": name, "chips": len(base),
            "power_W": round(sum(c.power_W for c in base), 1),
            "init_C": round(t0, 3),
            "traditional_C": round(t_old, 3), "traditional_s": round(sec_old, 1),
            "traditional_evals": ev_old,
            "new_C": round(t_new, 3), "new_s": round(sec_new, 1), "new_evals": ev_new,
        }
        with log.open("a") as fh:
            fh.write(json.dumps(row, ensure_ascii=False) + "\n")
        print(f"[{i}/6] {name}: 起点 {t0:.2f}  传统 {t_old:.2f} ({sec_old:.0f}s)  "
              f"最新 {t_new:.2f} ({sec_new:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
