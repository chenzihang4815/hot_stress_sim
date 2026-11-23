#!/usr/bin/env python3
"""布局维度上的温度-应力帕累托。

5 芯粒、环氧填充、芯片固定为 n5_epoxy 的初始选型，不换芯片、不调频率。
先用 Python 场扫 400 个随机布局，只看峰值温度和线长，挑出散布开的 18 个，
再对这 18 个各跑一次 ATSim，温度和应力都取 ATSim 的数。

输出写到本目录的 runs/，不碰 热仿真流水线/runs，也不碰 闭环全流程 的任何目录。
"""

from __future__ import annotations

import json
import random
import sys
from copy import deepcopy
from pathlib import Path

SRC = Path("/home/czh/任务仿真产生功耗文件工具/热仿真流水线/src")
sys.path.insert(0, str(SRC))

import numpy as np
from scipy.sparse.linalg import splu

import run_closed_loop_multidie as cl
import run_closed_loop_suite as suite
import run_fill_superposition as hs
import run_hw_area_field as hw
import run_sa_hw_smooth as sm
import run_wirelength as wire

HERE = Path("/home/czh/任务仿真产生功耗文件工具/论文证据/帕累托_布局")
OUT = HERE / "runs"
SEED = 20261117
N_SCAN = 400
N_KEEP = 18
FILL = "epoxy"


def scatter(dies: list[dict], rng: random.Random, tries: int = 4000) -> list[dict] | None:
    """把每颗芯粒随机放到不重叠的位置，40% 的概率顺便旋转。"""
    for _ in range(tries):
        trial = deepcopy(dies)
        for die in trial:
            if rng.random() < 0.4:
                die["w"], die["h"] = die["h"], die["w"]
            die["j0"] = rng.randrange(hs.N - die["w"] + 1)
            die["ib0"] = rng.randrange(hs.N - die["h"] + 1)
        if sm.legal(trial):
            return trial
    return None


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "layouts").mkdir(exist_ok=True)
    rng = random.Random(SEED)

    cat = cl.catalog()
    base = suite.members(cat, "n5")
    saved = json.loads(
        (SRC.parent / "runs").joinpath(
            "../../论文证据/闭环全流程_只收更优/runs/n5_epoxy/layouts/init.json"
        ).resolve().read_text(encoding="utf-8")
    )
    # members() 给出的是带目录信息的芯粒，位置用存档的初始布局，保证芯片选型一致。
    by_slot = {d["slot"]: d for d in saved}
    for die in base:
        src = by_slot[die["slot"]]
        die.update({k: src[k] for k in ("name", "part", "j0", "ib0", "w", "h") if k in src})
    cl.apply_power(base)
    matrix = wire.assumed_matrix(base)

    k_si = hw.paint(base, hs.K_SI)
    mat_si, g_si = hs.build_operator(k_si)
    fac_si = splu(mat_si.tocsc())

    pool = [("init", base)]
    seen = set()
    for i in range(N_SCAN):
        trial = scatter(base, rng)
        if trial is None:
            continue
        key = tuple(sorted((d["slot"], d["j0"], d["ib0"], d["w"], d["h"]) for d in trial))
        if key in seen:
            continue
        seen.add(key)
        pool.append((f"r{i:03d}", trial))
    print(f"pool {len(pool)}", flush=True)

    rows = []
    for tag, dies in pool:
        cl.apply_power(dies)
        temp = hw.solve(k_si, hw.power_map(dies), fac_si, g_si)[0]
        ev = cl.metrics(temp)
        wl = wire.wirelength_cells(dies, matrix, cl.CW_MM)
        rows.append({
            "tag": tag,
            "T_py_C": ev["Tmax_C"],
            "grad": ev["grad_max"],
            "wl_mm": wl,
            "dies": cl.pack(dies),
        })
    (OUT / "scan.json").write_text(
        json.dumps([{k: r[k] for k in ("tag", "T_py_C", "grad", "wl_mm")} for r in rows],
                   indent=1) + "\n", encoding="utf-8")

    pts = np.array([[r["T_py_C"], r["wl_mm"]] for r in rows])
    lo, hi = pts.min(0), pts.max(0)
    norm = (pts - lo) / np.maximum(hi - lo, 1e-9)
    keep = [0]  # 初始布局必留
    target = np.linspace(0.0, 1.0, N_KEEP)
    for f in target:
        # 沿温度从低到高的归一化轴均匀取点，线长只用来打破并列。
        order = np.argsort(np.abs(norm[:, 0] - f) + 0.15 * norm[:, 1])
        for idx in order:
            if int(idx) not in keep:
                keep.append(int(idx))
                break
    keep = keep[:N_KEEP]
    chosen = [rows[i] for i in keep]
    print("chosen:", flush=True)
    for r in chosen:
        print(f"  {r['tag']:8} T_py {r['T_py_C']:.2f}  wl {r['wl_mm']:.2f}", flush=True)
        (OUT / "layouts" / f"{r['tag']}.json").write_text(
            json.dumps(r["dies"], indent=1) + "\n", encoding="utf-8")
    (OUT / "chosen.json").write_text(
        json.dumps([r["tag"] for r in chosen]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
