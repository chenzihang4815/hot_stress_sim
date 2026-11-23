#!/usr/bin/env python3
"""多芯片方案的退火。不改生产 SA，也不用四芯片那一套。

芯片：2 颗 Xeon、3 张 GPU、1 颗 HBM3，加上宏数不少于 64 的 CIM。
单宏的 jia / colonnade / basic 不进这个封装。
CIM 的面积和功耗来自 AccelForge cim_macros_scale 的 summary_table.json。

代价只看温度和温度梯度，两项在起点都是 1：
  (T_python-45)/(T_init-45) + grad/grad_init
材料图是整层硅，矩阵只分解一次。布局变了只改功率。
Python 温度不是 HotSpot。HotSpot 只在起点和代价最低的那一步各算一次。
CPU 功耗是 380 K 模板合计。GPU 功耗是同一份 backprop 的周期加权。
不算应力。
"""

from __future__ import annotations

import json
import random
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import run_fill_superposition as hs  # noqa: E402
import run_hw_area_field as hw  # noqa: E402
import run_sa_hw_smooth as sm  # noqa: E402

OUT = ROOT / "热仿真流水线" / "runs" / "sa_many_macro"
N = hs.N
N_STEPS = 60
ALPHA = 0.90
SEED = 7
T_AMB_C = 45.0

# 2 CPU + 3 GPU。背景每通道只算一次。
HBM_W = (0.088 + 2 * (1.035 - 0.088) + 3 * (0.844 - 0.088)) * 16
XEON_W = 66.8485 + 35.4132 + 1.68429


def die(name: str, role: str, area: float, power: float, j0: int, ib0: int, n_macros: int | None) -> dict:
    d = hw.prepare({"name": name, "area_mm2": area, "power_W": power})
    d.update({"role": role, "j0": j0, "ib0": ib0, "n_macros": n_macros})
    return d


def package() -> list[dict]:
    """64 格、间隙至少 1 格。高宏 CIM 都占得下格子，不再把单宏贴在别人身上。"""
    return [
        die("Xeon_a", "CPU", 421.676, XEON_W, 0, 0, None),
        die("GPU_2060S", "GPU", 445.0, 71.12082725952361, 17, 0, None),
        die("GPU_TITANX", "GPU", 471.0, 67.4071250008286, 34, 0, None),
        die("CIM_jia_m64", "CIM", 273.99848050829274, 0.02051619509144024, 52, 0, 64),
        die("Xeon_b", "CPU", 421.676, XEON_W, 0, 23, None),
        die("GPU_3070", "GPU", 392.0, 66.88781634908364, 17, 23, None),
        die("HBM3", "HBM", 121.0, HBM_W, 32, 23, None),
        die("CIM_colonnade_m512", "CIM", 118.456264372224, 0.2967109271348416, 41, 23, 512),
        die("CIM_colonnade_m256", "CIM", 59.228132186112, 0.29677526103697566, 0, 46, 256),
        die("CIM_basic_m256", "CIM", 11.9601949587906, 0.0022112712579585407, 7, 46, 256),
        die("CIM_colonnade_m64", "CIM", 14.807033046528, 0.29716126444977986, 11, 46, 64),
    ]


def nudge(dies: list[dict], rng: random.Random) -> tuple[list[dict], str]:
    trial = deepcopy(dies)
    i = rng.randrange(len(trial))
    step = rng.choice([1, 1, 2])
    dj, di = rng.choice([(1, 0), (-1, 0), (0, 1), (0, -1)])
    trial[i]["j0"] += dj * step
    trial[i]["ib0"] += di * step
    return trial, f"nudge {trial[i]['name']} ({dj * step},{di * step})"


def rotate(dies: list[dict], rng: random.Random) -> tuple[list[dict], str]:
    trial = deepcopy(dies)
    i = rng.randrange(len(trial))
    trial[i]["w"], trial[i]["h"] = trial[i]["h"], trial[i]["w"]
    return trial, f"rotate {trial[i]['name']}"


def score(ev: dict, scales: dict) -> float:
    return (ev["Tmax_C"] - T_AMB_C) / scales["dT"] + ev["grad_max"] / scales["g0"]


def pack(dies: list[dict]) -> list[dict]:
    keys = ("name", "role", "j0", "ib0", "w", "h", "power_W", "area_mm2", "grid_area_mm2", "n_macros")
    return [{k: d.get(k) for k in keys} for d in dies]


def slim(ev: dict) -> dict:
    return {k: ev[k] for k in ("grad_max", "grad_rms", "Tmax_C")}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "layouts").mkdir(exist_ok=True)
    (OUT / "fields").mkdir(exist_ok=True)
    start = package()
    if not sm.legal(start):
        raise SystemExit("起点摆不下，间隙或边界不合法")
    k = np.full((N, N), hs.K_SI)
    mat, g_amb = hs.build_operator(k)
    fac = splu(mat.tocsc())

    rng_m = random.Random(SEED)
    rng_a = random.Random(SEED + 7919)
    current = deepcopy(start)
    cur = sm.evaluate(current, fac, g_amb)
    scales = {"dT": max(cur["Tmax_C"] - T_AMB_C, 1e-6), "g0": max(cur["grad_max"], 1e-9)}
    j_cur = score(cur, scales)
    best = deepcopy(current)
    best_ev = cur
    best_j = j_cur
    temp = 0.20
    log = [{
        "step": 0, "op": "init", "accept": True, "J": j_cur,
        "grad_max": cur["grad_max"], "grad_rms": cur["grad_rms"], "Tmax_C": cur["Tmax_C"],
    }]
    print(
        f"init J {j_cur:.3f}  grad {cur['grad_max']:.3f}  T {cur['Tmax_C']:.2f}  "
        f"chips {len(start)}",
        flush=True,
    )
    for step in range(1, N_STEPS + 1):
        if rng_m.random() < 0.75:
            trial, op = nudge(current, rng_m)
            move = "nudge"
        else:
            trial, op = rotate(current, rng_m)
            move = "rotate"
        if not sm.legal(trial):
            log.append({"step": step, "op": op, "move": move, "accept": False, "reason": "illegal"})
            temp *= ALPHA
            continue
        ev = sm.evaluate(trial, fac, g_amb)
        j_new = score(ev, scales)
        take, ap = sm.accept(j_cur, j_new, temp, rng_a)
        rec = {
            "step": step, "op": op, "move": move, "accept": take, "ap": ap,
            "J": j_new, "grad_max": ev["grad_max"], "grad_rms": ev["grad_rms"],
            "Tmax_C": ev["Tmax_C"], "anneal_T": temp,
        }
        if take:
            current = trial
            j_cur = j_new
            if j_new < best_j - 1e-12:
                best = deepcopy(trial)
                best_ev = ev
                best_j = j_new
                rec["new_best"] = True
                print(
                    f"  {step:02d} BEST {op}  J {j_new:.3f}  grad {ev['grad_max']:.3f}  "
                    f"T {ev['Tmax_C']:.2f}",
                    flush=True,
                )
        log.append(rec)
        temp *= ALPHA

    np.save(OUT / "fields" / "init.npy", cur["T"])
    np.save(OUT / "fields" / "best.npy", best_ev["T"])
    (OUT / "layouts" / "init.json").write_text(json.dumps(pack(start), indent=2) + "\n", encoding="utf-8")
    (OUT / "layouts" / "best.json").write_text(json.dumps(pack(best), indent=2) + "\n", encoding="utf-8")
    (OUT / "trace.json").write_text(json.dumps(log, indent=2) + "\n", encoding="utf-8")

    hs_init = sm.hotspot(start, OUT / "hotspot_init")
    hs_best = sm.hotspot(best, OUT / "hotspot_best")
    report = {
        "cost": "(T_python-45)/(T_init-45) + grad/grad_init",
        "material": "整层硅，k=130。不是环氧填充。矩阵只分解一次。",
        "steps": N_STEPS,
        "seed": SEED,
        "alpha": ALPHA,
        "n_chips": len(start),
        "hbm_W": HBM_W,
        "note": (
            "不用四芯片方案。CIM 只用 cim_macros_scale 里 n_macros>=64 的结果，"
            "不是 more_combos 里的单宏。CPU 功耗是 Xeon 380 K 合计，不是 matmul 重算。"
            "GPU 功耗是 backprop 周期加权。Python 温度不能当成 HotSpot。"
            "最优是代价最低的那一步，不是最后一步。没有算应力，没有改生产 SA。"
        ),
        "init": {**slim(cur), "J": score(cur, scales), "hotspot_C": hs_init},
        "best": {**slim(best_ev), "J": best_j, "hotspot_C": hs_best},
        "accepted": sum(1 for r in log if r.get("accept") and r["step"] > 0),
        "chips": pack(start),
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(
        f"best J {best_j:.3f}  grad {best_ev['grad_max']:.3f}  "
        f"python {best_ev['Tmax_C']:.2f}  HotSpot init {hs_init:.2f} best {hs_best:.2f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
