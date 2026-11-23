#!/usr/bin/env python3
"""每一步 HotSpot 之后立刻按 CPU 块温度修正 McPAT 漏电，并和钉死在 380 K 的走法比。

不跑 mcpat 二进制。官方模板是 device_type=0（HP），温度只进入
cacti/technology.cc 的 I_off_n[T-300]。栅漏、动态、面积都不随温度变。
一次模板求解要重做 CACTI 划分，超过 1 小时，乘出来的就是这张表。
所以 10 K 一档把 380 K 的亚阈乘 Ioff 比值，栅漏不动。380 K 档与官方 txt 一致。

修正用的是这颗 CPU 自己的 HotSpot 块温度，不是封装最高温度。
GPU 和 HBM 的功耗不改。s=1，不换芯片。不改生产 SA，不算应力。
"""

from __future__ import annotations

import json
import math
import random
import shutil
import sys
import time
from copy import deepcopy
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

HS_BIN = ROOT / "HotSpot-master" / "hotspot"
OUT = ROOT / "热仿真流水线" / "runs" / "leak_step_gap"
CASE = ROOT / "热仿真流水线" / "runs" / "swap_library" / "practice_case.json"

SIDE = 0.072
ROWS = 64
N_SA = 30
SEED = 22
T0 = sa.T0
ALPHA = sa.ALPHA
MAX_INNER = 6
T_MIN_K, T_MAX_K = 300, 400

OFFICIAL_380 = {
    "Penryn": {"dyn_W": 47.9315, "sub_W": 10.7217, "gate_W": 0.588125},
    "T2": {"dyn_W": 50.421, "sub_W": 15.0682, "gate_W": 3.50762},
    "Xeon": {"dyn_W": 66.8485, "sub_W": 35.4132, "gate_W": 1.68429},
}
IOFF_HP = {
    "Penryn": {
        0: 2.8e-7, 10: 3.28e-7, 20: 3.81e-7, 30: 4.39e-7, 40: 5.02e-7,
        50: 5.69e-7, 60: 6.42e-7, 70: 7.2e-7, 80: 8.03e-7, 90: 8.91e-7, 100: 9.84e-7,
    },
    "T2": {
        0: 1.96e-7, 10: 2.29e-7, 20: 2.66e-7, 30: 3.05e-7, 40: 3.49e-7,
        50: 3.95e-7, 60: 4.45e-7, 70: 4.97e-7, 80: 5.48e-7, 90: 5.94e-7, 100: 6.3e-7,
    },
}
IOFF_HP["Xeon"] = IOFF_HP["T2"]


def round10_k(T_C: float) -> int:
    k = T_C + 273.15
    r = int(round(k / 10.0) * 10)
    return max(T_MIN_K, min(T_MAX_K, r))


def leak_at(model: str, T_K: int) -> dict:
    if (T_K - 300) not in IOFF_HP[model]:
        raise ValueError(f"{model} {T_K} K off the 10 K grid")
    off = OFFICIAL_380[model]
    ratio = IOFF_HP[model][T_K - 300] / IOFF_HP[model][80]
    sub = off["sub_W"] * ratio
    leak = sub + off["gate_W"]
    return {
        "McPAT_T_K": T_K,
        "Ioff_ratio": ratio,
        "dyn_W": off["dyn_W"],
        "sub_W": sub,
        "gate_W": off["gate_W"],
        "leak_W": leak,
        "power_W": off["dyn_W"] + leak,
    }


def cpu_of(chips: list[sa.Chip]) -> sa.Chip:
    return next(c for c in chips if c.kind == "CPU")


def apply_bin(chips: list[sa.Chip], T_K: int) -> dict:
    info = leak_at(cpu_of(chips).model, T_K)
    cpu_of(chips).power_W = info["power_W"]
    return info


def hotspot(chips: list[sa.Chip], work: Path) -> dict:
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    write_hotspot_inputs([c.rect() for c in chips], work, rows=ROWS, cols=ROWS)
    thermal = run_hotspot(HS_BIN, work)
    cpu_name = cpu_of(chips).name
    temps = thermal["temperatures_K"]
    if cpu_name not in temps:
        raise RuntimeError(f"no block {cpu_name} in {sorted(temps)}")
    return {
        "Tmax_C": float(thermal["Tmax_C"]),
        "T_cpu_C": float(temps[cpu_name] - 273.15),
        "block_T_C": {k: float(v - 273.15) for k, v in temps.items()},
    }


def evaluate(chips: list[sa.Chip], mode: str, work: Path) -> dict:
    """mode=fixed 钉在 380 K。mode=dynamic 每步看到块温度就改漏电，直到档位不动。"""
    cur = deepcopy(chips)
    if mode == "fixed":
        info = apply_bin(cur, 380)
        th = hotspot(cur, work)
        return {
            "mode": mode,
            "converged": True,
            "inner_iters": 1,
            "McPAT_T_K": 380,
            "power_W": info["power_W"],
            "leak_W": info["leak_W"],
            "dyn_W": info["dyn_W"],
            "Ioff_ratio": info["Ioff_ratio"],
            **th,
            "trace": [{"McPAT_T_K": 380, **info, **th}],
        }
    if mode != "dynamic":
        raise ValueError(mode)
    T = 380
    trace = []
    th = None
    info = None
    for i in range(MAX_INNER):
        info = apply_bin(cur, T)
        th = hotspot(cur, work / f"i{i}")
        nxt = round10_k(th["T_cpu_C"])
        trace.append({"iter": i, "next_K": nxt, **info, **th})
        if nxt == T:
            break
        T = nxt
    converged = trace[-1]["next_K"] == trace[-1]["McPAT_T_K"]
    return {
        "mode": mode,
        "converged": converged,
        "inner_iters": len(trace),
        "McPAT_T_K": info["McPAT_T_K"],
        "power_W": info["power_W"],
        "leak_W": info["leak_W"],
        "dyn_W": info["dyn_W"],
        "Ioff_ratio": info["Ioff_ratio"],
        "Tmax_C": th["Tmax_C"],
        "T_cpu_C": th["T_cpu_C"],
        "block_T_C": th["block_T_C"],
        "trace": trace,
    }


def j_of(tmax_c: float) -> float:
    return (tmax_c - 45.0) / 25.0


def chips_from_layout(layout: dict) -> list[sa.Chip]:
    out = []
    for c in layout["chips"]:
        out.append(sa.Chip(
            c["name"], c["kind"], c["x"], c["y"], c["w"], c["h"],
            c["power_W"], c["area_mm2"], c.get("model") or "",
        ))
    return out


def pack(ev: dict) -> dict:
    keep = ("mode", "converged", "inner_iters", "McPAT_T_K", "power_W", "leak_W",
            "dyn_W", "Ioff_ratio", "Tmax_C", "T_cpu_C", "block_T_C")
    return {k: ev[k] for k in keep}


def static_pairs(case: dict) -> list[dict]:
    rows = []
    for key, layout in case["layouts"].items():
        chips = chips_from_layout(layout)
        fixed = evaluate(chips, "fixed", OUT / "static" / key / "fixed")
        dyn = evaluate(chips, "dynamic", OUT / "static" / key / "dynamic")
        row = {
            "layout": key,
            "cpu": layout["cpu_model"],
            "gpu": layout["gpu_id"],
            "fixed": pack(fixed),
            "dynamic": pack(dyn),
            "dTmax_C": fixed["Tmax_C"] - dyn["Tmax_C"],
            "dT_cpu_C": fixed["T_cpu_C"] - dyn["T_cpu_C"],
            "dP_cpu_W": fixed["power_W"] - dyn["power_W"],
            "dLeak_W": fixed["leak_W"] - dyn["leak_W"],
        }
        rows.append(row)
        print(
            f"static {key:22} fixed {fixed['Tmax_C']:.3f} °C  "
            f"dyn {dyn['Tmax_C']:.3f} °C  dT {row['dTmax_C']:+.3f}  "
            f"dP {row['dP_cpu_W']:+.2f} W  bin {dyn['McPAT_T_K']}  iters {dyn['inner_iters']}",
            flush=True,
        )
    return rows


def run_sa(chips0: list[sa.Chip], mode: str, work: Path) -> dict:
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    rng = random.Random(SEED)
    init = evaluate(chips0, mode, work / "init")
    current = deepcopy(chips0)
    apply_bin(current, init["McPAT_T_K"])
    best = deepcopy(current)
    j_cur = j_best = j_of(init["Tmax_C"])
    best_ev = pack(init)
    T = T0
    steps = []
    accepted = 0
    for step in range(1, N_SA + 1):
        trial, op = sa.propose(current, SIDE, SIDE, rng)
        legal = sa.legal(trial, SIDE, SIDE)
        rec = {"step": step, "op": op, "legal": legal, "anneal_T": T}
        if not legal:
            steps.append(rec)
            T *= ALPHA
            continue
        ev = evaluate(trial, mode, work / "try")
        j_new = j_of(ev["Tmax_C"])
        delta = j_cur - j_new
        ap = 1.0 if delta >= 0 else math.exp(delta / max(T, 1e-9))
        take = rng.random() < ap
        rec.update({
            "accepted": take,
            "Tmax_C": ev["Tmax_C"],
            "T_cpu_C": ev["T_cpu_C"],
            "power_W": ev["power_W"],
            "leak_W": ev["leak_W"],
            "McPAT_T_K": ev["McPAT_T_K"],
            "inner_iters": ev["inner_iters"],
        })
        steps.append(rec)
        if take:
            current = trial
            apply_bin(current, ev["McPAT_T_K"])
            j_cur = j_new
            accepted += 1
            if j_new < j_best - 1e-12:
                best = deepcopy(current)
                j_best = j_new
                best_ev = pack(ev)
        T *= ALPHA
    return {
        "mode": mode,
        "accepted": accepted,
        "n_steps": N_SA,
        "init": pack(init),
        "best_as_seen": best_ev,
        "trace": steps,
        "best_chips": [
            {"name": c.name, "kind": c.kind, "model": c.model, "x": c.x, "y": c.y, "w": c.w, "h": c.h,
             "power_W": c.power_W, "area_mm2": c.area_mm2}
            for c in best
        ],
    }


def cross_score(arm: dict, work: Path) -> dict:
    chips = [
        sa.Chip(c["name"], c["kind"], c["x"], c["y"], c["w"], c["h"], c["power_W"], c["area_mm2"], c["model"])
        for c in arm["best_chips"]
    ]
    fixed = evaluate(chips, "fixed", work / "fixed")
    dyn = evaluate(chips, "dynamic", work / "dynamic")
    return {"fixed": pack(fixed), "dynamic": pack(dyn),
            "dTmax_C": fixed["Tmax_C"] - dyn["Tmax_C"],
            "dP_cpu_W": fixed["power_W"] - dyn["power_W"]}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    case = json.loads(CASE.read_text(encoding="utf-8"))
    t0 = time.time()
    static = static_pairs(case)
    (OUT / "static.json").write_text(json.dumps(static, indent=2) + "\n", encoding="utf-8")

    sa_rows = []
    for key in ("Penryn__GPU_3070", "Xeon__GPU_3070"):
        chips0 = chips_from_layout(case["layouts"][key])
        for mode in ("fixed", "dynamic"):
            print(f"SA {key} {mode}", flush=True)
            arm = run_sa(chips0, mode, OUT / "sa" / key / mode)
            arm["layout"] = key
            scored = cross_score(arm, OUT / "sa" / key / f"{mode}_best_rescore")
            arm["best_rescored"] = scored
            # chips are raw data; keep them, drop nothing
            sa_rows.append(arm)
            (OUT / "sa.json").write_text(json.dumps(sa_rows, indent=2) + "\n", encoding="utf-8")
            print(
                f"  seen best {arm['best_as_seen']['Tmax_C']:.3f} °C  "
                f"rescore fixed {scored['fixed']['Tmax_C']:.3f}  dyn {scored['dynamic']['Tmax_C']:.3f}  "
                f"accepted {arm['accepted']}/{N_SA}",
                flush=True,
            )

    note = {
        "question": "每一步 HotSpot 后立刻修正漏电，和全程钉在 380 K，温度差多少",
        "method": "HP Ioff 表缩放 380 K 亚阈。不跑 mcpat 二进制。CPU 用自己的块温度圆整到 10 K，不是封装 Tmax。",
        "not": "不改生产 SA。GPU 和 HBM 功耗不动。s=1。应力不在这次。Python/HotSpot 的温度不要和 ATSim 混。",
        "sa": {"n_steps": N_SA, "seed": SEED, "T0": T0, "alpha": ALPHA, "cost": "(Tmax_C-45)/25"},
        "seconds": time.time() - t0,
    }
    (OUT / "result.json").write_text(
        json.dumps({"note": note, "static": static, "sa": sa_rows}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"done {note['seconds']:.1f} s", flush=True)


if __name__ == "__main__":
    main()
