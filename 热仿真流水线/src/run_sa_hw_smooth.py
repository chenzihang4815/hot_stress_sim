#!/usr/bin/env python3
"""退火的每一步除了挪位置，还看整层温度，决定要不要换仓库里的硬件。

仓库只有现成的硬件，不编功耗：
- CPU 任务：Penryn / T2 / Xeon，功耗是各自 McPAT 380 K 总功耗
- GPU 任务：3070 / 2060S / TITAN X，同一份 backprop
CIM 和 HBM 不换。

起点是最小的一套（Penryn + 3070），热点最尖。每一步做两件提案：
1. 把一颗芯片挪一格或两格。
2. 在当前温度场上，找足迹里梯度最尖的那颗 CPU 或 GPU，
   只试面积更大、而且放得下的仓库硬件（左下角不动，往外长）。
两件提案里，整层最大梯度更小的那件交给退火。代价就是这个最大梯度，
不是峰值温度。材料图用整层硅：前面已经量过，填充换成别的材料，
梯度几乎不动，而整层硅可以只分解一次矩阵。

对照是同一条随机序列、只挪位置、不许换硬件。
绝对温度不是 HotSpot。起点和两个终点另用 HotSpot 记一笔峰值。
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
from geometry import Rect  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

OUT = ROOT / "热仿真流水线" / "runs" / "sa_hw_smooth"
HS_BIN = ROOT / "HotSpot-master" / "hotspot"
N = hs.N
CW_MM = hs.CW * 1e3
GAP = 1  # 格，1.125 mm，大于流水线的 0.5 mm 最小间隙
N_STEPS = 30
ALPHA = 0.90
SEED = 7


def clonable(src: dict) -> dict:
    return {
        "name": src["name"], "role": src["role"],
        "w": src["w"], "h": src["h"],
        "power_W": src["power_W"], "area_mm2": src["area_mm2"],
        "grid_area_mm2": src["grid_area_mm2"],
        "j0": src.get("j0", 0), "ib0": src.get("ib0", 0),
    }


def warehouse() -> dict[str, list[dict]]:
    cpus, gpus, hbm, cim = hw.place_all()
    for c in cpus:
        c["role"] = "CPU"
    for g in gpus:
        g["role"] = "GPU"
    hbm["role"] = "HBM"
    cim["role"] = "CIM"
    return {
        "CPU": [clonable(c) for c in cpus],
        "GPU": [clonable(g) for g in gpus],
        "fixed": [clonable(hbm), clonable(cim)],
        "start": [clonable(cpus[0]), clonable(gpus[0]), clonable(hbm), clonable(cim)],
    }


def gap_cells(a: dict, b: dict) -> int:
    dx = max(b["j0"] - (a["j0"] + a["w"]), a["j0"] - (b["j0"] + b["w"]))
    dy = max(b["ib0"] - (a["ib0"] + a["h"]), a["ib0"] - (b["ib0"] + b["h"]))
    if dx < 0 and dy < 0:
        return max(dx, dy)
    if dx < 0:
        return dy
    if dy < 0:
        return dx
    return min(dx, dy)


def legal(dies: list[dict]) -> bool:
    for d in dies:
        if d["j0"] < 0 or d["ib0"] < 0:
            return False
        if d["j0"] + d["w"] > N or d["ib0"] + d["h"] > N:
            return False
    for i in range(len(dies)):
        for j in range(i + 1, len(dies)):
            if gap_cells(dies[i], dies[j]) < GAP:
                return False
    return True


def grad_max(T: np.ndarray) -> float:
    gx = np.abs(np.diff(T, axis=1)).max() / CW_MM
    gy = np.abs(np.diff(T, axis=0)).max() / CW_MM
    return float(max(gx, gy))


def grad_rms(T: np.ndarray) -> float:
    gx = np.diff(T, axis=1) / CW_MM
    gy = np.diff(T, axis=0) / CW_MM
    return float(np.sqrt(0.5 * (np.mean(gx ** 2) + np.mean(gy ** 2))))


def sharpness(T: np.ndarray, die: dict) -> float:
    """这颗芯片足迹边上的最大梯度，K/mm。"""
    cc = set(hw.cells(die))
    m = 0.0
    for i, j in cc:
        if j + 1 < N:
            m = max(m, abs(T[i, j + 1] - T[i, j]) / CW_MM)
        if i + 1 < N:
            m = max(m, abs(T[i + 1, j] - T[i, j]) / CW_MM)
    return float(m)


def evaluate(dies: list[dict], fac, g_amb) -> dict:
    k = np.full((N, N), hs.K_SI)
    state = hw.solve(k, hw.power_map(dies), fac, g_amb)
    T = state[0]
    per = []
    for d in dies:
        if d["role"] not in ("CPU", "GPU"):
            continue
        per.append({
            "name": d["name"], "role": d["role"],
            "sharp_K_per_mm": sharpness(T, d),
            "area_mm2": d["area_mm2"], "power_W": d["power_W"],
        })
    return {
        "T": T,
        "grad_max": grad_max(T),
        "grad_rms": grad_rms(T),
        "Tmax_C": float(T.max() - 273.15),
        "dies": per,
    }


def layout_move(dies: list[dict], rng: random.Random) -> tuple[list[dict], str]:
    trial = deepcopy(dies)
    movable = [i for i, d in enumerate(trial) if d["role"] in ("CPU", "GPU")]
    i = rng.choice(movable)
    step = rng.choice([1, 1, 2])
    dj, di = rng.choice([(1, 0), (-1, 0), (0, 1), (0, -1)])
    trial[i]["j0"] += dj * step
    trial[i]["ib0"] += di * step
    op = f"nudge {trial[i]['name']} ({dj * step},{di * step}) cells"
    return trial, op


def hardware_move(dies: list[dict], info: dict, catalog: dict[str, list[dict]]) -> tuple[list[dict] | None, str]:
    """最尖的那颗，换成放得下、而且整层梯度最低的更大硬件。"""
    ranked = sorted(info["dies"], key=lambda r: -r["sharp_K_per_mm"])
    for victim in ranked:
        idx = next(i for i, d in enumerate(dies) if d["name"] == victim["name"])
        cur = dies[idx]
        cands = [p for p in catalog[cur["role"]] if p["area_mm2"] > cur["area_mm2"] + 1e-6]
        best = None
        best_j = None
        for part in cands:
            trial = deepcopy(dies)
            trial[idx] = clonable(part)
            trial[idx]["j0"] = cur["j0"]
            trial[idx]["ib0"] = cur["ib0"]
            if not legal(trial):
                continue
            ev = evaluate(trial, info["fac"], info["g"])
            if best is None or ev["grad_max"] < best_j:
                best = (trial, part["name"], ev["grad_max"])
                best_j = ev["grad_max"]
        if best is not None:
            trial, name, _ = best
            return trial, f"swap {cur['name']} -> {name} (sharp {victim['sharp_K_per_mm']:.3f})"
    return None, "no larger part fits"


def accept(j_cur: float, j_new: float, temp: float, rng: random.Random) -> tuple[bool, float]:
    delta = j_cur - j_new
    if delta >= 0:
        return True, 1.0
    ap = float(np.exp(delta / max(temp, 1e-9)))
    return rng.random() < ap, ap


def run(mode: str, start: list[dict], catalog, fac, g_amb) -> dict:
    rng = random.Random(SEED)
    current = deepcopy(start)
    info = {"fac": fac, "g": g_amb}
    cur = evaluate(current, fac, g_amb)
    best = deepcopy(current)
    best_ev = cur
    temp = max(0.15, 0.25 * cur["grad_max"])
    log = [{
        "step": 0, "op": "init", "accept": True,
        "grad_max": cur["grad_max"], "grad_rms": cur["grad_rms"],
        "Tmax_C": cur["Tmax_C"], "hardware": [d["name"] for d in current],
    }]
    print(
        f"[{mode}] init grad {cur['grad_max']:.3f} K/mm  "
        f"Tmax {cur['Tmax_C']:.2f} °C  {[d['name'] for d in current]}",
        flush=True,
    )
    for step in range(1, N_STEPS + 1):
        layout, op_l = layout_move(current, rng)
        proposals = []
        if legal(layout):
            ev = evaluate(layout, fac, g_amb)
            proposals.append((layout, op_l, ev))
        if mode == "hardware":
            info_now = {"fac": fac, "g": g_amb, "dies": cur["dies"]}
            hw_trial, op_h = hardware_move(current, info_now, catalog)
            if hw_trial is not None:
                ev = evaluate(hw_trial, fac, g_amb)
                proposals.append((hw_trial, op_h, ev))
        if not proposals:
            log.append({"step": step, "op": "illegal", "accept": False})
            temp *= ALPHA
            continue
        trial, op, ev = min(proposals, key=lambda t: t[2]["grad_max"])
        take, ap = accept(cur["grad_max"], ev["grad_max"], temp, rng)
        rec = {
            "step": step, "op": op, "accept": take, "ap": ap,
            "grad_max": ev["grad_max"], "grad_rms": ev["grad_rms"],
            "Tmax_C": ev["Tmax_C"], "hardware": [d["name"] for d in trial],
            "T": temp,
        }
        if take:
            current = trial
            cur = ev
            if ev["grad_max"] < best_ev["grad_max"]:
                best = deepcopy(trial)
                best_ev = ev
                rec["new_best"] = True
            print(
                f"  {step:02d} {'BEST' if rec.get('new_best') else 'take'} {op}  "
                f"grad {ev['grad_max']:.3f}  Tmax {ev['Tmax_C']:.2f}",
                flush=True,
            )
        else:
            print(f"  {step:02d} skip {op}  grad {ev['grad_max']:.3f} ap {ap:.2f}", flush=True)
        log.append(rec)
        temp *= ALPHA
    return {"mode": mode, "best": best, "best_ev": {k: v for k, v in best_ev.items() if k != "T"},
            "log": log}


def to_rects(dies: list[dict]) -> list[Rect]:
    rects = []
    for d in dies:
        rects.append(Rect(
            d["name"], d["j0"] * hs.CW, d["ib0"] * hs.CW,
            d["w"] * hs.CW, d["h"] * hs.CW,
            d["power_W"], d["grid_area_mm2"], d["role"],
        ))
    return rects


def hotspot(dies: list[dict], work: Path) -> float:
    write_hotspot_inputs(to_rects(dies), work, rows=N, cols=N)
    out = run_hotspot(HS_BIN, work)
    temps = out["temperatures_K"]
    chips = [v for k, v in temps.items() if not k.startswith("anchor")]
    return max(chips) - 273.15


def dump_die(d: dict) -> dict:
    return {k: d[k] for k in ("name", "role", "j0", "ib0", "w", "h", "power_W", "area_mm2", "grid_area_mm2")}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cat = warehouse()
    start = cat["start"]
    if not legal(start):
        raise SystemExit("起点布局不合法")
    k = np.full((N, N), hs.K_SI)
    mat, g_amb = hs.build_operator(k)
    fac = splu(mat.tocsc())
    (OUT / "warehouse.json").write_text(json.dumps({
        "CPU": [dump_die(c) for c in cat["CPU"]],
        "GPU": [dump_die(c) for c in cat["GPU"]],
        "fixed": [dump_die(c) for c in cat["fixed"]],
        "note": "CPU 功耗是 McPAT 380 K 总功耗；GPU 功耗是同一份 backprop。只允许换成面积更大的。",
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    results = {}
    for mode in ("layout", "hardware"):
        results[mode] = run(mode, start, cat, fac, g_amb)
        (OUT / f"{mode}_trace.json").write_text(
            json.dumps(results[mode]["log"], indent=2) + "\n", encoding="utf-8")

    hs_init = hotspot(start, OUT / "hotspot_init")
    summary = {"hotspot_init_C": hs_init, "modes": {}}
    for mode, res in results.items():
        hs_best = hotspot(res["best"], OUT / f"hotspot_{mode}")
        summary["modes"][mode] = {
            "grad_max_init": res["log"][0]["grad_max"],
            "grad_max_best": res["best_ev"]["grad_max"],
            "grad_rms_best": res["best_ev"]["grad_rms"],
            "python_Tmax_init_C": res["log"][0]["Tmax_C"],
            "python_Tmax_best_C": res["best_ev"]["Tmax_C"],
            "hotspot_Tmax_best_C": hs_best,
            "hardware": [d["name"] for d in res["best"]],
            "dies": [dump_die(d) for d in res["best"]],
        }
        print(
            f"{mode}: grad {res['log'][0]['grad_max']:.3f} -> {res['best_ev']['grad_max']:.3f} K/mm  "
            f"HotSpot {hs_init:.2f} -> {hs_best:.2f} °C  hardware {summary['modes'][mode]['hardware']}",
            flush=True,
        )
    summary["note"] = (
        "代价是芯片层最大温度梯度，不是峰值。Python 场是整层硅、没有 spreader 外伸。"
        "HotSpot 峰值只在起点和两个终点各算一次。"
    )
    (OUT / "report.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"写出 {OUT / 'report.json'}", flush=True)


if __name__ == "__main__":
    main()
