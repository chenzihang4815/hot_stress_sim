#!/usr/bin/env python3
"""多组闭环。

内环代价是温度、梯度、线长各三分之一。定价用场叠加：填充和硅一样时
因子留着复用；填充不同时只在第一次重解，之后用残差迭代修正。
应力不进内环。布局定下来之后 ATSim 判一次，再调频率，应力仍超再修复。

外环布局冻住。向导按这一步的 ATSim 块平均来缩，不按 HotSpot，
因为 HotSpot 的芯片层不吃填充。倍率步长 0.1，在向导低于 90 C 的
前提下往上加，直到加不动。认门槛的仍是 ATSim：芯粒层最热一格 < 90 C，
正常芯片内部应力 < 200 MPa。加完若越出门槛，按动态功耗最大的那颗
每次降 0.1，最多 12 次。

温度对比三列都是块平均。应力是内部最大值，不是整包最大。
GPU 吞吐是 backprop 指令数 / 主频。CPU 没有这次任务的操作数。
"""

from __future__ import annotations

import gc
import json
import os
import subprocess
import sys
import time
import traceback
from copy import deepcopy
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import math
import random

from blocks import Block
from bstree import random_bstree_layout
from geometry import Rect
from lbfgs import random_lbfgs_layout
from mcts import random_mcts_layout

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from scipy.sparse.linalg import splu

import run_closed_loop_multidie as cl
import run_fill_superposition as hs
import run_hw_area_field as hw
import run_sa_hw_smooth as sm
import run_wirelength as wire

FONT = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
FENICS = Path("/home/czh/miniconda3/envs/fenicsx/bin/python")
SVG_PL = ROOT / "HotSpot-master" / "scripts" / "grid_thermal_map.pl"
SUITE = Path(os.environ.get(
    "SUITE_DIR",
    str(ROOT / "热仿真流水线" / "runs" / "closed_loop_suite"),
))
N = hs.N
SEED = 11
ALPHA = 0.90
N_INNER = int(os.environ.get("SUITE_STEPS", "40"))
SKIP_ATSIM = os.environ.get("SUITE_NO_ATSIM") == "1"
ONLY = {s for s in os.environ.get("SUITE_ONLY", "").split(",") if s}

# 导热是工程里的 0.8。弹性模量、泊松比、热膨胀材料表没有环氧卡，
# 这一次只写进案例 json，借用已有 organic 的力学数，不改 DEFAULT_MATERIALS。
EPOXY = {
    "kappa0_W_mK": 0.8,
    "T_ref_K": 300.0,
    "alpha": 0.0,
    "rho_kg_m3": 1800.0,
    "cp_J_kgK": 1000.0,
    "E_Pa": 2.5e10,
    "nu": 0.30,
    "alpha_cte_1_K": 1.5e-5,
}
SIM_INSTS = {
    "GPU_2060S": 370711,
    "GPU_3070": 370653,
    "GPU_TITANX": 370707,
}
N_COMPUTES = {
    "CIM_jia_m64": 589824,
    "CIM_colonnade_m512": 16384,
    "CIM_colonnade_m256": 16384,
    "CIM_basic_m256": 1024,
    "CIM_colonnade_m64": 16384,
}


def style() -> None:
    font_manager.fontManager.addfont(FONT)
    plt.rcParams["font.family"] = "Noto Sans CJK JP"
    plt.rcParams["axes.unicode_minus"] = False


def mem_mb() -> float:
    for line in open("/proc/meminfo", encoding="utf-8"):
        if line.startswith("MemAvailable"):
            return int(line.split()[1]) / 1024.0
    return 0.0


def wait_ram(need: float = 3000.0) -> None:
    for _ in range(12):
        free = mem_mb()
        print(f"  MemAvailable {free:.0f} MB", flush=True)
        if free >= need:
            return
        time.sleep(15)
    raise RuntimeError(f"可用内存低于 {need:.0f} MB，停止以免卡死")


def put(cat: dict, role: str, part: str, slot: str, j0: int, ib0: int) -> dict:
    die = deepcopy(cat["by"][(role, part)])
    die["name"] = f"{part}_{slot}"
    die["slot"] = slot
    die["j0"] = j0
    die["ib0"] = ib0
    die["s"] = 1.0
    return die


def cim_at(index: int, j0: int, ib0: int) -> dict:
    name, area, power, nmac = cl.CIMS[index]
    die = cl.fixed_die(name, "CIM", area, power, nmac)
    die["j0"], die["ib0"] = j0, ib0
    return die


def hbm_at(j0: int, ib0: int) -> dict:
    die = cl.fixed_die("HBM3", "HBM", 121.0, 0.0, None)
    die["j0"], die["ib0"] = j0, ib0
    return die


FILLS = {
    "silicon": {"k": 130.0, "mechanics": None},
    "epoxy": {"k": 0.8, "mechanics": "borrowed_organic"},
    "underfill": {"k": 1.0, "mechanics": "borrowed_organic"},
    "air": {"k": 0.028, "mechanics": None},
}


def fill_of(name: str) -> dict:
    """内置四种，也接受 "名字:导热" 这样的自定义。"""
    if ":" in name:
        label, k = name.split(":", 1)
        return {"name": label, "k": float(k), "mechanics": None, "custom": True}
    if name not in FILLS:
        raise ValueError(f"填充 {name} 不在 {sorted(FILLS)}，自定义写成 名字:导热")
    return {"name": name, "custom": False, **FILLS[name]}


def members(cat: dict, kind: str) -> list[dict]:
    """只定芯片清单，不定位置。位置交给三种初始布局。"""
    if kind == "n5":
        dies = [put(cat, "CPU", "Xeon", "Ca", 0, 0), put(cat, "GPU", "GPU_3070", "Ga", 0, 0),
                cim_at(0, 0, 0), hbm_at(0, 0), cim_at(1, 0, 0)]
    elif kind == "n6":
        dies = [put(cat, "CPU", "Xeon", "Ca", 0, 0), put(cat, "GPU", "GPU_3070", "Ga", 0, 0),
                put(cat, "GPU", "GPU_2060S", "Gb", 0, 0), cim_at(0, 0, 0), hbm_at(0, 0), cim_at(1, 0, 0)]
    elif kind == "n7":
        dies = [put(cat, "CPU", "Xeon", "Ca", 0, 0), put(cat, "CPU", "T2", "Cb", 0, 0),
                put(cat, "GPU", "GPU_3070", "Ga", 0, 0), cim_at(0, 0, 0), hbm_at(0, 0),
                cim_at(1, 0, 0), cim_at(2, 0, 0)]
    elif kind == "n8":
        dies = [put(cat, "CPU", "Xeon", "Ca", 0, 0), put(cat, "GPU", "GPU_2060S", "Ga", 0, 0),
                put(cat, "GPU", "GPU_3070", "Gb", 0, 0), cim_at(0, 0, 0),
                put(cat, "CPU", "Xeon", "Cb", 0, 0), hbm_at(0, 0), cim_at(1, 0, 0), cim_at(2, 0, 0)]
    elif kind == "n9":
        dies = [put(cat, "CPU", "Xeon", "Ca", 0, 0), put(cat, "GPU", "GPU_2060S", "Ga", 0, 0),
                put(cat, "GPU", "GPU_3070", "Gb", 0, 0), cim_at(0, 0, 0),
                put(cat, "CPU", "Xeon", "Cb", 0, 0), hbm_at(0, 0), cim_at(1, 0, 0),
                cim_at(2, 0, 0), cim_at(3, 0, 0)]
    elif kind == "n11":
        dies = cl.place(cat)
    else:
        raise ValueError(kind)
    cl.refresh_hbm(dies)
    cl.apply_power(dies)
    return dies


def shelf(cat: dict, kind: str) -> list[dict]:
    if kind == "n11":
        return cl.place(cat)
    if kind == "n5":
        dies = [
            put(cat, "CPU", "Xeon", "Ca", 0, 0),
            put(cat, "GPU", "GPU_3070", "Ga", 18, 0),
            cim_at(0, 34, 0),
            hbm_at(48, 0),
            cim_at(1, 0, 23),
        ]
    elif kind == "n6":
        dies = [
            put(cat, "CPU", "Xeon", "Ca", 0, 0),
            put(cat, "GPU", "GPU_3070", "Ga", 18, 0),
            put(cat, "GPU", "GPU_2060S", "Gb", 34, 0),
            cim_at(0, 52, 0),
            hbm_at(0, 23),
            cim_at(1, 18, 24),
        ]
    elif kind == "n7":
        dies = [
            put(cat, "CPU", "Xeon", "Ca", 0, 0),
            put(cat, "CPU", "T2", "Cb", 18, 0),
            put(cat, "GPU", "GPU_3070", "Ga", 32, 0),
            cim_at(0, 48, 0),
            hbm_at(0, 23),
            cim_at(1, 18, 24),
            cim_at(2, 32, 24),
        ]
    elif kind == "n8":
        dies = [
            put(cat, "CPU", "Xeon", "Ca", 0, 0),
            put(cat, "GPU", "GPU_2060S", "Ga", 18, 0),
            put(cat, "GPU", "GPU_3070", "Gb", 36, 0),
            cim_at(0, 52, 0),
            put(cat, "CPU", "Xeon", "Cb", 0, 24),
            hbm_at(18, 24),
            cim_at(1, 27, 24),
            cim_at(2, 38, 24),
        ]
    elif kind == "n9":
        dies = [
            put(cat, "CPU", "Xeon", "Ca", 0, 0),
            put(cat, "GPU", "GPU_2060S", "Ga", 18, 0),
            put(cat, "GPU", "GPU_3070", "Gb", 36, 0),
            cim_at(0, 52, 0),
            put(cat, "CPU", "Xeon", "Cb", 0, 24),
            hbm_at(18, 24),
            cim_at(1, 27, 24),
            cim_at(2, 38, 24),
            cim_at(3, 46, 24),
        ]
    else:
        raise ValueError(kind)
    cl.refresh_hbm(dies)
    cl.apply_power(dies)
    return dies


def _free_spot(dies: list[dict], mover: dict) -> tuple[int, int] | None:
    """框里第一个不与其余芯粒重叠、也不越界的位置。"""
    others = [d for d in dies if d is not mover]
    for ib0 in range(hs.N - mover["h"] + 1):
        for j0 in range(hs.N - mover["w"] + 1):
            trial = dict(mover)
            trial["j0"], trial["ib0"] = j0, ib0
            if all(sm.gap_cells(trial, d) >= sm.GAP for d in others):
                return j0, ib0
    return None


def _separate(dies: list[dict]) -> bool:
    """重叠的芯粒整颗挪到第一个空位，挪不开就转 90 度再找。"""
    for _ in range(len(dies) * 4):
        if sm.legal(dies):
            return True
        worst = None
        for i in range(len(dies)):
            for j in range(i + 1, len(dies)):
                gap = sm.gap_cells(dies[i], dies[j])
                if worst is None or gap < worst[0]:
                    worst = (gap, j)
        mover = dies[worst[1]]
        spot = _free_spot(dies, mover)
        if spot is None and mover["w"] != mover["h"]:
            mover["w"], mover["h"] = mover["h"], mover["w"]
            spot = _free_spot(dies, mover)
        if spot is None:
            return False
        mover["j0"], mover["ib0"] = spot
    return sm.legal(dies)


def snap_layout(dies: list[dict], rects: list[Rect]) -> list[dict]:
    """三种算法给的是米。落到 1.125 mm 的格子上，再拉开到合法间隙。"""
    by_name = {r.name: r for r in rects}
    out = deepcopy(dies)
    for die in out:
        rect = by_name[die["name"]]
        die["j0"] = int(round(rect.x / (cl.CW_MM * 1e-3)))
        die["ib0"] = int(round(rect.y / (cl.CW_MM * 1e-3)))
        die["j0"] = min(max(0, die["j0"]), hs.N - die["w"])
        die["ib0"] = min(max(0, die["ib0"]), hs.N - die["h"])
    if not _separate(out):
        raise RuntimeError("初始布局落格后仍不合法")
    return out


def initial_layout(dies: list[dict], k_fill: float) -> tuple[list[dict], str, dict]:
    """bstree、mcts、lbfgs 各跑一次，取得分最低的。这里要应力，所以用 ATSim。"""
    blocks = [Block(d["name"], d["role"], "catalog", d["area_mm2"], d["power_W"], "") for d in dies]
    scored = []
    for name, fn in (
        ("bstree", random_bstree_layout),
        ("mcts", random_mcts_layout),
        ("lbfgs", random_lbfgs_layout),
    ):
        rects, meta = fn(blocks, seed=SEED)
        placed = snap_layout(dies, rects)
        tag = "epoxy" if abs(k_fill - 0.8) < 1e-9 else "silicon"
        at = run_atsim(placed, tag, Path("/home/czh/.claude/jobs/367ec9f7/tmp") / f"init_{name}.json",
                       Path("/home/czh/.claude/jobs/367ec9f7/tmp") / f"init_{name}_result.json")
        score = at["die_cell_max_C"] / cl.T_GATE + at["vm_normal_chip_max_MPa"] / cl.S_GATE
        scored.append((score, name, placed, at))
        print(f"  init {name}  T {at['die_cell_max_C']:.1f}  vm {at['vm_normal_chip_max_MPa']:.1f}", flush=True)
    score, name, placed, at = min(scored, key=lambda r: r[0])
    return placed, name, {"score": score, "T": at["die_cell_max_C"], "vm": at["vm_normal_chip_max_MPa"]}


def iterate(dies: list[dict], k_fill: float, base: np.ndarray) -> np.ndarray:
    """填充和硅不同时只重解这一次，返回相对硅填充的残差。"""
    k = hw.paint(dies, k_fill)
    mat, g_amb = hs.build_operator(k)
    fac = splu(mat.tocsc())
    truth = hw.solve(k, hw.power_map(dies), fac, g_amb)[0]
    return truth - base


def joint(dies: list[dict], k_fill: float, check: bool) -> tuple[np.ndarray, float | None]:
    cl.apply_power(dies)
    k = hw.paint(dies, k_fill)
    mat, g_amb = hs.build_operator(k)
    fac = splu(mat.tocsc())
    truth = hw.solve(k, hw.power_map(dies), fac, g_amb)[0]
    err = None
    if check:
        acc = np.zeros_like(truth)
        for die in dies:
            one = hw.solve(k, hw.power_map([die]), fac, g_amb)[0]
            acc += one - hs.T_AMB
        acc += hs.T_AMB
        err = float(np.max(np.abs(acc - truth)))
    del fac, mat
    return truth, err


def rises_of(dies: list[dict], k_fill: float) -> dict[str, np.ndarray]:
    cl.apply_power(dies)
    k = hw.paint(dies, k_fill)
    mat, g_amb = hs.build_operator(k)
    fac = splu(mat.tocsc())
    out = {}
    for die in dies:
        one = hw.solve(k, hw.power_map([die]), fac, g_amb)[0]
        out[die["name"]] = one - hs.T_AMB
    del fac, mat
    return out


def field_at(dies: list[dict], rises: dict[str, np.ndarray], scales: dict[str, float]) -> np.ndarray:
    acc = np.zeros((N, N))
    for die in dies:
        if die["role"] in ("CPU", "GPU"):
            p0 = die["leak_W"] + die["dyn_W"]
            scale = float(scales[die["slot"]])
            power = die["leak_W"] + scale * die["dyn_W"]
        else:
            p0 = die["power_W"]
            power = p0
        if p0 <= 0:
            continue
        acc += rises[die["name"]] * (power / p0)
    return acc + hs.T_AMB


def py_blocks(temp: np.ndarray, dies: list[dict]) -> list[dict]:
    rows = []
    for die in dies:
        i0, i1 = hs.rows_of({"ib0": die["ib0"], "ib1": die["ib0"] + die["h"]})
        sub = temp[i0:i1, die["j0"]:die["j0"] + die["w"]]
        rows.append({
            "name": die["name"],
            "mean_C": float(sub.mean() - 273.15),
            "max_C": float(sub.max() - 273.15),
        })
    return rows


def hottest(rows: list[dict], key: str = "mean_C") -> dict:
    row = max(rows, key=lambda r: r[key])
    return {"name": row["name"], "C": row[key]}


def hs_blocks(temps: dict, dies: list[dict]) -> list[dict]:
    rows = []
    for die in dies:
        if die["name"] not in temps:
            continue
        rows.append({"name": die["name"], "mean_C": float(temps[die["name"]])})
    return rows


def at_blocks(atsim: dict) -> list[dict]:
    rows = []
    for row in atsim["per_chip_T"]:
        if row["T_column_mean_C"] is None:
            continue
        rows.append({
            "name": row["name"],
            "mean_C": float(row["T_column_mean_C"]),
            "max_C": float(row["T_column_max_C"]),
        })
    return rows


def cost_of(ev: dict, wl: float, scales: dict) -> float:
    """温度、梯度、线长各三分之一。梯度是温度匀不匀的指标，不并进温度项。

    三项都除以起点的值，所以量纲不同也能加。
    """
    return (
        (ev["Tmax_C"] - cl.T_AMB_C) / scales["dT"]
        + ev["grad_max"] / scales["g0"]
        + wl / scales["wl0"]
    ) / 3.0


def layout_mix(step: int, total: int) -> tuple[float, float, float, float]:
    """五项里改芯片固定 0.10，剩下 0.90 当作 1 分给布局四项。

    旋转在这 0.90 里固定 0.10。k 从 0.10 线性升到 0.70，
    剩下的 0.90 - k 由大跳和交换位置平分。
    返回的四项是占全部抽签的比例，加起来是 0.90。
    """
    done = step / max(total, 1)
    k = 0.10 + 0.60 * done
    rest = 0.90 - k
    share = 0.90
    return (share * k, share * rest / 2, share * rest / 2, share * 0.10)


def move_of(dies: list[dict], rng: random.Random, mix: tuple[float, float, float, float]) -> tuple[list[dict], str, str]:
    """四种动作：微调、大跳、交换位置、旋转。比例按 (微调, 大跳, 交换, 旋转)。"""
    trial = deepcopy(dies)
    roll = rng.random()
    b_nudge, b_jump, b_swap = mix[0], mix[0] + mix[1], mix[0] + mix[1] + mix[2]
    i = rng.randrange(len(trial))
    if roll < b_nudge:
        step = rng.choice([1, 1, 2, 3, 4])
        dj, di = rng.choice([(1, 0), (-1, 0), (0, 1), (0, -1)])
        trial[i]["j0"] += dj * step
        trial[i]["ib0"] += di * step
        return trial, f"nudge {trial[i]['name']} ({dj * step},{di * step})", "nudge"
    if roll < b_jump:
        trial[i]["j0"] = rng.randrange(hs.N - trial[i]["w"] + 1)
        trial[i]["ib0"] = rng.randrange(hs.N - trial[i]["h"] + 1)
        if rng.random() < 0.4:
            trial[i]["w"], trial[i]["h"] = trial[i]["h"], trial[i]["w"]
            return trial, f"jump+rot {trial[i]['name']}", "jump"
        return trial, f"jump {trial[i]['name']}", "jump"
    if roll < b_swap:
        j = rng.randrange(len(trial) - 1)
        if j >= i:
            j += 1
        trial[i]["j0"], trial[j]["j0"] = trial[j]["j0"], trial[i]["j0"]
        trial[i]["ib0"], trial[j]["ib0"] = trial[j]["ib0"], trial[i]["ib0"]
        return trial, f"swap {trial[i]['name']}<->{trial[j]['name']}", "swap"
    trial[i]["w"], trial[i]["h"] = trial[i]["h"], trial[i]["w"]
    return trial, f"rotate {trial[i]['name']}", "rotate"


def anneal(start: list[dict], cat: dict, k_fill: float, matrix: list[list[int]]) -> tuple[list[dict], list[dict], np.ndarray, np.ndarray]:
    rng_m = random.Random(SEED)
    rng_a = random.Random(SEED + 7919)
    current = deepcopy(start)
    k_si = hw.paint(current, hs.K_SI)
    mat_si, g_si = hs.build_operator(k_si)
    fac_si = splu(mat_si.tocsc())
    base0 = hw.solve(k_si, hw.power_map(current), fac_si, g_si)[0]
    residual = iterate(current, k_fill, base0) if abs(k_fill - hs.K_SI) > 1e-9 else None
    temp0 = base0 if residual is None else base0 + residual
    err0 = 0.0
    cur = cl.metrics(temp0)
    wl0 = wire.wirelength_cells(current, matrix, cl.CW_MM)
    scales = {
        "dT": max(cur["Tmax_C"] - cl.T_AMB_C, 1e-6),
        "g0": max(cur["grad_max"], 1e-9),
        "wl0": max(wl0, 1e-9),
    }
    j_cur = cost_of(cur, wl0, scales)
    best = deepcopy(current)
    best_j = j_cur
    best_temp = temp0
    log = [{
        "step": 0, "op": "init", "accept": True, "J": j_cur,
        "Tmax_C": cur["Tmax_C"], "grad_max": cur["grad_max"],
        "wirelength_mm_per_wire": wl0,
        "superposition_max_abs_K": err0,
        "hardware": [d["name"] for d in current],
    }]
    print(
        f"  init T {cur['Tmax_C']:.2f}  grad {cur['grad_max']:.3f}  "
        f"wl {wl0:.2f}  J {j_cur:.3f}",
        flush=True,
    )
    # 改芯片固定 0.10，其余 0.90 按步数分给布局四项。
    for step in range(1, N_INNER + 1):
        mix = layout_mix(step, N_INNER)
        if rng_m.random() < 0.10:
            trial, op = cl.swap_hw(current, cat, rng_m)
            move = "改芯片"
        else:
            trial, op, move = move_of(current, rng_m, mix)
        if not sm.legal(trial):
            log.append({"step": step, "op": op, "move": move, "accept": False, "reason": "illegal"})
            continue
        cl.apply_power(trial)
        trial_base = hw.solve(k_si, hw.power_map(trial), fac_si, g_si)[0]
        trial_t = trial_base if residual is None else trial_base + residual
        ev = cl.metrics(trial_t)
        wl_new = wire.wirelength_cells(trial, matrix, cl.CW_MM)
        j_new = cost_of(ev, wl_new, scales)
        # 更好的必收。较差的按步数收：第一步 0.20，最后一步 0，中间线性。
        # SUITE_WORSE_P 给定时改用固定概率，0 就是只收不差于当前的。
        if "SUITE_WORSE_P" in os.environ:
            worse_p = float(os.environ["SUITE_WORSE_P"])
        else:
            worse_p = 0.20 * (1.0 - step / N_INNER)
        if j_new <= j_cur:
            take, ap = True, 1.0
        else:
            take, ap = rng_a.random() < worse_p, worse_p
        rec = {
            "step": step, "op": op, "move": move, "accept": take, "ap": ap,
            "J": j_new, "Tmax_C": ev["Tmax_C"], "grad_max": ev["grad_max"],
            "wirelength_mm_per_wire": wl_new,
            "hardware": [d["name"] for d in trial],
        }
        if take:
            current = trial
            j_cur = j_new
            if j_new < best_j - 1e-12:
                best = deepcopy(trial)
                best_j = j_new
                best_temp = trial_t
                rec["new_best"] = True
                print(
                    f"  {step:02d} BEST {op}  T {ev['Tmax_C']:.2f}  "
                    f"wl {wl_new:.2f}  J {j_new:.3f}",
                    flush=True,
                )
        log.append(rec)
        del trial_t
    cl.apply_power(best)
    best_base = hw.solve(k_si, hw.power_map(best), fac_si, g_si)[0]
    best_temp = best_base if residual is None else best_base + residual
    err_b = 0.0
    if residual is not None:
        err_b = float(np.max(np.abs(iterate(best, k_fill, best_base) - residual)))
    log.append({"step": "best_check", "fill_residual_drift_K": err_b})
    print(f"  best fill residual drift {err_b:.3e} K", flush=True)
    return best, log, temp0, best_temp


def atsim_case(dies: list[dict], path: Path, fill: str) -> None:
    cl.apply_power(dies)
    payload = {
        "name": path.stem,
        "description": "closed_loop_suite",
        "chips": [
            {
                "name": d["name"],
                "x_mm": d["j0"] * cl.CW_MM,
                "y_mm": d["ib0"] * cl.CW_MM,
                "w_mm": d["w"] * cl.CW_MM,
                "h_mm": d["h"] * cl.CW_MM,
                "t_um": 150.0,
                "power_W": d["power_W"],
            }
            for d in dies
        ],
        "package": {
            "tim_um": 20, "t_spreader_mm": 1.0, "t_sink_mm": 6.9,
            "s_spreader_mm": 96.0, "s_sink_mm": 96.0, "frame_mm": 72.0,
            "tim_cover": "frame",
            "t_substrate_mm": 1.0, "t_interposer_um": 100.0,
            "r_convec_K_per_W": 0.0390625, "ambient_C": 45,
            "pitch_die_mm": 2.0, "pitch_coarse_mm": 8.0,
            "nz_die": 3, "nz_spreader": 3, "nz_sink": 4, "nz_tim": 1,
            "nz_substrate": 3, "nz_interposer": 2,
            "fill_die": True, "fill_xy": "chips",
            "fill_material": fill if fill in ("epoxy", "underfill", "air") else "silicon",
        },
    }
    if fill in ("epoxy", "underfill"):
        payload["materials"] = {"epoxy": EPOXY}
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def run_atsim(dies: list[dict], fill: str, case_path: Path, out_path: Path) -> dict:
    wait_ram()
    atsim_case(dies, case_path, fill)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "2"
    env["OPENBLAS_NUM_THREADS"] = "2"
    env["MKL_NUM_THREADS"] = "2"
    env["NUMEXPR_NUM_THREADS"] = "2"
    print(f"  ATSim {case_path.stem}", flush=True)
    with out_path.with_suffix(".log").open("w", encoding="utf-8") as log:
        proc = subprocess.run(
            [str(FENICS), str(SRC / "run_closed_loop_multidie.py"),
             "--atsim", str(case_path), str(out_path)],
            stdout=log, stderr=subprocess.STDOUT, env=env, check=False,
        )
    if proc.returncode != 0:
        tail = out_path.with_suffix(".log").read_text(encoding="utf-8", errors="replace")[-2000:]
        raise RuntimeError(f"ATSim {case_path.stem} rc={proc.returncode}\n{tail}")
    gc.collect()
    return json.loads(out_path.read_text(encoding="utf-8"))


def side(die: dict, scale: float) -> dict:
    freq = None if die.get("f_MHz") is None else float(die["f_MHz"]) * scale
    ops = SIM_INSTS.get(die["part"])
    ratio = None
    if ops is not None and freq:
        ratio = ops / freq
    return {
        "part": die["part"],
        "ops": ops,
        "f_MHz": freq,
        "ops_per_MHz": ratio,
        "cores": die.get("cores"),
        "s": scale,
    }


def throughput(init: list[dict], best: list[dict], freq: list[dict]) -> dict:
    best_by = {d["slot"]: d for d in best}
    freq_by = {d["slot"]: d for d in freq}
    slots = []
    for die in init:
        if die["role"] not in ("CPU", "GPU"):
            continue
        slots.append({
            "slot": die["slot"],
            "role": die["role"],
            "init": side(die, 1.0),
            "best": side(best_by[die["slot"]], 1.0),
            "freq": side(freq_by[die["slot"]], float(freq_by[die["slot"]]["s"])),
        })
    cims = []
    for die in init:
        if die["role"] != "CIM":
            continue
        cims.append({
            "name": die["part"],
            "n_computes": N_COMPUTES.get(die["part"]),
            "f_MHz": None,
            "note": "有计算次数，没有主频，不算操作数/频率",
        })
    return {
        "note": (
            "GPU 的操作数是 backprop 那一次 kernel 的指令数，频率是该卡主频乘倍率。"
            "操作数/频率的单位是次/MHz。同一条 kernel 指令数几乎不变，主频升高这个比值下降，"
            "它不是每秒完成的工作。CPU 的 McPAT 模板没有这次任务的操作数，不编。"
            "CIM、CPU、GPU 不加总。"
        ),
        "slots": slots,
        "cim": cims,
    }


def stage_row(py_rows, hs_rows, atsim, py_grad) -> dict:
    at_rows = at_blocks(atsim) if atsim else []
    return {
        "python": hottest(py_rows),
        "hotspot": hottest(hs_rows) if hs_rows else None,
        "atsim": hottest(at_rows) if at_rows else None,
        "python_grad_K_per_mm": py_grad,
        "atsim_grad_K_per_mm": None if atsim is None else atsim["grad_max_K_per_mm"],
        "stress_MPa": None if atsim is None else atsim["vm_normal_chip_max_MPa"],
        "stress_name": None if atsim is None else atsim["vm_normal_chip_name"],
        "gate_die_cell_max_C": None if atsim is None else atsim["die_cell_max_C"],
        "python_per_chip": py_rows,
        "hotspot_per_chip": hs_rows,
        "atsim_per_chip": at_rows,
    }


def climb(dies: list[dict], rises: dict[str, np.ndarray], kcal: float, mode: str) -> dict:
    slots = [d["slot"] for d in dies if d["role"] in ("CPU", "GPU")]
    scales = {slot: 1.0 for slot in slots}
    log = []

    def guide(sc: dict[str, float]) -> dict:
        ev = cl.metrics(field_at(dies, rises, sc))
        blocks = py_blocks(field_at(dies, rises, sc), dies)
        hot = hottest(blocks)
        est = cl.T_AMB_C + kcal * (hot["C"] - cl.T_AMB_C)
        return {"T_block_C": hot["C"], "T_est_C": est, "name": hot["name"],
                "grad_max": ev["grad_max"], "ok": est < cl.T_GATE}

    if mode == "outside" and guide(scales)["ok"]:
        log.append({"op": "stress_outside_hold", "scales": dict(scales), **guide(scales)})
        return {"scales": scales, "log": log, "guide": guide(scales), "held": True}

    for _ in range(30):
        g0 = guide(scales)
        if g0["ok"]:
            break
        best = None
        for slot in slots:
            if scales[slot] <= cl.S_LO + 1e-9:
                continue
            trial = dict(scales)
            trial[slot] = round(scales[slot] - cl.S_STEP, 1)
            g = guide(trial)
            if best is None or g["T_est_C"] < best[0]:
                best = (g["T_est_C"], slot, trial, g)
        if best is None:
            break
        scales = best[2]
        log.append({"op": f"drop {best[1]}", "scales": dict(scales), **best[3]})
        print(f"  freq drop {best[1]} -> {scales[best[1]]:.1f}  T_est {best[0]:.2f}", flush=True)

    bumped = 0
    while bumped < 80:
        pick = None
        for slot in slots:
            if scales[slot] >= cl.S_HI - 1e-9:
                continue
            trial = dict(scales)
            trial[slot] = round(scales[slot] + cl.S_STEP, 1)
            g = guide(trial)
            if not g["ok"]:
                continue
            if pick is None or g["T_est_C"] < pick[0]:
                pick = (g["T_est_C"], slot, trial, g)
        if pick is None:
            break
        scales = pick[2]
        bumped += 1
        log.append({"op": f"up {pick[1]}", "scales": dict(scales), **pick[3]})
    g = guide(scales)
    print(f"  freq sum_s {sum(scales.values()):.1f}  T_est {g['T_est_C']:.2f}  bumps {bumped}", flush=True)
    return {"scales": scales, "log": log, "guide": g, "held": False, "bumps": bumped}


STRESS_OVER = 20.0
STRESS_TRIES = 6
STRESS_WINS = 3


def stress_fix(dies: list[dict], cat: dict, k_fill: float, fill: str, atsim: dict, out: Path) -> dict:
    """频率已经调到极限。应力仍在 220 以内才做这一轮，超过 220 直接保留。

    最多 6 次尝试、3 次有效下降。还过不了就保留，标不达标。
    """
    log = []
    current = deepcopy(dies)
    result = atsim
    tries = 0
    wins = 0

    def over(vm: float) -> bool:
        return vm >= cl.S_GATE

    def record(op: str, take: bool) -> None:
        log.append({
            "try": tries, "op": op, "accept": take,
            "vm_MPa": result["vm_normal_chip_max_MPa"],
            "vm_name": result["vm_normal_chip_name"],
            "die_cell_max_C": result["die_cell_max_C"],
            "scales": cl.scales_of(current),
        })

    if result["vm_normal_chip_max_MPa"] > cl.S_GATE + STRESS_OVER:
        return {
            "dies": current, "atsim": result, "log": [], "tries": 0, "wins": 0,
            "passed": False, "status": "应力超过 220，不退火，保留作参考",
            "rule": "频率调到极限后应力仍超 220，不再动布局。",
        }

    rng_m = random.Random(SEED + 17)
    rng_a = random.Random(SEED + 23)
    while over(result["vm_normal_chip_max_MPa"]) and tries < STRESS_TRIES and wins < STRESS_WINS:
        roll = rng_m.random()
        if roll < 0.60:
            trial, op = cl.nudge(current, rng_m)
        elif roll < 0.80:
            trial, op = cl.rotate(current, rng_m)
        else:
            trial, op = cl.swap_hw(current, cat, rng_m)
        if not sm.legal(trial):
            continue
        tries += 1
        trial_at = run_atsim(trial, fill, out / f"fix_{tries}.json", out / f"fix_{tries}_result.json")
        take = trial_at["vm_normal_chip_max_MPa"] < result["vm_normal_chip_max_MPa"] - 1e-9
        if take:
            current = trial
            result = trial_at
            wins += 1
        record(op, take)
        print(f"  stress sa {tries} {op}  vm {trial_at['vm_normal_chip_max_MPa']:.1f}  take {take}", flush=True)

    passed = result["die_cell_max_C"] < cl.T_GATE and not over(result["vm_normal_chip_max_MPa"])
    return {
        "dies": current, "atsim": result, "log": log, "tries": tries, "wins": wins,
        "passed": passed,
        "status": "达标" if passed else "不达标，保留作参考",
        "rule": "先降频率，仍超再做应力退火。最多 6 次尝试、3 次有效下降。",
    }


def save_maps(out: Path, fields: list[tuple[str, np.ndarray, list[dict]]]) -> None:
    style()
    arrays = [t - 273.15 for _, t, _ in fields]
    vmin = min(float(a.min()) for a in arrays)
    vmax = max(float(a.max()) for a in arrays)
    fig, axes = plt.subplots(1, len(fields), figsize=(4.2 * len(fields), 4.4), squeeze=False)
    for ax, (title, _, dies), arr in zip(axes[0], fields, arrays):
        show = arr[::-1, :]
        im = ax.imshow(show, origin="lower", vmin=vmin, vmax=vmax, cmap="inferno")
        for die in dies:
            ax.add_patch(plt.Rectangle(
                (die["j0"] - 0.5, die["ib0"] - 0.5), die["w"], die["h"],
                fill=False, ec="white", lw=0.6,
            ))
        ax.set_title(title, loc="left")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(im, ax=axes[0].tolist(), fraction=0.03, label="C")
    fig.suptitle("Python 芯粒层，含本次填充和外伸环。左下角是 (0, 0)。", fontsize=10)
    fig.savefig(out, dpi=120)
    plt.close(fig)


def save_bars(path: Path, report: dict) -> None:
    style()
    labels = ["初始", "最优", "调频后"]
    keys = ["init", "best", "freq"]
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8))
    x = np.arange(3)
    w = 0.25
    for shift, solver, color in (
        (-w, "python", "#b85c38"),
        (0.0, "hotspot", "#1f4e79"),
        (w, "atsim", "#2f6f4e"),
    ):
        vals = []
        for key in keys:
            cell = report["compare"][key][solver]
            vals.append(np.nan if cell is None else cell["C"])
        axes[0].bar(x + shift, vals, w, label=solver, color=color)
    axes[0].set_xticks(x, labels)
    axes[0].set_ylabel("C")
    axes[0].set_title("最热芯片的块平均", loc="left")
    axes[0].legend(frameon=False, fontsize=8)
    stress = [report["compare"][k]["stress_MPa"] or np.nan for k in keys]
    axes[1].bar(labels, stress, color="#1f4e79")
    axes[1].axhline(200, color="#888", lw=0.8)
    axes[1].set_ylabel("MPa")
    axes[1].set_title("内部应力", loc="left")
    pg = [report["compare"][k]["python_grad_K_per_mm"] for k in keys]
    ag = [report["compare"][k]["atsim_grad_K_per_mm"] or np.nan for k in keys]
    axes[2].plot(labels, pg, "o-", color="#b85c38", label="Python")
    axes[2].plot(labels, ag, "o-", color="#2f6f4e", label="ATSim")
    axes[2].set_ylabel("K/mm")
    axes[2].set_title("温度梯度", loc="left")
    axes[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def save_ops(path: Path, slots: list[dict]) -> None:
    style()
    gpus = [s for s in slots if s["role"] == "GPU"]
    cpus = [s for s in slots if s["role"] == "CPU"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8))
    if gpus:
        x = np.arange(len(gpus))
        w = 0.25
        for shift, key, label, color in (
            (-w, "init", "初始", "#c5c1b9"),
            (0.0, "best", "最优", "#1f4e79"),
            (w, "freq", "调频后", "#b85c38"),
        ):
            vals = [s[key]["ops_per_MHz"] or np.nan for s in gpus]
            axes[0].bar(x + shift, vals, w, label=label, color=color)
        axes[0].set_xticks(x, [s["freq"]["part"].replace("GPU_", "") + " " + s["slot"] for s in gpus], rotation=15)
    axes[0].set_ylabel("次/MHz")
    axes[0].set_title("GPU 指令数/主频", loc="left")
    axes[0].legend(frameon=False, fontsize=8)
    if cpus:
        x = np.arange(len(cpus))
        w = 0.25
        for shift, key, label, color in (
            (-w, "init", "初始", "#c5c1b9"),
            (0.0, "best", "最优", "#1f4e79"),
            (w, "freq", "调频后", "#b85c38"),
        ):
            axes[1].bar(x + shift, [s[key]["f_MHz"] or np.nan for s in cpus], w, label=label, color=color)
        axes[1].set_xticks(x, [f"{s['freq']['part']} {s['slot']}" for s in cpus], rotation=15)
    axes[1].set_ylabel("MHz")
    axes[1].set_title("CPU 主频，无操作数", loc="left")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def write_svgs(out: Path) -> None:
    cl.OUT = out
    cl.write_svgs()


def run_case(name: str, fill: str, kind: str, cat: dict) -> None:
    out = SUITE / name
    if (out / "report.json").exists() and not os.environ.get("SUITE_FORCE"):
        print(f"skip {name}", flush=True)
        return
    out.mkdir(parents=True, exist_ok=True)
    (out / "layouts").mkdir(exist_ok=True)
    (out / "fields").mkdir(exist_ok=True)
    (out / "figures").mkdir(exist_ok=True)
    (out / "atsim").mkdir(exist_ok=True)
    print(f"\n== {name}  fill={fill}  steps={N_INNER} ==", flush=True)
    spec = fill_of(fill)
    k_fill = spec["k"]
    init_file = os.environ.get("SUITE_INIT_FILE")
    if os.environ.get("SUITE_SHELF"):
        start, init_name = shelf(cat, kind), "shelf"
    elif init_file:
        start = json.loads(Path(init_file).read_text(encoding="utf-8"))
        for die in start:
            die.setdefault("s", 1.0)
        cl.apply_power(start)
        init_name = f"file:{Path(init_file).stem}"
        print(f"  init from {init_file}", flush=True)
    else:
        start, init_name, init_score = initial_layout(members(cat, kind), k_fill)
        print(f"  chosen init {init_name}  score {init_score['score']:.3f}", flush=True)
    print(f"  init from {init_name}", flush=True)
    matrix = wire.assumed_matrix(start)
    t0 = time.time()
    best, log, t_init, t_best = anneal(start, cat, k_fill, matrix)
    (out / "trace.json").write_text(json.dumps(log, ensure_ascii=False) + "\n", encoding="utf-8")
    (out / "layouts" / "init.json").write_text(
        json.dumps(cl.pack(start), indent=2) + "\n", encoding="utf-8")
    (out / "layouts" / "best.json").write_text(
        json.dumps(cl.pack(best), indent=2) + "\n", encoding="utf-8")
    np.save(out / "fields" / "init.npy", t_init)
    np.save(out / "fields" / "best.npy", t_best)

    hs_init = cl.hotspot_block_max(start, out / "hotspot_init")
    hs_best = cl.hotspot_block_max(best, out / "hotspot_best")
    py_init = py_blocks(t_init, start)
    py_best = py_blocks(t_best, best)
    at_init = at_best = at_freq = None
    stress = None
    tuned = deepcopy(best)
    freq_info = {"scales": {d["slot"]: 1.0 for d in best if d["role"] in ("CPU", "GPU")},
                 "log": [], "guide": None, "held": True, "bumps": 0}
    walks = []
    if not SKIP_ATSIM:
        at_init = run_atsim(start, fill, out / "atsim" / "init_s1.json", out / "atsim" / "init_s1_result.json")
        at_best = run_atsim(best, fill, out / "atsim" / "best_s1.json", out / "atsim" / "best_s1_result.json")
        mode = cl.classify(at_best)
        py_hot = hottest(py_best)["C"]
        at_hot = hottest(at_blocks(at_best))["C"]
        kcal = (at_hot - cl.T_AMB_C) / max(py_hot - cl.T_AMB_C, 1e-6)
        print(f"  mode {mode}  kcal {kcal:.4f}", flush=True)
        rises = rises_of(best, k_fill)
        freq_info = climb(best, rises, kcal, mode)
        tuned = cl.with_scales(best, freq_info["scales"])
        for attempt in range(12):
            at_freq = run_atsim(
                tuned, fill,
                out / "atsim" / f"freq_{attempt}.json",
                out / "atsim" / f"freq_{attempt}_result.json",
            )
            walks.append({
                "attempt": attempt,
                "scales": cl.scales_of(tuned),
                "die_cell_max_C": at_freq["die_cell_max_C"],
                "vm_normal_chip_max_MPa": at_freq["vm_normal_chip_max_MPa"],
                "block_mean_C": hottest(at_blocks(at_freq))["C"],
            })
            if at_freq["die_cell_max_C"] < cl.T_GATE and at_freq["vm_normal_chip_max_MPa"] < cl.S_GATE:
                break
            if at_freq["die_cell_max_C"] >= cl.T_GATE:
                movable = [d for d in tuned if d["role"] in ("CPU", "GPU") and d["s"] > cl.S_LO + 1e-9]
                if not movable:
                    break
                victim = max(movable, key=lambda d: d["dyn_W"])
                victim["s"] = round(victim["s"] - cl.S_STEP, 1)
                cl.apply_power(tuned)
                print(f"  walkback {victim['name']} s {victim['s']:.1f}", flush=True)
                continue
            break
        del rises
        stress = None
        if at_freq["vm_normal_chip_max_MPa"] >= cl.S_GATE:
            stress = stress_fix(tuned, cat, k_fill, fill, at_freq, out / "atsim" / "stress")
            tuned = stress["dies"]
            at_freq = stress["atsim"]
            print(f"  stress {stress['status']}  tries {stress['tries']}  wins {stress['wins']}", flush=True)
    else:
        kcal = None
        mode = "skipped"
        stress = None
    t_freq = field_at(best, rises_of(best, k_fill), cl.scales_of(tuned))
    np.save(out / "fields" / "freq.npy", t_freq)
    hs_freq = cl.hotspot_block_max(tuned, out / "hotspot_freq")
    (out / "layouts" / "freq.json").write_text(
        json.dumps(cl.pack(tuned), indent=2) + "\n", encoding="utf-8")
    (out / "freq_trace.json").write_text(
        json.dumps(freq_info["log"], ensure_ascii=False) + "\n", encoding="utf-8")

    passed = None
    if at_freq is not None:
        passed = at_freq["die_cell_max_C"] < cl.T_GATE and at_freq["vm_normal_chip_max_MPa"] < cl.S_GATE
    report = {
        "case": name,
        "n_chips": len(start),
        "fill": fill,
        "k_fill": k_fill,
        "epoxy_mechanics": None if fill != "epoxy" else (
            "导热 0.8 是工程里的数。E=25 GPa，nu=0.30，CTE=1.5e-5 借用 ATSim 已有的 organic 卡，"
            "只写在本案例 json，没有改默认材料表。"
        ),
        "package": "72 mm 芯粒框，96 mm 热沉，r_convec 0.0390625",
        "inner_steps": N_INNER,
        "seed": SEED,
        "kcal_python_block_to_atsim_block": kcal,
        "kcal_note": "向导按最优点的 ATSim 块平均温升 / Python 块平均温升来缩。HotSpot 不吃填充，不拿来锚定。",
        "freq_mode_at_best_s1": mode,
        "freq_bumps": freq_info.get("bumps", 0),
        "freq_held_for_stress": freq_info.get("held", False),
        "gate_passed_after_freq": passed,
        "stress_in_anneal": False,
        "stress_after_layout": None if stress is None else {
            "status": stress["status"],
            "tries": stress["tries"],
            "wins": stress["wins"],
            "passed": stress["passed"],
            "rule": stress["rule"],
            "log": stress["log"],
        },
        "walkback": walks,
        "hardware_init": [d["name"] for d in start],
        "hardware_best": [d["name"] for d in best],
        "hbm_W_init": next(d["power_W"] for d in start if d["role"] == "HBM"),
        "hbm_W_best": next(d["power_W"] for d in best if d["role"] == "HBM"),
        "elapsed_s": round(time.time() - t0, 1),
        "compare": {
            "definition": (
                "温度是每颗芯片盖住的格子先求平均，再取最热那颗。"
                "Python 是 1.125 mm 芯粒层。HotSpot 是第 0 层块平均，而且芯片层导热仍是硅，不吃填充。"
                "ATSim 是 2 mm 网格上芯粒层柱的平均。应力是正常芯片内部 von Mises 最大，不是整包最大。"
                "梯度两列网格不同，不能并成一列。门槛仍看 ATSim 芯粒层最热一格和内部应力。"
            ),
            "init": stage_row(py_init, hs_blocks(hs_init["blocks_C"], start), at_init, cl.metrics(t_init)["grad_max"]),
            "best": stage_row(py_best, hs_blocks(hs_best["blocks_C"], best), at_best, cl.metrics(t_best)["grad_max"]),
            "freq": stage_row(py_blocks(t_freq, tuned), hs_blocks(hs_freq["blocks_C"], tuned), at_freq, cl.metrics(t_freq)["grad_max"]),
        },
        "throughput": throughput(start, best, tuned),
        "superposition_init_K": log[0].get("superposition_max_abs_K", 0.0),
        "fill_residual_drift_K": log[-1].get("fill_residual_drift_K", 0.0),
        "wirelength": {
            "formula": "TAP-2.5D 的中心曼哈顿，mm/wire。连接矩阵在退火前定死，初始和最优用同一张。",
            "netlist": "假设的：每颗 CPU、GPU 和 HBM 双向 1024，每颗 CIM 和第一颗 CPU 双向 128。不是仿真给的。",
            "connection_matrix": matrix,
            "order": [d["name"] for d in start],
            "init_mm_per_wire": log[0]["wirelength_mm_per_wire"],
            "best_mm_per_wire": wire.wirelength_cells(best, matrix, cl.CW_MM),
            "cost": "温度、梯度、线长各三分之一，都除以起点的值。梯度不并进温度项。",
        },
    }
    (out / "report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    save_maps(out / "figures" / "python_maps.png", [
        ("初始", t_init, start),
        ("最优", t_best, best),
        ("调频后", t_freq, tuned),
    ])
    save_bars(out / "figures" / "compare_bars.png", report)
    save_ops(out / "figures" / "ops_per_mhz.png", report["throughput"]["slots"])
    write_svgs(out)
    scale_path = out / "figures" / "svg_scale.json"
    if scale_path.exists():
        meta = json.loads(scale_path.read_text(encoding="utf-8"))
        meta["hotspot_fill"] = "HotSpot 芯片层导热仍是硅，这三张图不包含填充。填充在 Python 图和 ATSim 里。"
        meta["fill_this_case"] = fill
        scale_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"  wrote {out / 'report.json'}  {report['elapsed_s']} s", flush=True)
    del t_init, t_best, t_freq
    gc.collect()


def cell(report: dict, stage: str, solver: str) -> str:
    item = report["compare"][stage][solver]
    if item is None:
        return "-"
    return f"{item['name']} {item['C']:.1f}"


def summarize() -> None:
    style()
    rows = []
    old = ROOT / "热仿真流水线" / "runs" / "closed_loop_multidie" / "report.json"
    if old.exists():
        src = json.loads(old.read_text(encoding="utf-8"))
        ba = src["block_average"]
        made = {"case": "n11_silicon_已有", "fill": "silicon", "rerun": False, "compare": {}}
        for stage in ("init", "best", "freq"):
            hot = ba["hottest"][stage]
            at = src[stage]["atsim"]
            made["compare"][stage] = {
                "python": {"name": hot["python_hottest_block_mean"]["name"], "C": hot["python_hottest_block_mean"]["C"]},
                "hotspot": {"name": hot["hotspot_hottest_block_mean"]["name"], "C": hot["hotspot_hottest_block_mean"]["C"]},
                "atsim": {"name": hot["atsim_hottest_block_mean"]["name"], "C": hot["atsim_hottest_block_mean"]["C"]},
                "stress_MPa": at["vm_normal_chip_max_MPa"],
                "stress_name": at["vm_normal_chip_name"],
                "python_grad_K_per_mm": src[stage]["python_grad_K_per_mm"],
                "atsim_grad_K_per_mm": at["grad_max_K_per_mm"],
            }
        rows.append(made)
    for path in sorted(SUITE.glob("*/report.json")):
        rows.append(json.loads(path.read_text(encoding="utf-8")))
    (SUITE / "index.json").write_text(
        json.dumps([{
            "case": r["case"],
            "fill": r.get("fill"),
            "atsim_init_C": r["compare"]["init"]["atsim"]["C"] if r["compare"]["init"]["atsim"] else None,
            "atsim_best_C": r["compare"]["best"]["atsim"]["C"] if r["compare"]["best"]["atsim"] else None,
            "stress_init_MPa": r["compare"]["init"]["stress_MPa"],
            "stress_best_MPa": r["compare"]["best"]["stress_MPa"],
            "python_grad_init": r["compare"]["init"]["python_grad_K_per_mm"],
            "python_grad_best": r["compare"]["best"]["python_grad_K_per_mm"],
        } for r in rows], indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    def num(report: dict, stage: str, key: str):
        if key == "atsim":
            item = report["compare"][stage]["atsim"]
            return np.nan if item is None else item["C"]
        value = report["compare"][stage][key]
        return np.nan if value is None else value

    if len(rows) < 2:
        return
    labels = [r["case"].replace("_", "\n") for r in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2))
    axes[0].bar(x - 0.18, [num(r, "init", "atsim") for r in rows], 0.36, label="初始", color="#c5c1b9")
    axes[0].bar(x + 0.18, [num(r, "best", "atsim") for r in rows], 0.36, label="最优", color="#1f4e79")
    axes[0].axhline(90, color="#888", lw=0.8)
    axes[0].set_ylabel("C")
    axes[0].set_title("ATSim 最热芯片块平均", loc="left")
    axes[1].bar(x - 0.18, [num(r, "init", "stress_MPa") for r in rows], 0.36, label="初始", color="#c5c1b9")
    axes[1].bar(x + 0.18, [num(r, "best", "stress_MPa") for r in rows], 0.36, label="最优", color="#1f4e79")
    axes[1].axhline(200, color="#888", lw=0.8)
    axes[1].set_ylabel("MPa")
    axes[1].set_title("内部应力", loc="left")
    axes[2].bar(x - 0.18, [r["compare"]["init"]["python_grad_K_per_mm"] for r in rows], 0.36, label="初始", color="#c5c1b9")
    axes[2].bar(x + 0.18, [r["compare"]["best"]["python_grad_K_per_mm"] for r in rows], 0.36, label="最优", color="#b85c38")
    axes[2].set_ylabel("K/mm")
    axes[2].set_title("Python 梯度", loc="left")
    for ax in axes:
        ax.set_xticks(x, labels, fontsize=7)
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(SUITE / "suite_compare.png", dpi=120)
    plt.close(fig)
    print(f"suite index {SUITE / 'index.json'}", flush=True)


CASES = (
    ("n5_silicon", "silicon", "n5"),
    ("n5_epoxy", "epoxy", "n5"),
    ("n6_silicon", "silicon", "n6"),
    ("n6_epoxy", "epoxy", "n6"),
    ("n7_silicon", "silicon", "n7"),
    ("n7_epoxy", "epoxy", "n7"),
    ("n8_silicon", "silicon", "n8"),
    ("n8_epoxy", "epoxy", "n8"),
    ("n9_silicon", "silicon", "n9"),
    ("n9_epoxy", "epoxy", "n9"),
    ("n11_epoxy", "epoxy", "n11"),
    ("n11_silicon", "silicon", "n11"),
)


def main() -> None:
    SUITE.mkdir(parents=True, exist_ok=True)
    cat = cl.catalog()
    failed = []
    for name, fill, kind in CASES:
        if ONLY and name not in ONLY:
            continue
        try:
            run_case(name, fill, kind, cat)
        except Exception:
            failed.append(name)
            (SUITE / name).mkdir(parents=True, exist_ok=True)
            (SUITE / name / "error.txt").write_text(traceback.format_exc(), encoding="utf-8")
            print(traceback.format_exc(), flush=True)
        gc.collect()
    summarize()
    if failed:
        raise SystemExit(f"失败: {failed}")


if __name__ == "__main__":
    main()
