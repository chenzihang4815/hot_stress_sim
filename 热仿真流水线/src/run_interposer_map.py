#!/usr/bin/env python3
"""同一布局：芯片下面加不加 100 µm 硅中介层。

python 四层模型的底是绝热的。加上的中介层也只在底面绝热，
热仍然从热沉上表面离开。缝的填充不改。
ATSim 用它自己的封装：有机基板两边都在，只把中介层厚度从 100 µm 改成 0。
绝对温度不混用。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import run_factor_controls as fc  # noqa: E402
import run_fill_superposition as hs  # noqa: E402
import run_hw_area_field as hw  # noqa: E402

OUT = ROOT / "热仿真流水线" / "runs" / "interposer_map"
N = hs.N
T_AMB = hs.T_AMB
T_IP = 100e-6
K_IP = 130.0


def die_max(T: np.ndarray, die: dict) -> float:
    return float(max(T[i, j] for i, j in hw.cells(die)) - 273.15)


def solve_edge(k_chip: np.ndarray, power: np.ndarray) -> np.ndarray:
    """冷边、底面绝热的五层场，芯片层，单位 K。"""
    return solve_with_interposer(k_chip, power, k_tim=4.0, g_bottom=0.0, edge_to_ambient=True)[1]


def stale_area() -> None:
    """四颗里只换一颗的面积。其余三颗沿用旧材料图上的场，只重解被换的那颗。

    材料图是新足迹：旧硅格改回环氧，新足迹涂成硅。
    对照是同一张新材料图上四颗轮流加功耗再相加，应当只剩舍入。
    冷边、底面绝热、缝是环氧。python 场，不是 HotSpot，也不是 ATSim。
    CPU 功耗是各模板自己的 McPAT 380 K 总功耗，GPU 是各自的 backprop 平均。
    """
    from copy import deepcopy

    import run_sa_hw_smooth as sm

    OUT.mkdir(parents=True, exist_ok=True)
    dies, cat = fc.base_dies()
    if not sm.legal(dies):
        raise SystemExit("起点不合法")

    def swap_role(src, idx, name, pool):
        part = next(c for c in pool if c["name"] == name)
        out = deepcopy(src)
        j0, ib0 = out[idx]["j0"], out[idx]["ib0"]
        out[idx] = sm.clonable(part)
        out[idx]["j0"], out[idx]["ib0"] = j0, ib0
        return out

    # 更小的 CPU：面积缩到原来的一半和四分之一，功耗密度按 Penryn 的格子均分保留。
    pen = dies[0]
    small = []
    for frac, tag in ((0.5, "half"), (0.25, "quarter")):
        trial = deepcopy(dies)
        d = trial[0]
        cells_n = max(1, int(round(d["w"] * d["h"] * frac)))
        w = max(1, int(round(d["w"] * frac ** 0.5)))
        h = max(1, int(round(cells_n / w)))
        d["w"], d["h"] = w, h
        d["name"] = f"Penryn_{tag}"
        d["power_W"] = pen["power_W"] * (w * h) / (pen["w"] * pen["h"])
        d["area_mm2"] = pen["area_mm2"] * (w * h) / (pen["w"] * pen["h"])
        d["grid_area_mm2"] = w * h * (hs.CW * 1e3) ** 2
        small.append((f"cpu_{tag}", trial, 0))

    schemes = [
        ("cpu_T2", swap_role(dies, 0, "T2", cat["CPU"]), 0),
        ("cpu_Xeon", swap_role(dies, 0, "Xeon", cat["CPU"]), 0),
        ("gpu_2060S", swap_role(dies, 1, "GPU_2060S", cat["GPU"]), 1),
        ("gpu_TITANX", swap_role(dies, 1, "GPU_TITANX", cat["GPU"]), 1),
        *small,
    ]

    k_old = hw.paint(dies, 0.8)
    p_old = hw.power_map(dies)
    kept = []
    for i in range(len(dies)):
        one = np.zeros_like(p_old)
        cc = hw.cells(dies[i])
        for r, c in cc:
            one[r, c] = p_old[r, c]
        kept.append(solve_edge(k_old, one))
    check = sum(p - T_AMB for p in kept) + T_AMB
    truth0 = solve_edge(k_old, p_old)
    print(f"起点自检 {np.max(np.abs(check - truth0)):.2e} K  Tmax {truth0.max()-273.15:.3f}", flush=True)

    rows = []
    for name, trial, idx in schemes:
        ok = sm.legal(trial)
        if not ok:
            rows.append({"scheme": name, "legal": False})
            print(f"{name} 不合法", flush=True)
            continue
        k_new = hw.paint(trial, 0.8)
        p_new = hw.power_map(trial)
        truth = solve_edge(k_new, p_new)
        parts = list(kept)
        one = np.zeros_like(p_new)
        for r, c in hw.cells(trial[idx]):
            one[r, c] = p_new[r, c]
        parts[idx] = solve_edge(k_new, one)
        stale = sum(p - T_AMB for p in parts) + T_AMB
        fresh = []
        for i in range(len(trial)):
            q = np.zeros_like(p_new)
            for r, c in hw.cells(trial[i]):
                q[r, c] = p_new[r, c]
            fresh.append(solve_edge(k_new, q))
        frozen = sum(p - T_AMB for p in fresh) + T_AMB
        d_stale = stale - truth
        d_frozen = frozen - truth
        old_d, new_d = dies[idx], trial[idx]
        rec = {
            "scheme": name,
            "legal": True,
            "changed": old_d["name"] + " -> " + new_d["name"],
            "old_grid_mm2": old_d["w"] * old_d["h"] * (hs.CW * 1e3) ** 2,
            "new_grid_mm2": new_d["w"] * new_d["h"] * (hs.CW * 1e3) ** 2,
            "old_power_W": old_d["power_W"],
            "new_power_W": new_d["power_W"],
            "cells_changed": len(set(hw.cells(old_d)) ^ set(hw.cells(new_d))),
            "stale_max_K": float(np.max(np.abs(d_stale))),
            "stale_rmse_K": float(np.sqrt(np.mean(d_stale ** 2))),
            "frozen_max_K": float(np.max(np.abs(d_frozen))),
            "python_Tmax_C": float(truth.max() - 273.15),
            "die_C": {d["name"]: die_max(truth, d) for d in trial},
        }
        rows.append(rec)
        print(f"{name:12} {rec['changed']:28} {rec['old_grid_mm2']:7.1f} -> {rec['new_grid_mm2']:7.1f} mm2  "
              f"cells {rec['cells_changed']:4d}  stale {rec['stale_max_K']:.4f} K  "
              f"frozen {rec['frozen_max_K']:.1e}", flush=True)

    payload = {
        "note": (
            "默认：其余三颗的场不改，是在起点材料图上单独解出来的。"
            "被换的那一颗按新足迹重解：新硅、旧硅改回环氧 0.8。"
            "中介层 100 um 硅，最外一圈 45 C，底面绝热，TIM k=4。"
            "stale 是这个做法对真解的误差。frozen 是同一张新材料图上四颗分开再相加，应当接近 0。"
            "python 场，不是 HotSpot，也不是 ATSim。"
            "CPU 功耗是各模板 McPAT 380 K 总功耗，不是同一趟 matmul 重跑；"
            "GPU 功耗是各自的 backprop 平均。half 和 quarter 是把 Penryn 按格子等比缩小。"
        ),
        "start": [sm.dump_die(d) for d in dies],
        "rows": rows,
    }
    (OUT / "stale_area.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"写到 {OUT / 'stale_area.json'}", flush=True)


def stale_timing() -> None:
    """三种做法的误差和墙钟，同一批换面积方案。

    只换一颗：矩阵按新材料图分解一次，只解被换那颗的右端，其余三颗用旧场。
    四颗都换：同一张新材料图分解一次，四颗各解一次再相加。
    一次求解：新材料图分解一次，四颗功耗一起放进一个右端。
    分解只算一次，因为材料图定了矩阵就定了。时间不含画图、不含读盘。
    """
    import time
    from copy import deepcopy

    import run_sa_hw_smooth as sm

    dies, cat = fc.base_dies()

    def swap_role(src, idx, name, pool):
        part = next(c for c in pool if c["name"] == name)
        out = deepcopy(src)
        j0, ib0 = out[idx]["j0"], out[idx]["ib0"]
        out[idx] = sm.clonable(part)
        out[idx]["j0"], out[idx]["ib0"] = j0, ib0
        return out

    pen = dies[0]
    small = []
    for frac, tag in ((0.5, "half"), (0.25, "quarter")):
        trial = deepcopy(dies)
        d = trial[0]
        cells_n = max(1, int(round(d["w"] * d["h"] * frac)))
        w = max(1, int(round(d["w"] * frac ** 0.5)))
        h = max(1, int(round(cells_n / w)))
        d["w"], d["h"] = w, h
        d["name"] = f"Penryn_{tag}"
        d["power_W"] = pen["power_W"] * (w * h) / (pen["w"] * pen["h"])
        d["area_mm2"] = pen["area_mm2"] * (w * h) / (pen["w"] * pen["h"])
        small.append((f"cpu_{tag}", trial, 0))
    schemes = [
        ("cpu_T2", swap_role(dies, 0, "T2", cat["CPU"]), 0),
        ("cpu_Xeon", swap_role(dies, 0, "Xeon", cat["CPU"]), 0),
        ("gpu_2060S", swap_role(dies, 1, "GPU_2060S", cat["GPU"]), 1),
        ("gpu_TITANX", swap_role(dies, 1, "GPU_TITANX", cat["GPU"]), 1),
        *small,
    ]

    def one_power(layout, idx):
        p = np.zeros((N, N))
        share_map = hw.power_map(layout)
        for r, c in hw.cells(layout[idx]):
            p[r, c] = share_map[r, c]
        return p

    # 旧场先备好，不计入下面三种做法的时间：那是上一次布局已经付过的。
    fac_old, g_old, t_old = factor_edge(hw.paint(dies, 0.8))
    kept = [apply_edge(fac_old, g_old, one_power(dies, i)) for i in range(4)]
    del fac_old

    rows = []
    for name, trial, idx in schemes:
        if not sm.legal(trial):
            rows.append({"scheme": name, "legal": False})
            continue
        k_new = hw.paint(trial, 0.8)
        fac, g, t_factor = factor_edge(k_new)

        t0 = time.perf_counter()
        parts = list(kept)
        parts[idx] = apply_edge(fac, g, one_power(trial, idx))
        stale = sum(p - T_AMB for p in parts) + T_AMB
        t_stale = time.perf_counter() - t0

        t0 = time.perf_counter()
        fresh = [apply_edge(fac, g, one_power(trial, i)) for i in range(4)]
        frozen = sum(p - T_AMB for p in fresh) + T_AMB
        t_frozen = time.perf_counter() - t0

        t0 = time.perf_counter()
        direct = apply_edge(fac, g, hw.power_map(trial))
        t_direct = time.perf_counter() - t0

        d_stale = stale - direct
        d_frozen = frozen - direct
        old_d, new_d = dies[idx], trial[idx]
        rec = {
            "scheme": name,
            "changed": old_d["name"] + " -> " + new_d["name"],
            "old_grid_mm2": old_d["w"] * old_d["h"] * (hs.CW * 1e3) ** 2,
            "new_grid_mm2": new_d["w"] * new_d["h"] * (hs.CW * 1e3) ** 2,
            "factor_s": t_factor,
            "stale_solve_s": t_stale,
            "stale_total_s": t_factor + t_stale,
            "stale_max_K": float(np.max(np.abs(d_stale))),
            "frozen_solve_s": t_frozen,
            "frozen_total_s": t_factor + t_frozen,
            "frozen_max_K": float(np.max(np.abs(d_frozen))),
            "direct_solve_s": t_direct,
            "direct_total_s": t_factor + t_direct,
        }
        rows.append(rec)
        print(f"{name:12} stale {rec['stale_max_K']:.4f} K  {rec['stale_total_s']*1e3:7.1f} ms   "
              f"frozen {rec['frozen_max_K']:.1e} K  {rec['frozen_total_s']*1e3:7.1f} ms   "
              f"direct {rec['direct_total_s']*1e3:7.1f} ms   factor {t_factor*1e3:.1f} ms", flush=True)
        del fac

    payload = {
        "note": (
            "三种做法对着同一次换面积。direct 是新材料图上四颗功耗一次求解，当作真解，误差为 0。"
            "stale 只重解被换的那颗，其余三颗是上一次材料图上的场，时间不含准备旧场。"
            "frozen 是新材料图上四颗分开求解再相加。"
            "factor 是装配加分解，三种做法各付一次，因为材料图相同矩阵就相同。"
            "total 是 factor 加求解。墙钟是单次，不是平均。"
            "冷边 45 C，底面绝热，缝环氧 0.8，中介层 100 um。python，不是 HotSpot，也不是 ATSim。"
        ),
        "old_factor_s": t_old,
        "rows": rows,
    }
    (OUT / "stale_timing.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"写到 {OUT / 'stale_timing.json'}", flush=True)


def many_chips() -> None:
    """芯粒从 4 颗加到格子放满。每一步只换一颗的面积。

    三列都含「这一步的分解」。面积变了矩阵就变，分解省不掉。
    只重算一颗：分解 1 次，回代 1 次。
    每颗都重解：分解 1 次，回代 n 次。
    一次求解：分解 1 次，回代 1 次。
    """
    import time

    sizes = (4, 8, 16, 32, 64, 128)
    rows = []
    for n in sizes:
        old, new, idx = _chip_pair(n)
        k_old = _paint_chips(old)
        fac_old, g_old, _ = factor_edge(k_old)
        kept = [apply_edge(fac_old, g_old, _one(old, i)) for i in range(n)]
        del fac_old

        k_new = _paint_chips(new)
        t0 = time.perf_counter()
        fac, g, = _factor_only(k_new)
        t_factor = time.perf_counter() - t0

        def backsub(power):
            return apply_edge(fac, g, power)

        t0 = time.perf_counter()
        changed = backsub(_one(new, idx))
        t_one = time.perf_counter() - t0
        parts = list(kept)
        parts[idx] = changed
        stale = sum(p - T_AMB for p in parts) + T_AMB

        t0 = time.perf_counter()
        fresh = [backsub(_one(new, i)) for i in range(n)]
        t_each = time.perf_counter() - t0
        frozen = sum(p - T_AMB for p in fresh) + T_AMB

        full = np.zeros((N, N))
        for i in range(n):
            full += _one(new, i)
        t0 = time.perf_counter()
        direct = backsub(full)
        t_direct = time.perf_counter() - t0

        rec = {
            "n_chips": n,
            "changed_index": idx,
            "old_cells": len(hw.cells(old[idx])) if False else int((_one(old, idx) > 0).sum()),
            "new_cells": int((_one(new, idx) > 0).sum()),
            "factor_s": t_factor,
            "one_backsub_s": t_one,
            "each_backsub_s": t_each,
            "direct_backsub_s": t_direct,
            "stale_total_s": t_factor + t_one,
            "each_total_s": t_factor + t_each,
            "direct_total_s": t_factor + t_direct,
            "stale_max_K": float(np.max(np.abs(stale - direct))),
            "each_max_K": float(np.max(np.abs(frozen - direct))),
            "python_Tmax_C": float(direct.max() - 273.15),
        }
        rows.append(rec)
        print(f"n {n:4d}  factor {t_factor*1e3:7.0f} ms  "
              f"one {t_one*1e3:6.1f}  each {t_each*1e3:8.1f}  direct {t_direct*1e3:6.1f} ms  "
              f"stale {rec['stale_max_K']:.3f} K  each {rec['each_max_K']:.1e} K", flush=True)
        del fac

    payload = {
        "note": (
            "同一张 64x64 网格，芯粒数从 4 到 128，小芯粒 2x2，间隙一格。"
            "每套只把一颗 2x2 换成 4x4，功耗按格子加倍，其余不动。"
            "factor 是新材料图的装配加分解，三种做法都付一次。"
            "one 是只回代被换的那颗。each 是每颗各回代一次。direct 是全部功耗一次回代。"
            "total 是 factor 加对应的回代。each 对 direct 应当只剩舍入。"
            "冷边 45 C，底面绝热，缝环氧 0.8。python，不是 HotSpot，也不是 ATSim。"
        ),
        "rows": rows,
    }
    (OUT / "many_chips_timing.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"写到 {OUT / 'many_chips_timing.json'}", flush=True)


def _chip_pair(n: int):
    """放 n 颗 2x2，间距 5 格，所以一颗长成 4x4 仍留一格缝。"""
    pitch = 5
    cols = N // pitch
    chips = []
    for k in range(n):
        col, row = k % cols, k // cols
        j0, ib0 = 1 + col * pitch, 1 + row * pitch
        if j0 + 4 > N or ib0 + 4 > N:
            raise SystemExit(f"{n} 颗放不下")
        chips.append({
            "name": f"c{k}", "w": 2, "h": 2, "j0": j0, "ib0": ib0, "power_W": 2.0,
        })
    new = [dict(c) for c in chips]
    grown = dict(new[0])
    grown["w"], grown["h"] = 4, 4
    grown["power_W"] = 8.0
    new[0] = grown
    return chips, new, 0


def _paint_chips(chips) -> np.ndarray:
    k = np.full((N, N), 0.8)
    for d in chips:
        for r, c in hw.cells(d):
            k[r, c] = 130.0
    return k


def _one(chips, idx) -> np.ndarray:
    p = np.zeros((N, N))
    d = chips[idx]
    cc = hw.cells(d)
    share = d["power_W"] / len(cc)
    for r, c in cc:
        p[r, c] = share
    return p


def _factor_only(k_chip):
    fac, g, _ = factor_edge(k_chip)
    return fac, g


def stale_reuse() -> None:
    """矩阵已经按当前材料图分解过。下一次只换功耗、不换材料图时，三种做法还要付多少。

    这是退火里挪功耗但不改面积的那一步。面积变了的那一步，矩阵必须重分解，见 stale_timing。
    """
    import time

    dies, _cat = fc.base_dies()
    fac, g, t_factor = factor_edge(hw.paint(dies, 0.8))
    powers = [np.zeros((N, N)) for _ in range(4)]
    full = hw.power_map(dies)
    for i in range(4):
        for r, c in hw.cells(dies[i]):
            powers[i][r, c] = full[r, c]

    def bench(fn, n=5):
        fn()
        t0 = time.perf_counter()
        for _ in range(n):
            fn()
        return (time.perf_counter() - t0) / n

    t_one = bench(lambda: apply_edge(fac, g, powers[0]))
    t_four = bench(lambda: [apply_edge(fac, g, p) for p in powers])
    direct = apply_edge(fac, g, full)
    parts = [apply_edge(fac, g, p) for p in powers]
    summed = sum(p - T_AMB for p in parts) + T_AMB
    rec = {
        "factor_once_s": t_factor,
        "one_solve_s": t_one,
        "four_solves_s": t_four,
        "four_vs_direct_K": float(np.max(np.abs(summed - direct))),
        "note": (
            "材料图不变时矩阵不用重分解。只换一颗是 1 次回代，四颗分开是 4 次，"
            "一次求解是 1 次。四颗分开再相加和一次求解差在舍入。"
            "回代时间是 5 次平均。分解是这一张材料图付一次。"
        ),
    }
    print(f"分解 {t_factor*1e3:.1f} ms  回代一次 {t_one*1e3:.2f} ms  "
          f"回代四次 {t_four*1e3:.2f} ms  相加对一次求解 {rec['four_vs_direct_K']:.1e} K", flush=True)
    prev = json.loads((OUT / "stale_timing.json").read_text(encoding="utf-8"))
    prev["reuse"] = rec
    (OUT / "stale_timing.json").write_text(
        json.dumps(prev, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def solve_with_interposer(k_chip: np.ndarray, power: np.ndarray, k_tim: float = 4.0,
                           g_bottom: float = 0.0, edge_to_ambient: bool = False) -> np.ndarray:
    """五层：中介层、芯片、TIM、均热片、热沉。返回 (5, N, N)，单位 K。

    g_bottom 是中介层每一格底面到环境的电导。0 就是原来的绝热底。
    edge_to_ambient 把中介层最外一圈接到环境，用来试「热从板边导出去」。
    """
    layers = (
        ("ip", T_IP, K_IP),
        ("chip", 150e-6, None),
        ("tim", 20e-6, k_tim),
        ("sp", 1.0e-3, 400.0),
        ("sink", 6.9e-3, 400.0),
    )
    nl = len(layers)

    def idx(layer: int, i: int, j: int) -> int:
        return (layer * N + i) * N + j

    nunk = nl * N * N
    acc = [[] for _ in range(nunk)]

    def add(a: int, b: int, g: float) -> None:
        acc[a].append((b, -g))
        acc[b].append((a, -g))
        acc[a].append((a, g))
        acc[b].append((b, g))

    # 中介层整板都是硅，横向。
    g_ip = 1.0 / hs._rx(K_IP, T_IP)
    for i in range(N):
        for j in range(N):
            if j + 1 < N:
                add(idx(0, i, j), idx(0, i, j + 1), g_ip)
            if i + 1 < N:
                add(idx(0, i, j), idx(0, i + 1, j), g_ip)
    # 芯片层横向，缝是填充。
    t_chip = layers[1][1]
    for i in range(N):
        for j in range(N):
            if j + 1 < N:
                r = 0.5 * hs._rx(k_chip[i, j], t_chip) + 0.5 * hs._rx(k_chip[i, j + 1], t_chip)
                add(idx(1, i, j), idx(1, i, j + 1), 1.0 / r)
            if i + 1 < N:
                r = 0.5 * hs._rx(k_chip[i, j], t_chip) + 0.5 * hs._rx(k_chip[i + 1, j], t_chip)
                add(idx(1, i, j), idx(1, i + 1, j), 1.0 / r)
    for layer in (2, 3, 4):
        t, k = layers[layer][1], layers[layer][2]
        g = 1.0 / hs._rx(k, t)
        for i in range(N):
            for j in range(N):
                if j + 1 < N:
                    add(idx(layer, i, j), idx(layer, i, j + 1), g)
                if i + 1 < N:
                    add(idx(layer, i, j), idx(layer, i + 1, j), g)
    # 层间：用下一层的整层厚度，和原来的四层模型同一规矩。
    g_up = 1.0 / hs._rz(K_IP, T_IP)
    for i in range(N):
        for j in range(N):
            add(idx(0, i, j), idx(1, i, j), g_up)
            add(idx(1, i, j), idx(2, i, j), 1.0 / hs._rz(float(k_chip[i, j]), t_chip))
    for layer, (t, k) in ((2, (20e-6, k_tim)), (3, (1e-3, 400.0))):
        g = 1.0 / hs._rz(k, t)
        for i in range(N):
            for j in range(N):
                add(idx(layer, i, j), idx(layer + 1, i, j), g)
    area = hs.CW * hs.CW
    r_sink = hs._rz(400.0, 6.9e-3) + hs.R_CONVEC * (hs.S_SINK * hs.S_SINK) / area
    g_amb = 1.0 / r_sink
    g_edge = 1.0e4
    for i in range(N):
        for j in range(N):
            acc[idx(4, i, j)].append((idx(4, i, j), g_amb))
            g_out = g_bottom
            if edge_to_ambient and (i == 0 or j == 0 or i == N - 1 or j == N - 1):
                g_out += g_edge
            if g_out:
                acc[idx(0, i, j)].append((idx(0, i, j), g_out))

    rows, cols, data = [], [], []
    for a, terms in enumerate(acc):
        merged: dict[int, float] = {}
        for b, v in terms:
            merged[b] = merged.get(b, 0.0) + v
        for b, v in merged.items():
            rows.append(a)
            cols.append(b)
            data.append(v)
    mat = csr_matrix((data, (rows, cols)), shape=(nunk, nunk))
    rhs = np.zeros(nunk)
    for i in range(N):
        for j in range(N):
            rhs[idx(4, i, j)] = g_amb * T_AMB
            rhs[idx(1, i, j)] = power[i, j]
            g_out = g_bottom
            if edge_to_ambient and (i == 0 or j == 0 or i == N - 1 or j == N - 1):
                g_out += g_edge
            if g_out:
                rhs[idx(0, i, j)] += g_out * T_AMB
    temp = splu(mat.tocsc()).solve(rhs)
    return temp.reshape(nl, N, N)


def factor_edge(k_chip: np.ndarray):
    """冷边算子分解一次。换功耗只换右端，不再重装矩阵。"""
    import time

    t0 = time.perf_counter()
    mat, g_terms = _edge_system(k_chip)
    fac = splu(mat.tocsc())
    return fac, g_terms, time.perf_counter() - t0


def _edge_system(k_chip: np.ndarray):
    """和 solve_with_interposer 同一套五层、冷边、底面绝热，只拆出矩阵。"""
    k_tim = 4.0
    layers = (
        ("ip", T_IP, K_IP),
        ("chip", 150e-6, None),
        ("tim", 20e-6, k_tim),
        ("sp", 1.0e-3, 400.0),
        ("sink", 6.9e-3, 400.0),
    )
    nl = len(layers)

    def idx(layer, i, j):
        return (layer * N + i) * N + j

    nunk = nl * N * N
    acc = [[] for _ in range(nunk)]

    def add(a, b, g):
        acc[a].append((b, -g))
        acc[b].append((a, -g))
        acc[a].append((a, g))
        acc[b].append((b, g))

    g_ip = 1.0 / hs._rx(K_IP, T_IP)
    for i in range(N):
        for j in range(N):
            if j + 1 < N:
                add(idx(0, i, j), idx(0, i, j + 1), g_ip)
            if i + 1 < N:
                add(idx(0, i, j), idx(0, i + 1, j), g_ip)
    t_chip = layers[1][1]
    for i in range(N):
        for j in range(N):
            if j + 1 < N:
                r = 0.5 * hs._rx(k_chip[i, j], t_chip) + 0.5 * hs._rx(k_chip[i, j + 1], t_chip)
                add(idx(1, i, j), idx(1, i, j + 1), 1.0 / r)
            if i + 1 < N:
                r = 0.5 * hs._rx(k_chip[i, j], t_chip) + 0.5 * hs._rx(k_chip[i + 1, j], t_chip)
                add(idx(1, i, j), idx(1, i + 1, j), 1.0 / r)
    for layer in (2, 3, 4):
        t, k = layers[layer][1], layers[layer][2]
        g = 1.0 / hs._rx(k, t)
        for i in range(N):
            for j in range(N):
                if j + 1 < N:
                    add(idx(layer, i, j), idx(layer, i, j + 1), g)
                if i + 1 < N:
                    add(idx(layer, i, j), idx(layer, i + 1, j), g)
    g_up = 1.0 / hs._rz(K_IP, T_IP)
    for i in range(N):
        for j in range(N):
            add(idx(0, i, j), idx(1, i, j), g_up)
            add(idx(1, i, j), idx(2, i, j), 1.0 / hs._rz(float(k_chip[i, j]), t_chip))
    for layer, (t, k) in ((2, (20e-6, k_tim)), (3, (1e-3, 400.0))):
        g = 1.0 / hs._rz(k, t)
        for i in range(N):
            for j in range(N):
                add(idx(layer, i, j), idx(layer + 1, i, j), g)
    area = hs.CW * hs.CW
    r_sink = hs._rz(400.0, 6.9e-3) + hs.R_CONVEC * (hs.S_SINK * hs.S_SINK) / area
    g_amb = 1.0 / r_sink
    g_edge = 1.0e4
    edge = []
    for i in range(N):
        for j in range(N):
            acc[idx(4, i, j)].append((idx(4, i, j), g_amb))
            if i == 0 or j == 0 or i == N - 1 or j == N - 1:
                acc[idx(0, i, j)].append((idx(0, i, j), g_edge))
                edge.append(idx(0, i, j))
    rows, cols, data = [], [], []
    for a, terms in enumerate(acc):
        merged = {}
        for b, v in terms:
            merged[b] = merged.get(b, 0.0) + v
        for b, v in merged.items():
            rows.append(a)
            cols.append(b)
            data.append(v)
    mat = csr_matrix((data, (rows, cols)), shape=(nunk, nunk))
    return mat, (g_amb, g_edge, edge)


def apply_edge(fac, g_terms, power: np.ndarray) -> np.ndarray:
    g_amb, g_edge, edge = g_terms
    rhs = np.zeros(5 * N * N)
    rhs[4 * N * N:(5 * N * N)] = g_amb * T_AMB
    rhs[N * N:2 * N * N] = power.ravel()
    for p in edge:
        rhs[p] += g_edge * T_AMB
    return fac.solve(rhs).reshape(5, N, N)[1]


def summarize(T: np.ndarray, dies) -> dict:
    return {
        "python_Tmax_C": float(T.max() - 273.15),
        "dies_C": {d["name"]: die_max(T, d) for d in dies},
    }


def thermal() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dies, _ = fc.base_dies()
    power = hw.power_map(dies)
    report = {
        "layout": [fc.sm.dump_die(d) for d in dies],
        "note": (
            "python 没有中介层时就是原来的四层。加上的是整块 72 mm、100 µm、k=130 的硅，"
            "底面绝热。填充 0.8 是原来的环氧；0.5 是为了和 ATSim 的 organic 对齐。"
            "这些绝对温度不是 HotSpot，也不是 ATSim。"
        ),
        "cases": [],
    }
    fields = {}
    for k_fill, tag in ((0.8, "epoxy"), (0.5, "organic")):
        k_chip = hw.paint(dies, k_fill)
        no = fc.truth(dies, k_fill, 4.0, 400.0)[0]
        yes = solve_with_interposer(k_chip, power)[1]
        fields[tag] = (no, yes)
        np.save(OUT / f"{tag}_no.npy", no)
        np.save(OUT / f"{tag}_with.npy", yes)
        a, b = summarize(no, dies), summarize(yes, dies)
        dT = yes - no
        report["cases"].append({
            "fill": tag,
            "k_fill": k_fill,
            "without": a,
            "with": b,
            "dTmax_K": b["python_Tmax_C"] - a["python_Tmax_C"],
            "dies_dT_K": {n: b["dies_C"][n] - a["dies_C"][n] for n in a["dies_C"]},
            "field_dT_min_K": float(dT.min()),
            "field_dT_max_K": float(dT.max()),
        })
        print(f"{tag}  without {a['python_Tmax_C']:.3f}  with {b['python_Tmax_C']:.3f}  "
              f"dTmax {b['python_Tmax_C']-a['python_Tmax_C']:+.3f}  "
              f"field {dT.min():+.3f}..{dT.max():+.3f}", flush=True)
    (OUT / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_case(rects, path: Path, t_ip: float) -> None:
    """跟代理用的同一套粗网格，但不导入弹性模块。"""
    payload = {
        "name": path.stem,
        "description": "interposer on/off, thermal only",
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
            for c in rects
        ],
        "package": {
            "tim_um": 20, "t_spreader_mm": 1.0, "t_sink_mm": 6.9,
            "t_substrate_mm": 1.0, "t_interposer_um": t_ip,
            "r_convec_K_per_W": 0.1, "ambient_C": 45,
            "pitch_die_mm": 2.0, "pitch_coarse_mm": 8.0,
            "nz_die": 3, "nz_spreader": 3, "nz_sink": 4, "nz_tim": 1,
            "nz_substrate": 3, "nz_interposer": 2,
            "fill_die": True, "fill_xy": "chips", "fill_material": "organic",
        },
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def atsim() -> None:
    sys.path.insert(0, str(ROOT / "热应力python工程_最终"))
    from atsim.fvm3d import solve_fvm3d
    from atsim.geometry import load_case

    report = json.loads((OUT / "report.json").read_text(encoding="utf-8"))
    dies = json.loads(json.dumps(report["layout"]))
    rects = fc.sm.to_rects(dies)
    rows = []
    slices = {}
    for tag, t_ip in (("with", 100.0), ("no", 0.0)):
        jp = OUT / f"atsim_{tag}.json"
        write_case(rects, jp, t_ip)
        case = load_case(str(jp))
        th = solve_fvm3d(case, verbose=False)
        T, grid = th["T"], th["grid"]
        solid = grid["solid"]
        xc, yc, zc = grid["xc"], grid["yc"], grid["zc"]
        ox, oy = case.chip_ox, case.chip_oy
        dies_C = {}
        X = xc[:, None]
        Y = yc[None, :]
        for blk in case.blocks:
            mask_xy = blk.contains_xy(X - ox, Y - oy)
            z0 = 0.0 if blk.z0 is None else blk.z0
            z1 = blk.z1 if blk.z1 is not None else z0
            mz = (zc >= z0 - 1e-16) & (zc < z1 - 1e-16)
            if not np.any(mz):
                mz = (zc >= z0 - 1e-16) & (zc <= z1 + 1e-16)
            vals = T[:, :, mz][mask_xy]
            if vals.size == 0:
                raise RuntimeError(f"no cells in {blk.name}")
            dies_C[blk.name] = float(vals.max() - 273.15)
        # 芯片高度上的最高温度，用来画平面。
        z_die = [blk.z0 for blk in case.blocks if blk.z0 is not None]
        z1 = max(blk.z1 for blk in case.blocks)
        z0 = min(z_die) if z_die else 0.0
        mz = (zc >= z0) & (zc < z1)
        plane = np.full((len(xc), len(yc)), np.nan)
        sub = T[:, :, mz]
        msk = solid[:, :, mz].any(axis=2)
        plane[msk] = sub.max(axis=2)[msk] - 273.15
        slices[tag] = plane
        np.savez(OUT / f"atsim_{tag}_plane.npz", plane=plane, xc=xc, yc=yc)
        rec = {
            "tag": tag,
            "t_interposer_um": t_ip,
            "atsim_Tmax_C": float(T[solid].max() - 273.15),
            "dies_C": dies_C,
            "elapsed_s": th["elapsed_s"],
            "n": [int(grid["nx"]), int(grid["ny"]), int(grid["nz"])],
        }
        rows.append(rec)
        print(f"atsim {tag}  Tmax {rec['atsim_Tmax_C']:.3f}  {dies_C}  {rec['elapsed_s']:.1f}s", flush=True)
    a = next(r for r in rows if r["tag"] == "no")
    b = next(r for r in rows if r["tag"] == "with")
    report["atsim"] = {
        "fill_material": "organic",
        "note": "两边都有 1 mm 有机基板。对流热阻 0.1，不是 python 的 0.039。",
        "rows": rows,
        "dTmax_K": b["atsim_Tmax_C"] - a["atsim_Tmax_C"],
        "dies_dT_K": {n: b["dies_C"][n] - a["dies_C"][n] for n in a["dies_C"]},
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def plot() -> None:
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.patches import Rectangle

    font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    plt.rcParams["font.family"] = "Noto Sans CJK JP"
    plt.rcParams["axes.unicode_minus"] = False
    report = json.loads((OUT / "report.json").read_text(encoding="utf-8"))
    no = np.load(OUT / "epoxy_no.npy")
    yes = np.load(OUT / "epoxy_with.npy")
    dT = yes - no
    extent = (0, 72, 0, 72)
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.2), constrained_layout=True)
    lo = min(no.min(), yes.min()) - 273.15
    hi = max(no.max(), yes.max()) - 273.15
    for ax, field, title in (
        (axes[0], no, "没有中介层"),
        (axes[1], yes, "加 100 µm 硅中介层"),
    ):
        im = ax.imshow(field[::-1] - 273.15, origin="lower", extent=extent, vmin=lo, vmax=hi, cmap="inferno")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label("°C")
    lim = max(abs(dT.min()), abs(dT.max()))
    im = axes[2].imshow(dT[::-1], origin="lower", extent=extent, vmin=-lim, vmax=lim, cmap="coolwarm")
    axes[2].set_title("加上之后的温差")
    fig.colorbar(im, ax=axes[2], fraction=0.046, pad=0.04).set_label("K")
    cw = hs.CW * 1e3
    for ax in axes:
        for d in report["layout"]:
            ax.add_patch(Rectangle((d["j0"] * cw, d["ib0"] * cw), d["w"] * cw, d["h"] * cw,
                                   fill=False, ec="white", lw=0.6))
        ax.set_xlabel("mm")
        ax.set_ylabel("mm")
    fig.suptitle("填充是环氧 0.8。python 场，不是 HotSpot，也不是 ATSim。")
    fig_dir = OUT / "figures"
    fig_dir.mkdir(exist_ok=True)
    fig.savefig(fig_dir / "python_interposer.png", dpi=140)
    plt.close(fig)

    if (OUT / "atsim_no_plane.npz").exists():
        a = np.load(OUT / "atsim_no_plane.npz")
        b = np.load(OUT / "atsim_with_plane.npz")
        fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.2), constrained_layout=True)
        pa, pb = a["plane"], b["plane"]
        lo = np.nanmin([np.nanmin(pa), np.nanmin(pb)])
        hi = np.nanmax([np.nanmax(pa), np.nanmax(pb)])
        ext = (a["xc"][0] * 1e3, a["xc"][-1] * 1e3, a["yc"][0] * 1e3, a["yc"][-1] * 1e3)
        for ax, plane, title in ((axes[0], pa, "ATSim 无中介层"), (axes[1], pb, "ATSim 有中介层")):
            im = ax.imshow(plane.T, origin="lower", extent=ext, vmin=lo, vmax=hi, cmap="inferno")
            ax.set_title(title)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label("°C")
        # 网格可能因为中介层只差在 z 上，xy 应当一致。
        if pa.shape == pb.shape:
            dd = pb - pa
            lim = np.nanmax(np.abs(dd))
            im = axes[2].imshow(dd.T, origin="lower", extent=ext, vmin=-lim, vmax=lim, cmap="coolwarm")
            fig.colorbar(im, ax=axes[2], fraction=0.046, pad=0.04).set_label("K")
        axes[2].set_title("加上之后的温差")
        for ax in axes:
            ax.set_xlabel("mm")
            ax.set_ylabel("mm")
        fig.suptitle("ATSim，填充 organic。坐标是封装自己的，和左图的 72 mm 网格不是同一把尺。")
        fig.savefig(fig_dir / "atsim_interposer.png", dpi=140)
        plt.close(fig)
    print(f"图在 {fig_dir}", flush=True)


if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "thermal"
    if phase == "thermal":
        thermal()
    elif phase == "atsim":
        atsim()
    elif phase == "plot":
        plot()
    elif phase == "stale_area":
        stale_area()
    elif phase == "stale_timing":
        stale_timing()
    elif phase == "stale_reuse":
        stale_reuse()
    else:
        raise SystemExit(phase)
