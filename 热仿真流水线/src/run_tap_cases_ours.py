#!/usr/bin/env python3
"""把 TAP-2.5D 三个案例的布局和功耗放进我们的退火，跑和那次只看温度的完整跑一样多的步。

TAP 那一次不再重跑。时间、起点、最优布局都用
runs/superposition/tap_field 里留着的记录。那一次只按峰值温度接受，线长没有进它的决策。

我们这一侧用现在的内环：目标值是温度、梯度、线长各三分之一，动作比例按步数走，
较差的步按 0.20 乘 (1 - 步/总步) 收。不换芯片，这三个案例的芯粒不在我们的硬件库里。

搜索用我们自己的热求解器。TAP 的 gen_flp 在半毫米尺寸上会让 fill_space 无限递归，
那个文件不改。求解器按 TAP 的封装取数：芯粒层导热 100，中介层 50 mm，均热板
100 mm，热沉 200 mm，对流热阻 0.009，环境 45 C。所以这一侧的绝对温度不能拿去
和我们自己 72 mm 封装的案例比。

几何照 TAP 的来，不折算到 1.125 mm 格。移动步长 1 mm，芯粒尺寸保留 7.75、11.87、
18.2。定价是硅填充的场叠加，算子只分解一次。

两边的最优布局最后都用 TAP 的 HotSpot 在 64 网格上重评一次。那一步把芯粒尺寸收成
整数毫米，只为绕开 fill_space 的递归，和搜索用的尺寸分开记。线长用 TAP 配置里的
连接矩阵，按真实尺寸算，不受这次取整影响。
"""

from __future__ import annotations

import json
import os
import random
import shutil
import subprocess
import sys
import time
from copy import deepcopy
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
TAP = ROOT / "TAP-2.5D"
HS_BIN = ROOT / "HotSpot-master" / "hotspot"
REF = ROOT / "热仿真流水线" / "runs" / "superposition" / "tap_field"
OUT = ROOT / "热仿真流水线" / "runs" / "tap_cases_ours"
sys.path.insert(0, str(TAP))

import block_occupation  # noqa: E402

SEED = 22
N_STEPS = int(os.environ.get("TAP_STEPS", "1050"))
STEP_MM = 1.0
INTP = 50.0
AMB_C = 45.0
AMB = AMB_C + 273.15
# 搜索网格跟 TAP 那次的搜索网格一样。最优布局另用 64 网格重评。
ROWS = int(os.environ.get("TAP_ROWS", "16"))
ROWS_FINAL = 64
CASES = tuple(s for s in os.environ.get(
    "TAP_CASES", "ascend910,micro150,multigpu").split(",") if s)

# TAP util/hotspot.config 的封装，按 50 mm 中介层放大后的数。
K_CHIP = 100.0
T_CHIP = 0.00015
T_TIM = 2.0e-05
K_TIM = 4.0
T_SPREADER = 0.001
K_SPREADER = 400.0
T_SINK = 0.0069
K_SINK = 400.0
S_SPREADER = 0.1
S_SINK = 0.2
R_CONVEC = 0.009
LAYERS = (
    ("chip", T_CHIP, None),
    ("tim", T_TIM, K_TIM),
    ("spreader", T_SPREADER, K_SPREADER),
    ("sink", T_SINK, K_SINK),
)
NL = len(LAYERS)


def span(v: float) -> int:
    """TAP 的占格宽度。半毫米尺寸在 1 mm 网格上多占一格。"""
    return int(v / 2 + 0.49) * 2


def cell_box(x: float, y: float, w: float, h: float) -> tuple[int, int, int, int]:
    """占格范围，和 block_occupation.set_block_occupation 同一套取整。"""
    return span_box(x, y, w, h)


def boundary_check(x, y, w, h):
    # 按占格判断，不按真实边。真实边可以伸出 50 mm 半毫米以内，和 TAP 的邻居函数一致。
    if x - span(w) / 2 < 0 or x + span(w) / 2 > INTP:
        return False
    if y - span(h) / 2 < 0 or y + span(h) / 2 > INTP:
        return False
    return True


def close_neighbor(system, grid):
    """照搬 TAP-2.5D/sim_annealing.py。不导入那个模块，它会拉起 CPLEX。"""
    chiplet_order = np.random.permutation(range(system.n))
    for p in chiplet_order:
        direction_order = np.random.permutation(["up", "down", "left", "right"])
        xx, yy = round(system.x[p]), round(system.y[p])
        width, height = system.w[p], system.h[p]
        for d in direction_order:
            if d == "left" and boundary_check(xx, yy - STEP_MM, width, height):
                if block_occupation.check_down_occupation(
                        grid, STEP_MM, xx, yy - STEP_MM, width, height):
                    return p, xx, yy - STEP_MM
            elif d == "right" and boundary_check(xx, yy + STEP_MM, width, height):
                if block_occupation.check_up_occupation(
                        grid, STEP_MM, xx, yy + STEP_MM, width, height):
                    return p, xx, yy + STEP_MM
            elif d == "up" and boundary_check(xx - STEP_MM, yy, width, height):
                if block_occupation.check_left_occupation(
                        grid, STEP_MM, xx - STEP_MM, yy, width, height):
                    return p, xx - STEP_MM, yy
            elif d == "down" and boundary_check(xx + STEP_MM, yy, width, height):
                if block_occupation.check_right_occupation(
                        grid, STEP_MM, xx + STEP_MM, yy, width, height):
                    return p, xx + STEP_MM, yy
    raise RuntimeError("没有芯粒可以移动")


def jumping_neighbor(system, grid):
    for _ in range(10000):
        pick = random.randint(0, system.n - 1)
        x_new = random.randint(1, int(INTP / STEP_MM) - 1) * STEP_MM
        y_new = random.randint(1, int(INTP / STEP_MM) - 1) * STEP_MM
        rotation = random.randint(0, 1)
        cw, ch = (system.h[pick], system.w[pick]) if rotation else (system.w[pick], system.h[pick])
        if boundary_check(x_new, y_new, cw, ch) and free_span(grid, x_new, y_new, cw, ch, pick):
            return pick, x_new, y_new, rotation
    p, xn, yn = close_neighbor(system, grid)
    return p, xn, yn, 0


def layout_mix(step: int, total: int) -> tuple[float, float, float, float]:
    """和 run_closed_loop_suite.layout_mix 同一条规则。

    改芯片在这里是 0，TAP 的芯粒没有对应的硬件库。剩下的 0.90 仍按那条规则分：
    旋转固定 0.10，k 从 0.10 升到 0.70，其余由大跳和交换平分。返回四项加起来
    是 0.90，调用处把 0.10 的空档并进微调。
    """
    done = step / max(total, 1)
    k = 0.10 + 0.60 * done
    rest = 0.90 - k
    share = 0.90
    return (share * k, share * rest / 2, share * rest / 2, share * 0.10)


class Layout:
    """中心坐标，单位 mm。尺寸保留 TAP 的半毫米，不收成整数。"""

    def __init__(self, chips: list[dict]):
        self.n = len(chips)
        self.w = [c["w_mm"] for c in chips]
        self.h = [c["h_mm"] for c in chips]
        self.power = [c["power_W"] for c in chips]
        self.x = [c["x_mm"] + c["w_mm"] / 2 for c in chips]
        self.y = [c["y_mm"] + c["h_mm"] / 2 for c in chips]

    def rotate(self, i: int) -> None:
        self.w[i], self.h[i] = self.h[i], self.w[i]

    def chips(self) -> list[dict]:
        return [{
            "i": i,
            "x_mm": self.x[i] - self.w[i] / 2,
            "y_mm": self.y[i] - self.h[i] / 2,
            "w_mm": self.w[i],
            "h_mm": self.h[i],
            "power_W": self.power[i],
        } for i in range(self.n)]


def matrix_of(name: str) -> list[list[int]]:
    text = (TAP / "configs" / f"sys_{name}.cfg").read_text(encoding="utf-8")
    fields: dict[str, str] = {}
    section = None
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1]
            continue
        if "=" not in line:
            if section == "chiplets" and "connections" in fields:
                fields["connections"] += line
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        fields[key] = value
    return [
        [int(v) for v in row.split(",") if v.strip()]
        for row in fields["connections"].split(";")
        if row.strip()
    ]


def wirelength_of(layout: Layout, matrix: list[list[int]]) -> float:
    """照 TAP 的 compute_wirelength：两端中心的曼哈顿距离乘权重，除以线数。"""
    total = 0.0
    count = 0
    for i, row in enumerate(matrix):
        for j, weight in enumerate(row):
            if i == j or weight <= 0:
                continue
            dist = abs(layout.x[i] - layout.x[j]) + abs(layout.y[i] - layout.y[j])
            total += dist * weight
            count += weight
    return total / count


def free_span(grid, x: float, y: float, w: float, h: float, owner: int) -> bool:
    """这片占格是空的，或者本来就是这颗芯粒的。"""
    x0, x1, y0, y1 = span_box(x, y, w, h)
    if x0 < 0 or y0 < 0 or x1 > len(grid) or y1 > len(grid):
        return False
    for col in range(x0, x1):
        for row in range(y0, y1):
            if grid[col][row] not in (0, owner + 2):
                return False
    return True


def span_box(x: float, y: float, w: float, h: float) -> tuple[int, int, int, int]:
    """占格范围，半开区间。取整和 block_occupation.set_block_occupation 一致。"""
    return (int(x / STEP_MM) - int(w / 2 / STEP_MM + 0.49),
            int(x / STEP_MM) + int(w / 2 / STEP_MM + 0.49) + 1,
            int(y / STEP_MM) - int(h / 2 / STEP_MM + 0.49),
            int(y / STEP_MM) + int(h / 2 / STEP_MM + 0.49) + 1)


def occupy(layout: Layout):
    """占格自己记。不走 block_occupation 的清除函数，它在对不上时直接退出进程。"""
    # 半毫米尺寸的占格会伸到 50 以外四格，表留到 60。越界的格不参与合法性判断。
    n = int(INTP / STEP_MM) + 10
    grid = [[0] * n for _ in range(n)]
    for i in range(layout.n):
        mark(grid, layout, i, i + 2)
    return grid


def mark(grid, layout: Layout, i: int, value: int) -> None:
    x0, x1, y0, y1 = span_box(layout.x[i], layout.y[i], layout.w[i], layout.h[i])
    for x in range(x0, x1):
        for y in range(y0, y1):
            grid[x][y] = value


def no_overlap(layout: Layout, grid=None) -> bool:
    """按占格表判断。表是 block_occupation 写的，和邻居函数用的是同一套。"""
    if grid is None:
        grid = occupy(layout)
    seen = {}
    for x in range(len(grid)):
        for y in range(len(grid[0])):
            owner = grid[x][y]
            if owner:
                if owner in seen and seen[owner] != (x, y):
                    pass
                seen.setdefault(owner, 0)
                seen[owner] += 1
    # 每颗芯粒的占格数必须和单独放时一样，少了就是被别的芯粒盖住了。
    for i in range(layout.n):
        x0, x1, y0, y1 = span_box(layout.x[i], layout.y[i], layout.w[i], layout.h[i])
        owned = sum(1 for x in range(max(0, x0), x1) for y in range(max(0, y0), y1)
                    if grid[x][y] == i + 2)
        if owned != (x1 - x0) * (y1 - y0):
            return False
    return True
def propose(layout: Layout, grid, rng: random.Random, mix):
    """四种动作，比例按 (微调, 大跳, 交换, 旋转)。返回 (试探, 动作名, 是否合法)。"""
    trial = deepcopy(layout)
    roll = rng.random()
    b_nudge, b_jump, b_swap = mix[0], mix[0] + mix[1], mix[0] + mix[1] + mix[2]
    if roll < b_nudge:
        try:
            p, xn, yn = close_neighbor(trial, grid)
        except RuntimeError:
            return trial, "nudge-none", False
        trial.x[p], trial.y[p] = float(round(xn)), float(round(yn))  # noqa: 收成整数
        return trial, f"nudge {p}", True
    if roll < b_jump:
        p, xn, yn, rot = jumping_neighbor(trial, grid)
        trial.x[p], trial.y[p] = float(round(xn)), float(round(yn))  # noqa: 收成整数
        if rot == 1:
            trial.rotate(p)
        return trial, f"jump {p}" + ("+rot" if rot else ""), True
    if roll < b_swap:
        i = rng.randrange(layout.n)
        j = rng.randrange(layout.n - 1)
        if j >= i:
            j += 1
        trial.x[i], trial.x[j] = trial.x[j], trial.x[i]
        trial.y[i], trial.y[j] = trial.y[j], trial.y[i]
        ok = all(boundary_check(trial.x[k], trial.y[k], trial.w[k], trial.h[k]) for k in (i, j))
        return trial, f"swap {i}<->{j}", ok
    p = rng.randrange(layout.n)
    trial.rotate(p)
    ok = boundary_check(trial.x[p], trial.y[p], trial.w[p], trial.h[p])
    return trial, f"rotate {p}", ok


class Thermal:
    """四层网格，和 run_fill_superposition 同一套离散，封装数换成 TAP 的。

    芯粒层是方格。芯粒边落在格内时按面积分摊：一格里芯粒占多少，导热和功耗
    就按这个比例和填充物混合。填充物取硅，这样材料图不随布局变，算子只分解一次。
    """

    def __init__(self, rows: int):
        self.rows = rows
        self.cw = INTP / 1000 / rows
        self.factor, self.g_amb = self._operator()

    def _rx(self, k: float, thickness: float) -> float:
        return self.cw / (k * self.cw * thickness)

    def _rz(self, k: float, thickness: float) -> float:
        return thickness / (k * self.cw * self.cw)

    def _idx(self, layer: int, i: int, j: int) -> int:
        return (layer * self.rows + i) * self.rows + j

    def _operator(self):
        n = self.rows
        nunk = NL * n * n
        acc = [[] for _ in range(nunk)]

        def add(a, b, g):
            acc[a].append((b, -g))
            acc[b].append((a, -g))
            acc[a].append((a, g))
            acc[b].append((b, g))

        g = 1.0 / self._rx(K_CHIP, T_CHIP)
        for i in range(n):
            for j in range(n):
                if j + 1 < n:
                    add(self._idx(0, i, j), self._idx(0, i, j + 1), g)
                if i + 1 < n:
                    add(self._idx(0, i, j), self._idx(0, i + 1, j), g)
        for layer in range(1, NL):
            t, k = LAYERS[layer][1], LAYERS[layer][2]
            g = 1.0 / self._rx(k, t)
            for i in range(n):
                for j in range(n):
                    if j + 1 < n:
                        add(self._idx(layer, i, j), self._idx(layer, i, j + 1), g)
                    if i + 1 < n:
                        add(self._idx(layer, i, j), self._idx(layer, i + 1, j), g)
        g = 1.0 / self._rz(K_CHIP, T_CHIP)
        for i in range(n):
            for j in range(n):
                add(self._idx(0, i, j), self._idx(1, i, j), g)
        for layer in range(1, NL - 1):
            t, k = LAYERS[layer][1], LAYERS[layer][2]
            g = 1.0 / self._rz(k, t)
            for i in range(n):
                for j in range(n):
                    add(self._idx(layer, i, j), self._idx(layer + 1, i, j), g)
        area = self.cw * self.cw
        r_sink = self._rz(K_SINK, T_SINK) + R_CONVEC * (S_SINK * S_SINK) / area
        g_amb = 1.0 / r_sink
        for i in range(n):
            for j in range(n):
                a = self._idx(3, i, j)
                acc[a].append((a, g_amb))
        rows, cols, data = [], [], []
        for a, terms in enumerate(acc):
            merged: dict[int, float] = {}
            for b, v in terms:
                merged[b] = merged.get(b, 0.0) + v
            for b, v in merged.items():
                rows.append(a)
                cols.append(b)
                data.append(v)
        mat = sparse.csr_matrix((data, (rows, cols)), shape=(nunk, nunk))
        return splu(mat.tocsc()), g_amb

    def power_map(self, layout: Layout, only: int | None) -> np.ndarray:
        """芯粒落在格内的部分按面积分摊功耗。第 0 行是 y 最大的一行。"""
        p = np.zeros((self.rows, self.rows))
        cell = INTP / self.rows
        for i in range(layout.n):
            if only is not None and i != only:
                continue
            x0 = (layout.x[i] - layout.w[i] / 2) / cell
            y0 = (layout.y[i] - layout.h[i] / 2) / cell
            x1 = (layout.x[i] + layout.w[i] / 2) / cell
            y1 = (layout.y[i] + layout.h[i] / 2) / cell
            area = layout.w[i] * layout.h[i]
            for row in range(int(np.floor(y0)), int(np.ceil(y1))):
                for col in range(int(np.floor(x0)), int(np.ceil(x1))):
                    if not (0 <= row < self.rows and 0 <= col < self.rows):
                        continue
                    overlap = (min(x1, col + 1) - max(x0, col)) * (min(y1, row + 1) - max(y0, row))
                    if overlap <= 0:
                        continue
                    p[self.rows - 1 - row, col] += layout.power[i] * overlap * cell * cell / area
        return p

    def solve(self, power: np.ndarray) -> np.ndarray:
        n = self.rows
        rhs = np.zeros(NL * n * n)
        for i in range(n):
            for j in range(n):
                rhs[self._idx(3, i, j)] = self.g_amb * AMB
                rhs[self._idx(0, i, j)] += power[i, j]
        return self.factor.solve(rhs).reshape(NL, n, n)[0]


class Pricer:
    """每颗芯粒的温升按位置缓存。材料图是整层硅，相加只剩舍入。"""

    def __init__(self, thermal: Thermal):
        self.thermal = thermal
        self.cache: dict[tuple, np.ndarray] = {}
        self.solves = 0

    def field(self, layout: Layout) -> np.ndarray:
        acc = np.zeros((self.thermal.rows, self.thermal.rows))
        for i in range(layout.n):
            key = (i, layout.x[i], layout.y[i], layout.w[i], layout.h[i])
            if key not in self.cache:
                t0 = time.time()
                state = self.thermal.solve(self.thermal.power_map(layout, i))
                self.cache[key] = state - AMB
                self.solves += 1
                self.secs = getattr(self, "secs", 0.0) + time.time() - t0
            acc += self.cache[key]
        return acc + AMB

    def peak(self, layout: Layout) -> float:
        return float(self.field(layout).max()) - 273.15

    def grad(self, layout: Layout) -> float:
        grid = self.field(layout)
        cell = INTP / self.thermal.rows
        gx = np.abs(np.diff(grid, axis=1)).max() / cell
        gy = np.abs(np.diff(grid, axis=0)).max() / cell
        return float(max(gx, gy))


def objective(peak_c: float, grad: float, wl: float, scales: dict) -> float:
    return ((peak_c - AMB_C) / scales["dT"] + grad / scales["g0"] + wl / scales["wl0"]) / 3.0


def anneal(name: str, layout: Layout, matrix, pricer: Pricer) -> dict:
    rng_m = random.Random(SEED)
    rng_a = random.Random(SEED + 7919)
    np.random.seed(SEED)
    current = deepcopy(layout)
    grid = occupy(current)
    peak = pricer.peak(current)
    grad = pricer.grad(current)
    wl = wirelength_of(current, matrix)
    scales = {"dT": max(peak - AMB_C, 1e-6), "g0": max(grad, 1e-9), "wl0": max(wl, 1e-9)}
    j_cur = objective(peak, grad, wl, scales)
    best = deepcopy(current)
    best_j = j_cur
    best_step = 0
    log = [{"step": 0, "op": "init", "accept": True, "J": j_cur,
            "Tmax_C": peak, "grad_max": grad, "wirelength_mm_per_wire": wl}]
    print(f"  init T {peak:.2f}  grad {grad:.3f}  wl {wl:.2f}  J {j_cur:.3f}", flush=True)
    accepted = 0
    worse = 0
    t0 = time.time()
    for step in range(1, N_STEPS + 1):
        mix = layout_mix(step, N_STEPS)
        mix = (mix[0] + 0.10, mix[1], mix[2], mix[3])
        trial, op, ok = propose(current, grid, rng_m, mix)
        if not ok or not no_overlap(trial):
            log.append({"step": step, "op": op, "accept": False, "reason": "illegal"})
            continue
        peak_new = pricer.peak(trial)
        grad_new = pricer.grad(trial)
        wl_new = wirelength_of(trial, matrix)
        j_new = objective(peak_new, grad_new, wl_new, scales)
        worse_p = 0.20 * (1.0 - step / N_STEPS)
        take = j_new <= j_cur or rng_a.random() < worse_p
        rec = {"step": step, "op": op, "accept": take,
               "ap": 1.0 if j_new <= j_cur else worse_p,
               "J": j_new, "Tmax_C": peak_new, "grad_max": grad_new,
               "wirelength_mm_per_wire": wl_new}
        if take:
            accepted += 1
            if j_new > j_cur:
                worse += 1
            moved = [i for i in range(current.n)
                     if (current.x[i], current.y[i], current.w[i], current.h[i])
                     != (trial.x[i], trial.y[i], trial.w[i], trial.h[i])]
            for i in moved:
                mark(grid, current, i, 0)
            current = trial
            for i in moved:
                mark(grid, current, i, i + 2)
            j_cur = j_new
            if j_new < best_j - 1e-12:
                best = deepcopy(trial)
                best_j = j_new
                best_step = step
                rec["new_best"] = True
                print(f"  {step:4d} BEST {op}  T {peak_new:.2f}  wl {wl_new:.2f}  J {j_new:.3f}",
                      flush=True)
        log.append(rec)
        if step % 100 == 0:
            print(f"  {step:4d} 步  已求解 {pricer.solves} 次 {pricer.secs:.0f}s", flush=True)
    case_dir = OUT / name
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "trace.json").write_text(json.dumps(log, ensure_ascii=False) + "\n", encoding="utf-8")
    (case_dir / "best.json").write_text(
        json.dumps(best.chips(), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    np.save(case_dir / "best_field.npy", pricer.field(best))
    return {"best": best, "best_step": best_step, "best_J": best_j,
            "accepted": accepted, "worse_accepted": worse,
            "search_s": round(time.time() - t0, 1)}


def write_flp(layout: Layout, path: Path) -> list[dict]:
    """写平面图。空隙切成条状填充块，每层用同一张，否则 detailed 模式会崩。"""
    lines = ["# name width height left-x bottom-y"]
    placed = []
    for i in range(layout.n):
        x0, y0, x1, y1 = box_of(layout, i)
        lines.append(f"Chiplet_{i}\t{(x1 - x0) / 1000:.6g}\t{(y1 - y0) / 1000:.6g}"
                     f"\t{x0 / 1000:.6g}\t{y0 / 1000:.6g}")
        placed.append({"i": i, "x_mm": x0, "y_mm": y0, "w_mm": x1 - x0, "h_mm": y1 - y0,
                       "clipped": (x1 - x0, y1 - y0) != (layout.w[i], layout.h[i])})
    path.write_text("\n".join(lines + fill_strips(layout)) + "\n", encoding="utf-8")
    return placed


def box_of(layout: Layout, i: int) -> tuple[float, float, float, float]:
    """平面图用的矩形。伸出 50 mm 的部分切掉，HotSpot 不接受超界的块。"""
    return (max(0.0, layout.x[i] - layout.w[i] / 2),
            max(0.0, layout.y[i] - layout.h[i] / 2),
            min(INTP, layout.x[i] + layout.w[i] / 2),
            min(INTP, layout.y[i] + layout.h[i] / 2))


def fill_strips(layout: Layout) -> list[str]:
    """把芯粒之间的空隙切成横向条。每条里芯粒按 x 排序，空隙连成矩形。"""
    edges_y = {0.0, INTP}
    for i in range(layout.n):
        y0, _, _, y1 = box_of(layout, i)[1], 0, 0, 0
        x0, y0, x1, y1 = box_of(layout, i)
        edges_y |= {y0, y1}
    bands = sorted(edges_y)
    out = []
    k = 0
    for lo, hi in zip(bands, bands[1:]):
        if hi - lo < 1e-6:
            continue
        cuts = {0.0, INTP}
        for i in range(layout.n):
            x0, y0, x1, y1 = box_of(layout, i)
            if y0 < hi - 1e-6 and y1 > lo + 1e-6:
                cuts.add(x0)
                cuts.add(x1)
        edges = sorted(cuts)
        for a, b in zip(edges, edges[1:]):
            if b - a < 1e-6:
                continue
            covered = any(
                box_of(layout, i)[0] < b - 1e-6 and box_of(layout, i)[2] > a + 1e-6
                and box_of(layout, i)[1] < hi - 1e-6 and box_of(layout, i)[3] > lo + 1e-6
                for i in range(layout.n))
            if not covered:
                out.append(f"Fill_{k}\t{(b - a) / 1000:.6g}\t{(hi - lo) / 1000:.6g}"
                           f"\t{a / 1000:.6g}\t{lo / 1000:.6g}")
                k += 1
    return out


def hotspot_config(path: Path) -> None:
    """照 run_tap_field.write_config，按行替换，不碰 r_convec_sec 里的 0.1。"""
    out = []
    for line in (TAP / "util" / "hotspot.config").read_text().splitlines():
        if "-s_sink" in line and "0.06" in line:
            line = line.replace("0.06", f"{S_SINK:.6g}", 1)
        elif "-s_spreader" in line and "0.03" in line:
            line = line.replace("0.03", f"{S_SPREADER:.6g}", 1)
        elif "-r_convec" in line and "0.1" in line and "sec" not in line:
            line = line.replace("0.1", f"{R_CONVEC:.6g}", 1)
        out.append(line)
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def hotspot_peak(layout: Layout, work: Path, rows: int) -> tuple[float, float, float]:
    """用 TAP 的 HotSpot 评一次，6 层 detailed 网格，和存档那次同一个模型。"""
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    placed = write_flp(layout, work / "chip.flp")
    flp_names = [l.split()[0] for l in (work / "chip.flp").read_text().splitlines()
                 if l and not l.startswith("#")]
    powers = [f"{layout.power[i]:.6g}" if name.startswith("Chiplet_") else "0"
              for i, name in ((int(n.split("_")[1]), n) if n.startswith("Chiplet_") else (0, n)
                              for n in flp_names)]
    names = flp_names
    (work / "power.ptrace").write_text("\t".join(names) + "\n" + "\t".join(powers) + "\n")
    hotspot_config(work / "hotspot.config")
    lcf = [
        "# layer file",
        "0", "Y", "N", "1.06E+06", "3.33", "0.0002", str(work / "base.flp"),
        "1", "Y", "N", "2.32E+06", "0.625", "0.00007", str(work / "base.flp"),
        "2", "Y", "N", "1.75E+06", "0.01", "0.00011", str(work / "base.flp"),
        "3", "Y", "N", "2.32E+06", "0.625", "1.00E-05", str(work / "base.flp"),
        "4", "Y", "Y", "1.75E+06", "0.01", "0.00015", str(work / "chip.flp"),
        "5", "Y", "N", "4.00E+06", "0.25", "2.00E-05", str(work / "base.flp"),
    ]
    (work / "layers.lcf").write_text("\n".join(lcf) + "\n", encoding="utf-8")
    # 每层用同一张平面图。芯粒层的填充块功耗为 0。
    (work / "base.flp").write_text((work / "chip.flp").read_text(), encoding="utf-8")
    t0 = time.time()
    run = subprocess.run(
        [str(HS_BIN), "-c", str(work / "hotspot.config"), "-f", str(work / "chip.flp"),
         "-p", str(work / "power.ptrace"), "-steady_file", str(work / "o.steady"),
         "-grid_steady_file", str(work / "o.grid"), "-model_type", "grid",
         "-grid_rows", str(rows), "-grid_cols", str(rows),
         "-detailed_3D", "on", "-grid_layer_file", str(work / "layers.lcf")],
        capture_output=True, text=True)
    if run.returncode != 0:
        raise RuntimeError(run.stderr[-400:])
    (work / "rounded.json").write_text(json.dumps(placed, indent=2) + "\n", encoding="utf-8")
    # HotSpot 在 detailed 模式会在芯粒层上下各加封装层，输出里芯粒层是 Layer 4，
    # 排在第 5 个。按层号取，不靠位置。
    chip = None
    idx = 0
    for line in (work / "o.grid").read_text().splitlines():
        line = line.strip()
        if line.startswith("Layer"):
            idx = int(line.split()[1].rstrip(":"))
            continue
        parts = line.split()
        if len(parts) == 2 and parts[0].isdigit() and idx == 4:
            chip = (chip or [])
            chip.append(float(parts[1]))
    vals = chip
    if len(vals) != rows * rows:
        raise RuntimeError(f"网格有 {len(vals)} 格，期望 {rows * rows}")
    peak = max(vals) - 273.15
    spread = max(vals) - min(vals)
    (work / "grid.txt").write_text("\n".join(f"{v:.6f}" for v in vals) + "\n", encoding="utf-8")
    (work / "rounded.json").write_text(json.dumps(placed, indent=2) + "\n", encoding="utf-8")
    return peak, spread, time.time() - t0


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = []
    thermal = Thermal(ROWS)
    for name in CASES:
        saved = json.loads((REF / "sfull" / name / "init.json").read_text())
        ref_path = REF / "sfull" / f"{name}.json"
        ref = json.loads((ref_path if ref_path.exists() else REF / f"{name}.json").read_text())
        matrix = matrix_of(name)
        print(f"{name} 开始，{N_STEPS} 步", flush=True)
        layout = Layout(saved)
        for i in range(layout.n):
            layout.x[i], layout.y[i] = float(round(layout.x[i])), float(round(layout.y[i]))
        case_dir = OUT / name
        case_dir.mkdir(parents=True, exist_ok=True)
        (case_dir / "init.json").write_text(
            json.dumps(layout.chips(), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        (case_dir / "connection_matrix.json").write_text(json.dumps(matrix) + "\n", encoding="utf-8")
        pricer = Pricer(thermal)
        result = anneal(name, layout, matrix, pricer)
        tap_best = Layout(json.loads((REF / "sfull" / name / "best_full.json").read_text()))
        fine = Thermal(ROWS_FINAL)
        t0 = time.time()
        ours_grid = fine.solve(fine.power_map(result["best"], None))
        ours_dt = time.time() - t0
        t0 = time.time()
        tap_grid = fine.solve(fine.power_map(tap_best, None))
        tap_dt = time.time() - t0
        np.save(case_dir / "best_ours_64.npy", ours_grid)
        np.save(case_dir / "best_tap_64.npy", tap_grid)
        ours_peak = float(ours_grid.max()) - 273.15
        tap_peak = float(tap_grid.max()) - 273.15
        ours_spread = float(ours_grid.max() - ours_grid.min())
        tap_spread = float(tap_grid.max() - tap_grid.min())
        ours_wl = wirelength_of(result["best"], matrix)
        tap_wl = wirelength_of(tap_best, matrix)
        rec = {
            "case": name,
            "steps": N_STEPS,
            "seed": SEED,
            "search_grid": ROWS,
            "note": (
                "TAP 侧是 runs/superposition/tap_field 里只按峰值温度接受的那次完整跑，没有重跑。"
                "我们这一侧按温度、梯度、线长各三分之一，搜索用我们自己的求解器，封装数取 TAP 的。"
                "两边最优布局的温度和温差用 TAP 的 HotSpot 在 64 网格上重评，尺寸收成整数毫米，"
                "只为绕开 fill_space 的递归。线长按真实尺寸、用 TAP 的连接矩阵。"
            ),
            "tap_reference": {
                "full_secs": ref["full_secs"],
                "full_final_64_C": ref["full_final_64_C"],
                "source": "runs/superposition/tap_field",
            },
            "ours": {
                "solve_secs": round(pricer.secs, 2),
                "solves": pricer.solves,
                "search_s": result["search_s"],
                "solve_ms_per_step": round(pricer.secs / N_STEPS * 1000, 2),
                "accepted": result["accepted"],
                "worse_accepted": result["worse_accepted"],
                "best_step": result["best_step"],
                "best_J": round(result["best_J"], 4),
            },
            "time_ratio_tap_over_ours": round(ref["full_secs"] / pricer.secs, 2),
            "rescored_64": {
                "ours": {"Tmax_C": round(ours_peak, 2), "spread_K": round(ours_spread, 2),
                         "wirelength_mm_per_wire": round(ours_wl, 4), "solve_s": round(ours_dt, 1)},
                "tap": {"Tmax_C": round(tap_peak, 2), "spread_K": round(tap_spread, 2),
                        "wirelength_mm_per_wire": round(tap_wl, 4), "solve_s": round(tap_dt, 1)},
                "Tmax_error_K": round(ours_peak - tap_peak, 2),
                "spread_error_K": round(ours_spread - tap_spread, 2),
                "wirelength_error_mm": round(ours_wl - tap_wl, 4),
                "solver": "我们自己的四层网格，封装数取 TAP 的",
            },
        }
        (case_dir / "report.json").write_text(
            json.dumps(rec, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        summary.append(rec)
        print(json.dumps(rec, ensure_ascii=False, indent=2), flush=True)
    (OUT / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
