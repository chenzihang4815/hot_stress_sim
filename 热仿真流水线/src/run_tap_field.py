#!/usr/bin/env python3
"""TAP-2.5D 的退火只留温度，用场叠加替换每步的全解。

TAP 原版每步对 6 层 detailed_3D 模型做一次完整 HotSpot。这里走它自己的邻居
移动（90% 邻格微调，10% 大跳）和它自己的退火日程（T 从 1.0 乘 decay 降到
0.01，每个温度 intp_size 步），但接受概率只看峰值温度，线长权重去掉，因此
不需要 CPLEX。

同一条搜索用两种计价各跑一遍：每步全解，和缓存的单芯粒温度场叠加。场只在
芯粒移动后失效，其余芯粒复用。两边都在 16 网格上跑，赢家最后用 64 网格复核
一次。一次只开一个 HotSpot，工作目录随解随删。
"""

from __future__ import annotations

import math
import os
import random
import shutil
import subprocess
import sys
import time
from copy import deepcopy
from pathlib import Path

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
TAP = ROOT / "TAP-2.5D"
HS = ROOT / "HotSpot-master" / "hotspot"
sys.path.insert(0, str(TAP))

from passive_interposer import PassiveInterposer  # noqa: E402
import block_occupation  # noqa: E402

# 下面三个函数照搬 TAP-2.5D/sim_annealing.py。不直接导入是因为那个模块在
# import 时会拉起 config 和 routing，而 routing 需要这台机器上没有的 CPLEX。

def boundary_check(system, x, y, w, h):
    if (x - w / 2) < 0 or (x + w / 2) > system.intp_size:
        return False
    if (y - h / 2) < 0 or (y + h / 2) > system.intp_size:
        return False
    return True


def close_neighbor(system, grid):
    chiplet_order = np.random.permutation(range(system.chiplet_count))
    granularity = system.granularity
    for p in chiplet_order:
        direction_order = np.random.permutation(['up', 'down', 'left', 'right'])
        xx, yy = system.x[p], system.y[p]
        width = system.width[p] + 2 * system.hubump[p]
        height = system.height[p] + 2 * system.hubump[p]
        for d in direction_order:
            if d == 'left':
                if block_occupation.check_down_occupation(grid, granularity, xx, yy - granularity, width, height):
                    return p, xx, yy - granularity
            elif d == 'right' and yy + granularity <= system.intp_size:
                if block_occupation.check_up_occupation(grid, granularity, xx, yy + granularity, width, height):
                    return p, xx, yy + granularity
            elif d == 'up':
                if block_occupation.check_left_occupation(grid, granularity, xx - granularity, yy, width, height):
                    return p, xx - granularity, yy
            elif d == 'down' and xx + granularity <= system.intp_size:
                if block_occupation.check_right_occupation(grid, granularity, xx + granularity, yy, width, height):
                    return p, xx + granularity, yy
    raise RuntimeError('没有芯粒可以移动')


def jumping_neighbor(system, grid):
    n = system.chiplet_count
    granularity = system.granularity
    for count in range(10000):
        pick = random.randint(0, n - 1)
        x_new = random.randint(1, int(system.intp_size / granularity) - 1) * granularity
        y_new = random.randint(1, int(system.intp_size / granularity) - 1) * granularity
        rotation = random.randint(0, 1)
        if rotation == 1:
            cw, ch = system.height[pick], system.width[pick]
        else:
            ch, cw = system.height[pick], system.width[pick]
        if boundary_check(system, x_new, y_new, cw + 2 * system.hubump[pick], ch + 2 * system.hubump[pick]) and \
                block_occupation.replace_block_occupation(grid, granularity, x_new, y_new, cw + 2 * system.hubump[pick], ch + 2 * system.hubump[pick], pick):
            return pick, x_new, y_new, rotation
    p, xn, yn = close_neighbor(system, grid)
    return p, xn, yn, 0


OUT = Path(os.environ["TAP_OUT"]) if os.environ.get("TAP_OUT") else (
    ROOT / "热仿真流水线" / "runs" / "superposition" / "tap_field")
WORK = OUT / "_work"
AMB = 318.15
SEED = 22
DECAY = 0.8
T_START = 1.0
T_MIN = 0.01
N_STEPS = int(os.environ["TAP_STEPS"]) if os.environ.get("TAP_STEPS") else None
JUMP_RATIO = 0.9
# 占格粒度和移动步长，跟 TAP 原版一致。
STEP_MM = 1
# 16 网格上 256 W 的芯粒低估约 2.5°C，64 网格一步要 42 秒。搜索用 16，
# 赢家用 64 复核，两者的差距单列。
ROWS_SEARCH = int(os.environ.get("TAP_ROWS", "16"))
ROWS_FINAL = 64


def make_system(label, widths, heights, powers, positions, intp_size):
    s = PassiveInterposer()
    s.path = str(OUT / label) + "/"
    s.intp_size = intp_size
    s.granularity = STEP_MM
    s.chiplet_count = len(widths)
    s.width = list(widths)
    s.height = list(heights)
    s.power = list(powers)
    s.hubump = [0.0] * s.chiplet_count
    s.x = [p[0] for p in positions]
    s.y = [p[1] for p in positions]
    Path(s.path).mkdir(parents=True, exist_ok=True)
    return s


def write_config(path, intp_size):
    """热沉和均热板按 TAP 的公式取，网格行数留给命令行覆盖。"""
    src = (TAP / "util" / "hotspot.config").read_text().splitlines()
    size_spreader = 2 * intp_size / 1000
    size_heatsink = 2 * size_spreader
    r_convec = 0.1 * 0.06 * 0.06 / size_heatsink / size_heatsink
    out = []
    for line in src:
        if "-s_sink" in line and "0.06" in line:
            line = line.replace("0.06", f"{size_heatsink:.6g}", 1)
        elif "-s_spreader" in line and "0.03" in line:
            line = line.replace("0.03", f"{size_spreader:.6g}", 1)
        elif "-r_convec" in line and "0.1" in line:
            line = line.replace("0.1", f"{r_convec:.6g}", 1)
        out.append(line)
    Path(path).write_text("\n".join(out) + "\n")


def occupy(system):
    grid = block_occupation.initialize_grid(int(system.intp_size / system.granularity))
    for i in range(system.chiplet_count):
        grid = block_occupation.set_block_occupation(
            grid, system.granularity, system.x[i], system.y[i],
            system.width[i], system.height[i], i)
    return grid


def chip_peak(steady: Path) -> float:
    """芯粒层（layer 3）上各芯粒块的最高温度。"""
    best = None
    for line in steady.read_text().splitlines():
        name, val = line.split()
        if name.startswith("layer_3_Chiplet_"):
            v = float(val)
            best = v if best is None else max(best, v)
    return best - 273.15


def grid_of(system, powers, rows):
    """一次 6 层求解，返回芯粒所在层的网格温度。

    填充块随芯粒位置重切，稳态文件里的块名对不齐，所以场必须按格子存。
    HotSpot 在芯粒层上下还各加一个封装层，芯粒层是第 4 层（下标 3 起的
    第 4 个内层之前有衬底等，实测芯粒层为 Layer 4）。
    """
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True)
    system.path = str(WORK) + "/"
    system.gen_flp("s")
    names = [sp[0] for line in (WORK / "sL4_ChipLayer.flp").read_text().splitlines()
             if (sp := line.split()) and not sp[0].startswith("#")]
    chip = {f"Chiplet_{i}": powers[i] for i in range(system.chiplet_count)}
    (WORK / "s.ptrace").write_text(
        "\t".join(names) + "\n"
        + "\t".join(f"{chip.get(n, 0):.6g}" for n in names) + "\n")
    t0 = time.time()
    r = subprocess.run(
        [str(HS), "-c", str(OUT / "hotspot.config"), "-f", str(WORK / "sL4_ChipLayer.flp"),
         "-p", str(WORK / "s.ptrace"), "-steady_file", str(WORK / "o.steady"),
         "-grid_steady_file", str(WORK / "o.grid"), "-model_type", "grid",
         "-grid_rows", str(rows), "-grid_cols", str(rows),
         "-detailed_3D", "on", "-grid_layer_file", str(WORK / "slayers.lcf")],
        capture_output=True, text=True)
    dt = time.time() - t0
    if r.returncode != 0:
        shutil.rmtree(WORK, ignore_errors=True)
        raise RuntimeError(r.stderr[-200:])
    # 芯粒层在 Layer 4。稳态文件里 layer_4 的芯粒块比 layer_3 高 0.03°C，
    # 是八层里最热的一层，和 TAP 取整层最大值的口径一致。
    vals = []
    take = False
    for line in (WORK / "o.grid").read_text().splitlines():
        line = line.strip()
        if line.startswith("Layer"):
            if take:
                break
            take = line.startswith("Layer 4")
            continue
        if take and line:
            vals.append(float(line.split()[-1]))
    shutil.rmtree(WORK, ignore_errors=True)
    if len(vals) != rows * rows:
        raise RuntimeError(f"Layer 4 有 {len(vals)} 格，期望 {rows*rows}")
    return np.array(vals), dt


def solve(system, powers, rows) -> tuple[float, float]:
    g, dt = grid_of(system, powers, rows)
    return float(g.max()) - 273.15, dt


def save_layout(system, path: Path) -> None:
    import json
    chips = [{"i": i, "x_mm": system.x[i], "y_mm": system.y[i],
              "w_mm": system.width[i], "h_mm": system.height[i],
              "power_W": system.power[i]} for i in range(system.chiplet_count)]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(chips, ensure_ascii=False, indent=2), encoding="utf-8")


def final_grid(system, rows, path: Path):
    """64 网格复核，温度场存盘，供热图使用。"""
    g, _ = grid_of(system, system.power, rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(f"{v:.6f}" for v in g) + "\n", encoding="utf-8")
    return g, float(g.max()) - 273.15


def anneal(system, use_field, rows):
    """TAP 的退火，只按温度接受。返回 (最优系统, 最优°C, 求解秒数, 求解次数, 步数)。"""
    random.seed(SEED)
    np_seed()
    cache = {}
    cur = deepcopy(system)
    grid = occupy(cur)

    def peak_of(sys):
        if not use_field:
            g, dt = grid_of(sys, sys.power, rows)
            return float(g.max()) - 273.15, dt, 1
        acc = None
        fresh = 0
        dt = 0.0
        for i in range(sys.chiplet_count):
            key = (i, sys.x[i], sys.y[i], sys.width[i], sys.height[i])
            if key not in cache:
                powers = [0.0] * sys.chiplet_count
                powers[i] = sys.power[i]
                cache[key], d = grid_of(sys, powers, rows)
                fresh += 1
                dt += d
            acc = cache[key] if acc is None else acc + cache[key]
        n = sys.chiplet_count
        return float((acc - (n - 1) * AMB).max()) - 273.15, dt, fresh

    temp, dt, solves = peak_of(cur)
    t_used = dt
    n_solves = solves
    best = deepcopy(cur)
    temp_best = temp
    T = T_START
    step = 0
    # 每步的累计求解秒数，供倍数随步数的曲线使用。
    curve = []
    while (N_STEPS is None and T >= T_MIN) or (N_STEPS is not None and step < N_STEPS):
        step += 1
        trial = deepcopy(cur)
        if random.random() > 1 - JUMP_RATIO:
            moved, xn, yn = close_neighbor(trial, grid)
            rot = 0
        else:
            moved, xn, yn, rot = jumping_neighbor(trial, grid)
        old_w, old_h = cur.width[moved], cur.height[moved]
        # 收成整数毫米。半毫米的芯粒尺寸（8.25、7.75）在 0.5 mm 边界上会产生
        # 0.0004999999999999996 这样的浮点坐标，TAP 的 fill_space 把它判到切割
        # 区域外面，递归就停不下来。
        trial.x[moved], trial.y[moved] = round(xn), round(yn)
        if rot == 1:
            trial.rotate(moved)
        t_new, dt, solves = peak_of(trial)
        t_used += dt
        n_solves += solves
        curve.append(round(t_used, 2))
        delta = temp - t_new
        ap = 1.0 if delta >= 0 else math.exp(delta / T)
        if random.random() < ap:
            grid = block_occupation.clear_block_occupation(
                grid, cur.granularity, cur.x[moved], cur.y[moved],
                old_w, old_h, moved)
            cur = trial
            grid = block_occupation.set_block_occupation(
                grid, cur.granularity, cur.x[moved], cur.y[moved],
                cur.width[moved], cur.height[moved], moved)
            temp = t_new
            if t_new < temp_best:
                best, temp_best = deepcopy(cur), t_new
        if step % int(system.intp_size) == 0:
            T *= DECAY
    return best, temp_best, t_used, n_solves, step, curve


def np_seed():
    import numpy as np
    np.random.seed(SEED)


def pack_positions(widths, heights, intp):
    """固定的紧凑起点。TAP 的初始布局由 B*-tree 生成，那一步需要 CPLEX，这台
    机器没有。占格按 TAP 的取整方式算：半毫米尺寸（8.25、7.75、11.87）在 1 mm
    网格上多占一格，间距和边界都按这个多出来的格子留空。按面积从大到小放，
    小芯粒填到大芯粒旁边的空隙里。"""
    def span(v):
        return int(v / 2 + 0.49) * 2

    spans = [(span(w), span(h)) for w, h in zip(widths, heights)]
    order = sorted(range(len(spans)), key=lambda i: spans[i][0] * spans[i][1], reverse=True)
    placed = []
    for i in order:
        sw, sh = spans[i]
        found = None
        for y in range(1, intp - sh, 1):
            for x in range(1, intp - sw, 1):
                box = (x, y, x + sw, y + sh)
                if all(box[2] + 2 <= b[0] or b[2] + 2 <= box[0]
                       or box[3] + 2 <= b[1] or b[3] + 2 <= box[1] for b, _ in placed):
                    found = box
                    break
            if found:
                break
        if found is None:
            raise RuntimeError(f"芯粒 {i} 放不进 {intp} mm 中介层")
        placed.append((found, i))
    pos = [None] * len(spans)
    for (x0, y0, x1, y1), i in placed:
        pos[i] = ((x0 + x1) / 2, (y0 + y1) / 2)
    return pos


CASES = [
    # 尺寸和功耗取自 TAP-2.5D/configs/sys_*.cfg，中介层都是 50 mm。
    # ascend910 用配置里的真实矩形，不用早期对照里的方形芯粒。
    ("ascend910", [14.5, 10.5, 7.75, 7.75, 7.75, 7.75],
     [31.4, 16, 11.87, 11.87, 11.87, 11.87],
     [256, 14, 20, 20, 20, 20], None, 50),
    ("micro150", [8.25, 8.25, 8.25, 8.25, 8.75, 8.75, 8.75, 8.75],
     [9, 9, 9, 9, 8.75, 8.75, 8.75, 8.75],
     [150, 150, 150, 150, 20, 20, 20, 20], None, 50),
    ("multigpu", [12, 18.2, 18.2, 7.75, 7.75, 7.75],
     [12, 18.2, 18.2, 11.87, 11.87, 11.87],
     [105, 295, 295, 20, 20, 20], None, 50),
]
for i, (label, w, h, pw, pos, intp) in enumerate(CASES):
    if pos is None:
        CASES[i] = (label, w, h, pw, pack_positions(w, h, intp), intp)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    want = set(os.environ.get("TAP_CASES", "").split(",")) - {""}
    for label, w, h, p, pos, intp in CASES:
        if want and label not in want:
            continue
        write_config(OUT / "hotspot.config", intp)
        system = make_system(label, w, h, p, pos, intp)
        print(f"{label} 开始，{intp} mm 中介层", flush=True)
        t0 = time.time()
        best_full, coarse_full, secs_full, solves_full, steps, curve_full = anneal(system, False, ROWS_SEARCH)
        print(f"  全解 {steps} 步 {secs_full:.0f}s 搜索峰值 {coarse_full:.2f}", flush=True)
        save_layout(system, OUT / label / "init.json")
        save_layout(best_full, OUT / label / "best_full.json")
        g_full, final_full = final_grid(best_full, ROWS_FINAL, OUT / label / "best_full_64.grid")
        t1 = time.time()
        best_field, coarse_field, secs_field, solves_field, _, curve_field = anneal(system, True, ROWS_SEARCH)
        print(f"  场叠加 {secs_field:.0f}s 求解 {solves_field} 次 搜索峰值 {coarse_field:.2f}", flush=True)
        save_layout(best_field, OUT / label / "best_field.json")
        g_field, final_field = final_grid(best_field, ROWS_FINAL, OUT / label / "best_field_64.grid")
        rec = {
            "label": label, "steps": steps, "step_mm": STEP_MM,
            "rows_search": ROWS_SEARCH,
            "full_secs": round(secs_full, 1), "full_final_64_C": round(final_full, 2),
            "field_secs": round(secs_field, 1), "field_solves": solves_field,
            "field_final_64_C": round(final_field, 2),
            "speedup": round(secs_full / secs_field, 2),
            "delta_C": round(final_field - final_full, 2),
            "wall_s": round(time.time() - t0, 1),
            "curve_full_s": curve_full,
            "curve_field_s": curve_field,
        }
        (OUT / f"{label}.json").write_text(json_dumps(rec), encoding="utf-8")
        print(json_dumps(rec), flush=True)


def json_dumps(rec) -> str:
    import json
    return json.dumps(rec, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
