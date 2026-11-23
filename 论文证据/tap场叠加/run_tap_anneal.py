#!/usr/bin/env python3
"""在 TAP-2.5D 自己的退火上对比全解和场叠加。

退火、移动、初始布局都用 TAP 的代码：sim_annealing.py 的 close_neighbor、
jumping_neighbor，system.initial_placement 的 B*-tree。不改 TAP 目录里的任何文件。

温度只替换一个入口。TAP 每步在 passive_interposer.run_hotspot 里调一次
HotSpot，这里换成 run_hotspot_timed：全解是一次求解，场叠加是把每个芯粒单独
的温度场缓存起来再相加，没移动的芯粒不再求解。接受只按峰值温度，不按线长，
所以不调用 routing.solve_Cplex，这台机器没有 CPLEX。

网格行数由 TAP_ROWS 控制，写进 TAP 自己生成的 new_hotspot.config，不改
util/hotspot.config。赢家都回到 64 网格重评，热图用的是这一次重评的温度场。

环境变量：TAP_CASE、TAP_ROWS、TAP_STEPS、TAP_MODE（full 或 field）、TAP_OUT、TAP_SEED。
TAP_INIT 给出一份布局 json 时，跳过 B*-tree，从这份布局出发。
"""
import json
import os
import random
import subprocess
import sys
import time
from copy import deepcopy
from pathlib import Path

import numpy as np

# routing.py 第一行就 import cplex，config.py 又在顶层 import routing。这台机器
# 没有 CPLEX，而这一步只按温度接受，不调用 solve_Cplex。放一个空模块挡住导入，
# TAP 目录里的文件不动。
sys.modules.setdefault("cplex", type(sys)("cplex"))

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
TAP = ROOT / "TAP-2.5D"
sys.path.insert(0, str(TAP))

os.chdir(TAP)

# fastSA.py 第 73 行把浮点传给 random.randint，那是 Python 2 的写法，Python 3.14
# 直接报错。只在这一处取整，TAP 的文件不动。
_randint = random.randint
random.randint = lambda a, b, _r=_randint: _r(int(a), int(b))

import sim_annealing  # noqa: E402
import config  # noqa: E402

CASE = os.environ["TAP_CASE"]
ROWS = int(os.environ.get("TAP_ROWS", "64"))
STEPS = int(os.environ.get("TAP_STEPS", "120"))
MODE = os.environ["TAP_MODE"]
OUT = Path(os.environ["TAP_OUT"])
SEED = int(os.environ.get("TAP_SEED", "22"))
ROWS_FINAL = 64
AMBIENT = 318.15
JUMP_RATIO = 0.9
T_START = 1.0
DECAY = 0.8

HS = TAP / "util" / "hotspot"


def load_system():
    """TAP 的 read_config 从 sys.argv 读参数，这里临时替换。"""
    saved = sys.argv
    sys.argv = ["run_tap_anneal.py", "-c", f"configs/sys_{CASE}.cfg",
                "-d", str(OUT / "tap_work") + "/"]
    try:
        system = config.read_config()
    finally:
        sys.argv = saved
    return system


def set_rows(system, rows):
    """TAP 第一次 gen_flp 时把 util/hotspot.config 复制成 new_hotspot.config，
    之后不再覆盖。网格行数在这里改，原文件不动。"""
    cfg = Path(system.path) / "new_hotspot.config"
    lines = []
    for line in cfg.read_text().splitlines():
        if "-grid_rows" in line:
            line = line[:line.index("-grid_rows")] + f"-grid_rows\t\t\t{rows}"
        elif "-grid_cols" in line:
            line = line[:line.index("-grid_cols")] + f"-grid_cols\t\t\t{rows}"
        lines.append(line)
    cfg.write_text("\n".join(lines) + "\n")


def layer_grid(grid_file, rows):
    """TAP 自带的 util/hotspot 把网格写成 序号\\t温度，没有 Layer 分段。
    我们工程里的 HotSpot-master 才有 Layer 分段，两者不能混用。"""
    vals = []
    for line in Path(grid_file).read_text().splitlines():
        sp = line.split()
        if len(sp) == 2:
            vals.append(float(sp[1]))
    if len(vals) != rows * rows:
        raise RuntimeError(f"网格有 {len(vals)} 格，期望 {rows * rows}")
    return np.array(vals)


def solve_once(system, filename, powers, rows):
    """一次 HotSpot 求解。powers 为 None 时用 TAP 自己的 gen_ptrace。"""
    system.gen_flp(filename)
    if powers is None:
        system.gen_ptrace(filename)
    else:
        flp = Path(system.path) / f"{filename}L4_ChipLayer.flp"
        names = [sp[0] for line in flp.read_text().splitlines()
                 if (sp := line.split()) and not sp[0].startswith("#")]
        chip = {f"Chiplet_{i}": powers[i] for i in range(system.chiplet_count)}
        ptrace = Path(system.path) / f"{filename}.ptrace"
        ptrace.write_text("\t".join(names) + "\n"
                          + "\t".join(f"{chip.get(n, 0):.6g}" for n in names) + "\n")
    t0 = time.time()
    proc = subprocess.run(
        [str(HS), "-c", system.path + "new_hotspot.config",
         "-f", system.path + filename + "L4_ChipLayer.flp",
         "-p", system.path + filename + ".ptrace",
         "-steady_file", system.path + filename + ".steady",
         "-grid_steady_file", system.path + filename + ".grid.steady",
         "-model_type", "grid", "-detailed_3D", "on",
         "-grid_layer_file", system.path + filename + "layers.lcf"],
        capture_output=True)
    dt = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.decode()[-300:])
    grid = layer_grid(system.path + filename + ".grid.steady", rows)
    # TAP 的 clean_hotspot 用花括号通配，这个 shell 不展开，中间文件删不掉。
    # 自己删，不改 TAP 的函数。
    for f in Path(system.path).glob(filename + "*"):
        if f.is_file():
            f.unlink()
    return grid, dt


class Timer:
    """替换 TAP 的 run_hotspot。返回峰值 °C，把求解秒数和次数记在自己身上。"""

    def __init__(self, rows, use_field):
        self.rows = rows
        self.use_field = use_field
        self.cache = {}
        self.secs = 0.0
        self.solves = 0

    def __call__(self, system, filename):
        if not self.use_field:
            grid, dt = solve_once(system, filename, None, self.rows)
            self.secs += dt
            self.solves += 1
            return float(grid.max()) - 273.15
        acc = None
        n = system.chiplet_count
        for i in range(n):
            key = (i, system.x[i], system.y[i], system.width[i], system.height[i])
            if key not in self.cache:
                powers = [0.0] * n
                powers[i] = system.power[i]
                self.cache[key], dt = solve_once(system, filename, powers, self.rows)
                self.secs += dt
                self.solves += 1
            acc = self.cache[key] if acc is None else acc + self.cache[key]
        return float((acc - (n - 1) * AMBIENT).max()) - 273.15


def accept(old, new, T):
    """TAP 注释掉的只按温度的接受准则，sim_annealing.py 第 106 行。"""
    delta = old - new
    return 1.0 if delta > 0 else np.exp(delta / T)


def save_layout(system, path):
    chips = [{"i": i, "x_mm": system.x[i], "y_mm": system.y[i],
              "w_mm": system.width[i], "h_mm": system.height[i],
              "power_W": system.power[i]} for i in range(system.chiplet_count)]
    path.write_text(json.dumps(chips, ensure_ascii=False, indent=2), encoding="utf-8")


def anneal(system, timer):
    random.seed(SEED)
    np.random.seed(SEED)
    step = 1
    temp = timer(system, f"step_{step}")
    best = deepcopy(system)
    temp_best = temp
    grid = sim_annealing.block_occupation.initialize_grid(
        int(system.intp_size / system.granularity))
    for i in range(system.chiplet_count):
        grid = sim_annealing.block_occupation.set_block_occupation(
            grid, system.granularity, system.x[i], system.y[i],
            system.width[i] + 2 * system.hubump[i],
            system.height[i] + 2 * system.hubump[i], i)
    T = T_START
    trace = []
    while step < STEPS:
        step += 1
        if 1 - JUMP_RATIO > random.random():
            moved, xn, yn, rot = sim_annealing.jumping_neighbor(system, grid)
        else:
            moved, xn, yn = sim_annealing.close_neighbor(system, grid)
            rot = 0
        trial = deepcopy(system)
        trial.x[moved], trial.y[moved] = xn, yn
        if rot == 1:
            trial.rotate(moved)
        t_new = timer(trial, f"step_{step}")
        took = bool(random.random() < accept(temp, t_new, T))
        trace.append({"step": step, "moved": int(moved), "rot": int(rot),
                      "from": [system.x[moved], system.y[moved]],
                      "to": [xn, yn], "temp": round(t_new, 2), "accepted": took,
                      "solve_secs": round(timer.secs, 1), "solves": timer.solves})
        # 每步落盘。64 网格一步几十秒，中途断了也能从这里接着看。
        (OUT / "trace.json").write_text(
            json.dumps(trace, ensure_ascii=False), encoding="utf-8")
        if step % 10 == 0:
            print(f"  第 {step} 步 温度 {t_new:.2f} 累计求解 {timer.secs:.0f}s "
                  f"{timer.solves} 次", flush=True)
        if took:
            grid = sim_annealing.block_occupation.clear_block_occupation(
                grid, system.granularity, system.x[moved], system.y[moved],
                system.width[moved] + 2 * system.hubump[moved],
                system.height[moved] + 2 * system.hubump[moved], moved)
            system = trial
            grid = sim_annealing.block_occupation.set_block_occupation(
                grid, system.granularity, xn, yn,
                system.width[moved] + 2 * system.hubump[moved],
                system.height[moved] + 2 * system.hubump[moved], moved)
            temp = t_new
            if t_new < temp_best:
                best, temp_best = deepcopy(system), t_new
        if step % int(system.intp_size) == 0:
            T *= DECAY
    return best, temp_best, trace


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # B*-tree 起点用了随机数，种子必须在生成初始布局之前设，否则每次跑的起点
    # 不同，全解和场叠加就不是从同一个布局出发。
    random.seed(SEED)
    np.random.seed(SEED)
    system = load_system()
    system.gen_flp("init")
    init_path = os.environ.get("TAP_INIT")
    if init_path:
        given = json.loads(Path(init_path).read_text())
        for c in given:
            i = c["i"]
            system.x[i], system.y[i] = c["x_mm"], c["y_mm"]
            system.width[i], system.height[i] = c["w_mm"], c["h_mm"]
    set_rows(system, ROWS)
    save_layout(system, OUT / "init.json")
    # 初始布局也回到 64 网格重评，热图的左半幅用这一份。
    set_rows(system, ROWS_FINAL)
    init_grid, _ = solve_once(system, "init_eval", None, ROWS_FINAL)
    (OUT / "init_64.grid").write_text(
        "\n".join(f"{v:.6f}" for v in init_grid) + "\n", encoding="utf-8")
    set_rows(system, ROWS)

    timer = Timer(ROWS, MODE == "field")
    t0 = time.time()
    best, search_peak, trace = anneal(system, timer)
    (OUT / "trace.json").write_text(
        json.dumps(trace, ensure_ascii=False, indent=1), encoding="utf-8")
    wall = time.time() - t0
    save_layout(best, OUT / "best.json")

    # 赢家回到 64 网格重评，温度场存盘供热图使用。
    set_rows(best, ROWS_FINAL)
    grid, _ = solve_once(best, "final", None, ROWS_FINAL)
    (OUT / "best_64.grid").write_text(
        "\n".join(f"{v:.6f}" for v in grid) + "\n", encoding="utf-8")
    final = float(grid.max()) - 273.15

    summary = {
        "case": CASE, "rows_search": ROWS, "mode": MODE, "steps": STEPS, "seed": SEED,
        "init": os.environ.get("TAP_INIT"),
        "search_peak_C": round(search_peak, 2),
        "final_64_C": round(final, 2),
        "solve_secs": round(timer.secs, 1),
        "solves": timer.solves,
        "wall_s": round(wall, 1),
    }
    (OUT / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
