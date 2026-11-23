"""10 步退火。材料图停在起点，对每一步都按最新材料图重算。

2026-11-14 的原脚本，从 /tmp/rankfix/ten.py 收回来。输出改到
论文证据/填充场叠加/复现/，不覆盖 热仿真流水线/runs/interposer_map/ 的存档。

两条走同一批布局：接受由精确解决定。面积不变，位置变，材料图跟着变。
每一步存误差场，最后叠出每格的最大绝对误差。缝环氧。python，不是 HotSpot。
"""
import json
import random
import sys
import time
from copy import deepcopy
import pathlib
from pathlib import Path

import numpy as np

sys.path.insert(0, "/home/czh/任务仿真产生功耗文件工具/热仿真流水线/src")
import run_interposer_map as im  # noqa: E402
import run_hw_area_field as hw  # noqa: E402
import run_factor_controls as fc  # noqa: E402
import sa_thermal_feedback as sa  # noqa: E402

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
OUT = pathlib.Path(__file__).resolve().parent
STEPS = 10
SEED = 22
T0 = 2.0
ALPHA = 0.95


def to_chips(dies):
    cw = im.hs.CW
    return [sa.Chip(d["name"], d.get("role", ""), d["j0"] * cw, d["ib0"] * cw,
                    d["w"] * cw, d["h"] * cw, d["power_W"], d["area_mm2"]) for d in dies]


def to_dies(chips):
    cw = im.hs.CW
    return [{"name": c.name, "j0": int(round(c.x / cw)), "ib0": int(round(c.y / cw)),
             "w": max(1, int(round(c.w / cw))), "h": max(1, int(round(c.h / cw))),
             "power_W": c.power_W, "area_mm2": c.area_mm2} for c in chips]


def main():
    dies, _ = fc.base_dies()
    intp = im.hs.SIDE
    chips = to_chips(dies)
    if not sa.legal(chips, intp, intp):
        raise SystemExit("起点不合法")
    fac, g, t_fac = im.factor_edge(hw.paint(to_dies(chips), 0.8))

    rng = random.Random(SEED)
    cur = chips
    x0 = fac.solve(_rhs(hw.power_map(to_dies(cur)), g))
    t_cur = float(x0[im.N * im.N:2 * im.N * im.N].max() - 273.15)
    temp = T0
    rows, err_max = [], np.zeros((im.N, im.N))
    grids = {"step0_exact": x0[im.N * im.N:2 * im.N * im.N].reshape(im.N, im.N)}
    wall = {"exact": 0.0, "frozen": 0.0}
    layouts = []
    for step in range(1, STEPS + 1):
        trial, op = sa.propose(cur, intp, intp, rng)
        tries = 1
        while not sa.legal(trial, intp, intp) and tries < 40:
            trial, op = sa.propose(cur, intp, intp, rng)
            tries += 1
        if not sa.legal(trial, intp, intp):
            temp *= ALPHA
            print(f"step {step} 无合法候选", flush=True)
            continue
        d = to_dies(trial)
        k = hw.paint(d, 0.8)
        p = hw.power_map(d)
        t0 = time.perf_counter()
        mat, g_new = im._edge_system(k)
        exact = im.splu(mat.tocsc()).solve(_rhs(p, g_new))
        wall["exact"] += time.perf_counter() - t0
        t0 = time.perf_counter()
        frozen = fac.solve(_rhs(p, g))
        wall["frozen"] += time.perf_counter() - t0
        ef = exact[im.N * im.N:2 * im.N * im.N].reshape(im.N, im.N)
        ff = frozen[im.N * im.N:2 * im.N * im.N].reshape(im.N, im.N)
        err = ff - ef
        err_max = np.maximum(err_max, np.abs(err))
        tf = float(ef.max() - 273.15)
        tz = float(ff.max() - 273.15)
        took = tf <= t_cur or rng.random() < np.exp(-(tf - t_cur) / temp)
        rows.append({
            "step": step, "op": op, "accept": bool(took),
            "exact_Tmax_C": tf, "frozen_Tmax_C": tz,
            "peak_err_K": abs(tz - tf),
            "grid_max_K": float(np.max(np.abs(err))),
            "grid_mean_K": float(np.mean(np.abs(err))),
        })
        grids[f"step{step}_exact"] = ef
        grids[f"step{step}_frozen"] = ff
        grids[f"step{step}_err"] = err
        layouts.append([{"name": c.name, "j0": int(round(c.x / im.hs.CW)),
                         "ib0": int(round(c.y / im.hs.CW)),
                         "w": int(round(c.w / im.hs.CW)), "h": int(round(c.h / im.hs.CW))}
                        for c in trial])
        print(f"step {step:2d} {op:32} {'取' if took else '弃'}  "
              f"精确 {tf:6.2f}  旧图 {tz:6.2f}  整图 {rows[-1]['grid_max_K']:7.2f} K", flush=True)
        if took:
            cur, t_cur = trial, tf
        temp *= ALPHA

    grids["err_max_over_steps"] = err_max
    np.savez(OUT / "ten_frozen_grids.npz", **grids)
    (OUT / "ten_frozen_layouts.json").write_text(
        json.dumps(layouts, ensure_ascii=False), encoding="utf-8")
    (OUT / "ten_frozen.json").write_text(json.dumps({
        "note": (
            "10 步，种子 22，T0 2.0，每步乘 0.95。提议是生产退火的平移、跳动、对调、旋转，"
            "面积不变，位置变。接受由精确解决定，两条走同一批布局。"
            "frozen 用起点材料图的分解回代，不重分解。exact 每步按最新材料图重分解重解。"
            "err_max_over_steps 是 10 张误差场逐格取绝对值的最大。"
            "冷边 45 C，底面绝热，缝环氧 0.8。python，不是 HotSpot，也不是 ATSim。"
        ),
        "start_factor_s": t_fac,
        "wall_s": wall,
        "rows": rows,
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"精确 {wall['exact']:.2f} s   旧图 {wall['frozen']:.3f} s", flush=True)


def _rhs(power, g):
    g_amb, g_edge, edge = g
    rhs = np.zeros(5 * im.N * im.N)
    rhs[4 * im.N * im.N:5 * im.N * im.N] = g_amb * im.T_AMB
    rhs[im.N * im.N:2 * im.N * im.N] = np.asarray(power).ravel()
    for q in edge:
        rhs[q] += g_edge * im.T_AMB
    return rhs


if __name__ == "__main__":
    main()
