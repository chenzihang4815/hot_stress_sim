"""同一批候选上，三种做法对照，按芯粒数和填充展开。

5、7、9、11 芯粒，各跑空气 0.028、有机 0.5、环氧 0.8、底部填充 1.0、硅 130。
硅填充整张材料图是同一个数，芯片挪到哪里矩阵都不变，用来当零误差的对照。

三种做法估同一批候选，接受由全部分解决定，所以差的是估计准不准，
不是两条退火走运气走到了别处。
  full    每步按最新布局重画材料图、重分解、重解。
  iter    只在起点分解一次，接受时也不重分解。每步装配新矩阵后，
          用起点的分解做预条件，共轭梯度固定修正 8 次。
  frozen  只在起点分解一次，之后每步把新的功耗直接回代，材料图不再更新。

起点用 sa_heat_vs_coupled 的 lbfgs 布局，和材料场加速那一组是同一个封装。
11 芯粒直接用 11 颗的起点。9 芯粒从 10 芯粒的起点里去掉 CIM_basic，
和加速实验同一口径。间隙改成 0，否则生产间隙会把芯片推出 72 mm。

python 场，冷边 45 C，底面绝热，100 um 硅中介层。不是 HotSpot，也不是 ATSim。
种子 22，200 步，温度从 2.0 起每步乘 0.97。

已有结果的案例会跳过。中断后重跑只补缺的。
"""

import json
import os
import random
import sys
import time
from copy import deepcopy
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

sys.path.insert(0, "/tmp/rankfix")
sys.path.insert(0, "/home/czh/任务仿真产生功耗文件工具/热仿真流水线/src")
import run_interposer_map as im  # noqa: E402
import run_hw_area_field as hw  # noqa: E402
import sa_same  # noqa: E402

N = im.N
ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
RUNS = ROOT / "热仿真流水线" / "runs"
OUT = ROOT / "论文证据" / "填充场叠加" / "复现" / "same_path"
OUT.mkdir(parents=True, exist_ok=True)

# (芯粒数, 布局文件, 去掉的芯片名)
CASES = (
    (5, "sa_heat_vs_coupled_pkg_2cpu_1gpu_2cim/layout_start.json", ()),
    (7, "sa_heat_vs_coupled_pkg_3cpu_1gpu_3cim/layout_start.json", ()),
    (9, "sa_heat_vs_coupled_pkg_4cpu_2gpu_4cim/layout_start.json", ("CIM_basic",)),
    (11, "sa_heat_vs_coupled_pkg_4cpu_3gpu_4cim/layout_start.json", ()),
)
FILLS = (("air", 0.028), ("organic", 0.5), ("epoxy", 0.8),
         ("underfill", 1.0), ("silicon", 130.0))

STEPS = sa_same.STEPS
SEED = sa_same.SEED
N_ITER = sa_same.N_ITER


def load(rel, drop):
    import sa_thermal_feedback as sa
    raw = json.loads((RUNS / rel).read_text())
    chips = [sa.Chip(c["name"], c["kind"], c["x_m"], c["y_m"], c["w_m"], c["h_m"],
                     c["power_W"], c["area_mm2"]) for c in raw if c["name"] not in drop]
    intp = im.hs.SIDE
    sa.GAP_M = sa.GAP_GPU_M = 0.0
    chips = sa.spread_to_gap(chips, intp, intp)
    for c in chips:
        c.x = min(max(0.0, c.x), max(0.0, intp - c.w))
        c.y = min(max(0.0, c.y), max(0.0, intp - c.h))
    return chips, intp


def peak_of(x):
    return float(np.asarray(x[N * N:2 * N * N]).max() - 273.15)


def one(n_chip, rel, drop, fill_name, k_fill):
    import sa_thermal_feedback as sa
    chips0, intp = load(rel, drop)
    dies0 = sa_same.to_dies(chips0)
    mat0, g = im._edge_system(hw.paint(dies0, k_fill))
    fac_start = splu(mat0.tocsc())

    rng = random.Random(SEED)
    cur = deepcopy(chips0)
    t_cur = peak_of(fac_start.solve(sa_same.rhs_of(hw.power_map(dies0), g)))
    temp = sa_same.T0
    rows = []
    wall = {"full": 0.0, "iter": 0.0, "frozen": 0.0}
    n_legal = 0
    for step in range(1, STEPS + 1):
        trial, _ = sa.propose(cur, intp, intp, rng)
        tries = 1
        while not sa.legal(trial, intp, intp) and tries < 25:
            trial, _ = sa.propose(cur, intp, intp, rng)
            tries += 1
        if not sa.legal(trial, intp, intp):
            temp *= sa_same.ALPHA
            continue
        n_legal += 1
        dies = sa_same.to_dies(trial)
        mat, _ = im._edge_system(hw.paint(dies, k_fill))
        rhs = sa_same.rhs_of(hw.power_map(dies), g)

        t0 = time.perf_counter()
        x_full = splu(mat.tocsc()).solve(rhs)
        wall["full"] += time.perf_counter() - t0
        t0 = time.perf_counter()
        x_iter = sa_same.iterate(mat, rhs, fac_start)
        wall["iter"] += time.perf_counter() - t0
        t0 = time.perf_counter()
        x_frozen = fac_start.solve(rhs)
        wall["frozen"] += time.perf_counter() - t0

        tf = peak_of(x_full)
        ti = peak_of(x_iter)
        tz = peak_of(x_frozen)
        took = tf <= t_cur or rng.random() < np.exp(-(tf - t_cur) / temp)
        rows.append({"step": step, "full_C": tf, "iter_C": ti, "frozen_C": tz,
                     "accept": bool(took)})
        if took:
            cur, t_cur = trial, tf
        temp *= sa_same.ALPHA
        if step % 50 == 0:
            print(f"    step {step:3d}  全分解 {tf:.2f}  修正 {ti:.2f}  不更新 {tz:.2f}",
                  flush=True)
    return {"n_legal": n_legal, "wall_s": wall, "rows": rows,
            "n_chips": len(chips0), "power_W": sum(c.power_W for c in chips0)}


def summary_of(rows):
    full = np.array([r["full_C"] for r in rows])
    it = np.array([r["iter_C"] for r in rows])
    fr = np.array([r["frozen_C"] for r in rows])
    return {
        "legal_steps": len(rows),
        "best_full_C": float(full.min()),
        "iter_peak_err_max_K": float(np.max(np.abs(it - full))),
        "iter_peak_err_mean_K": float(np.mean(it - full)),
        "iter_peak_err_mean_abs_K": float(np.mean(np.abs(it - full))),
        "frozen_peak_err_max_K": float(np.max(np.abs(fr - full))),
        "frozen_peak_err_mean_K": float(np.mean(fr - full)),
        "frozen_under_count": int(np.sum(fr < full)),
    }


def main():
    only = {s for s in os.environ.get("SUITE_ONLY", "").split(",") if s}
    for n_chip, rel, drop in CASES:
        for fill_name, k_fill in FILLS:
            tag = f"n{n_chip}_{fill_name}"
            if only and tag not in only:
                continue
            dest = OUT / f"{tag}.json"
            if dest.exists():
                print(f"skip {tag}", flush=True)
                continue
            print(f"== {tag}  k={k_fill} ==", flush=True)
            t0 = time.perf_counter()
            try:
                res = one(n_chip, rel, drop, fill_name, k_fill)
            except Exception as exc:
                print(f"FAILED {tag}  {exc}", flush=True)
                import traceback
                traceback.print_exc()
                continue
            out = {
                "case": tag, "n_chips": res["n_chips"], "fill": fill_name,
                "k_fill": k_fill, "power_W": res["power_W"],
                "steps": STEPS, "seed": SEED, "n_iter": N_ITER,
                "layout": rel, "dropped": list(drop),
                "note": "三条方法估同一批候选，接受由全部分解决定。"
                        "修正 8 次只用起点的分解，接受时不重分解。"
                        "不更新是起点材料图的分解直接回代。"
                        "python 场，冷边 45 C，底面绝热，不是 HotSpot，也不是 ATSim。",
                "wall_s": res["wall_s"],
                "summary": summary_of(res["rows"]),
                "rows": res["rows"],
            }
            dest.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n",
                            encoding="utf-8")
            s = out["summary"]
            print(f"  {tag} 合法 {s['legal_steps']} 步  "
                  f"修正最大 {s['iter_peak_err_max_K']:.3f} K  "
                  f"不更新最大 {s['frozen_peak_err_max_K']:.2f} K  "
                  f"用时 全分解 {res['wall_s']['full']:.0f} s  "
                  f"修正 {res['wall_s']['iter']:.0f} s  "
                  f"共 {time.perf_counter() - t0:.0f} s", flush=True)
            del res
            import gc
            gc.collect()
    print("SAME_PATH_DONE", flush=True)


if __name__ == "__main__":
    main()
