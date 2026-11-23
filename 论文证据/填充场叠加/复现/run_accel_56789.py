"""5 到 9 芯粒，材料场加速对每步全部分解。

加速只在起点分解一次，每步装配新矩阵后用起点的分解固定修正 8 次，
接受时不重分解，也不按误差停。精确每步都重新分解。
两边种子相同，但接受由自己的峰值决定，所以布局会分开。
误差是加速选出的布局再用全部分解重算，重算的时间单独记，不算进退火。

5、6、7 芯粒用 sa_heat_vs_coupled 的 lbfgs 起点。
8、9 芯粒没有现成的封装，从 10、11 芯粒的起点里去掉功耗最小的 CIM。
间隙改成 0，和 sa_cases 同一口径，否则生产间隙会把芯片推出 72 mm。
填充环氧 k=0.8。python 场，冷边 45 C，底面绝热，不是 HotSpot，也不是 ATSim。
"""
import json
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
OUT = Path("/home/czh/任务仿真产生功耗文件工具/论文证据/填充场叠加/复现/accel_56789")

# (标签, 芯粒数, 布局文件, 去掉的芯片名)
CASES = (
    ("5 芯粒", 5, "sa_heat_vs_coupled_pkg_2cpu_1gpu_2cim/layout_start.json", ()),
    ("6 芯粒", 6, "sa_heat_vs_coupled_pkg_2cpu_2gpu_2cim/layout_start.json", ()),
    ("7 芯粒", 7, "sa_heat_vs_coupled_pkg_3cpu_1gpu_3cim/layout_start.json", ()),
    ("8 芯粒", 8, "sa_heat_vs_coupled_pkg_4cpu_2gpu_4cim/layout_start.json",
     ("CIM_basic", "CIM_wang")),
    ("9 芯粒", 9, "sa_heat_vs_coupled_pkg_4cpu_2gpu_4cim/layout_start.json",
     ("CIM_basic",)),
)


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


def chip_layer(x):
    return np.asarray(x[N * N:2 * N * N]).reshape(N, N).copy()


def peak(field):
    return float(field.max() - 273.15)


def anneal(chips0, intp, g, t_start, mode, fac0):
    import sa_thermal_feedback as sa
    rng = random.Random(sa_same.SEED)
    cur = deepcopy(chips0)
    best = deepcopy(chips0)
    t_cur = t_best = t_start
    best_field = None
    temp = sa_same.T0
    asm = fac_s = sol_s = 0.0
    n_accept = n_legal = 0
    for step in range(1, sa_same.STEPS + 1):
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
        t0 = time.perf_counter()
        mat, _ = im._edge_system(hw.paint(dies, 0.8))
        asm += time.perf_counter() - t0
        rhs = sa_same.rhs_of(hw.power_map(dies), g)
        t0 = time.perf_counter()
        if mode == "exact":
            fac = splu(mat.tocsc())
            fac_s += time.perf_counter() - t0
            t0 = time.perf_counter()
            field = chip_layer(fac.solve(rhs))
            sol_s += time.perf_counter() - t0
        else:
            field = chip_layer(sa_same.iterate(mat, rhs, fac0))
            sol_s += time.perf_counter() - t0
        tf = peak(field)
        took = tf <= t_cur or rng.random() < np.exp(-(tf - t_cur) / temp)
        if took:
            cur = trial
            t_cur = tf
            n_accept += 1
            if tf < t_best:
                best = deepcopy(trial)
                t_best = tf
                best_field = field
        temp *= sa_same.ALPHA
        if step % 50 == 0:
            print(f"  {mode:5s} step {step:3d}  当前 {t_cur:.2f}  最好 {t_best:.2f}",
                  flush=True)
    return {
        "assemble_s": asm, "factor_s": fac_s, "solve_s": sol_s,
        "accept": n_accept, "legal": n_legal, "priced_best_C": t_best,
        "dies": sa_same.to_dies(best), "field": best_field, "chips": best,
    }


def exact_field(chips, g):
    dies = sa_same.to_dies(chips)
    mat, _ = im._edge_system(hw.paint(dies, 0.8))
    rhs = sa_same.rhs_of(hw.power_map(dies), g)
    t0 = time.perf_counter()
    field = chip_layer(splu(mat.tocsc()).solve(rhs))
    return field, time.perf_counter() - t0


def grid_err(priced, exact):
    err = priced - exact
    return {
        "max_K": float(err.max()),
        "min_K": float(err.min()),
        "max_abs_K": float(np.max(np.abs(err))),
        "mean_abs_K": float(np.mean(np.abs(err))),
        "rmse_K": float(np.sqrt(np.mean(err ** 2))),
        "peak_price_minus_exact_K": float(priced.max() - exact.max()),
    }


def check():
    import sa_thermal_feedback as sa
    for label, n, rel, drop in CASES:
        chips, intp = load(rel, drop)
        dies = sa_same.to_dies(chips)
        overlap = 0
        seen = {}
        for d in dies:
            for i, j in hw.cells(d):
                overlap += (i, j) in seen
                seen[(i, j)] = d["name"]
        power = sum(c.power_W for c in chips)
        print(f"{label}  实际 {len(chips)} 颗  功耗 {power:.2f} W  "
              f"合法 {sa.legal(chips, intp, intp)}  足迹重叠格 {overlap}  "
              f"去掉 {','.join(drop) or '无'}", flush=True)


def dump(rows, arrays):
    note = (
        "5、6、7、8、9 芯粒。种子 22，200 步，T0=2，alpha=0.97，固定 8 次修正。"
        "间隙改成 0。填充环氧 k=0.8。加速只在起点分解一次，每步装配后固定 8 次修正，"
        "接受时不重分解。精确每步都重新分解。两边各自按自己的峰值接受。"
        "iter_exact_C 是加速选出的布局再用全部分解重算的峰值，时间在 check_s，不算进退火。"
        "5、6、7 用 sa_heat_vs_coupled 的 lbfgs 起点。8、9 没有现成封装，"
        "从 10 芯粒起点里去掉功耗最小的 CIM。python 场，冷边 45 C，底面绝热，"
        "100 um 硅中介层，TIM k=4，不是 HotSpot，也不是 ATSim。"
    )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "accel_56789.json").write_text(json.dumps({
        "note": note,
        "steps": sa_same.STEPS,
        "seed": sa_same.SEED,
        "n_iter": sa_same.N_ITER,
        "cases": rows,
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    np.savez(OUT / "accel_56789_fields.npz", **arrays)


def plot(rows, arrays):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.patches import Rectangle

    font_manager.fontManager.addfont(
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    plt.rcParams["font.family"] = "Noto Sans CJK JP"
    plt.rcParams["axes.unicode_minus"] = False
    n = len(rows)
    temps = []
    for i in range(n):
        temps.append(arrays[f"exact_{i}"] - 273.15)
        temps.append(arrays[f"rescored_{i}"] - 273.15)
    lo = min(float(t.min()) for t in temps)
    hi = max(float(t.max()) for t in temps)
    err_lim = max(float(np.max(np.abs(arrays[f"err_{i}"]))) for i in range(n))
    err_lim = max(err_lim, 1e-4)

    fig = plt.figure(figsize=(13.6, 2.45 * n))
    gs = fig.add_gridspec(n, 5, width_ratios=[1, 1, 0.045, 1, 0.045],
                          wspace=0.22, hspace=0.55)
    axes = np.empty((n, 3), dtype=object)
    for i in range(n):
        axes[i, 0] = fig.add_subplot(gs[i, 0])
        axes[i, 1] = fig.add_subplot(gs[i, 1])
        axes[i, 2] = fig.add_subplot(gs[i, 3])
    for i, rec in enumerate(rows):
        panels = (
            (arrays[f"exact_{i}"] - 273.15, rec["exact_dies"],
             f"每步全部分解\n{rec['exact_best_C']:.2f} C   {rec['exact_total_s']:.0f} s",
             "temp"),
            (arrays[f"rescored_{i}"] - 273.15, rec["iter_dies"],
             f"8 次迭代选出，精确重算\n{rec['iter_exact_C']:.2f} C"
             f"   报价 {rec['iter_priced_C']:.2f} C", "temp"),
            (arrays[f"err_{i}"], rec["iter_dies"],
             f"同一布局  迭代减精确\n网格最大 {rec['grid_err']['max_abs_K']:.3f} K"
             f"   峰值 {rec['grid_err']['peak_price_minus_exact_K']:+.3f} K", "err"),
        )
        for ax, (fld, dies, title, kind) in zip(axes[i], panels):
            if kind == "temp":
                ax.imshow(fld[::-1], origin="lower", cmap="viridis",
                          vmin=lo, vmax=hi, extent=[0, 72, 0, 72])
                ec = "white"
            else:
                ax.imshow(fld[::-1], origin="lower", cmap="coolwarm",
                          vmin=-err_lim, vmax=err_lim, extent=[0, 72, 0, 72])
                ec = "black"
            for c in dies:
                ax.add_patch(Rectangle((c["j0"] * 1.125, c["ib0"] * 1.125),
                                       c["w"] * 1.125, c["h"] * 1.125,
                                       fill=False, ec=ec, lw=0.45))
            ax.set_title(title, fontsize=9)
            ax.tick_params(labelsize=7)
        axes[i, 0].set_ylabel(rec["label"], fontsize=9)
    fig.colorbar(axes[0, 0].images[0], cax=fig.add_subplot(gs[:, 2])).set_label("C")
    fig.colorbar(axes[0, 2].images[0], cax=fig.add_subplot(gs[:, 4])).set_label("K")
    fig.subplots_adjust(top=0.93, left=0.06, right=0.97)
    fig.suptitle("材料场加速对每步全部分解，各自退火 200 步，种子 22。"
                 "加速只在开局分解一次，之后每步固定 8 次", fontsize=11)
    figdir = Path("/home/czh/任务仿真产生功耗文件工具/论文证据/填充场叠加/图")
    figdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figdir / "加速对比_温度场.png", dpi=130)

    fig, ax = plt.subplots(figsize=(8.2, 4.6), constrained_layout=True)
    y = np.arange(n)
    exact_s = [r["exact_total_s"] for r in rows]
    iter_s = [r["iter_total_s"] for r in rows]
    ax.barh(y + 0.16, exact_s, height=0.32, color="#4c72b0", label="每步全部分解")
    ax.barh(y - 0.16, iter_s, height=0.32, color="#dd8452",
            label="开局分解一次，每步 8 次")
    for i, rec in enumerate(rows):
        ax.text(rec["exact_total_s"] + 1.5, i + 0.16, f"{rec['exact_total_s']:.0f} s",
                va="center", fontsize=8)
        ax.text(rec["iter_total_s"] + 1.5, i - 0.16,
                f"{rec['iter_total_s']:.0f} s   {rec['speedup']:.2f} 倍",
                va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=9)
    ax.set_xlabel("200 步退火用时，s。含装配。复核那一次不算进去")
    ax.set_xlim(0, max(exact_s) * 1.32)
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title("材料场加速相对每步重新分解", fontsize=12)
    fig.savefig(figdir / "加速对比_用时.png", dpi=130)


def main():
    import sa_thermal_feedback as sa
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    arrays = {}
    for idx, (label, n, rel, drop) in enumerate(CASES):
        print(f"\n== {label}", flush=True)
        chips0, intp = load(rel, drop)
        if len(chips0) != n:
            raise SystemExit(f"{label} 实际 {len(chips0)} 颗")
        if not sa.legal(chips0, intp, intp):
            raise SystemExit(f"{label} 起点不合法")
        dies0 = sa_same.to_dies(chips0)
        t0 = time.perf_counter()
        mat0, g = im._edge_system(hw.paint(dies0, 0.8))
        assemble0 = time.perf_counter() - t0
        t0 = time.perf_counter()
        fac0 = splu(mat0.tocsc())
        factor0 = time.perf_counter() - t0
        start = chip_layer(fac0.solve(sa_same.rhs_of(hw.power_map(dies0), g)))
        t_start = peak(start)
        print(f"  起点 {t_start:.2f} C  分解 {factor0:.2f} s  "
              f"功耗 {sum(c.power_W for c in chips0):.2f} W", flush=True)

        full = anneal(chips0, intp, g, t_start, "exact", fac0)
        full_field = full.pop("field")
        full.pop("chips")
        if full_field is None:
            full_field = start
        it = anneal(chips0, intp, g, t_start, "iter", fac0)
        priced = it.pop("field")
        it_chips = it.pop("chips")
        if priced is None:
            priced = start
            it_chips = chips0
        rescored, check_s = exact_field(it_chips, g)
        err = grid_err(priced, rescored)
        total_exact = full["assemble_s"] + full["factor_s"] + full["solve_s"]
        total_iter = factor0 + it["assemble_s"] + it["solve_s"]
        rec = {
            "label": label,
            "n_chips": len(chips0),
            "source": rel,
            "dropped": list(drop),
            "chips": [{"name": c.name, "kind": c.kind, "power_W": c.power_W,
                       "area_mm2": c.area_mm2} for c in chips0],
            "power_W": float(sum(c.power_W for c in chips0)),
            "start_C": t_start,
            "factor_once_s": factor0,
            "assemble0_s": assemble0,
            "exact": {k: full[k] for k in
                      ("assemble_s", "factor_s", "solve_s", "accept", "legal",
                       "priced_best_C")},
            "exact_total_s": total_exact,
            "exact_best_C": peak(full_field),
            "exact_dies": full["dies"],
            "iter": {k: it[k] for k in
                     ("assemble_s", "solve_s", "accept", "legal", "priced_best_C")},
            "iter_total_s": total_iter,
            "iter_priced_C": it["priced_best_C"],
            "iter_exact_C": peak(rescored),
            "iter_dies": it["dies"],
            "check_s": check_s,
            "grid_err": err,
            "speedup": total_exact / total_iter,
        }
        rows.append(rec)
        arrays[f"start_{idx}"] = start
        arrays[f"exact_{idx}"] = full_field
        arrays[f"priced_{idx}"] = priced
        arrays[f"rescored_{idx}"] = rescored
        arrays[f"err_{idx}"] = priced - rescored
        dump(rows, arrays)
        print(f"  精确 {rec['exact_best_C']:.2f} C  {total_exact:.1f} s", flush=True)
        print(f"  加速报价 {rec['iter_priced_C']:.2f}  复核 {rec['iter_exact_C']:.2f}  "
              f"峰值差 {err['peak_price_minus_exact_K']:+.3f} K  "
              f"网格最大 {err['max_abs_K']:.3f} K  {total_iter:.1f} s  "
              f"{rec['speedup']:.2f} 倍", flush=True)
    plot(rows, arrays)
    print("写完", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--check":
        check()
    else:
        main()
