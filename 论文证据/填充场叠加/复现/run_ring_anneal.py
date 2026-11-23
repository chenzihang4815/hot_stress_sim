"""面积接近的芯粒，材料场加速对每步全部分解。求解换成带外伸环的四层算子。

芯粒取闭环流水线放大过面积的那一套，n5、n8、n11 三套。
最小的也有 11 mm2 以上，框按芯粒真实毫米尺寸画，不按格子取整。

退火只按峰值温度接受，只接受温度更低的布局，不算应力，不算线长。
种子 22，200 步，T0=2，alpha=0.97。加速只在起点分解一次，每步装配后
固定修正 8 次，接受时不重分解。精确每步都重新分解。两边各自按自己的
峰值接受，布局会分开。间隙 0。填充环氧 k=0.8。

算子是 ring_solver：芯片、TIM、铜均热片、铜热沉，网格到 72 mm，
72 mm 以外补 8 个外伸节点。没有中介层，没有冷边。不是 HotSpot 本身，
但和 HotSpot 同一套硅填充时峰值差 0.64 K。
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
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_hw_area_field as hw  # noqa: E402
import sa_same  # noqa: E402
import ring_solver as rs  # noqa: E402

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
HERE = Path("/home/czh/任务仿真产生功耗文件工具/论文证据/填充场叠加")
CASES = (
    ("n5", "热仿真流水线/runs/closed_loop_suite/n5_epoxy/hotspot_init/layout.json"),
    ("n8", "热仿真流水线/runs/closed_loop_suite/n8_epoxy/hotspot_init/layout.json"),
    ("n11", "热仿真流水线/runs/closed_loop_suite/n11_epoxy/hotspot_init/layout.json"),
)


def load(src):
    import sa_thermal_feedback as sa
    raw = json.loads((ROOT / src).read_text())
    chips = [sa.Chip(c["name"], c["kind"], c["x_m"], c["y_m"], c["w_m"], c["h_m"],
                     c["power_W"], c["area_mm2"]) for c in raw]
    side = rs.SIDE
    sa.GAP_M = sa.GAP_GPU_M = 0.0
    chips = sa.spread_to_gap(chips, side, side)
    for c in chips:
        c.x = min(max(0.0, c.x), max(0.0, side - c.w))
        c.y = min(max(0.0, c.y), max(0.0, side - c.h))
    return chips


def assemble(dies):
    return rs.build(hw.paint(dies, rs.K_FILL))


def solve_field(mat, dies, g_amb, fac):
    rhs = rs.rhs_of(hw.power_map(dies), g_amb)
    return rs.chip_layer(fac.solve(rhs))


def peak(field):
    return float(field.max() - 273.15)


def boxes_of(chips):
    return [{"name": c.name, "x_mm": c.x * 1e3, "y_mm": c.y * 1e3,
             "w_mm": c.w * 1e3, "h_mm": c.h * 1e3} for c in chips]


def anneal(chips0, g_amb, t_start, mode, fac0):
    import sa_thermal_feedback as sa
    rng = random.Random(sa_same.SEED)
    side = rs.SIDE
    cur = deepcopy(chips0)
    best = deepcopy(chips0)
    t_cur = t_best = t_start
    best_field = None
    temp = sa_same.T0
    asm = fac_s = sol_s = 0.0
    n_accept = n_legal = 0
    for step in range(1, sa_same.STEPS + 1):
        trial, _ = sa.propose(cur, side, side, rng)
        tries = 1
        while not sa.legal(trial, side, side) and tries < 25:
            trial, _ = sa.propose(cur, side, side, rng)
            tries += 1
        if not sa.legal(trial, side, side):
            temp *= sa_same.ALPHA
            continue
        n_legal += 1
        dies = sa_same.to_dies(trial)
        t0 = time.perf_counter()
        mat, _ = assemble(dies)
        asm += time.perf_counter() - t0
        rhs = rs.rhs_of(hw.power_map(dies), g_amb)
        t0 = time.perf_counter()
        if mode == "exact":
            fac = splu(mat.tocsc())
            fac_s += time.perf_counter() - t0
            t0 = time.perf_counter()
            field = rs.chip_layer(fac.solve(rhs))
            sol_s += time.perf_counter() - t0
        else:
            field = rs.chip_layer(sa_same.iterate(mat, rhs, fac0))
            sol_s += time.perf_counter() - t0
        tf = peak(field)
        took = tf <= t_cur
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
        "dies": sa_same.to_dies(best), "boxes": boxes_of(best),
        "chips": best, "field": best_field,
    }


def covered_mask(box_sets):
    cell = rs.SIDE * 1e3 / rs.N
    covered = np.zeros((rs.N, rs.N), dtype=bool)
    for boxes in box_sets:
        for c in boxes:
            j0 = int(round(c["x_mm"] / cell))
            ib0 = int(round(c["y_mm"] / cell))
            w = max(1, int(round(c["w_mm"] / cell)))
            h = max(1, int(round(c["h_mm"] / cell)))
            covered[ib0:ib0 + h, j0:j0 + w] = True
    return covered


def plot(label, start, full_field, iter_field, box_sets, rec):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.patches import Rectangle

    font_manager.fontManager.addfont(
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    plt.rcParams["font.family"] = "Noto Sans CJK JP"
    plt.rcParams["axes.unicode_minus"] = False

    temps = [start - 273.15, full_field - 273.15, iter_field - 273.15]
    inside = np.concatenate([t[covered_mask(box_sets)] for t in temps])
    lo = float(np.percentile(inside, 2))
    hi = max(float(t.max()) for t in temps)
    n = rec["n_chips"]
    panels = (
        (temps[0], box_sets[0], f"初始布局\n{rec['start_C']:.2f} C"),
        (temps[1], box_sets[1],
         f"每步全部分解，最优\n{rec['exact_best_C']:.2f} C   {rec['exact_total_s']:.0f} s"),
        (temps[2], box_sets[2],
         f"材料场加速，最优\n报价 {rec['iter_priced_C']:.2f} C"
         f"   精确重算 {rec['iter_exact_C']:.2f} C"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.8))
    for ax, (fld, boxes, title) in zip(axes, panels):
        ims = ax.imshow(fld[::-1], origin="lower", cmap="inferno",
                        vmin=lo, vmax=hi, extent=[0, 72, 0, 72])
        for c in boxes:
            ax.add_patch(Rectangle((c["x_mm"], c["y_mm"]), c["w_mm"], c["h_mm"],
                                   fill=False, ec="white", lw=1.1))
        ax.set_xlim(0, 72)
        ax.set_ylim(0, 72)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("mm")
    axes[0].set_ylabel("mm")
    fig.colorbar(ims, ax=list(axes), shrink=0.85, pad=0.02).set_label("C")
    areas = [c["area_mm2"] for c in rec["chips"]]
    fig.suptitle(f"{n} 芯粒，面积 {min(areas):.0f} 到 {max(areas):.0f} mm²。"
                 f"只接受温度更低的布局，种子 22，200 步", fontsize=12)
    out = HERE / "图" / f"外伸环_{n}芯粒_只降.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"写到 {out}", flush=True)


def one(label, src):
    import sa_thermal_feedback as sa
    out = HERE / "复现" / f"ring_{label}"
    out.mkdir(parents=True, exist_ok=True)
    chips0 = load(src)
    side = rs.SIDE
    if not sa.legal(chips0, side, side):
        raise SystemExit(f"{label} 起点不合法")
    dies0 = sa_same.to_dies(chips0)
    seen = {}
    overlap = 0
    for d in dies0:
        for i, j in hw.cells(d):
            overlap += (i, j) in seen
            seen[(i, j)] = d["name"]
    print(f"\n== {label}  {len(chips0)} 芯粒  "
          f"功耗 {sum(c.power_W for c in chips0):.2f} W  重叠格 {overlap}", flush=True)
    if overlap:
        raise SystemExit(f"{label} 起点有格子重叠")

    t0 = time.perf_counter()
    mat0, g_amb = assemble(dies0)
    print(f"  装配 {time.perf_counter() - t0:.2f} s", flush=True)
    t0 = time.perf_counter()
    fac0 = splu(mat0.tocsc())
    factor0 = time.perf_counter() - t0
    start = solve_field(mat0, dies0, g_amb, fac0)
    t_start = peak(start)
    print(f"  起点 {t_start:.2f} C  分解 {factor0:.2f} s", flush=True)

    full = anneal(chips0, g_amb, t_start, "exact", fac0)
    it = anneal(chips0, g_amb, t_start, "iter", fac0)
    full_field = full["field"] if full["field"] is not None else start
    priced = it["field"] if it["field"] is not None else start

    mat, _ = assemble(it["dies"])
    t0 = time.perf_counter()
    rescored = rs.chip_layer(splu(mat.tocsc()).solve(
        rs.rhs_of(hw.power_map(it["dies"]), g_amb)))
    check_s = time.perf_counter() - t0
    err = priced - rescored
    total_exact = full["assemble_s"] + full["factor_s"] + full["solve_s"]
    total_iter = factor0 + it["assemble_s"] + it["solve_s"]
    rec = {
        "label": label,
        "n_chips": len(chips0),
        "source": src,
        "solver": "四层网格加 8 个外伸节点，均热片 96 mm，热沉 96 mm，无冷边",
        "chips": [{"name": c.name, "kind": c.kind, "power_W": c.power_W,
                   "area_mm2": c.area_mm2} for c in chips0],
        "power_W": float(sum(c.power_W for c in chips0)),
        "start_C": t_start,
        "factor_once_s": factor0,
        "exact_total_s": total_exact,
        "exact_best_C": peak(full_field),
        "exact": {k: full[k] for k in
                  ("assemble_s", "factor_s", "solve_s", "accept", "legal")},
        "iter_total_s": total_iter,
        "iter_priced_C": it["priced_best_C"],
        "iter_exact_C": peak(rescored),
        "iter": {k: it[k] for k in ("assemble_s", "solve_s", "accept", "legal")},
        "check_s": check_s,
        "grid_max_abs_K": float(np.max(np.abs(err))),
        "peak_price_minus_exact_K": float(priced.max() - rescored.max()),
        "speedup": total_exact / total_iter,
        "steps": sa_same.STEPS,
        "seed": sa_same.SEED,
        "n_iter": sa_same.N_ITER,
        "acceptance": "只接受温度更低的布局",
    }
    (out / f"ring_{label}.json").write_text(
        json.dumps(rec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    boxes0 = boxes_of(chips0)
    for name, obj in (("start_dies", dies0), ("start_boxes", boxes0),
                      ("exact_best_dies", full["dies"]), ("exact_best_boxes", full["boxes"]),
                      ("iter_best_dies", it["dies"]), ("iter_best_boxes", it["boxes"])):
        (out / f"{name}.json").write_text(
            json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    np.savez(out / f"ring_{label}_fields.npz",
             start=start, exact_best=full_field, iter_priced=priced, iter_rescored=rescored)
    print(f"  精确 {rec['exact_best_C']:.2f} C  {total_exact:.1f} s", flush=True)
    print(f"  加速报价 {rec['iter_priced_C']:.2f}  复核 {rec['iter_exact_C']:.2f}  "
          f"网格最大 {rec['grid_max_abs_K']:.3f} K  {total_iter:.1f} s  "
          f"{rec['speedup']:.2f} 倍", flush=True)
    plot(label, start, full_field, rescored, (boxes0, full["boxes"], it["boxes"]), rec)
    return rec


def main():
    want = sys.argv[1:] or [c[0] for c in CASES]
    for label, src in CASES:
        if label in want:
            one(label, src)


if __name__ == "__main__":
    main()
