#!/usr/bin/env python3
"""同一张冻住的布局上，拆开各颗倍率，找温度和应力互相牺牲的点。

布局是 11 颗环氧的内环最优，不再移动、不再换芯片。倍率步长 0.5，
只动动态功耗大的三颗：Xeon_Ca、T2_Cb、GPU_3070_Ga。另外两颗 GPU
动态功耗小，固定在 1.0。每个取样点单独跑一次 ATSim。

认的数和门槛一样：芯粒层最热一格，正常芯片内部应力。
不是块平均，也不是整包最大应力。
"""

from __future__ import annotations

import json
import os
import sys
from itertools import product
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

import run_closed_loop_multidie as cl
import run_closed_loop_suite as suite

BEST = ROOT / "热仿真流水线" / "runs" / "closed_loop_suite" / "n11_epoxy" / "layouts" / "best.json"
CHECK = ROOT / "热仿真流水线" / "runs" / "closed_loop_suite" / "n11_epoxy" / "atsim" / "best_s1_result.json"
OUT = ROOT / "热仿真流水线" / "runs" / "pareto_freq"
FONT = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
# 动态功耗大的三颗才拆开。另外两颗 GPU 固定 1.0。
FREE = ("Ca", "Cb", "Ga")
LEVELS = (0.5, 1.0, 1.5, 2.0, 2.5)


def style() -> None:
    font_manager.fontManager.addfont(FONT)
    plt.rcParams["font.family"] = "Noto Sans CJK JP"
    plt.rcParams["axes.unicode_minus"] = False


def scales_of(combo: tuple[float, float, float]) -> dict[str, float]:
    out = {"Gb": 1.0, "Gc": 1.0}
    out.update(zip(FREE, combo))
    return out


def tag_of(combo: tuple[float, float, float]) -> str:
    return "s_" + "_".join(f"{slot}{value:.1f}" for slot, value in zip(FREE, combo))


def evaluate(base: list[dict], combo: tuple[float, float, float]) -> dict:
    tag = tag_of(combo)
    path = OUT / "atsim" / f"{tag}_result.json"
    scales = scales_of(combo)
    if not path.exists():
        tuned = cl.with_scales(base, scales)
        suite.run_atsim(tuned, "epoxy", OUT / "atsim" / f"{tag}.json", path)
    result = json.loads(path.read_text(encoding="utf-8"))
    tuned = cl.with_scales(base, scales)
    by = {d["slot"]: d for d in tuned}
    return {
        "tag": tag,
        "scales": scales,
        "die_cell_max_C": result["die_cell_max_C"],
        "vm_MPa": result["vm_normal_chip_max_MPa"],
        "vm_name": result["vm_normal_chip_name"],
        "power_W": round(sum(d["power_W"] for d in tuned), 3),
        "cpu_f_MHz": {slot: by[slot]["f_MHz"] * by[slot]["s"] for slot in ("Ca", "Cb")},
        "gpu_f_MHz": {slot: by[slot]["f_MHz"] * by[slot]["s"] for slot in ("Ga", "Gb", "Gc")},
    }


def dominates(a: dict, b: dict) -> bool:
    """温度低、应力低都是更好。"""
    better = a["die_cell_max_C"] <= b["die_cell_max_C"] and a["vm_MPa"] <= b["vm_MPa"]
    strict = a["die_cell_max_C"] < b["die_cell_max_C"] or a["vm_MPa"] < b["vm_MPa"]
    return better and strict


def draw(rows: list[dict]) -> None:
    style()
    front = [r for r in rows if r["on_front"]]
    front.sort(key=lambda r: r["die_cell_max_C"])
    fig, ax = plt.subplots(figsize=(8.6, 6.2))
    ax.scatter(
        [r["die_cell_max_C"] for r in rows if not r["on_front"]],
        [r["vm_MPa"] for r in rows if not r["on_front"]],
        s=18, color="#c5c1b9", label="被支配", zorder=2,
    )
    ax.plot(
        [r["die_cell_max_C"] for r in front],
        [r["vm_MPa"] for r in front],
        "-o", color="#1f4e79", ms=5, label="非支配", zorder=3,
    )
    marked = {front[0]["tag"], front[-1]["tag"]}
    knee = min(front, key=lambda r: r["die_cell_max_C"] / 90 + r["vm_MPa"] / 200)
    marked.add(knee["tag"])
    for row in rows:
        if row["tag"] not in marked:
            continue
        sc = row["scales"]
        text = f"Ca {sc['Ca']:.1f}  Cb {sc['Cb']:.1f}\nGa {sc['Ga']:.1f}"
        ax.annotate(
            text, (row["die_cell_max_C"], row["vm_MPa"]),
            textcoords="offset points", xytext=(8, 6), fontsize=8, color="#1f4e79",
        )
    ax.axvline(90, color="#888", lw=0.7)
    ax.axhline(200, color="#888", lw=0.7)
    ax.text(90.4, ax.get_ylim()[0], "", fontsize=1)
    ax.set_xlabel("芯粒层最热一格，C")
    ax.set_ylabel("正常芯片内部应力，MPa")
    ax.set_title("11 颗环氧，布局冻住，只拆三颗的倍率", loc="left")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "pareto.png", dpi=140)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "atsim").mkdir(exist_ok=True)
    base = json.loads(BEST.read_text(encoding="utf-8"))
    combos = list(product(LEVELS, repeat=len(FREE)))
    rows = []
    for i, combo in enumerate(combos, 1):
        row = evaluate(base, combo)
        rows.append(row)
        print(
            f"{i:03d}/{len(combos)} {row['tag']}  "
            f"T {row['die_cell_max_C']:.2f}  vm {row['vm_MPa']:.1f}",
            flush=True,
        )
    check = json.loads(CHECK.read_text(encoding="utf-8"))
    nominal = next(r for r in rows if r["tag"] == tag_of((1.0, 1.0, 1.0)))
    gap_t = abs(nominal["die_cell_max_C"] - check["die_cell_max_C"])
    gap_s = abs(nominal["vm_MPa"] - check["vm_normal_chip_max_MPa"])
    if gap_t > 0.05 or gap_s > 0.05:
        raise RuntimeError(f"s=1 对不上已有结果: {gap_t:.4f} C, {gap_s:.4f} MPa")
    for row in rows:
        row["on_front"] = not any(dominates(other, row) for other in rows if other is not row)
    (OUT / "samples.json").write_text(
        json.dumps({
            "layout": "n11_epoxy 内环最优，冻住",
            "fill": "epoxy",
            "free_slots": FREE,
            "fixed": {"Gb": 1.0, "Gc": 1.0},
            "levels": LEVELS,
            "objectives": "芯粒层最热一格越低越好，正常芯片内部应力越低越好",
            "nominal_matches_saved_run": {"dT_C": gap_t, "dMPa": gap_s},
            "n_front": sum(r["on_front"] for r in rows),
            "rows": rows,
        }, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    draw(rows)
    front = sorted((r for r in rows if r["on_front"]), key=lambda r: r["die_cell_max_C"])
    print(f"front {len(front)} / {len(rows)}", flush=True)
    for row in front:
        sc = row["scales"]
        print(
            f"  T {row['die_cell_max_C']:.2f}  vm {row['vm_MPa']:.1f}  "
            f"Ca {sc['Ca']:.1f} Cb {sc['Cb']:.1f} Ga {sc['Ga']:.1f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
