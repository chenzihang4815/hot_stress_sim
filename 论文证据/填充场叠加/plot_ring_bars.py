#!/usr/bin/env python3
"""外伸环算子上，材料场加速相对每步全部分解的用时和误差。"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

font_manager.fontManager.addfont(
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent


def load():
    rows = []
    for label in ("n5", "n8", "n11"):
        path = HERE / "复现" / f"ring_{label}" / f"ring_{label}.json"
        rows.append(json.loads(path.read_text()))
    return rows


def main() -> None:
    rows = load()
    labels = [f"{r['n_chips']} 芯粒" for r in rows]
    x = np.arange(len(rows))
    w = 0.34

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.8))

    ax = axes[0]
    exact = [r["exact_total_s"] for r in rows]
    accel = [r["iter_total_s"] for r in rows]
    ax.bar(x - w / 2, exact, width=w, color="#4c72b0", label="每步全部分解")
    ax.bar(x + w / 2, accel, width=w, color="#dd8452", label="材料场加速")
    for i, rec in enumerate(rows):
        ax.text(i - w / 2, rec["exact_total_s"] + 1.2, f"{rec['exact_total_s']:.0f} s",
                ha="center", fontsize=9)
        ax.text(i + w / 2, rec["iter_total_s"] + 1.2,
                f"{rec['iter_total_s']:.0f} s\n{rec['speedup']:.2f} 倍",
                ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("200 步退火用时，s")
    ax.set_ylim(0, max(exact) * 1.32)
    ax.legend(fontsize=9, loc="upper left")
    ax.set_title("用时", fontsize=12)

    ax = axes[1]
    grid = [r["grid_max_abs_K"] for r in rows]
    peak = [r["peak_price_minus_exact_K"] for r in rows]
    ax.bar(x - w / 2, grid, width=w, color="#542788", label="网格最大误差")
    ax.bar(x + w / 2, peak, width=w, color="#1b7837", label="峰值差，报价减精确重算")
    ax.axhline(0, color="black", lw=0.6)
    for i, (g, p) in enumerate(zip(grid, peak)):
        ax.text(i - w / 2, g + 0.002, f"{g:.3f}", ha="center", fontsize=9)
        ax.text(i + w / 2, 0.004, f"{p:+.4f}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("K")
    ax.set_ylim(-0.01, 0.1)
    ax.legend(fontsize=9, loc="upper left")
    ax.set_title("加速选出的布局，报价相对精确重算", fontsize=12)

    fig.suptitle("外伸环算子，材料场加速对每步全部分解。只接受温度更低的布局，种子 22，200 步",
                 fontsize=12)
    fig.tight_layout()
    out = HERE / "图" / "外伸环_加速对比.png"
    fig.savefig(out, dpi=140)
    print(f"写到 {out}")


if __name__ == "__main__":
    main()
