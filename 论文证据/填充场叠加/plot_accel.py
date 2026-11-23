#!/usr/bin/env python3
"""材料场加速选出的布局上，报价相对精确重算的误差。"""

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


def main() -> None:
    rows = json.loads((HERE / "复现/accel_56789/accel_56789.json").read_text())["cases"]
    rows += json.loads((HERE / "复现/accel_11/accel_56789.json").read_text())["cases"]
    labels = [r["label"] for r in rows]
    x = np.arange(len(rows))
    grid = [r["grid_err"]["max_abs_K"] for r in rows]
    peak = [r["grid_err"]["peak_price_minus_exact_K"] for r in rows]

    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    w = 0.36
    ax.bar(x - w / 2, grid, width=w, color="#542788", label="网格最大误差")
    colors = ["#b35806" if v < 0 else "#1b7837" for v in peak]
    ax.bar(x + w / 2, peak, width=w, color=colors, label="峰值差，报价减精确")
    ax.axhline(0, color="black", lw=0.6)
    for i, (g, p) in enumerate(zip(grid, peak)):
        ax.text(i - w / 2, g + 0.05, f"{g:.2f}", ha="center", fontsize=8)
        ax.text(i + w / 2, p + (0.08 if p >= 0 else -0.16), f"{p:+.2f}",
                ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("K")
    ax.set_ylim(-1.15, 3.1)
    ax.legend(fontsize=9)
    ax.set_title("加速选出的布局，用全部分解重算后的误差\n"
                 "峰值差为负是报价偏低", fontsize=12)
    fig.tight_layout()
    out = HERE / "图" / "加速对比_误差.png"
    fig.savefig(out, dpi=140)
    print(f"写到 {out}")


if __name__ == "__main__":
    main()
