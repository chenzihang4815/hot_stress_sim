#!/usr/bin/env python3
"""10 步的误差热图画在一张图里，每步一格。

每格是这一步的 frozen 减 exact，也就是材料图停在起点时多估出来的温度。
10 格共用一个色标。芯粒框是这一步的布局，被拒绝的那一步用虚线。

数据是 复现/ten_frozen_grids.npz 和 ten_frozen_layouts.json，
2026-11-17 重跑，和 2026-11-14 的存档逐字节相同。
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
import numpy as np

font_manager.fontManager.addfont(
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
DATA = HERE / "复现"
SIDE = 72.0
CELL = SIDE / 64


def main() -> None:
    z = np.load(DATA / "ten_frozen_grids.npz")
    layouts = json.loads((DATA / "ten_frozen_layouts.json").read_text())
    rows = json.loads((DATA / "ten_frozen.json").read_text())["rows"]
    err = [z[f"step{s}_err"][::-1] for s in range(1, 11)]
    # 90% 的格子误差在 0.6 K 以内，31 K 的峰值把色标拉满，浅的地方就看不见。
    # 色标收到 2 K，超出的格子画成两端的深色。
    span = 2.0

    fig, axes = plt.subplots(2, 5, figsize=(16, 7.2))
    im = None
    for ax, field, lay, rec in zip(axes.ravel(), err, layouts, rows):
        im = ax.imshow(field, origin="lower", extent=(0, SIDE, 0, SIDE),
                       cmap="coolwarm", vmin=-span, vmax=span)
        for d in lay:
            ax.add_patch(Rectangle(
                (d["j0"] * CELL, d["ib0"] * CELL), d["w"] * CELL, d["h"] * CELL,
                fill=False, ec="black", lw=0.6,
                ls="--" if not rec["accept"] else "-"))
        mark = "" if rec["accept"] else "，弃"
        ax.set_title(f"第 {rec['step']} 步{mark}  {rec['peak_err_K']:+.1f} K",
                     fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    bar = fig.colorbar(im, ax=axes, fraction=0.015, pad=0.01, extend="both")
    bar.set_label("frozen − exact (K)，超出 ±2 的画成两端深色")
    fig.suptitle("10 步，材料图停在起点的误差，每格是这一步的布局", fontsize=13)
    out = HERE / "图" / "十步每步误差.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"写到 {out}")


if __name__ == "__main__":
    main()
