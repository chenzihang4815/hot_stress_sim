#!/usr/bin/env python3
"""10 步的误差场用同一张图、同一套透明度叠起来。

每一步是 frozen 减 exact，也就是材料图停在起点时多估出来的温度。
10 张场各用 0.1 的透明度叠到同一张坐标上，颜色深的地方是多步都在估高的地方。
芯粒框按步从浅到深，标出芯片挪过的位置。第 6 步被拒绝，框用虚线。

数据是 复现/ten_frozen_grids.npz 和 ten_frozen_layouts.json，
2026-11-17 重跑，和 2026-11-14 的存档逐字节相同。
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import Normalize
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
    err = np.stack([z[f"step{s}_err"] for s in range(1, 11)])
    span = float(np.max(np.abs(err)))
    norm = Normalize(vmin=-span, vmax=span)
    cmap = plt.get_cmap("coolwarm")

    fig, ax = plt.subplots(figsize=(7.6, 6.8))
    blank = np.zeros_like(err[0])
    base = ax.imshow(blank, origin="lower", extent=(0, SIDE, 0, SIDE),
                     cmap=cmap, norm=norm)
    # 误差场的第 0 行是 HotSpot 的顶行，芯粒框的 ib0 从底边往上数，叠之前翻过来。
    for field in err:
        ax.imshow(cmap(norm(field[::-1])), origin="lower",
                  extent=(0, SIDE, 0, SIDE), alpha=0.1)

    for s, (lay, rec) in enumerate(zip(layouts, rows), start=1):
        shade = 0.25 + 0.7 * (s - 1) / 9
        for d in lay:
            ax.add_patch(Rectangle(
                (d["j0"] * CELL, d["ib0"] * CELL), d["w"] * CELL, d["h"] * CELL,
                fill=False, ec=(shade, shade, shade), lw=0.7,
                ls="--" if not rec["accept"] else "-"))

    ax.set_xlim(0, SIDE)
    ax.set_ylim(0, SIDE)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("10 步误差叠在一张上，每步透明度 0.1\n"
                 "红是材料图停在起点时估高了，框从浅到深是第 1 到第 10 步",
                 fontsize=11)
    fig.colorbar(base, ax=ax, fraction=0.046, pad=0.02).set_label(
        "frozen − exact (K)")
    fig.tight_layout()
    out = HERE / "图" / "十步叠加.png"
    fig.savefig(out, dpi=140)
    print(f"写到 {out}")
    print(f"10 步里整图最大 {err.max():+.2f} K，最小 {err.min():+.2f} K")


if __name__ == "__main__":
    main()
