#!/usr/bin/env python3
"""micro150 重跑的热图。左是 16 网格场叠加的初始布局，右是种子 23 重跑的最优。

温度场都是 64 网格重评。芯粒框按中心减半宽画，和 TAP 写 flp 的口径一致。
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
RUN = HERE / "runs" / "micro150_16_field_seed23"
ROWS = 64
SIDE = 60.0


def field(path):
    vals = [float(line) for line in Path(path).read_text().splitlines() if line.strip()]
    return np.array(vals).reshape(ROWS, ROWS) - 273.15


def main() -> None:
    init = field(RUN / "init_64.grid")
    best = field(RUN / "best_64.grid")
    start = json.loads((RUN / "init.json").read_text())
    winner = json.loads((RUN / "best.json").read_text())
    lo = min(float(init.min()), float(best.min()))
    hi = max(float(init.max()), float(best.max()))
    panels = (
        (init, start, f"初始布局\n{init.max():.2f} °C"),
        (best, winner, f"场叠加重跑，种子 23\n{best.max():.2f} °C"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.8))
    for ax, (fld, chips, title) in zip(axes, panels):
        im = ax.imshow(fld[::-1], origin="lower", cmap="inferno",
                       vmin=lo, vmax=hi, extent=[0, SIDE, 0, SIDE])
        for c in chips:
            # TAP 的 x、y 是芯粒中心，flp 里的左下角是中心减半宽。
            ax.add_patch(Rectangle(
                (c["x_mm"] - c["w_mm"] / 2, c["y_mm"] - c["h_mm"] / 2),
                c["w_mm"], c["h_mm"], fill=False, ec="white", lw=1.0))
        ax.set_xlim(0, SIDE)
        ax.set_ylim(0, SIDE)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("mm")
    axes[0].set_ylabel("mm")
    fig.colorbar(im, ax=list(axes), shrink=0.85, pad=0.02).set_label("°C")
    fig.suptitle("micro150，16 网格场叠加，初始布局不变，种子 23 重跑 120 步", fontsize=12)
    out = HERE / "图" / "micro150_16网格_重跑.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"写到 {out}")


if __name__ == "__main__":
    main()
