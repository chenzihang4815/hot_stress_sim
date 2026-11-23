#!/usr/bin/env python3
"""三条路径占总功耗的百分比。154.79 W，侧边和底面两个出口同时开。"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pathlib import Path

font_manager.fontManager.addfont(
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

# 百分比由 fill_with_both_exits.json 的热流除以 154.790098 W 得到，三条加起来是 100%。
FILLS = ["空气\n0.028", "有机\n0.5", "环氧\n0.8", "底部填充\n1.0", "硅\n130"]
TOP = [89.62, 77.11, 73.56, 71.65, 30.96]
BOTTOM = [7.33, 10.38, 10.34, 10.24, 5.04]
EDGE = [3.05, 12.50, 16.10, 18.11, 64.00]

HERE = Path(__file__).resolve().parent


def main() -> None:
    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    x = list(range(len(FILLS)))
    ax.bar(x, TOP, label="顶上热沉", color="#b35806")
    ax.bar(x, BOTTOM, bottom=TOP, label="底面", color="#f1a340")
    ax.bar(x, EDGE, bottom=[a + b for a, b in zip(TOP, BOTTOM)],
           label="侧边", color="#542788")
    for i in x:
        ax.annotate(f"{TOP[i]:.1f}", (i, TOP[i] / 2), ha="center", va="center",
                    color="white", fontsize=9)
        ax.annotate(f"{EDGE[i]:.1f}", (i, TOP[i] + BOTTOM[i] + EDGE[i] / 2),
                    ha="center", va="center", color="white", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(FILLS)
    ax.set_ylabel("占总功耗的百分比 (%)")
    ax.set_xlabel("缝的填充，导热系数 W/m·K")
    ax.set_ylim(0, 100)
    ax.set_title("154.79 W，三条路径的占比")
    ax.legend(loc="lower left")
    fig.tight_layout()
    out = HERE / "图" / "散热路径占比百分比.png"
    fig.savefig(out, dpi=140)
    print(f"写到 {out}")


if __name__ == "__main__":
    main()
