#!/usr/bin/env python3
"""10 步退火，材料图不更新对每步更新的误差。"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

font_manager.fontManager.addfont(
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
ROWS = json.loads((HERE / "复现" / "ten_frozen.json").read_text())["rows"]


def main() -> None:
    steps = [r["step"] for r in ROWS]
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.plot(steps, [r["peak_err_K"] for r in ROWS], "o-", color="#b35806",
            label="峰值误差")
    ax.plot(steps, [r["grid_max_K"] for r in ROWS], "s-", color="#542788",
            label="整图最大误差")
    for r in ROWS:
        if not r["accept"]:
            ax.axvline(r["step"], color="#999999", lw=0.6, ls="--")
    ax.annotate("第 6 步被拒绝", xy=(6, 29.15), xytext=(6.4, 26.2),
                fontsize=9, color="#666666")
    ax.set_xticks(steps)
    ax.set_xlabel("步")
    ax.set_ylabel("frozen − exact (K)")
    ax.set_title("10 步，材料图停在起点，缝环氧 0.8")
    ax.legend()
    ax.set_ylim(0, 36)
    fig.tight_layout()
    out = HERE / "图" / "十步误差.png"
    fig.savefig(out, dpi=140)
    print(f"写到 {out}")


if __name__ == "__main__":
    main()
