#!/usr/bin/env python3
"""120 步退火，场叠加的最优减全解的最优。正值是场叠加更热。"""

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
RUNS = HERE / "runs"
CASES = ("ascend910", "micro150", "multigpu")


def delta(case: str, rows: str) -> float:
    full = json.loads((RUNS / f"{case}_{rows}_full" / "summary.json").read_text())
    field = json.loads((RUNS / f"{case}_{rows}_field" / "summary.json").read_text())
    return field["final_64_C"] - full["final_64_C"]


def broken(ax_top, ax_bot) -> None:
    """纵轴左侧一对角标，标出断开的位置。"""
    ax_top.spines["bottom"].set_visible(False)
    ax_bot.spines["top"].set_visible(False)
    kwargs = dict(marker=[(-1, -0.5), (1, 0.5)], markersize=9,
                  linestyle="none", color="black", mec="black", mew=1,
                  clip_on=False)
    ax_top.plot([0], [0], transform=ax_top.transAxes, **kwargs)
    ax_bot.plot([0], [1], transform=ax_bot.transAxes, **kwargs)


def main() -> None:
    # micro150 的 64 网格差 2.49 °C，其余都在 0.04 以内。纵轴从中断开，
    # 否则小的五根贴在零线上看不见。
    fig, (top, bot) = plt.subplots(
        2, 1, sharex=True, figsize=(7.2, 5.4),
        gridspec_kw={"height_ratios": (1, 2.6), "hspace": 0.05})
    x = list(range(len(CASES)))
    width = 0.36
    d16 = [delta(c, "16") for c in CASES]
    d64 = [delta(c, "64") for c in CASES]
    for ax in (top, bot):
        ax.bar([i - width / 2 for i in x], d16, width, label="16 网格搜索",
               color="#b35806")
        ax.bar([i + width / 2 for i in x], d64, width, label="64 网格搜索",
               color="#542788")
        ax.axhline(0, color="black", lw=0.6)
    broken(top, bot)
    bot.set_ylim(-0.015, 0.075)
    top.set_ylim(2.35, 2.66)
    top.set_yticks([2.4, 2.5, 2.6])
    bot.set_xticks(x)
    bot.set_xticklabels(CASES)
    bot.set_ylabel("场叠加最优 − 全解最优 (°C)")
    bot.yaxis.set_label_coords(-0.09, 0.85)
    top.set_title("TAP 120 步，场叠加对全解的最优温度差")
    top.legend(loc="center left")
    for ax, lo, hi in ((bot, -0.015, 0.075), (top, 2.35, 2.66)):
        for bar in ax.patches:
            v = bar.get_height()
            if not lo <= v <= hi:
                continue
            ax.annotate(f"{v:+.2f}",
                        xy=(bar.get_x() + bar.get_width() / 2, v),
                        xytext=(0, 3), textcoords="offset points",
                        ha="center", fontsize=9)
    out = HERE / "图" / "误差柱状图.png"
    fig.savefig(out, dpi=140)
    print(f"写到 {out}")


if __name__ == "__main__":
    main()
