"""ascend910 上大芯片和小芯片的老化涨幅。

左图是老化涨幅。大芯片和 14 W 那颗各一条线，四颗 20 W 的小芯片画成色带，
带的上下沿是四颗里涨得最多和最少的。实线只有主路径，虚线加了次路径。

右图是大芯片比最冷小芯片热多少，四种情形各一条线。
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

ROOT = Path("/home/czh/任务仿真产生功耗文件工具/论文证据/填充场叠加")
DATA = json.loads((ROOT / "复现" / "ascend910_涨幅.json").read_text())

ORDER = ["air", "organic", "epoxy", "underfill", "silicon"]
XLABEL = ["空气", "有机", "环氧", "底填", "硅"]


def row(secondary, tim, fill):
    return next(r for r in DATA["rows"]
                if r["secondary"] == secondary and r["tim"] == tim and r["fill"] == fill)


def rise(secondary, fill, die):
    a = row(secondary, "x3", fill)["dies_C"][str(die)]
    b = row(secondary, "x1", fill)["dies_C"][str(die)]
    return a - b


def band(secondary, dies):
    """每颗芯粒一条涨幅，返回六种填充下的最低和最高。"""
    lo, hi = [], []
    for fill in ORDER:
        values = [rise(secondary, fill, d) for d in dies]
        lo.append(min(values))
        hi.append(max(values))
    return lo, hi


def gap(secondary, tim, fill):
    d = row(secondary, tim, fill)["dies_C"]
    return d["0"] - min(d[str(i)] for i in range(1, 6))


def main():
    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 5.6))
    xs = list(range(len(ORDER)))

    for sec, color, ls, tag in ((False, "#1f4e79", "-", "只有主路径"),
                                (True, "#2f6f4e", "--", "加次路径")):
        lo, hi = band(sec, [2, 3, 4, 5])
        left.fill_between(xs, lo, hi, color=color, alpha=0.15)
        left.plot(xs, lo, color=color, lw=0.6, ls=ls)
        left.plot(xs, hi, color=color, lw=0.6, ls=ls,
                  label=f"四颗 20 W 的范围，{tag}")

    left.plot(xs, [rise(False, f, 0) for f in ORDER], marker="o", color="#b85c38",
              lw=1.6, label="大芯片 256 W，只有主路径")
    left.plot(xs, [rise(True, f, 0) for f in ORDER], marker="s", color="#b85c38",
              lw=1.6, ls="--", label="大芯片 256 W，加次路径")
    left.plot(xs, [rise(False, f, 1) for f in ORDER], marker="o", color="#7a4e9a",
              lw=1.4, label="小芯片 14 W，只有主路径")
    left.plot(xs, [rise(True, f, 1) for f in ORDER], marker="s", color="#7a4e9a",
              lw=1.4, ls="--", label="小芯片 14 W，加次路径")
    left.set_ylabel("老化涨幅 °C")
    left.set_title("老化带来的温度涨幅")

    right.plot(xs, [gap(False, "x1", f) for f in ORDER], marker="o", color="#8a8a8a",
               lw=1.4, label="只有主路径，新")
    right.plot(xs, [gap(False, "x3", f) for f in ORDER], marker="o", color="#b85c38",
               lw=1.4, label="只有主路径，老化")
    right.plot(xs, [gap(True, "x1", f) for f in ORDER], marker="s", color="#1f4e79",
               lw=1.4, ls="--", label="加次路径，新")
    right.plot(xs, [gap(True, "x3", f) for f in ORDER], marker="s", color="#2f6f4e",
               lw=1.4, ls="--", label="加次路径，老化")
    right.set_ylabel("大芯片比最冷小芯片热 °C")
    right.set_title("最热和最冷的温差")

    for ax in (left, right):
        ax.set_xticks(xs)
        ax.set_xticklabels(XLABEL)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("ascend910，350 W。大芯片 256 W，功率密度 0.56 W/mm²；"
                 "14 W 那颗是 0.08 W/mm²。", fontsize=12)
    fig.tight_layout()
    fig.savefig(ROOT / "图" / "ascend910_涨幅.png", dpi=140)
    plt.close(fig)
    print("写到 图/ascend910_涨幅.png")


if __name__ == "__main__":
    main()
