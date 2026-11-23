"""ascend910 上两组相减之前的温度，拆成四张图。

两张是老化涨幅的两组：新的和老化的并排，六颗芯粒各自一根。
柱上是实际温度，右列是老化减新。
两张是最热减最冷的两组：大芯片和 14 W 并排，对上方是两者的差。
柱从环境温度 45 °C 量起。空心是新的，实心是老化的。
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch

font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

ROOT = Path("/home/czh/任务仿真产生功耗文件工具/论文证据/填充场叠加")
OUT = ROOT / "图"
DATA = json.loads((ROOT / "复现" / "ascend910_涨幅.json").read_text())

ORDER = ["air", "organic", "epoxy", "underfill", "silicon"]
XLABEL = ["空气", "有机", "环氧", "底填", "硅"]
DIES = [0, 2, 3, 4, 5, 1]
DIE_COLOR = {0: "#b85c38", 2: "#1f4e79", 3: "#4a90a4",
             4: "#2f6f4e", 5: "#c4a35a", 1: "#7a4e9a"}
DIE_NAME = {0: "大芯片 256 W", 2: "20 W (38, 22)", 3: "20 W (9, 25)",
            4: "20 W (9, 11)", 5: "20 W (35, 11)", 1: "14 W"}
FLOOR = 45.0


def row(secondary, tim, fill):
    return next(r for r in DATA["rows"]
                if r["secondary"] == secondary and r["tim"] == tim and r["fill"] == fill)


def temp(secondary, tim, fill, die):
    return row(secondary, tim, fill)["dies_C"][str(die)]


def legend_handles():
    handles = [Patch(facecolor=DIE_COLOR[d], edgecolor=DIE_COLOR[d], label=DIE_NAME[d])
               for d in DIES]
    handles.append(Patch(facecolor="white", edgecolor="#555555", label="空心是新的"))
    handles.append(Patch(facecolor="#555555", edgecolor="#555555", label="实心是老化的"))
    return handles


def rise_figure(secondary, title, name):
    fig, ax = plt.subplots(figsize=(13, 16))
    pitch = 0.62
    group = 9.2
    ax.set_xlim(FLOOR, 112)
    yticks, ylabels = [], []
    for fi, fill in enumerate(ORDER):
        for si, state in enumerate(("x1", "x3")):
            center = fi * group + (1 - si) * 4.2
            yticks.append(center)
            ylabels.append(f"{XLABEL[fi]}  {'新' if state == 'x1' else '老化'}")
            for di, die in enumerate(DIES):
                y = center + (2.5 - di) * pitch
                val = temp(secondary, state, fill, die)
                kw = {"height": pitch * 0.78, "linewidth": 1.1,
                      "edgecolor": DIE_COLOR[die]}
                kw["facecolor"] = "white" if state == "x1" else DIE_COLOR[die]
                ax.barh(y, val - FLOOR, left=FLOOR, **kw)
                ax.text(val + 0.5, y, f"{val:.2f}", va="center", ha="left",
                        fontsize=9, color="#222222")
            if state == "x3":
                for di, die in enumerate(DIES):
                    y = center + (2.5 - di) * pitch
                    delta = (temp(secondary, "x3", fill, die)
                             - temp(secondary, "x1", fill, die))
                    ax.text(100.5, y, f"{delta:+.2f}", va="center", ha="left",
                            fontsize=10, color=DIE_COLOR[die], fontweight="medium")
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=11)
    ax.set_xlabel("温度 °C", fontsize=12)
    ax.set_title(title, fontsize=14, pad=12)
    ax.axvline(90, color="#cccccc", lw=0.7, ls=":")
    ax.axvline(98, color="#dddddd", lw=0.6)
    ax.grid(axis="x", ls=":", lw=0.4, color="#dddddd")
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", labelsize=11)
    ax.text(100.5, 1.01, "老化减新", fontsize=11, color="#555555",
            transform=ax.get_xaxis_transform(), clip_on=False)
    fig.legend(handles=legend_handles(), loc="lower center", ncol=4,
               frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.045))
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"写到 图/{name}")


def gap_figure(secondary, title, name):
    fig, ax = plt.subplots(figsize=(13, 8))
    for fi, fill in enumerate(ORDER):
        specs = [(fi - 0.30, "x1", 0), (fi - 0.10, "x1", 1),
                 (fi + 0.10, "x3", 0), (fi + 0.30, "x3", 1)]
        vals = {}
        for x, state, die in specs:
            val = temp(secondary, state, fill, die)
            vals[(state, die)] = val
            kw = {"width": 0.18, "linewidth": 1.2, "edgecolor": DIE_COLOR[die]}
            kw["facecolor"] = "white" if state == "x1" else DIE_COLOR[die]
            ax.bar(x, val - FLOOR, bottom=FLOOR, **kw)
            ax.text(x, val + 0.5, f"{val:.2f}", ha="center", va="bottom",
                    rotation=90, fontsize=10, color="#222222")
        for state, x in (("x1", fi - 0.20), ("x3", fi + 0.20)):
            diff = vals[(state, 0)] - vals[(state, 1)]
            ax.text(x, 96.2, f"差 {diff:.2f}", ha="center", va="bottom",
                    fontsize=11, color="#222222")
    ax.set_xlim(-0.6, 4.6)
    ax.set_ylim(FLOOR, 102)
    ax.set_xticks(range(5))
    ax.set_xticklabels(XLABEL, fontsize=12)
    ax.set_ylabel("温度 °C", fontsize=12)
    ax.set_title(title, fontsize=14, pad=12)
    ax.axhline(90, color="#cccccc", lw=0.7, ls=":")
    ax.grid(axis="y", ls=":", lw=0.4, color="#dddddd")
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=11)
    handles = [Patch(facecolor=DIE_COLOR[0], edgecolor=DIE_COLOR[0], label="大芯片 256 W"),
               Patch(facecolor=DIE_COLOR[1], edgecolor=DIE_COLOR[1], label="小芯片 14 W"),
               Patch(facecolor="white", edgecolor="#555555", label="空心是新的"),
               Patch(facecolor="#555555", edgecolor="#555555", label="实心是老化的")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=11, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"写到 图/{name}")


def main():
    rise_figure(False, "只有主路径：新的和老化的温度，右列是老化减新",
                "ascend910_只有主路径_新旧.png")
    rise_figure(True, "加次路径：新的和老化的温度，右列是老化减新",
                "ascend910_加次路径_新旧.png")
    gap_figure(False, "只有主路径：大芯片减 14 W",
               "ascend910_只有主路径_大小芯片.png")
    gap_figure(True, "加次路径：大芯片减 14 W",
               "ascend910_加次路径_大小芯片.png")


if __name__ == "__main__":
    main()
