#!/usr/bin/env python3
"""120 步，16 网格和 64 网格，场叠加对全解的用时和最优温度差。"""

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
RUNS = HERE / "runs"
CASES = ("ascend910", "micro150", "multigpu")


def load(case: str, rows: str, mode: str) -> dict:
    return json.loads((RUNS / f"{case}_{rows}_{mode}" / "summary.json").read_text())


def main() -> None:
    labels = [f"{c}\n{rows} 网格" for rows in ("16", "64") for c in CASES]
    x = np.arange(len(labels))
    w = 0.36
    full = [load(c, rows, "full") for rows in ("16", "64") for c in CASES]
    field = [load(c, rows, "field") for rows in ("16", "64") for c in CASES]

    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.2))

    ax = axes[0]
    full_s = [r["solve_secs"] for r in full]
    field_s = [r["solve_secs"] for r in field]
    ax.bar(x - w / 2, full_s, width=w, color="#4c72b0", label="每步全解")
    ax.bar(x + w / 2, field_s, width=w, color="#dd8452", label="场叠加")
    ax.axvline(2.5, color="0.75", lw=0.6)
    ax.set_yscale("log")
    for i, (a, b) in enumerate(zip(full_s, field_s)):
        ax.text(i - w / 2, a * 1.12, f"{a:.0f} s", ha="center", fontsize=7.5)
        ax.text(i + w / 2, b * 1.12, f"{b:.0f} s\n{a / b:.2f} 倍",
                ha="center", fontsize=7.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("120 步求解用时，s，对数轴")
    ax.set_ylim(10, max(full_s) * 2.6)
    ax.legend(fontsize=9, loc="upper left")
    ax.set_title("用时", fontsize=12)

    ax = axes[1]
    delta = [b["final_64_C"] - a["final_64_C"] for a, b in zip(full, field)]
    ax.bar(x, delta, width=0.5, color="#542788")
    ax.axhline(0, color="black", lw=0.6)
    ax.axvline(2.5, color="0.75", lw=0.6)
    for i, v in enumerate(delta):
        ax.text(i, v + 0.06, f"{v:+.2f}", ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("场叠加最优 − 全解最优，°C")
    ax.set_ylim(-0.15, 3.0)
    ax.set_title("最优温度差，都是 64 网格重评", fontsize=12)

    fig.suptitle("TAP 120 步，场叠加对每步全解。左三组 16 网格搜索，右三组 64 网格搜索",
                 fontsize=12)
    fig.tight_layout()
    out = HERE / "图" / "用时与误差.png"
    fig.savefig(out, dpi=140)
    print(f"写到 {out}")


if __name__ == "__main__":
    main()
