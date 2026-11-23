#!/usr/bin/env python3
"""填充材料对三条散热路径占比的影响。数字来自重算，和存档逐位一致。"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pathlib import Path

font_manager.fontManager.addfont(
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

# 154.79 W，侧边接到 45 °C，底面先过 1 mm、k=0.5 的有机层再对流。
# 重算结果，和 热仿真流水线/runs/interposer_map/fill_with_both_exits.json 逐位一致。
FILLS = ["空气\n0.028", "有机\n0.5", "环氧\n0.8", "底部填充\n1.0", "硅\n130"]
TOP = [138.7173, 119.3590, 113.8646, 110.9103, 47.9294]
BOTTOM = [11.3525, 16.0748, 15.9997, 15.8434, 7.8025]
EDGE = [4.7202, 19.3563, 24.9258, 28.0364, 99.0581]

HERE = Path(__file__).resolve().parent


def main() -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    x = list(range(len(FILLS)))
    w = 0.26
    b1 = ax.bar([i - w for i in x], TOP, w, label="顶上热沉", color="#b35806")
    b2 = ax.bar(x, BOTTOM, w, label="底面", color="#f1a340")
    b3 = ax.bar([i + w for i in x], EDGE, w, label="侧边", color="#542788")
    ax.set_xticks(x)
    ax.set_xticklabels(FILLS)
    ax.set_ylabel("带走的热 (W)")
    ax.set_xlabel("缝的填充，导热系数 W/m·K")
    ax.set_title("154.79 W，三条路径各带走多少")
    ax.legend()
    for bars in (b1, b2, b3):
        for bar in bars:
            v = bar.get_height()
            ax.annotate(f"{v:.1f}",
                        xy=(bar.get_x() + bar.get_width() / 2, v),
                        xytext=(0, 3), textcoords="offset points",
                        ha="center", fontsize=8)
    ax.set_ylim(0, 165)
    fig.tight_layout()
    out = HERE / "图" / "散热路径占比.png"
    fig.savefig(out, dpi=140)
    print(f"写到 {out}")


if __name__ == "__main__":
    main()
