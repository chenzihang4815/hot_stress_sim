#!/usr/bin/env python3
"""46.3 倍对照图：TAP 原版（1 mm、每步 64 网格全解）对 16 网格场叠加。"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).parent
base = json.loads((HERE / "ascend910_64.json").read_text())
fast = json.loads((HERE / "ascend910_1mm.json").read_text())

fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))

ax = axes[0]
vals = [base["full_secs"] / 60, fast["field_secs"] / 60]
bars = ax.bar(["TAP 原版\n1 mm，每步 64 网格全解", "16 网格 + 场叠加"],
              vals, color=["#b0b0b0", "#3b6ea5"], width=0.55)
for b, v, s in zip(bars, vals, [base["full_secs"], fast["field_secs"]]):
    ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.1f} 分\n({s:.0f} s)",
            ha="center", va="bottom", fontsize=9)
ax.set_ylabel("求解用时")
ax.set_ylim(0, max(vals) * 1.3)
ax.set_title(f"用时  {base['full_secs'] / fast['field_secs']:.1f} 倍")

ax = axes[1]
temps = [base["full_final_64_C"], fast["field_final_64_C"]]
bars = ax.bar(["TAP 原版", "16 网格 + 场叠加"], temps,
              color=["#b0b0b0", "#3b6ea5"], width=0.55)
for b, t in zip(bars, temps):
    ax.text(b.get_x() + b.get_width() / 2, t, f"{t:.2f} °C",
            ha="center", va="bottom", fontsize=9)
ax.set_ylim(min(temps) - 1.0, max(temps) + 1.6)
ax.set_ylabel("最优布局的峰值温度（64 网格复核）")
ax.set_title(f"温差  {fast['field_final_64_C'] - base['full_final_64_C']:+.2f} °C")

fig.suptitle("ascend910，120 步，种子 22，1 mm 步长，线长项已去掉", fontsize=12)
fig.tight_layout()
fig.savefig(HERE / "tap46x.png", dpi=140)
print(HERE / "tap46x.png")
