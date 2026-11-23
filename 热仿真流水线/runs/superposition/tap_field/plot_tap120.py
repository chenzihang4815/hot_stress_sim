#!/usr/bin/env python3
"""120 步对照图。数据来自同一次运行的 ascend910.json。"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"

HERE = Path(__file__).parent
rec = json.loads((HERE / "ascend910.json").read_text())

fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))

ax = axes[0]
bars = ax.bar(["每步全解", "场叠加"], [rec["full_secs"], rec["field_secs"]],
              color=["#b0b0b0", "#3b6ea5"], width=0.55)
ax.bar_label(bars, fmt="%.1f s", padding=3)
ax.set_ylabel("求解用时 (秒)")
ax.set_ylim(0, rec["full_secs"] * 1.25)
ax.set_title(f"用时  {rec['speedup']:.2f} 倍")

ax = axes[1]
bars = ax.bar(["每步全解", "场叠加"], [rec["full_final_64_C"], rec["field_final_64_C"]],
              color=["#b0b0b0", "#3b6ea5"], width=0.55)
ax.bar_label(bars, fmt="%.2f °C", padding=3)
lo = min(rec["full_final_64_C"], rec["field_final_64_C"])
ax.set_ylim(lo - 1.0, lo + 2.2)
ax.set_ylabel("64 网格复核的峰值温度 (°C)")
ax.set_title(f"温差  {rec['delta_C']:+.2f} °C")

fig.suptitle("ascend910，120 步，16 网格搜索，同一条退火轨迹", fontsize=12)
fig.tight_layout()
fig.savefig(HERE / "tap120.png", dpi=140)
print(HERE / "tap120.png")
