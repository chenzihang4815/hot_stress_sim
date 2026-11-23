#!/usr/bin/env python3
"""三张热图：初始布局、TAP 原版 120 步、16 网格加场叠加 120 步。

温度场都是 64 网格、芯粒层（Layer 4），同一色标，单位 °C。
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).parent
ROWS = 64

def grid(path):
    vals = [float(v) for v in Path(path).read_text().split()]
    assert len(vals) == ROWS * ROWS, (path, len(vals))
    # HotSpot 按行从下往上写，imshow 原点在上方，翻转后 y 轴向上。
    return (np.array(vals).reshape(ROWS, ROWS) - 273.15)[::-1]


panels = [
    ("初始布局", HERE / "init_64.grid", HERE / "init.json"),
    ("TAP 原版，120 步", HERE / "replay64" / "best_full_64.grid",
     HERE / "replay64" / "best_full.json"),
    ("16 网格 + 场叠加，120 步", HERE / "replay16" / "best_field_64.grid",
     HERE / "replay16" / "best_field.json"),
]
fields = [grid(p) for _, p, _ in panels]
lo = min(f.min() for f in fields)
hi = max(f.max() for f in fields)

fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.6))
for ax, (title, _, layout), f in zip(axes, panels, fields):
    im = ax.imshow(f, extent=[0, 50, 0, 50], origin="lower",
                   vmin=lo, vmax=hi, cmap="inferno")
    for c in json.loads(Path(layout).read_text()):
        ax.add_patch(plt.Rectangle((c["x_mm"] - c["w_mm"] / 2, c["y_mm"] - c["h_mm"] / 2),
                                   c["w_mm"], c["h_mm"], fill=False,
                                   ec="white", lw=0.7))
    ax.set_title(f"{title}\n峰值 {f.max():.2f} °C", fontsize=11)
    ax.set_xlabel("x (mm)")
    ax.set_xlim(0, 50)
    ax.set_ylim(0, 50)
axes[0].set_ylabel("y (mm)")
fig.colorbar(im, ax=axes, fraction=0.046, pad=0.02, label="°C")
fig.suptitle("ascend910 芯粒层温度场，64 网格，同一色标", fontsize=12)
fig.savefig(HERE / "heatmaps.png", dpi=140, bbox_inches="tight")
print(HERE / "heatmaps.png")
