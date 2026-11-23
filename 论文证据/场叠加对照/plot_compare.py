#!/usr/bin/env python3
"""场叠加对照的两张图：时间，和 64 网格重评的峰值温度。

数字全部读自同目录的 table.json，不在这里重算。
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

HERE = Path(__file__).parent
table = json.loads((HERE / "table.json").read_text(encoding="utf-8"))
cases = table["cases"]
labels = [c["case"] for c in cases]
x = range(len(cases))
w = 0.36

FULL = "#b0b0b0"
FIELD = "#3b6ea5"

fig, ax = plt.subplots(figsize=(7.2, 4.4))
b1 = ax.bar([i - w / 2 for i in x], [c["full_solve_s"] for c in cases],
            w, color=FULL, label="全解")
b2 = ax.bar([i + w / 2 for i in x], [c["field_solve_s"] for c in cases],
            w, color=FIELD, label="场叠加")
for b in list(b1) + list(b2):
    v = b.get_height()
    ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.0f} s",
            ha="center", va="bottom", fontsize=8)
for i, c in enumerate(cases):
    ax.text(i, max(c["full_solve_s"], c["field_solve_s"]) * 1.09,
            f"{c['speedup']:.2f} 倍", ha="center", va="bottom", fontsize=9)
ax.set_xticks(list(x), labels)
ax.set_ylabel("求解用时 (s)")
ax.set_ylim(0, max(c["full_solve_s"] for c in cases) * 1.28)
ax.legend(frameon=False)
ax.set_title("1050 步，16 网格搜索，只按温度接受")
fig.tight_layout()
fig.savefig(HERE / "时间优化.png", dpi=160)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7.2, 4.4))
b1 = ax.bar([i - w / 2 for i in x], [c["best_full"]["Tmax_C"] for c in cases],
            w, color=FULL, label="全解最优")
b2 = ax.bar([i + w / 2 for i in x], [c["best_field"]["Tmax_C"] for c in cases],
            w, color=FIELD, label="场叠加最优")
base = 100.0
for b in list(b1) + list(b2):
    v = b.get_height()
    ax.text(b.get_x() + b.get_width() / 2, v + 0.6, f"{v:.2f}",
            ha="center", va="bottom", fontsize=8)
for i, c in enumerate(cases):
    top = max(c["best_full"]["Tmax_C"], c["best_field"]["Tmax_C"])
    ax.text(i, top + 4.2, f"{c['Tmax_error_K']:+.2f} °C",
            ha="center", va="bottom", fontsize=9)
hi = max(c["best_field"]["Tmax_C"] for c in cases)
ax.set_xticks(list(x), labels)
ax.set_ylabel("峰值温度 (°C)")
ax.set_ylim(base, hi + 12)
ax.legend(frameon=False, loc="upper left")
ax.set_title("最优布局的峰值温度，HotSpot 64 网格重评，纵轴从 100 °C 起")
fig.tight_layout()
fig.savefig(HERE / "峰值优化.png", dpi=160)
plt.close(fig)
print("written")
