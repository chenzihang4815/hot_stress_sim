"""ascend910 上，大芯片和小芯片的老化涨幅。

用论文那一轮的初始布局，六颗芯粒。芯粒 0 是 14.5×31.4 mm、256 W 的大芯片，
其余五颗是 14 到 20 W 的小芯片。总功耗 350 W。

四种情形：只有主路径或加次路径，TIM 新的或老化。五种填充。涨幅是老化减新。

python 场，64 网格，和 tim_aging_paths 同一套模型。中介层按 50 mm 取，
比四芯粒那一组的 96 mm 小，格子是 0.78125 mm。
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
sys.path.insert(0, str(ROOT / "热仿真流水线" / "src"))
sys.path.insert(0, str(ROOT / "论文证据" / "填充场叠加"))

import rerun_heat_paths as hp  # noqa: E402
import run_interposer_map as im  # noqa: E402

SIDE = 0.050
N = 64
CELL = SIDE / N
FILLS = [("air", 0.028), ("organic", 0.5), ("epoxy", 0.8),
         ("underfill", 1.0), ("silicon", 130.0)]
TIMS = (("x1", 4.0), ("x3", 4.0 / 3.0))
OUT = ROOT / "论文证据" / "填充场叠加" / "复现" / "ascend910_涨幅.json"


def grid_of(layout):
    """芯粒中心是 mm，宽高也是 mm。返回每格的填充导热和功耗。"""
    k = np.full((N, N), np.nan)
    power = np.zeros((N, N))
    for c in layout:
        x0 = (c["x_mm"] - c["w_mm"] / 2) / 1000
        y0 = (c["y_mm"] - c["h_mm"] / 2) / 1000
        x1 = x0 + c["w_mm"] / 1000
        y1 = y0 + c["h_mm"] / 1000
        cells = []
        for i in range(N):
            for j in range(N):
                cx = (j + 0.5) * CELL
                cy = (N - 1 - i + 0.5) * CELL
                if x0 <= cx < x1 and y0 <= cy < y1:
                    cells.append((i, j))
        if not cells:
            raise SystemExit(f"芯粒 {c['i']} 没有落到格子上")
        each = c["power_W"] / len(cells)
        for i, j in cells:
            k[i, j] = 130.0
            power[i, j] = each
        c["n_cells"] = len(cells)
    return k, power


def solve(k_chip, power, k_tim, secondary):
    if secondary:
        st = im.solve_with_interposer(k_chip, power, k_tim, g_bottom=hp.G_ORG,
                                       edge_to_ambient=True)
    else:
        st = im.solve_with_interposer(k_chip, power, k_tim, g_bottom=0.0,
                                       edge_to_ambient=False)
    return st[1]


def die_max(field, layout, chip):
    x0 = (chip["x_mm"] - chip["w_mm"] / 2) / 1000
    y0 = (chip["y_mm"] - chip["h_mm"] / 2) / 1000
    best = None
    for i in range(N):
        for j in range(N):
            cx = (j + 0.5) * CELL
            cy = (N - 1 - i + 0.5) * CELL
            if x0 <= cx < x0 + chip["w_mm"] / 1000 and y0 <= cy < y0 + chip["h_mm"] / 1000:
                v = float(field[i, j]) - 273.15
                best = v if best is None else max(best, v)
    return best


def main():
    layout = json.loads((ROOT / "论文证据" / "tap场叠加" / "runs"
                         / "ascend910_16_full" / "init.json").read_text())
    k_base, power = grid_of(layout)
    total = float(power.sum())
    print(f"总功耗 {total:.1f} W  格子 {CELL * 1000:.3f} mm", flush=True)
    for c in layout:
        print(f"  芯粒 {c['i']}  {c['w_mm']}×{c['h_mm']} mm  {c['power_W']:.0f} W"
              f"  {c['n_cells']} 格", flush=True)

    rows = []
    for secondary in (False, True):
        for tag, k_tim in TIMS:
            for fill, k_fill in FILLS:
                k_chip = np.where(np.isnan(k_base), k_fill, k_base)
                field = solve(k_chip, power, k_tim, secondary)
                dies = {c["i"]: round(die_max(field, layout, c), 3) for c in layout}
                rows.append({"fill": fill, "k_fill": k_fill, "tim": tag,
                             "secondary": secondary, "dies_C": dies,
                             "peak_C": round(float(field.max()) - 273.15, 3)})
                print(f"  {'次开' if secondary else '次关'} {tag} {fill:10s}"
                      f"  大芯片 {dies[0]:7.2f}  小芯片最高 {max(dies[i] for i in range(1, 6)):7.2f}",
                      flush=True)

    out = {
        "说明": "ascend910 初始布局，python 场，64 网格，中介层 50 mm。"
               "芯粒 0 是 256 W 的大芯片，其余五颗是小芯片。涨幅是老化减新。",
        "total_power_W": total,
        "cell_mm": CELL * 1000,
        "layout": layout,
        "rows": rows,
    }
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"写到 {OUT}", flush=True)


if __name__ == "__main__":
    main()
