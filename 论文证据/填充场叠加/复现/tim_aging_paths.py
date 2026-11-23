"""TIM 老化会不会改变主次散热路径的占比。

和 rerun_heat_paths.py 同一套模型、同一批芯片、同一组出口。五层：中介层
100 µm 硅、芯片层 150 µm、TIM 20 µm、铜均热片 1 mm、热沉 6.9 mm。芯片是
base_dies 的四颗，Penryn、GPU_3070、HBM3、CIM_jia。缝的填充走五种，芯片占的
格子固定是硅 130。

出口和 fill_with_both_exits 相同：底面先过 1 mm、k=0.5 的有机层再接到环境，
中介层最外一圈用 1e4 W/K 接到 45 °C。这一组才有次路径，封死底面的那组回答不了
占比的问题。

老化只改 TIM 的导热，厚度不变。新的是 k=4，老化是 k=4/3，热阻变成 3 倍。
热流从求解后的温度场积分，不另解。三条出口加起来等于总功耗。

python 场，不是 HotSpot，也不是 ATSim。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
sys.path.insert(0, str(ROOT / "热仿真流水线" / "src"))
sys.path.insert(0, str(ROOT / "论文证据" / "填充场叠加"))

import rerun_heat_paths as hp  # noqa: E402
import run_factor_controls as fc  # noqa: E402
import run_hw_area_field as hw  # noqa: E402
import run_interposer_map as im  # noqa: E402

FILLS = [("air", 0.028), ("organic", 0.5), ("epoxy", 0.8),
         ("underfill", 1.0), ("silicon", 130.0)]
TIMS = (("x1", 4.0), ("x3", 4.0 / 3.0))


def one(dies, power, silicon, fill, k_fill, tag, k_tim):
    k_chip = hw.paint(dies, k_fill)
    st = im.solve_with_interposer(k_chip, power, k_tim, g_bottom=hp.G_ORG,
                                  edge_to_ambient=True)
    rec = {
        "fill": fill, "k_fill": k_fill, "tim": tag, "k_tim": k_tim,
        "python_Tmax_C": float(st[1].max() - 273.15),
        "dies_C": {d["name"]: im.die_max(st[1], d) for d in dies},
    }
    rec.update(hp.flows(st, k_chip, hp.G_ORG, True, silicon))
    return rec


def main():
    dies, _ = fc.base_dies()
    power = hw.power_map(dies)
    silicon = np.zeros((im.N, im.N), dtype=bool)
    for d in dies:
        for i, j in hw.cells(d):
            silicon[i, j] = True
    total = float(power.sum())
    print(f"总功耗 {total:.6f} W", flush=True)

    rows = []
    for tag, k_tim in TIMS:
        for fill, k_fill in FILLS:
            rec = one(dies, power, silicon, fill, k_fill, tag, k_tim)
            gap = abs(rec["Q_sum_W"] - total)
            share_top = 100.0 * rec["Q_top_W"] / total
            share_bot = 100.0 * rec["Q_bottom_W"] / total
            share_edge = 100.0 * rec["Q_edge_W"] / total
            print(f"  {tag} {fill:10s} 峰值 {rec['python_Tmax_C']:6.2f}  "
                  f"顶上 {share_top:5.2f}%  底面 {share_bot:5.2f}%  "
                  f"侧边 {share_edge:5.2f}%  账差 {gap:.1e} W", flush=True)
            rows.append(rec)

    out = Path(__file__).resolve().parent
    (out / "tim_aging_paths.json").write_text(json.dumps({
        "note": "底面过 1 mm 有机层接到环境，侧边接到 45 °C，和 fill_with_both_exits "
                "同一组出口。老化是 TIM 导热从 4 降到 4/3，厚度不变。四芯粒，"
                "python 场，不是 HotSpot，也不是 ATSim。",
        "total_power_W": total,
        "rows": rows,
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"写到 {out / 'tim_aging_paths.json'}", flush=True)


if __name__ == "__main__":
    main()
