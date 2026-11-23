#!/usr/bin/env python3
"""同一趟任务，换三块真实芯片，把性能、功耗、温度、应力走完。

芯片只用 McPAT 官方模板，不发明 P(面积)：
- Penryn  45 nm  2 核  3.7 GHz  91.13 mm²
- T2      65 nm  8 核  1.4 GHz  273.52 mm²   （Niagara2）
- Xeon    65 nm  2 核  3.4 GHz  421.68 mm²

任务运算量用已有的 x86 192×192 matmul，N_op = 11_675_900，不随温度和 s 变。
IPC 各芯片相同、且不随 s 变，这是现有频率闭环的同一条假设。核数能不能用上，
matmul 本身没有在这三块芯片上重跑过，所以吞吐给两种，都写在结果里：

- serial：Θ ∝ f × s。任务只占一个核。
- parallel：Θ ∝ n_cores × f × s。任务铺满所有核。

功耗沿用闭环：P = s × 动态 + 漏电(温度档)。漏电不乘 s。动态和 380 K 官方表
的 sub/gate 来自 OFFICIAL_380，亚阈按 Ioff 表相对 80°C 缩放。
温度是 HotSpot 64 网格、72 mm 封装、芯片居中的方块。10 K 一档，从 350 K
迭代到档位不动。T < 90 °C 且芯片内 von Mises < 200 MPa 才算能用。
应力阶段用 ATSim，填充仍是硅，和现有代理一致。
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

from geometry import Rect  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402

HS = ROOT / "HotSpot-master" / "hotspot"
OUT = ROOT / "热仿真流水线" / "runs" / "task_chip_coupling"
N_OP = 11675900.0
TASK = "x86 192x192 matmul"
S_LIST = (0.5, 1.0, 1.5, 2.0)
ROWS = 64
SIDE = 0.072
T_LIM_C = 90.0
SIG_LIM_MPA = 200.0
MAX_INNER = 5
T_START_K = 350
T_MIN_K, T_MAX_K = 300, 400

# 与 run_full_OFFICIAL_380 / IOFF_HP 同一张表。热学阶段不导入那个模块，
# 因为它在 import 时会拉起 fenics。
OFFICIAL_380 = {
    "Penryn": {"area_mm2": 91.1348, "dyn_W": 47.9315, "sub_W": 10.7217, "gate_W": 0.588125},
    "T2": {"area_mm2": 273.523, "dyn_W": 50.421, "sub_W": 15.0682, "gate_W": 3.50762},
    "Xeon": {"area_mm2": 421.676, "dyn_W": 66.8485, "sub_W": 35.4132, "gate_W": 1.68429},
}
IOFF_HP = {
    "Penryn": {
        0: 2.8e-7, 10: 3.28e-7, 20: 3.81e-7, 30: 4.39e-7, 40: 5.02e-7, 50: 5.69e-7,
        60: 6.42e-7, 70: 7.2e-7, 80: 8.03e-7, 90: 8.91e-7, 100: 9.84e-7,
    },
    "T2": {
        0: 1.96e-7, 10: 2.29e-7, 20: 2.66e-7, 30: 3.05e-7, 40: 3.49e-7, 50: 3.95e-7,
        60: 4.45e-7, 70: 4.97e-7, 80: 5.48e-7, 90: 5.94e-7, 100: 6.3e-7,
    },
}
IOFF_HP["Xeon"] = IOFF_HP["T2"]


def leak_of(model: str, T_K: int) -> dict:
    off = OFFICIAL_380[model]
    tab = IOFF_HP[model]
    if (T_K - 300) not in tab:
        raise ValueError(f"{model} T={T_K} off the 10 K grid")
    ratio = tab[T_K - 300] / tab[80]
    sub = off["sub_W"] * ratio
    return {
        "model": model, "T_K": T_K, "dyn_W": off["dyn_W"], "sub_W": sub,
        "gate_W": off["gate_W"], "leak_W": sub + off["gate_W"], "Ioff_ratio": ratio,
        "area_mm2": off["area_mm2"],
    }


def round10_k(T_C: float) -> int:
    return max(T_MIN_K, min(T_MAX_K, int(round((T_C + 273.15) / 10.0) * 10)))

# 核数和主频来自官方 xml / T2 结果头，不是面积推出来的。
CHIPS = (
    {"model": "Penryn", "cores": 2, "f_MHz": 3700},
    {"model": "T2", "cores": 8, "f_MHz": 1400},
    {"model": "Xeon", "cores": 2, "f_MHz": 3400},
)


def die_rect(model: str, power_W: float) -> Rect:
    area = OFFICIAL_380[model]["area_mm2"]
    side = math.sqrt(area) * 1e-3
    x = (SIDE - side) / 2.0
    return Rect(model, x, x, side, side, power_W, area, "CPU")


def hotspot(model: str, power_W: float, work: Path) -> dict:
    write_hotspot_inputs([die_rect(model, power_W)], work, rows=ROWS, cols=ROWS)
    out = run_hotspot(HS, work)
    temps = out["temperatures_K"]
    # 锚点是零功耗的定位块，峰值只取芯片自己。
    t_k = temps[model]
    return {"T_K": t_k, "T_C": t_k - 273.15, "block_T_C": {k: v - 273.15 for k, v in temps.items()}}


def thermal_point(model: str, s: float, work: Path) -> dict:
    T = T_START_K
    trace = []
    solved = None
    for _ in range(MAX_INNER):
        m = leak_of(model, T)
        power = s * m["dyn_W"] + m["leak_W"]
        solved = hotspot(model, power, work / f"bin_{T}")
        nxt = round10_k(solved["T_C"])
        trace.append({"bin_K": T, "power_W": power, "dyn_W": s * m["dyn_W"],
                      "leak_W": m["leak_W"], "T_C": solved["T_C"], "next_bin_K": nxt})
        if nxt == T:
            break
        T = nxt
    m = leak_of(model, T)
    chip = next(c for c in CHIPS if c["model"] == model)
    serial = chip["f_MHz"] * s
    parallel = chip["cores"] * serial
    # 以 Penryn、s=1、单线程为 1。IPC 消掉，只剩频率和核数。
    ref = 3700.0
    return {
        "model": model,
        "s": s,
        "cores": chip["cores"],
        "f_MHz": chip["f_MHz"],
        "area_mm2": m["area_mm2"],
        "side_mm": math.sqrt(m["area_mm2"]),
        "McPAT_T_K": T,
        "converged": trace[-1]["next_bin_K"] == trace[-1]["bin_K"],
        "power_W": trace[-1]["power_W"],
        "dyn_W": trace[-1]["dyn_W"],
        "leak_W": trace[-1]["leak_W"],
        "leak_frac": trace[-1]["leak_W"] / trace[-1]["power_W"],
        "Tmax_C": trace[-1]["T_C"],
        "fail_T": trace[-1]["T_C"] >= T_LIM_C,
        "Theta_serial_rel": serial / ref,
        "Theta_parallel_rel": parallel / ref,
        "trace": trace,
    }


def run_thermal() -> list[dict]:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for chip in CHIPS:
        for s in S_LIST:
            rec = thermal_point(chip["model"], s, OUT / chip["model"] / f"s_{s:.2f}")
            rows.append(rec)
            (OUT / chip["model"] / f"s_{s:.2f}" / "thermal.json").write_text(
                json.dumps(rec, indent=2) + "\n", encoding="utf-8")
            print(
                f"{rec['model']:7s} s={s:.2f}  P={rec['power_W']:.2f} W "
                f"(dyn {rec['dyn_W']:.2f} leak {rec['leak_W']:.2f})  "
                f"T={rec['Tmax_C']:.2f} °C bin {rec['McPAT_T_K']}  "
                f"{'FAIL_T' if rec['fail_T'] else 'T_ok'}  "
                f"Θser {rec['Theta_serial_rel']:.3f}  Θpar {rec['Theta_parallel_rel']:.3f}",
                flush=True,
            )
    (OUT / "thermal.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    return rows


def run_stress(rows: list[dict]) -> None:
    """必须用 fenicsx 解释器。每个点算完就落盘。"""
    from types import SimpleNamespace

    from atsim_sa_proxy import eval_atsim_proxy, write_proxy_json

    for rec in rows:
        work = OUT / rec["model"] / f"s_{rec['s']:.2f}"
        chip = die_rect(rec["model"], rec["power_W"])
        proxy = SimpleNamespace(
            name=chip.name, x=chip.x, y=chip.y, w=chip.w, h=chip.h, power_W=chip.power_W)
        jp = work / "case.json"
        write_proxy_json([proxy], jp, f"{rec['model']} s={rec['s']}")
        st = eval_atsim_proxy(jp, verbose=False)
        keep = {k: st.get(k) for k in (
            "Tmax_C", "vm_normal_chip_max_MPa", "vm_normal_chip_name",
            "vm_normal_chip_p95_MPa", "w_pv_um", "elapsed_s")}
        rec["stress"] = keep
        rec["sigma_MPa"] = keep["vm_normal_chip_max_MPa"]
        rec["fail_sigma"] = (rec["sigma_MPa"] or 0) >= SIG_LIM_MPA
        rec["feasible"] = (not rec["fail_T"]) and (not rec["fail_sigma"])
        (work / "stress.json").write_text(
            json.dumps(keep, indent=2) + "\n", encoding="utf-8")
        (work / "thermal.json").write_text(
            json.dumps(rec, indent=2) + "\n", encoding="utf-8")
        print(
            f"stress {rec['model']:7s} s={rec['s']:.2f}  "
            f"σ={rec['sigma_MPa']:.2f} MPa  "
            f"{'FAIL_SIG' if rec['fail_sigma'] else 'SIG_ok'}  "
            f"{keep['elapsed_s']:.0f} s",
            flush=True,
        )
    (OUT / "thermal.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    summarize(rows)


def summarize(rows: list[dict]) -> None:
    feasible = [r for r in rows if r.get("feasible")]

    def best(key: str):
        if not feasible:
            return None
        return max(feasible, key=lambda r: r[key])

    report = {
        "task": TASK,
        "N_op": N_OP,
        "T_lim_C": T_LIM_C,
        "sigma_lim_MPa": SIG_LIM_MPA,
        "assumptions": [
            "P = s * McPAT_dyn + leak(T_bin); leak 不乘 s",
            "IPC 各芯片相同且不随 s 变",
            "serial 吞吐 ∝ f*s；parallel 吞吐 ∝ n_cores*f*s",
            "方芯片居中，72 mm 封装，HotSpot 64 网格",
            "ATSim 填充为硅，应力取芯片内部 von Mises",
        ],
        "winner_serial": best("Theta_serial_rel"),
        "winner_parallel": best("Theta_parallel_rel"),
        "cases": rows,
    }
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"写出 {OUT / 'report.json'}", flush=True)


def main() -> None:
    phase = sys.argv[1] if len(sys.argv) > 1 else "thermal"
    if phase == "thermal":
        run_thermal()
    elif phase == "stress":
        rows = json.loads((OUT / "thermal.json").read_text())
        run_stress(rows)
    else:
        raise SystemExit(f"unknown phase {phase}")


if __name__ == "__main__":
    main()
