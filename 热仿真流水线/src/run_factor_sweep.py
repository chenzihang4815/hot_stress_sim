#!/usr/bin/env python3
"""One-factor-at-a-time on the frozen 2+2+2 layout.

Baseline is big_balanced_2 at s=1, 350 K: two T2, two TITANX, two jia64.
Each case changes exactly one thing and re-solves HotSpot and real ATSim on
the same L* coordinates, so the temperature and stress deltas come only from
that one change. Ops come from the gem5 record and move only when the factor
is one that changes how much work a chip does per second.
"""

from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402
from run_full_chain import (  # noqa: E402
    GPU_CONST_W, GPU_STATIC_W, OPS, SIG_LIM_MPA, T_LIM_C, leak_of, real_stress,
)

BASE = ROOT / "热仿真流水线" / "runs" / "big_balanced_2"
OUT = ROOT / "热仿真流水线" / "runs" / "factor_sweep"
HS = ROOT / "HotSpot-master" / "hotspot"

N_OP = OPS["CPU_matmul"]["N_op"]
T_SIM = OPS["CPU_matmul"]["t_sim_s"]


def load_layout() -> list[sa.Chip]:
    chips = []
    for c in json.loads((BASE / "Lstar.json").read_text()):
        chips.append(sa.Chip(c["name"], c["kind"], c["x_m"], c["y_m"], c["w_m"],
                             c["h_m"], c["power_W"], c["area_mm2"]))
    return chips


def power(chips, s_cpu, s_gpu, v_scale, t_k, gated, cores_scale) -> dict:
    """v_scale squares CPU dynamic power. cores_scale multiplies one CPU's
    dynamic power and its ops, standing in for a core count McPAT was never
    re-run for."""
    cpu_dyn = cpu_leak = gpu_dyn = gpu_leak = gpu_const = cim = n_op = 0.0
    for c in chips:
        if c.kind == "CPU":
            m = leak_of("T2", t_k)
            scale = cores_scale if c.name == "CPU_T2_0" else 1.0
            if c.name in gated:
                dyn, leak = 0.0, 0.0
            else:
                dyn, leak = s_cpu * m["dyn_W"] * v_scale * scale, m["leak_W"]
            c.power_W = dyn + leak
            cpu_dyn += dyn
            cpu_leak += leak
            n_op += 0.0 if c.name in gated else N_OP * scale
        elif c.kind == "GPU":
            base = 67.41
            dyn0 = max(0.0, base - GPU_STATIC_W["GPU_TITANX"] - GPU_CONST_W["GPU_TITANX"])
            if c.name in gated:
                c.power_W = 0.0
            else:
                c.power_W = s_gpu * dyn0 + GPU_STATIC_W["GPU_TITANX"] + GPU_CONST_W["GPU_TITANX"]
                gpu_dyn += s_gpu * dyn0
                gpu_leak += GPU_STATIC_W["GPU_TITANX"]
                gpu_const += GPU_CONST_W["GPU_TITANX"]
        else:
            cim += c.power_W
    tot = cpu_dyn + cpu_leak
    return {
        "cpu_dyn_W": cpu_dyn, "cpu_leak_W": cpu_leak,
        "cpu_leak_frac": cpu_leak / tot if tot else None,
        "gpu_dyn_W": gpu_dyn, "pkg_W": tot + gpu_dyn + gpu_leak + gpu_const + cim,
        "N_op": n_op, "ops_per_s": n_op / T_SIM,
    }


def hotspot(chips, tag) -> float:
    d = OUT / tag
    write_hotspot_inputs(sa.to_rects(chips), d)
    return run_hotspot(HS, d)["Tmax_C"]


CASES = [
    ("基线 s=1 350K", dict()),
    ("频率 只降CPU s=0.5", dict(s_cpu=0.5)),
    ("频率 CPU和GPU都 s=0.5", dict(s_cpu=0.5, s_gpu=0.5)),
    ("频率 升到 s=1.5", dict(s_cpu=1.5, s_gpu=1.5)),
    ("电压 降到 0.9", dict(v_scale=0.9 ** 2)),
    ("电压 升到 1.1", dict(v_scale=1.1 ** 2)),
    ("门控 关一颗CPU", dict(gated={"CPU_T2_1"})),
    ("门控 关一颗GPU", dict(gated={"GPU_TITANX_1"})),
    ("核数 一颗CPU翻倍", dict(cores_scale=2.0)),
    ("温度档 360K", dict(t_k=360)),
    ("温度档 340K", dict(t_k=340)),
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    base = load_layout()
    rows = []
    for name, kw in CASES:
        chips = deepcopy(base)
        br = power(chips, kw.get("s_cpu", 1.0), kw.get("s_gpu", 1.0),
                   kw.get("v_scale", 1.0), kw.get("t_k", 350),
                   kw.get("gated", set()), kw.get("cores_scale", 1.0))
        t_hs = hotspot(chips, "hs_" + str(len(rows)))
        st = real_stress(chips, OUT / ("at_" + str(len(rows))), "f")
        ok = st["Tmax_C"] < T_LIM_C and (st["vm_normal_chip_max_MPa"] or 1e9) < SIG_LIM_MPA
        row = {"case": name, **br, "Tmax_hotspot_C": t_hs, "Tmax_atsim_C": st["Tmax_C"],
               "sigma_MPa": st["vm_normal_chip_max_MPa"],
               "sigma_chip": st["vm_normal_chip_name"], "in_window": ok}
        rows.append(row)
        (OUT / "rows.json").write_text(json.dumps(rows, indent=1, ensure_ascii=False) + "\n")
        print("%-26s pkg %6.1fW leak %5.2f ops/s %6.2e HS %5.2f AT %5.2f sig %6.1f %s" % (
            name, br["pkg_W"], br["cpu_leak_W"], br["ops_per_s"], t_hs,
            st["Tmax_C"], st["vm_normal_chip_max_MPa"], "窗内" if ok else "窗外"), flush=True)
    print("WROTE", OUT)


if __name__ == "__main__":
    main()
