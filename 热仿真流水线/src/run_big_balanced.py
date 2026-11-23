#!/usr/bin/env python3
"""Big balanced packages: n CPUs + n GPUs + n large CIMs, full chain.

CIM is jia at 64 macros (274 mm^2, 20.5 mW), so all three kinds are large.
CPU is T2 (274 mm^2), GPU is TITANX (471 mm^2). One n per run, chosen by the
first argument. The chain is the same as run_full_chain: 350 K start, heat-only
SA to L*, real ATSim, frequency sweep, re-SA from L*.
"""

from __future__ import annotations

import json
import shutil
import sys
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
SRC = ROOT / "热仿真流水线" / "src"
sys.path.insert(0, str(SRC))

import sa_thermal_feedback as sa  # noqa: E402
from blocks import Block  # noqa: E402
from hotspot_io import run_hotspot, write_hotspot_inputs  # noqa: E402
from run_full_chain import (  # noqa: E402
    COLORS, GPU_AREA, GPU_CONST_W, GPU_POWER, GPU_STATIC_W, N_SA, OFFICIAL_380,
    SEED, T_START_K, Ledger, apply_power, lbfgs_clamped, leak_of, log_md,
    real_stress, round10_k, stash,
)
import run_full_chain as C  # noqa: E402

JIA64 = (273.9985, 0.020516)  # 64 macros, from cim_macros_scale
GPU = "GPU_TITANX"
CPU = "T2"


def build_balanced(n: int):
    """n of each kind, every chip uniquely named so the floorplan is valid."""
    blocks, cpu_xml, gpu_model = [], {}, {}
    for i in range(n):
        name = f"CPU_T2_{i}"
        off = OFFICIAL_380[CPU]
        blocks.append(Block(name, "CPU", f"McPAT {CPU}", off["area_mm2"],
                            off["dyn_W"] + off["sub_W"] + off["gate_W"], CPU))
        cpu_xml[name] = CPU
    for i in range(n):
        name = f"GPU_TITANX_{i}"
        blocks.append(Block(name, "GPU", "AccelWattch", GPU_AREA[GPU], GPU_POWER[GPU], GPU))
        gpu_model[name] = GPU
    for i in range(n):
        blocks.append(Block(f"CIM_jia64_{i}", "CIM", "AccelForge", JIA64[0], JIA64[1], "jia x64"))
    rects, _ = lbfgs_clamped(blocks, seed=SEED)
    chips = [sa.Chip(r.name, r.kind, r.x, r.y, r.w, r.h, r.power_W, r.area_mm2) for r in rects]
    for c in chips:
        if c.kind == "GPU":
            c.model = gpu_model[c.name]
    return chips, cpu_xml


def stash_by_model(chips) -> None:
    for c in chips:
        if c.kind == "GPU":
            m = c.model
            c._dyn0 = max(0.0, c.power_W - GPU_STATIC_W[m] - GPU_CONST_W[m])
            c._static = GPU_STATIC_W[m]
            c._const = GPU_CONST_W[m]
        else:
            c._p0 = c.power_W


def main() -> None:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    out = ROOT / "热仿真流水线" / "runs" / f"big_balanced_{n}"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    C.OUT = out  # real_stress and the chain write under OUT

    # Point the chain's power functions at model-aware versions.
    C.stash = stash_by_model
    orig_apply = C.apply_power

    def apply_named(chips, s, cpu_xml, T_K):
        # apply_power keys GPU static by chip name; give it the model name.
        saved = {}
        for c in chips:
            if c.kind == "GPU":
                saved[id(c)] = c.name
                c.name = c.model
        try:
            return orig_apply(chips, s, cpu_xml, T_K)
        finally:
            for c in chips:
                if id(c) in saved:
                    c.name = saved[id(c)]

    C.apply_power = apply_named

    spec = {"id": f"big_{n}", "label": f"{n}×T2+{n}×TITANX+{n}×jia64",
            "cpus": [CPU] * n, "gpus": [GPU] * n}
    md = out / "EXPERIMENT.md"
    log_md(md, f"# 大芯粒均衡方案：各 {n} 颗")
    log_md(md, "")
    log_md(md, f"- started (UTC): {datetime.now(timezone.utc).isoformat()}")
    log_md(md, f"- CPU {CPU} {OFFICIAL_380[CPU]['area_mm2']:.0f} mm²，"
               f"GPU {GPU} {GPU_AREA[GPU]:.0f} mm²，CIM jia×64宏 {JIA64[0]:.0f} mm² / {JIA64[1]*1e3:.1f} mW。")
    log_md(md, f"- 各 {n} 颗，共 {3*n} 颗。CIM 面积与 T2 同级，但功耗只有几十毫瓦。")
    log_md(md, f"- 链路与全链相同：350 K 起步，只热 SA {N_SA} 步，真 ATSim 应力，s 扫 0.5/1/1.5/2。")
    log_md(md, "")

    chips, cpu_xml = build_balanced(n)
    stash_by_model(chips)
    (out / "lbfgs_init.json").write_text(
        json.dumps(sa.snapshot(chips), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    # Reuse the chain's per-scheme driver by handing it a spec it understands,
    # but with our chips already built. run_scheme rebuilds chips itself, so
    # drive the stages here instead.
    from run_full_chain import close_loop, theta, S_LIST, T_LIM_C, SIG_LIM_MPA, MAX_INNER

    slot = out
    ledger = Ledger(slot, md)
    log_md(md, "| 步 | 阶段 | 档K | s | Tmax | 应力MPa | CPU漏电W | dyn/leak | 漏电占比 | 封装W | 操作数 | 轮廓mm | 相对上步 |")
    log_md(md, "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|")

    log_md(md, "")
    log_md(md, "### A 架构回馈（起点 350 K，只热 SA）")
    chips, T_star, br, tmax = close_loop(
        chips, 1.0, cpu_xml, T_START_K, slot / "A", md, ledger, "A", with_sa=True)
    (slot / "Lstar.json").write_text(
        json.dumps(sa.snapshot(chips), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    Lstar = deepcopy(chips)

    log_md(md, "")
    log_md(md, "### B 收敛布局上的真 ATSim")
    st = real_stress(chips, slot / "B", "Lstar")
    ledger.add("B 真ATSim", chips, br, st["Tmax_C"], st)
    log_md(md, f"- 应力 {st['vm_normal_chip_max_MPa']:.1f} MPa @ {st['vm_normal_chip_name']}，"
               f"翘曲 {st['w_pv_um']:.2f} µm，ndof {st['ndof']}，用时 {st['elapsed_s']:.0f} s")

    log_md(md, "")
    log_md(md, "### C 选频（L* 冻住，每档真 ATSim）")
    for s in S_LIST:
        T = T_star
        cur = deepcopy(Lstar)
        for i in range(MAX_INNER):
            b = C.apply_power(cur, s, cpu_xml, T)
            t_hs = sa.eval_hotspot(cur, slot / f"C_s{s}_i{i}")
            nxt = round10_k(t_hs)
            if nxt == T:
                break
            T = nxt
        st_s = real_stress(cur, slot / "C", f"s{s}")
        ok = st_s["Tmax_C"] < T_LIM_C and (st_s["vm_normal_chip_max_MPa"] or 1e9) < SIG_LIM_MPA
        ledger.add(f"C s={s:g}", cur, b, st_s["Tmax_C"], st_s)
        log_md(md, f"- s={s:g} 档 {T} K  Tmax {st_s['Tmax_C']:.2f}  "
                   f"应力 {st_s['vm_normal_chip_max_MPa']:.1f}  Θ {theta(s):.3f}  "
                   f"{'窗内' if ok else '窗外'}")

    (slot / "ledger.json").write_text(
        json.dumps(ledger.rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    log_md(md, "")
    log_md(md, f"日志：`{out}`")
    print("WROTE", out)


if __name__ == "__main__":
    main()
