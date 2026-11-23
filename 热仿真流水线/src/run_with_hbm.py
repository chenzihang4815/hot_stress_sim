#!/usr/bin/env python3
"""big_balanced_2 plus one HBM3 stack, full chain.

Same 2x T2 + 2x TITANX + 2x jia64 as run_big_balanced.py, with one HBM die
added as a fourth chiplet kind. HBM power is a fixed stack power and is never
scaled by frequency, because apply_power only scales CPU and GPU. Writes to
runs/with_hbm_2 and does not touch runs/big_balanced_2.

HBM assumptions:
- area 11 x 11 mm, the Samsung 8-Hi HBM3 cube. McPAT has no DRAM die model.
- power is one Ramulator channel at peak, times the 16 channels of an HBM3
  stack. Both channel powers were measured on HBM3_16Gb_8hi / HBM3_6400Mbps
  with the Ayna model, single channel, requests sent back to back, so this is
  an upper bound, not the average over a task.
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
from run_full_chain import (  # noqa: E402
    GPU_AREA, GPU_POWER, MAX_INNER, N_SA, OFFICIAL_380, S_LIST, SEED,
    SIG_LIM_MPA, T_LIM_C, T_START_K, Ledger, close_loop, lbfgs_clamped,
    log_md, real_stress, round10_k, theta,
)
import run_full_chain as C  # noqa: E402
from run_big_balanced import GPU, JIA64, stash_by_model  # noqa: E402

# One channel, back to back, measured 2026-11-02.
# Background (refresh + standby, ~0.088 W) is paid once per stack; the
# traffic-dependent part scales with who is driving the stack.
CHAN_BG_W = 0.088
CHAN_W = {
    "gem5_matmul192": 1.035,   # CPU task alone
    "gpu_backprop": 0.844,     # GPU task alone
    # Both at once: one background plus each task's traffic.
    "cpu_plus_gpu": CHAN_BG_W + (1.035 - CHAN_BG_W) + (0.844 - CHAN_BG_W),
}
N_CHAN = 16  # HBM3: 16 channels, 32 pseudo-channels
HBM_AREA_MM2 = 11.0 * 11.0


def build(hbm_w: float, which: str):
    blocks, cpu_xml, gpu_model = [], {}, {}
    for i in range(2):
        name = f"CPU_T2_{i}"
        off = OFFICIAL_380["T2"]
        blocks.append(Block(name, "CPU", "McPAT T2", off["area_mm2"],
                            off["dyn_W"] + off["sub_W"] + off["gate_W"], "T2"))
        cpu_xml[name] = "T2"
    for i in range(2):
        name = f"GPU_TITANX_{i}"
        blocks.append(Block(name, "GPU", "AccelWattch", GPU_AREA[GPU],
                            GPU_POWER[GPU], GPU))
        gpu_model[name] = GPU
    for i in range(2):
        blocks.append(Block(f"CIM_jia64_{i}", "CIM", "AccelForge",
                            JIA64[0], JIA64[1], "jia x64"))
    blocks.append(Block("HBM3_stack", "HBM", "Ramulator2.1",
                        HBM_AREA_MM2, hbm_w, which))
    rects, _ = lbfgs_clamped(blocks, seed=SEED)
    chips = [sa.Chip(r.name, r.kind, r.x, r.y, r.w, r.h, r.power_W, r.area_mm2)
             for r in rects]
    for c in chips:
        if c.kind == "GPU":
            c.model = gpu_model[c.name]
    return chips, cpu_xml


def main() -> None:
    which = sys.argv[1] if len(sys.argv) > 1 else "cpu_plus_gpu"
    if which not in CHAN_W:
        raise SystemExit(f"unknown trace {which}, choose from {list(CHAN_W)}")
    hbm_w = CHAN_W[which] * N_CHAN

    tag = "" if which == "gem5_matmul192" else f"_{which}"
    out = ROOT / "热仿真流水线" / "runs" / f"with_hbm_2{tag}"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    C.OUT = out
    C.stash = stash_by_model
    base_apply = C.apply_power

    def apply_named(chips, s, cpu_xml, T_K):
        saved = {}
        for c in chips:
            if c.kind == "GPU":
                saved[id(c)] = c.name
                c.name = c.model
        try:
            return base_apply(chips, s, cpu_xml, T_K)
        finally:
            for c in chips:
                if id(c) in saved:
                    c.name = saved[id(c)]
    C.apply_power = apply_named

    md = out / "EXPERIMENT.md"
    log_md(md, "# 大芯粒均衡方案加一颗 HBM3：2×T2 + 2×TITANX + 2×jia64 + 1×HBM3")
    log_md(md, "")
    log_md(md, f"- started (UTC): {datetime.now(timezone.utc).isoformat()}")
    log_md(md, f"- HBM3 16Gb 8hi 6400Mbps，面积 11×11 mm = {HBM_AREA_MM2:.0f} mm²。"
               f"McPAT 没有 DRAM 颗粒模型，面积取三星 8 层 HBM3 的商用尺寸。")
    log_md(md, f"- 功耗：Ramulator 单通道峰值 {CHAN_W[which]:.3f} W × {N_CHAN} 通道"
               f" = {hbm_w:.2f} W。trace 为 {which}，请求背靠背发送，Ayna 用默认电气参数。"
               f"16 通道按同时打满计，是上界，不是任务全程的平均。")
    log_md(md, "- HBM 功耗不随频率缩放，频率环只动 CPU 和 GPU。")
    log_md(md, f"- 其余与全链相同：350 K 起步，只热 SA {N_SA} 步，真 ATSim 应力，s 扫 0.5/1/1.5/2。")
    log_md(md, "")

    chips, cpu_xml = build(hbm_w, which)
    stash_by_model(chips)
    (out / "lbfgs_init.json").write_text(
        json.dumps(sa.snapshot(chips), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    ledger = Ledger(out, md)
    log_md(md, "| 步 | 阶段 | 档K | s | Tmax | 应力MPa | CPU漏电W | dyn/leak | 漏电占比 | 封装W | 操作数 | 轮廓mm | 相对上步 |")
    log_md(md, "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|")

    log_md(md, "")
    log_md(md, "### A 架构回馈（起点 350 K，只热 SA）")
    chips, T_star, br, tmax = close_loop(
        chips, 1.0, cpu_xml, T_START_K, out / "A", md, ledger, "A", with_sa=True)
    (out / "Lstar.json").write_text(
        json.dumps(sa.snapshot(chips), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    Lstar = deepcopy(chips)

    log_md(md, "")
    log_md(md, "### B 收敛布局上的真 ATSim")
    st = real_stress(chips, out / "B", "Lstar")
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
            t_hs = sa.eval_hotspot(cur, out / f"C_s{s}_i{i}")
            nxt = round10_k(t_hs)
            if nxt == T:
                break
            T = nxt
        st_s = real_stress(cur, out / "C", f"s{s}")
        ok = st_s["Tmax_C"] < T_LIM_C and (st_s["vm_normal_chip_max_MPa"] or 1e9) < SIG_LIM_MPA
        ledger.add(f"C s={s:g}", cur, b, st_s["Tmax_C"], st_s)
        log_md(md, f"- s={s:g} 档 {T} K  Tmax {st_s['Tmax_C']:.2f}  "
                   f"应力 {st_s['vm_normal_chip_max_MPa']:.1f}  Θ {theta(s):.3f}  "
                   f"{'窗内' if ok else '窗外'}")

    (out / "ledger.json").write_text(
        json.dumps(ledger.rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    log_md(md, "")
    log_md(md, f"日志：`{out}`")
    print("WROTE", out, "HBM", f"{hbm_w:.2f} W")


if __name__ == "__main__":
    main()
