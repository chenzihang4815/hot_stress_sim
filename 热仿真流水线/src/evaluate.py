"""对一个设计点给出温度、应力和吞吐。

这是工具的入口。输入是 design.Design 加一个芯粒清单，输出是一次评估：
64 网格复核的峰值温度、最热芯粒的应力、相对吞吐。搜索在内部用 16 网格
场叠加完成，调用方拿不到 16 网格的温度，报出去的只有复核值。

和 run_full_chain.py 的差别只在参数从哪来。那里每个方案的封装、材料、
频率都写在模块顶部，这里全部从 Design 读。旧脚本不动，已有的实验还能复现。

电压还没进功耗。McPAT 的电压是逐部件的，而现有漏电表是在标称电压下取的，
直接按 v² 缩放会和那张表重复计算。design.ChipOp.v 先留着，等接上按电压
重跑的漏电表再启用。
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import run_full_chain as chain
import sa_thermal_feedback as sa
from atsim_sa_proxy import write_proxy_json, eval_atsim_proxy
from design import ChipOp, Design
from hotspot_io import run_hotspot, write_hotspot_inputs

T_LIM_C = 90.0
SIG_LIM_MPA = 200.0


def _scale(c, design: Design) -> float:
    return design.operating.get(c.name, ChipOp()).s


def apply(chips, design: Design, T_K: int) -> None:
    """按设计点重算每颗芯粒的功耗。频率只放大动态部分，漏电只看温度档。"""
    chain.stash(chips)
    for c in chips:
        s = _scale(c, design)
        if c.kind == "CPU":
            m = chain.leak_of(c.model, T_K)
            c.power_W = s * m["dyn_W"] + m["leak_W"]
        elif c.kind == "GPU":
            c.power_W = s * c._dyn0 + c._static + c._const
        else:
            c.power_W = c._p0


def evaluate(chips, design: Design, work: Path, search: bool = True) -> dict:
    """温度闭环、可选的布局搜索、64 网格复核、一次应力。"""
    work.mkdir(parents=True, exist_ok=True)

    T = chain.T_START_K
    tmax = None
    for _ in range(chain.MAX_INNER):
        apply(chips, design, T)
        tmax = _tmax(chips, design, work / "bin", 64)
        nxt = chain.round10_k(tmax)
        if nxt == T:
            break
        T = nxt

    info = None
    if search:
        chips, _, info = chain.heat_sa(chips, work / "search", lambda m: None)
        apply(chips, design, T)
        tmax = _tmax(chips, design, work / "final64", 64)

    stress = _stress(chips, design, work / "stress")
    groups = {"CPU": [], "GPU": [], "CIM": []}
    for c in chips:
        if c.kind in groups:
            groups[c.kind].append(_scale(c, design) if c.kind != "CIM" else 1.0)
    n = sum(len(v) for v in groups.values())
    theta = sum(sum(v) for v in groups.values()) / n if n else 0.0
    sigma = stress.get("vm_normal_chip_max_MPa")

    rec = {
        "Tmax_C": tmax,
        "sigma_MPa": sigma,
        "sigma_chip": stress.get("vm_normal_chip_name"),
        "Theta": theta,
        "McPAT_T_K": T,
        "feasible": tmax is not None and tmax < T_LIM_C and (sigma or 0) < SIG_LIM_MPA,
        "search": info,
        "design": design.to_dict(),
    }
    (work / "eval.json").write_text(
        json.dumps(rec, indent=2, ensure_ascii=False, default=str) + "\n")
    return rec


def _tmax(chips, design, work, rows) -> float:
    if work.exists():
        shutil.rmtree(work)
    write_hotspot_inputs(sa.to_rects(chips), work, rows=rows, cols=rows,
                         package=design.package)
    return float(run_hotspot(chain.HS_BIN, work)["Tmax_C"])


def _stress(chips, design, work) -> dict:
    work.mkdir(parents=True, exist_ok=True)
    jp = work / "case.json"
    write_proxy_json(chips, jp, "eval", design=design)
    st = eval_atsim_proxy(jp, verbose=False)
    return {k: st.get(k) for k in (
        "Tmax_C", "vm_normal_chip_max_MPa", "vm_normal_chip_name", "w_pv_um", "elapsed_s")}
