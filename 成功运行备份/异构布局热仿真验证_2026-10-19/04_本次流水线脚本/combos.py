"""任务 × 硬件搭配。CPU 一律是 gem5+McPAT，不再把官方 A9 当成第二颗 CPU。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from blocks import (
    BACKUP,
    Block,
    accelforge_from_summary,
    cim_basic_analog,
    cpu_from_gem5_mcpat,
    gpu_backprop,
    load_cpu_macro_blocks,
)

GEM5 = BACKUP / "gem5"


@dataclass(frozen=True)
class Mix:
    id: str
    title: str
    why: str
    blocks: tuple[Block, ...]


def _cpu(rel: str, name: str, notes: str) -> Block:
    return cpu_from_gem5_mcpat(GEM5 / rel, name=name, notes=notes)


CPU_L2 = _cpu(
    "output/mcpat_from_gem5.txt",
    "CPU",
    "gem5 TimingSimpleCPU + L1/L2，hello",
)
CPU_TIMING_NOCACHE = _cpu(
    "多案例/01_se_timing_无cache/mcpat.txt",
    "CPU",
    "gem5 TimingSimpleCPU 无 cache，hello",
)
CPU_ATOMIC = _cpu(
    "多案例/02_se_atomic_无cache/mcpat.txt",
    "CPU",
    "gem5 AtomicSimpleCPU，hello",
)
CPU_STARTER = _cpu(
    "多案例/03_starter_se_atomic/mcpat.txt",
    "CPU",
    "gem5 ARM starter_se atomic，hello",
)


def _cim(combo: str, name: str) -> Block:
    return accelforge_from_summary(combo, kind="CIM", name=name)


def _eyeriss() -> Block:
    return accelforge_from_summary("eyeriss+matmul32", kind="ACCEL", name="ACCEL_Eyeriss")


def all_mixes() -> list[Mix]:
    macros = tuple(load_cpu_macro_blocks())
    return [
        Mix(
            "cpuL2_cimBasic",
            "CPU L1+L2 + CIM basic",
            "面积接近的小系统，布局算法差最容易看出来",
            (CPU_L2, cim_basic_analog()),
        ),
        Mix(
            "cpuL2_cimJia",
            "CPU L1+L2 + CIM jia",
            "CIM 面积 4.3 mm²，和 CPU 核同量级",
            (CPU_L2, _cim("jia_jssc_2020+builtin", "CIM_jia")),
        ),
        Mix(
            "cpuStarter_cimColonnade",
            "CPU starter + CIM colonnade",
            "CPU 动态功耗最高的 gem5 案例 × 功耗较高的 CIM",
            (CPU_STARTER, _cim("colonnade_jssc_2021+builtin", "CIM_colonnade")),
        ),
        Mix(
            "cpuAtomic_cimWang",
            "CPU atomic + CIM wang",
            "无 L2 的小 CPU × 中等 CIM",
            (CPU_ATOMIC, _cim("wang_vlsi_2022+builtin", "CIM_wang")),
        ),
        Mix(
            "cpuTiming_eyeriss_cimBasic",
            "CPU no-cache + Eyeriss + CIM",
            "数字加速器代替 GPU，三块面积都在 mm² 量级",
            (CPU_TIMING_NOCACHE, _eyeriss(), cim_basic_analog()),
        ),
        Mix(
            "cpuL2_gpu_cimBasic",
            "CPU L1+L2 + GPU + CIM",
            "完整异构；GPU 815 mm² / 73 W 会压过另外两块",
            (CPU_L2, gpu_backprop(), cim_basic_analog()),
        ),
        Mix(
            "cpuStarter_gpu_cimJia",
            "CPU starter + GPU + CIM jia",
            "同样有 GPU，CIM 换成面积更大的 jia",
            (CPU_STARTER, gpu_backprop(), _cim("jia_jssc_2020+builtin", "CIM_jia")),
        ),
        Mix(
            "cpuMacros_cimJia",
            "CPU macros + CIM jia",
            "CPU 拆成 4 块再加 CIM，给 B*-tree/MCTS 更多模块",
            macros + (_cim("jia_jssc_2020+builtin", "CIM_jia"),),
        ),
    ]
