"""每个案例固定三颗：CPU、GPU、CIM。

CPU = gem5 + McPAT（一种 CPU，换任务/机器配置）。
GPU = Accel-Sim AccelWattch（本机目前只有 backprop 一份功耗）。
CIM = AccelForge（换宏/任务组合）。
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from blocks import (
    BACKUP,
    Block,
    accelforge_from_summary,
    cim_basic_analog,
    cpu_from_gem5_mcpat,
    gpu_backprop,
)

GEM5 = BACKUP / "gem5"


@dataclass(frozen=True)
class Mix:
    id: str
    title: str
    why: str
    blocks: tuple[Block, ...]


def _cpu(rel: str, notes: str) -> Block:
    b = cpu_from_gem5_mcpat(GEM5 / rel, name="CPU", notes=notes)
    return replace(b, name="CPU", kind="CPU")


def _gpu() -> Block:
    b = gpu_backprop()
    return replace(b, name="GPU", kind="GPU")


def _cim_basic() -> Block:
    b = cim_basic_analog()
    return replace(b, name="CIM", kind="CIM")


def _cim(combo: str, notes: str) -> Block:
    b = accelforge_from_summary(combo, kind="CIM", name="CIM")
    return replace(b, name="CIM", kind="CIM", notes=notes + "；" + b.notes)


CPU_VARIANTS = [
    (
        "L2",
        "CPU Timing L1+L2",
        _cpu("output/mcpat_from_gem5.txt", "gem5 TimingSimpleCPU + L1/L2，hello"),
    ),
    (
        "noCache",
        "CPU Timing no-cache",
        _cpu("多案例/01_se_timing_无cache/mcpat.txt", "gem5 TimingSimpleCPU 无 cache，hello"),
    ),
    (
        "atomic",
        "CPU Atomic",
        _cpu("多案例/02_se_atomic_无cache/mcpat.txt", "gem5 AtomicSimpleCPU，hello"),
    ),
    (
        "starter",
        "CPU starter",
        _cpu("多案例/03_starter_se_atomic/mcpat.txt", "gem5 ARM starter_se atomic，hello"),
    ),
]

CIM_VARIANTS = [
    ("basic", "CIM basic_analog", _cim_basic()),
    (
        "jia",
        "CIM jia",
        _cim("jia_jssc_2020+builtin", "AccelForge jia_jssc_2020 自带任务"),
    ),
    (
        "colonnade",
        "CIM colonnade",
        _cim("colonnade_jssc_2021+builtin", "AccelForge colonnade 自带任务"),
    ),
    (
        "wang",
        "CIM wang",
        _cim("wang_vlsi_2022+builtin", "AccelForge wang_vlsi_2022 自带任务"),
    ),
]


def all_mixes() -> list[Mix]:
    gpu = _gpu()
    mixes: list[Mix] = []
    for cpu_id, cpu_title, cpu in CPU_VARIANTS:
        for cim_id, cim_title, cim in CIM_VARIANTS:
            blocks = (cpu, gpu, cim)
            kinds = [b.kind for b in blocks]
            names = [b.name for b in blocks]
            if kinds != ["CPU", "GPU", "CIM"] or names != ["CPU", "GPU", "CIM"]:
                raise ValueError(f"每个案例必须正好是 CPU/GPU/CIM 三颗，得到 {names} {kinds}")
            mixes.append(
                Mix(
                    id=f"{cpu_id}__{cim_id}",
                    title=f"{cpu_title} + GPU backprop + {cim_title}",
                    why="三颗 chiplet：换 CPU 的 gem5 配置和 CIM 宏，GPU 目前只有 backprop",
                    blocks=blocks,
                )
            )
    return mixes
