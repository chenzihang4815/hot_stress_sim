"""把四套已跑通的功耗/面积结果收成 HotSpot 可用的矩形块。"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

ROOT = Path("/home/czh/任务仿真产生功耗文件工具")
BACKUP = ROOT / "成功运行备份"

# AccelWattch 不输出面积。QV100 配置对应 Volta / Quadro GV100 一类，
# 公开裸片面积约 815 mm²（NVIDIA GV100）。这是文献尺寸，不是仿真器算出来的。
GPU_DIE_MM2 = 815.0


@dataclass
class Block:
    name: str
    kind: str
    source: str
    area_mm2: float
    power_W: float
    notes: str

    @property
    def area_m2(self) -> float:
        return self.area_mm2 * 1e-6


def _mcpat_processor(txt: Path, name: str, kind: str, source: str, notes: str) -> Block:
    area = power_dyn = leakage = None
    in_processor = False
    for raw in txt.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if line.startswith("Processor:"):
            in_processor = True
            continue
        if in_processor and (
            line.startswith("Total Cores")
            or line.startswith("Total L2")
            or line.startswith("Total First")
            or line.startswith("Total NoC")
            or line.startswith("****")
        ):
            break
        if not in_processor:
            continue
        if line.startswith("Area"):
            area = float(line.split("=")[1].split("mm")[0])
        elif line.startswith("Runtime Dynamic"):
            power_dyn = float(line.split("=")[1].split("W")[0])
        elif line.startswith("Total Leakage"):
            leakage = float(line.split("=")[1].split("W")[0])
    if area is None or power_dyn is None or leakage is None:
        raise ValueError(f"未能从 {txt} 解析 Processor 级面积/功耗")
    return Block(
        name=name,
        kind=kind,
        source=source,
        area_mm2=area,
        power_W=power_dyn + leakage,
        notes=notes + f"；P=RuntimeDynamic({power_dyn:.6g})+Leakage({leakage:.6g})",
    )


def _gpu_from_accelwattch(log: Path) -> Block:
    avg = None
    kernel = ""
    for raw in log.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.startswith("kernel_name"):
            kernel = raw.split("=", 1)[1].strip()
        if raw.startswith("kernel_avg_power"):
            avg = float(raw.split("=")[1].strip())
            break
    if avg is None:
        raise ValueError(f"未能从 {log} 解析 kernel_avg_power")
    return Block(
        name="GPU_AccelSim",
        kind="GPU",
        source=str(log),
        area_mm2=GPU_DIE_MM2,
        power_W=avg,
        notes=(
            f"backprop kernel={kernel}；功率来自 AccelWattch kernel_avg_power；"
            f"面积不是仿真器输出，采用 GV100 裸片 {GPU_DIE_MM2} mm²"
        ),
    )


def _cim_from_accelforge(result: Path) -> Block:
    data = json.loads(result.read_text(encoding="utf-8"))
    energy = float(data["energy_J"])
    latency = float(data["latency_s"])
    area_m2 = sum(float(v) for v in data["area_m2"].values())
    if latency <= 0 or area_m2 <= 0:
        raise ValueError(f"CIM 结果能量/延迟/面积异常: {result}")
    return Block(
        name="CIM_AccelForge",
        kind="CIM",
        source=str(result),
        area_mm2=area_m2 * 1e6,
        power_W=energy / latency,
        notes="功率=energy_J/latency_s；面积=各组件 area_m2 之和（MainMemory=0 不占面积）",
    )


def cpu_from_gem5_mcpat(txt: Path, name: str, notes: str) -> Block:
    """gem5 不算功耗：stats → McPAT XML → 这一份 txt 才是 CPU chiplet。"""
    return _mcpat_processor(txt, name=name, kind="CPU", source=str(txt), notes=notes)


def gpu_backprop() -> Block:
    return _gpu_from_accelwattch(
        BACKUP / "Accel-Sim/官方下载/run_backprop/accelwattch_power_report.log"
    )


def cim_basic_analog() -> Block:
    return _cim_from_accelforge(BACKUP / "AccelForge/cim_basic_analog/output/result.json")


def accelforge_from_summary(combo: str, kind: str, name: str) -> Block:
    rows = json.loads((BACKUP / "AccelForge/more_combos/output/summary_table.json").read_text(encoding="utf-8"))
    hit = next((r for r in rows if r["combo"] == combo), None)
    if hit is None or hit.get("area_total_mm2") in (None, 0) or hit.get("avg_power_W") in (None, 0):
        raise ValueError(f"AccelForge 组合不可用: {combo}")
    return Block(
        name=name,
        kind=kind,
        source=str(BACKUP / "AccelForge/more_combos/output/summary_table.json") + f"::{combo}",
        area_mm2=float(hit["area_total_mm2"]),
        power_W=float(hit["avg_power_W"]),
        notes=f"AccelForge {combo}；P=energy/latency={hit['avg_power_W']:.6g} W",
    )


def load_four_chiplets() -> list[Block]:
    """旧入口：一种 CPU（gem5+McPAT）+ GPU + CIM。不再把官方 McPAT A9 当成第二颗 CPU。"""
    return [
        cpu_from_gem5_mcpat(
            BACKUP / "gem5/output/mcpat_from_gem5.txt",
            name="CPU",
            notes="gem5 TimingSimpleCPU + L1/L2，hello；stats 转 McPAT",
        ),
        gpu_backprop(),
        cim_basic_analog(),
    ]


def load_cpu_macro_blocks(txt: Path | None = None) -> list[Block]:
    """gem5+McPAT 粗粒度宏块（Core / L2 / NoC / MC），不是官方 A9 第二颗 CPU。"""
    txt = txt or (BACKUP / "gem5/output/mcpat_from_gem5.txt")
    cores = _section_metrics(txt, "Total Cores:")
    l2 = _section_metrics(txt, "Total L2s:")
    noc = _section_metrics(txt, "Total NoCs")
    mc = _section_metrics(txt, "Total MCs:")
    return [
        Block("CPU_Core", "CPU", str(txt), cores["area"], cores["dyn"] + cores["leak"], "Total Cores"),
        Block("CPU_L2", "CPU", str(txt), l2["area"], l2["dyn"] + l2["leak"], "L2"),
        Block("CPU_NoC", "CPU", str(txt), noc["area"], noc["dyn"] + noc["leak"], "NoC/Bus"),
        Block("CPU_MC", "CPU", str(txt), mc["area"], mc["dyn"] + mc["leak"], "Memory Controller"),
    ]


def _section_metrics(txt: Path, header: str) -> dict[str, float]:
    lines = txt.read_text(encoding="utf-8", errors="replace").splitlines()
    start = None
    for i, raw in enumerate(lines):
        if header in raw:
            start = i
            break
    if start is None:
        raise ValueError(f"{txt} 中没有 {header}")
    area = dyn = leak_s = leak_g = None
    for raw in lines[start + 1 : start + 16]:
        line = raw.strip()
        if not line:
            if area is not None and dyn is not None:
                break
            continue
        if line.startswith("****") or (
            area is not None
            and (
                line.startswith("Total ")
                or (line.endswith(":") and not line[0].isspace() and "=" not in line)
            )
        ):
            break
        if line.startswith("Area"):
            area = float(line.split("=")[1].split("mm")[0])
        elif line.startswith("Runtime Dynamic"):
            dyn = float(line.split("=")[1].split("W")[0])
        elif line.startswith("Subthreshold Leakage"):
            leak_s = float(line.split("=")[1].split("W")[0])
        elif line.startswith("Gate Leakage"):
            leak_g = float(line.split("=")[1].split("W")[0])
    if None in (area, dyn, leak_s, leak_g):
        raise ValueError(f"解析 {header} 失败")
    return {"area": area, "dyn": dyn, "leak": leak_s + leak_g}


def dump_blocks(blocks: list[Block], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = []
    for b in blocks:
        row = asdict(b)
        row["side_if_square_mm"] = math.sqrt(b.area_mm2)
        row["power_density_W_per_mm2"] = b.power_W / b.area_mm2 if b.area_mm2 else 0.0
        payload.append(row)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
