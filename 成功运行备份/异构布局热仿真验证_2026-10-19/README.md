# 异构布局热仿真验证备份（2026-10-19）

CPU = gem5 + McPAT（一种 CPU，不再把官方 ARM_A9 当第二颗）。
GPU = Accel-Sim AccelWattch backprop。
CIM / 数字加速器 = AccelForge。
布局 = bstree / mcts / lbfgs 随机初始布局。
热仿真 = 本目录上级的 HotSpot-master 副本。

## 目录

- `00_总表/`  8 组搭配 × 3 算法的 CSV / JSON
- `01_算法对比图/`  Tmax、轮廓面积、空隙比、ΔT
- `02_HotSpot每次运行/`  24 次：flp / ptrace / block.steady / floorplan.png
- `03_本次用到的功耗源案例/`  实际喂给布局的 CPU/GPU/CIM 输出
- `04_本次流水线脚本/`  当时用的 Python

## 8 组搭配

| id | 硬件/任务 | 为什么选 |
|---|---|---|
| cpuL2_cimBasic | gem5 Timing+L1/L2 hello + CIM basic_analog | 小系统，算法差最容易看 |
| cpuL2_cimJia | 同上 CPU + CIM jia_jssc_2020 | CIM 4.3 mm²，和 CPU 同量级 |
| cpuStarter_cimColonnade | gem5 starter atomic hello + CIM colonnade | CPU 动态功耗最高 × 较热 CIM |
| cpuAtomic_cimWang | gem5 Atomic hello + CIM wang | 无 L2 小 CPU |
| cpuTiming_eyeriss_cimBasic | gem5 Timing 无 cache + Eyeriss matmul32 + CIM basic | 用数字加速器代替 GPU |
| cpuL2_gpu_cimBasic | Timing+L1/L2 + GPU backprop + CIM basic | 完整异构；GPU 815 mm² / 73 W 会压过其余 |
| cpuStarter_gpu_cimJia | starter CPU + GPU + CIM jia | 同样有 GPU，CIM 更大 |
| cpuMacros_cimJia | CPU 拆 Core/L2/NoC/MC + CIM jia | 给 B*-tree/MCTS 更多模块 |

## 功耗源文件

- CPU Timing+L1/L2: `成功运行备份/gem5/output/mcpat_from_gem5.txt`（任务 hello）
- CPU Timing 无 cache: `gem5/多案例/01_se_timing_无cache/mcpat.txt`
- CPU Atomic: `gem5/多案例/02_se_atomic_无cache/mcpat.txt`
- CPU starter: `gem5/多案例/03_starter_se_atomic/mcpat.txt`
- GPU: `Accel-Sim/官方下载/run_backprop/accelwattch_power_report.log`
- CIM basic: `AccelForge/cim_basic_analog/output/result.json`
- 其余 CIM / Eyeriss: `AccelForge/more_combos/output/summary_table.json`

GPU 面积不是仿真器给的，暂用 GV100 裸片 815 mm²。
