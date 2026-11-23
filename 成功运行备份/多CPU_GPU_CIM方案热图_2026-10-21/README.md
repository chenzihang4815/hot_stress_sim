# 多 CPU / GPU / CIM 封装方案（给后续热应力回馈）

单组「1 CPU + 1 GPU + 1 CIM」回馈旋钮太少：GPU 几乎决定 Tmax。下面三套都是 **同一封装里多颗不同硬件/任务**，布局 L-BFGS，热图 = HotSpot `split_grid_steady.py` + `grid_thermal_map.py`。

ATSim JSON：`热应力python工程/cases/multi_chiplet_schemes/`

## 方案 A  `pkg_2cpu_1gpu_2cim`（建议当回馈基线）

5 颗。总功耗 68.9 W，轮廓约 612 mm²，块 Tmax 55.6 °C，网格峰值 329.53 K ≈ 56.4 °C。

| 颗 | 硬件 / 任务 | 面积 | 功耗 | 功率密度 |
|---|---|---:|---:|---:|
| CPU_matmul | x86 Timing+L1/L2，64×64 矩阵乘 | 13.63 mm² | 1.00 W | 7.3 W/cm² |
| CPU_threads | x86 Timing+L1/L2，gem5 threads | 13.63 mm² | 0.96 W | 7.0 W/cm² |
| GPU_3070 | RTX3070，backprop 按该卡重跑 | **392 mm²** | 66.9 W | 17.1 W/cm² |
| CIM_jia | jia JSSC 2020 | 4.28 mm² | 2.8 mW | 0.065 W/cm² |
| CIM_colonnade | colonnade JSSC 2021 | 0.23 mm² | 41 mW | 17.8 W/cm² |

回馈可动：两颗 CPU 的频率/任务；jia vs colonnade 的位置（冷大块 vs 热小块）。

## 方案 B  `pkg_3cpu_1gpu_3cim`

7 颗。总功耗 71.1 W，块 Tmax 56.0 °C，网格峰值 330.1 K ≈ 56.9 °C。

| 颗 | 硬件 / 任务 | 面积 | 功耗 | 功率密度 |
|---|---|---:|---:|---:|
| CPU_matmul | x86 矩阵乘 | 13.63 mm² | 1.00 W | 7.3 W/cm² |
| CPU_threads | x86 threads | 13.63 mm² | 0.96 W | 7.0 W/cm² |
| CPU_starter | ARM starter atomic hello | **2.90 mm²** | **1.65 W** | **56.9 W/cm²** |
| GPU_TITANX | Pascal Titan X，backprop 重跑 | 471 mm² | 67.4 W | 14.3 W/cm² |
| CIM_jia / colonnade / wang | 三种 CIM | 4.28 / 0.23 / 0.12 mm² | 2.8 / 41 / 26 mW | 0.065 / 17.8 / 21.3 W/cm² |

starter 是高密度小核，热图右下角更暗。三种 CIM 给架构/映射旋钮。

## 方案 C  `pkg_2cpu_2gpu_2cim`

6 颗、**两颗 GPU**。总功耗 **140.7 W**，块 Tmax **63.7 °C**，网格峰值 338.0 K ≈ **64.9 °C**。

| 颗 | 硬件 / 任务 | 面积 | 功耗 |
|---|---|---:|---:|
| CPU_matmul | x86 矩阵乘 | 13.63 mm² | 1.00 W |
| CPU_starter | ARM starter | 2.90 mm² | 1.65 W |
| GPU_3070 | RTX3070 | 392 mm² | 66.9 W |
| GPU_2060S | RTX2060 Super | 445 mm² | 71.1 W |
| CIM_jia | jia | 4.28 mm² | 2.8 mW |
| CIM_wang | wang | 0.12 mm² | 26 mW |

两颗 GPU 并排，热点最高，适合试「把两块 GPU 拉开 / 换边」。

## 热图

色标为开尔文（官方脚本）。PNG：本目录 `pkg_*.png`。
