# 5 CPU + ≥3 GPU + 5 CIM

热图：HotSpot 官方 `split_grid_steady.py` + `grid_thermal_map.py`（开尔文）。
布局：L-BFGS。ATSim JSON：`热应力python工程/cases/large_package/`。

## pkg_5cpu_3gpu_5cim（13 颗）

总功耗 **211 W**，块 Tmax **72.2 °C**。

## pkg_5cpu_4gpu_5cim（14 颗，多一颗同规格 3070）

总功耗 **278 W**，块 Tmax **76.6 °C**。

### 5 颗 CPU

| 名字 | 任务 / 机器 | 面积 | 功耗 |
|---|---|---:|---:|
| CPU_matmul | x86 Timing+L1/L2，64×64 矩阵乘 | 13.63 mm² | 1.00 W |
| CPU_threads | x86 Timing+L1/L2，gem5 threads | 13.63 mm² | 0.96 W |
| CPU_armL2 | ARM Timing+L1/L2 hello | 13.63 mm² | 0.87 W |
| CPU_atomic | ARM Atomic hello | 2.90 mm² | 1.14 W |
| CPU_starter | ARM starter atomic hello | 2.90 mm² | 1.65 W |

### GPU

| 名字 | 配置 | 面积 | 功耗（该卡重跑 backprop） |
|---|---|---:|---:|
| GPU_3070 | RTX3070 | 392 mm² | 66.9 W |
| GPU_2060S | RTX2060 Super | 445 mm² | 71.1 W |
| GPU_TITANX | Pascal Titan X | 471 mm² | 67.4 W |
| GPU_3070b | 仅 14 颗方案：第二颗 3070 | 392 mm² | 66.9 W |

### 5 颗 CIM

| 名字 | 面积 | 功耗 |
|---|---:|---:|
| CIM_jia | 4.28 mm² | 2.8 mW |
| CIM_colonnade | 0.23 mm² | 41 mW |
| CIM_wang | 0.12 mm² | 26 mW |
| CIM_basic | 0.047 mm² | 2.6 mW |
| CIM_sinangil | 0.0026 mm² | 0.79 mW |
