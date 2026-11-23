# 三颗 CPU / GPU / CIM 随机布局 + HotSpot（2026-10-19）

**每个案例固定 3 颗 chiplet：CPU、GPU、CIM。** 没有两颗、没有宏块拆分、没有 Eyeriss 顶替 GPU。

- CPU = gem5 + McPAT（一种 CPU；换 4 种 gem5 机器/任务配置）
- GPU = Accel-Sim AccelWattch backprop（本机目前只有这一份 GPU 功耗）
- CIM = AccelForge（换 4 种宏）
- 布局 = bstree / mcts / lbfgs
- 热仿真 = `任务仿真产生功耗文件工具/HotSpot-master/hotspot`

16 组搭配 × 3 算法 = 48 次 HotSpot。

## 目录

- `00_总表/` CSV + JSON
- `01_算法对比图/` Tmax、轮廓面积、空隙比、ΔT
- `02_HotSpot每次运行/` 每次的 flp / ptrace / block.steady / floorplan.png
- `03_本次用到的功耗源案例/`
- `04_本次流水线脚本/`

## CPU 四种（任务都是 ARM hello）

| id | gem5 | 面积 | 功耗 |
|---|---|---|---|
| L2 | TimingSimpleCPU + L1/L2 | 13.63 mm² | 0.868 W |
| noCache | TimingSimpleCPU 无 cache | 2.90 mm² | 0.676 W |
| atomic | AtomicSimpleCPU | 2.90 mm² | 1.144 W |
| starter | ARM starter_se atomic | 2.90 mm² | 1.651 W |

## GPU（所有案例同一颗）

AccelWattch `kernel_avg_power` = 72.955 W；面积暂用 GV100 裸片 815 mm²（仿真器不算面积）。

## CIM 四种

| id | AccelForge | 面积 | 功耗 |
|---|---|---|---|
| basic | basic_analog | 0.047 mm² | 2.6 mW |
| jia | jia_jssc_2020 | 4.28 mm² | 2.8 mW |
| colonnade | colonnade_jssc_2021 | 0.231 mm² | 41.2 mW |
| wang | wang_vlsi_2022 | 0.122 mm² | 26.0 mW |

案例 id = `{CPU}__{CIM}`，例如 `starter__jia` = starter CPU + GPU backprop + CIM jia。
