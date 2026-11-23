# 13芯 McPAT 温度闭环（冻结 SA 最优布局）

- started (UTC): 2026-11-10T10:59:38.560277+00:00
- 封装：`pkg_5cpu_3gpu_5cim`
- 布局：**冻结** heat-only SA 最优（Tmax 原 70.04 °C，功耗来自默认 350 K McPAT）
- 每轮：用本轮 HotSpot **封装 Tmax** 圆整到 McPAT 10 K 档，改 5 颗 CPU 的 XML `temperature`，重跑 McPAT
- 用新的 Runtime Dynamic + Total Leakage 替换 CPU 功耗，**GPU/CIM 不变**，再跑 HotSpot
- 停止：McPAT 温度档与上一轮相同，或达到 8 轮
- gem5 操作数 / 仿真时间 **不随 McPAT T 变**

## Round 0: McPAT T = 350 K
- CPU_matmul: dyn=0.217977 W  leak=0.78389 W  tot=1.00187 W  leak_frac=0.782  N_op=1297322 (unchanged)
- CPU_threads: dyn=0.172267 W  leak=0.78389 W  tot=0.956157 W  leak_frac=0.820  N_op=1808180 (unchanged)
- CPU_armL2: dyn=0.056212 W  leak=0.6197 W  tot=0.675912 W  leak_frac=0.917  N_op=4994 (unchanged)
- CPU_atomic: dyn=0.523847 W  leak=0.6197 W  tot=1.14355 W  leak_frac=0.542  N_op=5020 (unchanged)
- CPU_starter: dyn=1.03098 W  leak=0.6197 W  tot=1.65068 W  leak_frac=0.375  N_op=5020 (unchanged)
- 封装功耗 CPU=5.4282 W (dyn 2.0013+leak 3.4269)  GPU=205.413 W  CIM=0.0734 W  tot=210.915 W
- HotSpot Tmax = **58.460 °C** (331.61 K) → 下一档 McPAT **330 K**

## Round 1: McPAT T = 330 K
- CPU_matmul: dyn=0.217977 W  leak=0.566811 W  tot=0.784788 W  leak_frac=0.722  N_op=1297322 (unchanged)
- CPU_threads: dyn=0.172267 W  leak=0.566811 W  tot=0.739078 W  leak_frac=0.767  N_op=1808180 (unchanged)
- CPU_armL2: dyn=0.056212 W  leak=0.448838 W  tot=0.50505 W  leak_frac=0.889  N_op=4994 (unchanged)
- CPU_atomic: dyn=0.523847 W  leak=0.448838 W  tot=0.972685 W  leak_frac=0.461  N_op=5020 (unchanged)
- CPU_starter: dyn=1.03098 W  leak=0.448838 W  tot=1.47982 W  leak_frac=0.303  N_op=5020 (unchanged)
- 封装功耗 CPU=4.4814 W (dyn 2.0013+leak 2.4801)  GPU=205.413 W  CIM=0.0734 W  tot=209.968 W
- HotSpot Tmax = **58.410 °C** (331.56 K) → 下一档 McPAT **330 K**
- **收敛**：下一档仍是 330 K，停止。

## 各轮对照

| 轮 | McPAT T | CPU dyn | CPU leak | CPU 合计 | 封装合计 | HotSpot Tmax | 下一档 | 操作数 |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0 | 350 K | 2.0013 | 3.4269 | 5.4282 | 210.915 | 58.460 °C | 330 K | 五颗均不变 |
| 1 | 330 K | 2.0013 | 2.4801 | 4.4814 | 209.968 | 58.410 °C | 330 K | 五颗均不变 |

## 规律（写给下一轮查阅）
1. 相对官方 350 K：末轮 CPU 漏电 -0.9467 W，封装 Tmax -0.0500 K。
2. 操作数来自 gem5，McPAT 温度改不了 `simInsts` / `simSeconds`。
3. Runtime Dynamic 若各轮相同，说明 McPAT 动态功耗不吃 temperature。
4. GPU/CIM 功耗未改；封装 Tmax 主要由 GPU 决定，CPU 漏电变化对 Tmax 可能很小。

日志目录：`/home/czh/任务仿真产生功耗文件工具/热仿真流水线/runs/mcpat_13chip_tmax_loop_10K_anchored`
