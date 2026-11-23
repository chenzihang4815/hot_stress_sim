# L-BFGS 三颗布局：HotSpot vs ATSim 热对比（2026-10-20）

同一套 **lbfgs** 初始布局（CPU / GPU / CIM），先写成 ATSim `chips[]`，再：

- ATSim：`scripts/run_case3.py`（3D FVM，含衬底/中介层）
- HotSpot：用 ATSim 导出的 `.flp/.ptrace/.config`（封装边长与对流电阻对齐）

只比温度。GPU 面积仍是 GV100 815 mm²。

相对误差 = (HotSpot − ATSim) / ATSim × 100%。负号表示 HotSpot 更低。

## Tmax（全域峰值）

| 案例 | 功耗 W | HotSpot Tmax °C | ATSim Tmax °C | Δ °C | 相对 |
|---|---:|---:|---:|---:|---:|
| L2 + CIM basic | 73.83 | 55.31 | 56.61 | −1.30 | −2.29% |
| L2 + CIM jia | 73.83 | 55.38 | 56.64 | −1.26 | −2.23% |
| noCache + CIM jia | 73.63 | 55.35 | 56.58 | −1.23 | −2.17% |
| atomic + CIM wang | 74.12 | 55.42 | 56.60 | −1.18 | −2.08% |
| starter + CIM jia | 74.61 | 56.51 | 57.77 | −1.26 | −2.18% |
| starter + CIM colonnade | 74.65 | 56.32 | 57.96 | −1.64 | −2.83% |

Tmax 相对误差约 **−2.1% ~ −2.8%**（约 1.2–1.6 °C）。HotSpot 始终略低，符合块平均 vs 连续体峰值。

## 分块：HotSpot 块温 vs ATSim 芯片体平均

| 案例 | CPU Δ °C (相对) | GPU Δ °C (相对) | CIM Δ °C (相对) |
|---|---|---|---|
| L2 + basic | −0.13 (−0.25%) | −0.39 (−0.71%) | −0.17 (−0.31%) |
| L2 + jia | −0.17 (−0.32%) | −0.38 (−0.69%) | +0.10 (+0.18%) |
| noCache + jia | −0.42 (−0.77%) | −0.38 (−0.69%) | +0.11 (+0.20%) |
| atomic + wang | −1.01 (−1.79%) | −0.35 (−0.62%) | −0.40 (−0.74%) |
| starter + jia | −1.10 (−1.90%) | −0.38 (−0.68%) | +0.12 (+0.22%) |
| starter + colonnade | −1.47 (−2.54%) | −0.34 (−0.60%) | −0.26 (−0.47%) |

GPU（功耗大头）体平均误差约 **0.6–0.7%**。小而热的 starter CPU 差到约 2.5%，因为 HotSpot 是块平均，ATSim 热点更尖。

## 注意

- ATSim 有有机衬底 + 硅中介层；HotSpot 没有这两层，只作趋势/块温对照。
- 环境 45 °C，R_convec=0.1 K/W，芯片 150 μm，TIM 20 μm。
- 完整数字：`compare_table.csv`；每次运行在 `runs/`。
