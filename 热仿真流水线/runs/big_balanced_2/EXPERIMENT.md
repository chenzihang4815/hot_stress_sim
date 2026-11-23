# 大芯粒均衡方案：各 2 颗

- started (UTC): 2026-10-31T05:04:49.987166+00:00
- CPU T2 274 mm²，GPU GPU_TITANX 471 mm²，CIM jia×64宏 274 mm² / 20.5 mW。
- 各 2 颗，共 6 颗。CIM 面积与 T2 同级，但功耗只有几十毫瓦。
- 链路与全链相同：350 K 起步，只热 SA 120 步，真 ATSim 应力，s 扫 0.5/1/1.5/2。

| 步 | 阶段 | 档K | s | Tmax | 应力MPa | CPU漏电W | dyn/leak | 漏电占比 | 封装W | 操作数 | 轮廓mm | 相对上步 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|

### A 架构回馈（起点 350 K，只热 SA）
| 0 | A 内环0 | 350 | 1 | 77.21 | — | 28.74 | 3.509 | 0.222 | 264.4 | 2.34e+07 | 59.2x48.1 |
- A 外环0: 只热 SA 120 步，起点为当前布局
    SA step  1 BEST rotate CPU_T2_1 Tmax=77.160
    SA step  9 BEST rotate CIM_jia64_1 Tmax=76.910
    SA step 10 BEST nudge GPU_TITANX_0 (1.0,0.0)mm Tmax=76.810
    SA step 17 BEST nudge CPU_T2_1 (0.0,-0.5)mm Tmax=76.730
    SA step 36 BEST swap CPU_T2_1<->GPU_TITANX_1 Tmax=76.520
    SA step 37 BEST nudge GPU_TITANX_0 (0.0,-2.0)mm Tmax=76.400
    SA step 71 BEST swap CIM_jia64_1<->CPU_T2_1 Tmax=75.940
    SA step 72 BEST nudge GPU_TITANX_1 (0.5,0.0)mm Tmax=75.900
    SA step 80 BEST nudge GPU_TITANX_1 (1.0,0.0)mm Tmax=75.880
    SA step 113 BEST nudge GPU_TITANX_1 (0.0,1.0)mm Tmax=75.860
    SA step 115 BEST rotate GPU_TITANX_0 Tmax=75.820
    SA step 117 BEST nudge CPU_T2_0 (-1.0,0.0)mm Tmax=75.790
    SA step 119 BEST swap GPU_TITANX_0<->GPU_TITANX_1 Tmax=75.590
| 1 | A SA0 | 350 | 1 | 75.59 | — | 28.74 | 3.509 | 0.222 | 264.4 | 2.34e+07 | 60.1x58.0 |  ΔT -1.62  Δleak +0.00 W
    SA 77.20→75.59 °C  accepted 62/120
- **收敛** A：温度档 350 K 与布局一致

### B 收敛布局上的真 ATSim
| 2 | B 真ATSim | 350 | 1 | 80.12 | 175.3 | 28.74 | 3.509 | 0.222 | 264.4 | 2.34e+07 | 60.1x58.0 |  ΔT +4.53  Δleak +0.00 W
- 应力 175.3 MPa @ GPU_TITANX_1，翘曲 4.88 µm，ndof 76332，用时 71 s

### C 选频（L* 冻住，每档真 ATSim）
| 3 | C s=0.5 | 340 | 0.5 | 72.22 | 140.8 | 26.21 | 1.924 | 0.342 | 205.8 | 2.34e+07 | 60.1x58.0 |  ΔT -7.90  Δleak -2.53 W
- s=0.5 档 340 K  Tmax 72.22  应力 140.8  Θ 0.667  窗内
| 4 | C s=1 | 350 | 1 | 80.12 | 175.3 | 28.74 | 3.509 | 0.222 | 264.4 | 2.34e+07 | 60.1x58.0 |  ΔT +7.90  Δleak +2.53 W
- s=1 档 350 K  Tmax 80.12  应力 175.3  Θ 1.000  窗内
| 5 | C s=1.5 | 360 | 1.5 | 89.23 | 210.1 | 31.49 | 4.804 | 0.172 | 323.3 | 2.34e+07 | 60.1x58.0 |  ΔT +9.11  Δleak +2.75 W
- s=1.5 档 360 K  Tmax 89.23  应力 210.1  Θ 1.333  窗外
| 6 | C s=2 | 360 | 2 | 97.98 | 243.2 | 31.49 | 6.405 | 0.135 | 379.5 | 2.34e+07 | 60.1x58.0 |  ΔT +8.75  Δleak +0.00 W
- s=2 档 360 K  Tmax 97.98  应力 243.2  Θ 1.667  窗外

日志：`/home/czh/任务仿真产生功耗文件工具/热仿真流水线/runs/big_balanced_2`
