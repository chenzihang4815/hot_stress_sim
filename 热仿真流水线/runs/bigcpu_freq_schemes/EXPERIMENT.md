# 大 CPU 多方案 × 升/降频

- started (UTC): 2026-10-26T02:38:00.160749+00:00
- GPU/CIM 与 13 芯相同；CPU 换成 McPAT 官方大核模板（不是 gem5 任务）
- Penryn 91.1 mm²，P=47.93+11.31=59.24 W；Niagara2/T2 273.5 mm²，P=50.42+18.58=69.00 W；Xeon 421.7 mm²，P=66.85+37.10=103.95 W
- 没有真正 500+ mm² 的官方模板；混搭里最大的一颗用 Xeon 422 mm²
- s=f/f0，CPU/GPU **动态** ×s，漏电与 GPU CONST/STATIC、CIM 不乘
- 布局冻结后只重跑 HotSpot（12 GB，本对照不做 ATSim）
- 窗口记录：Tmax < 90.0 °C（仅标注，不作为唯一目标）

- 建成 `pkg_5penryn_3gpu_5cim`  CPU面积=456 mm²  CPU功耗=296.2 W  封装 501.7 W / 1768 mm²  Tmax=108.09 °C
- 建成 `pkg_5t2_3gpu_5cim`  CPU面积=1368 mm²  CPU功耗=345.0 W  封装 550.5 W / 2680 mm²  Tmax=110.54 °C
- 建成 `pkg_2penryn_2t2_1xeon_3gpu_5cim`  CPU面积=1151 mm²  CPU功耗=360.4 W  封装 565.9 W / 2464 mm²  Tmax=114.14 °C
- 复用 `pkg_5xeon_3gpu_5cim`  CPU面积=422 mm²  Tmax=128.74 °C

## 频率扫描

### pkg_5penryn_3gpu_5cim
- s=0.5  Tmax=90.41 °C  P_pkg=365.8 W  CPU dyn/leak=2.119  leak_frac_cpu=0.321  T=90.4>90
- s=0.75  Tmax=99.25 °C  P_pkg=433.8 W  CPU dyn/leak=3.179  leak_frac_cpu=0.239  T=99.2>90
- s=1  Tmax=108.09 °C  P_pkg=501.7 W  CPU dyn/leak=4.238  leak_frac_cpu=0.191  T=108.1>90
- s=1.25  Tmax=116.93 °C  P_pkg=569.6 W  CPU dyn/leak=5.298  leak_frac_cpu=0.159  T=116.9>90
- s=1.5  Tmax=125.77 °C  P_pkg=637.5 W  CPU dyn/leak=6.357  leak_frac_cpu=0.136  T=125.8>90
- s=2  Tmax=143.44 °C  P_pkg=773.4 W  CPU dyn/leak=8.476  leak_frac_cpu=0.106  T=143.4>90
- 本扫点 **没有** s 落在 Tmax<90（s=0.5 仍超）

### pkg_5t2_3gpu_5cim
- s=0.5  Tmax=93.68 °C  P_pkg=408.4 W  CPU dyn/leak=1.357  leak_frac_cpu=0.424  T=93.7>90
- s=0.75  Tmax=102.11 °C  P_pkg=479.4 W  CPU dyn/leak=2.036  leak_frac_cpu=0.329  T=102.1>90
- s=1  Tmax=110.54 °C  P_pkg=550.5 W  CPU dyn/leak=2.714  leak_frac_cpu=0.269  T=110.5>90
- s=1.25  Tmax=118.97 °C  P_pkg=621.5 W  CPU dyn/leak=3.393  leak_frac_cpu=0.228  T=119.0>90
- s=1.5  Tmax=127.39 °C  P_pkg=692.5 W  CPU dyn/leak=4.072  leak_frac_cpu=0.197  T=127.4>90
- s=2  Tmax=144.25 °C  P_pkg=834.6 W  CPU dyn/leak=5.429  leak_frac_cpu=0.156  T=144.2>90
- 本扫点 **没有** s 落在 Tmax<90（s=0.5 仍超）

### pkg_2penryn_2t2_1xeon_3gpu_5cim
- s=0.5  Tmax=95.50 °C  P_pkg=418.1 W  CPU dyn/leak=1.360  leak_frac_cpu=0.424  T=95.5>90
- s=0.75  Tmax=104.82 °C  P_pkg=492.0 W  CPU dyn/leak=2.041  leak_frac_cpu=0.329  T=104.8>90
- s=1  Tmax=114.14 °C  P_pkg=565.9 W  CPU dyn/leak=2.721  leak_frac_cpu=0.269  T=114.1>90
- s=1.25  Tmax=123.46 °C  P_pkg=639.8 W  CPU dyn/leak=3.401  leak_frac_cpu=0.227  T=123.5>90
- s=1.5  Tmax=132.78 °C  P_pkg=713.7 W  CPU dyn/leak=4.081  leak_frac_cpu=0.197  T=132.8>90
- s=2  Tmax=151.43 °C  P_pkg=861.5 W  CPU dyn/leak=5.441  leak_frac_cpu=0.155  T=151.4>90
- 本扫点 **没有** s 落在 Tmax<90（s=0.5 仍超）

### pkg_5xeon_3gpu_5cim
- s=0.5  Tmax=107.79 °C  P_pkg=542.1 W  CPU dyn/leak=0.901  leak_frac_cpu=0.526  T=107.8>90
- s=0.75  Tmax=118.26 °C  P_pkg=633.6 W  CPU dyn/leak=1.351  leak_frac_cpu=0.425  T=118.3>90
- s=1  Tmax=128.74 °C  P_pkg=725.2 W  CPU dyn/leak=1.802  leak_frac_cpu=0.357  T=128.7>90
- s=1.25  Tmax=139.22 °C  P_pkg=816.8 W  CPU dyn/leak=2.252  leak_frac_cpu=0.307  T=139.2>90
- s=1.5  Tmax=149.70 °C  P_pkg=908.4 W  CPU dyn/leak=2.703  leak_frac_cpu=0.270  T=149.7>90
- s=2  Tmax=170.66 °C  P_pkg=1091.5 W  CPU dyn/leak=3.604  leak_frac_cpu=0.217  T=170.7>90
- 本扫点 **没有** s 落在 Tmax<90（s=0.5 仍超）

## 总表（s=1）

| 方案 | CPU 模型 | CPU 面积 | CPU 功耗 s=1 | 封装功耗 | Tmax s=1 | 窗内最高 s |
|---|---|---:|---:|---:|---:|---:|
| 5×Penryn | 5×Penryn | 456 | 296.2 | 501.7 | 108.1 | 无 |
| 5×T2 | 5×Niagara2 | 1368 | 345.0 | 550.5 | 110.5 | 无 |
| 2P+2T+1X | 2×Penryn+2×T2+1×Xeon | 1151 | 360.4 | 565.9 | 114.1 | 无 |
| 5×Xeon | 5×Xeon | 2108 | 519.7 | 725.2 | 128.7 | 无 |

图：`/home/czh/任务仿真产生功耗文件工具/热仿真流水线/runs/bigcpu_freq_schemes/freq_tmax_and_ratio.png`  `/home/czh/任务仿真产生功耗文件工具/热仿真流水线/runs/bigcpu_freq_schemes/s1_tmax.png`
