# 大 CPU 闭环：频率 → 热 → McPAT 漏电 → SA 布局

- started (UTC): 2026-10-26T14:09:51.245196+00:00
- 每套、每个 s：McPAT 温度档 ↔ HotSpot 直到 10 K 档稳定，然后 **只热 SA 16 步**，再回写 T；外环最多 3 次
- CPU 动态 = s × McPAT Runtime Dynamic(T)；漏电 = McPAT Total Leakage(T)，不随 s
- GPU 动态 ×s，STATIC 当漏电不乘；CIM 不乘
- 本闭环 **不跑 ATSim**（12 GB）；布局只服务 Tmax/漏电
- 漏电 **不是** 现场跑 mcpat。官方模板 device_type=0（HP），McPAT 里温度只改 I_off（cacti/technology.cc），HP 栅漏电不随温度，面积和动态功耗也不随温度。一档模板要跑 >1 h，重复的是与温度无关的 CACTI 分区搜索。
- 所以每档用官方 380 K ExampleResults：亚阈值漏电 × ITRS HP I_off 表比值，栅漏电固定。380 K 处与官方 Total Leakage 逐位一致（脚本启动时核对）。这就是二进制会打印的数。
- 窗：Tmax < 90.0 °C。窗内赢家 = 最低 CPU leak_frac（任务占比最高）

## 5×Penryn  (`pkg_5penryn_3gpu_5cim`)
### s=0.5
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=90.410 °C → 360 K  CPU dyn/leak=2.119 leak_frac=0.321
    inner 1: McPAT 360 K → Tmax=89.000 °C → 360 K  CPU dyn/leak=2.616 leak_frac=0.277
- outer 0: 只热 SA 16 步
    SA Tmax 88.580→88.580 °C  next McPAT 360 K  accepted 6/16
- **收敛** T 档 360 K，布局+漏电一致
- 结果 s=0.5  McPAT=360 K  Tmax=88.58 °C  CPU dyn/leak=2.616  leak_frac_cpu=0.277  P_pkg=355.1 W  OK

### s=1
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=108.090 °C → 380 K  CPU dyn/leak=4.238 leak_frac=0.191
- outer 0: 只热 SA 16 步
    SA Tmax 107.360→107.360 °C  next McPAT 380 K  accepted 6/16
- **收敛** T 档 380 K，布局+漏电一致
- 结果 s=1  McPAT=380 K  Tmax=107.36 °C  CPU dyn/leak=4.238  leak_frac_cpu=0.191  P_pkg=501.7 W  T=107.4>90

### s=1.5
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=125.770 °C → 400 K  CPU dyn/leak=6.357 leak_frac=0.136
    inner 1: McPAT 400 K → Tmax=127.350 °C → 400 K  CPU dyn/leak=5.238 leak_frac=0.160
- outer 0: 只热 SA 16 步
    SA Tmax 126.340→126.340 °C  next McPAT 400 K  accepted 6/16
- **收敛** T 档 400 K，布局+漏电一致
- 结果 s=1.5  McPAT=400 K  Tmax=126.34 °C  CPU dyn/leak=5.238  leak_frac_cpu=0.160  P_pkg=649.6 W  T=126.3>90

### s=2
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=143.440 °C → 400 K  CPU dyn/leak=8.476 leak_frac=0.106
    inner 1: McPAT 400 K → Tmax=145.030 °C → 400 K  CPU dyn/leak=6.984 leak_frac=0.125
- outer 0: 只热 SA 16 步
    SA step  1 BEST rotate CPU_armL2 Tmax=143.830
    SA Tmax 143.860→143.830 °C  next McPAT 400 K  accepted 6/16
- **收敛** T 档 400 K，布局+漏电一致
- 结果 s=2  McPAT=400 K  Tmax=143.83 °C  CPU dyn/leak=6.984  leak_frac_cpu=0.125  P_pkg=785.5 W  T=143.8>90

## 5×T2  (`pkg_5t2_3gpu_5cim`)
### s=0.5
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=93.680 °C → 370 K  CPU dyn/leak=1.357 leak_frac=0.424
    inner 1: McPAT 370 K → Tmax=92.860 °C → 370 K  CPU dyn/leak=1.468 leak_frac=0.405
- outer 0: 只热 SA 16 步
    SA step 13 BEST jump GPU_TITANX Tmax=91.490
    SA step 14 BEST nudge GPU_3070 (0.0,-0.5)mm Tmax=91.400
    SA Tmax 92.850→91.400 °C  next McPAT 360 K  accepted 6/16
- outer 1: McPAT 闭环（布局暂冻）
    inner 0: McPAT 360 K → Tmax=90.550 °C → 360 K  CPU dyn/leak=1.601 leak_frac=0.384
- outer 1: 只热 SA 16 步
    SA Tmax 90.550→90.550 °C  next McPAT 360 K  accepted 10/16
- **收敛** T 档 360 K，布局+漏电一致
- 结果 s=0.5  McPAT=360 K  Tmax=90.55 °C  CPU dyn/leak=1.601  leak_frac_cpu=0.384  P_pkg=394.2 W  T=90.6>90

### s=1
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=110.540 °C → 380 K  CPU dyn/leak=2.714 leak_frac=0.269
- outer 0: 只热 SA 16 步
    SA step 13 BEST jump GPU_TITANX Tmax=109.190
    SA step 14 BEST nudge GPU_3070 (0.0,-0.5)mm Tmax=109.050
    SA Tmax 110.520→109.050 °C  next McPAT 380 K  accepted 6/16
- **收敛** T 档 380 K，布局+漏电一致
- 结果 s=1  McPAT=380 K  Tmax=109.05 °C  CPU dyn/leak=2.714  leak_frac_cpu=0.269  P_pkg=550.5 W  T=109.1>90

### s=1.5
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=127.390 °C → 400 K  CPU dyn/leak=4.072 leak_frac=0.197
    inner 1: McPAT 400 K → Tmax=128.720 °C → 400 K  CPU dyn/leak=3.631 leak_frac=0.216
- outer 0: 只热 SA 16 步
    SA step 13 BEST jump GPU_TITANX Tmax=127.400
    SA step 14 BEST nudge GPU_3070 (0.0,-0.5)mm Tmax=127.210
    SA Tmax 128.700→127.210 °C  next McPAT 400 K  accepted 6/16
- **收敛** T 档 400 K，布局+漏电一致
- 结果 s=1.5  McPAT=400 K  Tmax=127.21 °C  CPU dyn/leak=3.631  leak_frac_cpu=0.216  P_pkg=703.8 W  T=127.2>90

### s=2
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=144.250 °C → 400 K  CPU dyn/leak=5.429 leak_frac=0.156
    inner 1: McPAT 400 K → Tmax=145.580 °C → 400 K  CPU dyn/leak=4.841 leak_frac=0.171
- outer 0: 只热 SA 16 步
    SA step 13 BEST jump GPU_TITANX Tmax=144.260
    SA step 14 BEST nudge GPU_3070 (0.0,-0.5)mm Tmax=144.030
    SA Tmax 145.540→144.030 °C  next McPAT 400 K  accepted 6/16
- **收敛** T 档 400 K，布局+漏电一致
- 结果 s=2  McPAT=400 K  Tmax=144.03 °C  CPU dyn/leak=4.841  leak_frac_cpu=0.171  P_pkg=845.9 W  T=144.0>90

## 2P+2T+1X  (`pkg_2penryn_2t2_1xeon_3gpu_5cim`)
### s=0.5
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=95.500 °C → 370 K  CPU dyn/leak=1.360 leak_frac=0.424
    inner 1: McPAT 370 K → Tmax=94.510 °C → 370 K  CPU dyn/leak=1.488 leak_frac=0.402
- outer 0: 只热 SA 16 步
    SA step  8 BEST jump GPU_2060S Tmax=91.850
    SA Tmax 94.360→91.850 °C  next McPAT 360 K  accepted 9/16
- outer 1: McPAT 闭环（布局暂冻）
    inner 0: McPAT 360 K → Tmax=90.880 °C → 360 K  CPU dyn/leak=1.642 leak_frac=0.378
- outer 1: 只热 SA 16 步
    SA Tmax 90.880→90.880 °C  next McPAT 360 K  accepted 6/16
- **收敛** T 档 360 K，布局+漏电一致
- 结果 s=0.5  McPAT=360 K  Tmax=90.88 °C  CPU dyn/leak=1.642  leak_frac_cpu=0.378  P_pkg=401.5 W  T=90.9>90

### s=1
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=114.140 °C → 390 K  CPU dyn/leak=2.721 leak_frac=0.269
    inner 1: McPAT 390 K → Tmax=115.100 °C → 390 K  CPU dyn/leak=2.517 leak_frac=0.284
- outer 0: 只热 SA 16 步
    SA step  2 BEST nudge GPU_TITANX (0.5,0.0)mm Tmax=114.890
    SA step  8 BEST jump GPU_2060S Tmax=112.060
    SA Tmax 114.900→112.060 °C  next McPAT 390 K  accepted 9/16
- **收敛** T 档 390 K，布局+漏电一致
- 结果 s=1  McPAT=390 K  Tmax=112.06 °C  CPU dyn/leak=2.517  leak_frac_cpu=0.284  P_pkg=573.8 W  T=112.1>90

### s=1.5
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=132.780 °C → 400 K  CPU dyn/leak=4.081 leak_frac=0.197
    inner 1: McPAT 400 K → Tmax=134.600 °C → 400 K  CPU dyn/leak=3.545 leak_frac=0.220
- outer 0: 只热 SA 16 步
    SA step  8 BEST jump GPU_2060S Tmax=131.190
    SA Tmax 134.350→131.190 °C  next McPAT 400 K  accepted 9/16
- **收敛** T 档 400 K，布局+漏电一致
- 结果 s=1.5  McPAT=400 K  Tmax=131.19 °C  CPU dyn/leak=3.545  leak_frac_cpu=0.220  P_pkg=728.3 W  T=131.2>90

### s=2
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=151.430 °C → 400 K  CPU dyn/leak=5.441 leak_frac=0.155
    inner 1: McPAT 400 K → Tmax=153.240 °C → 400 K  CPU dyn/leak=4.727 leak_frac=0.175
- outer 0: 只热 SA 16 步
    SA step  8 BEST jump GPU_2060S Tmax=149.460
    SA Tmax 152.940→149.460 °C  next McPAT 400 K  accepted 9/16
- **收敛** T 档 400 K，布局+漏电一致
- 结果 s=2  McPAT=400 K  Tmax=149.46 °C  CPU dyn/leak=4.727  leak_frac_cpu=0.175  P_pkg=876.1 W  T=149.5>90

## 5×Xeon  (`pkg_5xeon_3gpu_5cim`)
### s=0.5
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=107.790 °C → 380 K  CPU dyn/leak=0.901 leak_frac=0.526
- outer 0: 只热 SA 16 步
    SA Tmax 107.760→107.760 °C  next McPAT 380 K  accepted 8/16
- **收敛** T 档 380 K，布局+漏电一致
- 结果 s=0.5  McPAT=380 K  Tmax=107.76 °C  CPU dyn/leak=0.901  leak_frac_cpu=0.526  P_pkg=542.1 W  T=107.8>90

### s=1
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=128.740 °C → 400 K  CPU dyn/leak=1.802 leak_frac=0.357
    inner 1: McPAT 400 K → Tmax=131.760 °C → 400 K  CPU dyn/leak=1.577 leak_frac=0.388
- outer 0: 只热 SA 16 步
    SA Tmax 131.730→131.730 °C  next McPAT 400 K  accepted 8/16
- **收敛** T 档 400 K，布局+漏电一致
- 结果 s=1  McPAT=400 K  Tmax=131.73 °C  CPU dyn/leak=1.577  leak_frac_cpu=0.388  P_pkg=751.7 W  T=131.7>90

### s=1.5
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=149.700 °C → 400 K  CPU dyn/leak=2.703 leak_frac=0.270
    inner 1: McPAT 400 K → Tmax=152.710 °C → 400 K  CPU dyn/leak=2.365 leak_frac=0.297
- outer 0: 只热 SA 16 步
    SA Tmax 152.690→152.690 °C  next McPAT 400 K  accepted 8/16
- **收敛** T 档 400 K，布局+漏电一致
- 结果 s=1.5  McPAT=400 K  Tmax=152.69 °C  CPU dyn/leak=2.365  leak_frac_cpu=0.297  P_pkg=934.8 W  T=152.7>90

### s=2
- outer 0: McPAT 闭环（布局暂冻）
    inner 0: McPAT 380 K → Tmax=170.660 °C → 400 K  CPU dyn/leak=3.604 leak_frac=0.217
    inner 1: McPAT 400 K → Tmax=173.670 °C → 400 K  CPU dyn/leak=3.153 leak_frac=0.241
- outer 0: 只热 SA 16 步
    SA Tmax 173.640→173.640 °C  next McPAT 400 K  accepted 8/16
- **收敛** T 档 400 K，布局+漏电一致
- 结果 s=2  McPAT=400 K  Tmax=173.64 °C  CPU dyn/leak=3.153  leak_frac_cpu=0.241  P_pkg=1118.0 W  T=173.6>90

## 总表

| 方案 | s | McPAT T | Tmax | CPU dyn/leak | CPU 漏电占比 | 封装 W | 窗 |
|---|---:|---:|---:|---:|---:|---:|---|
| 5×Penryn | 0.5 | 360 K | 88.58 | 2.616 | 0.277 | 355.1 | 是 |
| 5×Penryn | 1 | 380 K | 107.36 | 4.238 | 0.191 | 501.7 | 否 |
| 5×Penryn | 1.5 | 400 K | 126.34 | 5.238 | 0.160 | 649.6 | 否 |
| 5×Penryn | 2 | 400 K | 143.83 | 6.984 | 0.125 | 785.5 | 否 |
| 5×T2 | 0.5 | 360 K | 90.55 | 1.601 | 0.384 | 394.2 | 否 |
| 5×T2 | 1 | 380 K | 109.05 | 2.714 | 0.269 | 550.5 | 否 |
| 5×T2 | 1.5 | 400 K | 127.21 | 3.631 | 0.216 | 703.8 | 否 |
| 5×T2 | 2 | 400 K | 144.03 | 4.841 | 0.171 | 845.9 | 否 |
| 2P+2T+1X | 0.5 | 360 K | 90.88 | 1.642 | 0.378 | 401.5 | 否 |
| 2P+2T+1X | 1 | 390 K | 112.06 | 2.517 | 0.284 | 573.8 | 否 |
| 2P+2T+1X | 1.5 | 400 K | 131.19 | 3.545 | 0.220 | 728.3 | 否 |
| 2P+2T+1X | 2 | 400 K | 149.46 | 4.727 | 0.175 | 876.1 | 否 |
| 5×Xeon | 0.5 | 380 K | 107.76 | 0.901 | 0.526 | 542.1 | 否 |
| 5×Xeon | 1 | 400 K | 131.73 | 1.577 | 0.388 | 751.7 | 否 |
| 5×Xeon | 1.5 | 400 K | 152.69 | 2.365 | 0.297 | 934.8 | 否 |
| 5×Xeon | 2 | 400 K | 173.64 | 3.153 | 0.241 | 1118.0 | 否 |

## 窗内更好的方案

- 5×Penryn: s=0.5  Tmax=88.58  leak_frac=0.277  dyn/leak=2.616 ← 赢家
