# 13芯：频率 s × McPAT 温度 → 任务/漏电占比

- started (UTC): 2026-10-25T14:57:42.173034+00:00
- 布局冻结：`pkg_5cpu_3gpu_5cim` heat-only SA 最优
- s = f/f0，CPU 与 GPU **动态**同乘 s；CIM 不乘
- CPU 漏电：McPAT Total Leakage（改 XML `temperature`）
- GPU 漏电：AccelWattch STATICP，不随 s、也不进 McPAT
- GPU CONSTP、CIM 功耗固定
- 两条臂：`frozen350` 漏电锁 350 K；`mcpatT` 用 HotSpot Tmax 回写 McPAT 直到 10 K 档稳定
- 看的量：任务/漏电 = P_dyn/P_leak；漏电占比 = P_leak/(P_dyn+P_leak)（CPU）或 /P_pkg

## s = 0.5
- frozen350: Tmax=68.040C  CPU dyn/leak=0.292  leak_frac_cpu=0.774  pkg dyn/leak=0.213
- mcpatT T=340K (2 inner): Tmax=67.980C  CPU dyn/leak=0.343  leak_frac_cpu=0.745  CPU leak=2.9186W

## s = 0.75
- frozen350: Tmax=68.980C  CPU dyn/leak=0.438  leak_frac_cpu=0.695  pkg dyn/leak=0.320
- mcpatT T=340K (1 inner): Tmax=68.930C  CPU dyn/leak=0.514  leak_frac_cpu=0.660  CPU leak=2.9186W

## s = 1.0
- frozen350: Tmax=70.020C  CPU dyn/leak=0.584  leak_frac_cpu=0.631  pkg dyn/leak=0.426
- mcpatT T=340K (1 inner): Tmax=69.870C  CPU dyn/leak=0.686  leak_frac_cpu=0.593  CPU leak=2.9186W

## s = 1.25
- frozen350: Tmax=71.220C  CPU dyn/leak=0.730  leak_frac_cpu=0.578  pkg dyn/leak=0.533
- mcpatT T=340K (1 inner): Tmax=71.070C  CPU dyn/leak=0.857  leak_frac_cpu=0.538  CPU leak=2.9186W

## s = 1.5
- frozen350: Tmax=72.430C  CPU dyn/leak=0.876  leak_frac_cpu=0.533  pkg dyn/leak=0.639
- mcpatT T=350K (2 inner): Tmax=72.430C  CPU dyn/leak=0.876  leak_frac_cpu=0.533  CPU leak=3.4269W

## s = 1.75
- frozen350: Tmax=73.640C  CPU dyn/leak=1.022  leak_frac_cpu=0.495  pkg dyn/leak=0.746
- mcpatT T=350K (1 inner): Tmax=73.640C  CPU dyn/leak=1.022  leak_frac_cpu=0.495  CPU leak=3.4269W

## s = 2.0
- frozen350: Tmax=74.840C  CPU dyn/leak=1.168  leak_frac_cpu=0.461  pkg dyn/leak=0.853
- mcpatT T=350K (1 inner): Tmax=74.840C  CPU dyn/leak=1.168  leak_frac_cpu=0.461  CPU leak=3.4269W

## s = 2.25
- frozen350: Tmax=76.050C  CPU dyn/leak=1.314  leak_frac_cpu=0.432  pkg dyn/leak=0.959
- mcpatT T=350K (1 inner): Tmax=76.050C  CPU dyn/leak=1.314  leak_frac_cpu=0.432  CPU leak=3.4269W

## s = 2.5
- frozen350: Tmax=77.250C  CPU dyn/leak=1.460  leak_frac_cpu=0.407  pkg dyn/leak=1.066
- mcpatT T=350K (1 inner): Tmax=77.250C  CPU dyn/leak=1.460  leak_frac_cpu=0.407  CPU leak=3.4269W

## s = 3.0
- frozen350: Tmax=79.660C  CPU dyn/leak=1.752  leak_frac_cpu=0.363  pkg dyn/leak=1.279
- mcpatT T=350K (1 inner): Tmax=79.660C  CPU dyn/leak=1.752  leak_frac_cpu=0.363  CPU leak=3.4269W

## s = 4.0
- frozen350: Tmax=84.490C  CPU dyn/leak=2.336  leak_frac_cpu=0.300  pkg dyn/leak=1.705
- mcpatT T=360K (2 inner): Tmax=84.650C  CPU dyn/leak=2.019  leak_frac_cpu=0.331  CPU leak=3.9645W

## frozen350（漏电不随温度）

| s | McPAT T | Tmax | CPU dyn | CPU leak | CPU dyn/leak | CPU 漏电占比 | 封装 dyn | 封装 leak | 封装 dyn/leak | 封装漏电占比 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 350 K | 68.04 °C | 1.001 | 3.427 | 0.292 | 0.774 | 17.02 | 79.83 | 0.213 | 0.412 |
| 0.75 | 350 K | 68.98 °C | 1.501 | 3.427 | 0.438 | 0.695 | 25.52 | 79.83 | 0.320 | 0.394 |
| 1 | 350 K | 70.02 °C | 2.001 | 3.427 | 0.584 | 0.631 | 34.03 | 79.83 | 0.426 | 0.379 |
| 1.25 | 350 K | 71.22 °C | 2.502 | 3.427 | 0.730 | 0.578 | 42.54 | 79.83 | 0.533 | 0.364 |
| 1.5 | 350 K | 72.43 °C | 3.002 | 3.427 | 0.876 | 0.533 | 51.05 | 79.83 | 0.639 | 0.350 |
| 1.75 | 350 K | 73.64 °C | 3.502 | 3.427 | 1.022 | 0.495 | 59.56 | 79.83 | 0.746 | 0.338 |
| 2 | 350 K | 74.84 °C | 4.003 | 3.427 | 1.168 | 0.461 | 68.07 | 79.83 | 0.853 | 0.326 |
| 2.25 | 350 K | 76.05 °C | 4.503 | 3.427 | 1.314 | 0.432 | 76.57 | 79.83 | 0.959 | 0.315 |
| 2.5 | 350 K | 77.25 °C | 5.003 | 3.427 | 1.460 | 0.407 | 85.08 | 79.83 | 1.066 | 0.305 |
| 3 | 350 K | 79.66 °C | 6.004 | 3.427 | 1.752 | 0.363 | 102.10 | 79.83 | 1.279 | 0.286 |
| 4 | 350 K | 84.49 °C | 8.005 | 3.427 | 2.336 | 0.300 | 136.13 | 79.83 | 1.705 | 0.255 |

## mcpatT（Tmax 回写 350K）

| s | McPAT T | Tmax | CPU dyn | CPU leak | CPU dyn/leak | CPU 漏电占比 | 封装 dyn | 封装 leak | 封装 dyn/leak | 封装漏电占比 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 340 K | 67.98 °C | 1.001 | 2.919 | 0.343 | 0.745 | 17.02 | 79.32 | 0.215 | 0.410 |
| 0.75 | 340 K | 68.93 °C | 1.501 | 2.919 | 0.514 | 0.660 | 25.52 | 79.32 | 0.322 | 0.393 |
| 1 | 340 K | 69.87 °C | 2.001 | 2.919 | 0.686 | 0.593 | 34.03 | 79.32 | 0.429 | 0.377 |
| 1.25 | 340 K | 71.07 °C | 2.502 | 2.919 | 0.857 | 0.538 | 42.54 | 79.32 | 0.536 | 0.362 |
| 1.5 | 350 K | 72.43 °C | 3.002 | 3.427 | 0.876 | 0.533 | 51.05 | 79.83 | 0.639 | 0.350 |
| 1.75 | 350 K | 73.64 °C | 3.502 | 3.427 | 1.022 | 0.495 | 59.56 | 79.83 | 0.746 | 0.338 |
| 2 | 350 K | 74.84 °C | 4.003 | 3.427 | 1.168 | 0.461 | 68.07 | 79.83 | 0.853 | 0.326 |
| 2.25 | 350 K | 76.05 °C | 4.503 | 3.427 | 1.314 | 0.432 | 76.57 | 79.83 | 0.959 | 0.315 |
| 2.5 | 350 K | 77.25 °C | 5.003 | 3.427 | 1.460 | 0.407 | 85.08 | 79.83 | 1.066 | 0.305 |
| 3 | 350 K | 79.66 °C | 6.004 | 3.427 | 1.752 | 0.363 | 102.10 | 79.83 | 1.279 | 0.286 |
| 4 | 360 K | 84.65 °C | 8.005 | 3.964 | 2.019 | 0.331 | 136.13 | 80.37 | 1.694 | 0.256 |

## 升频会不会让任务/漏电比下降？

在本扫点上，CPU dyn/leak **没有** 因升频下降：漏电台阶相对动态的 s 倍增不够大。
封装 dyn/leak 本扫点上也是升频单调升（GPU 动态跟着 s 涨，STATIC 漏电固定）。

对照：`frozen350` 里漏电常数，dyn∝s，**dyn/leak 必随 s 单调升**。
`mcpatT` 里漏电按 10 K 档跳，理论上可以在档边界把比压下去；是否发生取决于 Tmax 有没有跨档。

