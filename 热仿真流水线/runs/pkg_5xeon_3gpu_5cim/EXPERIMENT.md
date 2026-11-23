# 13芯：5 颗 CPU 全部换成 McPAT 官方 Xeon

- started (UTC): 2026-10-26T02:16:09.350709+00:00
- 旧方案 `pkg_5cpu_3gpu_5cim` **不覆盖**（仍是 TimingSimple 13.6 mm²）
- 新方案 id：`pkg_5xeon_3gpu_5cim`
- 每颗 CPU：Xeon.xml 官方结果 Area=421.676 mm²，P=RuntimeDynamic 66.8485 + Leakage 37.0975 = **103.946 W**
- GPU / CIM 与 13 芯原方案相同
- 初始布局：L-BFGS seed=13
- 名字仍用 CPU_matmul 等槽位，便于和旧布局对照，硬件已不是 gem5 小核

## 芯粒

| 名称 | 种类 | 面积 mm² | 功耗 W |
|---|---|---:|---:|
| CPU_matmul | CPU | 421.7 | 103.9 |
| CPU_threads | CPU | 421.7 | 103.9 |
| CPU_armL2 | CPU | 421.7 | 103.9 |
| CPU_atomic | CPU | 421.7 | 103.9 |
| CPU_starter | CPU | 421.7 | 103.9 |
| GPU_3070 | GPU | 392 | 66.89 |
| GPU_2060S | GPU | 445 | 71.12 |
| GPU_TITANX | GPU | 471 | 67.41 |
| CIM_jia | CIM | 4.281 | 0.002773 |
| CIM_colonnade | CIM | 0.2314 | 0.0412 |
| CIM_wang | CIM | 0.122 | 0.026 |
| CIM_basic | CIM | 0.04672 | 0.002625 |
| CIM_sinangil | CIM | 0.002577 | 0.0007948 |

## L-BFGS + HotSpot

- 硅面积合计 **3421 mm²**（旧 13 芯约 1340 mm² 量级，CPU 只有 ~47 mm²）
- 封装功耗合计 **725.2 W**（旧 13 芯约 211 W，CPU 只有 ~5 W）
- 轮廓 5885 mm²  deadspace=0.419
- HotSpot Tmax = **128.740 °C**  Tmin=108.600 °C
- 热图：`/home/czh/任务仿真产生功耗文件工具/热仿真流水线/runs/pkg_5xeon_3gpu_5cim/pkg_5xeon_3gpu_5cim_hotspot.png`
- 布局：`/home/czh/任务仿真产生功耗文件工具/热仿真流水线/runs/pkg_5xeon_3gpu_5cim/layout.json`
- 已登记到 `sa_compare_schemes/pkg_5xeon_3gpu_5cim` 供后续 SA
