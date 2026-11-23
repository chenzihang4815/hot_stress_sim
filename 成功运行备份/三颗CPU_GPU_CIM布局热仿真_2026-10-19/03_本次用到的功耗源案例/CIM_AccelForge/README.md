# AccelForge 多组合（CIM + 数字加速器）

官方 `examples/arches/` **没有 CPU、也没有 GPU**。两类硬件：

1. **CIM 宏**：`compute_in_memory/*.yaml`（阵列 + ADC + 驱动）
2. **数字 DNN 加速器**：Eyeriss / Simba / NVDLA / TPU-v4i / Snowcat / simple  
   是卷积/矩阵加速器（PE 阵列、缓冲），不是通用 CPU，也不是 NVIDIA GPU。

任务都是 arch 自带的小 matmul，或 `workloads/basic` 里 32×32 一次 matmul。

---

## 已跑通

| 类型 | 组合 | 面积 mm² | 平均功率 | 说明 |
|---|---|---|---|---|
| CIM | basic_analog + 自带任务 | 0.047 | **2.63 mW** | 先前案例 |
| CIM | colonnade + 自带 | 0.231 | **41.2 mW** | |
| CIM | jia_jssc_2020 + 自带 | 4.28 | **2.77 mW** | 面积最大 |
| CIM | sinangil + 自带 | 0.0026 | **0.79 mW** | 很小的宏 |
| CIM | wang_vlsi_2022 + 自带 | 0.122 | **26.0 mW** | |
| 数字 | eyeriss + matmul32 | 7.79 | E/t ≈ 39.7 mW | 能量/延迟为物理单位 |
| 数字 | simba + matmul32 | 2.22 | E/t ≈ 30.1 mW | 同上 |

CIM 功率 = `energy_J / latency_s`。  
Eyeriss/Simba 的 energy、latency 也是 J / s，功率可同样相除。

---

## 没跑通（本机 pip 包 vs 仓库 YAML 不完全匹配）

| 组合 | 原因 |
|---|---|
| isaac_isca_2016 | 表达式求值失败（`throughput_scale`） |
| tpu_v4i、nvdla | YAML 字段 `separate_read_write_ports` 当前安装的 accelforge 不认 |
| snowcat | 跑完但 latency=0、面积空，能量是相对单位，不当芯片结果 |

完整数字在 `output/summary_table.json`、`all_results.json`。
