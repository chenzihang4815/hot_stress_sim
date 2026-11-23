# AccelForge 基础案例备份（来自完整 accelforge-main）

完整源码：`/home/czh/任务仿真产生功耗文件工具/accelforge-main/accelforge-main`  
官方：[GitHub AccelForge](https://github.com/Accelergy-Project/accelforge) · [文档](https://accelergy-project.github.io/accelforge/)

这不是周期级仿真。输入是 **YAML（硬件 + 任务 + 可选 mapping）**，输出是 **能量、延迟、面积**。

---

## 运行指令（两例都已跑通）

```bash
source /home/czh/miniconda3/etc/profile.d/conda.sh
conda activate accelforge

# 案例 1：最简单数字加速器 toy
python /home/czh/任务仿真产生功耗文件工具/成功运行备份/AccelForge/simple_matmul/scripts/run.py

# 案例 2：CIM basic_analog
python /home/czh/任务仿真产生功耗文件工具/成功运行备份/AccelForge/cim_basic_analog/scripts/run.py
```

环境：conda `accelforge`，Python 3.12，`pip install accelforge`。

---

## 输入是什么、格式如何

**不是一个文件。** 用 Python 把多份 YAML 合成 `Spec`：

| 输入 | 格式 | 含义 |
|---|---|---|
| **arch YAML** | 顶层键 `arch:`，下面 `nodes:` 列表 | 硬件：Memory / Compute / CIM 的 ADC、驱动、阵列 |
| **workload YAML** | 顶层键 `workload:` | 任务：einsum 名、张量形状、位宽 |
| **mapping YAML**（可省） | 顶层键 `mapping:` | 循环怎么铺到哪个存储/计算单元 |

YAML 可用 Jinja（`{% set M = 64 %}`）。CIM 的 arch 文件里往往还带 `workload:` / `variables:`。

### 案例 1 三个输入文件（备份在 `simple_matmul/input/`）

- `simple.yaml`：MainMemory、GlobalBuffer、MAC。能量故意写成 1，面积 0（学接口用）。
- `matmuls.yaml`：一串 matmul；运行时 `N_EINSUMS=1, M=64, KN=32`。
- `unfused_matmuls_to_simple.yaml`：数据放哪、三重循环、在 MAC 上 compute。

### 案例 2 CIM（备份在 `cim_basic_analog/input/`）

- `basic_analog.yaml`：模拟 CIM 宏（ADC、行列驱动、CimUnit）+ 自带 workload。
- `_include.yaml`、`_include_functions.py`、`memory_cells/`：工艺、器件、直方图。

---

## 输出是什么、怎么读

没有 `stats.txt`。Python 对象 `result`，备份成 `output/result.json`。

### 案例 1 `simple_matmul/output/result.json`

| 字段 | 这次的值 | 含义 |
|---|---|---|
| energy | 2228224 | 总能量（YAML 里单位=1 的相对量，**不是焦耳**） |
| latency | 65536 | 延迟（这里 = MAC 次数，因为 MAC throughput=1） |
| n_computes | 65536 | MAC 次数 = 64×32×32 |
| energy_per_component | GlobalBuffer 最大 | 各硬件块能量 |

面积全是 0，因为 toy arch 写死了 `area: 0`。

### 案例 2 `cim_basic_analog/output/result.json`

| 字段 | 这次的值 | 含义 |
|---|---|---|
| energy_J | 2.63e-10 | 总能量 **焦耳** |
| latency_s | 1.00e-7 | 延迟 **秒**（100 ns） |
| n_computes | 1024 | MAC 次数（占满约 32×32 阵列） |
| energy_per_MAC_J | 2.56e-13 | 约 **0.256 pJ/MAC** |
| energy_per_component_J | ADC 最大 | ADC / RowDrivers / CimUnit |
| area_m2 | ADC ≈ 4.56e-8 | **平方米**（ADC ≈ 0.046 mm²） |

和 McPAT 类似：面积看硬件，能量随任务/映射变。CIM 里 ADC 往往占能量大头。

---

## 备份目录

```
成功运行备份/AccelForge/
  README.md
  simple_matmul/input/     三份 YAML
  simple_matmul/output/result.json
  simple_matmul/scripts/run.py
  cim_basic_analog/input/  basic_analog.yaml + include + memory_cells
  cim_basic_analog/output/result.json
  cim_basic_analog/scripts/run.py
```
