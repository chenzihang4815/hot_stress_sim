# CPU + GPU + CIM 三芯片案例（L-BFGS 初始布局）

每个 JSON 都是 **正好 3 颗**：`CPU`、`GPU`、`CIM`。

- CPU：gem5 + McPAT（4 种机器，任务都是 ARM hello）
- GPU：Accel-Sim AccelWattch backprop（16 组同一颗；面积 815 mm² 为 GV100 文献值）
- CIM：AccelForge 4 种宏
- 坐标：`defaults.py` 里拍板的 **lbfgs**，seed=13

默认入口（给 ATSim 直接跑）：

```bash
python scripts/run_case3.py cases/case3_cpu_gpu_cim.json
```

这是 `L2__jia`：CPU 13.6 mm² / 0.87 W，GPU 815 mm² / 73.0 W，CIM 4.28 mm² / 2.8 mW。

其余 16 组在本目录，文件名 = `{CPU配置}__{CIM宏}.json`。

| 文件 | CPU | CIM |
|---|---|---|
| `L2__basic.json` | CPU Timing L1+L2 | CIM basic |
| `L2__jia.json` | CPU Timing L1+L2 | CIM jia |
| `L2__colonnade.json` | CPU Timing L1+L2 | CIM colonnade |
| `L2__wang.json` | CPU Timing L1+L2 | CIM wang |
| `noCache__basic.json` | CPU Timing no-cache | CIM basic |
| `noCache__jia.json` | CPU Timing no-cache | CIM jia |
| `noCache__colonnade.json` | CPU Timing no-cache | CIM colonnade |
| `noCache__wang.json` | CPU Timing no-cache | CIM wang |
| `atomic__basic.json` | CPU Atomic | CIM basic |
| `atomic__jia.json` | CPU Atomic | CIM jia |
| `atomic__colonnade.json` | CPU Atomic | CIM colonnade |
| `atomic__wang.json` | CPU Atomic | CIM wang |
| `starter__basic.json` | CPU starter | CIM basic |
| `starter__jia.json` | CPU starter | CIM jia |
| `starter__colonnade.json` | CPU starter | CIM colonnade |
| `starter__wang.json` | CPU starter | CIM wang |

跑某一组：

```bash
python scripts/run_case3.py cases/cpu_gpu_cim_lbfgs/starter__jia.json
```
