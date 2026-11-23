# Accel-Sim 原工程缺件清单（2026-10-07 本机核对）

原路径：`/home/czh/任务仿真产生功耗文件工具/accel-sim-framework-dev/accel-sim-framework-dev`

结论：**不能跑。** 这是 Accel-Sim 前端 git 仓库，仿真后端和输入按官方设计不进 git。

## 已有

- SASS 前端：`gpu-simulator/accel-sim.cc`、`trace-driven/`、`trace-parser/`
- ISA：`ISA_Def/{kepler,pascal,volta,turing,ampere,hopper}_opcode.h`
- 功耗操作码映射：`ISA_Def/accelwattch_component_mapping.h`（只有映射，没有电容/XML）
- 轨迹配置：`gpu-simulator/configs/tested-cfgs/SM7_QV100/trace.config` 等（只有延迟，没有 SM 几何）
- 作业配置：`util/job_launching/configs/define-standard-cfgs.yml`（引用 `$GPGPUSIM_ROOT/.../gpgpusim.config`）
- `get-accel-sim-traces.py`、`setup_environment.sh`
- 本机编译器：g++、cmake、libxml2-dev、libzstd-dev、部分 boost

## 缺少

| 项 | 本机状态 |
|---|---|
| `gpu-simulator/gpgpu-sim` | 目录不存在（gitignore） |
| `accel-sim.out` | 不存在 |
| 任意 `accelwattch_*.xml` | 0 个 |
| 任意 `gpgpusim.config` | 0 个（它们在 GPGPU-Sim 树里） |
| `hw_run/` 轨迹 | 不存在 |
| `kernelslist.g` | 原工程中没有 |
| `gpu-app-collection` | 不存在 |
| CUDA toolkit / `nvcc` | 没有 |
| `/usr/local/cuda` | 没有 |
| NVIDIA GPU / `nvidia-smi` | 没有 |

## 补齐顺序

1. 安装 CUDA toolkit（即使没卡，setup 脚本也要 `CUDA_INSTALL_PATH`）。
2. `source gpu-simulator/setup_environment.sh` → clone GPGPU-Sim。
3. cmake 编译 → `gpu-simulator/bin/release/accel-sim.out`。
4. `./get-accel-sim-traces.py` 下官方 traces（没卡时的唯一轨迹来源）。
5. 功耗用 `QV100-Accelwattch_SASS_SIM` 或 `GV100-Accelwattch_SASS_SIM`，不要用 H100 开 AccelWattch。

不要使用 `最简单案例/Accel-sim/sample_traces/`，那是格式示范。
