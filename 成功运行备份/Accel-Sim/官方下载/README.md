# Accel-Sim：用官方指令拿到的输入（2026-10-08）

本机 **没有 CUDA / nvcc**，编不出 `accel-sim.out`，所以 **还没有仿真输出**。  
下面是按官方 README / `get-accel-sim-traces.py` 已经下好、可观察的输入。

## 用了哪些官方指令

轨迹清单和最小套件（Purdue FTP，脚本默认地址）：

```bash
# 等价于 ./get-accel-sim-traces.py -a tesla-v100/rodinia_2.0-ft
wget ftp://ftp.ecn.purdue.edu/tgrogers/accel-sim/traces/1.1.0.trace.summary.txt
wget ftp://ftp.ecn.purdue.edu/tgrogers/accel-sim/traces/tesla-v100/1.1.0.latest/rodinia_2.0-ft.tgz
tar -xzvf rodinia_2.0-ft.tgz
```

`rodinia_2.0-ft.tgz` 压缩 **7.2 MB**（清单里最小的一套；cutlass 要 77 GB）。  
解压后整套在原工程：

`accel-sim-framework-dev/.../hw_run/tesla-v100/rodinia_2.0-ft/`

GPGPU-Sim 后端（gitignore，README 说 setup 时 clone）：

```bash
git clone -b dev https://github.com/accel-sim/gpgpu-sim_distribution.git \
  gpu-simulator/gpgpu-sim
```

这里才有 `gpgpusim.config` 和 `accelwattch_*.xml`。

## 输入是什么（对应 gem5 的「配置 + 任务」）

| 角色 | 文件 | 说明 |
|---|---|---|
| 任务（SASS 轨迹） | `traces/backprop-rodinia-2.0-ft/kernelslist.g` | 剧本：Memcpy + kernel-1/2.traceg |
| 任务本体 | `kernel-1.traceg`（13MB）、`kernel-2.traceg`（4.7MB） | 真机 NVBit 采的 Volta SASS |
| 机器：性能 | `configs/SM7_QV100/gpgpusim.config` | Quadro V100 |
| 机器：轨迹延迟 | `configs/SM7_QV100/trace.config` | 来自 Accel-Sim 前端 |
| 机器：功耗 | `configs/SM7_QV100/accelwattch_sass_sim.xml` | AccelWattch 标定 XML |

`kernelslist.g` 内容：

```
MemcpyHtoD,...
kernel-1.traceg
MemcpyHtoD,...
kernel-2.traceg
```

kernel 头（真实，不是假案例）：

```
-kernel name = _Z22bpnn_layerforward_CUDAPfS_S_S_ii
-grid dim = (1,256,1)
-block dim = (16,16,1)
-binary version = 70
-nvbit version = 1.4
```

## 本机还缺什么才能「跑出输出」

- CUDA toolkit（`CUDA_INSTALL_PATH`，cmake 要 `cuda.h`）
- 编出 `gpu-simulator/bin/release/accel-sim.out`

有 CUDA 之后的命令见同目录 `若本机有CUDA如何跑.txt`。

输出（现在没有）会是仿真目录里的周期统计，打开功耗后还有 `accelwattch_power_report.log`（整卡平均瓦特 + INT/FPU/L2/DRAM 等，**没有 XY**）。

## 备份目录

```
成功运行备份/Accel-Sim/官方下载/
  1.1.0.trace.summary.txt     官方轨迹目录
  rodinia_2.0-ft.tgz          原包
  traces/backprop-.../        一个完整 app 的 kernelslist + two kernels
  configs/SM7_QV100/          性能+功耗+trace 配置
  configs/SM7_GV100/          功耗标定常用的 GV100
```

完整 10 个 rodinia 小应用仍在原工程 `hw_run/tesla-v100/`，备份只拷了 backprop 方便看格式。
