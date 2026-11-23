# 三个原工程试跑备份（2026-10-07，本机 Ubuntu）

本目录只放**原工程里真正跑过或明确缺件**的材料，不含「最简单案例」里的假 stats / 假 trace。

三条链互相独立：

```
CPU 性能：gem5        输入=配置脚本+程序二进制     输出=stats.txt（活动计数，不是功耗）
CPU 功耗：McPAT       输入=XML（架构 param + 活动 stat） 输出=分层面积 mm² + 功耗 W
GPU 仿真：Accel-Sim   输入=SASS 轨迹 + GPGPU-Sim 配置  输出=周期 + AccelWattch 功耗
```

不要把 gem5 的 `stats.txt` 喂给 Accel-Sim。

---

## 本机结论一览

| 工具 | 原工程路径 | 本机结果 |
|---|---|---|
| **McPAT 1.0** | `McPAT_v1.0.tar/McPAT_v1.0/mcpat` | **成功**。官方 `ARM_A9_2000.xml` 与官方 ExampleResults 数字一致。 |
| **gem5** | `gem5-stable/gem5-stable` | **未完成**。源码齐，ARM 目标编了一半；第一次因找不到 conda 的 `libpython3.11.so` 失败，已设 `LD_LIBRARY_PATH` 续编。GCC 15.2 超出官方支持范围（官方 GCC 11–14.2）。 |
| **Accel-Sim** | `accel-sim-framework-dev/accel-sim-framework-dev` | **无法跑**。前端源码齐，缺 GPGPU-Sim 后端、功耗 XML、轨迹、CUDA、可执行文件。 |

下面按工具写：已有 / 缺少 / 名词。

---

## 1. McPAT —— 已成功

### 做了什么

原版 `makefile` 写死 `g++ -m32`（32 位）。本机是 64 位 Ubuntu，没有 32 位工具链，所以用：

```bash
cd McPAT_v1.0.tar/McPAT_v1.0/mcpat
make -j8 CXX="g++" CC="gcc"
./mcpat -infile ProcessorDescriptionFiles/ARM_A9_2000.xml -print_level 5 -opt_for_clk 1
```

警告若干（空括号、缺 return），**链接成功**。跑官方 ARM Cortex-A9 例题，stdout 与 `ExampleResults/ARM_A9_2000` 一致。

另开了官方 `Penryn.xml`（Intel 酷睿类 OOO 大核）。McPAT 会对目标频率做设计空间搜索，Penryn 比 A9 慢一个数量级；跑了约 20 分钟仍停在 “is computing…”，为给 gem5 腾核已中止。**不是缺输入。** 需要时用备份里的 `input/Penryn.xml` 自己再跑。成功案例以 ARM_A9 为准。

### 备份里有什么

```
McPAT/
  bin/mcpat                          64 位可执行文件
  input/ARM_A9_2000.xml              官方输入
  input/Penryn.xml                   官方输入（更大）
  output/ARM_A9_2000.txt             本次跑出的完整功耗树
  output/ARM_A9_2000.official_expected.txt
  output/ARM_A9_2000_blocks.csv      从 txt 抽出的块级表（给热仿真用）
  scripts/run_mcpat.sh               复跑脚本
  scripts/read_mcpat_output.py       抽 Processor 级数字
```

复跑：

```bash
./McPAT/scripts/run_mcpat.sh
# 或指定 XML：
./McPAT/scripts/run_mcpat.sh /绝对路径/ARM_A9_2000.xml 5
```

### 输入长什么样

XML 分两类标签：

- `<param>`：**芯片长什么样**。核数、工艺节点 nm、主频 MHz、cache 容量、流水线宽度。热仿真里的「块的个数和面积」主要来自这些。
- `<stat>`：**这段程序有多忙**。周期数、提交指令数、cache 访问/缺失。同一块芯片、不同任务，改的是这些。

官方 A9 模板是 2 核、40 nm、2 GHz、in-order 风格参数（`machine_type` 等），stat 是一份固定的活动计数，不是 gem5 现跑出来的。

### 输出（Processor 级，已与官方对齐）

| 量 | 值 | 含义 |
|---|---|---|
| Area | 5.51463 mm² | 整处理器面积 |
| Peak Power | 1.77235 W | 架构最坏情况，可大于实际任务功耗 |
| Total Leakage | 0.109185 W | 漏电（亚阈 + 栅） |
| Peak Dynamic | 1.66316 W | 峰值动态 |
| Runtime Dynamic | **2.98062 W** | **这段 stat 对应的平均动态功耗**，热仿真稳态热源优先用这个 |
| Core clock | 2000 MHz | |
| Technology | 40 nm | |

分层（热仿真可贴 floorplan 的块，完整表见 csv）：

- Processor 5.51 mm²，Runtime Dynamic 2.98 W
- Total Cores（2 个）面积 4.96 mm²
- 单核 Core 2.48 mm²，Runtime Dynamic 1.08 W
- Instruction Cache / Data Cache / Execution Unit / NoC 等子块都有 Area + 四种功耗

注意：McPAT 打印树里，子块 Runtime Dynamic **加起来不必等于父块**（有共享/复制/目录项），热仿真用父块或自己选定的一层，不要把所有叶子再加一遍当芯片总功耗。

### 专业名词

- **McPAT**：HP Labs 的 CPU **架构级**功耗/面积/时序模型（MICRO 2009）。输入是微结构参数 + 性能计数器，不是门级网表。
- **CACTI**：McPAT 内嵌的 cache/SRAM/互连面积延迟模型。
- **Peak Power vs Runtime Dynamic**：Peak 是「所有单元按峰值开关」的上界；Runtime Dynamic 用 XML 里的访问计数，随任务变。
- **Leakage / 漏电**：哪怕不翻转也耗的静态电。McPAT 用 XML 的 `temperature`（开尔文，官方 A9 是 340 K ≈ 67 °C）算亚阈漏电和栅漏电。温度升高漏电升得很快，所以热闭环里可以只改 temperature 重跑 McPAT，不必重跑 gem5。
- **ITRS LOP device**：Low Operating Power 工艺档。XML `device_type=2`。
- **print_level 0–5**：5 打到 cache/ALU 这一层。
- **opt_for_clk**：1=尽量满足 XML 目标频率（可能换更大的电路）；0=只优化能量延迟积。

---

## 2. gem5 —— 源码齐，二进制未完成

### 已有

- 完整 gem5 源码：`gem5-stable/gem5-stable`
- 官方 ARM hello 配置：`configs/example/gem5_library/arm-hello.py`
- 仓库自带 ARM hello 二进制：`tests/test-progs/hello/bin/arm/linux/hello`（静态 ARM ELF）
- 经典 SE 脚本：`configs/deprecated/example/se.py`（可加 `--caches --l2cache`）
- `build/ARM/` 下已有约 1200+ 个 `.o`，说明有人/先前会话已经开始编 ARM 目标
- 本机：g++ 15.2、SCons 4.11、Python 3.11（conda）、m4、zlib、conda 的 protobuf 35.1（`pkg-config` 在 `~/miniconda3/lib/pkgconfig`）

### 缺 / 卡住

1. **还没有 `build/ARM/gem5.opt`**。续编命令：

   ```bash
   export PATH="$HOME/miniconda3/bin:$PATH"
   export LD_LIBRARY_PATH="$HOME/miniconda3/lib:$LD_LIBRARY_PATH"
   export PKG_CONFIG_PATH="$HOME/miniconda3/lib/pkgconfig:$PKG_CONFIG_PATH"
   cd gem5-stable/gem5-stable
   python3 -m SCons build/ARM/gem5.opt -j8
   ```

   第一次失败原因：configure 测 Python 时动态链接器找不到 `libpython3.11.so.1.0`（文件在 conda 的 `lib/`，默认 `LD_LIBRARY_PATH` 没有）。加上面第二行即可。

2. **GCC 15.2 不在官方支持列表**（gem5 写的是 GCC 11–14.2）。可能编过、也可能在后期 C++ 报错。若失败，需要系统装 GCC 13/14，或 conda 里装旧编译器。本机 **没有 sudo 免密**，不能直接 `apt install g++-14`。

3. **系统 apt 没有 `libprotobuf-dev`**，用的是 conda protobuf。必须带着 `PKG_CONFIG_PATH` 编，否则 scons 会说找不到 protobuf。

4. `arm-hello.py` 会从 gem5 Resources **网上拉** `arm-hello64-static`。若断网，改用仓库自带：

   ```bash
   ./build/ARM/gem5.opt configs/deprecated/example/se.py \
       --cpu-type TimingSimpleCPU --caches --l2cache \
       --cmd tests/test-progs/hello/bin/arm/linux/hello
   ```

   hello 极短，cache 计数很少，这是正常的。McPAT 要「有东西吃」需要更长的程序或至少带 cache 的配置。

### 编出来之后期望的输入 / 输出

- **输入**：Python 配置（CPU 类型、频率、cache）+ 客户程序（ARM 二进制）。
- **输出目录**（`-d m5out`）：
  - `stats.txt`：周期、指令、cache 命中/缺失。这是给 McPAT 的活动计数来源。
  - `config.ini`：实例化后的机器。
- gem5 **不算功耗**。仓库里的 `MathExprPowerModel` 是示意公式，论文级不要用。

### 专业名词

- **SE vs FS**：Syscall Emulation 只仿真用户态，系统调用由 gem5 代理，跑 hello 用这个。Full System 要内核和磁盘映像，重得多。
- **TimingSimpleCPU**：按时序、每拍大约一条指令的简单核，不是乱序 O3。O3 更慢、stats 字段也更多。
- **ISA 目标 `build/ARM`**：这个二进制里只链了 ARM 译码器。跑 x86 程序要 `build/X86` 或 `build/ALL`。
- **scons**：gem5 的构建系统（类似 make）。
- **stats.txt**：模拟结束时的计数器快照。McPAT 要的是 `numCycles`、`committedInsts`、`icache.overallAccesses` 这类，不是温度。

---

## 3. Accel-Sim —— 无法在本机跑

原仓库按设计就不包含仿真后端和轨迹。`.gitignore` 明确排除：

```
gpu-simulator/gpgpu-sim
hw_run/
gpu-app-collection/
```

### 已有（前端）

| 文件 | 作用 |
|---|---|
| `gpu-simulator/accel-sim.cc` | SASS 轨迹驱动的仿真封装 |
| `gpu-simulator/trace-parser/` | 解析 `kernelslist.g` + `.traceg` |
| `gpu-simulator/ISA_Def/*_opcode.h` | Kepler→Hopper 的 SASS 操作码 |
| `gpu-simulator/ISA_Def/accelwattch_component_mapping.h` | SASS 操作码 → 功耗单元（INT/FPU/…） |
| `gpu-simulator/configs/tested-cfgs/*/trace.config` | 轨迹侧延迟（QV100/GV100/H100…） |
| `util/job_launching/configs/define-standard-cfgs.yml` | 如何拼 `QV100-Accelwattch_SASS_SIM` |
| `get-accel-sim-traces.py` | 从 Purdue FTP 下官方预采集轨迹 |
| `gpu-simulator/setup_environment.sh` | **会 git clone GPGPU-Sim** |
| cmake / g++ / libxml2-dev / libzstd-dev / 一批 boost 头 | 本机编译器侧大体够 |

### 缺少（缺一就不能真跑）

| 缺什么 | 为什么致命 | 怎么补 |
|---|---|---|
| `gpu-simulator/gpgpu-sim/` | 周期级性能模型 + AccelWattch 引擎 + `gpgpusim.config` + `accelwattch_*.xml` 全在这里 | `source gpu-simulator/setup_environment.sh`（clone `gpgpu-sim_distribution` 的 `dev` 分支） |
| `accelwattch_sass_sim.xml` 等 | 功耗模型参数（电容、缩放） | 随 GPGPU-Sim 来，本树里 **0 个 xml** |
| `gpgpusim.config` | SM 数、L2、时钟域、存储 | 同样在 GPGPU-Sim 的 `configs/tested-cfgs/SM7_QV100/` 等 |
| `accel-sim.out` | 仿真器二进制 | clone 后 cmake 编译 |
| `kernelslist.g` + `kernel-*.traceg` | 要重放的 SASS 轨迹 | 有 GPU：NVBit 采集；没 GPU：`./get-accel-sim-traces.py` |
| CUDA toolkit（`nvcc`、`/usr/local/cuda`） | setup 脚本默认要 `CUDA_INSTALL_PATH`；编 tracer 更要 | 本机没有，也 **没有 NVIDIA GPU**（无 `nvidia-smi`） |
| `gpu-app-collection` | PTX 模式要的 CUDA 程序 | 另 clone；SASS 模式用现成 traces 可暂不需要 |

`define-standard-cfgs.yml` 里功耗配置是：

```
QV100-Accelwattch_SASS_SIM
  = $GPGPUSIM_ROOT/.../SM7_QV100/gpgpusim.config
  + -power_simulation_enabled 1
  + -power_simulation_mode 0          # 0=仿真器活动；1=硬件计数；2=混合
  + -accelwattch_xml_file accelwattch_sass_sim.xml
```

没有 GPGPU-Sim，这条路径是空的。

### 本机最短补齐（仍缺 CUDA）

```bash
cd accel-sim-framework-dev/accel-sim-framework-dev
export CUDA_INSTALL_PATH=/usr/local/cuda   # 本机目前不存在
source ./gpu-simulator/setup_environment.sh
cmake -S ./gpu-simulator/ -B ./gpu-simulator/build
cmake --build ./gpu-simulator/build -j$(nproc)
cmake --install ./gpu-simulator/build
./get-accel-sim-traces.py
```

功耗请用 **Volta QV100/GV100 + Accelwattch_SASS_SIM**。H100 的性能模型有，AccelWattch XML 是按 Volta 标定的，不要一上来 H100 开功耗。

没有 GPU 也可以：装 CUDA **toolkit**（不必有卡）+ 官方 traces。采自己的程序轨迹则必须有 NVIDIA 卡和 NVBit。

### 专业名词

- **SASS**：NVIDIA GPU 的机器指令（真机上跑的那种），不是 CUDA C，也不是 PTX。Accel-Sim 轨迹驱动模式重放的是 SASS。
- **PTX**：CUDA 的中间表示。GPGPU-Sim 还能功能仿真 PTX（不需要轨迹，但要程序）。精度/覆盖和 SASS 重放不同。
- **NVBit**：NVIDIA 的二进制插桩工具，在真机上把每条 SASS 记成 `.traceg`。
- **kernelslist.g**：一次运行的「剧本」：先 `MemcpyHtoD`，再列出每个 `kernel-*.traceg`，再 `MemcpyDtoH`。仿真器入口是这个文件。
- **GPGPU-Sim**：周期级 GPU 性能模型（SM 流水线、L1/L2、NoC、DRAM）。Accel-Sim 是它的轨迹前端。
- **AccelWattch**：挂在 GPGPU-Sim 上的 GPU 功耗模型（MICRO 2021），标定对象是 **Volta GV100**。输出 `accelwattch_power_report.log`：整卡平均瓦特 + INT/FPU/L2/DRAM 等功能单元瓦特。
- **没有 XY 坐标**：AccelWattch 不给 per-SM 的空间图。热仿真要自己做 floorplan，把功能单元功率摊到块上。
- **SM / cluster**：Streaming Multiprocessor，GPU 上的计算簇。改 `-gpgpu_n_clusters` 等于改 SM 数，块的个数会变。

---

## 和热应力工具的接口（跑通之后怎么用）

1. **现在就能用的**：`McPAT/output/ARM_A9_2000_blocks.csv` 和 `ARM_A9_2000.txt`。把选定层级的 Area、Runtime Dynamic 填进 ATSim 的 `chips[]`。这是架构块，不是 Innovus 宏。
2. **gem5 编出来之后**：真 `stats.txt` → 填 McPAT `<stat>` → 再出一套与任务相关的 P_i。面积（`<param>`）不变则 floorplan 可复用。
3. **Accel-Sim**：本机先不要指望。GPU 功耗链放到有 CUDA 的机器；标定用 QV100。

环的建议仍然是：布局循环只用 `{块, P, A}`；周期级仿真不要进内环。
