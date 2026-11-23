# gem5 官方输入归纳（配置脚本 + 被仿真程序）

原工程：`/home/czh/任务仿真产生功耗文件工具/gem5-stable/gem5-stable`  
本目录是原样拷贝，方便离线翻，**不是**再编一份 gem5。跑仿真仍要用原工程里的 `build/ARM/gem5.opt`（脚本会 `import` 同目录的 `common/` 等）。

gem5 一次运行需要两类输入：

```
配置脚本  = 机器长什么样（几个核、有没有 cache）
被仿真程序 = 任务（ARM/x86/... 的 ELF）
```

---

## 一、配置脚本一共多少、在哪

原路径：`configs/`  
本备份：`官方输入清单/configs/`（已去掉 `__pycache__`）

**Python 文件共 159 个**（完整相对路径见 `configs_py_清单.txt`）。

不要把 159 都理解成「官方案例」。大约分成：

| 角色 | 大约数量 | 是什么 |
|---|---|---|
| **可直接丢给 gem5.opt 的入口脚本** | 约 50～70 | 里面会 `simulate()` / 搭 System 再跑 |
| **辅助库（不能单独当唯一输入）** | 其余 | `common/`、`ruby/`、`topologies/` 等，被入口脚本 import |

按文件夹：

| 备份内路径 | 原工程路径 | py 个数（约） | 用途 |
|---|---|---|---|
| `configs/example/` | 同上 | 78 | **主示例区**：stdlib、ARM starter、SE/FS、Ruby 测试 |
| `configs/deprecated/` | 同上 | 5 | 旧示例；**真正能用的 se.py 在这里** |
| `configs/learning_gem5/` | 同上 | 14 | Learning gem5 教程 |
| `configs/common/` | 同上 | 25 | 选项解析、cache 类、CPU 列表（库） |
| `configs/ruby/` | 同上 | 17 | 一致性协议配置（库） |
| `configs/topologies/` | 同上 | 10 | 片上网络拓扑（库） |
| `configs/dram/` `nvm/` `network/` `splash2/` `dist/` | 同上 | 若干 | DRAM/NVM 扫描、Splash-2、分布式 |
| `configs/boot/` | 同上 | 无 py，若干 `.rcS` | 全系统启动后跑的 guest 脚本 |

`configs/example/se.py` **不能跑**：打开只有一行 fatal，让你改用 deprecated 那份。

---

## 二、建议当作「官方入口案例」的脚本

下面这些才是「机器图纸」。ISA 列的是脚本面向的指令集；本机目前只有 `build/ARM/gem5.opt`，**只有 ARM 行能直接跑**。

### A. 用户态 SE（只要 ELF，不要内核）— 最接近已跑通的 hello

| 脚本（相对 configs/） | ISA | 程序从哪来 | 本机 |
|---|---|---|---|
| `deprecated/example/se.py` | 由 `--cpu-type` 和 gem5 二进制决定 | `--cmd` 指定。仓库 ARM hello 已有 | **已跑通** |
| `example/arm/starter_se.py` | ARM | 命令行里的二进制 | 可跑 |
| `example/gem5_library/arm-hello.py` | ARM | 联网下 `arm-hello64-static` | 可跑，无 cache |
| `example/arm/etrace_se.py` | ARM | 弹性指令轨迹 | 要额外 trace |
| `learning_gem5/part1/simple-arm.py` | ARM | 脚本内指定 binary | 可跑，无 cache |
| `learning_gem5/part1/two_level.py` | 教程默认偏 x86 | 命令行 binary | 需改 CPU 才配 ARM gem5 |
| `example/gem5_library/riscvmatched-hello.py` | RISC-V | Resources | 要 `build/RISCV` |
| `example/gem5_library/power-hello.py` | POWER | Resources | 要 POWER 包 |
| `example/gem5_library/fdp-hello.py` | 看脚本 | Resources | 看 ISA |

`se.py` 的开关（机器侧）：`--cpu-type`、`--caches`、`--l2cache`、`--l1d_size`、`--num-cpus`、`--cpu-clock`、`--mem-size` …  
任务侧：`--cmd`。

`starter_se.py` 的 CPU 四档：`atomic` / `minor` / `hpi` / `o3`（不是 TimingSimpleCPU 那张长表）。

### B. 全系统 FS（要内核 + 磁盘，Resources 下载）

| 脚本 | ISA | 本机 |
|---|---|---|
| `example/gem5_library/arm-ubuntu-run.py` | ARM | 镜像大、内存紧，不建议当下一例 |
| `example/gem5_library/arm-ubuntu-run-with-kvm.py` | ARM | 还要 KVM |
| `example/arm/starter_fs.py`、`fs_bigLITTLE.py`、`baremetal.py` | ARM | 要 firmware/kernel |
| `example/fs.py` | 通用旧 FS | deprecated 思路 |
| `x86-ubuntu-run.py` 等 | x86 | 要 X86 gem5 |
| `riscv-ubuntu-run.py`、`riscvmatched-fs.py` | RISC-V | 要 RISCV gem5 |

### C. 基准 / GPU / 其它（几乎都要 X86 或额外资源）

`x86-npb-benchmarks.py`、`x86-spec-cpu2017-benchmarks.py`、`x86-gapbs-benchmarks.py`、`x86-mi300x-gpu.py`、`multisim-*.py`、`dram/sweep.py`、`example/garnet_synth_traffic.py` 等。

完整入口文件名可用：

```bash
grep -E "entry" 本 README 不够时，看 configs_py_清单.txt
```

或打开备份里的 `configs/example/gem5_library/` 逐个看文件头注释。

---

## 三、被仿真程序（仓库里真正带了二进制的）

原路径：`tests/test-progs/`  
本备份：`官方输入清单/test-progs/`

**git 里编好的 ELF 只有 8 个**，外加源码和 Makefile。大基准、Ubuntu 镜像 **不在仓库里**。

### 二进制

| 备份内路径 | 格式 | 配哪份 gem5 |
|---|---|---|
| `test-progs/hello/bin/arm/linux/hello` | ARM 32 位静态 | **本机这份 ARM gem5.opt**（已跑通） |
| `test-progs/hello/bin/x86/linux/hello` | x86-64 静态 | 要 `build/X86` |
| `test-progs/hello/bin/x86/linux/hello32` | i386 静态 | 要 X86 |
| `test-progs/hello/bin/riscv/linux/hello` | RISC-V 64 静态 | 要 RISCV |
| `test-progs/hello/bin/power/linux/hello32` | PowerPC 32 | 要 POWER |
| `test-progs/hello/bin/power/linux/hello64le` | PPC64 LE | 要 POWER |
| `test-progs/hello/bin/power/linux/hello64be` | PPC64 BE | 要 POWER |
| `test-progs/threads/bin/x86/linux/threads` | x86-64 **动态**链接 | 要 X86；SE 动态库更麻烦 |

### 源码（官方准备了、二进制未放进 git 的）

| 路径 | 说明 |
|---|---|
| `hello/src/hello.c` | `printf("Hello world!\n");` |
| `hello/src/Makefile.arm` 等 | 用 **dockcross Docker** 交叉编译，仓库 **不含** ARM gcc |
| `m5-exit/src/m5-exit.c` | 调用 gem5 的 m5_exit，**没有现成 bin/** |
| `threads/src/threads.cpp` | 多线程示例 |

### 官方准备了、但不在 git 里的程序

- `arm-hello.py` 用的 `arm-hello64-static` → [resources.gem5.org](https://resources.gem5.org)  
- Ubuntu 内核 / 磁盘  
- NPB、SPEC、GAPBS 等  

---

## 四、和「已跑通那一次」的关系

真正跑过的输入仍在上一级，不要和这份总清单混淆：

```
成功运行备份/gem5/
  input/hello + hello.c     ← 那一次的任务拷贝
  scripts/se.py             ← 那一次的配置拷贝
  output/stats.txt          ← 那一次的输出
  官方输入清单/             ← 官方 configs + test-progs 全集（本目录）
```

本机当前能完整离线复现的组合只有：

```
build/ARM/gem5.opt
+ configs/deprecated/example/se.py
+ tests/test-progs/hello/bin/arm/linux/hello
```

对应备份：`../scripts/run_arm_hello.sh`。

---

## 五、怎么用这份备份观察

- **看有哪些入口**：翻 `configs/example/gem5_library/`、`configs/example/arm/`、`configs/deprecated/example/`、`configs/learning_gem5/`。  
- **看辅助库**：`configs/common/Options.py`（`--cpu-type` / `--caches` 定义处）。  
- **看任务**：`test-progs/hello/`。  
- **不要**在本备份目录里直接 `gem5.opt 官方输入清单/configs/...`：相对 `import common` 会找不到，除非把 `PYTHONPATH` 指回原工程 `configs/`。观察用备份，执行为：

```bash
cd /home/czh/任务仿真产生功耗文件工具/gem5-stable/gem5-stable
./build/ARM/gem5.opt configs/deprecated/example/se.py --help
```
