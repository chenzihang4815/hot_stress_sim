# 又跑通的 4 个官方 ARM hello（同一程序，不同「机器」）

任务全部是仓库自带：

`tests/test-progs/hello/bin/arm/linux/hello`

仿真器：

`/home/czh/任务仿真产生功耗文件工具/gem5-stable/gem5-stable/build/ARM/gem5.opt`

都要先：

```bash
export LD_LIBRARY_PATH="$HOME/miniconda3/lib:$LD_LIBRARY_PATH"
cd /home/czh/任务仿真产生功耗文件工具/gem5-stable/gem5-stable
```

之前已经跑过的第 0 个（带 L1/L2）在上一级 `../output/`。

## 命令和结果

### 0）se.py + Timing + cache（上次）

```bash
./build/ARM/gem5.opt -d m5out \
  configs/deprecated/example/se.py \
  --cpu-type TimingSimpleCPU --caches --l2cache \
  --cmd tests/test-progs/hello/bin/arm/linux/hello
```

Hello world。约 4994 指令、69599 周期，**有** icache/dcache/l2 计数。

### 1）同一 se.py，Timing，但不要 cache

```bash
./build/ARM/gem5.opt -d 多案例/01_se_timing_无cache \
  configs/deprecated/example/se.py \
  --cpu-type TimingSimpleCPU \
  --cmd tests/test-progs/hello/bin/arm/linux/hello
```

Hello world。指令仍约 4994，但周期变成 **667118**（没 cache，访存更慢）。stats 里 **没有** icache。

### 2）se.py + Atomic（最快、时序最粗）

```bash
./build/ARM/gem5.opt -d 多案例/02_se_atomic_无cache \
  configs/deprecated/example/se.py \
  --cpu-type AtomicSimpleCPU \
  --cmd tests/test-progs/hello/bin/arm/linux/hello
```

Hello world。周期只有 **5841**。Atomic 几乎不模拟访存延迟。

### 3）换配置脚本：ARM starter_se.py

```bash
./build/ARM/gem5.opt -d 多案例/03_starter_se_atomic \
  configs/example/arm/starter_se.py \
  --cpu atomic -- \
  tests/test-progs/hello/bin/arm/linux/hello
```

注意：程序写在 `--` 后面，不是 `--cmd`。`--cpu atomic` 是 starter 自己的四档之一。Hello world。

### 4）再换脚本：教程 simple-arm.py（程序写死在脚本里）

```bash
./build/ARM/gem5.opt -d 多案例/04_simple-arm教程 \
  configs/learning_gem5/part1/simple-arm.py
```

没有 `--cmd`。脚本内部指向同一个 ARM hello。Hello world。无 cache。

## 对照（同一 hello）

| 案例 | 配置脚本 | 核 | cache | 指令 | 周期 |
|---|---|---|---|---|---|
| 0 | se.py | Timing | L1+L2 | ~4994 | 69599 |
| 1 | se.py | Timing | 无 | ~4994 | 667118 |
| 2 | se.py | Atomic | 无 | ~5020 | 5841 |
| 3 | starter_se.py | atomic | 无（atomic 档） | ~5020 | 字段名不同 |
| 4 | simple-arm.py | Timing | 无 | ~5914 | 388138 |

要看的文件：每个子目录里的 `command.txt`、`sim.stdout`、`stats.txt`。
