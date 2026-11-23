# gem5 备份目录说明

## 已真跑通的一次（ARM hello + se.py）

见 `input/`、`scripts/`、`output/`。说明仍见下方「那一次」。

## 官方输入全集（配置 + 仓库内程序）

见 **`官方输入清单/`**：

- `configs/`：原工程 159 个 Python 配置（含入口脚本和辅助库）
- `test-progs/`：仓库里自带的 hello / threads 二进制和源码
- `README.md`：分类、数量、哪份能在本机 ARM gem5 上跑
- `configs_py_清单.txt`：159 个 `.py` 的相对路径

原工程根：`/home/czh/任务仿真产生功耗文件工具/gem5-stable/gem5-stable`

---

# 已跑通：官方 ARM hello（2026-10-07）

仿真器：原工程 `build/ARM/gem5.opt`  
配置：`configs/deprecated/example/se.py`  
程序：`tests/test-progs/hello/bin/arm/linux/hello`（源码 `hello.c`）

```bash
./build/ARM/gem5.opt -d output \
    configs/deprecated/example/se.py \
    --cpu-type TimingSimpleCPU --caches --l2cache \
    --cmd tests/test-progs/hello/bin/arm/linux/hello
```

终端出现 `Hello world!`。`output/stats.txt` 为真活动计数。复跑：`scripts/run_arm_hello.sh`。
