# 功耗 → 初始布局 → HotSpot

在 `/home/czh/任务仿真产生功耗文件工具` 下操作。

## 已拍板

**初始布局策略 = L-BFGS**（代码里叫 `lbfgs`，对应你写的 L-FBGS）。

依据：2026-10-19 在 16 组「固定 CPU / GPU / CIM 三颗」上对比 bstree / mcts / lbfgs，只看这一轮 HotSpot Tmax 时 lbfgs 15/16 独赢。常量写在 `src/defaults.py`。

bstree / mcts 仍留着，只给对照实验（`sweep_combos.py`），新跑默认不再走它们。

每个案例固定三颗：

| 颗 | 来源 |
|---|---|
| CPU | gem5 + McPAT（一种 CPU，换 gem5 配置） |
| GPU | Accel-Sim AccelWattch（目前只有 backprop；面积暂用 GV100 815 mm²） |
| CIM | AccelForge（换宏） |

## 怎么跑

默认（只 L-BFGS）：

```bash
python3 /home/czh/任务仿真产生功耗文件工具/热仿真流水线/src/run_pipeline.py
```

三种算法对照（已经跑过，不必当日常入口）：

```bash
python3 /home/czh/任务仿真产生功耗文件工具/热仿真流水线/src/sweep_combos.py
```

HotSpot：`任务仿真产生功耗文件工具/HotSpot-master/hotspot`
工程全貌与提速结论见 [SUMMARY.md](SUMMARY.md)。
