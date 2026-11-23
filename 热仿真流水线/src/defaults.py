"""已拍板的流水线默认值。

2026-10-19 在 16 组「固定 CPU/GPU/CIM 三颗」上对比 bstree / mcts / lbfgs。
只看这一轮 HotSpot Tmax：lbfgs 15/16 独赢。用户随后定为初始布局策略。
bstree / mcts 仍留在仓库里，只给对照实验用，新跑默认不再走它们。
"""

INITIAL_LAYOUT = "lbfgs"
INITIAL_LAYOUT_SEED = 13
