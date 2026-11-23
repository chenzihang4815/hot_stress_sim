# Python 热仿真工具说明

对应王启盼论文第 2 章：**芯片级稳态热仿真（ATSim3D）**。  
当前工程实现的是其中两块：**热阻网络（粗粒度）** 和 **热传导方程的有限体积离散（细粒度）**。  
第 2 章的 Kirchhoff 变换、全局–局部分层并行、TSV 子域，代码里有骨架，Case 1–3 的主路径还没用上。

根目录：`E:\给claude操作文件临时存放\芯片热仿真comsol和hotspot和代码`

---

## 1. 论文第 2 章 vs 本工程

论文稳态方程：

\[
\nabla\cdot\big(\kappa(\boldsymbol r,T)\nabla T\big)+P(\boldsymbol r,T)=0
\]

| 论文 ATSim3D | 本工程现在实际跑的 |
|---|---|
| 全局粗网格热阻网络，给局部 Dirichlet / \(T_g\) | `atsim/global_local.py` 的 `solve_global_rc`。Case 1 会算一遍作对照，**封装 Case 2/3 主求解不用它** |
| 改进 Kirchhoff，把 \(\kappa(T)\) 变成对 \(\psi\) 的线性泊松 | `atsim/materials.py` 有正/逆变换。JSON 里 `solver.nonlinear=false`，**主路径未启用** |
| 局部子域 FVM 解 \(\psi\)，外层 Picard 更新漏电 | **主路径直接对 \(T\) 做 3D FVM**（`atsim/fvm3d.py`），常数 \(\kappa\)，一次线性求解 |
| Robin 对流放在全局层 | Case 1：硅顶面；Case 2/3：散热器顶面。都在 FVM 里直接施加 |

一句话：算法家族是论文第 2 章的热阻网络 + 热传导 FVM；工程上为了和 COMSOL/HotSpot 对齐封装几何，主求解器是 **整包连续体 FVM**，不是论文里的「粗网格 + 多子域 Kirchhoff」。

---

## 2. 热传导方程在哪个文件

**唯一组装并求解 \(\nabla\cdot(\kappa\nabla T)+P=0\) 的文件：`atsim/fvm3d.py`**

| 函数 | 行号量级 | 做什么 |
|---|---|---|
| `_assemble` | ~40 行起 | 单元中心 FVM：六个面热流 + 热源体积积分 + 顶面 Robin，组装 \(G T = b\) |
| `_solve_linear` | 文件后部 | `scipy.sparse.linalg.spsolve`（CG 失败则回退直接法） |
| `solve_fvm3d` | 文件末尾 | 建网格、铺功耗、铺 \(\kappa\)、调用组装与求解 |

界面热导（论文 FVM「调和平均」的非均匀网格写法）在 `_assemble` 里：

\[
R=\frac{d_p}{2\kappa_p A}+\frac{d_n}{2\kappa_n A},\quad g=1/R
\]

不要和「把 \(\kappa\) 调和平均后再除中心距」搞混；TIM（薄、k=4）接到铜（厚、k=400）必须用串联热阻，否则热点会虚高约 8 K。

**热阻网络（论文全局粗网格）在：`atsim/global_local.py` → `solve_global_rc`**

把每个功耗块当成一个节点：块顶面对流 \(g=hA\)，相邻块侧向 \(g_{ij}=\kappa\cdot L_\mathrm{shared}\cdot t / L_{ij}\)。Case 1 用它估块温；Case 3 主路径不走这里。

Kirchhoff：`atsim/materials.py`（`kirchhoff_forward` / `kirchhoff_inverse`）。

---

## 3. 工程结构（数据流）

```
cases/*.json                    规格（Case 3 只填 chips[]）
        │
        ▼
atsim/geometry.py               load_case
        │  若有 chips[] → atsim/layout.py 展开封装/网格/探针
        ▼
atsim/package_grid.py           叠层 + 空腔 mask + 功耗只写在芯片硅体积
        ▼
atsim/fvm3d.py                  热传导方程 FVM  ← 计算核心
        ▼
atsim/visualize.py              中面热图、探针
atsim/hotspot_io.py             写出 HotSpot .flp/.ptrace/.config
```

| 文件 | 角色 | 是不是「算热」 |
|---|---|---|
| `atsim/fvm3d.py` | 热传导方程离散 + 求解 | **是，主求解器** |
| `atsim/global_local.py` | 论文热阻网络 + Kirchhoff 局部框架 | 辅；Case 1 对照 |
| `atsim/materials.py` | \(\kappa(T)\)、Kirchhoff 积分 | 物性，不求解场 |
| `atsim/package_grid.py` | 几何网格、固体/空腔、功耗密度 | 前处理 |
| `atsim/layout.py` | `chips[]` → 完整封装 Case | 前处理，Case 3 入口 |
| `atsim/geometry.py` | JSON → Case 对象 | 前处理 |
| `atsim/hotspot_io.py` | 导出 HotSpot | 后处理 |
| `atsim/visualize.py` | 图和 CSV | 后处理 |
| `atsim/elasticity.py` | 单向热弹性（hex，`div σ=0`） | 应力 |
| `atsim/warpage.py` | 中介层中面 \(w=u_z\) 峰谷 | 封装级翘曲，布局指标 |
| `atsim/compare.py` | 和 COMSOL 点云比 | 后处理 |
| `scripts/run_case1.py` | 仅芯片 | 入口 |
| `scripts/run_case2.py` | 单硅片 + 封装 | 入口 |
| `scripts/run_case3.py` | 多芯片，只改 JSON | 入口 |
| `scripts/run_case3_stress.py` | 多芯片热+应力+翘曲 | 入口 |

---

## 4. 三个算例怎么用

| | 规格 | 入口 | 输出 |
|---|---|---|---|
| Case 1 仅芯片 | `cases/case1_2d_ic.json` | `python scripts/run_case1.py` | `results/python/` |
| Case 2 单片+封装 | `cases/case2_package.json` | `python scripts/run_case2.py` | `results/python_case2/` |
| Case 3 多芯片热 | `cases/case3_chips.json` | `python scripts/run_case3.py` | `results/python_case3_chips/` |
| Case 3 多芯片应力 | 同一 JSON | `python scripts/run_case3_stress.py` | `results/python_case3_chips_stress/` |

Case 3 你只改 `chips[]`：名字、`x_mm,y_mm,w_mm,h_mm,t_um,power_W`。  
封装、网格、探针由 `layout.py` 生成。说明见 `cases/README_chips.md`。

芯片间隙默认是 **空腔**（空气 / 互连），中面热图只画芯片。不要填硅。TIM / 衬底 / 中介层平面默认等于芯片包围盒。应力：`python scripts/run_case3_stress.py`，同一 `chips[]`。

Case 2/3 叠层默认 2.5D：有机衬底 → 硅中介层（无 TSV）→ 芯片 → TIM → 均热板 → 散热器。中介层下无 C4，界面按理想粘结。力学默认固定衬底底面。

---

## 5. 物理边界（Case 2/3 封装）

自下而上：有机衬底 → 硅中介层（无 TSV）→ 硅芯片（发热，间隙空气）→ TIM（k=4，盖满均热板底）→ 铜均热板 → 铜散热器。  
对流只加在散热器顶：\(h=1/(R_\mathrm{convec} A_\mathrm{sink})\)，默认 \(R=0.1\,\mathrm{K/W}\)。侧面和底面绝热。应力约束在衬底底面。

---

## 6. 和 COMSOL / HotSpot

三家共用同一套几何和功耗。Python 与 COMSOL 比的是连续体场；HotSpot 报块平均，峰值会略低。  
Case 1 硅面已对齐（误差 < 1 K）。Case 2/3 已换成 2.5D，旧 COMSOL CSV 对不上，需按新手册重建。HotSpot 没有衬底/中介层，只作趋势。
