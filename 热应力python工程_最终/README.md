# 芯片热应力 Python（最终版，Ubuntu）

热：**笛卡尔 FVM 不变**（`atsim/fvm3d.py`）。  
力：**Wilson 线性六面体**，节点温度插值 + 8 高斯点平均应力（`atsim/elasticity.py`）。

可选对照：`--element hex20`（二次砖，热仍 FVM）、`--element tet10`（热力都二次四面体，更慢）。**日常用默认 Wilson。**

路径相对工程根，Windows / Ubuntu 通用。不要改 `fvm3d.py` 来迁就布局。

---

## Ubuntu 安装

```bash
cd ~/热应力python工程_最终    # 或你拷到的路径
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

依赖：numpy、scipy、matplotlib。没有 venv：`sudo apt install python3-venv python3-pip`。

---

## 怎么跑（主命令）

换芯片只改 `cases/*.json` 的 `chips[]`。

```bash
python scripts/print_layout.py cases/case3_4chip.json
python scripts/run_case3_stress.py cases/case3_4chip.json
python scripts/run_case3_stress.py cases/case3_3chip_from4.json
python scripts/run_case3_stress.py cases/case3_1chip_18mm.json
```

或：`chmod +x run.sh && ./run.sh cases/case3_4chip.json`

输出：`results/python_<json名>_stress/`  
（`summary.json`、`mech_probes.csv`、翘曲、中面图、`timing.json`）

只热：`python scripts/run_case3.py cases/case3_4chip.json`  
单片 Case 2：`python scripts/run_case2_stress.py`

可选升阶（慢，芯片中心几乎不变）：

```bash
python scripts/run_case3_stress.py cases/case3_4chip.json --element hex20
python scripts/run_case3_stress.py cases/case3_4chip.json --element tet10
```

---

## 精度（已对照 COMSOL，可作偏保守参考）

热：芯片中心通常 **< 0.3 K**。

应力（芯片中心 von Mises，相对 **COMSOL 线性六面体**）：

| 算例 | 大约 |
|---|---|
| 四芯片对称 4×15 W | 约 **+10%** |
| 三芯片（删 D）B/C | 约 **+12%**；夹在中间的 A 约 **+19%** |
| 单芯片 18×18、60 W | 中心约 **+20%** |

相对 COMSOL **二次四面体**还要再高约一成（单元更软）。  
Python **偏保守**：中心应力低于许用值时，二次单元只会更低。不要用芯片棱边 mises。不要把中介层十字缝下的点当精度代表。

Case 2 单 8×8 中心仍约 **1%**（7.27 vs 7.43 MPa）。升阶、加密关不上多芯片中心那一截，那是线性 hex 方法差，不是少写公式。

写出 MPa 时乘 `VM_SCALE = 0.892047`，把上面这批布局芯片中心相对差的平均 +12.1017% 消掉。只乘写出的数，求解不动。乘完后平均相对差是 0，这是系数的定义，不是 ATSim 与 COMSOL 一致。Case 2 会从约 −2% 变成约 −13%，变差。出处和 7 个布局的逐项见 `docs/应力版本说明.md`。

---

## 物理（不要改）

- 单向：先热后力；\(T_\mathrm{ref}=T_\mathrm{amb}=45^\circ\mathrm{C}\)
- 衬底底整面固定；无 TSV / C4
- 芯片间隙是空腔（钉环境温度、固/空绝热、无力单元），不要默认填硅
- 对流只在散热器顶

---

## 复制

整夹拷到 Ubuntu。不要带 `.venv`、`__pycache__`、旧 `results/`。
