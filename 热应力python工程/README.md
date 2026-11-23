# 芯片热应力 Python 工程（Ubuntu 可独立运行）

从原工作区抽出的 **热 + 单向热弹性 + 封装翘曲** 代码。不含 COMSOL `.mph`、HotSpot 源码、计算结果。

路径全部相对工程根目录，Windows / Ubuntu 都能跑。

---

## 目录

```
热应力python工程/
  atsim/                 求解器（不要为了改布局去改 PDE）
    fvm3d.py             唯一热传导方程  ∇·(κ∇T)+P=0
    elasticity.py        单向热弹性  ∇·σ=0, σ=C:(ε-αΔT I)
    layout.py            chips[] → 2.5D 封装叠层
    warpage.py           中介层中面 uz 峰谷
  cases/                 算例 JSON（换芯片只改这里）
  scripts/               入口
  docs/                  公式与几何手册
  results/               运行后自动生成
  requirements.txt
```

主路径（泛化能力最好）：**只改 `cases/*.json` 的 `chips[]`，跑 `scripts/run_case3_stress.py`。**

---

## Ubuntu 安装

Python 3.8+（3.10 / 3.12 均可）。建议虚拟环境：

```bash
cd ~/热应力python工程          # 或你复制到的路径
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

依赖：`numpy`、`scipy`、`matplotlib`。

若系统没有 `python3-venv`：

```bash
sudo apt update
sudo apt install python3-venv python3-pip
```

---

## 怎么跑

先看布局（不求解，秒出）：

```bash
python scripts/print_layout.py cases/case3_4chip.json
```

多芯片热 + 应力 + 翘曲（主命令）：

```bash
python scripts/run_case3_stress.py cases/case3_chips.json    # 3 芯片基线
python scripts/run_case3_stress.py cases/case3_4chip.json    # 2×2
python scripts/run_case3_stress.py cases/你的.json
```

只热、不算应力：

```bash
python scripts/run_case3.py cases/case3_4chip.json
```

单片封装对照（Case 2，几何写在 JSON 的 stack 里，不是 chips[]）：

```bash
python scripts/run_case2_stress.py
python scripts/run_case2_stress.py cases/case2_package_hot.json
```

输出在 `results/python_<json名>_stress/`：

| 文件 | 内容 |
|---|---|
| `summary.json` | 探针、翘曲峰谷、计时 |
| `mech_probes.csv` | 关键点 T / uz / vm / sxx / syy |
| `warpage_interposer.csv` `.png` | 中介层中面 \(w=u_z\) |
| `overview_chips.png` | 芯片中面 T / von Mises / uz / \|u\| |
| `timing.json` | 热/力组装与求解时间 |

四芯片默认网格大约 1–3 分钟（视机器而定）。`case3_chips_fine.json`（0.25 mm）很慢，不要当日常算例。

---

## 换芯片方案

每个芯片 6 个数：名字、左下角 (x,y)、长、宽、厚度、功耗。单位 **mm / μm / W**。

```json
{
  "name": "my_layout",
  "chips": [
    {"name": "CPU", "x_mm": 8,  "y_mm": 8,  "w_mm": 10, "h_mm": 10, "t_um": 150, "power_W": 25},
    {"name": "GPU", "x_mm": 22, "y_mm": 8,  "w_mm": 12, "h_mm": 10, "t_um": 150, "power_W": 40}
  ]
}
```

程序会：芯片包围盒 → 衬底/中介层/TIM 同平面 → 均热板、散热器按边距放大 → 整簇居中。间隙是空气，不填硅。无 TSV / C4。

**不要改** `atsim/fvm3d.py` 来迁就布局。

现成 JSON：

| 文件 | 说明 |
|---|---|
| `case3_chips.json` | 3 芯片 25/40/8 W |
| `case3_2chip.json` | 两颗 |
| `case3_4chip.json` | 2×2，各 15 W |
| `case3_6chip.json` | GPU+CPU+四颗 HBM |
| `case3_10chip.json` | 十芯粒 |
| `case3_4chip_10mm.json` | 四颗 10×10 紧贴 |
| `case2_package.json` | 单 8×8 硅片封装 |

更细的字段见 `cases/README_chips.md`。

---

## 物理（主路径）

1. **热**：单元中心 3D FVM，串联面热阻，散热器顶 Robin，空腔钉在环境温度。
2. **力**：温度算完后单向加载，\(T_\mathrm{ref}=T_\mathrm{amb}=45^\circ\mathrm{C}\)。默认 Wilson/QM6 六面体。衬底底面固定。
3. **翘曲**：中介层中面 \(u_z\) 的峰谷（peak-to-valley），单位 μm。

公式见 `docs/05_工程成果与公式.md`。

---

## 可选：和 COMSOL 比

比较脚本**不进**主求解。有 CSV 再单独跑：

```bash
python scripts/compare_case3_stress.py path/to/die_points.csv
python scripts/compare_case3_interposer.py path/to/interposer.csv
```

本包不含 COMSOL 模型和点云。

---

## 复制到 Ubuntu

整夹拷过去即可，例如：

```bash
# 在 Windows 上打成 zip，或用 scp / U盘
scp -r 热应力python工程 user@ubuntu:~/
```

不要带 `__pycache__`、`.venv`、`results/`（空目录可以保留）。
