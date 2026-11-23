# 多芯片布局入口（Case 3）

求解算法仍是 Case 2 的 3D FVM，**不要改 `atsim/fvm3d.py`**。  
换方案只改 `chips[]` JSON，不要改 `atsim/fvm3d.py`。对照表：

| JSON | 用途 |
|---|---|
| `case3_cpu_gpu_cim.json` | **CPU + GPU + CIM 三颗**（L-BFGS 布局，默认入口） |
| `cpu_gpu_cim_lbfgs/*.json` | 同上三颗，16 组 gem5 CPU × AccelForge CIM（GPU 都是 backprop） |
| `case3_chips.json` | 3 芯片基线 CPU/GPU/HBM 25/40/8 W（旧，第三颗是 HBM 不是 CIM） |
| `case3_chips_hot.json` | 同一布局，功率 ×2 |
| `case3_chips_fine.json` | 同一布局，平面 0.25 mm（慢，仅网格试验） |
| `case3_2chip.json` | 两颗，少一条自由边 |
| `case3_4chip.json` | 2×2 对称 15 W |
| `case3_6chip.json` | GPU+CPU+四颗 HBM |
| `case3_3chip_spread.json` | 三芯片拉大间距 |
| `case3_3chip_tim_pads.json` | 三芯片，TIM 每颗一块垫（COMSOL 更稳） |
| `case3_10chip.json` | 十芯粒 |

## 你需要填的

每个芯片 6 个数：名字、左下角 (x,y)、长、宽、厚度、功耗。单位默认 **mm / μm / W**。

```json
{
  "name": "case3_3chip",
  "chips": [
    {"name": "CPU", "x_mm": 8,  "y_mm": 8,  "w_mm": 10, "h_mm": 10, "t_um": 150, "power_W": 25},
    {"name": "GPU", "x_mm": 22, "y_mm": 8,  "w_mm": 12, "h_mm": 10, "t_um": 150, "power_W": 40},
    {"name": "HBM", "x_mm": 8,  "y_mm": 22, "w_mm": 8,  "h_mm": 8,  "t_um": 150, "power_W": 8}
  ]
}
```

坐标系：`x_mm/y_mm` 是布局相对坐标，**不是**散热器左下角。程序会：

1. 算芯片包围盒
2. **衬底 / 中介层 / TIM 平面 = 该包围盒**（同长同宽）
3. 均热板 = 包围盒 + 一圈 `spreader_margin_mm`（默认 11 mm）
4. 散热器 = 均热板 + 一圈 `sink_margin_mm`（默认 15 mm）
5. 把整簇芯片平移到封装正中

芯片 **不重叠**。换布局时不用手改衬底/中介层/TIM 尺寸。

可选 `package`（不写则用下面默认）：

- `tim_um`、`t_spreader_mm`、`t_sink_mm`
- `t_substrate_mm`（默认 1.0）、`t_interposer_um`（默认 100）
- `r_convec_K_per_W`（默认 0.1）、`ambient_C`（默认 45）
- `fill_die` 默认 false；`tim_cover` 默认 `"chips"`（TIM 平面 = 芯片包围盒）

## 自动生成

```
衬底（有机，平面 = 芯片包围盒，底面力学固定）
  → 中介层（整块硅，无 TSV，平面同包围盒）
    → 各芯片（间隙空气，不填硅）
      → TIM（平面同包围盒）
        → 均热板 → 散热器
```

若 TIM 要铺满均热板：`"tim_cover": "spreader"`。  
若只要每颗芯片顶一块小垫：`"tim_cover": "none"`。

## 运行

热：

```bat
python scripts/print_layout.py cases/case3_4chip.json
python scripts/run_case3.py
python scripts/run_case3.py cases/case3_4chip.json
```

应力（同一份 `chips[]` JSON，单向：先热后力）：

```bat
python scripts/run_case3_stress.py
python scripts/run_case3_stress.py cases/case3_10chip.json
```

热输出：`results/python_case3_chips/`  
应力输出：`results/python_case3_chips_stress/`（`mech_probes.csv`、`summary.json`、`timing.json`、中面图、**封装翘曲** `warpage_interposer.png` / `warpage_interposer.csv`）

芯片内部 von Mises 仍可看可靠性；**布局更该看封装级翘曲**：中介层中面 \(w=u_z\) 的峰谷（peak-to-valley），对应 ATMPlace 的整包弯曲，不是芯片内部应力。提取在 `atsim/warpage.py`（不新解薄板方程，从已有 3D 弹性 `uz` 取样）。衬底底固定时 \(u_z\) 已含约束，峰谷才是弯曲量。棱边 mises max 不当判据。

`atsim/fvm3d.py` 和 `atsim/elasticity.py` 都不用改。约束默认衬底底固定。

和 COMSOL 比是可选脚本，**不要**放进上面两条主命令。有 CSV 时再单独跑：

```bat
python scripts/compare_case3_stress.py path/to/die_points.csv
python scripts/compare_case3_interposer.py path/to/interposer.csv
```

公式、误差、如何换芯片方案见 `05_工程成果与公式.md`。

