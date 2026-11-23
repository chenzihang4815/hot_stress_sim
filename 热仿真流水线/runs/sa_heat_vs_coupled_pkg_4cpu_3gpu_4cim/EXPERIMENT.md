# Heat-only vs heat+stress SA — pkg_4cpu_3gpu_4cim

- started (UTC): 2026-10-25T12:58:08.741975+00:00
- title: 4 CPU + 3 GPU + 4 CIM（11 颗）
- why: 夹在 10 与 13 之间
- n_chiplets: 11
- seed=22  N_STEPS=16  T0=2.0  alpha=0.95
- gaps: any 0.5 mm, GPU–GPU 2.0 mm
- coupled J = (Tmax-45)/25.0 + 1.0*vm_chip_max/100.0 + 0.5*w_pv/4.0
- vm = worst interior von Mises among normal chiplets (CPU/GPU, CIM area≥1 mm²)
- ATSim package overrides: {"pitch_die_mm": 2.0, "pitch_coarse_mm": 8.0}
- one scheme at a time (12 GB); CG only, no spsolve

## Run log

- start outline 59.12 x 44.48 mm  min_gap=0.56 mm  GPU-GPU=2.62 mm
- `heat_only` init  Tmax=70.550C  vm_chipmax=149.00 MPa (GPU_TITANX)  w_pv=4.221 um  J=1.0220  stress=138.8s
  `heat_only` step   1 accept rotate CPU_armL2  Tmax=70.550C  J=1.0220
  `heat_only` step   3 accept nudge GPU_3070 (-1.0,0.0)mm  Tmax=70.580C  J=1.0232
  `heat_only` step   6 accept swap CIM_wang<->CIM_colonnade  Tmax=70.580C  J=1.0232
  `heat_only` step   9 accept nudge GPU_3070 (-2.0,0.0)mm  Tmax=70.680C  J=1.0272
  `heat_only` step  10 accept nudge GPU_TITANX (0.0,0.5)mm  Tmax=70.730C  J=1.0292
  `heat_only` step  12 accept nudge CPU_armL2 (0.0,1.0)mm  Tmax=70.730C  J=1.0292
  `heat_only` step  15 accept nudge GPU_2060S (0.0,-0.5)mm  Tmax=70.740C  J=1.0296
  `heat_only` step  16 accept nudge CPU_starter (-1.5,0.0)mm  Tmax=70.750C  J=1.0300
- `heat_only` rescoring best with silicon-fill stress …
- `heat_only` done  Tmax 70.550→70.550C  vm 149.00→149.00 MPa  wpv 4.221→4.221 um  accepted 8/16
- `coupled` init  Tmax=70.550C  vm_chipmax=149.00 MPa (GPU_TITANX)  w_pv=4.221 um  J=3.0397  stress=105.7s
  `coupled` step   1 BEST rotate CPU_armL2  Tmax=70.550C  vm=149.00  wpv=4.221  J=3.0397
  `coupled` step   3 BEST nudge GPU_3070 (-1.0,0.0)mm  Tmax=70.580C  vm=147.64  wpv=4.231  J=3.0285
  `coupled` step   6 accept swap CIM_wang<->CIM_colonnade  Tmax=70.580C  J=3.0286
  `coupled` step   9 BEST nudge GPU_3070 (-2.0,0.0)mm  Tmax=70.680C  vm=120.76  wpv=3.781  J=2.7074
  `coupled` step  10 accept nudge GPU_TITANX (0.0,0.5)mm  Tmax=70.730C  J=2.7591
  `coupled` step  12 accept nudge CPU_armL2 (0.0,1.0)mm  Tmax=70.730C  J=2.7591
  `coupled` step  15 accept nudge GPU_2060S (0.0,-0.5)mm  Tmax=70.740C  J=2.8109
  `coupled` step  16 accept nudge CPU_starter (-1.5,0.0)mm  Tmax=70.750C  J=2.8121
- `coupled` done  Tmax 70.550→70.680C  vm 149.00→120.76 MPa  wpv 4.221→3.781 um  accepted 8/16

## Table
| arm | Tmax init→best | vm_chip_max init→best | w_pv init→best | accepted |
|---|---|---|---|---|
| heat_only | 70.550→70.550 | 149.00→149.00 (GPU_TITANX) | 4.221→4.221 | 8/16 |
| coupled | 70.550→70.680 | 149.00→120.76 (GPU_TITANX) | 4.221→3.781 | 8/16 |

coupled − heat_only at best: dTmax=0.130 C, d_vm=-28.24 MPa, d_wpv=-0.440 um

wall 1231.1s  finished 2026-10-25T13:18:39.839884+00:00
