# Heat-only vs heat+stress SA — pkg_4cpu_2gpu_4cim

- started (UTC): 2026-10-23T08:26:04.984641+00:00
- title: 4 CPU + 2 GPU + 4 CIM（10 颗，新方案）
- why: 中等规模：两颗不同 GPU、四颗 CPU、四颗 CIM，比 13/14 芯省内存，仍有异构旋钮
- n_chiplets: 10
- seed=22  N_STEPS=16  T0=2.0  alpha=0.95
- gaps: any 0.5 mm, GPU–GPU 2.0 mm
- coupled J = (Tmax-45)/25.0 + 1.0*vm_chip_max/100.0 + 0.5*w_pv/4.0
- vm = worst interior von Mises among normal chiplets (CPU/GPU, CIM area≥1 mm²)
- ATSim package overrides: {"pitch_die_mm": 2.0, "pitch_coarse_mm": 8.0}
- one scheme at a time (12 GB); CG only, no spsolve

## Run log

- start outline 36.62 x 45.98 mm  min_gap=0.58 mm  GPU-GPU=2.36 mm
- `heat_only` init  Tmax=62.530C  vm_chipmax=76.03 MPa (CPU_starter)  w_pv=2.221 um  J=0.7012  stress=64.9s
  `heat_only` step   1 accept rotate CPU_armL2  Tmax=62.530C  J=0.7012
  `heat_only` step   3 BEST nudge GPU_3070 (-1.0,0.0)mm  Tmax=62.480C  vm=76.03  wpv=2.221  J=0.6992
  `heat_only` step   6 accept swap CIM_basic<->CIM_wang  Tmax=62.480C  J=0.6992
  `heat_only` step   9 BEST nudge GPU_3070 (-2.0,0.0)mm  Tmax=62.440C  vm=76.03  wpv=2.221  J=0.6976
  `heat_only` step  10 accept nudge CIM_jia (0.0,0.5)mm  Tmax=62.440C  J=0.6976
  `heat_only` step  12 accept nudge CPU_armL2 (0.0,1.0)mm  Tmax=62.440C  J=0.6976
  `heat_only` step  15 accept nudge GPU_TITANX (0.0,-0.5)mm  Tmax=62.440C  J=0.6976
  `heat_only` step  16 accept nudge CPU_starter (-1.5,0.0)mm  Tmax=62.440C  J=0.6976
- `heat_only` rescoring best with silicon-fill stress …
- `heat_only` done  Tmax 62.530→62.440C  vm 76.03→70.60 MPa  wpv 2.221→2.372 um  accepted 8/16
- `coupled` init  Tmax=62.530C  vm_chipmax=76.03 MPa (CPU_starter)  w_pv=2.221 um  J=1.7391  stress=64.2s
  `coupled` step   1 BEST rotate CPU_armL2  Tmax=62.530C  vm=76.03  wpv=2.221  J=1.7391
  `coupled` step   3 accept nudge GPU_3070 (-1.0,0.0)mm  Tmax=62.480C  J=1.7520
  `coupled` step   6 accept swap CIM_basic<->CIM_wang  Tmax=62.480C  J=1.7530
  `coupled` step   9 BEST nudge GPU_3070 (-2.0,0.0)mm  Tmax=62.440C  vm=70.60  wpv=2.372  J=1.7001
  `coupled` step  10 accept nudge CIM_jia (0.0,0.5)mm  Tmax=62.440C  J=1.7001
  `coupled` step  12 accept nudge CPU_armL2 (0.0,1.0)mm  Tmax=62.440C  J=1.7002
  `coupled` step  15 BEST nudge GPU_TITANX (0.0,-0.5)mm  Tmax=62.440C  vm=70.76  wpv=2.329  J=1.6964
  `coupled` step  16 accept nudge CPU_starter (-1.5,0.0)mm  Tmax=62.440C  J=1.7810
- `coupled` done  Tmax 62.530→62.440C  vm 76.03→70.76 MPa  wpv 2.221→2.329 um  accepted 8/16

## Table
| arm | Tmax init→best | vm_chip_max init→best | w_pv init→best | accepted |
|---|---|---|---|---|
| heat_only | 62.530→62.440 | 76.03→70.60 (CPU_starter) | 2.221→2.372 | 8/16 |
| coupled | 62.530→62.440 | 76.03→70.76 (CPU_starter) | 2.221→2.329 | 8/16 |

coupled − heat_only at best: dTmax=0.000 C, d_vm=0.16 MPa, d_wpv=-0.042 um

wall 747.8s  finished 2026-10-23T08:38:32.717182+00:00
