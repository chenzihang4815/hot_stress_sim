# Heat-only vs heat+stress SA (13 chip, routing gaps)

- started (UTC): 2026-10-22T14:27:45.199031+00:00
- seed=22  N_STEPS=16  T0=2.0  alpha=0.95
- gaps: any 0.5 mm, GPU–GPU 2.0 mm
- coupled J = (Tmax-45)/25.0 + 1.0*vm_mean/70.0 + 0.5*w_pv/4.0
- stress proxy: ATSim fill_die silicon, chip-center von Mises, CG no spsolve
- STRESS_EVERY=1 (coupled trials)
- does **not** overwrite `runs/sa_gapped_13chip`

## Run log

- start outline 57.12 x 40.55 mm  min_gap=0.54 mm  GPU-GPU=2.62 mm
- `heat_only` init  Tmax=70.990C  vm_mean=70.40 MPa  w_pv=4.457 um  J=1.0396  stress=171.9s
  `heat_only` step   1 accept rotate CPU_armL2  Tmax=70.990C  J=1.0396
  `heat_only` step   4 accept rotate CIM_wang  Tmax=70.990C  J=1.0396
  `heat_only` step   6 BEST jump GPU_3070  Tmax=70.040C  vm=70.40  wpv=4.457  J=1.0016
  `heat_only` step   8 accept nudge GPU_2060S (-1.0,0.0)mm  Tmax=70.620C  J=1.0248
  `heat_only` step   9 accept nudge CIM_wang (-1.5,0.0)mm  Tmax=70.630C  J=1.0252
  `heat_only` step  10 accept nudge GPU_2060S (0.0,0.5)mm  Tmax=70.330C  J=1.0132
  `heat_only` step  11 accept nudge CIM_jia (0.0,-1.0)mm  Tmax=70.330C  J=1.0132
  `heat_only` step  12 accept rotate CPU_armL2  Tmax=70.330C  J=1.0132
  `heat_only` step  15 accept nudge GPU_3070 (0.0,-0.5)mm  Tmax=70.780C  J=1.0312
  `heat_only` step  16 accept nudge CPU_atomic (-1.5,0.0)mm  Tmax=70.730C  J=1.0292
- `heat_only` rescoring best with silicon-fill stress …
- `heat_only` done  Tmax 70.990→70.040C  vm 70.40→56.89 MPa  wpv 4.457→4.479 um  accepted 10/16
- `coupled` init  Tmax=70.990C  vm_mean=70.40 MPa  w_pv=4.457 um  J=2.6024  stress=142.6s
  `coupled` step   1 accept rotate CPU_armL2  Tmax=70.990C  J=2.6168
  `coupled` step   4 accept rotate CIM_wang  Tmax=70.990C  J=2.6175
  `coupled` step   6 BEST jump GPU_3070  Tmax=70.040C  vm=56.89  wpv=4.479  J=2.3743
  `coupled` step   8 BEST nudge GPU_2060S (-1.0,0.0)mm  Tmax=70.620C  vm=53.93  wpv=4.480  J=2.3553
  `coupled` step   9 accept nudge CIM_wang (-1.5,0.0)mm  Tmax=70.630C  J=2.3746
  `coupled` step  10 accept nudge GPU_2060S (0.0,0.5)mm  Tmax=70.330C  J=2.3880
  `coupled` step  11 accept nudge CIM_jia (0.0,-1.0)mm  Tmax=70.330C  J=2.3894
  `coupled` step  12 accept rotate CPU_armL2  Tmax=70.330C  J=2.3894
  `coupled` step  15 accept nudge GPU_3070 (0.0,-0.5)mm  Tmax=70.780C  J=2.4088
  `coupled` step  16 BEST nudge CPU_atomic (-1.5,0.0)mm  Tmax=70.730C  vm=52.56  wpv=4.472  J=2.3390
- `coupled` done  Tmax 70.990→70.730C  vm 70.40→52.56 MPa  wpv 4.457→4.472 um  accepted 10/16

## Table
| arm | Tmax init→best | vm_mean init→best | w_pv init→best | accepted |
|---|---|---|---|---|
| heat_only | 70.990→70.040 | 70.40→56.89 | 4.457→4.479 | 10/16 |
| coupled | 70.990→70.730 | 70.40→52.56 | 4.457→4.472 | 10/16 |

coupled − heat_only at best: dTmax=0.690 C, d_vm=-4.34 MPa, d_wpv=-0.007 um

wall 1981.8s  finished 2026-10-22T15:00:46.994211+00:00
