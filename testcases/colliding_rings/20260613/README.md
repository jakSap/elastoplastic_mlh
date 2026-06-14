# 20260613 — GIZMO tensile-correction sweep

Sweep over the GIZMO tensile-instability measures ported to the demonstrator,
on the colliding rubber rings test. IC (`rings_hopkins.h5`) and all config
values are identical to `../20260611/run_1` (Hopkins relaxed-glass IC, N=5932,
NNNTarget=32, timeEnd=300). Only compile-time flags in
`demonstrator/include/parameter.h` differ between variants; each variant
directory contains the exact `parameter.h` snapshot it was built with.

## Flags (GIZMO source mapping)

| Flag | GIZMO source | Effect |
|---|---|---|
| `KERNEL_FUNCTION` | `kernel.h` (same numbering: 3/6/7/9) | Kernel: cubic spline / Wendland C2 / C4 / C6 |
| `TENSILE_CORRECTION` | `elastic_physics.c` `get_negative_pressure_tensilecorrfac` | Master switch: computes the Monaghan 2000 factor f = 0.2 (W(q)/W(Δp))⁴ with the active kernel |
| `TENSILE_CORRECTION_1` | `hydro/hydro_core_meshless.h` dummy_pressure | Negative-pressure shift before the Riemann problem, re-subtracted scaled by (1−f) |
| `TENSILE_CORRECTION_2` | `solids/elastic_stress_tensor_force.h` `#if 0` branch | Damp the whole deviatoric stress of a negative-pressure side by (1−f) |
| `TENSILE_CORRECTION_3` | `solids/elastic_stress_tensor_force.h` default branch | Damp only tensile (positive) principal stress components by (1−f), always |

`TENSILE_CORRECTION_2` and `_3` are mutually exclusive (GIZMO's `#if 0/#else`).
GIZMO's commented-out approach-limiter block (elastic_stress_tensor_force.h:52-73)
is not ported: it is disabled in GIZMO itself and requires the split
hydro/stress flux structure the demonstrator does not have.

## Variants

| Directory | Kernel | TC1 | TC2 | TC3 | Note |
|---|---|---|---|---|---|
| cubic_tc_none | cubic | – | – | – | no correction baseline |
| cubic_tc1 | cubic | x | – | – | = 20260611/run_1 setup |
| cubic_tc2 | cubic | – | x | – | |
| cubic_tc3 | cubic | – | – | x | |
| cubic_tc1_tc2 | cubic | x | x | – | |
| cubic_tc1_tc3 | cubic | x | – | x | |
| wc2_tc_none | Wendland C2 | – | – | – | |
| wc2_tc1 | Wendland C2 | x | – | – | |
| wc2_tc2 | Wendland C2 | – | x | – | |
| wc2_tc3 | Wendland C2 | – | – | x | |
| wc2_tc1_tc2 | Wendland C2 | x | x | – | = GIZMO simple_stressflux variant |
| wc2_tc1_tc3 | Wendland C2 | x | – | x | = GIZMO reference (KERNEL_FUNCTION=6, default branches) |
| wc4_tc1_tc3 | Wendland C4 | x | – | x | GIZMO defaults with C4 |

Reproduce with `./run_all.sh`. Per-variant: `config.info`, `parameter.h`,
`build.log`, `log.txt`, `exit_code`, `output/*.h5`; sweep result in `summary.txt`.
