# Basalt impact — failure-model sweep

miluphcuda's canonical brittle-fragmentation benchmark is the basalt impact
(Nakamura & Fujiwara 1991; Benz & Asphaug 1994). Here a 2D basalt **target disk**
(R=5) is struck by a small fast **projectile disk** (R=0.8) from the left. Single
material (matId=0; the demonstrator has no multi-material ICs yet), normalized
units (rho0 = A = c = 1), so the impact runs in the gentle Mach<~0.3 regime the
MFM scheme is stable in — a hypervelocity (Mach>1) impact collapses the CFL
timestep.

The **same IC** (`basalt_impact-2D.h5`, ~8k particles, with Weibull flaws) is run
through every failure model, simple → complex, one result folder each:

| folder | flags | result (final svm = sqrt(2 J2)) |
|---|---|---|
| `01_elastic` | none | baseline, svm up to 0.077 |
| `02_vonmises` | `VON_MISES_PLASTICITY` | svm capped at Y0 = 0.010 |
| `03_mohrcoulomb` | `MOHR_COULOMB_PLASTICITY` | svm capped at Y_M = 0.030 |
| `04_druckerprager` | `DRUCKER_PRAGER_PLASTICITY` | svm capped at Y_M = 0.030 |
| `05_fragmentation` | `FRAGMENTATION` | damage on (neg.) pressure; fracture on impact side, interior intact |
| `06_fragmentation_damageS` | `FRAGMENTATION DAMAGE_ACTS_ON_S` | also damages deviatoric stress (sensitive: stress spikes at fracture faces) |
| `07_collins_damage` | `COLLINS_PLASTICITY FRAGMENTATION DAMAGE_ACTS_ON_S` | damage-blended yield; least fracture (plasticity absorbs energy) |

## Reproduce
```
python3 generateIC_basalt.py          # writes basalt_impact-2D.h5 (with flaws)
./run_sweep.sh 10                      # builds 7 binaries, runs all in parallel to t=10
python3 plot_sweep.py                  # sweep_stress.png, sweep_damage.png
```
Each combo is a separate compile-time build (flags are in `parameter.h`); the
script sed-toggles them, builds, and copies the binary into the result folder.

## Calibration notes
- `Y0=0.01` puts the yield surface just above the bulk stress (svm median ~0.002,
  p90 ~0.013) so the impact region yields while the bulk stays elastic.
- `DAMAGE_START_TIME=0.6` (set for the fragmentation builds): the kernel-sum
  density is truncated (rho ~ 0.9 -> spurious tension everywhere) for the first
  step until the smoothing length converges. Without the gate that startup
  transient latches **every** Weibull flaw at once and damage saturates the whole
  body. With it, damage is correctly confined to the impact-driven tension.
- Weibull `k=5e12`, `m=9.5`. With `MAX_NUM_FLAWS=1` only the weakest flaw per
  particle is stored — coarse; raise both for a graded crack-growth response.
