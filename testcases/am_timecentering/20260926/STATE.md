# am_timecentering 2026-09-26 — state for a resumed turn

Task (user, 2026-09-26 ~12:25): AM conservation from the time-centering angle; analyse RECON_FACE_DV
from it; build/test fixes; goal = 2nd-order L_z error convergence on the rotating disc (smooth field).
Second angle (user, mid-turn): elastic wave structure smoothed by the HLL of the GIZMO stress flux,
which does not discriminate between wave types.

## DONE
- v_int sign fix committed: newMFM `6b89f87` (boost check runs/inv_ref vs inv_boost: max|du| 2e-14).
- Findings (all measured on disc d0.05 omega 0.2, one rotation unless noted):
  1. DRIFT_KDK 1 (standard) -> dL/L = +omega^2 dt/2 per unit time (kick at x^n, drift with v^{n+1/2});
     7 testplan disc runs ratio meas/pred 0.89-0.97. Endpoint drift or TIME_DKD kill it (8.2e-3 -> ~5e-4).
  2. Same cause: rhoE/rhoK drifts +omega^2 dt per unit time (kick B uses grad v^{n+1} on x^n), disc inflates.
  3. RECON_FACE_DV jump carries the MUSCL predictor, which misses dt/2 (x_f-x_i).d_t grad v -> spurious
     jump dt/2 Gdot x_ij -> transverse torque first order in dt (TQDIAG wtT 3.2e-4 -> 1.7e-4 at CFL/2).
     New flag RECON_FACE_DV_SPATIAL 1 (uncommitted) -> wtT 12x smaller and dt-independent.
  4. RECON_FACE_DV 2 with predicted jump explodes (C_dkd_dv2_tq: rmax 1382); with spatial jump stable.
  5. New flag TIME_KDK 1 (GIZMO kick-drift-kick, needs DRIFT_KDK 0, RECON_TIME_PRED 0; uncommitted):
     rhoE/rhoK 3.7e-6, dE 1e-14, max|dL/L| 2.3e-4 (dev channel, spatial).
  6. Default build with new flags off is bit-identical (runs/inv_defcheck vs inv_ref).
  7. Waves N64: RECON_FACE_DV 2 cuts damping 13x (P), 15x (shear), 17x (sound).
- Remaining disc error = deviatoric traction channel (SPH area), concentrated in the 3d surface band,
  net ~1% of sum|tau_i|; grows d0.05 -> d0.035 in TIME_DKD+sp. ELASTIC_FACE_MODE 2 ladders queued.

## RUNNING (queues under setsid; queue.log; summary.txt)
- jobs_ladder.txt  : std / dkd_sp / dkd_sp_fm2 at d 0.1..0.0177 (+ dkd_sp d0.0125)
- jobs_ladder2.txt : gkdk / gkdk_dv2
- jobs_ladder3.txt : gkdk_fm2 / gkdk_fm2_dv2
- jobs_waves.txt, jobs_waves2.txt : W1 (tillotson, ideal), W2, W3 at N64/128/256 for std, sp, sp_dv2,
  dkd_sp_dv2, gkdk, gkdk_dv2, e0_* (ideal gas).

## ANALYSIS
- `python3 ladder.py [variants]` (orders), `python3 tqparse.py run`, `python3 tqmap.py run`,
  `python3 wave_table.py`, `python3 analyze.py run`.
- Theory draft: scratchpad theory_draft.md (copy in report).

## LEFT
- Fit orders; decide scheme (TIME_KDK vs TIME_DKD+SPATIAL; face mode; RECON_FACE_DV 2).
- Movies (disc_movie.py) of the top rungs + the exploding dv2 run; send via SendUserFile.
- Report (German, verbose + summary + daily note), errors first. Commit flags (ask user about defaults).

## Update ~14:30
- New flags (uncommitted): ELASTIC_FACE_MODE 3 (split MFM half faces, dev torque = round-off with TC 0),
  RECON_FACE_DV_MINMOD (radial channel only; transverse minmod brakes spin -2.5e-3, see X_gkdk_dv2_mmboth_d0.025).
- Standard ladder (square lattice): p = 0.98, R2 1.000 (5.9e-3 / 4.2e-3 / 3.0e-3 at d .035/.025/.0177).
- TIME_DKD + SPATIAL + RECON_FACE_DV 2 (face 0, TC on): stable, 5.1e-4 .. 5.3e-5 (d .1 .. .025), p ~1.5 (pairwise 2.0/0.35/1.7/2.5).
- TIME_DKD + SPATIAL + face 3 + TC 0: 8.8e-5 .. 1.7e-5 (d .07 .. .025), pairwise 1.9/1.0/2.0.
- Square-lattice cut is NOT a smooth refinement family (Lz/Lrig 0.985/1.00006/0.994/...): new concentric-ring
  IC gen/disc_rings_gen.py (ics/discR_*), ladder queued: jobs_rings.txt (std, dkd_sp_dv2, dkd_sp_fm3_tc0, dkd_sp_fm3).
- ERROR RUNS: TIME_KDK + RECON_FACE_DV 2 disrupts at t~10 (d.025: 4.3M negE, u -133; d.0177: 161k) -> killed,
  logs trimmed. TIME_KDK + face 2/3 cools one outermost surface particle (u_min 0.07/0.50 at d.035).
  TIME_KDK deprioritised; TIME_DKD route is the candidate.
- Killed by me (exit 143): L_dkd_sp_d0.0125, L_gkdk_dv2_d0.0177/.025, L_gkdk_fm2_d0.0177, L_dkd_sp_fm2_d0.0177,
  L_gkdk_fm2_tc0_d0.025, L_gkdk_fm3_tc0_d0.0177, L_gkdk_dv2_mm_d0.0177, W1_sound/gkdk_fm2_dv2_N128,
  W3_shear/gkdk_fm2_dv2_tc0_N128. skip_lowprio dirs = never run.

## Update ~14:45
- Waves: RECON_FACE_DV 2 + SPATIAL + ELASTIC_FACE_MODE 2/3 on the standard drift = 2nd order in W1/W2/W3
  (sound N256 p 2.05). TIME_DKD version: velocity 2nd, but density-coupled (L1rho sound p 1.0, P-wave defect p 1.0,
  ideal gas p 1.2-1.35) FIRST order -> TIME_DKD costs wave density order.
- Disc rings: TIME_DKD+SP+DV2+FM3 (TC on) p 1.85 over .0714-.0357; u_min sinks at the ring-layout centre
  (ring k~2) -> radial minmod variant (_mmr) on rings queued (jobs_rings3).
- TIME_KDK: surface-particle cooling only with TC OFF; with TC on u_min 0.87-0.92. KDK+DV2 needs radial minmod
  (L_gkdk_dv2_mmr_d0.025 stable, 9.8e-5). Final KDK candidate gkdk_fm3_dv2_mmr (TC on) queued on rings + all waves
  (queue_mixed.sh jobs_kdkfinal.txt).
- figs/timecentering.png done (fig_timecentering.py). Failure movies sent (movies/fail_*.mp4).

## FINAL ~16:45
Reports written: MA-Obsidian/2026-09-26-am-zeitzentrierung{,-summary}.md + daily note; figs copied (fig1, leiter),
failure movies in vault and sent. Memory: project_am_timecentering, feedback_no_self_kill.
Still running (not needed for the conclusions): L_Rdkd_sp_fm3_tc0_d0.0179 and leftover queue items.
Nothing committed except 6b89f87 (sign fix). New flags live uncommitted in demonstrator/ (defaults 0, bit-identical).
