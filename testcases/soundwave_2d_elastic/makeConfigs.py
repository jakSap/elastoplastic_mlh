#!/usr/bin/env python3
"""Write the run configs and create the output directories.

n in {64, 128, 256}, N_NN = 32, Tillotson with mu > 0, masses corrected for the
kernel-sum density offset -> config_elastic_N<n>_corr.info.

mu, c_l and the period are read from the IC file's attributes, so changing
generateIC.py's --mu automatically retimes the runs: PASSINGS wave periods
lambda/c_l, with SNAPS_PER_PASSING snapshots each, so a snapshot sits exactly
on every period boundary and must reproduce the IC.

Dump-index caveat (MeshlessScheme.cpp:193): with ADAPTIVE_TIMESTEP the next
dump is scheduled against dumpTimes[dumpStep+1] although dumpStep has already
been incremented by the dump, so dumpTimes[1] is skipped and file k > 0 is
written at t = (k+1)*timeStep, not k*timeStep. The period boundary k*T is
therefore file k*SNAPS_PER_PASSING - 1, and timeEnd carries TAIL_DUMPS extra
dump intervals so the last boundary is actually reached before the loop exits.
Always select snapshots by their stored /time, never by file index.
"""

import math
import os
import h5py

RESOLUTIONS = (64, 128, 256)
NNN = 32
BOX_Y = 0.25
PASSINGS = 3
SNAPS_PER_PASSING = 50
# extra dump intervals past the last period boundary, see the dump-index caveat
TAIL_DUMPS = 3

TEMPLATE = """\
; 2D linear travelling ELASTIC wave (Hopkins 2015, sec. 4.2.1, with mu > 0),
; Tillotson EOS, {n} x {ny} particles (n = {n} per wavelength), N_NN = {nnn}.
;
; MASSES CORRECTED for the kernel-sum density offset at N_NN = {nnn}
; ({offset:+.8e}, IC {ic}): <rho> = rho_bar = 1, so P_bar = 0 and the wave has
; genuine tensile half-cycles.
;
; The IC also carries /Sxx, /Sxy, /Syy: with mu > 0 the travelling-wave
; eigenmode has S_xx = -4/3 mu ln(rho/rho0), and starting from S = 0 would split
; the wave. See generateIC.py.
;
; Build: build_elastic/ -- EOS 2, ELASTIC 1 + GIZMO_ELASTIC_FLUX 0 (xSplitElasticHLLC),
;        mu = {mu:.9g} (Poisson 0.3 at K = 1), TENSILE_CORRECTION 1 + _LOCAL 1,
;        SURFACE_VOLCORR 0, EXPLICIT_VOL_INTEGRATION 1, PERIODIC_BOUNDARIES 1,
;        Wendland C2. See makeBuilds.sh.
; Run from this directory:  build_elastic/bin/mlh -c {cfg}   (or ./runAll.sh)

initFile {ic}
outDir   {out}

; Longitudinal wave speed c_l = sqrt((K + 4/3 mu)/rho0) = {cl:.9f}, so one wave
; passing takes lambda/c_l = {period:.9f} -- NOT 1 as in the mu = 0 case.
; {passings} passings, {snaps} snapshots each. Because of the dump-index
; off-by-one in MeshlessScheme.cpp:193 (file k > 0 is written at t = (k+1)*dt),
; the period boundaries are files {b1}, {b2} and {b3}; timeEnd runs {tail} dump
; intervals past {passings} T so the last one is actually written.
timeStep {dt:.12g}
timeEnd  {tend:.12g}
h5DumpInterval 1

periodicBoxLimits {{
    lowerX 0.0
    lowerY 0.0
    lowerZ 0.0
    upperX 1.0
    upperY {boxy}
    upperZ 0.0
}}

; ~ kernel support for N_NN = {nnn} on the lattice: sqrt(N_NN/pi) * dx, dx = 1/{n}
kernelSize {ks:.6g}

sml {{
    NNNTarget {nnn}
    tol       1e-3
    maxIter   200
    maxFactor 1.2
    hMinFactor 0.1
    hMaxFactor 3.0
    warnFraction  0.05
    panicFraction 0.25
}}

; read only for EOS 0, kept for config-parser compatibility
hydro_gamma 1.6666666666666667

; Hopkins-canonical Tillotson (as testcases/soundwave_2d): a = b = B = 0,
; alpha = beta = 0 -> P = A (rho/rho0 - 1) in every regime, K = A/rho0 = 1.
; IC u = 0.01 > 0 (u <= 0 logs an ERROR every step), below u_iv so the cold
; branch is used.
;
; mu > 0 is the difference to the mu = 0 baseline: the deviatoric stress is
; live, the normal traction is P - S_xx, and the wave runs at c_l, not at
; sqrt(K/rho0) = 1.
materials {{
    mat0 {{
        matId 0
        rho0  1.0
        A     1.0
        B     0.0
        u0    1.0
        a     0.0
        b     0.0
        alpha 0.0
        beta  0.0
        u_iv  0.05
        u_cv  0.18
        mu    {mu:.17g}
    }}
}}
"""

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    for n in RESOLUTIONS:
        stem = "elastic_N{}_corr".format(n)
        ic = "soundwave_{}.h5".format(stem)
        if not os.path.exists(ic):
            raise SystemExit("missing IC {} -- run ./generateIC.py --all".format(ic))
        with h5py.File(ic, "r") as f:
            mu, cl = float(f.attrs["mu"]), float(f.attrs["c_l"])
            period, offset = float(f.attrs["period"]), float(f.attrs["kernel_offset"])
        cfg, out = "config_{}.info".format(stem), "output_{}".format(stem)
        with open(cfg, "w") as f:
            f.write(TEMPLATE.format(
                n=n, ny=int(round(BOX_Y * n)), nnn=NNN, cfg=cfg, out=out, ic=ic,
                boxy=BOX_Y, ks=math.sqrt(NNN / math.pi) / n, mu=mu, cl=cl,
                period=period, offset=offset, passings=PASSINGS,
                snaps=SNAPS_PER_PASSING,
                b1=SNAPS_PER_PASSING - 1, b2=2 * SNAPS_PER_PASSING - 1,
                b3=3 * SNAPS_PER_PASSING - 1, tail=TAIL_DUMPS,
                dt=period / SNAPS_PER_PASSING,
                tend=(PASSINGS * SNAPS_PER_PASSING + TAIL_DUMPS) * period / SNAPS_PER_PASSING))
        os.makedirs(out, exist_ok=True)
        print("{} -> {}  (t_end = {:.6f} = {} x {:.6f} + {} dumps, dt_dump = {:.6f})".format(
            cfg, out, (PASSINGS * SNAPS_PER_PASSING + TAIL_DUMPS) * period / SNAPS_PER_PASSING,
            PASSINGS, period, TAIL_DUMPS, period / SNAPS_PER_PASSING))
