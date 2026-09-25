#!/usr/bin/env python3
"""Write the run configs and create the output directories for the 2D shear wave.

n in {64, 128, 256}, N_NN = 32, Tillotson with mu > 0, masses corrected for the kernel-sum
density offset -> config_shear_N<n>_corr.info.

c_s and the period are read from the IC file's attributes, so changing generateIC.py's --mu
automatically retimes the runs: PASSINGS wave periods lambda/c_s with SNAPS_PER_PASSING
snapshots each, so a snapshot sits on every period boundary and must reproduce the IC.

Note the period is 1.471960 (c_s = 0.679366), NOT 0.786796 as in the elastic P-wave case:
at the same dt (set by the *longitudinal* signal speed c_l = 1.270978) a three-passing shear
run therefore takes ~1.87x as many steps as the P-wave run at the same resolution.

Dump-index caveat (MeshlessScheme.cpp:193): with ADAPTIVE_TIMESTEP the next dump is scheduled
against dumpTimes[dumpStep+1] although dumpStep has already been incremented, so dumpTimes[1]
is skipped and file k > 0 is written at t = (k+1)*timeStep. The period boundary k*T is
therefore file k*SNAPS_PER_PASSING - 1, and timeEnd carries TAIL_DUMPS extra dump intervals so
the last boundary is reached before the loop exits. Always select snapshots by stored /time.

CFL is a compile-time flag (parameter.h:151), not a config key: the CFL control run uses
its own build directory (cflTest.sh) and only needs a distinct run stem, hence --suffix.

Usage: ./makeConfigs.py  [--suffix _cfl05]
"""

import argparse
import math
import os

import h5py

RESOLUTIONS = (64, 128, 256)
NNN = 32
BOX_Y = 0.25
PASSINGS = 3
SNAPS_PER_PASSING = 50
TAIL_DUMPS = 3

TEMPLATE = """\
; 2D linear travelling SHEAR wave, Tillotson EOS with mu > 0, {n} x {ny} particles
; (n = {n} per wavelength), N_NN = {nnn}.
;
; MASSES CORRECTED for the kernel-sum density offset at N_NN = {nnn}
; ({offset:+.8e}, IC {ic}): <rho> = rho_bar = 1, so P_bar = 0.
;
; The IC carries /Sxy = -rho c_s v_y: a travelling shear wave is only an eigenmode with that
; stress (S = 0 would give a standing wave). rho is uniform -- a shear wave perturbs the
; density only at O(A^2) -- so EOS, kernel-density offset and tensile correction drop out of
; the restoring force, which is 100 % deviatoric. See generateIC.py.
;
; Build: build_shear/ -- EOS 2, ELASTIC 1 + GIZMO_ELASTIC_FLUX 0 (xSplitElasticHLLC),
;        mu = {mu:.9g} (Poisson 0.3 at K = 1), TENSILE_CORRECTION 1 + _LOCAL 1,
;        SURFACE_VOLCORR 0, EXPLICIT_VOL_INTEGRATION 1, PERIODIC_BOUNDARIES 1,
;        Wendland C2. See makeBuilds.sh.
; Run from this directory:  build_shear/bin/mlh -c {cfg}   (or ./runAll.sh)

initFile {ic}
outDir   {out}

; Shear wave speed c_s = sqrt(mu/rho0) = {cs:.9f}, so one passing takes lambda/c_s =
; {period:.9f}. {passings} passings, {snaps} snapshots each. Because of the dump-index
; off-by-one in MeshlessScheme.cpp:193 (file k > 0 is written at t = (k+1)*dt), the period
; boundaries are files {b1}, {b2} and {b3}; timeEnd runs {tail} dump intervals past
; {passings} T so the last one is actually written.
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

; Hopkins-canonical Tillotson: a = b = B = 0, alpha = beta = 0 -> P = A (rho/rho0 - 1),
; K = A/rho0 = 1. IC u = 0.01 > 0 (u <= 0 logs an ERROR every step), below u_iv so the cold
; branch is used. With uniform rho and corrected masses P = 0 everywhere at t = 0.
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
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--resolutions", type=int, nargs="+", default=RESOLUTIONS)
    ap.add_argument("--suffix", default="", help="append to the run stem, e.g. _cfl05")
    args = ap.parse_args()

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    for n in args.resolutions:
        ic = "shearwave_N{}_corr.h5".format(n)
        if not os.path.exists(ic):
            raise SystemExit("missing IC {} -- run ./generateIC.py --all".format(ic))
        with h5py.File(ic, "r") as f:
            mu, cs = float(f.attrs["mu"]), float(f.attrs["c_wave"])
            period, offset = float(f.attrs["period"]), float(f.attrs["kernel_offset"])
        stem = "shear_N{}_corr{}".format(n, args.suffix)
        cfg, out = "config_{}.info".format(stem), "output_{}".format(stem)
        with open(cfg, "w") as f:
            f.write(TEMPLATE.format(
                n=n, ny=int(round(BOX_Y * n)), nnn=NNN, cfg=cfg, out=out, ic=ic,
                boxy=BOX_Y, ks=math.sqrt(NNN / math.pi) / n, mu=mu, cs=cs,
                period=period, offset=offset, passings=PASSINGS, snaps=SNAPS_PER_PASSING,
                b1=SNAPS_PER_PASSING - 1, b2=2 * SNAPS_PER_PASSING - 1,
                b3=3 * SNAPS_PER_PASSING - 1, tail=TAIL_DUMPS,
                dt=period / SNAPS_PER_PASSING,
                tend=(PASSINGS * SNAPS_PER_PASSING + TAIL_DUMPS) * period / SNAPS_PER_PASSING))
        os.makedirs(out, exist_ok=True)
        print("{} -> {}  (t_end = {:.6f} = {} x {:.6f} + {} dumps, dt_dump = {:.6f})".format(
            cfg, out, (PASSINGS * SNAPS_PER_PASSING + TAIL_DUMPS) * period / SNAPS_PER_PASSING,
            PASSINGS, period, TAIL_DUMPS, period / SNAPS_PER_PASSING))
