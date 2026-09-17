#!/usr/bin/env python3
"""Write the run configs and create the output directories.

Main set:   3 EOS x n in {64, 128, 256, 512}, N_NN = 32  -> config_<eos>_N<n>.info
NNN set:    3 EOS x n in {64, 128, 256, 512}, N_NN = 64  -> config_<eos>_N<n>_nnn64.info
Corrected:  Murnaghan x n in {64, 128, 256}, N_NN = 32, masses corrected for the
            kernel-sum density offset (generateIC.py --correct-nnn 32)
                                                        -> config_murnaghan_N<n>_corr.info

N_NN is a runtime parameter (sml.NNNTarget), so both sets use the same builds.
Box: periodic [0, 1) x [0, 1/4); see generateIC.py.
"""

import math
import os

EOS_ID = {"ideal": 0, "murnaghan": 1, "tillotson": 2}
BOX_Y = 0.25

HEADER = """\
; 2D linear travelling soundwave (Hopkins 2015, sec. 4.2.1), {eos} EOS,
; {n} x {ny} particles (n = {n} per wavelength), N_NN = {nnn}.{corrnote}
; Build: build_{eos}/ -- EOS {eosid}, ELASTIC 1 + GIZMO_ELASTIC_FLUX 0 (xSplitElasticHLLC),
;        mu = 0 (no solid dynamics), TENSILE_CORRECTION 1, SURFACE_VOLCORR 0,
;        EXPLICIT_VOL_INTEGRATION 1, PERIODIC_BOUNDARIES 1, Wendland C2. See makeBuilds.sh.
; Run from this directory:  build_{eos}/bin/mlh -c {cfg}   (or ./runAll.sh)

initFile {ic}
outDir   {out}

; One wave period, lambda/c = 1. ADAPTIVE_TIMESTEP: snapshots every
; timeStep*h5DumpInterval = 0.01.
timeStep 0.01
timeEnd  1.0
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

"""

EOS_BLOCK = {
    "ideal": """\
; Hopkins: gamma = 5/3, rho = 1, P = 3/5 -> c = 1 (u = 0.9 in the IC).
; c^2 = gamma (gamma - 1) u does not depend on rho, so the kernel-sum density
; offset only shifts P_bar.
hydro_gamma 1.6666666666666667
""",
    "murnaghan": """\
; read only for EOS 0, kept for config-parser compatibility
hydro_gamma 1.6666666666666667

; Repo-standard Murnaghan (testcases/rip_test): P = K0/n ((rho/rho0)^n - 1),
; c^2 = K0/rho0 (rho/rho0)^(n-1) = 1 at rho = rho0 = 1. The kernel-sum density
; offset (+4.2e-3 at N_NN = 32) gives P_bar = +4.2e-3 and c = 1.0021.
; mu = 0: no deviatoric stress.
materials {
    mat0 {
        matId 0
        K0    1.0
        n     2.0
        rho0  1.0
        mu    0.0
    }
}
""",
    "tillotson": """\
; read only for EOS 0, kept for config-parser compatibility
hydro_gamma 1.6666666666666667

; Hopkins-canonical Tillotson (testcases/periodic_wall_bounce): a = b = B = 0,
; alpha = beta = 0 -> P = A (rho/rho0 - 1) in every regime, c^2 = A/rho0 = 1.
; The kernel-sum density offset (+4.2e-3 at N_NN = 32) gives P_bar = +4.2e-3.
; IC u = 0.01 > 0 (u <= 0 logs an ERROR every step), below u_iv so the cold
; branch is used. mu = 0: no deviatoric stress.
materials {
    mat0 {
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
        mu    0.0
    }
}
""",
}

RUNS = [(eos, n, 32, False) for eos in EOS_ID for n in (64, 128, 256, 512)] \
     + [(eos, n, 64, False) for eos in EOS_ID for n in (64, 128, 256, 512)] \
     + [("murnaghan", n, 32, True) for n in (64, 128, 256)]

CORRNOTE = """
; MASSES CORRECTED for the kernel-sum density offset at N_NN = {nnn} (IC {ic}):
; <rho> = rho_bar = 1, so P_bar ~ 0 (the wave has tensile half-cycles) and c = 1.
; The offset remarks in the EOS block below refer to nominal masses."""

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    for eos, n, nnn, corr in RUNS:
        stem = "{}_N{}".format(eos, n) + ("" if nnn == 32 else "_nnn{}".format(nnn)) \
             + ("_corr" if corr else "")
        cfg, out = "config_{}.info".format(stem), "output_{}".format(stem)
        ic = "soundwave_{}.h5".format(stem) if corr else "soundwave_{}_N{}.h5".format(eos, n)
        if not os.path.exists(ic):
            raise SystemExit("missing IC " + ic)
        with open(cfg, "w") as f:
            f.write(HEADER.format(eos=eos, eosid=EOS_ID[eos], n=n, ny=int(round(BOX_Y * n)),
                                  nnn=nnn, cfg=cfg, out=out, boxy=BOX_Y, ic=ic,
                                  corrnote=CORRNOTE.format(nnn=nnn, ic=ic) if corr else "",
                                  ks=math.sqrt(nnn / math.pi) / n))
            f.write(EOS_BLOCK[eos])
        os.makedirs(out, exist_ok=True)
        print(cfg, "->", out)
