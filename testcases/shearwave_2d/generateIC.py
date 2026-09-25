#!/usr/bin/env python3
"""2D linear travelling SHEAR wave -- the elastically dominated companion to
testcases/soundwave_2d_elastic.

Same geometry as the soundwave cases (periodic [0,1) x [0,1/4), square lattice of n x n/4
equal-mass particles, wave in +x, uniform in y, A = 1e-6, Wendland C2, N_NN = 32) and the
same material (Tillotson in Hopkins-canonical form, K = 1, mu = 6/13 from Poisson's ratio
0.3), but the perturbation is transverse:

    rho  = rho_bar                        (a shear wave perturbs rho only at O(A^2))
    v_y  = c_s A sin(k x),  v_x = 0       c_s = sqrt(mu/rho_0) = 0.679366220
    S_xy = -rho c_s v_y,    S_xx = S_yy = 0

Why this test exists. The elastic P-wave of 2026-09-18 converged at order 2, but its
restoring traction is only 38 % deviatoric -- and that is structural, not a parameter choice:
for a longitudinal wave f = (4/3)mu/(K + (4/3)mu) never exceeds 2/3, even at nu = 0. The
shear wave is the mode whose restoring force is **100 % deviatoric**, and which to leading
order carries no density perturbation, so the EOS, the kernel-sum density offset and the
tensile correction all drop out of it. What is left is pure solid dynamics: the stress
integrator, the S <-> v coupling, the deviatoric flux and the elastic wave fan.

The period is lambda/c_s = 1.471960144, i.e. 1.87x the P-wave period (c_l = 1.270978), so at
the same dt a three-passing run takes 1.87x as many steps.

A travelling shear wave needs S_xy in the IC (S = 0 would give a standing wave); this is
possible since fd88e21, which made InitialDistribution read /Sxx, /Sxy, /Syy.

The mode, the geometry and the kernel-offset mass correction all come from
testcases/wave_validation/waveic.py; this script only fixes the defaults of this testcase.

Usage:
    ./generateIC.py -n 64
    ./generateIC.py --all                 # n = 64, 128, 256
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "wave_validation"))
import waveic  # noqa: E402

RESOLUTIONS = (64, 128, 256)
CELLS = (4, 1)          # box [0,1) x [0,1/4) in lattice cells
MODES = (1, 0)          # one wavelength along x, uniform in y


def make_ic(n, amplitude, mu, nnn, out=None):
    if n % 4:
        raise SystemExit("n must be a multiple of 4 (box height 1/4), got {}".format(n))
    fname = out or "shearwave_N{}{}.h5".format(n, "" if nnn is None else "_corr")
    return waveic.build(CELLS, n // 4, MODES, mode="s", eos="tillotson", mu=mu,
                        amplitude=amplitude, nnn=nnn, out=fname)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", "--resolution", type=int, help="particles per wavelength")
    ap.add_argument("-A", "--amplitude", type=float, default=1e-6,
                    help="v_y amplitude in units of c_s (default: %(default)s)")
    ap.add_argument("--mu", type=float, default=waveic.MU_NU03,
                    help="shear modulus (default: %(default).9g, i.e. Poisson 0.3 at K = 1)")
    ap.add_argument("--correct-nnn", type=int, default=32, metavar="N_NN",
                    help="divide masses by 1 + kernel-sum offset for this N_NN (default: %(default)s)")
    ap.add_argument("--no-correct", action="store_true", help="nominal masses")
    ap.add_argument("-o", "--out", help="output file name (single resolution only)")
    ap.add_argument("--all", action="store_true",
                    help="generate all resolutions {}".format(RESOLUTIONS))
    args = ap.parse_args()

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    nnn = None if args.no_correct else args.correct_nnn
    if args.all:
        for n in RESOLUTIONS:
            make_ic(n, args.amplitude, args.mu, nnn)
    else:
        if args.resolution is None:
            ap.error("need -n, or --all")
        make_ic(args.resolution, args.amplitude, args.mu, nnn, args.out)
