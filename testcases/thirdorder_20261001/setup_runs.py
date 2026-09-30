#!/usr/bin/env python3
"""Create run dirs run_<tag>_N<n>/ (config.info + BUILD) for the soundwave ladder.
Usage: ./setup_runs.py tag build n1 n2 ... [--cfl C] [--tend T] [--nnn N]
  dt = C / n (c = 1, dx = 1/n), default C = 1/3 (the 09-25 fixed-CFL ladder, M = 3n steps per period).
  Dumps every M/4 steps (h5DumpInterval = M/4 + 1: the demonstrator's non-adaptive counter advances
  twice on a dump step), so t = 1 is a dump."""
import argparse, os, math
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("tag"); ap.add_argument("build"); ap.add_argument("ns", type=int, nargs="+")
ap.add_argument("--cfl", type=float, default=1/3); ap.add_argument("--tend", type=float, default=1.0)
ap.add_argument("--nnn", type=int, default=32); ap.add_argument("--ic", default="soundwave_ideal_N{n}.h5")
a = ap.parse_args()
for n in a.ns:
    M = round(1/(a.cfl/n))            # steps per unit time
    assert abs(M*a.cfl/n - 1) < 1e-12, "C/n must divide 1"
    assert M % 4 == 0
    dt = 1.0/M
    d = os.path.join(HERE, "run_%s_N%d" % (a.tag, n)); os.makedirs(d, exist_ok=True)
    ks = math.sqrt(a.nnn/math.pi)/n
    cfg = f"""; thirdorder study 2026-10-01: 2D linear soundwave (Hopkins 2015 sec. 4.2.1), ideal gas gamma 5/3,
; P = 3/5, c = 1, A = 1e-6, box 1 x 0.25, n = {n} per wavelength. Build {a.build} (see makeBuilds.sh).
; fixed dt = {a.cfl:.6g}/n (M = {M} steps per unit time), dumps every M/4 steps.
initFile ics/{a.ic.format(n=n)}
outDir run_{a.tag}_N{n}/output

timeStep {dt!r}
timeEnd {a.tend + 0.1}
h5DumpInterval {M//4 + 1}

periodicBoxLimits {{
    lowerX 0.0
    lowerY 0.0
    lowerZ 0.0
    upperX 1.0
    upperY 0.25
    upperZ 0.0
}}

kernelSize {ks!r}

sml {{
    NNNTarget {a.nnn}
    tol       1e-3
    maxIter   200
    maxFactor 1.2
    hMinFactor 0.1
    hMaxFactor 3.0
    warnFraction  0.05
    panicFraction 0.25
}}

hydro_gamma 1.6666666666666667
"""
    open(os.path.join(d, "config.info"), "w").write(cfg)
    open(os.path.join(d, "BUILD"), "w").write("build_%s\n" % a.build)
    os.makedirs(os.path.join(d, "output"), exist_ok=True)
    print(d)
