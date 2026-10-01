#!/usr/bin/env python3
"""Run dirs disc/run_<tag>_d<d>/ for the rotating disc (test plan A.1: R = 1, omega = 0.2, Tillotson rho0 = A = 1,
mu = 0.22, u = 1) on the thirdorder builds. Fixed dt (TIME_RK3/RK4 need it; the second-order builds get the same dt):
M steps per rotation with c dt/d = CFL, c = c_l + omega R (default CFL 1/3), 32 dumps per rotation.
Paths are relative to the parent study dir, where run_queue.sh runs.
Usage: disc/setup_disc.py tag build d [--rot 3.2] [--cfl 0.3333]"""
import argparse, os, math, h5py
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("tag"); ap.add_argument("build"); ap.add_argument("d", type=float)
ap.add_argument("--rot", type=float, default=3.2); ap.add_argument("--cfl", type=float, default=1/3)
ap.add_argument("--nnn", type=int, default=32)
a = ap.parse_args()
ic = "disc/ics/disc_d%g_w0.2.h5" % a.d
with h5py.File(os.path.join(HERE, "..", ic)) as h:
    T = float(h.attrs["period"]); mu = float(h.attrs["mu"]); om = float(h.attrs["omega"]); R = float(h.attrs["radius"])
c = math.sqrt(1 + 4/3*mu) + om*R
M = 32*math.ceil(T*c/(a.cfl*a.d)/32)          # steps per rotation, multiple of 32
dt = T/M
name = "run_%s_d%g" % (a.tag, a.d)
d = os.path.join(HERE, name); os.makedirs(os.path.join(d, "output"), exist_ok=True)
cfg = f"""; thirdorder study 2026-10-01, rotating disc (test plan A.1), build {a.build} (see ../makeBuilds.sh).
; d = {a.d}, omega = {om}, T_rot = {T!r}; fixed dt = T/{M} (c dt/d = {c*dt/a.d:.4f}, c = c_l + omega R), {a.rot} rotations,
; 32 dumps per rotation.
initFile {ic}
outDir disc/{name}/output
timeStep {dt!r}
timeEnd {a.rot*T + 0.5*dt!r}
h5DumpInterval {M//32 + 1}
periodicBoxLimits {{
    lowerX -3.0
    lowerY -3.0
    lowerZ 0.0
    upperX 3.0
    upperY 3.0
    upperZ 0.0
}}
kernelSize {math.sqrt(a.nnn/math.pi)*a.d!r}
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
materials {{
    mat0 {{ matId 0  rho0 1.0  A 1.0  B 0.0  u0 1.0  a 0.0  b 0.0  alpha 0.0  beta 0.0
           u_iv 0.05  u_cv 0.18  mu {mu!r} }}
}}
"""
open(os.path.join(d, "config.info"), "w").write(cfg)
open(os.path.join(d, "BUILD"), "w").write("build_%s\n" % a.build)
print("disc/" + name, "dt = %.6g, M = %d per rotation" % (dt, M))
