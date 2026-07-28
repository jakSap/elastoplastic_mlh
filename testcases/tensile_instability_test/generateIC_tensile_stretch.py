#!/usr/bin/env python3
"""
Homologous tensile-stretch IC for the newMFM demonstrator.

Unlike generateIC_tensile.py, which imposes a *step* velocity field
(vx = -v for x<0, +v for x>=0), this script imposes a *linear* (homologous)
velocity ramp

    vx = strain_rate * x,   vy = 0.

Why the ramp and not the step:

  The step field has div(v) = 2v * delta(x): all the strain is concentrated in
  a single plane at x=0, the bulk has div(v)=0, so the interior never goes into
  tension -- it simply translates rigidly and the bar tears at the center.
  (The step is a perfectly well-posed Riemann problem, but its exact solution
  is a centered double-rarefaction / cavitation: separation, not stretching.)

  The linear ramp has div(v) = strain_rate, *uniform* everywhere. Every face
  sees a tiny diverging jump ~ strain_rate*dx, the whole bar expands
  homologously, rho drops below rho0 everywhere, and the Murnaghan EOS

      P = K0/n * ((rho/rho0)^n - 1)

  returns P < 0 (uniform, sustained tension) throughout the bulk. That is the
  state in which the SPH/meshless tensile (pairing) instability can actually
  develop, so the cubic-spline vs Wendland-C2 comparison becomes meaningful.

IMPORTANT -- boundary conditions:

  A linear ramp vx = strain_rate*x is INCOMPATIBLE with x-periodicity (it has a
  jump of 2*strain_rate*L across the seam). Run this IC with a FREE-ended bar:
  set PERIODIC_BOUNDARIES 0 (or at least non-periodic in x) in parameter.h.
  The ends pull apart, the center stays put, the bulk stretches uniformly.

Optional --prestretch:
  Seed the lattice already expanded by a factor (1 + strain) with ZERO velocity,
  i.e. a static, uniformly pre-tensioned bar. This is the cleanest pure
  pairing-instability probe (no kinematic transient at all). Combine with a
  small/zero --strain-rate to hold the tension roughly static.

Single material (matId=0), normalized units (rho0 = K0 = 1, c = 1 at rho0).
u is set to 0: with a=b=alpha=beta=0, u is inert for pressure (P depends only
on rho), and u=0 keeps HLLC energy fluxes minimal.

Run:  python3 generateIC_tensile_stretch.py [--strain-rate 0.1] [--delta-p 0.05] [--prestretch 0.0] [--plot]
"""

import os
import argparse
import numpy as np
import h5py

# ----------------------------------------------------------------------------
# Parameters (normalized units; must match config materials + box)
# ----------------------------------------------------------------------------
x_lo, x_hi = -8.0, 8.0  # bar extent in x (FREE ends -- not periodic in x)
y_lo, y_hi = -1.5, 1.5    # transverse extent
dp_default = 0.05          # default particle spacing
sr_default = 0.1           # default strain rate; |vx| at the ends = sr * x_hi
rho0 = 1.0                 # reference density (P = 0 at rho0)
u0   = 1.0                 # Tillotson reference specific energy (u0 in config)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--strain-rate", type=float, default=sr_default,
                    help="homologous strain rate a in vx = a*x (default: %(default)s)")
    ap.add_argument("--delta-p", type=float, default=dp_default,
                    help="particle spacing (default: %(default)s)")
    ap.add_argument("--prestretch", type=float, default=0.0,
                    help="seed lattice expanded by (1+strain) for static "
                         "pre-tension (default: %(default)s -> none)")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()
    a = args.strain_rate
    delta_p = args.delta_p
    pre = args.prestretch

    out_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            f"tensile_stretch_ic_sr{a:g}_dp{delta_p:g}.h5")

    # Cell-centred lattice on the reference (unstretched) box.
    nx = int(round((x_hi - x_lo) / delta_p))
    ny = int(round((y_hi - y_lo) / delta_p))
    xs = x_lo + (np.arange(nx) + 0.5) * delta_p
    ys = y_lo + (np.arange(ny) + 0.5) * delta_p
    gx, gy = np.meshgrid(xs, ys)
    gx = gx.ravel()
    gy = gy.ravel()
    N = gx.size

    # Mass is fixed by the *reference* spacing so that, on the reference lattice,
    # the SPH density equals rho0 (P = 0). Any expansion then lowers rho below
    # rho0 and produces tension.
    mass = delta_p ** 2 * rho0
    m = np.ones(N) * mass

    # Optional static pre-stretch: expand positions in x by (1 + pre). The
    # spacing grows to delta_p*(1+pre), so the SPH density drops to
    # rho0/(1+pre) < rho0 -> uniform tension already at t=0.
    if pre != 0.0:
        gx = gx * (1.0 + pre)

    x = np.zeros((N, 2))
    x[:, 0] = gx
    x[:, 1] = gy

    # Homologous velocity ramp: uniform div(v) = a -> uniform stretching.
    vel = np.zeros((N, 2))
    vel[:, 0] = a * gx

    u = np.ones(N)             # a=b=alpha=beta=0: u inert for pressure; 0 keeps HLLC energy fluxes minimal
    materialId = np.zeros(N, dtype=np.int8)

    with h5py.File(out_file, "w") as h5f:
        h5f.create_dataset("x", data=x)
        h5f.create_dataset("v", data=vel)
        h5f.create_dataset("m", data=m)
        h5f.create_dataset("u", data=u)
        h5f.create_dataset("materialId", data=materialId)
        h5f.create_dataset("time", data=0.0 * m)

    rho_seed = rho0 / (1.0 + pre) if pre != 0.0 else rho0
    print(f"N = {N}  ({nx} x {ny}),  delta_p = {delta_p}")
    print(f"bar x in [{x_lo},{x_hi}] (FREE ends), y in [{y_lo},{y_hi}]")
    print(f"strain_rate a = {a}  ->  div(v) = {a} uniform,  |vx|(ends) = {a*x_hi}")
    print(f"prestretch = {pre}  ->  seed rho approx {rho_seed:.4f}  "
          f"(P {'<' if pre != 0.0 else '='} 0 at t=0)")
    print(f"mass/particle = {mass}")
    print(f"Saved {out_file}")
    print("NOTE: run with PERIODIC_BOUNDARIES 0 (free ends); a linear ramp is "
          "incompatible with x-periodicity.")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 2), dpi=130)
        sc = ax.scatter(x[:, 0], x[:, 1], s=1, c=vel[:, 0], cmap="coolwarm")
        ax.set_aspect("equal")
        fig.colorbar(sc, label="vx")
        ax.set_title(f"homologous tensile-stretch IC  N={N}  a={a}")
        png = out_file.replace(".h5", ".png")
        plt.savefig(png, bbox_inches="tight")
        print("Saved", png)


if __name__ == "__main__":
    main()
