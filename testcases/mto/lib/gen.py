"""Initial conditions for the multi-testcase-optimizing (mto) testplan.

Every generator writes the demonstrator's HDF5 IC format (x, v, m, u, materialId,
Sxx/Sxy/Syy) and stores its parameters as attributes, so a run directory is
self-describing. All bodies are square lattices, uniform density rho0 = 1,
u = 1, S = 0, Hopkins-canonical Tillotson (P = A (rho/rho0 - 1), K = 1).

  rings   two annuli (r 3..4) head-on along x                 -> test T1
  blocks  two 0.5 x 0.5 squares, offset in y, head-on along x  -> test T2
  disc    rigidly rotating disc R = 1                          -> test A1
  rod     two rods head-on along x (Taylor impact)             -> test P1
  wave    periodic travelling wave (waveic.py)                 -> tests W1-W3

The head-on pairs are mirror-symmetric about x = 0, so each is equivalent to one
body hitting a rigid, frictionless wall; a periodic one-body variant would need
the same particle number (a body cannot collide with its own periodic image).
"""

import numpy as np
import h5py

from waveic import build as wave_build

RHO0 = 1.0
U_BODY = 1.0


def _lattice(xmin, xmax, ymin, ymax, d):
    """Cell-centred square lattice filling [xmin,xmax] x [ymin,ymax]."""
    nx = int(round((xmax - xmin) / d))
    ny = int(round((ymax - ymin) / d))
    xs = xmin + (np.arange(nx) + 0.5) * d
    ys = ymin + (np.arange(ny) + 0.5) * d
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    return np.column_stack([X.ravel(), Y.ravel()])


def _write(out, pos, vel, d, attrs, S=None):
    N = len(pos)
    with h5py.File(out, "w") as f:
        f.create_dataset("x", data=pos)
        f.create_dataset("v", data=vel)
        f.create_dataset("m", data=np.full(N, RHO0 * d * d))
        f.create_dataset("u", data=np.full(N, U_BODY))
        f.create_dataset("materialId", data=np.zeros(N, dtype=np.int8))
        for name in ("Sxx", "Sxy", "Syy"):
            f.create_dataset(name, data=np.zeros(N) if S is None else S[name])
        attrs = dict(attrs, N=N, spacing=d, rho0=RHO0, u=U_BODY)
        for k, v in attrs.items():
            f.attrs[k] = v
    return N


def rings(out, d, r_in=3.0, r_out=4.0, shift=5.0, v=0.059):
    """Two annuli centred at (-shift, 0) and (+shift, 0), moving towards each other.

    Lattice centred on each ring centre (same construction as
    colliding_rings/generateIC.py), so both rings are exact mirror images.
    """
    n = int(np.ceil(r_out / d)) + 1
    i = (np.arange(-n, n) + 0.5) * d
    X, Y = np.meshgrid(i, i, indexing="ij")
    r = np.hypot(X, Y)
    keep = (r >= r_in) & (r <= r_out)
    ring = np.column_stack([X[keep], Y[keep]])
    pos = np.vstack([ring + [-shift, 0.0], ring + [shift, 0.0]])
    vel = np.zeros_like(pos)
    vel[: len(ring), 0] = v
    vel[len(ring):, 0] = -v
    return _write(out, pos, vel, d, dict(test="rings", r_in=r_in, r_out=r_out,
                                         shift=shift, velocity=v,
                                         t_contact=(2 * shift - 2 * r_out) / (2 * v)))


def blocks(out, d, size=0.5, gap=0.1, yoff=0.25, v=0.059):
    """Two size x size squares, gap apart in x, offset by yoff in y, closing at 2v.

    Offset contact: the collision carries angular momentum and drives a corner
    into the other block's face (the 09-26 T.2/A.2 setup, ic_blocks.py -v 0.059
    -g 0.1 -y 0.25).
    """
    a = _lattice(-gap / 2 - size, -gap / 2, -yoff / 2 - size / 2, -yoff / 2 + size / 2, d)
    b = _lattice(gap / 2, gap / 2 + size, yoff / 2 - size / 2, yoff / 2 + size / 2, d)
    pos = np.vstack([a, b])
    vel = np.zeros_like(pos)
    vel[: len(a), 0] = v
    vel[len(a):, 0] = -v
    return _write(out, pos, vel, d, dict(test="blocks", size=size, gap=gap, yoff=yoff,
                                         velocity=v, t_contact=gap / (2 * v)))


def disc(out, d, R=1.0, omega=0.2):
    """Rigidly rotating disc, v = omega x r, S = 0 (breathes; L_z must hold)."""
    n = int(np.ceil(R / d)) + 2
    i = np.arange(-n, n + 1) * d
    X, Y = np.meshgrid(i, i, indexing="ij")
    keep = X * X + Y * Y <= R * R
    pos = np.column_stack([X[keep], Y[keep]])
    vel = np.column_stack([-omega * pos[:, 1], omega * pos[:, 0]])
    return _write(out, pos, vel, d, dict(test="disc", radius=R, omega=omega,
                                         period=2 * np.pi / omega))


def rod(out, d, length=1.0, width=0.25, gap=0.1, v=0.2):
    """Two rods along x, mirror-symmetric about x = 0, closing at 2v (Taylor impact)."""
    a = _lattice(-gap / 2 - length, -gap / 2, -width / 2, width / 2, d)
    b = _lattice(gap / 2, gap / 2 + length, -width / 2, width / 2, d)
    pos = np.vstack([a, b])
    vel = np.zeros_like(pos)
    vel[: len(a), 0] = v
    vel[len(a):, 0] = -v
    return _write(out, pos, vel, d, dict(test="rod", length=length, width=width, gap=gap,
                                         velocity=v, t_contact=gap / (2 * v)))


def wave(out, n, mode, eos, mu, amplitude=1e-6):
    """Travelling wave along +x in the [0,1) x [0,1/4) box, n particles per wavelength."""
    if n % 4:
        raise ValueError("n must be a multiple of 4 (box is 4 x 1 cells)")
    wave_build((4, 1), n // 4, (1, 0), mode=mode, eos=eos, mu=mu, amplitude=amplitude,
               nnn=32, out=out, quiet=True)
    with h5py.File(out, "r") as f:
        return len(f["m"])
