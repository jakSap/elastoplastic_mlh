"""Initial conditions for the multi-testcase-optimizing (mto) testplan.

Every generator writes the demonstrator's HDF5 IC format (x, v, m, u, materialId,
Sxx/Sxy/Syy) and stores its parameters as attributes, so a run directory is
self-describing. All bodies are square lattices, uniform density rho0 = 1,
u = 1, S = 0 (except the seeded disc), Hopkins-canonical Tillotson (P = A (rho/rho0 - 1), K = 1).

  rings   two annuli (r 3..4) head-on along x                 -> test T1
  blocks  two 0.5 x 0.5 squares, offset in y, head-on along x  -> test T2
  disc    rigidly rotating disc R = 1, in spinning equilibrium -> test A1
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


def circles(out, d, R=0.25, gap=0.12, yoff=0.25, v=0.25):
    """Two discs of radius R, gap apart in x, centres offset by yoff in y, closing at 2v (P2).

    Offset (glancing) plastic impact: the collision carries angular momentum, so the discs
    flatten, shear and spin up. Disc B is the point mirror of disc A, (x, y) -> (-x, -y), so
    the IC is exactly point-symmetric (lattice centred on each disc centre, as in rings).
    """
    n = int(np.ceil(R / d)) + 1
    i = (np.arange(-n, n) + 0.5) * d
    X, Y = np.meshgrid(i, i, indexing="ij")
    keep = np.hypot(X, Y) < R
    a = np.column_stack([X[keep], Y[keep]]) + [-(R + gap / 2), -yoff / 2]
    pos = np.vstack([a, -a])
    vel = np.zeros_like(pos)
    vel[: len(a), 0] = v
    vel[len(a):, 0] = -v
    return _write(out, pos, vel, d, dict(test="circles", R=R, gap=gap, yoff=yoff, velocity=v,
                                         t_contact=gap / (2 * v)))


def disc_equilibrium(Rb, omega, mu, dev_dim, K=1.0, rho0=RHO0, R_end=None):
    """Radial equilibrium of a disc spinning at omega, free surface at reference radius Rb.

    The code's own constitutive law, no linearisation: P = K (rho/rho0 - 1) (the Tillotson
    of the body tests, A = K, a = b = B = 0) and the hypoelastic deviator integrated along
    the coaxial spin-up path, S = 2 mu dev(ln V) (Hencky; Jaumann + rigid rotation keeps the
    principal axes e_r, e_theta corotating). dev_dim = 2: 2D-traceless (STRESS_2D_GIZMO 1),
    3: in-plane part of the 3D deviator (plane strain, STRESS_2D_GIZMO 0). With
    lr = dr/dR, lt = r/R, the Lagrangian radial balance is
        d sigma_rr/dR = -rho0 omega^2 R - lr (sigma_rr - sigma_tt)/r,  sigma_rr(Rb) = 0,
    shot from the centre (lr = lt = l0) for l0. Returns r(R) and sigma_rr, sigma_tt, S_rr,
    S_tt, P as functions of the reference radius R <= max(Rb, R_end) (R_end: outermost
    particle, which may lie beyond the area-equivalent surface Rb).
    """
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq

    c = 2.0 * mu / dev_dim

    def stress(lr, lt):
        P = K * (1.0 / (lr * lt) - 1.0)                      # rho/rho0 = 1/(lr lt)
        Srr = 2.0 * mu * np.log(lr) - c * np.log(lr * lt)
        Stt = 2.0 * mu * np.log(lt) - c * np.log(lr * lt)
        return P, Srr, Stt

    def rhs(R, y):
        r, lr = y
        lt = r / R
        dlt = (lr - lt) / R
        dsrr = -rho0 * omega ** 2 * R - lr * 2.0 * mu * np.log(lr / lt) / r
        ds_dlr = K / (lr * lr * lt) + (2.0 * mu - c) / lr      # d sigma_rr / d lr
        ds_dlt = K / (lr * lt * lt) - c / lt                   # d sigma_rr / d lt
        return [lr, (dsrr - ds_dlt * dlt) / ds_dlr]

    R0 = 1e-6 * Rb                                             # lr - lt = O(R^2) at the centre
    R1 = max(Rb, R_end or Rb)

    def shoot(l0):
        return solve_ivp(rhs, (R0, R1), [l0 * R0, l0], rtol=1e-12, atol=1e-14, dense_output=True)

    def srr_edge(l0):
        r, lr = shoot(l0).sol(Rb)
        P, Srr, _ = stress(lr, r / Rb)
        return -P + Srr

    l0 = brentq(srr_edge, 1.0, 1.5, xtol=1e-15)
    sol = shoot(l0).sol

    def at(R):
        R = np.maximum(np.asarray(R, dtype=float), R0)
        r, lr = sol(R)
        P, Srr, Stt = stress(lr, r / R)
        return dict(r=r, P=P, Srr=Srr, Stt=Stt, srr=-P + Srr, stt=-P + Stt)

    return at, l0


def disc(out, d, R=1.0, omega=0.2, mu=None, dev_dim=3, K=1.0):
    """Rigidly rotating disc, v = omega x r (L_z must hold).

    mu = None: unstressed lattice, S = 0 (breathes). Otherwise the disc starts in its
    continuum spinning equilibrium (disc_equilibrium; surface effects of the SPH density
    ignored): every lattice point is displaced radially r(R) -- the density, hence P, is
    set by the stretched lattice -- and carries S = S_rr e_r e_r + S_tt e_t e_t. The free
    surface is at the area-equivalent reference radius sqrt(N d^2 / pi).
    """
    n = int(np.ceil(R / d)) + 2
    i = np.arange(-n, n + 1) * d
    X, Y = np.meshgrid(i, i, indexing="ij")
    keep = X * X + Y * Y <= R * R
    pos = np.column_stack([X[keep], Y[keep]])
    attrs = dict(test="disc", radius=R, omega=omega, period=2 * np.pi / omega)
    S = None
    if mu is not None:
        Rb = np.sqrt(len(pos) * d * d / np.pi)
        R0 = np.hypot(pos[:, 0], pos[:, 1])
        eq, l0 = disc_equilibrium(Rb, omega, mu, dev_dim, K=K, R_end=R0.max())
        f = eq(R0)
        scale = np.where(R0 > 0, f["r"] / np.maximum(R0, 1e-300), l0)
        pos = pos * scale[:, None]
        r = np.hypot(pos[:, 0], pos[:, 1])
        cs = np.where(r > 0, pos[:, 0] / np.maximum(r, 1e-300), 1.0)
        sn = np.where(r > 0, pos[:, 1] / np.maximum(r, 1e-300), 0.0)
        Srr, Stt = f["Srr"], f["Stt"]
        S = dict(Sxx=Srr * cs * cs + Stt * sn * sn, Syy=Srr * sn * sn + Stt * cs * cs,
                 Sxy=(Srr - Stt) * cs * sn)
        edge = eq(Rb)
        attrs.update(seeded=1, mu=mu, dev_dim=dev_dim, K=K, R_ref=Rb, R_eq=float(edge["r"]),
                     lambda0=l0, P_centre=float(f["P"][np.argmin(R0)]))
    vel = np.column_stack([-omega * pos[:, 1], omega * pos[:, 0]])
    return _write(out, pos, vel, d, attrs, S=S)


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


def wave(out, n, mode, eos, mu, amplitude=1e-6, nnn=32, kernel=6, pmod=4.0 / 3.0):
    """Travelling wave along +x in the [0,1) x [0,1/4) box, n particles per wavelength."""
    if n % 4:
        raise ValueError("n must be a multiple of 4 (box is 4 x 1 cells)")
    wave_build((4, 1), n // 4, (1, 0), mode=mode, eos=eos, mu=mu, amplitude=amplitude,
               nnn=nnn, out=out, quiet=True, kernel=kernel, pmod=pmod)
    with h5py.File(out, "r") as f:
        return len(f["m"])
