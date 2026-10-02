#!/usr/bin/env python3
"""Initial conditions for linear travelling waves on a periodic lattice, in 2D and 3D,
with an arbitrary (commensurate) wave-vector direction.

Shared by testcases/soundwave_2d, testcases/soundwave_2d_elastic, testcases/shearwave_2d
and testcases/wave_rotated_2d; their generateIC.py are thin wrappers that fix defaults and
file names. See testcases/wave_validation/PLAN.md for the study this serves.

Geometry
--------
The wave must be periodic in the box AND the box must be an integer number of lattice
spacings in every direction. With a box of cells*scale lattice cells per axis and mode
numbers modes = (p, q[, r]),

    L_i = cells_i * scale * d,   k = 2 pi (p/Lx, q/Ly[, r/Lz]),   tan(theta) = (q Lx)/(p Ly)

which is always rational: an exactly 60 deg wave vector does not exist on a periodic square
lattice. d is fixed by normalising Lx = 1. The resolution is ppw = lambda/d particles per
wavelength, lambda = 2 pi/|k|.

Examples
    cells=(4,1),  scale=n/4, modes=(1,0)  -> the legacy [0,1)x[0,1/4) box, wave along +x,
                                            n particles per wavelength (ppw = n)
    cells=(7,4),  scale=s,   modes=(1,1)  -> theta = 60.255 deg, ppw = 3.4729 s, N = 28 s^2
    cells=(7,4),  scale=s,   modes=(2,0)  -> its axis-aligned twin, ppw = 3.5 s, same N

Modes
-----
"p" - longitudinal (sound wave at mu = 0, elastic P-wave at mu > 0), direction n_hat:

    rho  = rho_bar (1 + A sin(k.x))
    v    = c_l A sin(k.x) n_hat,        c_l = sqrt((K + 4/3 mu)/rho_0)
    S_ij = -2 mu ln(rho/rho_0) (n_i n_j - delta_ij/3)
    u    = isentropic profile of the EOS

  The stress is not optional at mu > 0: with S = 0 the deviator references each particle's
  own initial density while the pressure references rho_0, the state is not an eigenmode and
  the wave splits (see the 2026-09-18 report, sec. 0).

"s" - shear, propagation n_hat, polarisation e_hat perpendicular to n_hat:

    rho  = rho_bar                      (a shear wave perturbs the density only at O(A^2))
    v    = c_s A sin(k.x) e_hat,        c_s = sqrt(mu/rho_0)
    S_ij = -rho_0 c_s (n_i e_j + e_i n_j) (c_s A sin(k.x))

  For n_hat = e_x, e_hat = e_y this is S_xy = -rho c_s v_y. The restoring force is 100 %
  deviatoric and, to leading order, there is no density perturbation, so EOS, kernel-sum
  density offset and tensile correction drop out of it.

Sign conventions checked against integrateStressTensor (Particles.cpp:1150ff) and
xSplitElasticHLLC (HLLC.cpp:640). The code's 2D deviator uses the 3D trace convention with
dvz_dz = 0 (visc_xx = (4/3) mu dvx_dx - (2/3) mu dvy_dy = 2 mu (eps_xx - tr eps/3)), so
S_xx + S_yy = -S_zz and the formulae above hold unchanged in 2D and 3D.

Lattice displacement
--------------------
Both modes are one-dimensional in the wave coordinate sigma = n_hat.x, so the density
profile of the "p" mode is imposed exactly by displacing each particle along n_hat:
solve sigma + (A/k)(1 - cos k sigma) = sigma_q (Newton, to machine precision) and move the
particle by (sigma - sigma_q) n_hat. The Jacobian is 1/(1 + A sin k sigma), i.e. the density
is exact, and periodicity survives because k is a reciprocal-lattice vector of the box.

Kernel-sum density offset
-------------------------
rho_i = m_i sum_j W_ij on a lattice sits a constant above m/d^DIM that depends only on N_NN.
h/d is solved from the solver's own definition (V sum_j W = N_NN with V = pi h^2 in 2D and
(4/3) pi h^3 in 3D, Particles.cpp:2616-2619), then offset = d^DIM sum_j W - 1 and the masses
are divided by 1 + offset, so <rho> = rho_bar, P_bar = 0 and the wave has genuine tensile
half-cycles. 2D at N_NN = 32: +4.237934e-3 (reproduces the solver to 3e-11).

Usage as a CLI (mostly for spot checks; the testcases call make_ic()):
    ./waveic.py --cells 7 4 --scale 18 --modes 1 1 --mode p --eos tillotson --mu 0.4615384615
    ./waveic.py --cells 4 1 --scale 16 --modes 1 0 --mode s --mu 0.4615384615 -o shear_N64.h5
"""

import numpy as np
import h5py

# --- material normalisations, identical to testcases/soundwave_2d -------------------------

RHO_BAR = 1.0

GAMMA = 5.0 / 3.0                       # ideal gas (Hopkins)
P_BAR_IDEAL = 3.0 / 5.0
MURN_K0, MURN_N, MURN_RHO0 = 1.0, 2.0, 1.0      # Murnaghan (repo standard, rip_test)
TILL_RHO0, TILL_A = 1.0, 1.0                    # Tillotson, Hopkins-canonical
U_BAR_SOLID = 1.0                               # u > 0 well clear of 0 (09-25: solids at u_bar = 1)

EOS_NAMES = ("ideal", "murnaghan", "tillotson")
# Poisson's ratio 0.3 at K = 1: mu = 3 K (1 - 2 nu)/(2 (1 + nu)) = 6/13
MU_NU03 = 3.0 * TILL_A * (1.0 - 2.0 * 0.3) / (2.0 * (1.0 + 0.3))


def eos_state(eos, rho):
    """Pressure, sound speed and isentropic u(rho), relative to RHO_BAR."""
    if eos == "ideal":
        u_bar = P_BAR_IDEAL / ((GAMMA - 1.0) * RHO_BAR)
        u = u_bar * (rho / RHO_BAR) ** (GAMMA - 1.0)
        P = (GAMMA - 1.0) * rho * u
        c = np.sqrt(GAMMA * P / rho)
    elif eos == "murnaghan":
        K0, n, r0 = MURN_K0, MURN_N, MURN_RHO0

        def e(r):  # int P/rho^2 drho, as EquationOfState::EOSInternalEnergy
            return (K0 / (n * (n - 1.0) * r0**n) * (r**(n - 1.0) - r0**(n - 1.0))
                    + K0 / n * (1.0 / r - 1.0 / r0))
        u = U_BAR_SOLID + e(rho) - e(RHO_BAR)
        P = K0 / n * ((rho / r0) ** n - 1.0)
        c = np.sqrt(K0 / r0 * (rho / r0) ** (n - 1.0))
    elif eos == "tillotson":
        A, r0 = TILL_A, TILL_RHO0

        def e(r):  # int A (r/r0 - 1)/r^2 dr
            return A / r0 * np.log(r) + A / r
        u = U_BAR_SOLID + e(rho) - e(RHO_BAR)
        P = A * (rho / r0 - 1.0)
        c = np.sqrt(A / r0) * np.ones_like(rho)
    else:
        raise ValueError("unknown EOS {!r}, expected one of {}".format(eos, EOS_NAMES))
    return P, c, u


def bulk_modulus(eos):
    """K = rho dP/drho at rho = RHO_BAR, the K entering c_l."""
    _, c, _ = eos_state(eos, np.array([RHO_BAR]))
    return float(c[0]) ** 2 * RHO_BAR


def longitudinal_speed(eos, mu):
    """c_l = sqrt((K + 4/3 mu)/rho_0), what xSplitElasticHLLC uses (HLLC.cpp:640)."""
    return np.sqrt((bulk_modulus(eos) + 4.0 / 3.0 * mu) / RHO_BAR)


def shear_speed(mu):
    """c_s = sqrt(mu/rho_0)."""
    return np.sqrt(mu / RHO_BAR)


# --- kernel-sum density offset ------------------------------------------------------------

def lattice_kernel_sum(h_over_d, dim=2, R=12):
    """d^DIM * sum_j W(r_j, h) over the infinite cubic lattice (Wendland C2, GIZMO norms)."""
    i = np.arange(-R, R + 1)
    if dim == 2:
        X, Y = np.meshgrid(i, i)
        u = np.hypot(X, Y) / h_over_d
        t = np.clip(1.0 - u, 0.0, None)
        return (7.0 / np.pi / h_over_d**2 * t**4 * (1.0 + 4.0 * u)).sum()
    X, Y, Z = np.meshgrid(i, i, i)
    u = np.sqrt(X * X + Y * Y + Z * Z) / h_over_d
    t = np.clip(1.0 - u, 0.0, None)
    return (21.0 / (2.0 * np.pi) / h_over_d**3 * t**4 * (1.0 + 4.0 * u)).sum()


def kernel_offset(nnn, dim=2):
    """Kernel-sum density offset for target neighbour number nnn (bisection on h/d).

    The solver's neighbour count is V(h) sum_j W with V = pi h^2 (2D) / (4/3) pi h^3 (3D),
    Particles.cpp:2616-2619.
    """
    vcoef = np.pi if dim == 2 else 4.0 / 3.0 * np.pi
    lo, hi = 1.0, 10.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if vcoef * mid**dim * lattice_kernel_sum(mid, dim) < nnn:
            lo = mid
        else:
            hi = mid
    return lattice_kernel_sum(0.5 * (lo + hi), dim) - 1.0


# --- geometry -----------------------------------------------------------------------------

class Geometry:
    """Periodic box of cells[i]*scale lattice cells per axis, Lx normalised to 1."""

    def __init__(self, cells, scale):
        self.cells = tuple(int(c) for c in cells)
        self.dim = len(self.cells)
        if self.dim not in (2, 3):
            raise ValueError("need 2 or 3 cell counts, got {}".format(self.cells))
        if any(c < 1 for c in self.cells):
            raise ValueError("cell counts must be positive: {}".format(self.cells))
        self.scale = int(scale)
        self.counts = tuple(c * self.scale for c in self.cells)   # particles per axis
        self.d = 1.0 / self.counts[0]                             # lattice spacing, Lx = 1
        self.L = tuple(n * self.d for n in self.counts)
        self.N = int(np.prod(self.counts))
        self.volume = float(np.prod(self.L))

    def lattice(self):
        """Cell-centred lattice positions, x-major (indexing='ij', C order)."""
        axes = [(np.arange(n) + 0.5) * self.d for n in self.counts]
        grids = np.meshgrid(*axes, indexing="ij")
        return np.column_stack([g.ravel() for g in grids])

    def __repr__(self):
        return "Geometry(cells={}, scale={}, counts={}, L={}, N={})".format(
            self.cells, self.scale, self.counts, tuple(round(l, 12) for l in self.L), self.N)


class Wave:
    """Reciprocal-lattice wave vector k = 2 pi (modes_i / L_i) on a Geometry."""

    def __init__(self, geometry, modes):
        self.geo = geometry
        self.modes = tuple(int(m) for m in modes)
        if len(self.modes) != geometry.dim:
            raise ValueError("need {} mode numbers, got {}".format(geometry.dim, self.modes))
        if all(m == 0 for m in self.modes):
            raise ValueError("at least one mode number must be non-zero")
        self.k = 2.0 * np.pi * np.array([m / L for m, L in zip(self.modes, geometry.L)])
        self.kmag = float(np.linalg.norm(self.k))
        self.n_hat = self.k / self.kmag
        self.wavelength = 2.0 * np.pi / self.kmag
        self.ppw = self.wavelength / geometry.d

    @property
    def theta(self):
        """Angle to the x axis in degrees (in the xy plane)."""
        return float(np.degrees(np.arctan2(self.n_hat[1], self.n_hat[0])))

    @property
    def phi(self):
        """Polar angle to the z axis in degrees (3D only, else 90)."""
        if self.geo.dim == 2:
            return 90.0
        return float(np.degrees(np.arccos(self.n_hat[2])))

    def phase(self, pos):
        """k.x for every particle."""
        return pos @ self.k

    def perpendicular(self, hint=None):
        """A unit vector perpendicular to n_hat (the default shear polarisation).

        2D: the left normal (-n_y, n_x). 3D: `hint` (default: the axis least aligned with
        n_hat) projected out of n_hat and normalised -- recorded in the IC attributes.
        """
        if self.geo.dim == 2:
            return np.array([-self.n_hat[1], self.n_hat[0]])
        if hint is None:
            hint = np.zeros(3)
            hint[int(np.argmin(np.abs(self.n_hat)))] = 1.0
        hint = np.asarray(hint, dtype=float)
        e = hint - np.dot(hint, self.n_hat) * self.n_hat
        norm = np.linalg.norm(e)
        if norm < 1e-12:
            raise ValueError("polarisation hint {} is parallel to n_hat".format(hint))
        return e / norm

    def __repr__(self):
        return ("Wave(modes={}, lambda={:.9g}, ppw={:.6g}, theta={:.4f} deg"
                .format(self.modes, self.wavelength, self.ppw, self.theta)
                + (", phi={:.4f} deg)".format(self.phi) if self.geo.dim == 3 else ")"))


def lagrangian_to_euler(q, A, k):
    """Solve sigma + (A/k)(1 - cos k sigma) = q for sigma (Newton, machine precision)."""
    x = q.copy()
    for _ in range(50):
        f = x + A / k * (1.0 - np.cos(k * x)) - q
        dx = f / (1.0 + A * np.sin(k * x))
        x -= dx
        if np.max(np.abs(dx)) < 1e-17:
            break
    return x


# --- the IC -------------------------------------------------------------------------------

def make_ic(geometry, wave, mode="p", eos="tillotson", mu=0.0, amplitude=1e-6,
            nnn=32, polarisation=None, out=None, quiet=False):
    """Build a travelling-wave IC and write it to HDF5. Returns the file name.

    mode  "p" longitudinal (sound / P-wave), "s" shear (needs mu > 0)
    nnn   target neighbour number for the mass correction; None = nominal masses
    """
    dim, A, k = geometry.dim, float(amplitude), wave.kmag
    if mode not in ("p", "s"):
        raise ValueError("mode must be 'p' or 's', got {!r}".format(mode))
    if mode == "s" and mu <= 0.0:
        raise ValueError("a shear wave needs mu > 0")

    pos = geometry.lattice()
    e_hat = wave.perpendicular(polarisation) if mode == "s" else None

    if mode == "p":
        # displace along n_hat so that the kernel-sum density reproduces rho(sigma)
        sigma_q = pos @ wave.n_hat
        sigma = lagrangian_to_euler(sigma_q, A, k)
        pos = pos + np.outer(sigma - sigma_q, wave.n_hat)
        phase = sigma                      # = n_hat.x after the displacement
        # An oblique displacement can push a particle out of the box; the configuration is
        # periodic, so wrapping is exact. (No-op for the axis-aligned legacy cases.)
        pos = np.mod(pos, np.array(geometry.L))
        rho = RHO_BAR * (1.0 + A * np.sin(k * phase))
        P, _, u = eos_state(eos, rho)
        c_wave = longitudinal_speed(eos, mu)
        vel = np.outer(c_wave * A * np.sin(k * phase), wave.n_hat)
        # S_ij = -2 mu ln(rho/rho_0) (n_i n_j - delta_ij/3)
        lnr = np.log(rho / TILL_RHO0)
        shape = np.outer(wave.n_hat, wave.n_hat) - np.eye(dim) / 3.0
        S = (-2.0 * mu * lnr)[:, None, None] * shape[None, :, :]
    else:
        rho = np.full(geometry.N, RHO_BAR)
        P, _, u = eos_state(eos, rho)
        c_wave = shear_speed(mu)
        phase = pos @ wave.n_hat
        v_amp = c_wave * A * np.sin(k * phase)
        vel = np.outer(v_amp, e_hat)
        # S_ij = -rho_0 c_s (n_i e_j + e_i n_j) v
        shape = np.outer(wave.n_hat, e_hat) + np.outer(e_hat, wave.n_hat)
        S = (-RHO_BAR * c_wave * v_amp)[:, None, None] * shape[None, :, :]

    offset = 0.0 if nnn is None else kernel_offset(nnn, dim)
    m = np.full(geometry.N, RHO_BAR * geometry.volume / geometry.N / (1.0 + offset))
    matId = np.zeros(geometry.N, dtype=np.int8)
    P_bar = float(eos_state(eos, np.array([RHO_BAR]))[0][0])
    period = wave.wavelength / c_wave

    fname = out or "wave_{}_{}_s{}.h5".format(mode, eos, geometry.scale)
    with h5py.File(fname, "w") as f:
        f.create_dataset("x", data=pos)
        f.create_dataset("v", data=vel)
        f.create_dataset("m", data=m)
        f.create_dataset("u", data=u)
        f.create_dataset("materialId", data=matId)
        # same dataset names as dump2file, so a snapshot can be fed back in (fd88e21)
        names = (("Sxx", 0, 0), ("Sxy", 0, 1), ("Syy", 1, 1)) if dim == 2 else \
                (("Sxx", 0, 0), ("Sxy", 0, 1), ("Sxz", 0, 2),
                 ("Syy", 1, 1), ("Syz", 1, 2), ("Szz", 2, 2))
        for name, a, b in names:
            f.create_dataset(name, data=np.ascontiguousarray(S[:, a, b]))
        for key, val in (
                ("dim", dim), ("mode", mode), ("eos", eos),
                ("cells", np.array(geometry.cells)), ("scale", geometry.scale),
                ("counts", np.array(geometry.counts)), ("d", geometry.d),
                ("box", np.array(geometry.L)), ("modes", np.array(wave.modes)),
                ("k", wave.k), ("n_hat", wave.n_hat),
                ("e_hat", e_hat if e_hat is not None else np.zeros(dim)),
                ("theta", wave.theta), ("phi", wave.phi),
                ("wavelength", wave.wavelength), ("ppw", wave.ppw),
                ("amplitude", A), ("rho_bar", RHO_BAR), ("P_bar", P_bar),
                ("mu", mu), ("K0", bulk_modulus(eos)),
                ("c_l", longitudinal_speed(eos, mu)), ("c_s", shear_speed(mu)),
                ("c_wave", c_wave), ("period", period),
                ("kernel_offset", offset), ("correct_nnn", 0 if nnn is None else nnn)):
            f.attrs[key] = val

    if not quiet:
        print("Wrote {}: {} = {} particles, d = {:g}, lambda = {:.9g} ({:.4g} ppw), "
              "theta = {:.4f} deg, mode {}, c = {:.9f}, T = {:.9f}, "
              "|S|max = {:.4e}, mass factor 1/(1 {:+.8e})".format(
                  fname, "x".join(str(c) for c in geometry.counts), geometry.N, geometry.d,
                  wave.wavelength, wave.ppw, wave.theta, mode, c_wave, period,
                  np.abs(S).max(), offset))
    return fname


def build(cells, scale, modes, **kwargs):
    """Geometry + Wave + make_ic in one call."""
    geo = Geometry(cells, scale)
    wave = Wave(geo, modes)
    return make_ic(geo, wave, **kwargs)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", type=int, nargs="+", required=True,
                    help="box aspect in lattice cells, e.g. 7 4 (2D) or 3 2 2 (3D)")
    ap.add_argument("--scale", type=int, required=True, help="lattice scale s")
    ap.add_argument("--modes", type=int, nargs="+", required=True,
                    help="wavelengths along each box axis, e.g. 1 1")
    ap.add_argument("--mode", choices=("p", "s"), default="p", help="wave type")
    ap.add_argument("-e", "--eos", choices=EOS_NAMES, default="tillotson")
    ap.add_argument("--mu", type=float, default=0.0, help="shear modulus (default 0)")
    ap.add_argument("-A", "--amplitude", type=float, default=1e-6)
    ap.add_argument("--correct-nnn", type=int, default=32, metavar="N_NN")
    ap.add_argument("--no-correct", action="store_true", help="nominal masses")
    ap.add_argument("--polarisation", type=float, nargs="+", default=None,
                    help="3D shear polarisation hint (projected onto the plane of n_hat)")
    ap.add_argument("-o", "--out")
    args = ap.parse_args()

    build(args.cells, args.scale, args.modes, mode=args.mode, eos=args.eos, mu=args.mu,
          amplitude=args.amplitude, nnn=None if args.no_correct else args.correct_nnn,
          polarisation=args.polarisation, out=args.out)
