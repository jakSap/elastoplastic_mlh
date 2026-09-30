"""Initial conditions for the plasticity tests (demonstrator .h5 layout: x, v, m, u, materialId
[, Sxx, Sxy, Syy]). One IC per (test case, rung), shared by every variant and by GIZMO
(converted, never regenerated). ICs are cached in <study>/ics/ and listed in IC_SOURCES.txt."""
import hashlib
import math
import os
import shutil

import h5py
import numpy as np


def _md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()[:12]


def _write(path, x, v, m, u, S=None, attrs=None):
    with h5py.File(path, "w") as f:
        f.create_dataset("x", data=np.asarray(x, dtype=np.float64))
        f.create_dataset("v", data=np.asarray(v, dtype=np.float64))
        f.create_dataset("m", data=np.asarray(m, dtype=np.float64))
        f.create_dataset("u", data=np.asarray(u, dtype=np.float64))
        f.create_dataset("materialId", data=np.zeros(len(m), dtype=np.int8))
        if S is not None:
            for k, a in zip(("Sxx", "Sxy", "Syy"), S):
                f.create_dataset(k, data=np.asarray(a, dtype=np.float64))
        for k, val in (attrs or {}).items():
            f.attrs[k] = val


def _lattice(x0, x1, y0, y1, d):
    nx, ny = int(round((x1 - x0) / d)), int(round((y1 - y0) / d))
    if abs(nx * d - (x1 - x0)) > 1e-9 * d or abs(ny * d - (y1 - y0)) > 1e-9 * d:
        raise ValueError("extent (%g, %g) is not a multiple of d = %g" % (x1 - x0, y1 - y0, d))
    X, Y = np.meshgrid(x0 + (np.arange(nx) + 0.5) * d, y0 + (np.arange(ny) + 0.5) * d, indexing="ij")
    return np.column_stack([X.ravel(), Y.ravel()])


# ------------------------------------------------------------------------ P0b block
def drive_matrix(drive, rate):
    """Velocity gradient L0 (v = L0 x) of the P0b drives; rate = engineering strain rate."""
    return {"uni": np.array([[-rate, 0.], [0., 0.]]),              # uniaxial compression
            "psh": np.array([[0., rate / 2], [rate / 2, 0.]]),     # pure shear, gamma_dot = rate
            "ssh": np.array([[0., rate], [0., 0.]])}[drive]        # simple shear v_x = rate y


def gen_block(path, d, mat, drive, rate, half):
    x = _lattice(-half, half, -half, half, d)
    L0 = drive_matrix(drive, rate)
    v = x @ L0.T
    N = len(x)
    _write(path, x, v, np.full(N, mat["rho0"] * d * d), np.full(N, mat["u0"]),
           attrs=dict(drive=drive, rate=rate, half=half, d=d))
    return N


# ---------------------------------------------------------------------- P2 / P3 / P4
def gen_slabs(path, d, mat, L, H, up):
    xl, xr = _lattice(-L, 0., 0., H, d), _lattice(0., L, 0., H, d)
    x = np.vstack([xl, xr])
    v = np.zeros_like(x)
    v[:len(xl), 0], v[len(xl):, 0] = up, -up
    N = len(x)
    _write(path, x, v, np.full(N, mat["rho0"] * d * d), np.full(N, mat["u0"]),
           attrs=dict(L=L, H=H, up=up, d=d, n_left=len(xl)))
    return N


def gen_flyer(path, d, mat, Lf, Lt, H, vf):
    xf, xt = _lattice(-Lf, 0., 0., H, d), _lattice(0., Lt, 0., H, d)
    x = np.vstack([xf, xt])
    v = np.zeros_like(x)
    v[:len(xf), 0] = vf
    N = len(x)
    _write(path, x, v, np.full(N, mat["rho0"] * d * d), np.full(N, mat["u0"]),
           attrs=dict(Lf=Lf, Lt=Lt, H=H, vf=vf, d=d, n_flyer=len(xf)))
    return N


def gen_bars(path, d, mat, L, W, v):
    """Two bars L x W, mirror-symmetric about x = 0, facing surfaces a gap d apart."""
    xr = _lattice(0.5 * d, 0.5 * d + L, -W / 2, W / 2, d)
    xl = xr.copy()
    xl[:, 0] *= -1
    x = np.vstack([xl, xr])
    vel = np.zeros_like(x)
    vel[:len(xl), 0], vel[len(xl):, 0] = v, -v
    N = len(x)
    _write(path, x, vel, np.full(N, mat["rho0"] * d * d), np.full(N, mat["u0"]),
           attrs=dict(L=L, W=W, v=v, d=d, n_left=len(xl)))
    return N


# ------------------------------------------------------------------ P1d true-2D waves
# Box [0,1) x [0,1/4), n x n/4 lattice, wave along +x, A = 1e-6, masses corrected for the
# kernel-sum density offset at N_NN (as testcases/soundwave_2d_elastic/generateIC.py).
def _lattice_kernel_sum(h_over_dx, R=12):
    i = np.arange(-R, R + 1)
    X, Y = np.meshgrid(i, i)
    u = np.hypot(X, Y) / h_over_dx
    t = np.clip(1.0 - u, 0.0, None)
    return (7.0 / np.pi / h_over_dx**2 * t**4 * (1.0 + 4.0 * u)).sum()


def kernel_offset(nnn):
    lo, hi = 1.0, 10.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if np.pi * mid**2 * _lattice_kernel_sum(mid) < nnn:
            lo = mid
        else:
            hi = mid
    return _lattice_kernel_sum(0.5 * (lo + hi)) - 1.0


def gen_wave(path, n, mat, mode, nnn, A=1e-6):
    """True-2D eigenmodes (SOLID_2D_MODEL 1):
       p: rho = 1 + A sin kx, v_x = c_l A sin kx, c_l = sqrt((K + mu)/rho0),
          S_ij = -2 mu ln(rho/rho0) (n_i n_j - delta_ij/2) -> S_xx = -S_yy = -mu ln(rho/rho0)
       s: v_y = c_s A sin kx, S_xy = -rho0 c_s v_y, rho uniform."""
    K, r0, mu = mat["A"], mat["rho0"], mat["mu"]
    k, dx = 2 * np.pi, 1.0 / n
    ny = n // 4
    q = (np.arange(n) + 0.5) * dx
    xs = q.copy()
    if mode == "p":
        for _ in range(60):                       # x + (A/k)(1 - cos kx) = q
            f = xs + A / k * (1 - np.cos(k * xs)) - q
            xs -= f / (1 + A * np.sin(k * xs))
    X, Y = np.meshgrid(xs, (np.arange(ny) + 0.5) * dx, indexing="ij")
    x = np.column_stack([X.ravel(), Y.ravel()])
    N = len(x)
    off = kernel_offset(nnn)
    m = np.full(N, 1.0 * 1.0 * 0.25 / N / (1 + off))
    v = np.zeros_like(x)
    if mode == "p":
        rho = 1.0 * (1 + A * np.sin(k * x[:, 0]))
        c = math.sqrt((K + mu) / r0)
        v[:, 0] = c * A * np.sin(k * x[:, 0])
        e = lambda r: K / r0 * np.log(r) + K / r        # isentropic u of P = K(r/r0 - 1)
        u = mat["u0"] + e(rho) - e(1.0)
        lnr = np.log(rho / r0)
        S = (-mu * lnr, np.zeros(N), mu * lnr)
    else:
        c = math.sqrt(mu / r0)
        v[:, 1] = c * A * np.sin(k * x[:, 0])
        u = np.full(N, mat["u0"])
        S = (np.zeros(N), -r0 * c * v[:, 1], np.zeros(N))
    _write(path, x, v, m, u, S=S, attrs=dict(mode=mode, n=n, amplitude=A, c=c, mu=mu, K0=K,
                                             kernel_offset=off, box_x=1.0, box_y=0.25,
                                             model="true 2D (SOLID_2D_MODEL 1)"))
    return N


# --------------------------------------------------------------------------- driver
def ic_name(test, T, rung):
    kind, a = T["ic"]
    if kind == "copy":
        return os.path.basename(a["src"])
    # generated ICs carry a hash of their arguments, so a changed manifest never reuses a stale file
    h = hashlib.md5(repr(sorted(a.items())).encode()).hexdigest()[:6]
    if kind == "block":
        return "block_%s_d%s_%s.h5" % (a["drive"], rung, h)
    if kind == "wave":
        return "wave%s_true2d_n%s_%s.h5" % (a["mode"], rung, h)
    return "%s_d%s_%s.h5" % (kind, rung, h)


def ensure_ic(S, test, rung):
    """Create (once) the IC for (test, rung). Returns (study-relative path, kernelSize, meta)."""
    T = S.M.TESTS[test]
    kind, a = T["ic"]
    name = ic_name(test, T, rung)
    path = S.p("ics", name)
    mat = {**S.M.MATERIAL, **T.get("material", {})}
    KS = math.sqrt(S.M.NNN / math.pi)
    ks = KS / float(rung) if T.get("ks_scale") == "n" else KS * float(rung)
    if not os.path.exists(path):
        if kind == "copy":
            shutil.copy(a["src"], path)
            src = "copy of %s" % a["src"]
        else:
            fn = {"block": gen_block, "slabs": gen_slabs, "flyer": gen_flyer, "bars": gen_bars}.get(kind)
            if kind == "wave":
                gen_wave(path, int(rung), mat, a["mode"], S.M.NNN)
            else:
                fn(path, float(rung), mat, **a)
            src = "ptlib/gen.py %s %s rung %s" % (kind, a, rung)
        with open(S.p("ics", "IC_SOURCES.txt"), "a") as f:
            f.write("%s  md5 %s  %s\n" % (name, _md5(path), src))
    with h5py.File(path, "r") as f:
        N = len(f["m"])
    return "ics/" + name, ks, dict(N=N, name=name)
