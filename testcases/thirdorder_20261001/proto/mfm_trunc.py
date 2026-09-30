#!/usr/bin/env python3
"""Truncation-error prototype of the demonstrator's periodic MFM right-hand side (2D, ideal gas, mu = 0).

Evaluates the semi-discrete rates on the EXACT travelling-wave state at t = 0 (point values at the displaced
lattice positions of generateIC.py) and compares them with the exact rates:
    dv_x/dt   = -c^2 A k cos(k x)            (Lagrangian, linear)
    du/dt     = -(P/rho) c A k cos(k x)
    dln rho/dt= -c A k cos(k x)
Mirrors the demonstrator (Particles.cpp, periodic path): Wendland C2, h from N_NN = pi h^2 omega (self term
included), psi_j(x_i) = W(r_ij, h_i)/omega_i, first-order psiTilde for the faces,
A_ij = V_i psiTilde_j(x_i) - V_j psiTilde_i(x_j) with V = m/rho (rho = explicit density or kernel density),
midpoint quadrature point, face frame (v_i + v_j)/2, reconstructed (rho, P, v); linearised acoustic Riemann
solver (the HLLC star state in the small-amplitude limit), MFM fluxes F_mom = P* A_ij, F_E = P* v*_n |A_ij|.
Explicit density rate = -div v from the (higher-order) WLS gradient, as integrateExplicitVolumeHalfStep.

Options select the ingredients under test (reconstruction order + quadratic coefficient, gradient order,
which density sets V, background-pressure handling)."""
import numpy as np
from scipy.spatial import cKDTree

LX, LY = 1.0, 0.25
GAMMA = 5/3
A_AMP = 1e-6
K = 2*np.pi


def W(r, h):
    q = r/h
    t = np.clip(1.0 - q, 0.0, None)
    return 7.0/(np.pi*h*h)*t**4*(1.0 + 4.0*q)


def lattice(n, A=A_AMP):
    dx = LX/n; ny = int(round(LY/dx))
    q = (np.arange(n) + .5)*dx
    x = q.copy()
    for _ in range(60):
        f = x + A/K*(1 - np.cos(K*x)) - q
        x -= f/(1 + A*np.sin(K*x))
    ys = (np.arange(ny) + .5)*dx
    X, Y = np.meshgrid(x, ys, indexing="ij")
    return np.column_stack([X.ravel(), Y.ravel()]), dx


class Geometry:
    def __init__(self, pos, dx, nnn=32):
        self.pos = pos; N = len(pos); self.N = N
        h0 = np.sqrt(nnn/np.pi)*dx
        R = 1.3*h0
        tree = cKDTree(pos, boxsize=[LX, LY])
        pr = tree.query_pairs(R, output_type="ndarray")
        i = np.concatenate([pr[:, 0], pr[:, 1]]); j = np.concatenate([pr[:, 1], pr[:, 0]])
        d = pos[j] - pos[i]
        d[:, 0] -= LX*np.round(d[:, 0]/LX); d[:, 1] -= LY*np.round(d[:, 1]/LY)
        self.np_ = len(pr)
        self.i, self.j, self.d = i, j, d                   # d = x_j - x_i
        self.r = np.hypot(d[:, 0], d[:, 1])
        # h: pi h^2 (W(0,h) + sum_j W(r_ij, h)) = nnn, bisection per particle
        lo = np.full(N, 0.5*h0); hi = np.full(N, 1.29*h0)
        for _ in range(80):
            h = .5*(lo + hi)
            om = W(0., h) + np.bincount(i, W(self.r, h[i]), N)
            f = np.pi*h*h*om - nnn
            lo = np.where(f < 0, h, lo); hi = np.where(f >= 0, h, hi)
        self.h = .5*(lo + hi)
        assert self.h.max() < R
        self.w = W(self.r, self.h[i])                      # W(r_ij, h_i) (gather)
        self.omega = W(0., self.h) + np.bincount(i, self.w, N)
        # first-order psiTilde_j(x_i) = E_i^-1 (x_j - x_i) W_ij/omega_i
        wn = self.w/self.omega[i]
        E = np.zeros((N, 2, 2))
        for a in range(2):
            for b in range(2):
                E[:, a, b] = np.bincount(i, d[:, a]*d[:, b]*wn, N)
        Einv = np.linalg.inv(E)
        self.psiT = np.einsum("pab,pb->pa", Einv[i], d)*wn[:, None]
        self.rev = np.concatenate([np.arange(self.np_, 2*self.np_), np.arange(self.np_)])  # index of (j,i)
        self._wls = {}

    def wls(self, order):
        """Per-pair weights G (npairs, P) with coefficient_c(i) = sum_j G[p, c] (f_j - f_i); coefficients
        ordered as derivatives: fx, fy, fxx, fxy, fyy, fxxx, fxxy, fxyy, fyyy (Taylor basis with factorials)."""
        if order in self._wls: return self._wls[order]
        i, d, N = self.i, self.d, self.N
        s = 1.0/self.h[i]
        X, Y = d[:, 0]*s, d[:, 1]*s
        cols = [X, Y]
        if order >= 2: cols += [X*X/2, X*Y, Y*Y/2]
        if order >= 3: cols += [X**3/6, X*X*Y/2, X*Y*Y/2, Y**3/6]
        if order >= 4: cols += [X**4/24, X**3*Y/6, X*X*Y*Y/4, X*Y**3/6, Y**4/24]
        b = np.column_stack(cols); P = b.shape[1]
        wn = self.w/self.omega[i]
        M = np.zeros((N, P, P))
        for a in range(P):
            for c in range(P):
                M[:, a, c] = np.bincount(i, b[:, a]*b[:, c]*wn, N)
        Minv = np.linalg.inv(M)
        G = np.einsum("pab,pb->pa", Minv[i], b)*wn[:, None]
        # undo the 1/h scaling: coefficient of a degree-q monomial scales with h^-q
        deg = np.array([1, 1] + [2]*3*(order >= 2) + [3]*4*(order >= 3) + [4]*5*(order >= 4))[:P]
        G = G*(s[:, None]**deg[None, :])
        self._wls[order] = G
        return G

    def coeffs(self, f, order):
        G = self.wls(order)
        return np.stack([np.bincount(self.i, G[:, c]*(f[self.j] - f[self.i]), self.N) for c in range(G.shape[1])], 1)


def face_value(g, f, C, side_d, quad):
    """Reconstruct f at x_i + side_d from coefficients C (per particle, evaluated at pair source g).
    quad: factor on the quadratic Taylor term (1 = Taylor, 2/3 = point-value/deconvolution corrected, 0 = linear)."""
    dxv, dyv = side_d[:, 0], side_d[:, 1]
    val = f[g] + C[g, 0]*dxv + C[g, 1]*dyv
    if C.shape[1] >= 5 and quad != 0:
        val = val + quad*(C[g, 2]*dxv*dxv/2 + C[g, 3]*dxv*dyv + C[g, 4]*dyv*dyv/2)
    return val


def rates(n, recon=1, quad=1.0, grad_order=1, divv_order=None, vol="exact", pbg=False, riemann="acoustic",
          quad3=0.0, nnn=32, face_order=1):
    """recon: 1 linear / 2 quadratic reconstruction of the face states (coefficients from WLS of order
    max(recon, grad_order)); quad: factor on the quadratic term; grad_order: WLS order used for the gradient in
    the reconstruction; divv_order: WLS order for the explicit-density divergence (default grad_order);
    vol: 'exact' (V = m/rho_exact, explicit volume seeded exactly) or 'kernel' (V = 1/omega);
    pbg: subtract the (global) background pressure from the momentum flux."""
    pos, dx = lattice(n)
    g = Geometry(pos, dx, nnn)
    N = g.N; x = pos[:, 0]
    rho_bar = 1.0; P_bar = 0.6; c = 1.0
    rho = rho_bar*(1 + A_AMP*np.sin(K*x))
    u = P_bar/((GAMMA-1)*rho_bar)*(rho/rho_bar)**(GAMMA-1)
    P = (GAMMA-1)*rho*u
    vx = c*A_AMP*np.sin(K*x); vy = np.zeros(N)
    m = rho_bar*LX*LY/N
    rhoK = m*g.omega
    if vol == "exact":
        V = m/rho
    elif vol == "kernel":
        V = 1.0/g.omega
    else:
        raise ValueError(vol)
    rhoV = m/V                       # density the solver sees
    PV = (GAMMA-1)*rhoV*u
    i, j, d = g.i, g.j, g.d
    psiF = g.psiT if face_order == 1 else g.wls(face_order)[:, :2]
    Aij = V[i, None]*psiF - V[j, None]*psiF[g.rev]
    Aabs = np.hypot(Aij[:, 0], Aij[:, 1]); nh = Aij/np.where(Aabs > 0, Aabs, 1.0)[:, None]
    ordc = max(recon, grad_order)
    Cr = g.coeffs(rhoV, ordc); Cp = g.coeffs(PV, ordc); Cu = g.coeffs(vx, ordc); Cv = g.coeffs(vy, ordc)
    q = quad if recon >= 2 else 0.0
    half = d/2
    rL = face_value(i, rhoV, Cr, half, q); rR = face_value(j, rhoV, Cr, -half, q)
    pL = face_value(i, PV, Cp, half, q);   pR = face_value(j, PV, Cp, -half, q)
    uL = face_value(i, vx, Cu, half, q);   uR = face_value(j, vx, Cu, -half, q)
    wL = face_value(i, vy, Cv, half, q);   wR = face_value(j, vy, Cv, -half, q)
    vnL = uL*nh[:, 0] + wL*nh[:, 1]; vnR = uR*nh[:, 0] + wR*nh[:, 1]
    rbar = .5*(rL + rR); cs = c
    Z = rbar*cs
    if riemann == "acoustic":
        Ps = .5*(pL + pR) + .5*Z*(vnL - vnR)
        Vs = .5*(vnL + vnR) + .5*(pL - pR)/Z
    elif riemann == "central":
        Ps = .5*(pL + pR); Vs = .5*(vnL + vnR)
    else:
        raise ValueError(riemann)
    Pm = Ps - P_bar if pbg else Ps
    mom = -np.stack([np.bincount(i, Aij[:, a]*Pm, N) for a in range(2)], 1)   # m dv/dt
    en = -np.bincount(i, Aabs*Ps*Vs, N)                                         # dE/dt (v frame included in Vs)
    dvdt = mom/m
    dudt = en/m - (vx*dvdt[:, 0] + vy*dvdt[:, 1])
    do = divv_order or grad_order
    Cu2 = g.coeffs(vx, do); Cv2 = g.coeffs(vy, do)
    dlnrho = -(Cu2[:, 0] + Cv2[:, 1])
    ex_a = -c*c*A_AMP*K*np.cos(K*x)
    ex_u = -(P/rho)*c*A_AMP*K*np.cos(K*x)
    ex_r = -c*A_AMP*K*np.cos(K*x)
    sc = A_AMP*K
    out = dict(ev=np.mean(np.abs(dvdt[:, 0] - ex_a))/sc, evy=np.mean(np.abs(dvdt[:, 1]))/sc,
               eu=np.mean(np.abs(dudt - ex_u))/(sc*P_bar/rho_bar), er=np.mean(np.abs(dlnrho - ex_r))/sc,
               ek=np.mean(np.abs(rhoK/np.mean(rhoK) - rho/np.mean(rho)))/A_AMP)
    # geometric channel alone: -P_bar sum_j A_ij / m, vs the V grad ln(rho_true/rho_E) prediction (= 0 for exact)
    geo = -P_bar*np.stack([np.bincount(i, Aij[:, a], N) for a in range(2)], 1)/m
    out["geo"] = np.mean(np.abs(geo[:, 0]))/sc
    return out


if __name__ == "__main__":
    import sys
    ns = [32, 64, 128, 256]
    cases = [
        ("lin GO3 face3", dict(recon=1, grad_order=3, face_order=3)),
        ("lin GO3 face2", dict(recon=1, grad_order=3, face_order=2)),
        ("quad2/3 GO3 face3", dict(recon=2, grad_order=3, quad=2/3, face_order=3)),
        ("lin GO3 face3 kernelV", dict(recon=1, grad_order=3, face_order=3, vol="kernel")),
    ]
    cases_old = [
        ("lin GO1 exactV", dict(recon=1, grad_order=1)),
        ("lin GO3 exactV", dict(recon=1, grad_order=3)),
        ("quad Taylor GO3", dict(recon=2, grad_order=3, quad=1.0)),
        ("quad 2/3 GO3", dict(recon=2, grad_order=3, quad=2/3)),
        ("quad 2/3 GO3 pbg", dict(recon=2, grad_order=3, quad=2/3, pbg=True)),
        ("quad 2/3 GO3 kernelV", dict(recon=2, grad_order=3, quad=2/3, vol="kernel")),
    ]
    for name, kw in cases:
        prev = None
        for n in ns:
            r = rates(n, **kw)
            s = "%-22s n=%4d  ev %.3e  eu %.3e  er %.3e  geo %.3e  evy %.1e  ek %.2e" % (
                name, n, r["ev"], r["eu"], r["er"], r["geo"], r["evy"], r["ek"])
            if prev:
                s += "   p: ev %.2f eu %.2f er %.2f geo %.2f" % tuple(np.log2(prev[k]/r[k]) for k in ("ev", "eu", "er", "geo"))
            print(s, flush=True); prev = r
        print()
