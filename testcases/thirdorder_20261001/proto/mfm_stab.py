#!/usr/bin/env python3
"""Linear stability of the demonstrator's periodic MFM right-hand side (2D, ideal gas, mu = 0).

State U = (x, y, vx, vy, rhoE, u) per particle; RHS as in mfm_trunc.py but for an arbitrary state:
geometry (h from N_NN, omega, psiTilde, WLS weights, faces A_ij) rebuilt from the positions,
V = m/rhoE, face states reconstructed linearly from WLS gradients, linear acoustic Riemann solver,
d ln rhoE/dt = -div v (WLS). The Jacobian at the uniform lattice at rest (rhoE = 1, P = P0) is
built by central differences; its eigenvalues say whether the semi-discrete scheme is stable.
Usage: ./mfm_stab.py n [variant ...]"""
import sys
import numpy as np
from scipy.spatial import cKDTree
from mfm_trunc import W, LX, LY, GAMMA

P0 = 0.6


class Scheme:
    def __init__(self, n, face_order=1, grad_order=1, swap=False, P0=P0, nnn=32, divv_order=None,
                 riemann="acoustic", quad=0.0, ny=12, pbg=False, divv_face=False, beta=0.0, hp_passes=1,
                 qmode="ratio"):
        self.n = n; self.dx = LX/n; self.LY = ny*self.dx
        q = (np.arange(n) + .5)*self.dx; ys = (np.arange(ny) + .5)*self.dx
        X, Y = np.meshgrid(q, ys, indexing="ij")
        self.X0 = np.column_stack([X.ravel(), Y.ravel()]); self.N = len(self.X0)
        self.m = LX*self.LY/self.N
        self.face_order, self.grad_order, self.swap, self.P0, self.nnn = face_order, grad_order, swap, P0, nnn
        self.divv_order = divv_order or grad_order
        self.riemann = riemann; self.quad = quad; self.pbg = pbg; self.divv_face = divv_face
        self.beta = beta; self.hp_passes = hp_passes; self.qmode = qmode
        h0 = np.sqrt(nnn/np.pi)*self.dx
        tree = cKDTree(self.X0, boxsize=[LX, self.LY])
        pr = tree.query_pairs(1.35*h0, output_type="ndarray")
        self.i = np.concatenate([pr[:, 0], pr[:, 1]]); self.j = np.concatenate([pr[:, 1], pr[:, 0]])
        self.np_ = len(pr)
        self.rev = np.concatenate([np.arange(self.np_, 2*self.np_), np.arange(self.np_)])
        self.h0 = h0
        u0 = P0/((GAMMA - 1)*1.0) if P0 > 0 else 0.9
        self.U0 = np.concatenate([self.X0[:, 0], self.X0[:, 1], np.zeros(self.N), np.zeros(self.N),
                                  np.ones(self.N), np.full(self.N, u0)])

    def wls(self, d, w, h, order):
        i, N = self.i, self.N
        s = 1.0/h[i]; X, Y = d[:, 0]*s, d[:, 1]*s
        cols = [X, Y]
        if order >= 2: cols += [X*X/2, X*Y, Y*Y/2]
        if order >= 3: cols += [X**3/6, X*X*Y/2, X*Y*Y/2, Y**3/6]
        b = np.column_stack(cols); P = b.shape[1]
        M = np.zeros((N, P, P))
        for a in range(P):
            for c in range(P):
                M[:, a, c] = np.bincount(i, b[:, a]*b[:, c]*w, N)
        G = np.einsum("pab,pb->pa", np.linalg.inv(M)[i], b)*w[:, None]
        deg = np.array([1, 1] + [2]*3 + [3]*4)[:P]
        return G*(s[:, None]**deg[None, :])

    def rhs(self, U):
        N, i, j, m = self.N, self.i, self.j, self.m
        x, y, vx, vy, rhoE, u = (U[k*N:(k+1)*N] for k in range(6))
        d = np.column_stack([x[j] - x[i], y[j] - y[i]])
        d[:, 0] -= LX*np.round(d[:, 0]/LX); d[:, 1] -= self.LY*np.round(d[:, 1]/self.LY)
        r = np.hypot(d[:, 0], d[:, 1])
        lo = np.full(N, 0.6*self.h0); hi = np.full(N, 1.3*self.h0)
        for _ in range(60):
            h = .5*(lo + hi)
            om = W(0., h) + np.bincount(i, W(r, h[i]), N)
            f = np.pi*h*h*om - self.nnn
            lo = np.where(f < 0, h, lo); hi = np.where(f >= 0, h, hi)
        h = .5*(lo + hi)
        w = W(r, h[i])
        om_final = W(0., h) + np.bincount(i, w, N)
        E = np.zeros((N, 2, 2))
        for a in range(2):
            for b in range(2):
                E[:, a, b] = np.bincount(i, d[:, a]*d[:, b]*w, N)
        psi1 = np.einsum("pab,pb->pa", np.linalg.inv(E)[i], d)*w[:, None]
        Gg = self.wls(d, w, h, max(self.grad_order, 2 if self.quad else 1)) if (self.grad_order > 1 or self.quad) else None
        psiF = psi1 if self.face_order == 1 else self.wls(d, w, h, self.face_order)[:, :2]
        V = m/rhoE
        Aij = V[i, None]*psiF - V[j, None]*psiF[self.rev]
        Aabs = np.hypot(Aij[:, 0], Aij[:, 1]); nh = Aij/np.where(Aabs > 0, Aabs, 1.)[:, None]
        P = (GAMMA - 1)*rhoE*u
        def coeffs(fv, order):
            G = psi1 if order == 1 else self.wls(d, w, h, order)
            return np.stack([np.bincount(i, G[:, c]*(fv[j] - fv[i]), N) for c in range(G.shape[1])], 1)
        go = max(self.grad_order, 2) if self.quad else self.grad_order
        Cp, Cu, Cv = coeffs(P, go), coeffs(vx, go), coeffs(vy, go)
        def face(f, C, g, dd):
            val = f[g] + C[g, 0]*dd[:, 0] + C[g, 1]*dd[:, 1]
            if self.quad and C.shape[1] >= 5:
                val = val + self.quad*(C[g, 2]*dd[:, 0]**2/2 + C[g, 3]*dd[:, 0]*dd[:, 1] + C[g, 4]*dd[:, 1]**2/2)
            return val
        pL = face(P, Cp, i, d/2); pR = face(P, Cp, j, -d/2)
        uL = face(vx, Cu, i, d/2); uR = face(vx, Cu, j, -d/2)
        wL = face(vy, Cv, i, d/2); wR = face(vy, Cv, j, -d/2)
        vnL = uL*nh[:, 0] + wL*nh[:, 1]; vnR = uR*nh[:, 0] + wR*nh[:, 1]
        if self.swap:
            back = np.sum(Aij*d, 1) < 0
            pL, pR = np.where(back, pR, pL), np.where(back, pL, pR)
            vnL, vnR = np.where(back, vnR, vnL), np.where(back, vnL, vnR)
        rho_f = .5*(rhoE[i] + rhoE[j]); cs = np.sqrt(GAMMA*P0/1.0) if self.P0 > 0 else 1.0
        Z = rho_f*cs
        if self.riemann == "acoustic":
            Ps = .5*(pL + pR) + .5*Z*(vnL - vnR); Vs = .5*(vnL + vnR) + .5*(pL - pR)/Z
        else:
            Ps = .5*(pL + pR); Vs = .5*(vnL + vnR)
        Pm = Ps - self.P0 if self.pbg else Ps
        mom = -np.stack([np.bincount(i, Aij[:, a]*Pm, N) for a in range(2)], 1)
        self._A1 = V[i, None]*psi1 - V[j, None]*psi1[self.rev]
        if self.beta:
            # high-pass restoring force: q = ln(rho_K/rho_E), q_hp = q - S^p[q] (Shepard kernel average),
            # pair "pressure" beta * Pbar_ij * (q_hp_i + q_hp_j)/2 on the FIRST-ORDER faces
            S = lambda f: (W(0., h)*f + np.bincount(i, w*f[j], N))/om_final
            if self.qmode == "ratio":
                q = np.log(m*om_final/rhoE)
                qs = q.copy()
                for _ in range(self.hp_passes):
                    qs = S(qs)
                qhp = q - qs
            else:                               # positions only, twicing high-pass (I - S)^2
                q = np.log(m*om_final)
                r1 = q - S(q)
                qhp = r1 - S(r1)
            A1 = self._A1
            Pb = .5*(P[i] + P[j])
            fr = self.beta*Pb*.5*(qhp[i] + qhp[j])
            mom -= np.stack([np.bincount(i, A1[:, a]*fr, N) for a in range(2)], 1)
        en = -np.bincount(i, Aabs*Ps*Vs, N)
        ax, ay = mom[:, 0]/m, mom[:, 1]/m
        dudt = en/m - (vx*ax + vy*ay)
        Cd = coeffs(vx, self.divv_order), coeffs(vy, self.divv_order)
        drho = -rhoE*(Cd[0][:, 0] + Cd[1][:, 1])
        if self.divv_face:   # MFM-consistent volume rate: dV/dt = sum_j A_ij . v*_ij
            dV = np.bincount(i, Aabs*Vs, N)
            drho = -rhoE*dV/V
        return np.concatenate([vx, vy, ax, ay, drho, dudt])

    def jacobian(self):
        N = self.N; U0 = self.U0; n6 = 6*N
        scale = np.concatenate([np.full(2*N, self.dx), np.full(2*N, 1.0), np.ones(N), np.ones(N)])
        J = np.zeros((n6, n6))
        eps = 1e-6
        for k in range(n6):
            e = np.zeros(n6); e[k] = eps*scale[k]
            J[:, k] = (self.rhs(U0 + e) - self.rhs(U0 - e))/(2*e[k])
        return J


VARIANTS = {
    "base": dict(face_order=1, grad_order=1),
    "go3": dict(face_order=1, grad_order=3),
    "f3": dict(face_order=3, grad_order=3),
    "f3swap": dict(face_order=3, grad_order=3, swap=True),
    "f3pbg": dict(face_order=3, grad_order=3, pbg=True),
    "f3pbgqT": dict(face_order=3, grad_order=3, pbg=True, quad=1.0),
    "f3pbgq23": dict(face_order=3, grad_order=3, pbg=True, quad=2/3),
    "f1pbgq23": dict(face_order=1, grad_order=3, pbg=True, quad=2/3),
    "f3b1": dict(face_order=3, grad_order=3, beta=1.0),
    "f3t1": dict(face_order=3, grad_order=3, beta=1.0, qmode="pos"),
    "f3t4": dict(face_order=3, grad_order=3, beta=4.0, qmode="pos"),
    "f3t16": dict(face_order=3, grad_order=3, beta=16.0, qmode="pos"),
    "f3tm4": dict(face_order=3, grad_order=3, beta=-4.0, qmode="pos"),
    "f3b2": dict(face_order=3, grad_order=3, beta=2.0),
    "f3b4": dict(face_order=3, grad_order=3, beta=4.0),
    "f3b2p2": dict(face_order=3, grad_order=3, beta=2.0, hp_passes=2),
    "f3b1p0": dict(face_order=3, grad_order=3, beta=1.0, hp_passes=0),
    "basepbg": dict(face_order=1, grad_order=1, pbg=True),
    "f3dvf": dict(face_order=3, grad_order=3, divv_face=True),
    "f3dvfpbg": dict(face_order=3, grad_order=3, divv_face=True, pbg=True),
    "basedvf": dict(face_order=1, grad_order=1, divv_face=True),
    "f3P0": dict(face_order=3, grad_order=3, P0=0.0),
    "f2": dict(face_order=2, grad_order=3),
    "f3central": dict(face_order=3, grad_order=3, riemann="central"),
    "basecentral": dict(face_order=1, grad_order=1, riemann="central"),
}

if __name__ == "__main__":
    n = int(sys.argv[1])
    for name in sys.argv[2:]:
        S = Scheme(n, **VARIANTS[name])
        J = S.jacobian()
        lam = np.linalg.eigvals(J)
        o = np.argsort(-lam.real)
        print("%-12s n=%d N=%d  max Re %.4e  (|lam| max %.3e)  top: %s" % (
            name, n, S.N, lam.real.max(), np.abs(lam).max(),
            " ".join("%.3g%+.3gi" % (lam[k].real, lam[k].imag) for k in o[:6])), flush=True)
        np.save("jac_%s_n%d.npy" % (name, n), J)
