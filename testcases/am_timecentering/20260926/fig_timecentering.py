#!/usr/bin/env python3
"""Figure 1: time-centring evidence on the rotating disc (d0.05, omega 0.2, one rotation).
(a) dL/L(t): standard (DRIFT_KDK 1) vs endpoint drift vs TIME_DKD, with the kick/drift prediction +omega^2 dt/2 t.
(b) rho_E/rho_K - 1 (bulk r<0.7): same runs, prediction +omega^2 dt t.
(c) integrated torque of the transverse RECON_FACE_DV channel: predicted jump vs spatial jump, CFL 0.2 / 0.1."""
import os, glob, re, numpy as np, h5py, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from analyze import dts
from tqparse import parse
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]; INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
W = 0.2
def series(run):
    T, L, R = [], [], []
    for f in sorted(glob.glob(f"runs/{run}/output/*.h5")):
        with h5py.File(f) as h:
            x = h["x"][()]; v = h["v"][()]; m = h["m"][()]; T.append(h["time"][()][0])
            L.append(np.sum(m*(x[:, 0]*v[:, 1]-x[:, 1]*v[:, 0]))); b = np.linalg.norm(x, axis=1) < 0.7
            rk = h["rhoKernel"][()]; R.append(np.mean(h["rho"][()][b]/rk[b]) - 1 if rk[b].min() > 0 else np.nan)
    T = np.array(T); L = np.array(L); return T, (L-L[0])/L[0], np.array(R)
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(1, 3, figsize=(13, 3.9), facecolor=SURF)
for a in ax: a.set_facecolor(SURF); a.grid(True, color=GRID, lw=0.6); a.set_xlabel("t")
runs = [("A_vfix_d0.05", "Standard (DRIFT_KDK 1)"), ("A_ep_d0.05", "Endpunkt-Drift"), ("A_dkd_d0.05", "TIME_DKD")]
dtm = dts("A_vfix_d0.05")[0].mean()
for k, (r, lab) in enumerate(runs):
    T, rel, R = series(r)
    ax[0].plot(T, rel, color=C[k], lw=2, label=lab); ax[1].plot(T[2:], R[2:], color=C[k], lw=2, label=lab)
tt = np.linspace(0, 31.4, 50)
ax[0].plot(tt, W*W*dtm/2*tt, "--", color=INK2, lw=1.2, label=r"Vorhersage $\omega^2\,dt\,t/2$")
ax[1].plot(tt, W*W*dtm*tt, "--", color=INK2, lw=1.2, label=r"Vorhersage $\omega^2\,dt\,t$")
ax[0].set_ylabel(r"$\Delta L_z/L_z$"); ax[0].set_title("(a) Drehimpuls, d0.05, CFL 0.2", loc="left", color=INK)
ax[1].set_ylabel(r"$\rho_E/\rho_K-1$ (Kern $r<0.7$)"); ax[1].set_title("(b) explizite gegen Kerndichte", loc="left", color=INK)
for k, (r, lab) in enumerate([("B_dkd_tq_d0.05", "Sprung mit Prädiktor, CFL 0.2"), ("B_dkd_tq_c01_d0.05", "Sprung mit Prädiktor, CFL 0.1"),
                              ("C_dkd_sp_tq_d0.05", "räumlicher Sprung, CFL 0.2"), ("C_dkd_sp_tq_c01_d0.05", "räumlicher Sprung, CFL 0.1")]):
    t, I = parse(r); L0 = 0.310390
    ax[2].plot(t, I[:, 2]/L0, color=C[k], lw=2, label=lab)
ax[2].set_ylabel(r"$\int$ Drehmoment Querkanal $/L_z$"); ax[2].set_title("(c) RECON_FACE_DV-Querkanal (TIME_DKD)", loc="left", color=INK)
for a in ax: a.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig("figs/timecentering.png", dpi=130, facecolor=SURF); print("figs/timecentering.png")
