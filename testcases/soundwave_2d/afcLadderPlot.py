#!/usr/bin/env python3
"""Ideal gas: does the adiabatic face correction restore second order in dt?
Data: Mars fixed-dt sweep, 2026-09-23/24."""
import numpy as np, matplotlib as mpl, matplotlib.pyplot as plt

SURF, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
OFF, ON = "#2a78d6", "#eb6834"

DT64 = np.array([5.208333e-03,3.906250e-03,2.604167e-03,1.953125e-03,1.302083e-03,
                 9.765625e-04,6.510417e-04,4.882812e-04,3.255208e-04,2.441406e-04,
                 1.627604e-04,1.220703e-04,8.138021e-05,6.103516e-05])
DT128 = np.array([2.604167e-03,1.302083e-03,6.510417e-04,3.255208e-04,1.627604e-04,
                  8.138021e-05,4.069010e-05])
D = {
 "ideal":        (DT64,  np.array([3.677817e-02,2.954357e-02,2.233587e-02,1.875613e-02,1.497209e-02,1.309291e-02,1.117858e-02,1.019465e-02,9.261952e-03,8.792764e-03,8.316562e-03,8.076339e-03,7.836371e-03,7.719676e-03])),
 "idealAFC":     (DT64,  np.array([1.692194e-02,1.462170e-02,1.157805e-02,1.012349e-02,8.765403e-03,8.094344e-03,7.439658e-03,7.123119e-03,6.787064e-03,6.620957e-03,6.450541e-03,6.365031e-03,6.277670e-03,6.235233e-03])),
 "idealLimN":    (DT64,  np.array([3.549090e-02,2.793252e-02,2.041606e-02,1.667705e-02,1.294780e-02,1.108573e-02,9.227710e-03,8.299866e-03,7.373151e-03,6.909654e-03,6.444107e-03,6.211997e-03,5.980978e-03,5.865676e-03])),
 "idealLimNAFC": (DT64,  np.array([1.570959e-02,1.280828e-02,9.976473e-03,8.586611e-03,7.214116e-03,6.534290e-03,5.858790e-03,5.522618e-03,5.187526e-03,5.020407e-03,4.853536e-03,4.770193e-03,4.686926e-03,4.645329e-03])),
 "idealN128":    (DT128, np.array([1.675347e-02,9.168044e-03,5.335547e-03,3.373255e-03,2.414629e-03,1.933632e-03,1.692484e-03])),
 "idealAFCN128": (DT128, np.array([6.664464e-03,3.903415e-03,2.541725e-03,1.875866e-03,1.543016e-03,1.378400e-03,1.293930e-03])),
}
def fit_p(E, dt):
    best = None
    for p in np.arange(.2, 4.001, .001):
        X = np.vstack([dt**p, np.ones_like(dt)]).T
        c, *_ = np.linalg.lstsq(X, E, rcond=None)
        ss = float(np.sum((E - X@c)**2))
        if best is None or ss < best[0]: best = (ss, p, c[0], c[1])
    return best[1], best[2], best[3]
def fit_lin(E, dt):
    X = np.vstack([dt, dt**2, np.ones_like(dt)]).T
    c, *_ = np.linalg.lstsq(X, E, rcond=None); return c

mpl.rcParams.update({"figure.facecolor":SURF,"axes.facecolor":SURF,"savefig.facecolor":SURF,
 "text.color":INK,"axes.labelcolor":INK2,"axes.edgecolor":MUTED,"xtick.color":INK2,
 "ytick.color":INK2,"font.size":9,"axes.linewidth":.8,"axes.grid":True,
 "grid.color":"#e4e3de","grid.linewidth":.7,"legend.frameon":False})
fig,((axA,axB),(axC,axD)) = plt.subplots(2,2,figsize=(12.6,9.4))
F = {k: fit_p(*D[k][::-1]) for k in D}

def panel(ax, keys, title, xlab=True):
    for k, col, ls, lab in keys:
        dt, E = D[k]; p,a,b = F[k]
        g = np.logspace(np.log10(dt.min())-.3, np.log10(dt.max())+.25, 200)
        ax.plot(g, a*g**p + b, ls, color=col, lw=2, zorder=2)
        ax.axhline(b, color=col, lw=1, ls=(0,(1,3)), alpha=.8, zorder=1)
        ax.plot(dt, E, "o", color=col, ms=7.5, mec=SURF, mew=1.3, zorder=3,
                label="%s   $p=%.2f$" % (lab, p))
    ax.set(xscale="log", yscale="log", ylabel="$1-a_\\rho$")
    if xlab: ax.set_xlabel("Zeitschritt  $\\Delta t$")
    ax.set_title(title, color=INK, fontsize=11, loc="left", pad=9)
    ax.legend(loc="lower right", fontsize=8.5, labelcolor=INK2)

panel(axA, [("ideal",OFF,"-","Korrektur aus"),("idealAFC",ON,"-","Korrektur AN")],
      "(a)  ideales Gas, $n=64$ — Dämpfung")
panel(axB, [("idealN128",OFF,"-","Korrektur aus"),("idealAFCN128",ON,"-","Korrektur AN")],
      "(b)  ideales Gas, $n=128$")

# (c) floor removed -> slope is p; the two curves are parallel
gl = np.logspace(np.log10(4e-5), np.log10(3e-3), 40)
for sl,st in ((1,(0,(6,3))),(2,(0,(2,2)))):
    yg = 3.2e-2*(gl/3e-3)**sl
    axC.plot(gl, yg, color=MUTED, lw=1.1, ls=st, zorder=1)
    axC.annotate("Steigung %d"%sl,(gl[0],yg[0]),xytext=(-3,0),textcoords="offset points",
                 fontsize=8.5,color=MUTED,style="italic",ha="right",va="center")
for k,col,ls,lab in (("ideal",OFF,"-","aus, m. Limiter"),("idealAFC",ON,"-","AN, m. Limiter"),
                     ("idealLimN",OFF,"--","aus, o. Limiter"),("idealLimNAFC",ON,"--","AN, o. Limiter")):
    dt,E = D[k]; p,a,b = F[k]
    axC.plot(dt, E-b, ls, marker="o", color=col, lw=2, ms=6.5, mec=SURF, mew=1.2,
             zorder=3, label="%s  $p=%.2f$"%(lab,p))
axC.set(xscale="log", yscale="log", xlabel="Zeitschritt  $\\Delta t$",
        ylabel="$1-a_\\rho-b$   (Boden abgezogen)")
axC.set_title("(c)  Boden abgezogen — beide bleiben Steigung 1", color=INK, fontsize=11, loc="left", pad=9)
axC.legend(loc="lower right", fontsize=8, labelcolor=INK2)

# (d) the linear coefficient from E = a dt + c dt^2 + b
rows = [("$n{=}64$, m. Limiter","ideal","idealAFC"),("$n{=}64$, o. Limiter","idealLimN","idealLimNAFC"),
        ("$n{=}128$","idealN128","idealAFCN128")]
y = np.arange(len(rows))[::-1]; h=.33
for i,(lab,ko,kn) in zip(y,rows):
    for off,k,col,tag in ((h/2,ko,OFF,"aus"),(-h/2,kn,ON,"AN")):
        a = fit_lin(*D[k][::-1])[0]
        axD.barh(i+off, a, height=h-.05, color=col, zorder=2)
        axD.annotate("%.2f  (%s)"%(a,tag),(a,i+off),xytext=(5,0),textcoords="offset points",
                     va="center",fontsize=8.5,color=INK2)
axD.axvline(0, color=MUTED, lw=1.2)
axD.set_yticks(y,[r[0] for r in rows],fontsize=9.5)
axD.set(xlabel="linearer Koeffizient $a$ in  $E=a\\,\\Delta t+c\\,\\Delta t^2+b$", xlim=(0,9.0))
axD.set_title("(d)  Der $O(\\Delta t)$-Term verschwindet nicht", color=INK, fontsize=11, loc="left", pad=9)
axD.grid(axis="y", visible=False)
axD.annotate("wäre die Ordnung auf 2 gehoben,\nmüsste $a\\to0$ gehen", (.40,.46),
             xycoords="axes fraction", fontsize=9.5, color=INK2, style="italic")

fig.suptitle("Adiabatische Flächenkorrektur: senkt die Dämpfung um Faktor 2.8, "
             "ändert die Ordnung nicht", fontsize=13.5, color=INK, x=.012, ha="left", y=.985)
fig.text(.012,.947,"Fixed-dt-Leiter (ADAPTIVE_TIMESTEP 0), ideales Gas, $t=1$, 14 bzw. 7 Sprossen "
  "über Faktor 85 in $\\Delta t$ — überbestimmter Fit statt der alten drei Punkte.\n"
  "Der räumliche Boden $b$ (gepunktet) fällt von $n=64$ auf $128$ um Faktor 5.1 (aus) bzw. 5.0 (AN), "
  "ist also $\\gtrsim O(h^2)$: der Defekt sitzt allein in der Zeit.",
  fontsize=8.6,color=INK2,ha="left",va="top")
fig.tight_layout(rect=(0,0,1,.917))
out="/private/tmp/claude-503/-Users-jsappler-Documents-MA-all/aae1b576-4ff2-48dc-8834-792701dcba0f/scratchpad/afc_result.png"
fig.savefig(out,dpi=170); print("wrote",out)
for k in D: print("%-14s p=%.3f  a_lin=%.3f" % (k, F[k][0], fit_lin(*D[k][::-1])[0]))
