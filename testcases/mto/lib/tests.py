"""Test registry of the mto testplan (vault: 2026-10-02-mto-testplan.md).

A *test* is a physical testcase with a resolution ladder of *rungs*; every rung
is one run. A test names the build *variant* it needs (compile-time flags that
belong to the testcase -- periodicity, EOS, plasticity -- not to the candidate).

Costs: the estimate is N * steps * C_PS with dt = CFL * h / (2 c), h = d sqrt(32/pi),
i.e. dt ~ 0.32 d / c at CFL 0.2. C_PS is calibrated per machine in calib.json
(`mto.py calib`); until then the Mars value 3.5e-5 s is assumed.
"""

import math

SQ = math.sqrt(32.0 / math.pi)          # h / d on a lattice at N_NN = 32
MU_BODY = 0.22                          # rings / blocks / disc material (09-26 plan)
MU_WAVE = 6.0 / 13.0                    # nu = 0.3 at K = 1 (wave studies)

# rod material, Taylor impact: rho c_l v = 0.26 > Y0 = 0.1 -> classic mushroom, rebound at t~3
# (GIZMO Y0 scan 2026-10-02: 0.05 stays welded, 0.2 barely yields; gizmo/y0_scan.png)
ROD = dict(mu=0.5, Y0=0.1, v=0.2, length=1.0, width=0.25, gap=0.1)


def till(mu, Y0=None):
    m = dict(matId=0, rho0=1.0, A=1.0, B=0.0, u0=1.0, a=0.0, b=0.0, alpha=0.0, beta=0.0,
             u_iv=0.05, u_cv=0.18, mu=mu)
    if Y0 is not None:
        m["Y0"] = Y0
    return m


MURN = dict(matId=0, K0=1.0, n=2.0, rho0=1.0, mu=0.0)

# variant -> flags forced on top of base + candidate (the testcase owns these)
VARIANTS = {
    # SURFACE_VOLCORR is a scheme choice (base.flags = 1), not owned by the body tests;
    # the periodic wave boxes have no free surface, where it is a no-op (09-25).
    "body":    dict(PERIODIC_BOUNDARIES=0, EOS=2),
    "plastic": dict(PERIODIC_BOUNDARIES=0, EOS=2,
                    VON_MISES_PLASTICITY=1, PLASTICITY_YIELD_GIZMO=1),
    "wave_e0": dict(PERIODIC_BOUNDARIES=1, EOS=0, SURFACE_VOLCORR=0),
    "wave_e1": dict(PERIODIC_BOUNDARIES=1, EOS=1, SURFACE_VOLCORR=0),
    "wave_e2": dict(PERIODIC_BOUNDARIES=1, EOS=2, SURFACE_VOLCORR=0),
}

EOS_VARIANT = {"ideal": "wave_e0", "murnaghan": "wave_e1", "tillotson": "wave_e2"}


def _wave_speed(eos, mu, mode):
    if mode == "s":
        return math.sqrt(mu)
    K = 5.0 / 3.0 * 0.6 if eos == "ideal" else 1.0   # gamma P_bar for the ideal gas
    return math.sqrt(K + 4.0 / 3.0 * mu)


def _body(test, gen, d, t_end, dump, box, c, mat, **genargs):
    n_est = None   # filled by mto.py from the IC
    return dict(test=test, gen=gen, d=d, genargs=genargs, t_end=t_end, dump=dump,
                box=box, kernel=d * SQ, c=c, materials=[mat], N=n_est)


def rungs_T1():   # colliding rings, contact at t ~ 7.5 (shift 5), rebound done by 180
    c = math.sqrt(1 + 4 / 3 * MU_BODY)
    return {"d{:g}".format(d): _body("T1", "rings", d, 180.0, 1.0, 20.0, c, till(MU_BODY))
            for d in (0.2, 0.1, 0.08, 0.06)}


def rungs_T2():   # offset blocks, contact t ~ 0.85
    c = math.sqrt(1 + 4 / 3 * MU_BODY)
    return {"d{:g}".format(d): _body("T2", "blocks", d, 20.0, 0.2, 3.0, c, till(MU_BODY))
            for d in (0.05, 0.025, 0.0177, 0.0125)}


def rungs_A1():   # rotating disc, one rotation at omega = 0.2
    c = math.sqrt(1 + 4 / 3 * MU_BODY)
    T = 2 * math.pi / 0.2
    return {"d{:g}".format(d): _body("A1", "disc", d, T, T / 20, 3.0, c, till(MU_BODY),
                                     omega=0.2)
            for d in (0.05, 0.035, 0.025, 0.0177)}


def rungs_P1():   # Taylor rods
    c = math.sqrt(1 + 4 / 3 * ROD["mu"])
    w = ROD["width"]
    out = {}
    for cells in (10, 16, 24):
        d = w / cells
        out["w{}".format(cells)] = _body("P1", "rod", d, 6.0, 0.1, 4.0, c,
                                         till(ROD["mu"], ROD["Y0"]),
                                         length=ROD["length"], width=w, gap=ROD["gap"],
                                         v=ROD["v"])
    return out


def rungs_wave(test, mode, eos, mu, ns=(32, 64, 128, 256)):
    c = _wave_speed(eos, mu, mode)
    T = 1.0 / c
    out = {}
    for n in ns:
        d = 1.0 / n
        mat = None if eos == "ideal" else (dict(MURN, mu=mu) if eos == "murnaghan"
                                           else till(mu))
        out["n{}".format(n)] = dict(test=test, gen="wave", d=d,
                                    genargs=dict(n=n, mode=mode, eos=eos, mu=mu),
                                    t_end=T, dump=T / 10, box=None, kernel=d * SQ,
                                    c=_wave_speed(eos, mu, "p"),
                                    materials=[mat] if mat else [], N=n * n // 4,
                                    eos=eos)
    return out


TESTS = {
    # name: (variant, rungs, what)
    "T1": ("body", rungs_T1, "colliding rings: no rip, no clumps"),
    "T2": ("body", rungs_T2, "offset blocks: no rip, no interior u/rho lines"),
    "P1": ("plastic", rungs_P1, "Taylor rods: contact-plane clumping, fold, flow vs GIZMO"),
    "A1": ("body", rungs_A1, "rotating disc: |dLz/Lz| after one rotation ~ d^2"),
    "W1i": ("wave_e0", lambda: rungs_wave("W1i", "p", "ideal", 0.0), "sound, ideal gas"),
    "W1m": ("wave_e1", lambda: rungs_wave("W1m", "p", "murnaghan", 0.0), "sound, Murnaghan"),
    "W1t": ("wave_e2", lambda: rungs_wave("W1t", "p", "tillotson", 0.0), "sound, Tillotson"),
    "W2": ("wave_e2", lambda: rungs_wave("W2", "p", "tillotson", MU_WAVE), "P-wave, mu = 6/13"),
    "W3": ("wave_e2", lambda: rungs_wave("W3", "s", "tillotson", MU_WAVE), "shear wave"),
}

# Tiers = which (test, rung) pairs a stage runs. Order inside a tier does not matter:
# the queue sorts by estimated cost, longest first.
TIERS = {
    # S: smoke, ~5 min wall at -j 4. Catches crashes, gross rips, broken waves.
    "S": [("T1", "d0.2"), ("T2", "d0.025"), ("A1", "d0.05"), ("P1", "w10"),
          ("W1t", "n64"), ("W2", "n64"), ("W3", "n64")],
    # T: tensile gate. Blocks ladder (fast stand-in) + the three ring rungs.
    "T": [("T2", "d0.05"), ("T2", "d0.025"), ("T2", "d0.0177"), ("T2", "d0.0125"),
          ("T1", "d0.1"), ("T1", "d0.08"), ("T1", "d0.06")],
    # P: plastic gate (compare with the GIZMO control in gizmo/).
    "P": [("P1", "w10"), ("P1", "w16"), ("P1", "w24")],
    # A: angular momentum order.
    "A": [("A1", "d0.05"), ("A1", "d0.035"), ("A1", "d0.025"), ("A1", "d0.0177")],
    # W: wave order; n = 32 is diagnostic only, the order is judged on 64/128/256.
    "W": [(t, "n{}".format(n)) for t in ("W1i", "W1m", "W1t", "W2", "W3")
          for n in (64, 128, 256)],
    "Wq": [(t, "n{}".format(n)) for t in ("W1t", "W2", "W3") for n in (32, 64, 128)],
}


def rung(test, tag):
    variant, fn, _ = TESTS[test]
    r = dict(fn()[tag])
    r["variant"] = variant
    r["tag"] = tag
    return r


def all_rungs(test):
    return list(TESTS[test][1]().keys())
