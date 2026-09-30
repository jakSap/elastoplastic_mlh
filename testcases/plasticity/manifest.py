"""Plasticity test plan -- the ONLY file that defines what is tested.

Plan: MA-Obsidian/2026-09-30-testplan-plastisch.md (Appendix P to the elastic test plan).
Edit this file, then start a NEW study with `./pt new` (a study freezes a copy of it).
Never edit a study's frozen copy (STUDY.json / manifest_frozen.py).

Rules that keep it reliable:
  * Code version: COMMIT only. Change it with `./pt pin <hash>`, not by hand-editing
    builds. Builds are extracted from git (`git archive`), never copied from a working tree.
  * A variant is a flag diff against the demonstrator defaults AT ITS COMMIT. Add a new
    variant here with a comment saying why; never change the flags of an existing one --
    add a new name instead (old studies refer to the old meaning).
  * A test is one case x rungs x variants + one checker. Numbers in PASS come from the plan.
"""
import math

# --------------------------------------------------------------------------- code
REPO = "/home/sappler/Documents/newMFM"                          # git repo (shared by all worktrees)
BRANCH = "plasticity-tests"
COMMIT = "30a6b3b"         # code under test (on BRANCH). Change with: ./pt pin <hash>
BASE_COMMIT = "3835f4a"    # am-timecentering state the branch started from (P1a only)
HOLIB = "/home/sappler/Documents/newMFM/demonstrator/holib"      # untracked header libs (HighFive)
GIZMO_SRC = "/home/sappler/Documents/GIZMO/ring_collision/master"
HIRES_ICS = "/home/sappler/Documents/newMFM/testcases/hires_20260929/ics"   # elastic-plan ICs (P1)
STUDY_ROOT = "/home/sappler/Documents/newMFM/testcases"           # studies: plasticity_<date>_<commit7>/

# ----------------------------------------------------------------------- material
# Canonical Tillotson of the elastic plan (P = A (rho/rho0 - 1), K = 1), plan P.1.2.
MU = 6.0 / 13.0
Y0 = 0.01
MATERIAL = dict(rho0=1.0, A=1.0, B=0.0, u0=1.0, a=0.0, b=0.0, alpha=0.0, beta=0.0,
                u_iv=0.05, u_cv=0.18, mu=MU, Y0=Y0)
NNN = 32
CFL_NOTE = "CFL 0.2 (demonstrator default), GIZMO CourantFac 0.2"

# ----------------------------------------------------------------------- variants
# flags: diff against the demonstrator defaults at the variant's commit (plan P.1.1 / P.2).
PREF = {
    "SOLID_2D_MODEL": 1,          # GIZMO true 2D (user decision 2026-09-30)
    "VON_MISES_PLASTICITY": 1,    # constant Y0 from the material block
    "YIELD_RETURN": 1,            # GIZMO return S *= Y^2/(d J2)
    "OUTPUT_PLASTIC_STRAIN": 1,   # /plasticStrain (P-R3)
    "OUTPUT_FAB": 1,              # TC factor per particle (P-R3)
    "EXPLICIT_VOL_SEED_RHO0": 1,  # t = 0 surface density deficit (plan P.1.1)
}
VARIANTS = {
    "pref":   dict(flags=PREF, doc="reference"),
    "pel":    dict(flags={**PREF, "VON_MISES_PLASTICITY": 0, "OUTPUT_PLASTIC_STRAIN": 0},
                   doc="elastic baseline of the true-2D model"),
    "prr":    dict(flags={**PREF, "YIELD_RETURN": 0}, doc="exact radial return"),
    "ptc0":   dict(flags={**PREF, "TENSILE_CORRECTION": 0, "TENSILE_CORRECTION_PARTICLE": 0,
                          "TENSILE_CORRECTION_DAMP_LIMIT": 0, "TENSILE_CORRECTION_TRACTION_GATE": 0,
                          "OUTPUT_FAB": 0},       # = tgate studies' tc0 + gate off (parameter.h #errors)
                   doc="tensile correction off (P-R4)"),
    "pps":    dict(flags={**PREF, "SOLID_2D_MODEL": 0}, doc="old plane-strain rate"),
    "pns":    dict(flags={**PREF, "EXPLICIT_VOL_SEED_RHO0": 0},   # 2026-09-30: is the t=0 rho0 seed behind the P4 negE?
                   doc="pref without the rho0 density seed (kernel-sum start)"),
    "legacy": dict(flags={}, doc="branch defaults (P1a)"),
    "base":   dict(flags={}, commit=BASE_COMMIT, doc="3835f4a itself (P1a)"),
    "gz":     dict(code="gizmo", doc="GIZMO MFM reference (P-R1)"),
}

# Build groups: extra flags a test needs on top of its variant (plan P.2, second table).
BUILD_GROUPS = {
    "std":  {},
    "per":  {"PERIODIC_BOUNDARIES": 1},                                   # P2, P3 (periodic y)
    "wave": {"PERIODIC_BOUNDARIES": 1, "SURFACE_VOLCORR": 0, "EXPLICIT_VOL_SEED_RHO0": 0},  # Tier W
}

# ------------------------------------------------------------------------- tests
# Keys of a test:
#   tier      plan section           ic      generator name (ptlib/gen.py) + args
#   rungs     list of resolution tags (strings; "0.025" = lattice spacing d, "64" = n)
#   t_end, dump_dt                   box     (x0, y0, x1, y1) of periodicBoxLimits
#   group     build group            material  overrides of MATERIAL
#   variants  newMFM variants (+ "gz" for the GIZMO arm)
#   check     checker name (ptlib/check.py) + its args
#   timeout   seconds per run (default 4 h)
KS = math.sqrt(NNN / math.pi)
P1_RUNGS_NOTE = "P1 uses the elastic plan's smoke configs verbatim (tgate2 setup_runs.py)"

def _p1(tier, case, variants, check, material=None, group="std", **extra):
    base = {
        "rings":  dict(ic=("copy", {"src": HIRES_ICS + "/rings_deltap0.2-2D.h5"}), rungs=["0.2"],
                       dt=1.0, t_end=200.0, dump_dt=1.0, box=(-20., -10., 20., 10.), material={"mu": 0.22}),
        "blocks": dict(ic=("copy", {"src": HIRES_ICS + "/two_blocks_offset_deltap0.05-2D.h5"}), rungs=["0.05"],
                       t_end=20.0, dump_dt=0.2, box=(-3., -3., 3., 3.), material={"mu": 0.22}),
        "disc":   dict(ic=("copy", {"src": HIRES_ICS + "/disc_d0.025_w0.2.h5"}), rungs=["0.025"],
                       t_end=100.531, dump_dt=2.0, box=(-3., -3., 3., 3.), material={"mu": 0.22}),
        "wave":   dict(ic=("copy", {"src": HIRES_ICS + "/elastic_N64.h5"}), rungs=["64"],
                       t_end=3.06 / math.sqrt(1 + 4 / 3 * MU), dump_dt=5 / math.sqrt(1 + 4 / 3 * MU) / 50,
                       box=(0., 0., 1., 0.25), material={}, group="wave", ks_scale="n"),
    }[case]
    t = dict(tier=tier, group=group, variants=variants, check=check)
    t.update(base)
    t.update(extra)
    if material is not None:
        t["material"] = {**t["material"], **material}
    t.setdefault("ks_scale", "d")
    return t

TESTS = {
    # ---- P1a: branch inert (legacy == base, bitwise), elastic smoke configs ---------------
    "P1a_rings":  _p1("P1a", "rings",  ["legacy", "base"], ("identical", {"a": "legacy", "b": "base"})),
    "P1a_blocks": _p1("P1a", "blocks", ["legacy", "base"], ("identical", {"a": "legacy", "b": "base"})),
    "P1a_disc":   _p1("P1a", "disc",   ["legacy", "base"], ("identical", {"a": "legacy", "b": "base"})),
    "P1a_wave":   _p1("P1a", "wave",   ["legacy", "base"], ("identical", {"a": "legacy", "b": "base"})),
    # ---- P1b: plasticity inert at Y0 = 1e30 (pref == pel bitwise, gamma_p == 0) ------------
    "P1b_rings":  _p1("P1b", "rings",  ["pref", "pel"], ("identical", {"a": "pref", "b": "pel"}),
                      material={"Y0": 1e30}),
    "P1b_blocks": _p1("P1b", "blocks", ["pref", "pel"], ("identical", {"a": "pref", "b": "pel"}),
                      material={"Y0": 1e30}),
    # ---- P1c: rigid rotation does not yield at the reference Y0 --------------------------
    "P1c_disc":   _p1("P1c", "disc",   ["pref", "pel"], ("identical", {"a": "pref", "b": "pel"})),
    # ---- P1d: elastic baseline of SOLID_2D_MODEL 1 (Tier S smoke + W.2/W.3) --------------
    "P1d_rings":  _p1("P1d", "rings",  ["pel"], ("smoke", {})),
    "P1d_blocks": _p1("P1d", "blocks", ["pel"], ("smoke", {})),
    "P1d_w2": dict(tier="P1d", ic=("wave", {"mode": "p"}), rungs=["64", "128", "256"],
                   t_end=3 * 1.0 / math.sqrt(1 + MU), dump_dt=1.0 / math.sqrt(1 + MU) / 20,
                   box=(0., 0., 1., 0.25), group="wave", material={}, ks_scale="n",
                   variants=["pel", "gz"],
                   check=("wavespeed", {"mode": "p", "c": math.sqrt(1 + MU), "tol": 1e-3, "at": "256"})),
    "P1d_w3": dict(tier="P1d", ic=("wave", {"mode": "s"}), rungs=["64", "128", "256"],
                   t_end=3 * 1.0 / math.sqrt(MU), dump_dt=1.0 / math.sqrt(MU) / 20,
                   box=(0., 0., 1., 0.25), group="wave", material={}, ks_scale="n",
                   variants=["pel", "gz"],
                   check=("wavespeed", {"mode": "s", "c": math.sqrt(MU), "tol": 1e-3, "at": "256"})),
    # ---- P0b: homogeneous deformation of a free block ------------------------------------
    **{f"P0b_{drv}": dict(tier="P0b", ic=("block", {"drive": drv, "rate": 0.1, "half": 1.0}),
                          rungs=["0.05", "0.025", "0.0125"], t_end=0.55, dump_dt=0.02,
                          box=(-3., -3., 3., 3.), group="std", material={},
                          variants=["pref", "pns", "prr", "ptc0", "pps", "gz"],
                          check=("p0b", {"drive": drv, "rate": 0.1, "inner": 0.25,
                                         "tol_elastic": 0.01, "tol_plateau": 0.02, "tol_plateau_prr": 0.01,
                                         "tol_scatter": 0.02, "tol_gz": 0.02, "tol_tc": 0.005}))
       for drv in ("uni", "psh", "ssh")},
    # ---- P2: plate impact, elastic precursor + plastic shock ------------------------------
    "P2_plate": dict(tier="P2", ic=("slabs", {"L": 4.0, "H": 0.2, "up": 0.1}),
                     rungs=["0.02", "0.01", "0.005"], t_end=3.0, dump_dt=0.1,
                     box=(-8., 0., 8., 0.2), group="per", material={},
                     variants=["pref", "pns", "prr", "ptc0", "pps", "gz"],
                     check=("p2", {"up": 0.1, "t_fit": (1.0, 3.0), "tol_U": 0.01, "tol_hel": 0.02,
                                   "tol_sigma": 0.01, "margin_gz": 0.005})),
    # ---- P3: thin flyer, release + attenuation (1D reference code) -----------------------
    # Lf = 0.2: in the 1D reference the flyer's release catches the plastic shock at X ~ 1.3-1.8 and
    # attenuates it to ~0.07 by X ~ 2.9 within t = 3 (Lf = 0.5 would need t ~ 6: no attenuation in the domain).
    "P3_flyer": dict(tier="P3", ic=("flyer", {"Lf": 0.2, "Lt": 4.0, "H": 0.2, "vf": 0.2}),
                     rungs=["0.02", "0.01", "0.005"], t_end=3.0, dump_dt=0.05,
                     box=(-6., 0., 10., 0.2), group="per", material={},   # same aspect as P2: one GIZMO build
                     variants=["pref", "pns", "prr", "ptc0", "gz"],
                     check=("p3", {"vf": 0.2, "Lf": 0.2, "Lt": 4.0, "gauges": (0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
                                   "arrival_min": 1.0,
                                   "tol_peak": 0.03, "tol_arrival": 0.01, "margin_gz": 0.005,
                                   "ref_cells": 20000})),
    # ---- P4: Taylor bars (symmetric bar impact) ----------------------------------------
    "P4_taylor": dict(tier="P4", ic=("bars", {"L": 2.0, "W": 0.5, "v": 0.1}),
                      rungs=["0.05", "0.025", "0.0125", "0.00625"], t_end=12.0, dump_dt=0.25,
                      box=(-4., -2., 4., 2.), group="std", material={},
                      variants=["pref", "pns", "prr", "ptc0", "gz"],
                      check=("p4", {"L": 2.0, "W": 0.5, "v": 0.1, "tol_Lconv": 0.01, "tol_Wconv": 0.03,
                                    "tol_Lgz": 0.01, "tol_Wgz": 0.03})),
}

# Gate order (plan P.3). `./pt status` shows it; `./pt next` names the first open gate.
ORDER = ["P0a", "P1a", "P1d", "P1b", "P1c", "P0b", "P2", "P3", "P4"]
GATES = {"P0a", "P1a", "P1d", "P0b"}
