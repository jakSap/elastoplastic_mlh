"""GIZMO MFM arm (plan P.5): build, IC conversion, params, run, snapshot conversion.

One GIZMO build per geometry: non-periodic ("gizmo_std") or periodic with the test's box aspect
("gizmo_per_<Ly/Lx>"; BOX_LONG_Y is compile-time). Snapshots are converted into the demonstrator
layout (output/NNNNNN.h5: time, x, v, m, rho, u, P, sml, Sxx, Sxy, Syy, ids), sorted by ParticleIDs
and shifted back to the demonstrator frame, so every checker reads both codes the same way."""
import bisect
import glob
import json
import os
import re
import shutil
import subprocess

import h5py
import numpy as np

CONFIG = """HYDRO_MESHLESS_FINITE_MASS
BOX_SPATIAL_DIMENSION=2
KERNEL_FUNCTION=6
EOS_TILLOTSON
EOS_ELASTIC
SELFGRAVITY_OFF
PREVENT_PARTICLE_MERGE_SPLIT
INPUT_IN_DOUBLEPRECISION
OUTPUT_IN_DOUBLEPRECISION
"""


def geometry(S, test):
    T = S.M.TESTS[test]
    if S.M.BUILD_GROUPS[T["group"]].get("PERIODIC_BOUNDARIES"):
        x0, y0, x1, y1 = T["box"]
        return "per_%g" % ((y1 - y0) / (x1 - x0)), dict(periodic=True, Lx=x1 - x0, aspect=(y1 - y0) / (x1 - x0),
                                                        shift=(-x0, -y0))
    return "std", dict(periodic=False, shift=(0.0, 0.0))


def build_dir(S, key):
    return S.p("gizmo_" + key)


def _config(geo):
    c = CONFIG
    if geo["periodic"]:
        c += "BOX_PERIODIC\nBOX_LONG_X=1\nBOX_LONG_Y=%r\n" % geo["aspect"]
    return c


def build_ok(S, key_or_group, test=None):
    B = build_dir(S, key_or_group)
    if not os.path.exists(os.path.join(B, "GIZMO_BUILD.json")):
        return False, "missing"
    if not os.path.exists(os.path.join(B, "src", "GIZMO")):
        return False, "no binary (make failed?)"
    return True, "ok"


def build(S, key, geo, force=False):
    B = build_dir(S, key)
    ok, _ = build_ok(S, key)
    if ok and not force:
        print("  %-22s up to date" % ("gizmo_" + key))
        return True
    shutil.rmtree(B, ignore_errors=True)
    os.makedirs(B)
    subprocess.run(["rsync", "-a", "--exclude", "*.o", "--exclude", "GIZMO", "--exclude", "GIZMO_SPH",
                    "--exclude", "*.hdf5", "--exclude", "output*", "--exclude", "*.params*",
                    S.M.GIZMO_SRC + "/", os.path.join(B, "src") + "/"], check=True)
    open(os.path.join(B, "src", "Config.sh"), "w").write(_config(geo))
    r = subprocess.run("make clean > /dev/null 2>&1; make -j6 > ../build.log 2>&1", shell=True,
                       cwd=os.path.join(B, "src"))
    if r.returncode or not os.path.exists(os.path.join(B, "src", "GIZMO")):
        print("  %-22s FAILED (see %s/build.log)" % ("gizmo_" + key, B))
        return False
    head = subprocess.run(["git", "-C", S.M.GIZMO_SRC, "rev-parse", "--short", "HEAD"], text=True,
                          stdout=subprocess.PIPE).stdout.strip()
    json.dump(dict(source=S.M.GIZMO_SRC, git=head, geometry=geo, config=_config(geo)),
              open(os.path.join(B, "GIZMO_BUILD.json"), "w"), indent=1)
    shutil.copy(os.path.join(B, "src", "Config.sh"), os.path.join(B, "Config.sh"))
    print("  %-22s built  (GIZMO %s)" % ("gizmo_" + key, head))
    return True


def demo_to_gizmo(src, dst, ks, geo, u_floor=None):
    with h5py.File(src, "r") as f:
        x, v, m, u = f["x"][()], f["v"][()], f["m"][()], f["u"][()]
        S = [f[k][()] if k in f else np.zeros(len(m)) for k in ("Sxx", "Sxy", "Syy")]
    N = len(m)
    pos = np.zeros((N, 3))
    pos[:, :2] = x + np.array(geo["shift"])
    vel = np.zeros((N, 3))
    vel[:, :2] = v
    st = np.zeros((N, 9))
    st[:, 0], st[:, 1], st[:, 3], st[:, 4] = S[0], S[1], S[1], S[2]
    with h5py.File(dst, "w") as o:
        h = o.create_group("Header")
        npart = np.array([N, 0, 0, 0, 0, 0], dtype=np.int32)
        for k, val in (("NumPart_ThisFile", npart), ("NumPart_Total", npart),
                       ("NumPart_Total_HighWord", np.zeros(6, dtype=np.int32)), ("MassTable", np.zeros(6)),
                       ("Time", 0.0), ("Redshift", 0.0), ("BoxSize", float(geo.get("Lx", 0.0))),
                       ("NumFilesPerSnapshot", 1), ("Omega0", 1.0), ("OmegaLambda", 0.0), ("HubbleParam", 1.0),
                       ("Flag_DoublePrecision", 1)):
            h.attrs[k] = val
        for k in ("Flag_Sfr", "Flag_Cooling", "Flag_StellarAge", "Flag_Metals", "Flag_Feedback",
                  "Flag_IC_Info", "Flag_AgeTracers"):
            h.attrs[k] = 0
        p = o.create_group("PartType0")
        p["Coordinates"], p["Velocities"], p["Masses"] = pos, vel, m.astype(np.float64)
        p["InternalEnergy"] = u.astype(np.float64)
        p["Density"] = np.ones(N)
        p["SmoothingLength"] = np.full(N, ks)
        p["ParticleIDs"] = np.arange(1, N + 1, dtype=np.uint32)
        p["CompositionType"] = np.zeros(N, dtype=np.uint32)
        p["StressTensor"] = st
        p["ParticleChildIDsNumber"] = np.zeros(N, dtype=np.uint32)
        p["ParticleIDGenerationNumber"] = np.zeros(N, dtype=np.uint32)


def setup_run(S, test, rung, rel, ic_rel, ks, meta):
    from .core import snap_t_end
    T = S.M.TESTS[test]
    key, geo = geometry(S, test)
    d = S.p("runs", rel)
    open(os.path.join(d, "BUILD"), "w").write("gizmo_%s\n" % key)
    demo_to_gizmo(S.p(ic_rel), os.path.join(d, "ic.hdf5"), ks, geo)
    json.dump(geo, open(os.path.join(d, "gizmo_geometry.json"), "w"))
    mat = {**S.M.MATERIAL, **T.get("material", {})}
    te, dt, n = snap_t_end(T["t_end"], T["dump_dt"])
    par = [("InitCondFile", "ic"), ("OutputDir", "gizmo_out"), ("TimeBegin", 0), ("TimeMax", te),
           ("TimeBetSnapshot", dt), ("TimeOfFirstSnapshot", 0), ("MaxSizeTimestep", dt),
           ("CourantFac", 0.2), ("DesNumNgb", S.M.NNN)]
    if geo["periodic"]:
        par.append(("BoxSize", geo["Lx"]))
    par += [("Tillotson_EOS_params_rho_0", mat["rho0"]), ("Tillotson_EOS_params_A", mat["A"]),
            ("Tillotson_EOS_params_B", mat["B"]), ("Tillotson_EOS_params_a", mat["a"]),
            ("Tillotson_EOS_params_b", mat["b"]), ("Tillotson_EOS_params_u_0", mat["u0"]),
            ("Tillotson_EOS_params_u_s", mat["u_iv"]), ("Tillotson_EOS_params_u_s_prime", mat["u_cv"]),
            ("Tillotson_EOS_params_alpha", mat["alpha"]), ("Tillotson_EOS_params_beta", mat["beta"]),
            ("Tillotson_EOS_params_mu", mat["mu"]), ("Tillotson_EOS_params_Y0", mat["Y0"]),
            ("UnitLength_in_cm", 1), ("UnitMass_in_g", 1), ("UnitVelocity_in_cm_per_s", 1)]
    par += [("Softening_Type%d" % i, 0.001) for i in range(6)]
    open(os.path.join(d, "gizmo.params"), "w").write("".join("%-32s %r\n" % (k, v) if not isinstance(v, str)
                                                             else "%-32s %s\n" % (k, v) for k, v in par))
    open(os.path.join(d, "NP"), "w").write("%d\n" % (2 if meta.get("N", 0) > 20000 else 1))


def run_command(S, rel):
    d = S.p("runs", rel)
    tag = open(os.path.join(d, "BUILD")).read().strip()
    np_ = int(open(os.path.join(d, "NP")).read())
    return ["mpirun", "-np", str(np_), S.p(tag, "src", "GIZMO"), "gizmo.params"]


def convert_output(S, rel):
    d = S.p("runs", rel)
    geo = json.load(open(os.path.join(d, "gizmo_geometry.json")))
    files = sorted(glob.glob(os.path.join(d, "gizmo_out", "snapshot_*.hdf5")))
    if not files:
        return 3
    shift = np.array(geo["shift"])
    # GIZMO writes a snapshot inside find_next_sync_point_and_drift (run.c): positions are drifted to the
    # output time, but the kicked fields (Elastic_Stress_Tensor, InternalEnergy) only had the opening half
    # kick of the step [T_k, T_k+1] that contains the output time, so they belong to (T_k + T_k+1)/2.
    # Recorded as /time_S; the checkers compare stresses at time_S.
    sync = sorted(float(m.group(1)) for m in re.finditer(r"Sync-Point \d+, Time: ([-0-9.eE+]+)",
                                                         open(os.path.join(d, "run.log"), errors="replace").read()))

    def stress_time(t):
        # the first snapshot (t = T_0) is also written after the opening half kick of step 1
        k = max(1, bisect.bisect_left(sync, t))    # sync[k-1] < t <= sync[k]
        if k >= len(sync):
            return t
        return 0.5 * (sync[k - 1] + sync[k])
    for i, f in enumerate(files):
        with h5py.File(f, "r") as g:
            t = float(g["Header"].attrs["Time"])
            p = g["PartType0"]
            o = np.argsort(p["ParticleIDs"][()])
            st = p["StressTensor"][()][o] if "StressTensor" in p else None
            out = dict(time=np.array([t]), time_S=np.array([stress_time(t)]),
                       x=p["Coordinates"][()][o][:, :2] - shift,
                       v=p["Velocities"][()][o][:, :2], m=p["Masses"][()][o], rho=p["Density"][()][o],
                       u=p["InternalEnergy"][()][o], sml=p["SmoothingLength"][()][o],
                       ids=p["ParticleIDs"][()][o])
            if "Pressure" in p:
                out["P"] = p["Pressure"][()][o]
            if st is not None:
                st = st.reshape(-1, 9)
                out["Sxx"], out["Sxy"], out["Syy"] = st[:, 0], st[:, 1], st[:, 4]
        with h5py.File(os.path.join(d, "output", "%06d.h5" % i), "w") as h:
            for k, a in out.items():
                h.create_dataset(k, data=a)
    for f in files:
        os.remove(f)
    return 0
