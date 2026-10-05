#!/usr/bin/env python3
"""GIZMO control for the rotating disc A1 (angular momentum), half or full rotation, optionally MPI.

    ./disc.py d0.005 [--np 3] [--frac 0.5]      IC from ../ics/A1_<rung>.h5, run, convert

Same GIZMO build, parameter file and snapshot conversion as rod.py (it imports it). The disc is purely
elastic (Tillotson, mu = MU_BODY, Y0 = 1e20 -> no yielding, no plasticity) in the mto 'body' setup;
t_end = frac * t_end(A1) (frac 0.5 = half a rotation), snapshot spacing as in A1.
Run directory gizmo/runs/disc_<rung>_f<frac> (params file named rod.params so the queue counts it).
"""
import argparse, glob, os, subprocess, sys
import h5py
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rod
import tests as T


def setup(tag, frac):
    r = T.rung("A1", tag)
    src = os.path.join(HERE, "..", "ics", "A1_{}.h5".format(tag))
    if not os.path.exists(src):   # unstressed lattice: GIZMO zeroes S at start, so no seeded equilibrium (mto.py's A1_<tag>_eq*)
        import gen
        gen.disc(src, r["d"], **r["genargs"])
    D = os.path.join(HERE, "runs", "disc_{}_f{:g}".format(tag, frac))
    os.makedirs(os.path.join(D, "output"), exist_ok=True)
    with h5py.File(src, "r") as f:
        x, v, m, u = f["x"][:], f["v"][:], f["m"][:], f["u"][:]
    N = len(m)
    pos = np.zeros((N, 3)); pos[:, :2] = x + rod.OFFSET
    vel = np.zeros((N, 3)); vel[:, :2] = v
    with h5py.File(os.path.join(D, "ic.hdf5"), "w") as f:
        h = f.create_group("Header")
        npart = np.array([N, 0, 0, 0, 0, 0], dtype=np.int32)
        h.attrs["NumPart_ThisFile"] = npart
        h.attrs["NumPart_Total"] = npart.astype(np.uint32)
        h.attrs["NumPart_Total_HighWord"] = np.zeros(6, dtype=np.uint32)
        h.attrs["MassTable"] = np.zeros(6)
        for k, val in (("Time", 0.0), ("Redshift", 0.0), ("BoxSize", 2 * rod.OFFSET), ("NumFilesPerSnapshot", 1),
                       ("Omega0", 0.0), ("OmegaLambda", 0.0), ("HubbleParam", 1.0), ("Flag_Sfr", 0), ("Flag_Cooling", 0),
                       ("Flag_StellarAge", 0), ("Flag_Metals", 0), ("Flag_Feedback", 0), ("Flag_DoublePrecision", 1)):
            h.attrs[k] = val
        g = f.create_group("PartType0")
        g.create_dataset("Coordinates", data=pos)
        g.create_dataset("Velocities", data=vel)
        g.create_dataset("Masses", data=m)
        g.create_dataset("InternalEnergy", data=u)
        g.create_dataset("ParticleIDs", data=np.arange(1, N + 1, dtype=np.uint32))
        g.create_dataset("CompositionType", data=np.zeros(N, dtype=np.int32))
    mat = r["materials"][0]
    text = rod.PARAMS.format(tag="A1 " + tag, dump=r["dump"] * min(1.0, 2 * frac), t_end=frac * r["t_end"], box=2 * rod.OFFSET,
                             mu=mat["mu"], Y0=1e20)
    text = text.replace("TimeLimitCPU                100000", "TimeLimitCPU                1e9")
    open(os.path.join(D, "rod.params"), "w").write(text)
    return D, N


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tag")
    ap.add_argument("--np", type=int, default=1, help="MPI ranks")
    ap.add_argument("--frac", type=float, default=0.5, help="fraction of the A1 end time (0.5 = half a rotation)")
    ap.add_argument("--convert-only", action="store_true")
    a = ap.parse_args()
    if a.convert_only:
        print(rod.convert(os.path.join(HERE, "runs", "disc_{}_f{:g}".format(a.tag, a.frac))))
        sys.exit()
    D, N = setup(a.tag, a.frac)
    cmd = (["mpirun", "-np", str(a.np)] if a.np > 1 else []) + [rod.GIZMO, "rod.params"]
    with open(os.path.join(D, "run.log"), "w") as lf:
        rc = subprocess.run(cmd, cwd=D, stdout=lf, stderr=subprocess.STDOUT).returncode
    out = rod.convert(D)
    print("disc {} (N={}, np={}): exit {}, {} snapshots".format(a.tag, N, a.np, rc, len(glob.glob(os.path.join(out, "*.h5")))))
