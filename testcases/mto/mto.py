#!/usr/bin/env python3
"""mto -- multi-testcase-optimizing harness.

One command line for the whole testplan (vault: 2026-10-02-mto-testplan.md):

  ./mto.py build  CAND [-v VARIANT ...]     build the variants a candidate needs
  ./mto.py run    CAND -t S [-t T ...] [-j 4] [--stop-on-fail] [--force]
                  CAND --only T1:d0.1 W2:n64  run single rungs
  ./mto.py status [CAND]                    progress of every run (ETA from the log)
  ./mto.py analyze CAND [--tests T1 W2]     (re)analyse finished runs -> result.json
  ./mto.py score  [CAND ...]                scoreboard, one verdict per test
  ./mto.py frames CAND TEST [RUNG]          PNG montage of a body run (or all rungs)
  ./mto.py calib                            seconds per particle-step from finished runs
  ./mto.py plan   CAND -t S -t T            dry run: jobs, N, estimated cost

Candidates are flag files, candidates/<name>.flags, applied on top of flags/base.flags;
the testcase variant (periodicity, EOS, plasticity; lib/tests.py) is applied last.
All code changes live behind flags in ONE source tree (demonstrator/), so a candidate
is fully described by its flag file plus the source checksum recorded in every run.

Layout:
  builds/<cand>/<variant>/   fresh source copy + make; STAMP = md5(source + flags)
  ics/<test>_<rung>.h5       cached ICs (generated on first use)
  runs/<cand>/<test>/<rung>/ config.info, run.log, output/, exit_code, meta.json, result.json
  runs/<cand>/queue.log      one line per finished job
"""

import argparse
import fcntl
from contextlib import contextmanager
import concurrent.futures as cf
import glob
import hashlib
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "lib"))

import tests as T            # noqa: E402

DEMO = os.path.normpath(os.path.join(HERE, "..", "..", "demonstrator"))
BUILDS, ICS, RUNS = (os.path.join(HERE, d) for d in ("builds", "ics", "runs"))
CALIB = os.path.join(HERE, "calib.json")
LOG_MAX = 1.5e9              # kill a run whose log exceeds this (error floods)
C_PS_DEFAULT = 3.5e-5        # s per particle-step, Mars calibration (09-22)


# ---------------------------------------------------------------- flags / builds

def read_flags(path):
    out = {}
    if not os.path.exists(path):
        raise SystemExit("no such flag file: " + path)
    for line in open(path):
        line = line.split("#", 1)[0].strip()
        if line and not line.startswith("@"):
            k, v = line.split(None, 1)
            out[k] = v.strip()
    return out


BODY_TESTS = ("T1", "T2", "P1", "A1")


def mat_overrides(cand):
    """Runtime material overrides of a candidate: '@mat KEY VALUE [TEST ...]' lines in its
    flag file (default tests: the body tests, never the waves, which stay inviscid/elastic).
    They change the config text, hence the config hash, so affected rungs rerun."""
    out = []
    path = os.path.join(HERE, "candidates", cand + ".flags")
    for line in open(path):
        line = line.split("#", 1)[0].split()
        if line and line[0] == "@mat":
            out.append((line[1], line[2], tuple(line[3:]) or BODY_TESTS))
    return out


def rung_for(cand, test, tag):
    """The rung as the candidate runs it: lib/tests.py + the candidate's @mat overrides."""
    r = T.rung(test, tag)
    ov = [(k, v) for k, v, tests in mat_overrides(cand) if test in tests]
    if ov:
        r["materials"] = [dict(m, **{k: float(v) for k, v in ov}) for m in r["materials"]]
    return r


def cand_flags(cand, variant):
    f = read_flags(os.path.join(HERE, "flags", "base.flags"))
    f.update(read_flags(os.path.join(HERE, "candidates", cand + ".flags")))
    f.update({k: str(v) for k, v in T.VARIANTS[variant].items()})
    return f


def source_hash():
    h = hashlib.md5()
    for sub in ("src", "include"):
        for name in sorted(os.listdir(os.path.join(DEMO, sub))):
            if name.endswith((".cpp", ".h")) and name != "parameter.h":
                h.update(open(os.path.join(DEMO, sub, name), "rb").read())
    h.update(open(os.path.join(DEMO, "Makefile"), "rb").read())
    return h.hexdigest()


def set_flag(text, name, value):
    pat = re.compile(r"^#define {} .*$".format(re.escape(name)), re.M)
    if not pat.search(text):
        raise SystemExit("flag {} not found in parameter.h".format(name))
    return pat.sub("#define {} {}".format(name, value), text, count=1)


@contextmanager
def _lock(path):
    """Exclusive inter-process lock (parallel queue jobs share builds and ICs)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def _built(B, stamp):
    st = os.path.join(B, "STAMP")
    return (os.path.exists(os.path.join(B, "bin", "mlh")) and os.path.exists(st)
            and open(st).read().strip() == stamp)


def build(cand, variant, quiet=False):
    flags = cand_flags(cand, variant)
    stamp = hashlib.md5((source_hash() + json.dumps(flags, sort_keys=True)
                         + open(os.path.join(DEMO, "include", "parameter.h")).read())
                        .encode()).hexdigest()
    B = os.path.join(BUILDS, cand, variant)
    if _built(B, stamp):
        return B, stamp
    with _lock(os.path.join(BUILDS, cand, variant + ".lock")):
        if _built(B, stamp):                     # another job built it while we waited
            return B, stamp
        return _build_locked(cand, variant, flags, stamp, B, quiet)


def _build_locked(cand, variant, flags, stamp, B, quiet):
    binp = os.path.join(B, "bin", "mlh")
    if not quiet:
        print("build {}/{} ...".format(cand, variant), flush=True)
    shutil.rmtree(B, ignore_errors=True)
    os.makedirs(B)
    for item in ("Makefile", "include", "src", "tests"):
        s = os.path.join(DEMO, item)
        if os.path.isdir(s):
            shutil.copytree(s, os.path.join(B, item))
        elif os.path.exists(s):
            shutil.copy(s, B)
    os.symlink(os.path.join(DEMO, "holib"), os.path.join(B, "holib"))
    P = os.path.join(B, "include", "parameter.h")
    text = open(P).read()
    for k, v in flags.items():
        text = set_flag(text, k, v)
    open(P, "w").write(text)
    with open(os.path.join(B, "SNAPSHOT.txt"), "w") as f:
        f.write("candidate {} variant {}\nsource md5 {}\ngit {}\n{}\n\nflags:\n".format(
            cand, variant, source_hash(),
            subprocess.run(["git", "-C", DEMO, "rev-parse", "HEAD"], capture_output=True,
                           text=True).stdout.strip(), time.strftime("%Y-%m-%d %H:%M")))
        for k, v in flags.items():
            f.write("  {} {}\n".format(k, v))
    r = subprocess.run(["make", "-j{}".format(os.cpu_count())], cwd=B,
                       stdout=open(os.path.join(B, "make.log"), "w"), stderr=subprocess.STDOUT)
    if r.returncode or not os.path.exists(binp):
        errs = [l for l in open(os.path.join(B, "make.log")) if "error" in l][:10]
        raise SystemExit("build {}/{} failed:\n{}".format(cand, variant, "".join(errs)))
    open(os.path.join(B, "STAMP"), "w").write(stamp)
    return B, stamp


# ---------------------------------------------------------------- ICs / configs

def ic_path(test, tag):
    return os.path.join(ICS, "{}_{}.h5".format(test, tag))


def make_ic(test, tag):
    import gen
    import h5py
    r = T.rung(test, tag)
    path = ic_path(test, tag)
    os.makedirs(ICS, exist_ok=True)
    with _lock(os.path.join(ICS, ".lock")):
        return _make_ic_locked(test, tag, r, path, gen, h5py)


def _make_ic_locked(test, tag, r, path, gen, h5py):
    if not os.path.exists(path):
        if r["gen"] == "wave":
            gen.wave(path, **r["genargs"])
        else:
            getattr(gen, r["gen"])(path, r["d"], **r["genargs"])
    with h5py.File(path, "r") as f:
        return len(f["m"])


SML = """sml {
    NNNTarget 32
    tol       1e-3
    maxIter   200
    maxFactor 1.2
    hMinFactor 0.1
    hMaxFactor 3.0
    warnFraction  0.05
    panicFraction 0.25
}
"""


def config_text(r, ic, outdir):
    if r["box"] is None:            # periodic wave box
        lo, hi = (0.0, 0.0), (1.0, 0.25)
    else:
        lo, hi = (-r["box"], -r["box"]), (r["box"], r["box"])
    s = ["; mto {} {} (generated by mto.py)".format(r["test"], r["tag"]),
         "initFile {}".format(ic), "outDir   {}".format(outdir),
         "timeStep {:.12g}".format(r["dump"]), "timeEnd  {:.12g}".format(r["t_end"]),
         "h5DumpInterval 1",
         "periodicBoxLimits {",
         "    lowerX {}".format(lo[0]), "    lowerY {}".format(lo[1]), "    lowerZ 0.0",
         "    upperX {}".format(hi[0]), "    upperY {}".format(hi[1]), "    upperZ 0.0", "}",
         "kernelSize {:.9g}".format(r["kernel"]), SML,
         "hydro_gamma 1.6666666666666667"]
    if r["materials"]:
        s.append("materials {")
        for k, m in enumerate(r["materials"]):
            s.append("    mat{} {{".format(k))
            s += ["        {} {}".format(a, b) for a, b in m.items()]
            s.append("    }")
        s.append("}")
    return "\n".join(s) + "\n"


# ---------------------------------------------------------------- cost model

def calib():
    try:
        return json.load(open(CALIB))
    except (OSError, ValueError):
        return {}


def est_steps(r):
    return r["t_end"] / (0.2 * r["kernel"] / (2.0 * r["c"]))


def est_cost(r, N):
    return N * est_steps(r) * calib().get(r["variant"], calib().get("all", C_PS_DEFAULT))


# ---------------------------------------------------------------- running

def cfg_hash(r, test, tag, D):
    """Identity of a rung's runtime setup: the config text (material, t_end, box, ...)."""
    return hashlib.md5(config_text(r, ic_path(test, tag), os.path.join(D, "output"))
                       .encode()).hexdigest()


def run_dir(cand, test, tag):
    return os.path.join(RUNS, cand, test, tag)


def jobs_for(args):
    pairs = []
    for tier in args.tier or []:
        pairs += T.TIERS[tier]
    for spec in args.only or []:
        test, tag = spec.split(":")
        pairs.append((test, tag))
    seen, out = set(), []
    for p in pairs:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def last_time(log):
    """Last 'TIME: t, STEP: n' of a run log (reads the tail only)."""
    try:
        with open(log, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - 20000))
            tail = f.read().decode(errors="replace")
    except OSError:
        return None, None
    m = re.findall(r"TIME: ([0-9.eE+-]+), STEP: (\d+)", tail)
    return (float(m[-1][0]), int(m[-1][1])) if m else (None, None)


def run_one(cand, test, tag, stamp, B, stop_flag):
    import analysis
    r = rung_for(cand, test, tag)
    D = run_dir(cand, test, tag)
    if stop_flag.is_set():
        return test, tag, "skipped", None
    shutil.rmtree(D, ignore_errors=True)
    os.makedirs(os.path.join(D, "output"))
    N = make_ic(test, tag)
    open(os.path.join(D, "config.info"), "w").write(
        config_text(r, ic_path(test, tag), os.path.join(D, "output")))
    meta = dict(cand=cand, test=test, tag=tag, stamp=stamp, source=source_hash(), N=N,
                t_end=r["t_end"], d=r["d"], variant=r["variant"], est_s=est_cost(r, N),
                cfg=cfg_hash(r, test, tag, D), materials=r["materials"],
                start=time.time())
    json.dump(meta, open(os.path.join(D, "meta.json"), "w"), indent=1)
    log = os.path.join(D, "run.log")
    with open(log, "w") as lf:
        p = subprocess.Popen([os.path.join(B, "bin", "mlh"), "-c", "config.info"], cwd=D,
                             stdout=lf, stderr=subprocess.STDOUT, start_new_session=True)
        open(os.path.join(D, "pid"), "w").write(str(p.pid))
        killed = None
        while p.poll() is None:
            time.sleep(5)
            if os.path.getsize(log) > LOG_MAX:
                killed = "log>{:.1g}B".format(LOG_MAX)
            if killed:
                os.killpg(p.pid, signal.SIGTERM)
                p.wait()
    meta.update(end=time.time(), exit=p.returncode, killed=killed)
    meta["wall_s"] = meta["end"] - meta["start"]
    meta["t_last"], meta["steps"] = last_time(log)
    json.dump(meta, open(os.path.join(D, "meta.json"), "w"), indent=1)
    open(os.path.join(D, "exit_code"), "w").write(str(p.returncode))
    try:
        res = analysis.analyse_run(D)
        verdict = res.get("verdict", "?")
    except Exception as e:                      # analysis must never kill the queue
        verdict = "analysis-error: {}".format(e)
    with open(os.path.join(RUNS, cand, "queue.log"), "a") as q:
        q.write("{} {} {} {} exit={} wall={:.0f}s est={:.0f}s {}\n".format(
            time.strftime("%H:%M:%S"), cand, test, tag, p.returncode, meta["wall_s"],
            meta["est_s"], verdict))
    return test, tag, verdict, meta


def cmd_run(args):
    pairs = jobs_for(args)
    if not pairs:
        raise SystemExit("nothing to run: give -t TIER or --only TEST:RUNG")
    variants = sorted({T.rung(t, g)["variant"] for t, g in pairs})
    built = {v: build(args.cand, v) for v in variants}
    todo = []
    for test, tag in pairs:
        r = rung_for(args.cand, test, tag)
        B, stamp = built[r["variant"]]
        D = run_dir(args.cand, test, tag)
        mp = os.path.join(D, "meta.json")
        if not args.force and os.path.exists(os.path.join(D, "exit_code")):
            m = json.load(open(mp))
            same_build = m.get("stamp") == stamp or args.keep_finished
            if (same_build and m.get("exit") == 0
                    and m.get("cfg") == cfg_hash(r, test, tag, D)):
                continue
        N = make_ic(test, tag)
        todo.append((est_cost(r, N), test, tag, B, stamp))
    todo.sort(reverse=True)                      # longest first: critical path starts early
    os.makedirs(os.path.join(RUNS, args.cand), exist_ok=True)
    tot = sum(c for c, *_ in todo)
    print("{} jobs, est. {:.0f} min serial, -j {}".format(len(todo), tot / 60, args.jobs))
    stop = threading.Event()
    with cf.ThreadPoolExecutor(args.jobs) as ex:
        futs = [ex.submit(run_one, args.cand, t, g, s, B, stop) for _, t, g, B, s in todo]
        for fu in cf.as_completed(futs):
            test, tag, verdict, meta = fu.result()
            print("{} {}:{} -> {}".format(time.strftime("%H:%M:%S"), test, tag, verdict),
                  flush=True)
            if args.stop_on_fail and verdict == "FAIL":
                stop.set()
    if stop.is_set():
        print("stopped launching after a FAIL (--stop-on-fail)")


def cmd_build(args):
    for v in args.variant or sorted(T.VARIANTS):
        B, stamp = build(args.cand, v)
        print("{}/{}: {} ({})".format(args.cand, v, B, stamp[:8]))


def cmd_plan(args):
    tot = 0
    for test, tag in jobs_for(args):
        r = rung_for(args.cand, test, tag)
        N = make_ic(test, tag)
        c = est_cost(r, N)
        tot += c
        print("{:5s} {:8s} {:8s} N={:6d} steps~{:7.0f} est {:7.0f} s".format(
            test, tag, r["variant"], N, est_steps(r), c))
    print("total {:.0f} min serial".format(tot / 60))


def candidates_with_runs():
    """Candidate directories under runs/ (logs and other files living there are skipped)."""
    if not os.path.isdir(RUNS):
        return []
    return sorted(d for d in os.listdir(RUNS) if os.path.isdir(os.path.join(RUNS, d)))


def cmd_status(args):
    cands = [args.cand] if args.cand else candidates_with_runs()
    for cand in cands:
        base = os.path.join(RUNS, cand)
        for test in sorted(os.listdir(base)):
            if not os.path.isdir(os.path.join(base, test)):
                continue
            for tag in sorted(os.listdir(os.path.join(base, test))):
                D = os.path.join(base, test, tag)
                try:
                    m = json.load(open(os.path.join(D, "meta.json")))
                except (OSError, ValueError):
                    continue
                t, step = last_time(os.path.join(D, "run.log"))
                if "exit" in m:
                    res = os.path.join(D, "result.json")
                    v = json.load(open(res)).get("verdict", "?") if os.path.exists(res) else "-"
                    st = "done exit={} {:.0f}s {}".format(m["exit"], m["wall_s"], v)
                else:
                    el = time.time() - m["start"]
                    frac = (t or 0) / m["t_end"]
                    eta = el / frac - el if frac > 0 else float("nan")
                    st = "run  t={:.4g}/{:.4g} ({:.0%}) {:.0f}s, ETA {:.0f}s".format(
                        t or 0, m["t_end"], frac, el, eta)
                print("{:10s} {:5s} {:8s} {}".format(cand, test, tag, st))


def cmd_analyze(args):
    import analysis
    base = os.path.join(RUNS, args.cand)
    for test in sorted(os.listdir(base)):
        if args.tests and test not in args.tests:
            continue
        D0 = os.path.join(base, test)
        if not os.path.isdir(D0):
            continue
        for tag in sorted(os.listdir(D0)):
            D = os.path.join(D0, tag)
            if os.path.exists(os.path.join(D, "exit_code")):
                res = analysis.analyse_run(D)
                print("{:5s} {:8s} {}  {}".format(test, tag, res.get("verdict"),
                                                  res.get("summary", "")))


def cmd_score(args):
    import analysis
    cands = args.cands or candidates_with_runs()
    print(analysis.scoreboard(RUNS, cands))


def cmd_frames(args):
    import analysis
    tags = [args.rung] if args.rung else sorted(os.listdir(os.path.join(RUNS, args.cand,
                                                                        args.test)))
    for tag in tags:
        png = analysis.frames(run_dir(args.cand, args.test, tag), n=args.n,
                              field=args.field)
        print(png)


def cmd_calib(args):
    # runs/<cand>/calib_exclude.txt: "TEST RUNG" lines of runs whose wall time is not a
    # clean measurement (e.g. paused with SIGSTOP while running)
    excl = set()
    for f in glob.glob(os.path.join(RUNS, "*", "calib_exclude.txt")):
        cand = os.path.basename(os.path.dirname(f))
        excl |= {(cand,) + tuple(l.split()) for l in open(f) if l.strip()}
    rates = {}
    for root, dirs, files in os.walk(RUNS):
        if "meta.json" in files:
            m = json.load(open(os.path.join(root, "meta.json")))
            if (m.get("cand"), m.get("test"), m.get("tag")) in excl:
                continue
            if m.get("exit") == 0 and m.get("steps") and m["wall_s"] > 20:
                rates.setdefault(m["variant"], []).append(m["wall_s"] / (m["N"] * m["steps"]))
    out = {v: sorted(r)[len(r) // 2] for v, r in rates.items()}
    if out:
        out["all"] = sorted(sum(rates.values(), []))[len(sum(rates.values(), [])) // 2]
    json.dump(out, open(CALIB, "w"), indent=1)
    for v, c in sorted(out.items()):
        print("{:8s} {:.3e} s/particle-step ({} runs)".format(v, c, len(rates.get(v, [0]))))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("build"); p.add_argument("cand"); p.add_argument("-v", "--variant", action="append")
    p.set_defaults(fn=cmd_build)
    for name, fn in (("run", cmd_run), ("plan", cmd_plan)):
        p = sp.add_parser(name); p.add_argument("cand")
        p.add_argument("-t", "--tier", action="append", choices=sorted(T.TIERS))
        p.add_argument("--only", nargs="+", metavar="TEST:RUNG")
        p.add_argument("-j", "--jobs", type=int, default=4)
        p.add_argument("--force", action="store_true", help="rerun finished jobs")
        p.add_argument("--stop-on-fail", action="store_true")
        p.add_argument("--keep-finished", action="store_true",
                       help="skip finished runs with the same config even if the build stamp "
                            "changed -- only after a source change that is default-off and "
                            "verified bit-identical (e.g. candidate ref_chk)")
        p.set_defaults(fn=fn)
    p = sp.add_parser("status"); p.add_argument("cand", nargs="?"); p.set_defaults(fn=cmd_status)
    p = sp.add_parser("analyze"); p.add_argument("cand"); p.add_argument("--tests", nargs="+")
    p.set_defaults(fn=cmd_analyze)
    p = sp.add_parser("score"); p.add_argument("cands", nargs="*"); p.set_defaults(fn=cmd_score)
    p = sp.add_parser("frames"); p.add_argument("cand"); p.add_argument("test")
    p.add_argument("rung", nargs="?"); p.add_argument("-n", type=int, default=6)
    p.add_argument("--field", default="rho"); p.set_defaults(fn=cmd_frames)
    p = sp.add_parser("calib"); p.set_defaults(fn=cmd_calib)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
