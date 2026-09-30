"""Study bookkeeping for the plasticity test plan: studies, builds, run dirs, queue, status.

Layout of a study (compatible with the elastic test-plan studies and binac-cpu-runs):
    <STUDY_ROOT>/plasticity_<YYYYMMDD>_<commit7>/
        STUDY.json, manifest_frozen.py      what this study is (never edited)
        build_<variant>[_<group>]/          newMFM build: SNAPSHOT.txt, bin/mlh, make.log
        gizmo_<group>/                       GIZMO build: Config.sh, GIZMO, build.log
        ics/                                 IC files (+ IC_SOURCES.txt)
        runs/<test>/<variant>_<rung>/        config.info | gizmo.params, BUILD, parameter.h, output/,
                                             run.log, exit_code
        jobs_<test>.txt                      newMFM jobs (binac_stage.sh input)
        summary.txt                          one line per finished run
        checks/<test>.txt, .json             checker verdicts
        movies/                              rendered movies
"""
import datetime
import glob
import hashlib
import importlib.util
import json
import math
import os
import shutil
import signal
import subprocess
import sys
import time

HARNESS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def die(msg):
    print("ERROR: " + msg, file=sys.stderr)
    sys.exit(2)


def sh(cmd, cwd=None, check=True, capture=True):
    r = subprocess.run(cmd, cwd=cwd, shell=isinstance(cmd, str), text=True,
                       stdout=subprocess.PIPE if capture else None,
                       stderr=subprocess.STDOUT if capture else None)
    if check and r.returncode != 0:
        die("command failed (%d): %s\n%s" % (r.returncode, cmd, r.stdout or ""))
    return (r.stdout or "").strip()


# ------------------------------------------------------------------ manifest / study
def load_manifest(path):
    spec = importlib.util.spec_from_file_location("manifest_%d" % abs(hash(path)), path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def live_manifest():
    return load_manifest(os.path.join(HARNESS, "manifest.py"))


def full_commit(M, c):
    return sh(["git", "-C", M.REPO, "rev-parse", "--verify", c + "^{commit}"])


class Study:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        self.meta = json.load(open(os.path.join(self.path, "STUDY.json")))
        self.M = load_manifest(os.path.join(self.path, "manifest_frozen.py"))

    def p(self, *a):
        return os.path.join(self.path, *a)

    # ---- tests / runs
    def tests(self, sel=None):
        """Test names matching sel: None/'all' -> all; a tier ('P2'); a test name; comma list."""
        names = list(self.M.TESTS)
        if sel in (None, "all"):
            return names
        out = []
        for s in sel.split(","):
            hit = [n for n in names if n == s or self.M.TESTS[n]["tier"] == s]
            if not hit:
                die("no test or tier '%s' (have: %s)" % (s, ", ".join(names)))
            out += [h for h in hit if h not in out]
        return out

    def runs(self, test, variants=None, rungs=None):
        T = self.M.TESTS[test]
        for v in T["variants"]:
            if variants and v not in variants:
                continue
            for r in T["rungs"]:
                if rungs and r not in rungs:
                    continue
                yield v, r, "%s/%s_%s" % (test, v, ("n" if T.get("ks_scale") == "n" else "d") + r)

    def code(self, variant):
        return self.M.VARIANTS[variant].get("code", "newmfm")

    def build_tag(self, variant, group, test=None):
        if self.code(variant) == "gizmo":
            from . import gizmo
            return "gizmo_" + gizmo.geometry(self, test)[0]
        return "build_" + (variant if group == "std" else "%s_%s" % (variant, group))


def study_name(M, date=None):
    return "plasticity_%s_%s" % (date or datetime.date.today().strftime("%Y%m%d"), M.COMMIT[:7])


def find_study(arg=None):
    """-s argument, else the newest study of the live manifest's COMMIT."""
    if arg:
        if not os.path.exists(os.path.join(arg, "STUDY.json")):
            die("%s is not a study (no STUDY.json)" % arg)
        return Study(arg)
    M = live_manifest()
    cands = sorted(glob.glob(os.path.join(M.STUDY_ROOT, "plasticity_*_%s" % M.COMMIT[:7])))
    cands = [c for c in cands if os.path.exists(os.path.join(c, "STUDY.json"))]
    if not cands:
        die("no study for COMMIT %s yet -- run: ./pt new" % M.COMMIT)
    return Study(cands[-1])


def harness_version():
    head = sh(["git", "-C", HARNESS, "rev-parse", "--short", "HEAD"], check=False)
    dirty = sh(["git", "-C", HARNESS, "status", "--short", "--", "."], check=False)
    return head + ("+dirty" if dirty else "")


def cmd_new(args):
    M = live_manifest()
    full = full_commit(M, M.COMMIT)
    anc = subprocess.run(["git", "-C", M.REPO, "merge-base", "--is-ancestor", full, M.BRANCH]).returncode
    if anc != 0:
        die("COMMIT %s is not on branch %s" % (M.COMMIT, M.BRANCH))
    path = os.path.join(M.STUDY_ROOT, study_name(M, args.date))
    if os.path.exists(path):
        die("%s exists already (one study per date and commit)" % path)
    os.makedirs(os.path.join(path, "ics"))
    shutil.copy(os.path.join(HARNESS, "manifest.py"), os.path.join(path, "manifest_frozen.py"))
    meta = dict(commit=full, commit7=M.COMMIT[:7], base_commit=full_commit(M, M.BASE_COMMIT),
                branch=M.BRANCH, created=datetime.datetime.now().isoformat(timespec="seconds"),
                harness=harness_version(), harness_dir=HARNESS,
                plan="MA-Obsidian/2026-09-30-testplan-plastisch.md")
    json.dump(meta, open(os.path.join(path, "STUDY.json"), "w"), indent=1)
    print("new study:", path)
    print("next: ./pt build   (then ./pt next)")


def cmd_refreeze(args):
    """Copy the live manifest into the study if it only ADDS tests/variants."""
    S = find_study(args.study)
    M = live_manifest()
    old = S.M
    if (M.COMMIT, M.BASE_COMMIT) != (old.COMMIT, old.BASE_COMMIT):
        die("COMMIT changed: start a new study (./pt new)")
    for k, v in old.VARIANTS.items():
        if M.VARIANTS.get(k) != v:
            die("variant '%s' changed or removed: start a new study" % k)
    for k, v in old.TESTS.items():
        if repr(M.TESTS.get(k)) != repr(v):
            die("test '%s' changed or removed: start a new study" % k)
    for k in ("MATERIAL", "BUILD_GROUPS"):
        if getattr(M, k) != getattr(old, k):
            die("%s changed: start a new study" % k)
    shutil.copy(os.path.join(HARNESS, "manifest.py"), S.p("manifest_frozen.py"))
    print("refrozen: new tests/variants are now part of", S.path)


def cmd_pin(args):
    M = live_manifest()
    full = full_commit(M, args.commit)
    if subprocess.run(["git", "-C", M.REPO, "merge-base", "--is-ancestor", full, M.BRANCH]).returncode:
        die("%s is not on branch %s" % (args.commit, M.BRANCH))
    path = os.path.join(HARNESS, "manifest.py")
    src = open(path).read()
    import re
    new, n = re.subn(r'^COMMIT = "[0-9a-f]+"', 'COMMIT = "%s"' % full[:7], src, flags=re.M)
    if n != 1:
        die("could not find the COMMIT line in manifest.py")
    open(path, "w").write(new)
    print("manifest COMMIT -> %s (%s)" % (full[:7], sh(["git", "-C", M.REPO, "log", "-1", "--format=%s", full])))
    print("next: commit manifest.py, then ./pt new")


# ------------------------------------------------------------------------- builds
def flags_for(S, variant, group):
    return {**S.M.VARIANTS[variant].get("flags", {}), **S.M.BUILD_GROUPS[group]}


def needed_builds(S, tests=None):
    """[(variant, group, test)] -- one entry per distinct build dir."""
    out, seen = [], set()
    for t in S.tests(tests):
        T = S.M.TESTS[t]
        for v in T["variants"]:
            tag = S.build_tag(v, T["group"], t)
            if tag not in seen:
                seen.add(tag)
                out.append((v, T["group"], t))
    return out


def build_status(S, v, g, t):
    if S.code(v) == "gizmo":
        from . import gizmo
        return gizmo.build_ok(S, gizmo.geometry(S, t)[0])
    return snapshot_ok(S, v, g)


def set_flag(ph, name, val):
    lines = open(ph).read().split("\n")
    hits = [i for i, l in enumerate(lines) if l.split()[:2] == ["#define", name]]
    if len(hits) != 1:
        die("flag %s: %d '#define' lines in %s (need exactly 1)" % (name, len(hits), ph))
    lines[hits[0]] = "#define %s %s" % (name, val)
    open(ph, "w").write("\n".join(lines))


def snapshot_ok(S, variant, group):
    """(ok, reason) -- build exists, binary exists, SNAPSHOT matches commit + flags + md5."""
    B = S.p(S.build_tag(variant, group))
    snap = os.path.join(B, "SNAPSHOT.txt")
    if not os.path.exists(snap):
        return False, "missing"
    info = json.loads(open(snap).read().split("\n#JSON ", 1)[1].split("\n")[0])
    commit = S.meta["base_commit"] if S.M.VARIANTS[variant].get("commit") else S.meta["commit"]
    if info["commit"] != commit:
        return False, "commit %s != %s" % (info["commit"][:7], commit[:7])
    if info["flags"] != {k: str(v) for k, v in flags_for(S, variant, group).items()}:
        return False, "flags differ from manifest"
    have = sh("cat src/* include/* | md5sum | cut -c1-12", cwd=B)
    if have != info["md5"]:
        return False, "sources modified after build (md5)"
    if not os.path.exists(os.path.join(B, "bin", "mlh")):
        return False, "no binary (make failed?)"
    return True, "ok"


def build_newmfm(S, variant, group, jobs=4, force=False):
    tag = S.build_tag(variant, group)
    ok, why = snapshot_ok(S, variant, group)
    if ok and not force:
        print("  %-22s up to date" % tag)
        return True
    B = S.p(tag)
    commit = S.meta["base_commit"] if S.M.VARIANTS[variant].get("commit") else S.meta["commit"]
    flags = flags_for(S, variant, group)
    shutil.rmtree(B, ignore_errors=True)
    os.makedirs(B)
    sh("git -C %s archive %s demonstrator/Makefile demonstrator/include demonstrator/src demonstrator/tests"
       " | tar -x --strip-components=1 -C %s" % (S.M.REPO, commit, B))
    os.symlink(S.M.HOLIB, os.path.join(B, "holib"))
    ph = os.path.join(B, "include", "parameter.h")
    default = open(ph).read()
    for k, v in flags.items():
        set_flag(ph, k, v)
    md5 = sh("cat src/* include/* | md5sum | cut -c1-12", cwd=B)
    diff = subprocess.run(["diff", "-", ph], input=default, text=True, stdout=subprocess.PIPE).stdout
    info = dict(commit=commit, variant=variant, group=group, md5=md5,
                flags={k: str(v) for k, v in flags.items()})
    with open(os.path.join(B, "SNAPSHOT.txt"), "w") as f:
        f.write("%s  %s  flags: %s\n" % (tag, datetime.datetime.now().strftime("%F %T"),
                                          " ".join("%s:%s" % kv for kv in flags.items())))
        f.write(commit + "\n(git archive of the commit -- no working tree involved)\n")
        f.write("md5(src+include): %s\n\n" % md5)
        f.write("parameter.h diff vs the commit's default parameter.h:\n%s\n" % diff)
        f.write("#JSON %s\n" % json.dumps(info, sort_keys=True))
    r = subprocess.run("make -j%d > make.log 2>&1" % jobs, shell=True, cwd=B)
    if r.returncode != 0 or not os.path.exists(os.path.join(B, "bin", "mlh")):
        print("  %-22s FAILED (see %s/make.log):" % (tag, tag))
        print("    " + sh("grep -m5 -iE 'error' make.log", cwd=B, check=False).replace("\n", "\n    "))
        return False
    print("  %-22s built  (%s, %d flags)" % (tag, commit[:7], len(flags)))
    return True


def cmd_build(args):
    S = find_study(args.study)
    from . import gizmo
    todo = needed_builds(S, args.tests)
    if args.variants:
        todo = [b for b in todo if b[0] in args.variants.split(",")]
    print("builds for", S.path)
    bad = 0
    for v, g, t in todo:
        if S.code(v) == "gizmo":
            key, geo = gizmo.geometry(S, t)
            bad += not gizmo.build(S, key, geo, force=args.force)
        else:
            bad += not build_newmfm(S, v, g, force=args.force)
    if bad:
        die("%d build(s) failed" % bad)


# -------------------------------------------------------------------------- setup
def snap_t_end(t_end, dt):
    """(t_end', dt', n) with t_end' == n*dt' and (int)(t_end'/dt') == n in double arithmetic.
    MeshlessScheme.cpp builds its dump times as i*timeStep*h5DumpInterval with
    (int)(timeEnd/timeStep)/h5DumpInterval + 1 entries and stops at t >= timeEnd; a t_end that is
    not an exact multiple reads past the dumpTimes array."""
    n = max(1, int(round(t_end / dt)))
    for _ in range(64):
        te = n * dt * 1
        if int(te / dt) == n:
            return te, dt, n
        dt = math.nextafter(dt, math.inf)
    die("cannot make t_end = %r a multiple of dump_dt = %r" % (t_end, dt))


def material_block(mat):
    keys = ["rho0", "A", "B", "u0", "a", "b", "alpha", "beta", "u_iv", "u_cv", "mu", "Y0"]
    body = "".join("        %-6s %r\n" % (k, float(mat[k])) for k in keys)
    return "materials {\n    mat0 {\n        matId 0\n" + body + "    }\n}"


def write_newmfm_config(S, test, variant, rung, rel, ic_rel, ks):
    T = S.M.TESTS[test]
    mat = {**S.M.MATERIAL, **T.get("material", {})}
    te, dt, n = snap_t_end(T["t_end"], T["dump_dt"])
    x0, y0, x1, y1 = T["box"]
    cfg = ["; plasticity test plan %s / %s  (%s)" % (test, os.path.basename(rel), S.build_tag(variant, T["group"])),
           "initFile %s" % ic_rel, "outDir runs/%s/output" % rel,
           "timeStep %r" % dt, "timeEnd %r" % te, "h5DumpInterval 1",
           "periodicBoxLimits {\n    lowerX %r\n    lowerY %r\n    lowerZ 0.0\n    upperX %r\n    upperY %r\n"
           "    upperZ 0.0\n}" % (x0, y0, x1, y1),
           "kernelSize %r" % ks,
           "sml {\n    NNNTarget %d\n    tol       1e-3\n    maxIter   200\n    maxFactor 1.2\n    hMinFactor 0.1\n"
           "    hMaxFactor 3.0\n    warnFraction  0.05\n    panicFraction 0.25\n}" % S.M.NNN,
           "hydro_gamma 1.6666666666666667", material_block(mat)]
    open(S.p("runs", rel, "config.info"), "w").write("\n".join(cfg) + "\n")
    return te


def cmd_setup(args):
    S = find_study(args.study)
    from . import gen, gizmo
    variants = args.variants.split(",") if args.variants else None
    rungs = args.rungs.split(",") if args.rungs else None
    for test in S.tests(args.tests):
        T = S.M.TESTS[test]
        new = []
        for v, r, rel in S.runs(test, variants, rungs):
            d = S.p("runs", rel)
            if os.path.exists(os.path.join(d, "exit_code")):
                continue                                   # finished runs are never touched
            ic_rel, ks, meta = gen.ensure_ic(S, test, r)
            os.makedirs(os.path.join(d, "output"), exist_ok=True)
            tag = S.build_tag(v, T["group"], test)
            open(os.path.join(d, "BUILD"), "w").write(tag + "\n")
            if S.code(v) == "gizmo":
                os.makedirs(os.path.join(d, "gizmo_out"), exist_ok=True)
                gizmo.setup_run(S, test, r, rel, ic_rel, ks, meta)
            else:
                write_newmfm_config(S, test, v, r, rel, ic_rel, ks)
                ph = S.p(tag, "include", "parameter.h")
                if not os.path.exists(ph):
                    die("build %s missing -- run ./pt build first" % tag)
                shutil.copy(ph, d)
            new.append((meta.get("N", 0), rel, v))
        # newMFM job list for binac-cpu-runs (cost-ordered, largest first)
        jl = [rel for n, rel, v in sorted(new, key=lambda x: -x[0]) if S.code(v) != "gizmo"]
        if jl:
            with open(S.p("jobs_%s.txt" % test), "a") as f:
                f.write("".join(j + "\n" for j in jl))
        print("  %-12s %d new run dir(s)" % (test, len(new)))


# ---------------------------------------------------------------------------- run
def run_state(S, rel):
    d = S.p("runs", rel)
    if not os.path.exists(d):
        return "absent"
    if os.path.exists(os.path.join(d, "exit_code")):
        return "done" if open(os.path.join(d, "exit_code")).read().strip() == "0" else "FAILED"
    if os.path.exists(os.path.join(d, "pid")):
        try:
            os.kill(int(open(os.path.join(d, "pid")).read()), 0)
            return "running"
        except (OSError, ValueError):
            return "DIED"
    return "pending"


def launch(S, rel):
    """Start one run in its own process group; returns (Popen, t0)."""
    d = S.p("runs", rel)
    tag = open(os.path.join(d, "BUILD")).read().strip()
    log = open(os.path.join(d, "run.log"), "w")
    if tag.startswith("gizmo_"):
        from . import gizmo
        cmd = gizmo.run_command(S, rel)
        cwd = d
    else:
        cmd = ["./%s/bin/mlh" % tag, "-c", "runs/%s/config.info" % rel]
        cwd = S.path
    p = subprocess.Popen(cmd, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
    open(os.path.join(d, "pid"), "w").write(str(p.pid))
    return p, time.time()


def finish(S, rel, rc, t0):
    d = S.p("runs", rel)
    tag = open(os.path.join(d, "BUILD")).read().strip()
    if tag.startswith("gizmo_") and rc == 0:
        from . import gizmo
        rc = gizmo.convert_output(S, rel)
    open(os.path.join(d, "exit_code"), "w").write("%d\n" % rc)
    try:
        os.remove(os.path.join(d, "pid"))
    except OSError:
        pass
    log = os.path.join(d, "run.log")
    txt = open(log, errors="replace").read() if os.path.getsize(log) < 2e8 else ""
    nerr = sum(1 for l in txt.split("\n") if "ERROR" in l)
    nneg = txt.count("Negative internal")
    line = "%s exit=%d snaps=%d wall=%ds errors=%d negE_lines=%d local\n" % (
        rel, rc, len(glob.glob(os.path.join(d, "output", "*.h5"))), time.time() - t0, nerr, nneg)
    open(S.p("summary.txt"), "a").write(line)
    print("  finished", line.strip())


def cmd_run(args):
    S = find_study(args.study)
    if args.detach:
        log = S.p("queue_%s.log" % datetime.datetime.now().strftime("%H%M%S"))
        argv = [a for a in sys.argv if a != "--detach"]
        subprocess.Popen([sys.executable] + argv, stdout=open(log, "w"), stderr=subprocess.STDOUT,
                         stdin=subprocess.DEVNULL, start_new_session=True, cwd=os.getcwd())
        print("queue detached, log:", log)
        return
    variants = args.variants.split(",") if args.variants else None
    rungs = args.rungs.split(",") if args.rungs else None
    todo = []
    for test in S.tests(args.tests):
        for v, r, rel in S.runs(test, variants, rungs):
            st = run_state(S, rel)
            if st == "pending":
                ok, why = build_status(S, v, S.M.TESTS[test]["group"], test)
                if not ok:
                    print("  skip %s: build %s" % (rel, why))
                if ok:
                    todo.append(rel)
            elif st == "absent":
                print("  skip %s: not set up (./pt setup %s)" % (rel, test))
            elif st == "DIED":
                print("  skip %s: a previous run died without exit_code -- inspect, then ./pt reset %s"
                      % (rel, rel))
    if not todo:
        print("nothing to run")
        return
    maxj = args.jobs
    timeout = args.timeout
    print("running %d run(s), %d at a time (GIZMO runs count as 2 slots)" % (len(todo), maxj))
    active = {}
    slots = lambda: sum(2 if open(S.p("runs", r, "BUILD")).read().startswith("gizmo_") else 1 for r in active)
    while todo or active:
        while todo and slots() < maxj:
            rel = todo.pop(0)
            active[rel] = launch(S, rel)
            print("  started", rel)
        time.sleep(3)
        for rel, (p, t0) in list(active.items()):
            rc = p.poll()
            logsz = os.path.getsize(S.p("runs", rel, "run.log"))
            if rc is None and (time.time() - t0 > timeout or logsz > 2e9):
                os.killpg(p.pid, signal.SIGTERM)
                p.wait()
                rc = 124 if time.time() - t0 > timeout else 125
                print("  KILLED %s (%s)" % (rel, "timeout" if rc == 124 else "log > 2 GB"))
            if rc is not None:
                finish(S, rel, rc, t0)
                del active[rel]


def cmd_reset(args):
    S = find_study(args.study)
    d = S.p("runs", args.run)
    if not os.path.isdir(d):
        die("no run dir " + d)
    for f in ("exit_code", "pid", "run.log"):
        try:
            os.remove(os.path.join(d, f))
        except OSError:
            pass
    shutil.rmtree(os.path.join(d, "output"), ignore_errors=True)
    os.makedirs(os.path.join(d, "output"))
    print("reset", args.run, "(output deleted; summary.txt keeps the old line)")


# ------------------------------------------------------------------------- status
def cmd_status(args):
    S = find_study(args.study)
    M = S.M
    live = live_manifest()
    print("study    ", S.path)
    print("commit   ", S.meta["commit"][:7], "on", S.meta["branch"], "  harness", S.meta["harness"])
    head = sh(["git", "-C", M.REPO, "rev-parse", "--short", M.BRANCH], check=False)
    if head != S.meta["commit"][:7]:
        print("NOTE      branch %s is now at %s (study pins %s)" % (M.BRANCH, head, S.meta["commit"][:7]))
    if live.COMMIT != M.COMMIT:
        print("NOTE      live manifest pins %s -> this study is not the current one" % live.COMMIT)
    # errors first (elastic plan rule 6)
    errs = []
    for test in S.tests(args.tests):
        for v, r, rel in S.runs(test):
            st = run_state(S, rel)
            if st in ("FAILED", "DIED"):
                errs.append("%-40s %s (exit %s)" % (rel, st, open(S.p("runs", rel, "exit_code")).read().strip()
                                                    if st == "FAILED" else "-"))
    if os.path.exists(S.p("summary.txt")):
        for l in open(S.p("summary.txt")):
            f = dict(kv.split("=", 1) for kv in l.split()[1:] if "=" in kv)
            if f.get("exit") == "0" and (f.get("errors", "0") != "0" or f.get("negE_lines", "0") != "0"):
                errs.append("%-40s exit 0 but errors=%s negE_lines=%s" % (l.split()[0], f.get("errors"),
                                                                          f.get("negE_lines")))
    print("\nERRORS   " + ("none" if not errs else "\n  " + "\n  ".join(errs)))
    # builds
    print("\nBUILDS")
    from . import gizmo
    for v, g, t in needed_builds(S, args.tests):
        ok, why = build_status(S, v, g, t)
        print("  %-24s %s" % (S.build_tag(v, g, t), "ok" if ok else why))
    # runs + checks per test, in gate order
    print("\nTESTS (gate order)  runs: done/total  [running, pending, FAILED]   check")
    for tier in M.ORDER:
        if tier == "P0a":
            print("  %-4s %-12s %-28s %s" % (tier, "(unit)", "", check_verdict(S, "P0a")))
            continue
        for test in [t for t in S.tests(args.tests) if M.TESTS[t]["tier"] == tier]:
            sts = [run_state(S, rel) for v, r, rel in S.runs(test)]
            c = {k: sts.count(k) for k in set(sts)}
            print("  %-4s %-12s %2d/%-2d [%d run, %d pend, %d FAIL]   %s" % (
                tier, test, c.get("done", 0), len(sts), c.get("running", 0), c.get("pending", 0) + c.get("absent", 0),
                c.get("FAILED", 0) + c.get("DIED", 0), check_verdict(S, test)))


def check_verdict(S, test):
    j = S.p("checks", test + ".json")
    if not os.path.exists(j):
        return "-"
    v = json.load(open(j))
    return "%s (%s)" % (v["verdict"], v.get("when", ""))


def cmd_next(args):
    S = find_study(args.study)
    for tier in S.M.ORDER:
        tests = ["P0a"] if tier == "P0a" else [t for t in S.M.TESTS if S.M.TESTS[t]["tier"] == tier]
        verdicts = [check_verdict(S, t).split(" ")[0] for t in tests]
        if all(v == "PASS" for v in verdicts):
            continue
        gate = " (GATE)" if tier in S.M.GATES else ""
        print("next open tier: %s%s -- tests %s, verdicts %s" % (tier, gate, ", ".join(tests), verdicts))
        if tier == "P0a":
            print("do: ./pt unit")
        else:
            print("do: ./pt setup %s && ./pt run %s --detach   (when finished) ./pt check %s" % (tier, tier, tier))
        return
    print("all tiers PASS")


# ------------------------------------------------------------------------- verify
def cmd_verify(args):
    S = find_study(args.study)
    bad = 0
    from . import gizmo
    for v, g, t in needed_builds(S):
        tag = S.build_tag(v, g, t)
        if not os.path.exists(S.p(tag)):
            continue
        ok, why = build_status(S, v, g, t)
        if not ok:
            bad += 1
            print("BUILD %-22s %s" % (tag, why))
    for test in S.tests():
        for v, r, rel in S.runs(test):
            d = S.p("runs", rel)
            if not os.path.exists(os.path.join(d, "BUILD")):
                continue
            tag = open(os.path.join(d, "BUILD")).read().strip()
            if S.code(v) != "gizmo":
                a, b = os.path.join(d, "parameter.h"), S.p(tag, "include", "parameter.h")
                if not (os.path.exists(a) and os.path.exists(b) and open(a).read() == open(b).read()):
                    bad += 1
                    print("RUN   %-40s parameter.h differs from %s" % (rel, tag))
    print("verify: %s" % ("OK -- every build matches the manifest, every run its build" if not bad
                          else "%d problem(s)" % bad))
    if bad:
        sys.exit(1)


# ---------------------------------------------------------------------------- unit
def cmd_unit(args):
    """P0a: extract the pinned commit, run make test-plasticity, record the verdict."""
    S = find_study(args.study)
    U = S.p("unit_P0a")
    shutil.rmtree(U, ignore_errors=True)
    os.makedirs(U)
    sh("git -C %s archive %s demonstrator/Makefile demonstrator/include demonstrator/src demonstrator/tests"
       " | tar -x --strip-components=1 -C %s" % (S.M.REPO, S.meta["commit"], U))
    os.symlink(S.M.HOLIB, os.path.join(U, "holib"))
    r = subprocess.run("make test-plasticity", shell=True, cwd=U, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out = r.stdout
    os.makedirs(S.p("checks"), exist_ok=True)
    open(S.p("checks", "P0a.txt"), "w").write(out)
    verdict = "PASS" if r.returncode == 0 and "PASSED" in out else "FAIL"
    json.dump(dict(verdict=verdict, when=datetime.datetime.now().strftime("%m-%d %H:%M"),
                   commit=S.meta["commit"][:7]), open(S.p("checks", "P0a.json"), "w"))
    print("\n".join(l for l in out.split("\n") if l.startswith("  ") or "PASS" in l or "FAIL" in l))
    print("P0a:", verdict)
