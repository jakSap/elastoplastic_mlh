#!/usr/bin/env python3
"""Fixed-dt configs + job queue for the time-order sweep (2026-09-24).

Same dt path as the 2026-09-23 sweep: ADAPTIVE_TIMESTEP 0, dt = t_target/M
exactly.  Two differences, both cost, neither method:

  * timeEnd = 1.15 * t_target instead of 2.2.  The run only has to REACH the
    target -- the loop breaks at the top of the first iteration with
    t >= timeEnd, after the dump block, so iteration M (t = t_target) is passed
    and dumped well before the break.  Halves the CPU per run.
  * the dump interval is unchanged from the old sweep, because it is the part
    that is easy to get wrong: `step` is incremented once per iteration AND once
    more inside the dump block, so with interval I the j-th dump sits at
    t_j = j (I-1) dt.  Landing a dump exactly on t_target = M dt therefore needs
    (I-1) | M -- hence I = d+1 with d the largest divisor of M below M/20.
"""
import os, re, sys

ROOT = os.path.expanduser(os.environ.get("ROOT", "~/Documents/tord_sweep"))
SW = os.path.join(ROOT, "testcases", "soundwave_2d")
END_FAC = 1.15

LAD_N64 = [192, 256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096, 6144, 8192, 12288, 16384]
LAD_N128 = [384, 768, 1536, 3072, 6144, 12288, 24576]

# tier, variant tag, base config, n, ladder, cost weight per step
CASES = [
    # tier 0 -- is the off-centred kick the remaining O(dt)?
    (0, "w2dkd",      "config_ideal_N64.info",  64, LAD_N64, 1.6),
    (0, "w2kdk",      "config_ideal_N64.info",  64, LAD_N64, 1.0),
    (0, "w2dkdpred",  "config_ideal_N64.info",  64, LAD_N64, 1.6),
    # tier 1 -- the two cheaper answers to the same question, plus best-of
    (1, "w2nopred",   "config_ideal_N64.info",  64, LAD_N64, 1.0),
    (1, "w2frame",    "config_ideal_N64.info",  64, LAD_N64, 1.0),
    (1, "w2dkdfull",  "config_ideal_N64.info",  64, LAD_N64, 1.6),
    (1, "w2dkdLimN",  "config_ideal_N64.info",  64, LAD_N64, 1.6),
    # tier 2 -- n = 128 for the three that decide it
    (2, "w2dkd",      "config_ideal_N128.info", 128, LAD_N128, 6.4),
    (2, "w2dkdpred",  "config_ideal_N128.info", 128, LAD_N128, 6.4),
    (2, "w2kdk",      "config_ideal_N128.info", 128, LAD_N128, 4.0),
]



def rewrite(base, out_cfg, dt, t_end, interval, outdir, M):
    src = open(base).read()

    def sub(key, val, s):
        if re.search(r"(?m)^%s\s" % key, s):
            return re.sub(r"(?m)^%s\s+\S+.*$" % key, "%s %s" % (key, val), s)
        return s + "\n%s %s\n" % (key, val)

    src = sub("timeStep", repr(dt), src)
    src = sub("timeEnd", repr(t_end), src)
    src = sub("h5DumpInterval", str(interval), src)
    src = sub("outDir", outdir, src)
    src = ("; time-order sweep 2026-09-24: ADAPTIVE_TIMESTEP 0, dt = t_target/%d,\n"
           "; dump interval %d -> dumps at t = j*%d*dt, so j = %d lands on t_target.\n"
           % (M, interval, interval - 1, M // (interval - 1))) + src
    open(out_cfg, "w").write(src)


jobs = []
for tier, tag, basecfg, n, ladder, w in CASES:
    base = os.path.join(SW, basecfg)
    build = "build_t_%s" % tag
    if not os.path.exists(base):
        print("MISSING base config:", base, file=sys.stderr); sys.exit(1)
    if not os.path.exists(os.path.join(SW, build, "bin", "mlh")):
        print("MISSING build:", build, file=sys.stderr); sys.exit(1)
    for M in ladder:
        name = "%s_N%d_M%d" % (tag, n, M)
        cfg = os.path.join(SW, "config_tord_%s.info" % name)
        out = "output_tord_%s" % name
        d = max(k for k in range(1, M // 20 + 1) if M % k == 0)
        rewrite(base, cfg, 1.0 / M, END_FAC * 1.0, d + 1, out, M)
        os.makedirs(os.path.join(SW, out), exist_ok=True)
        jobs.append((tier, M * w, SW, build, os.path.basename(cfg), out, name))

# longest-processing-time first inside each tier: better packing across 6 workers
jobs.sort(key=lambda j: (j[0], -j[1]))
with open(os.path.join(ROOT, "jobs.txt"), "w") as f:
    for tier, cost, dd, build, cfg, out, name in jobs:
        f.write("\t".join([dd, build, cfg, out, name]) + "\n")
print("%d jobs, %d builds, sum of steps = %d"
      % (len(jobs), len({j[3] for j in jobs}), sum(int(j[1]) for j in jobs)))
