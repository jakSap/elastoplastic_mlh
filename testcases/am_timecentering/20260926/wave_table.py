#!/usr/bin/env python3
"""wave_table.py: W1-W3 metrics (testplan_20260926/analysis/wave_analysis.py, unchanged) for runsw/<test>/<variant>_N<n>,
next to the testplan's own ref rows (standard scheme, night 2026-09-26). Orders p between consecutive n."""
import os, sys, glob, re, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); TP = os.path.join(HERE, "..", "..", "testplan_20260926")
sys.path.insert(0, os.path.join(TP, "analysis")); import wave_analysis as wa
res = []
for d in sorted(glob.glob(os.path.join(HERE, "runsw", "*", "*"))):
    if not os.path.exists(os.path.join(d, "exit_code")) or not glob.glob(os.path.join(d, "output", "*.h5")): continue
    test = d.split(os.sep)[-2]; n = int(re.search(r"_N(\d+)", d).group(1))
    r = wa.analyse("W1_sound" if test == "W1_ideal" else test, d, n)
    if r: r["exit"] = open(os.path.join(d, "exit_code")).read().strip(); r["grp"] = test + "/" + re.sub(r"_N\d+$", "", os.path.basename(d)); res.append(r)
for r in json.load(open(os.path.join(TP, "analysis", "wave_results.json"))):
    m = re.match(r"runs/(W[123]_\w+?)/ref(_tillotson|_ideal)?_N(\d+)$", r["run"])
    if m and (m.group(1) != "W1_sound" or m.group(2)):
        r["grp"] = ("W1_ideal" if m.group(2) == "_ideal" else m.group(1)) + "/REF(testplan)"; res.append(r)
keys = ["L1rho", "L1v", "damp", "dc", "L1S", "defect", "contam_vx", "dE"]
groups = {}
for r in res: groups.setdefault(r["grp"], []).append(r)
for g in sorted(groups):
    rs = sorted(groups[g], key=lambda r: r["n"]); ks = [k for k in keys if k in rs[0]]
    print("## %s" % g); print("%5s %5s " % ("n", "exit") + " ".join("%10s" % k for k in ks))
    for i, r in enumerate(rs):
        print("%5d %5s " % (r["n"], r["exit"]) + " ".join("%+10.3e" % r[k] for k in ks))
        if i: print("%11s " % "p" + " ".join("%10.2f" % wa.order(rs[i-1][k], r[k]) if k != "dE" else "%10s" % "" for k in ks))
    print()
