#!/usr/bin/env python3
"""Dynamic mto queue: keeps up to 4 jobs running, always starting the FIRST not-yet-started line
of the queue file, which is re-read every time a slot frees -- lines can be inserted or appended
while it runs. State (started lines) in runs/queue_state.txt; jobs via queue_job.sh.

    ./runqueue.py queue/2026-10-03.txt [-j 4]
Runs on battery too (user 10-03); waits only while Low Power Mode is on.
"""
import os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
Q = sys.argv[1] if len(sys.argv) > 1 else "queue/2026-10-03.txt"
J = int(sys.argv[sys.argv.index("-j") + 1]) if "-j" in sys.argv else 4
STATE, LOG = "runs/queue_state.txt", "runs/queue_jobs.log"

def log(msg):
    with open(LOG, "a") as f: f.write(time.strftime("%H:%M:%S ") + msg + "\n")

def power_ok():
    """Normal power mode (Low Power Mode off); battery is fine."""
    p = subprocess.run(["pmset", "-g"], capture_output=True, text=True).stdout
    return any(l.split()[:2] == ["lowpowermode", "0"] for l in p.splitlines())

def entries():
    return [" ".join(l.split()) for l in open(Q) if l.strip() and not l.lstrip().startswith("#")]

def busy(running):
    """Slots in use: every own job holds one slot from its start (also while it builds or makes its
    IC), plus solvers that are not ours (orphans of an earlier runner, manual runs)."""
    kids = {}
    for l in subprocess.run(["ps", "-axo", "pid=,ppid="], capture_output=True, text=True).stdout.splitlines():
        pid, ppid = map(int, l.split()); kids.setdefault(ppid, []).append(pid)
    mine, stack = set(), [p.pid for p in running.values()]
    while stack:
        q = stack.pop()
        if q not in mine: mine.add(q); stack += kids.get(q, [])
    foreign = 0
    for pat in ("bin/mlh -c config.info", "GIZMO rod.params"):
        foreign += sum(int(q) not in mine for q in
                       subprocess.run(["pgrep", "-f", pat], capture_output=True, text=True).stdout.split())
    return len(running) + foreign

def started():
    return set(open(STATE).read().splitlines()) if os.path.exists(STATE) else set()

while not power_ok():
    log("waiting for Low Power Mode off"); time.sleep(120)
log("dynamic queue {} started (-j {})".format(Q, J))
running = {}
while True:
    for e, p in list(running.items()):
        if p.poll() is not None: del running[e]
    done = started()
    todo = [e for e in entries() if e not in done]
    while todo and busy(running) < J:
        e = todo.pop(0)
        with open(STATE, "a") as f: f.write(e + "\n")
        running[e] = subprocess.Popen(["./queue_job.sh"] + e.split())
        time.sleep(2)
    if not todo and not running: break
    time.sleep(10)
log("dynamic queue {} finished".format(Q))
