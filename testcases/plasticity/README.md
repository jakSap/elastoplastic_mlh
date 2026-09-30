# Plasticity test harness (`pt`)

Runs the plasticity test plan — `MA-Obsidian/2026-09-30-testplan-plastisch.md` (Appendix P to the
elastic test plan) — for the demonstrator (newMFM) against GIZMO's MFM. Written so that a model with
little context can operate it without ever touching flags or source trees by hand.

Everything lives on branch `plasticity-tests` (worktree `~/Documents/newMFM-plasticity`). Run all
commands from this directory: `cd ~/Documents/newMFM-plasticity/testcases/plasticity`.

## The three rules

1. **Never edit flags, `parameter.h`, build dirs or configs by hand.** The only file that defines
   what is tested is `manifest.py`. Builds are extracted from git (`git archive <COMMIT>`), never
   copied from a working tree, and every build records its commit, flag set and source md5.
2. **The code version is `COMMIT` in `manifest.py`.** Change it only with `./pt pin <hash>` (it
   refuses commits that are not on the branch), commit `manifest.py`, then `./pt new` — a new code
   version is a new study. An old study keeps its own frozen copy of the manifest and stays valid.
3. **Gate order is binding** (`./pt next` tells you where you are). Report errors first (exit ≠ 0,
   `ERROR` lines, negative internal energy) — `./pt status` lists them at the top.

## Typical session

```bash
./pt status                 # where are we? errors first, builds, runs, verdicts per test
./pt next                   # first tier that is not PASS + the exact commands for it
./pt setup P0b              # ICs + run dirs for a tier (or one test: P0b_psh)
./pt run P0b --detach       # background queue (log: <study>/queue_HHMMSS.log); -j 4 slots default
./pt status P0b             # progress
./pt check P0b              # verdicts -> <study>/checks/P0b_*.txt|json
./pt movie P0b_psh --field gp --variants pref --rungs 0.0125   # look at <study>/movies/*_last.png
```

Start of a code version: `./pt new`, `./pt unit` (P0a), `./pt build`.

## Commands

| Command | What it does |
|---|---|
| `new [--date]` | create `newMFM/testcases/plasticity_<date>_<commit7>/` for the manifest's COMMIT |
| `pin <hash>` | set COMMIT in `manifest.py` (on the branch only) |
| `refreeze` | add *new* tests/variants from `manifest.py` to the current study (refuses changes to existing ones) |
| `unit` | P0a: `make test-plasticity` on a git-archive of COMMIT |
| `build [tests] [--variants] [--force]` | every build the tests need; up-to-date builds are skipped |
| `setup [tests] [--variants] [--rungs]` | ICs (cached in `ics/`, listed in `ics/IC_SOURCES.txt`) and run dirs; finished runs are never touched |
| `run [tests] [-j N] [--timeout s] [--detach]` | run pending runs; GIZMO runs take 2 slots; a run is killed at the timeout (exit 124) or a 2 GB log (125) |
| `reset <test>/<run>` | delete one run's result so it runs again |
| `status [tests]` / `next` / `verify` | overview / what to do next / provenance check |
| `check [tests]` | run the checkers: PASS, FAIL, INCOMPLETE, LOOK (+ FINDING lines) |
| `accept <test> "<note>"` | LOOK → PASS after you inspected the movie |
| `movie [tests] --field gp|j2|sigma|rho|fab` | mp4 + last-frame PNG per run in `movies/` |

`[tests]` is a test name (`P2_plate`), a tier (`P2`) or a comma list; default all.

## What the names mean

- **Variants** (`manifest.py`, `VARIANTS`): `pref` reference, `pel` elastic baseline, `prr` exact
  radial return, `ptc0` tensile correction off, `pps` old plane-strain rate, `legacy`/`base` branch
  defaults vs. `3835f4a` (P1a only), `gz` GIZMO.
- **Build dirs**: `build_<variant>` or `build_<variant>_<group>` (`per` = periodic in y, `wave` =
  Tier-W flags), `gizmo_std` / `gizmo_per_<aspect>`.
- **Run dirs**: `runs/<test>/<variant>_d<spacing>` or `_n<particles per wavelength>`.

## When something goes wrong

- A build FAILED: read `<study>/build_<tag>/make.log`. Do not patch the build dir — fix the code on the
  branch, commit, `./pt pin`, `./pt new`.
- A run FAILED: `./pt status` shows it first. Read its `run.log` tail. Report it (elastic plan rule 6)
  before looking at any verdict. Rerun only after understanding it (`./pt reset`).
- `DIED` (no exit_code, process gone): the machine or the queue was interrupted. `./pt reset` it.
- `status` says the branch moved past the study's commit: that is information, not an error. Results
  stay tied to the study's commit. Pin the new commit only when told to.
- Big queues: the newMFM jobs of a test are listed in `<study>/jobs_<test>.txt` in the layout of the
  `binac-cpu-runs` skill (stage with `binac_stage.sh <study> jobs_<test>.txt`). GIZMO runs stay local.

## Files

`pt` (CLI) · `manifest.py` (the definition) · `ptlib/core.py` (study, builds, queue, status) ·
`ptlib/gen.py` (ICs) · `ptlib/gizmo.py` (GIZMO build, params, conversion) · `ptlib/check.py`
(checkers, material-point ODE, exact plate solution, 1D reference code) · `ptlib/movie.py`.
