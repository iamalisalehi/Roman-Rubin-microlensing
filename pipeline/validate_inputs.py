#!/usr/bin/env python3
"""Cheap sanity checks on the two big downloaded inputs, run before the pipeline spends minutes on them.

    python3 pipeline/validate_inputs.py besancon FILE     # header names the builder needs; row widths
    python3 pipeline/validate_inputs.py opsim DB          # `observations` table, its columns, time span

Exit 0 = usable (warnings may be printed), 1 = not usable (reason on stderr). Reads only the first few
thousand lines of the Besancon file and one aggregate query of the OpSim database.
"""
import ast
import os
import sqlite3
import sys
import urllib.parse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# The columns CMD/BolometricCorrection.py (DEFAULT_COLUMNS) and CMD/lens_ml_table.py (COLS) read.
BESANCON_NEEDED = ["Teff", "logg", "Pop", "Age", "Mass", "Mbol", "[M/H]", "[a/Fe]", "CL", "Typ"]
SAMPLE_ROWS = 5000


def fail(msg):
    print(f"validate_inputs: {msg}", file=sys.stderr)
    sys.exit(1)


def besancon(path):
    if not os.path.isfile(path):
        fail(f"no such file: {path}")
    with open(path, errors="replace") as fh:
        first = fh.readline()
        if not first.startswith("#"):
            fail(f"{path}: first line is not a '#' column header (is this a Besancon catalogue?)")
        names = first.lstrip("#").split()
        missing = [c for c in BESANCON_NEEDED if c not in names]
        if missing:
            fail(f"{path}: header lacks column(s) {missing}; has {len(names)}: {' '.join(names)}")
        good = bad = 0
        for _ in range(SAMPLE_ROWS):
            line = fh.readline()
            if not line:
                break
            if line.startswith("#") or not line.strip():
                continue
            if len(line.split()) == len(names):
                good += 1
            else:
                bad += 1
    if good == 0 or bad > 0.01 * (good + bad):
        fail(f"{path}: {bad} of the first {good + bad} data rows do not have {len(names)} fields")
    note = f" ({bad} damaged rows in the first {good + bad}; the builder skips them)" if bad else ""
    print(f"besancon: header ok, {len(names)} columns, all {len(BESANCON_NEEDED)} needed ones present; "
          f"{os.path.getsize(path) / 1e9:.2f} GB{note}")
    if len(names) != 38:
        print(f"  warning: bos10 has 38 columns, this file has {len(names)}")


def literal_from_source(path, name):
    """The literal value assigned to a top-level `name = ...` in a Python source file."""
    for node in ast.parse(open(path).read()).body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == name for t in node.targets):
            return ast.literal_eval(node.value)
    raise KeyError(name)


def opsim(path):
    if not os.path.isfile(path):
        fail(f"no such file: {path}")
    script = os.path.join(ROOT, "Baseline", "readbaselineBulge.py")
    nam0, idx = literal_from_source(script, "nam0"), literal_from_source(script, "idx")
    need = [nam0[i] for i in idx]
    t0 = literal_from_source(script, "TIME0_MJD")
    uri = "file:" + urllib.parse.quote(os.path.abspath(path)) + "?mode=ro"
    try:
        con = sqlite3.connect(uri, uri=True)
        cols = [r[1] for r in con.execute("PRAGMA table_info(observations)")]
    except sqlite3.DatabaseError as e:
        fail(f"{path}: not a usable sqlite database ({e})")
    if not cols:
        fail(f"{path}: no table `observations` (is this an OpSim baseline database?)")
    missing = [c for c in need if c not in cols]
    if missing:
        fail(f"{path}: `observations` lacks column(s) {missing} that readbaselineBulge.py selects")
    n, lo, hi = con.execute("SELECT count(*), min(observationStartMJD), max(observationStartMJD) "
                            "FROM observations").fetchone()
    con.close()
    sys.path.insert(0, os.path.join(ROOT, "analysis"))
    try:
        from cparams import P
        tobs = float(P.Tobs)
    except Exception as e:                                     # numpy missing, header moved, ...
        tobs = None
        print(f"  (could not read Tobs from config/parameters.h: {e})")
    print(f"opsim: {n:,} visits, MJD {lo:.3f} .. {hi:.3f}; the simulation's day 0 is MJD {t0:.3f}")
    d0 = lo - t0
    if d0 > 1.0:
        print(f"  WARNING: this baseline starts {d0:.1f} d AFTER day 0, so Rubin has no visits at all in the "
              f"first {d0:.1f} d of the simulated clock (Roman's mission is placed on that clock).")
    elif d0 < -1.0:
        print(f"  note: this baseline starts {-d0:.1f} d BEFORE day 0; visits before it are dropped.")
    if tobs is not None:
        end = t0 + tobs
        if hi < end - 1.0:
            print(f"  note: this baseline ends {end - hi:.1f} d BEFORE the clock does (day {tobs:.1f} = "
                  f"MJD {end:.3f}): Rubin has no visits in the last {end - hi:.1f} d of the window "
                  "(a 10-yr baseline that starts before day 0 always does this).")
        elif hi > end + 1.0:
            print(f"  note: visits after day {tobs:.1f} (MJD {end:.3f}) are dropped.")


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in ("besancon", "opsim"):
        sys.exit(__doc__)
    {"besancon": besancon, "opsim": opsim}[sys.argv[1]](sys.argv[2])


if __name__ == "__main__":
    main()
