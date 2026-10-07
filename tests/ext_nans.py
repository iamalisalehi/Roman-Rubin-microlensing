"""QC of the simulator's extinction tables, files/ext/ext_tables.dat (Deviation 70).

Run from tests/ (`python ext_nans.py`) or the repo root (`python tests/ext_nans.py`). Checks what
src/galaxy/extinction.cpp readExtinction() also refuses -- a non-finite value, a decreasing profile, a row count or
a distance grid that disagrees with the header -- and prints a summary of A_V at 8 kpc. The tables
before Deviation 70 (one Bayestar/DECaPS file per pointing, 78 of them all-NaN) are archived in
files/ext_bayestar_v1/ and then deleted (2026-10-01); if regenerated, pass that directory to check
them the old way.
"""
import glob
import os
import sys

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
arg = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "files", "ext", "ext_tables.dat")

if os.path.isdir(arg):                                   # the archived one-file-per-table layout
    total, partial = [], []
    for fname in glob.glob(os.path.join(arg, "*.txt")):
        ext = np.loadtxt(fname)[:, 3]
        n = np.isnan(ext).sum()
        if n == len(ext):
            total.append(fname)
        elif n:
            partial.append((fname, n, len(ext)))
    print(f"{len(total)} sightlines entirely NaN; {len(partial)} with isolated NaN gaps")
    sys.exit(1 if total or partial else 0)

head = {}
with open(arg) as f:
    for line in f:
        if not line.startswith("#"):
            break
        p = line[1:].split()
        if p and p[0] in ("n_tables", "n_dist", "k"):
            head[p[0]] = float(p[1])
        elif p and p[0] == "dist":
            head["dist"] = np.array(p[1:], float)
t = np.loadtxt(arg, comments="#", dtype=np.float32)
lb, av = t[:, :2], t[:, 2:]
errs = []
if t.shape[0] != head["n_tables"]:
    errs.append(f"{t.shape[0]} rows, header says {int(head['n_tables'])}")
if av.shape[1] != head["n_dist"] or head["dist"].size != head["n_dist"]:
    errs.append("distance grid disagrees with header")
if not np.isfinite(av).all():
    errs.append(f"{(~np.isfinite(av)).any(axis=1).sum()} rows with non-finite values")
if (np.diff(av, axis=1) < -1e-4).any():
    errs.append(f"{(np.diff(av, axis=1) < -1e-4).any(axis=1).sum()} rows decreasing with distance")
i8 = int(np.argmin(np.abs(head["dist"] - 8.0)))
b = lb[:, 1]
print(f"{t.shape[0]} tables x {av.shape[1]} distances, k = {head['k']}")
for lo, hi in ((0, 0.5), (0.5, 1), (1, 1.5), (1.5, 2.5), (2.5, 6)):
    m = (np.abs(b) >= lo) & (np.abs(b) < hi)
    print(f"  |b| {lo:.1f}-{hi:.1f}: {m.sum():5d} tables, median A_V(8 kpc) {np.median(av[m, i8]):6.2f}")
print("ERRORS: " + "; ".join(errs) if errs else "all checks held")
sys.exit(1 if errs else 0)
