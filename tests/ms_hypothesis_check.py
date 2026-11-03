"""Check CMD/components/ms_hypothesis.dat (lens_ml_table.py --hypothesis).

    python tests/ms_hypothesis_check.py PATH [LENS_ML]

LENS_ML defaults to CMD/components/lens_ml.dat next to this repository (pass the main checkout's
file from a worktree). Checks: masses strictly increasing; Mab_r and Mab_F146 decrease with mass
(reversals up to 0.05 mag between neighbours are allowed and reported); over 0.10-1.00 Msun within
0.5 mag of lens_ml.dat component 0 (the thin disc) in r and F146. Prints the lowest and highest
mass and the worst difference; exits non-zero on failure.
"""
import os
import sys

import numpy as np

REVERSAL_MAX = 0.05     # mag, allowed between neighbouring rows
LENS_ML_TOL = 0.5       # mag, against lens_ml.dat component 0
M_LO, M_HI = 0.10, 1.00
fails = []


def check(ok, msg):
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        fails.append(msg)


def main():
    if len(sys.argv) not in (2, 3):
        raise SystemExit(__doc__)
    path = sys.argv[1]
    here = os.path.dirname(os.path.abspath(__file__))
    lens_ml = sys.argv[2] if len(sys.argv) == 3 else os.path.join(here, "..", "CMD", "components", "lens_ml.dat")
    for p in (path, lens_ml):
        if not os.path.isfile(p):
            print(f"FAIL missing file {p}")
            return 1
    t = np.loadtxt(path, comments="#")
    check(t.ndim == 2 and t.shape[1] == 8, f"{path}: 8 columns (mass + 7 bands), {len(t)} rows")
    if fails:
        return 1
    m, r, f = t[:, 0], t[:, 3], t[:, 7]
    check(bool(np.isfinite(t).all()), "every value finite")
    check(bool((np.diff(m) > 0).all()), "masses strictly increasing")
    check(m[0] >= 0.08, f"lowest mass {m[0]:.4f} >= 0.08")
    for name, col in (("Mab_r", r), ("Mab_F146", f)):
        rev = np.diff(col)
        bad = np.flatnonzero(rev > 0)
        for i in bad:
            print(f"  reversal in {name}: {m[i]:.3f} -> {m[i + 1]:.3f} Msun, +{rev[i]:.3f} mag")
        worst_rev = rev[bad].max() if bad.size else 0.0
        check(worst_rev <= REVERSAL_MAX,
              f"{name} decreases with mass ({bad.size} reversals, worst {worst_rev:.3f} mag <= {REVERSAL_MAX})")
    L = np.loadtxt(lens_ml, comments="#")
    L = L[L[:, 0] == 0]
    mc = 0.5 * (L[:, 1] + L[:, 2])
    sel = (m >= M_LO) & (m <= M_HI)
    worst = 0.0
    for name, col, j in (("r", r, 6), ("F146", f, 10)):
        d = col[sel] - np.interp(m[sel], mc, L[:, j])
        k = np.argmax(np.abs(d))
        print(f"  worst {name} difference to lens_ml.dat comp 0: {d[k]:+.3f} mag at {m[sel][k]:.3f} Msun")
        worst = max(worst, abs(d[k]))
    check(worst <= LENS_ML_TOL, f"within {LENS_ML_TOL} mag of lens_ml.dat over {M_LO}-{M_HI} Msun (worst {worst:.3f})")
    print(f"lowest mass {m[0]:.4f} Msun, highest {m[-1]:.4f} Msun, {len(m)} rows")
    print("PASS" if not fails else "FAIL")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
