"""Check CMD/components/ms_hypothesis.dat (lens_ml_table.py --hypothesis).

    python tests/ms_hypothesis_check.py PATH [--ref LENS_ML] [--deployed LENS_ML]

Checks: masses strictly increasing; every band (u g r i z y F146) decreases with mass (reversals up to
0.05 mag between neighbours are allowed; every reversal is reported); over 0.10-1.00 Msun within 0.5 mag
of component 0 (the thin disc) of the reference lens table in r and F146. --ref defaults to the empirical
lens table, CMD/components_staging/bos10_empirical_fillnone/lens_ml.dat of this repository
(the thin disc on the same Pecaut & Mamajek relation). The difference to the deployed lens-light table
(--deployed, default CMD/components/lens_ml.dat; skipped with a printed line when absent) is printed for
information only.
Prints the lowest and highest mass; exits non-zero on failure.
"""
import argparse
import os
import sys

import numpy as np

REVERSAL_MAX = 0.05     # mag, allowed between neighbouring rows
LENS_ML_TOL = 0.5       # mag, against lens_ml.dat component 0
M_LO, M_HI = 0.10, 1.00
CMD = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "CMD")
BANDS = ["u", "g", "r", "i", "z", "y", "F146"]
fails = []


def check(ok, msg):
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        fails.append(msg)


def compare(m, t, lens_ml):
    """Worst (difference, mass) of the table to component 0 of lens_ml over M_LO-M_HI in r and F146."""
    L = np.loadtxt(lens_ml, comments="#")
    L = L[L[:, 0] == 0]
    mc = 0.5 * (L[:, 1] + L[:, 2])
    sel = (m >= M_LO) & (m <= M_HI)
    out = []
    for name, col, j in (("r", 3, 6), ("F146", 7, 10)):
        d = t[sel, col] - np.interp(m[sel], mc, L[:, j])
        k = np.argmax(np.abs(d))
        out.append((name, d[k], m[sel][k]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path")
    ap.add_argument("--ref", default=os.path.join(CMD, "components_staging", "bos10_empirical_fillnone", "lens_ml.dat"))
    ap.add_argument("--deployed", default=os.path.join(CMD, "components", "lens_ml.dat"))
    a = ap.parse_args()
    for p in (a.path, a.ref):
        if not os.path.isfile(p):
            print(f"FAIL missing file {p}")
            return 1
    t = np.loadtxt(a.path, comments="#")
    check(t.ndim == 2 and t.shape[1] == 8, f"{a.path}: 8 columns (mass + 7 bands), {len(t)} rows")
    if fails:
        return 1
    m = t[:, 0]
    check(bool(np.isfinite(t).all()), "every value finite")
    check(bool((np.diff(m) > 0).all()), "masses strictly increasing")
    check(m[0] >= 0.08, f"lowest mass {m[0]:.4f} >= 0.08")
    for j, name in enumerate(BANDS, start=1):
        rev = np.diff(t[:, j])
        bad = np.flatnonzero(rev > 0)
        for i in bad:
            print(f"  reversal in Mab_{name}: {m[i]:.3f} -> {m[i + 1]:.3f} Msun, +{rev[i]:.3f} mag")
        worst_rev = rev[bad].max() if bad.size else 0.0
        check(worst_rev <= REVERSAL_MAX,
              f"Mab_{name} decreases with mass ({bad.size} reversals, worst {worst_rev:.3f} mag <= {REVERSAL_MAX})")
    worst = 0.0
    for name, d, mm in compare(m, t, a.ref):
        print(f"  worst {name} difference to the reference lens table: {d:+.3f} mag at {mm:.3f} Msun")
        worst = max(worst, abs(d))
    check(worst <= LENS_ML_TOL, f"within {LENS_ML_TOL} mag of {os.path.basename(os.path.dirname(a.ref))}/lens_ml.dat "
                                f"over {M_LO}-{M_HI} Msun (worst {worst:.3f})")
    if os.path.isfile(a.deployed):
        for name, d, mm in compare(m, t, a.deployed):
            print(f"  INFO deployed lens-light table (besancon dwarfs), not asserted: worst {name} {d:+.3f} mag at {mm:.3f} Msun")
    else:
        print(f"deployed table not found: skipped ({a.deployed})")
    print(f"lowest mass {m[0]:.4f} Msun, highest {m[-1]:.4f} Msun, {len(m)} rows")
    print("PASS" if not fails else "FAIL")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
