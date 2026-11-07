"""Check CMD/components/ms_hypothesis.dat (lens_ml_table.py --hypothesis).

    python tests/ms_hypothesis_check.py PATH [--lens-table LENS_ML]

The table holds, per Galactic component (column comp, 0-3), the lens table's rows below 1 Msun and the
Pecaut & Mamajek rows above. Checks: 9 columns (comp mass + 7 bands), all finite; per component masses strictly
increasing and every band decreasing with mass (reversals up to 0.05 mag between neighbours are allowed and
reported); below 1 Msun the rows equal --lens-table (default CMD/components/lens_ml.dat) at its bin centres,
per component, EXACTLY (to the 4 printed decimals); above 1 Msun the rows are the one Pecaut & Mamajek
sequence (a component may lack the rows its own faintest-so-far magnitudes make non-monotone). Prints
the step in each band at 1 Msun (last lens-table row to first sequence row), per component, and the lowest and
highest mass. Exits non-zero on failure.
"""
import argparse
import os
import sys

import numpy as np

REVERSAL_MAX = 0.05     # mag, allowed between neighbouring rows
JOIN = 1.0              # Msun
CMD = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "CMD")
BANDS = ["u", "g", "r", "i", "z", "y", "F146"]
fails = []


def check(ok, msg):
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        fails.append(msg)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path")
    ap.add_argument("--lens-table", default=os.path.join(CMD, "components", "lens_ml.dat"))
    a = ap.parse_args()
    for p in (a.path, a.lens_table):
        if not os.path.isfile(p):
            print(f"FAIL missing file {p}")
            return 1
    t = np.loadtxt(a.path, comments="#", ndmin=2)
    check(t.shape[1] == 9, f"{a.path}: 9 columns (comp, mass, 7 bands), {len(t)} rows")
    if fails:
        return 1
    check(bool(np.isfinite(t).all()), "every value finite")
    L = np.loadtxt(a.lens_table, comments="#")
    sub = {}
    for c in range(4):
        s = t[t[:, 0] == c]
        check(len(s) > 0, f"comp {c}: {len(s)} rows")
        if not len(s):
            continue
        m = s[:, 1]
        sub[c] = s
        check(bool((np.diff(m) > 0).all()), f"comp {c}: masses strictly increasing ({m[0]:.4f} to {m[-1]:.4f} Msun)")
        for j, name in enumerate(BANDS, start=2):
            rev = np.diff(s[:, j])
            bad = np.flatnonzero(rev > 0)
            for i in bad:
                print(f"  comp {c} reversal in Mab_{name}: {m[i]:.3f} -> {m[i + 1]:.3f} Msun, +{rev[i]:.3f} mag"
                      + ("  (lens table's own)" if m[i + 1] < JOIN else ""))
            hi = bad[m[bad + 1] > JOIN]
            worst = rev[hi].max() if hi.size else 0.0
            check(worst <= REVERSAL_MAX, f"comp {c}: Mab_{name} decreases with mass above {JOIN} Msun "
                                         f"({hi.size} reversals, worst {worst:.3f})")
        lc = L[L[:, 0] == c]
        want = np.column_stack([0.5 * (lc[:, 1] + lc[:, 2]), lc[:, 4:11]])
        got = s[s[:, 1] < JOIN][:, 1:]
        same = got.shape == want.shape and np.allclose(got, want, rtol=0, atol=6e-5)
        check(same, f"comp {c}: rows below {JOIN} Msun equal the lens table's ({len(got)} vs {len(want)} bin centres, "
                    f"max |diff| {np.abs(got - want).max() if got.shape == want.shape else float('nan'):.5f} mag)")
        first = s[s[:, 1] > JOIN][0] if (s[:, 1] > JOIN).any() else None
        if first is not None:
            step = first[2:] - got[-1, 1:]
            print(f"  comp {c}: step at {JOIN} Msun ({got[-1, 0]:.2f} -> {first[1]:.3f} Msun), M_AB: "
                  + " ".join(f"{b} {d:+.3f}" for b, d in zip(BANDS, step)))
    if len(sub) == 4:
        hi = [sub[c][sub[c][:, 1] > JOIN][:, 1:] for c in range(4)]
        seq = {}
        ok = True
        for h in hi:
            for r in h:
                ok &= bool(np.array_equal(seq.setdefault(round(r[0], 4), r), r))
        check(ok, f"above {JOIN} Msun every component's rows are the same sequence ({len(seq)} distinct masses; "
                  f"rows per component {[len(h) for h in hi]})")
    print(f"lowest mass {t[:, 1].min():.4f} Msun, highest {t[:, 1].max():.4f} Msun, {len(t)} rows")
    print("PASS" if not fails else "FAIL")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
