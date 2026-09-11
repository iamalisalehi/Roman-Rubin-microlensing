#!/usr/bin/env python3
"""Independent recomputation of Step H5's astrometric numbers, as a check on h5_astrometric_shift.py.

WHY A SECOND SCRIPT
-------------------
Every number in PROGRESS.md 5f is a claim about the astrometric deflection, and one code path
agreeing with itself is not evidence. This script reads the same table and recomputes the same
quantities from the physics, written independently of h5_astrometric_shift.py. Where the two
disagree, one of them is wrong and the disagreement is the finding.

It also carries the two checks that are statements about physics rather than statistics, and
that a result has to pass before it is quotable:

  1. sigma_joint <= sigma_single, EVENT BY EVENT. F[SJOINT] == F[SRUBIN] + F[SROMAN] holds
     element by element (Step C5), so the joint matrix sees strictly more information and its
     marginalised error cannot exceed either single-survey error. A violation is a bug.

  2. PRECISION IMPROVES WITH THE SIZE OF THE DEFLECTION. If sigma(theta_E)/theta_E did not fall
     as the peak centroid shift grows, the astrometric matrix would not be responding to the
     deflection at all, and the forecast would be measuring something else.

THE PHYSICS
-----------
A point lens displaces the centroid of the source's (unresolved) images from the true source
position by

    delta(u) = theta_E * u / (u^2 + 2)          [mas]

which is NOT monotonic: it rises from zero, peaks at u = sqrt(2), and falls again. The
trajectory sweeps u from large values down to u0 and back, so the largest deflection the event
ever reaches is

    theta_E / sqrt(8)              if u0 <= sqrt(2)   -- the trajectory crosses the peak
    theta_E * u0 / (u0^2 + 2)      if u0 >  sqrt(2)   -- it never gets there; closest is best

Two consequences worth holding on to: a high-magnification event (small u0) is a POOR
astrometric event, because the shift vanishes as u0 -> 0; and the lens mass
Ml = theta_E/(kappa piE) needs one observable from each matrix, which is why it is the only
quantity where a change in either shows up as a change in the science.

INPUT
-----
A whitespace column extract of test5.dat. The column rule is the one every script here has to
respect: the header line begins with a lone '#', so a data row's awk index = header position - 1
(verified empirically, PROGRESS.md). Produce the extract with:

    awk '/^#/ {next} ($38==1 || $39==1 || $40==1) {
      print $3,$5,$6,$8,$9,$11,$12,$36,$37,$38,$39,$40,$41,$42,$43,$47,$48,$49,$50,$51,$52,
            $78,$79,$80,$81,$82,$83,$84,$85,$86,$89,$91,$92 }' test5.dat

with this header prepended:

    tE piE tetE u0 Ml Dl Ds ndwL ndwR detL detR detJ okA_J okA_L okA_R sigpiE_J sigpiE_L
    sigpiE_R sigtetE_J sigtetE_L sigtetE_R relMl_J relMl_L relMl_R okB_J okB_L okB_R condB_J
    condB_L condB_R w_area nepL_pk nepR_pk magb_F146

Add `$70` (magb_F146) to the awk print list to get the per-source precision; without it the
script says so and falls back to the bright-source floor rather than quietly using the wrong
number.

Needs the project venv (.roman/bin/python) -- /usr/bin/python3 has numpy but no pandas.
"""
import argparse
import sys

import numpy as np
import pandas as pd

# errRomanA() from Bulge.h / helper.cpp, Step H4. Per EXPOSURE, in mas.
#
# The 1.1 mas that gets quoted as "Roman's astrometric precision" is the BRIGHT-SOURCE
# centroiding floor -- 1% of the 110 mas pixel -- and it does not apply to a GBTDS microlensing
# source. Those are faint bulge main-sequence stars, median F146 ~ 23, where the background
# term dominates and the real per-exposure precision is several mas. Using the floor as if it
# were the typical value understates the noise by a factor of six; this script did exactly that
# in its first version.
ROMAN_AST_FLOOR = 1.1      # mas, the floor itself
ROMAN_AST_MFLR = 20.62     # mag, below which the floor binds
ROMAN_AST_MBKG = 23.5      # mag, above which background dominates
ROMAN_AST_SLOPE_SRC = 0.3329
ROMAN_AST_SBKG = 10.0
ROMAN_AST_SLOPE_BKG = 0.4


def roman_ast_error(mag):
    """Per-exposure astrometric error at F146 magnitude `mag`, in mas."""
    mag = np.asarray(mag, dtype=float)
    out = np.full(mag.shape, ROMAN_AST_FLOOR)
    mid = (mag > ROMAN_AST_MFLR) & (mag <= ROMAN_AST_MBKG)
    out[mid] = ROMAN_AST_FLOOR * 10.0 ** (ROMAN_AST_SLOPE_SRC * (mag[mid] - ROMAN_AST_MFLR))
    hi = mag > ROMAN_AST_MBKG
    out[hi] = ROMAN_AST_SBKG * 10.0 ** (ROMAN_AST_SLOPE_BKG * (mag[hi] - ROMAN_AST_MBKG))
    return out


def load(path):
    d = pd.read_csv(path, sep=r'\s+')
    r2 = np.sqrt(2.0)
    u0 = d.u0.to_numpy()
    tetE = d.tetE.to_numpy()
    # The peak deflection reachable on this event's trajectory -- see the docstring.
    d['shift_max'] = np.where(u0 <= r2, tetE / np.sqrt(8.0), tetE * u0 / (u0 * u0 + 2.0))
    return d


def report_signal(cov):
    print("\n--- the centroid deflection, over Roman-observed detections ---")
    for q in (0.05, 0.25, 0.50, 0.75, 0.95):
        print(f"  {q:4.0%}  theta_E {np.quantile(cov.tetE, q):8.4f} mas"
              f"   max shift {np.quantile(cov.shift_max, q):8.4f} mas")
    r2 = np.sqrt(2.0)
    print(f"  fraction with u0 <= sqrt(2), i.e. reaching the astrometric peak: "
          f"{(cov.u0 <= r2).mean():.1%}")
    n = int(cov.ndwR.median())
    if "magb_F146" in cov.columns:
        sig1 = roman_ast_error(cov.magb_F146.to_numpy())
        med1 = float(np.median(sig1))
        above = float((cov.shift_max.to_numpy() > sig1).mean())
        src = "per source, errRomanA at its own F146"
    else:
        med1 = ROMAN_AST_FLOOR
        above = float('nan')
        src = ("NO magb_F146 COLUMN -- falling back to the bright-source FLOOR, which is NOT "
               "the typical value for these faint sources; re-extract with $70 to fix")
    print(f"\n  Roman per-exposure astrometric precision: {med1:.4f} mas   ({src})")
    print(f"  median peak signal / single-exposure precision: {cov.shift_max.median()/med1:.4f}")
    if above == above:  # not NaN
        print(f"  fraction of events whose peak shift exceeds ONE exposure: {above:.4%}")
    print(f"  Roman exposures per event (ndwR): median {n:,}")
    print(f"  an INDEPENDENT error averages to {med1/np.sqrt(n):.4f} mas over that many,")
    print(f"  putting the median signal at {cov.shift_max.median()/(med1/np.sqrt(n)):.2f} sigma")
    print("  -- that independence is an assumption, and the dominant caveat here. OPEN_ITEMS.md")


def report_tetE(cov):
    print("\n--- fractional precision on theta_E, where the astrometric matrix inverted ---")
    for lab, sc, oc in (("joint", "sigtetE_J", "okB_J"),
                        ("Rubin only", "sigtetE_L", "okB_L"),
                        ("Roman only", "sigtetE_R", "okB_R")):
        m = (cov[oc] == 1) & (cov[sc] > 0)
        if not m.sum():
            print(f"  {lab:11s} none")
            continue
        r = (cov.loc[m, sc] / cov.loc[m, "tetE"]).to_numpy()
        print(f"  {lab:11s} n={m.sum():7,}  median {np.median(r):8.4f}"
              f"   <10%: {(r < 0.10).mean():6.1%}   <1%: {(r < 0.01).mean():6.1%}")

    print("\n  theta_E measured better than THR -- and by whom")
    for thr in (0.30, 0.10, 0.03):
        j = (cov.okB_J == 1) & (cov.sigtetE_J > 0) & ((cov.sigtetE_J / cov.tetE) < thr)
        l = (cov.okB_L == 1) & (cov.sigtetE_L > 0) & ((cov.sigtetE_L / cov.tetE) < thr)
        r = (cov.okB_R == 1) & (cov.sigtetE_R > 0) & ((cov.sigtetE_R / cov.tetE) < thr)
        print(f"    <{thr:4.0%}: joint {int(j.sum()):6,} ({j.mean():5.1%})"
              f"   Rubin {int(l.sum()):6,} ({l.mean():5.1%})"
              f"   Roman {int(r.sum()):6,} ({r.mean():5.1%})"
              f"   joint-only {int((j & ~l & ~r).sum()):5,}")


def report_mass(cov):
    # A finite sigma is NOT a measurement. Counting non-sentinel relMl values found zero masses
    # exclusive to the joint fit, which was meaningless: relMl_L exists for 99.9% of these
    # events with a median fractional error of 15. Thresholds, always.
    print("\n--- lens mass Ml = theta_E/(kappa piE): needs BOTH matrices ---")
    for thr in (1.00, 0.30, 0.10):
        j = (cov.relMl_J > 0) & (cov.relMl_J < thr)
        l = (cov.relMl_L > 0) & (cov.relMl_L < thr)
        r = (cov.relMl_R > 0) & (cov.relMl_R < thr)
        print(f"  <{thr:4.0%}: joint {int(j.sum()):6,} ({j.mean():5.1%})"
              f"   Rubin {int(l.sum()):6,} ({l.mean():5.1%})"
              f"   Roman {int(r.sum()):6,} ({r.mean():5.1%})"
              f"   joint gets it, NEITHER alone {int((j & ~l & ~r).sum()):5,}")


def check_monotone_information(cov):
    """sigma_joint <= sigma_single, event by event. A violation is a bug, not a result."""
    print("\n--- CHECK 1: sigma_joint <= sigma_single, event by event ---")
    bad = 0
    for sc, oc in (("sigtetE_L", "okB_L"), ("sigtetE_R", "okB_R")):
        m = (cov.okB_J == 1) & (cov[oc] == 1) & (cov.sigtetE_J > 0) & (cov[sc] > 0)
        viol = int((cov.loc[m, "sigtetE_J"] > cov.loc[m, sc] * (1.0 + 1e-9)).sum())
        bad += viol
        print(f"  sigtetE_J vs {sc}: n={int(m.sum()):7,}  violations {viol}")
    print("  PASS" if bad == 0 else f"  *** FAILED: {bad} violations ***")
    return bad == 0


def check_responds_to_deflection(cov):
    """A bigger deflection must be measured better, or the matrix is not seeing it."""
    print("\n--- CHECK 2: precision improves with the size of the deflection ---")
    m = (cov.okB_J == 1) & (cov.sigtetE_J > 0)
    c = cov[m]
    rel = (c.sigtetE_J / c.tetE).to_numpy()
    q = pd.qcut(c.shift_max, 5, labels=False, duplicates='drop')
    for b in range(int(q.max()) + 1):
        s = q == b
        print(f"  quintile {b+1}: median shift {np.median(c.shift_max[s]):7.4f} mas"
              f"   median sigma/theta_E {np.median(rel[s]):8.4f}")
    r = float(np.corrcoef(np.log10(c.shift_max), np.log10(rel))[0, 1])
    print(f"  corr(log delta_max, log sigma/theta_E) = {r:+.3f}  (must be NEGATIVE)")
    print("  PASS" if r < 0 else "  *** FAILED ***")
    return r < 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("extract", help="column extract of test5.dat (see docstring)")
    a = ap.parse_args()

    d = load(a.extract)
    print(f"detected events: {len(d):,}")
    cov = d[d.ndwR > 0].copy()
    print(f"  of which Roman observed: {len(cov):,} ({len(cov)/len(d):.1%})")
    print("  Roman is the only astrometric instrument here worth the name, so everything")
    print("  below is over those events.")

    report_signal(cov)
    report_tetE(cov)
    report_mass(cov)
    ok = check_monotone_information(cov)
    ok = check_responds_to_deflection(cov) and ok
    print("\n" + ("ALL CHECKS PASSED" if ok else "*** A CHECK FAILED -- see above ***"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
