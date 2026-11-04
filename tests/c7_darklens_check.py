"""Check the follow-up labels and the dark-lens probability in analysis/romanlib.py.

    .roman/bin/python tests/c7_darklens_check.py [--ms PATH]

Synthetic one-row frames and a fixed seed. --ms (optional) is a real ms_hypothesis.dat; when given it
is also loaded and checked for shape and monotonicity.
Checks: a 10 Msun lens against a faint blend is dark with probability 1; a lens whose blend is 1.5x
its own main-sequence flux is not (P < 0.75); below the table's lowest mass nothing is excluded; -1
mass error gives NaN; absurd errors stay finite in [0, 1]; the Rubin flux sums the bands of
`ref_bands` (a blend 1% below the two-band sum is excluded, 1% above is not, and r alone would
differ); the joint verdict is "excluded in any tested band"; blend_class on one row per class.
Exits non-zero on failure.
"""
import argparse
import os
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "analysis"))
import romanlib as R          # noqa: E402

BANDS = ["u", "g", "r", "i", "z", "y", "F146"]
fails = []


def check(ok, msg):
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        fails.append(msg)


def synthetic_table():
    """Monotone toy main sequence, 0.08-60 Msun: M_b = 14 - 0.7 k - 4 log10(M/0.08), k the band index."""
    m = np.geomspace(0.08, 60.0, 40)
    cols = [m] + [14.0 - 0.7 * k - 4.0 * np.log10(m / 0.08) for k in range(len(BANDS))]
    f = tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False)
    f.write("# toy main sequence\n# mass " + " ".join(f"Mab_{b}" for b in BANDS) + "\n")
    np.savetxt(f, np.column_stack(cols), fmt="%.6f")
    f.close()
    return f.name


def frame(**kw):
    base = dict(Ml=1.0, Dl=4.0, mbs0=20.0, fb0=0.5, mbs1=20.0, fb1=0.5, ndw_R=50,
                relMl_J=0.1, relMl_L=0.1, relMl_R=0.1,
                sigfb0_J=-1.0, sigfb0_L=-1.0, sigfb1_J=-1.0, sigfb1_R=-1.0)
    for b in BANDS:
        base[f"Alens_{b}"] = 0.0
    base.update(kw)
    return pd.DataFrame([base])


def ms_flux(ms, idx, dl, alens, bands):
    """Hand-computed main-sequence flux at table row `idx` (an exact grid mass); returns (flux, mass)."""
    row = ms.iloc[idx]
    mu = 5.0 * np.log10(dl * 100.0)
    return sum(10 ** (-0.4 * (row[f"Mab_{b}"] + mu + alens.get(b, 0.0))) for b in bands), float(row["mass"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ms", default=None, help="a real ms_hypothesis.dat to load-check")
    args = ap.parse_args()

    path = synthetic_table()
    ms = R.load_ms_hypothesis(path)
    check(list(ms.columns) == ["mass"] + [f"Mab_{b}" for b in BANDS] and len(ms) == 40
          and bool(np.all(np.diff(ms["mass"]) > 0)),
          "load_ms_hypothesis: columns mass, Mab_<band>; 40 rows increasing")
    if args.ms:
        real = R.load_ms_hypothesis(args.ms)
        check(real["mass"].min() < 0.1 and real["mass"].max() > 50 and bool(np.all(np.diff(real["mass"]) > 0)),
              f"real table: {len(real)} rows, {real['mass'].min()}-{real['mass'].max()} Msun")

    def pdark(df, part, **kw):
        return float(R.dark_lens_probability(df, part, ms, **kw).iloc[0])

    # 1. a 10 Msun lens against a faint blend
    df = frame(Ml=10.0, Dl=4.0, Alens_F146=0.5, mbs1=20.0, fb1=0.95, sigfb1_R=0.01, relMl_R=0.1)
    check(pdark(df, "roman") == 1.0, "10 Msun at 4 kpc, faint blend: roman P(dark) == 1")

    # 2. blend flux = 1.5 x the lens's own main-sequence flux
    fl, mass = ms_flux(ms, 17, 4.0, {}, ["F146"])
    fb1 = 0.4
    p = pdark(frame(Ml=mass, Dl=4.0, mbs1=-2.5 * np.log10(1.5 * fl / (1 - fb1)), fb1=fb1,
                    sigfb1_R=0.02, relMl_R=0.05), "roman")
    check(p < 0.75, f"blend 1.5x the main-sequence flux: P(dark) = {p:.3f} < 0.75")
    p = pdark(frame(Ml=mass, Dl=4.0, mbs1=-2.5 * np.log10(0.6 * fl / (1 - fb1)), fb1=fb1,
                    sigfb1_R=0.02, relMl_R=0.05), "roman")
    check(p > 0.9, f"blend 0.6x the main-sequence flux: P(dark) = {p:.3f} > 0.9 (the test discriminates)")

    # 3. below the table
    p = pdark(frame(Ml=0.05, mbs1=12.0, fb1=0.99, sigfb1_R=0.001, relMl_R=0.05), "roman")
    check(p == 0.0, f"0.05 Msun, below the table: P(dark) = {p} == 0")

    # 4. sentinel and absurd errors
    check(np.isnan(pdark(frame(sigfb1_R=0.02, relMl_R=-1.0), "roman")), "relMl_R = -1: NaN")
    check(np.isnan(pdark(frame(sigfb1_R=0.02, relMl_R=0.0), "roman")), "relMl_R = 0: NaN")
    p = pdark(frame(sigfb1_R=10.0, relMl_R=5.0), "roman")
    check(np.isfinite(p) and 0.0 <= p <= 1.0, f"relMl_R = 5, sigfb1 = 10: P(dark) = {p:.3f} finite in [0, 1]")
    check(np.isnan(pdark(frame(relMl_R=0.1), "roman")), "sigfb1_R = -1 (no tested band): NaN")

    # 5. Rubin band sums the reference bands
    alens = {"r": 0.4, "i": 1.1}
    f2, m2 = ms_flux(ms, 20, 6.0, alens, ["r", "i"])
    f_r, _ = ms_flux(ms, 20, 6.0, alens, ["r"])
    check(f_r < 0.99 * f2, "toy table: r alone is more than 1% below r + i (the check can tell them apart)")
    fb0 = 0.3

    def rubin_row(scale):
        return frame(Ml=m2, Dl=6.0, mbs0=-2.5 * np.log10(scale * f2 / (1 - fb0)), fb0=fb0,
                     sigfb0_L=1e-9, relMl_L=1e-6, Alens_r=0.4, Alens_i=1.1)

    got = [pdark(rubin_row(s), "rubin", ref_bands=[2, 3]) for s in (0.99, 1.01)]
    check(got == [1.0, 0.0], f"ref_bands=[2, 3]: Rubin flux is the r + i hand sum (0.99, 1.01 -> {got})")
    check(pdark(rubin_row(0.99), "rubin", ref_bands=[2]) == 0.0,
          "ref_bands=[2] on the 0.99 row: r alone cannot exclude (a different answer)")
    check(pdark(rubin_row(0.99), "rubin") == pdark(rubin_row(0.99), "rubin", ref_bands=list(R.P.RUBIN_REF_BANDS)),
          "default ref_bands = P.RUBIN_REF_BANDS")

    # 6. joint: excluded in ANY tested band; a band with sigma <= 0 is not tested
    f1, m1 = ms_flux(ms, 25, 3.0, {}, ["F146"])
    fr, _ = ms_flux(ms, 25, 3.0, {}, ["r"])

    def mb(f, ratio):                  # baseline magnitude of a blend that is `ratio` x flux f (fb = 0.4)
        return -2.5 * np.log10(ratio * f / 0.6)

    def joint(r0, r1, s0, s1):
        return pdark(frame(Ml=m1, Dl=3.0, fb0=0.4, fb1=0.4, relMl_J=1e-6, mbs0=mb(fr, r0), mbs1=mb(f1, r1),
                           sigfb0_J=s0, sigfb1_J=s1), "joint")

    got = (joint(3.0, 3.0, 1e-9, 1e-9), joint(3.0, 0.3, 1e-9, 1e-9), joint(0.3, 3.0, -1.0, 1e-9),
           joint(0.3, 3.0, 1e-9, -1.0))
    check(got == (0.0, 1.0, 0.0, 1.0), f"joint: any-band exclusion, untested band ignored: {got} == (0, 1, 0, 1)")

    # 7. blend_class: one row per class
    df = pd.DataFrame([(0.9, 0.9, 80), (0.3, 0.9, 80), (0.9, 0.3, 80), (0.3, 0.3, 80), (0.3, -1.0, 0)],
                      columns=["fb0", "fb1", "ndw_R"])
    want = ["clean-both", "resolved-by-roman", "blended-roman-only", "blended-both", "rubin-blended-no-roman"]
    got = list(R.blend_class(df))
    check(got == want, f"blend_class one row per class: {got}")
    edge = R.blend_class(pd.DataFrame({"fb0": [0.6, 0.3], "fb1": [0.6, 0.9], "ndw_R": [3, 2]}))
    check(list(edge) == ["clean-both", "rubin-blended-no-roman"],
          "blend_class: f_S = 0.6 is clean; ndw_R = 2 (< kMinTeleEpochs) is no Roman data")
    check(R.BLEND_CLEAN_FS == 0.6 and R.DARK_P_MIN == 0.75, "constants 0.6 and 0.75")

    os.unlink(path)
    print()
    if fails:
        print(f"FAIL: {len(fails)} check(s)")
        sys.exit(1)
    print("PASS: all checks")


if __name__ == "__main__":
    main()
