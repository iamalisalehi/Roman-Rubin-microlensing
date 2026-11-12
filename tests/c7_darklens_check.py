"""Check the follow-up labels and the dark-lens probability in analysis/romanlib.py.

    .roman/bin/python tests/c7_darklens_check.py [--ms PATH]

Synthetic one-row frames and a fixed seed. --ms (optional) is a real ms_hypothesis.dat; when given it
is also loaded and checked for shape and monotonicity. The per-component layout (comp mass bands) and the
older layout without comp are both exercised: the component of each event (`struc`) selects its relation.
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


def synthetic_table(comps=None, shift=0.0):
    """Monotone toy main sequence, 0.08-60 Msun: M_b = 14 - 0.7 k - 4 log10(M/0.08) + shift*comp, k the band
    index. comps=None writes the older layout (no comp column); else one block per listed component."""
    m = np.geomspace(0.08, 60.0, 40)
    f = tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False)
    f.write("# toy main sequence\n# " + ("comp " if comps else "") + "mass " + " ".join(f"Mab_{b}" for b in BANDS) + "\n")
    for c in (comps or [None]):
        cols = [m] + [14.0 - 0.7 * k - 4.0 * np.log10(m / 0.08) + (shift * c if c else 0.0) for k in range(len(BANDS))]
        if c is not None:
            cols = [np.full_like(m, c)] + cols
        np.savetxt(f, np.column_stack(cols), fmt="%.6f")
    f.close()
    return f.name


def frame(**kw):
    base = dict(struc=0, Ml=1.0, Dl=4.0, mbs0=20.0, fb0=0.5, mbs1=20.0, fb1=0.5, ndw_R=50,
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


# The dark-lens test as it was before the i and z bands (a verbatim copy of analysis/romanlib.py's
# dark_lens_probability at 1efbdda; only the module prefix R. added to the names it calls). The new function
# must reproduce it exactly on a table without the i/z columns, and never fall below it with them.
def _dark_lens_probability_ref(df, partition, ms, n_samples=2000, seed=12345, ref_bands=None):
    """P(dark) per event: the fraction of posterior samples in which the blend cannot hide the lens.

    partition "joint" tests Rubin's reference band(s) (sigfb0_J > 0) and F146 (sigfb1_J > 0) with the
    joint mass error relMl_J; "rubin" the reference band(s) with sigfb0_L and relMl_L; "roman" F146
    with sigfb1_R and relMl_R. `ref_bands` are indices into R.FILTERS (default R.P.RUBIN_REF_BANDS);
    Rubin's baseline mbs0 and the lens flux are both summed over those bands together.

    Per sample: ln M ~ N(ln Ml, relMl); f_S ~ N(fb, sigfb) TRUNCATED to [0, 1] (the posterior under
    the physical prior; clipping would pile the mass of a poorly measured f_S on f_S = 1, zero blend
    light, and call every lens dark), independently per band (the joint fit's f_S correlation is
    not carried); a sigfb that is effectively zero gives f_S = fb; blend flux = 10^(-0.4 mbs) (1 - f_S); a
    main-sequence star of mass M at D_L has flux 10^(-0.4 (Mab(M) + 5 log10(D_L 100) + Alens)), Mab
    linear in ln M, the top row above the table, and a sample below the table's lowest mass is NOT
    excluded (a brown dwarf can hide anywhere). A sample is excluded when that flux exceeds the blend
    flux in ANY tested band. D_L is held at its true value (an approximation: the forecast has no
    distance error, and the mass error is the dominant one); the mass centre is the true Ml.

    The relation is the event's own Galactic component's (column `struc`, 0 thin disc, 1 bulge, 2 thick
    disc, 3 halo; the lens's component for every population, remnants included); an old-layout table
    has one relation for all.

    Float in [0, 1]; NaN (no verdict) when relMl <= 0 (the -1 sentinel), no band is tested, or a
    tested band's inputs are not finite. Vectorised over events in chunks; deterministic for a seed.
    """
    spec = {"joint": ("relMl_J", "sigfb0_J", "sigfb1_J"), "rubin": ("relMl_L", "sigfb0_L", None),
            "roman": ("relMl_R", None, "sigfb1_R")}
    if partition not in spec:
        raise ValueError(f"partition {partition!r}; expected one of {sorted(spec)}")
    relcol, s0col, s1col = spec[partition]
    rubin = list(R.P.RUBIN_REF_BANDS if ref_bands is None else ref_bands)
    col = lambda c: df[c].to_numpy(float)
    n = len(df)
    rel = col(relcol)
    lnM0 = np.log(col("Ml"))
    mu = 5.0 * np.log10(col("Dl") * 100.0)
    struc = col("struc")
    # per tested band: (fb, sigma, baseline mag, [(table column, lens-extinction column)])
    bands = []
    if s0col:
        bands.append((col("fb0"), col(s0col), col("mbs0"),
                      [(f"Mab_{R.FILTERS[b]}", f"Alens_{R.FILTERS[b]}") for b in rubin]))
    if s1col:
        bands.append((col("fb1"), col(s1col), col("mbs1"), [("Mab_F146", "Alens_F146")]))
    active = [b[1] > 0 for b in bands]
    ok = (rel > 0) & np.isfinite(rel) & np.isfinite(lnM0) & np.isfinite(mu) & np.any(active, axis=0)
    for (fb, sg, mbs, cols), act in zip(bands, active):
        good = np.isfinite(fb) & np.isfinite(sg) & np.isfinite(mbs)
        for _, ac in cols:
            good &= np.isfinite(col(ac))
        ok &= ~act | good
    ok &= np.isin(struc, (0, 1, 2, 3))
    rel_of = {c: R._ms_relation(ms, c) for c in (0, 1, 2, 3) if (ms["comp"] == -1).all() or (ms["comp"] == c).any()}
    ok &= np.isin(struc, list(rel_of))
    lnm_of = {c: np.log(g["mass"].to_numpy(float)) for c, g in rel_of.items()}
    out = np.full(n, np.nan)
    rng = np.random.default_rng(seed)
    idx = np.flatnonzero(ok)
    chunk = max(1, 4_000_000 // n_samples)
    for lo in range(0, idx.size, chunk):
        i = idx[lo:lo + chunk]
        lnM = lnM0[i, None] + rel[i, None] * rng.standard_normal((i.size, n_samples))
        excl = np.zeros(lnM.shape, bool)
        for (fb, sg, mbs, cols), act in zip(bands, active):
            t = act[i]
            if not t.any():
                continue
            fs = R._truncnorm01(fb[i, None], np.where(t, sg[i], 0.0)[:, None], rng.random(lnM.shape))
            blend = 10.0 ** (-0.4 * mbs[i, None]) * (1.0 - fs)
            lens = 0.0
            low = np.zeros(lnM.shape, bool)
            for mc, ac in cols:
                M = np.empty(lnM.shape)
                for c, g in rel_of.items():       # each event's own component relation
                    r = np.flatnonzero(struc[i] == c)
                    if r.size:
                        M[r] = np.interp(lnM[r], lnm_of[c], g[mc].to_numpy(float))
                        low[r] = lnM[r] < lnm_of[c][0]
                lens = lens + 10.0 ** (-0.4 * (M + mu[i, None] + col(ac)[i, None]))
            excl |= t[:, None] & ~low & (lens > blend)
        out[i] = excl.mean(axis=1)
    return pd.Series(out, index=df.index, name=f"pdark_{partition}")


IZ_COLS = ["sigfb_i_J", "sigfb_i_L", "sigfb_z_J", "sigfb_z_L"]


def with_iz(df, si_J=-1.0, si_L=-1.0, sz_J=-1.0, sz_L=-1.0, mbs_i=21.0, mbs_z=21.0, blend_i=0.4, blend_z=0.4):
    """A copy of `df` with the seven new columns (mbs_i, mbs_z, sigfb_{i,z}_{J,L}, nbandL) and the true i and z
    source flux fractions (blend_i, blend_z); scalars or per-row arrays."""
    d = df.copy()
    d["mbs_i"], d["mbs_z"] = mbs_i, mbs_z
    d["sigfb_i_J"], d["sigfb_i_L"], d["sigfb_z_J"], d["sigfb_z_L"] = si_J, si_L, sz_J, sz_L
    d["nbandL"] = 2
    d["blend_i"], d["blend_z"] = blend_i, blend_z
    return d


def mixed_frame():
    """Several rows with mixed sentinels, one per situation (old-layout columns only)."""
    rows = [dict(Ml=1.0, Dl=4.0, fb0=0.5, sigfb0_J=0.1, sigfb0_L=0.1, sigfb1_J=0.1, sigfb1_R=0.1),
            dict(Ml=3.0, Dl=5.0, fb0=0.3, sigfb0_J=0.2, sigfb0_L=-1.0, sigfb1_J=-1.0, sigfb1_R=0.2),
            dict(Ml=0.5, Dl=3.0, fb0=0.8, sigfb0_J=-1.0, sigfb0_L=0.05, sigfb1_J=0.05, sigfb1_R=-1.0),
            dict(Ml=2.0, Dl=6.0, fb0=0.4, sigfb0_J=-1.0, sigfb0_L=-1.0, sigfb1_J=-1.0, sigfb1_R=-1.0),
            dict(Ml=8.0, Dl=2.0, fb0=0.6, sigfb0_J=0.3, sigfb0_L=0.3, sigfb1_J=0.3, sigfb1_R=0.3, relMl_J=-1.0),
            dict(Ml=1.5, Dl=7.0, fb0=0.2, sigfb0_J=5.0, sigfb0_L=5.0, sigfb1_J=5.0, sigfb1_R=5.0, relMl_L=5.0)]
    return pd.concat([frame(**r) for r in rows], ignore_index=True)


def same_series(a, b):
    a, b = a.to_numpy(float), b.to_numpy(float)
    return bool(np.array_equal(a, b, equal_nan=True))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ms", default=None, help="a real ms_hypothesis.dat to load-check")
    args = ap.parse_args()

    # the older layout (no comp column) is one relation for every component; the tests below run on it
    path = synthetic_table()
    ms = R.load_ms_hypothesis(path)
    check(list(ms.columns) == ["comp", "mass"] + [f"Mab_{b}" for b in BANDS] and len(ms) == 40
          and bool((ms["comp"] == -1).all()) and bool(np.all(np.diff(ms["mass"]) > 0)),
          "load_ms_hypothesis: old layout loads (comp = -1), 40 rows increasing")
    if args.ms:
        real = R.load_ms_hypothesis(args.ms)
        ok = real["mass"].min() < 0.1 and real["mass"].max() > 50
        if (real["comp"] >= 0).all():
            ok &= sorted(real["comp"].unique()) == [0, 1, 2, 3]
        check(bool(ok), f"real table: {len(real)} rows, {real['mass'].min()}-{real['mass'].max()} Msun")

    # per-component layout: component c is c * 1 mag fainter; each event uses its own relation
    pc = synthetic_table(comps=[0, 1, 2, 3], shift=1.0)
    msc = R.load_ms_hypothesis(pc)
    check(len(msc) == 160 and sorted(msc["comp"].unique()) == [0, 1, 2, 3], "per-component layout: 4 x 40 rows")
    mass_c = float(msc["mass"].iloc[17])
    row = lambda c: frame(struc=c, Ml=mass_c, Dl=4.0, mbs1=20.0, fb1=0.5, sigfb1_R=1e-9, relMl_R=1e-6)
    fl0 = 10 ** (-0.4 * (float(msc[(msc.comp == 0)].iloc[17]["Mab_F146"]) + 5 * np.log10(400.0)))
    for c in (0, 3):
        d = frame(struc=c, Ml=mass_c, Dl=4.0, fb1=0.5, sigfb1_R=1e-9, relMl_R=1e-6,
                  mbs1=-2.5 * np.log10(1.5 * fl0 * 10 ** (-0.4 * c) / 0.5))
        # blend = 1.5 x comp-c lens flux => not excluded; with the comp-0 flux (brighter than comp 3's) it would be at c=3
        p = float(R.dark_lens_probability(d, "roman", msc).iloc[0])
        check(p == 0.0, f"per-component: struc {c}, blend 1.5x its own relation's flux: P(dark) = {p}")
    d3 = frame(struc=3, Ml=mass_c, Dl=4.0, fb1=0.5, sigfb1_R=1e-9, relMl_R=1e-6,
               mbs1=-2.5 * np.log10(0.9 * fl0 / 0.5))
    d0 = d3.assign(struc=0)
    p3, p0 = (float(R.dark_lens_probability(d, "roman", msc).iloc[0]) for d in (d3, d0))
    check(p3 == 0.0 and p0 == 1.0, f"same blend (0.9x comp-0 flux): comp 3 (3 mag fainter) P = {p3}, comp 0 P = {p0}")
    dn = d0.assign(struc=7)
    check(np.isnan(float(R.dark_lens_probability(dn, "roman", msc).iloc[0])), "struc outside 0-3: NaN")
    os.unlink(pc)

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

    # 4b. an unmeasured f_S must not make every lens dark (truncated, not clipped, f_S).
    # A 0.5 Msun lens whose blend is 100x its own main-sequence light; sigma(f_S) = 10 in each band.
    mass_h = float(ms["mass"].iloc[8])
    fl_r, _ = ms_flux(ms, 8, 4.0, {}, ["r"])
    fl_f, _ = ms_flux(ms, 8, 4.0, {}, ["F146"])
    for part, s0, s1 in (("roman", -1.0, 10.0), ("rubin", 10.0, -1.0), ("joint", 10.0, 10.0)):
        d = frame(Ml=mass_h, Dl=4.0, fb0=0.5, fb1=0.5, relMl_J=1e-6, relMl_L=1e-6, relMl_R=1e-6,
                  mbs0=-2.5 * np.log10(100 * fl_r / 0.5), mbs1=-2.5 * np.log10(100 * fl_f / 0.5),
                  sigfb0_J=s0, sigfb0_L=s0, sigfb1_J=s1, sigfb1_R=s1)
        p = pdark(d, part, ref_bands=[2])
        check(p < 0.1, f"{part}: blend 100x the lens light, f_S unmeasured (sigma 10): P(dark) = {p:.3f} < 0.1")
    check(R._truncnorm01(0.5, 0.0, 0.3) == 0.5 and 0 <= R._truncnorm01(1.0, 1e-300, 0.9) <= 1,
          "truncated f_S: sigma -> 0 gives f_S = fb, no NaN")

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
    df = pd.DataFrame([(0.9, 0.9, 80), (0.3, 0.9, 80), (0.9, 0.3, 80), (0.3, 0.3, 80), (0.3, 0.9, 0),
                       (0.9, 0.3, 0)],
                      columns=["fb0", "fb1", "ndw_R"])
    want = ["clean-both", "resolved-by-roman", "blended-roman-only", "blended-both", "rubin-blended-no-roman",
            "clean-rubin-no-roman"]
    got = list(R.blend_class(df))
    check(got == want, f"blend_class one row per class: {got}")
    edge = R.blend_class(pd.DataFrame({"fb0": [0.6, 0.3], "fb1": [0.6, 0.9], "ndw_R": [3, 2]}))
    check(list(edge) == ["clean-both", "rubin-blended-no-roman"],
          "blend_class: f_S = 0.6 is clean; ndw_R = 2 (< kMinTeleEpochs) is no Roman data")
    nr = R.blend_class(pd.DataFrame({"fb0": [0.9, 0.9], "fb1": [0.9, 0.1], "ndw_R": [0, 1]}))
    check(list(nr) == ["clean-rubin-no-roman"] * 2, "blend_class: Rubin clean without Roman data ignores fb1")
    check(R.BLEND_CLEAN_FS == 0.6 and R.DARK_P_MIN == 0.75, "constants 0.6 and 0.75")

    # 8. i and z bands in the dark-lens test (rulings R-a to R-d)
    mf = mixed_frame()
    for part in ("joint", "rubin", "roman"):
        check(same_series(R.dark_lens_probability(mf, part, ms), _dark_lens_probability_ref(mf, part, ms)),
              f"(a) {part}: a table without the i/z columns reproduces the old function exactly")
        check(same_series(R.dark_lens_probability(with_iz(mf), part, ms), _dark_lens_probability_ref(mf, part, ms)),
              f"(d) {part}: the i/z columns present, all four sigmas -1: bit-identical to the old function")

    # (b) r not fitted, i fitted: a verdict where the old function had none
    fi, mi = ms_flux(ms, 20, 5.0, {}, ["i"])
    fbi = 0.9

    def i_row(ratio, **kw):
        # blend flux in i = `ratio` x the lens's main-sequence flux in i at the event's mass (f_S,i = 0.9)
        d = frame(Ml=mi, Dl=5.0, fb0=0.5, relMl_J=1e-3, relMl_L=1e-3, sigfb0_J=-1.0, sigfb0_L=-1.0, **kw)
        return with_iz(d, si_J=1e-9, si_L=1e-9, mbs_i=-2.5 * np.log10(ratio * fi / (1 - fbi)), blend_i=fbi)

    for part in ("rubin", "joint"):
        d = i_row(0.5)
        p, p_old = float(R.dark_lens_probability(d, part, ms).iloc[0]), float(_dark_lens_probability_ref(d, part, ms).iloc[0])
        check(np.isfinite(p) and p == 1.0, f"(b) {part}: r unfitted, i fitted, lens twice the blend: P(dark) = {p} (finite, 1)")
        check(np.isnan(p_old), f"(b) {part}: the old function has no verdict on this row ({p_old})")
        p = float(R.dark_lens_probability(i_row(2.0), part, ms).iloc[0])
        check(p == 0.0, f"(b) {part}: same row, blend twice the lens: P(dark) = {p} (the i band discriminates)")

    # (c) r loose (sigma 10: f_S nearly uniform), i tight and too faint to hide the lens
    fr20, _ = ms_flux(ms, 20, 5.0, {}, ["r"])
    dc = frame(Ml=mi, Dl=5.0, fb0=0.5, relMl_J=1e-3, relMl_L=1e-3, sigfb0_J=10.0, sigfb0_L=10.0,
               mbs0=-2.5 * np.log10(100 * fr20 / 0.5))
    dci = with_iz(dc, si_J=1e-9, si_L=1e-9, mbs_i=-2.5 * np.log10(0.5 * fi / (1 - fbi)), blend_i=fbi)
    for part in ("rubin", "joint"):
        p_r = float(R.dark_lens_probability(dc, part, ms, ref_bands=[2]).iloc[0])
        p_ri = float(R.dark_lens_probability(dci, part, ms, ref_bands=[2]).iloc[0])
        check(p_ri - p_r > 0.9, f"(c) {part}: r loose, i tight and faint: P(dark) {p_r:.3f} (r alone) -> {p_ri:.3f} (r + i), gain > 0.9")

    # (e) a random multi-row frame, i/z active on some rows: never below the old function, exactly
    rg = np.random.default_rng(7)
    n = 80
    rnd = pd.concat([frame(Ml=float(np.exp(rg.uniform(np.log(0.15), np.log(20)))), Dl=float(rg.uniform(1.5, 8)),
                           fb0=float(rg.uniform(0.05, 0.95)), fb1=float(rg.uniform(0.05, 0.95)),
                           mbs0=float(rg.uniform(18, 24)), mbs1=float(rg.uniform(18, 24)),
                           relMl_J=float(rg.choice([-1.0, 0.05, 0.3])), relMl_L=float(rg.choice([-1.0, 0.1, 0.4])),
                           relMl_R=float(rg.choice([0.1, 0.5])),
                           sigfb0_J=float(rg.choice([-1.0, 0.02, 0.2])), sigfb0_L=float(rg.choice([-1.0, 0.05, 0.3])),
                           sigfb1_J=float(rg.choice([-1.0, 0.05])), sigfb1_R=float(rg.choice([-1.0, 0.1])))
                     for _ in range(n)], ignore_index=True)
    pick = lambda: rg.choice([-1.0, 0.01, 0.1, 0.5], size=n)
    rz = with_iz(rnd, si_J=pick(), si_L=pick(), sz_J=pick(), sz_L=pick(), mbs_i=rg.uniform(18, 24, n),
                 mbs_z=rg.uniform(18, 24, n), blend_i=rg.uniform(0.05, 0.95, n), blend_z=rg.uniform(0.05, 0.95, n))
    for part in ("joint", "rubin", "roman"):
        new = R.dark_lens_probability(rz, part, ms).to_numpy(float)
        old = _dark_lens_probability_ref(rnd, part, ms).to_numpy(float)
        fin = np.isfinite(old)
        check(bool(np.isfinite(new[fin]).all()) and bool((new[fin] >= old[fin]).all()),
              f"(e) {part}: {fin.sum()} rows with an old verdict, P_new >= P_old on every one, exactly")
        if part != "roman":
            check(int((new[fin] > old[fin]).sum()) > 3,
                  f"(e) {part}: the i/z bands raise P(dark) on {(new[fin] > old[fin]).sum()} rows (the check is not vacuous)")
        else:
            check(same_series(R.dark_lens_probability(rz, part, ms), _dark_lens_probability_ref(rnd, part, ms)),
                  "(e) roman: unchanged by the i/z columns")

    os.unlink(path)
    print()
    if fails:
        print(f"FAIL: {len(fails)} check(s)")
        sys.exit(1)
    print("PASS: all checks")


if __name__ == "__main__":
    main()
