"""Main-sequence mass -> absolute magnitude table for luminous lenses.

The simulator draws a lens's mass from its population's mass function, with no star attached. For
the ordinary-star (bulge) population a lens below the main-sequence turnoff is a living star, and its
light dilutes both the magnification and the astrometric centroid shift. The per-component CMD files
in components/ hold stars, not a mass-magnitude relation, so this script goes back to the raw
Besancon catalogue (no visibility filter), keeps dwarfs (CL = 5, Typ < 9), and computes their
absolute magnitudes with the same MIST bolometric corrections (AB) and BC-grid clamp as
BolometricCorrection.py, so lens and source light are on one system.

Output <out>: one row per (component, mass bin):
    comp  m_lo  m_hi  n  Mab_u Mab_g Mab_r Mab_i Mab_z Mab_y Mab_F146      (AB)
comp: 0 thin disk (Besancon Pop 1-7), 1 bulge (Pop 10), 2 thick disk (Pop 8, 11), 3 halo (Pop 9),
the GalacticComponent order of include/common.h. Mass bins 0.08-1.00 Msun in 0.02 steps, in the
star's true (catalogue) mass; below 0.08 Msun (brown dwarfs) the lens is treated as dark. n is the
number of the component's dwarfs in the bin. Comment lines ('#') above the table list what was done.

    cd CMD && ../.roman/bin/python lens_ml_table.py [--input Besancon/bos10] [--dwarfs MODE] [--out PATH]

  --input PATH      the Besancon catalogue (default Besancon/bos10). Rows whose field count differs from
                    the header's are skipped and counted, as in BolometricCorrection.py.
  --dwarfs MODE     empirical (default) or besancon. The thick disc (comp 2) and the halo (comp 3) are
                    always made the besancon way:
        besancon    per component and bin, the medians of that component's own dwarfs (Teff, logg and
                    Mbol over all distances; bos10 has no magnitude limit) and the MIST bolometric
                    corrections at the component's median [M/H] and [a/Fe] (clamped onto the BC grid,
                    as for the stars).
        empirical   for the thin disc and the bulge, the relation the sources use after
                    BolometricCorrection.py --dwarfs empirical: Mbol_E and Teff_E of the Pecaut &
                    Mamajek dwarf sequence at the bin-centre mass, logg from the builder's formula,
                    BCs at the component's median metallicity. Between 0.6 and 0.7 Msun the builder's
                    taper blends it with the component's own medians, (1 - t) own + t empirical,
                    t = clip((0.7 - m) / 0.1, 0, 1); from 0.7 Msun up the own medians are used.
  --out PATH        default components_staging/<input>_<dwarfs>_fillnone/lens_ml.dat. Writing
                    components/lens_ml.dat needs the explicit path. The fill-floor variants of the
                    builder share the lens table of their --dwarfs mode: pass --out for them.
  --dmax KPC        only for --input Besancon/bos9.dat (default 1.5), whose magnitude limit (V <= 29)
                    leaves only anomalously bright dwarfs far away: the table is then made from nearby
                    disc dwarfs alone, pooled over the components, and needs --dwarfs besancon.
                    Ignored for bos10.

A bin in which a component has no stars (the bulge below its mass floor of 0.156 Msun; the thick disc
and the halo below ~0.155; the halo's high-mass bins) takes the empirical relation at the component's
median metallicity, with n = 0; the output lists those bins in a comment line.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

ORIG_CWD = os.getcwd()
here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)
sys.path.insert(0, here)
import BolometricCorrection as BC      # noqa: E402  (importing builds nothing)

COLS = ["Teff", "logg", "Pop", "Mass", "Mbol", "[M/H]", "[a/Fe]", "CL", "Typ"]
COMP_OF_POP = {**{p: 0 for p in range(1, 8)}, 10: 1, 8: 2, 11: 2, 9: 3}   # as in BolometricCorrection.py
BINS = np.round(np.arange(0.08, 1.0001, 0.02), 4)
BANDS = ["LSST_u", "LSST_g", "LSST_r", "LSST_i", "LSST_z", "LSST_y", "Roman_F146"]
COMP_NAMES = ["thin_disk", "bulge", "thick_disk", "halo"]
EMPIRICAL_COMPS = (0, 1)               # the components whose dwarfs the empirical shift moves


def resolve(path):
    """A user path: relative to where the command was typed if it exists there, else to CMD/."""
    if os.path.isabs(path) or os.path.exists(os.path.join(ORIG_CWD, path)):
        return os.path.abspath(os.path.join(ORIG_CWD, path))
    return os.path.join(here, path)


def table_modern(a):
    """bos10-style table: own medians, the empirical relation where asked, BCs at the median metallicity."""
    cat, rep = BC.read_catalogue(a.input, COLS)
    comp = BC.component_codes(cat["Pop"])
    mq = BC.qmass(cat["Mass"])
    dwarf = (cat["CL"] == 5) & (cat["Typ"] < BC.DARK_TYP[0])
    in_range = (mq >= BINS[0] - BC.MASS_EPS) & (mq < BINS[-1] - BC.MASS_EPS)
    emp = BC.EmpiricalDwarfs()
    nb = len(BINS) - 1
    centres = 0.5 * (BINS[:-1] + BINS[1:])
    t = BC.taper(centres)
    rows, comments = [], []
    for c in range(4):
        sel = (comp == c) & dwarf & in_range
        mh = float(np.median(cat["M_H"][sel].astype(np.float64)))
        afe = float(np.median(cat["a_Fe"][sel].astype(np.float64)))
        mh_c, afe_c = BC.clamp_to_grid(np.array([mh]), np.array([afe]))
        ib = np.digitize(mq[sel], BINS) - 1
        g = pd.DataFrame({"b": ib, "Teff": cat["Teff"][sel].astype(np.float64), "logg": cat["logg"][sel].astype(np.float64),
                          "Mbol": cat["Mbol"][sel].astype(np.float64)}).groupby("b").agg(
            n=("Mbol", "size"), Teff=("Teff", "median"), logg=("logg", "median"), Mbol=("Mbol", "median")).reindex(range(nb))
        n = g["n"].fillna(0).to_numpy().astype(int)
        have = n > 0
        logT = np.log10(g["Teff"].to_numpy())
        logg = g["logg"].to_numpy().copy()
        mbol = g["Mbol"].to_numpy().copy()
        if a.dwarfs == "empirical" and c in EMPIRICAL_COMPS:
            logT = (1.0 - t) * logT + t * emp.logT(centres)
            mbol = (1.0 - t) * mbol + t * emp.mbol(centres)
            logg = np.where(t > 0, BC.logg_formula(centres, logT, mbol), logg)
        empty = ~have
        logT[empty], mbol[empty] = emp.logT(centres[empty]), emp.mbol(centres[empty])
        logg[empty] = BC.logg_formula(centres[empty], logT[empty], mbol[empty])
        rows.append(dict(c=c, n=n, logT=logT, logg=logg, mbol=mbol, mh=mh, afe=afe, mh_c=float(mh_c[0]),
                         afe_c=float(afe_c[0]), empty=empty, n_dwarfs=int(sel.sum())))
        if empty.any():
            e = np.flatnonzero(empty)
            comments.append(f"comp {c} ({COMP_NAMES[c]}): {len(e)} bins with no stars, n = 0, empirical relation: "
                            + _ranges(BINS, e))
    bcs = BC.BCSet()
    for r in rows:
        r["mags"] = bcs.mags(r["logT"], r["logg"], r["mbol"], np.full(nb, r["mh"]), np.full(nb, r["afe"]), clamp=True)
        if not np.isfinite(r["mags"]).all():
            raise RuntimeError(f"comp {r['c']}: a lens-light bin has no MIST BC after the grid clamp")
    return rows, comments, rep


def _ranges(bins, idx):
    """'0.08-0.16' style description of consecutive bin indices."""
    out, start = [], idx[0]
    for i, j in zip(idx, list(idx[1:]) + [None]):
        if j is None or j != i + 1:
            out.append(f"{bins[start]:.2f}-{bins[i + 1]:.2f}")
            start = j
    return ", ".join(out)


def write_modern(a, rows, comments, rep):
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    nb = len(BINS) - 1
    with open(a.out, "w") as f:
        f.write(f"# lens mass -> absolute magnitude (AB): median Teff, logg, Mbol of each component's own dwarfs "
                f"(CL = 5, Typ < 9; all distances) per TRUE-mass bin in {BC.shown(a.input)}, MIST BCs "
                "at each component's median [M/H], [a/Fe] (clamped onto the BC grid); lens_ml_table.py\n")
        f.write(f"# dwarfs = {a.dwarfs}: " + (
            "comps 0 (thin) and 1 (bulge) on the Pecaut & Mamajek (2013) dwarf sequence at the bin-centre mass "
            "below 0.6 Msun, blended with the component's own medians between 0.6 and 0.7 (the builder's taper), "
            "own medians above; comps 2, 3 own medians" if a.dwarfs == "empirical" else "all components: own medians")
            + f"; input rows skipped (field count != {rep['n_fields']}): {rep['n_skipped']}\n")
        f.write("# component median [M/H] [a/Fe] used for the BCs (clamped values after the arrow): "
                + "; ".join(f"comp {r['c']}: {r['mh']:.3f} {r['afe']:.3f} -> {r['mh_c']:.3f} {r['afe_c']:.3f}" for r in rows) + "\n")
        for line in comments:
            f.write(f"# {line}\n")
        f.write("# comp m_lo m_hi n Mab_u Mab_g Mab_r Mab_i Mab_z Mab_y Mab_F146\n")
        for r in rows:
            for b in range(nb):
                f.write(f"{r['c']} {BINS[b]:.2f} {BINS[b + 1]:.2f} {int(r['n'][b])} "
                        + " ".join(f"{v:.4f}" for v in r["mags"][b]) + "\n")
    print(f"wrote {a.out}: {4 * nb} rows")


def table_legacy(a):
    """Legacy bos9 table: nearby dwarfs only, pooled over the components."""
    near, meta = [], []
    rep = {}
    for chunk in BC.iter_catalogue_chunks(a.input, COLS + ["Dist"], 3_000_000, rep):
        c = chunk[(chunk.CL == 5) & chunk.Pop.isin(list(COMP_OF_POP))]
        meta.append(c.assign(comp=c.Pop.map(COMP_OF_POP))[["comp", "[M/H]", "[a/Fe]"]]
                    .sample(min(len(c), 200_000), random_state=1))
        near.append(c[(c.Dist < a.dmax) & (c.Mass >= BINS[0]) & (c.Mass < BINS[-1])])
    near = pd.concat(near, ignore_index=True)
    meta = pd.concat(meta, ignore_index=True)
    met = meta.groupby("comp")[["[M/H]", "[a/Fe]"]].median()
    print(f"{len(near)} dwarfs within {a.dmax} kpc; component median [M/H], [a/Fe]:\n{met}")
    near["bin"] = np.digitize(near.Mass, BINS) - 1
    rel = near.groupby("bin").agg(Teff=("Teff", "median"), logg=("logg", "median"),
                                  Mbol=("Mbol", "median"), n=("Mass", "size")).reset_index()
    rows = []
    for comp in range(4):
        r = rel.copy()
        r["[M/H]"], r["[a/Fe]"], r["comp"] = met.loc[comp, "[M/H]"], met.loc[comp, "[a/Fe]"], comp
        rows.append(r)
    df = pd.concat(rows, ignore_index=True)

    # Same BC machinery as BolometricCorrection.py (AV = 0; the simulator applies extinction).
    rubin = BC.MISTBolometricCorrection("Rubin"); rubin.input_data = df.reset_index(drop=True); rubin.interp()
    roman = BC.MISTBolometricCorrection("F146");  roman.input_data = rubin.input_data; roman.interp()
    out = roman.input_data
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    nrow = 0
    with open(a.out, "w") as f:
        f.write(f"# lens mass -> absolute magnitude (AB): median Teff, logg, Mbol of dwarfs (CL=5) within "
                f"{a.dmax} kpc in {BC.shown(a.input)}, MIST BCs at each component's median [M/H]; "
                "lens_ml_table.py\n")
        f.write("# comp m_lo m_hi n Mab_u Mab_g Mab_r Mab_i Mab_z Mab_y Mab_F146\n")
        for r in out.sort_values(["comp", "bin"]).itertuples(index=False):
            mags = [getattr(r, x) for x in BANDS]
            if not np.all(np.isfinite(mags)):
                continue
            f.write(f"{int(r.comp)} {BINS[r.bin]:.2f} {BINS[r.bin + 1]:.2f} {int(r.n)} "
                    + " ".join(f"{v:.4f}" for v in mags) + "\n")
            nrow += 1
    print(f"wrote {a.out}: {nrow} rows")


def summarise(path):
    t = pd.read_csv(path, sep=r"\s+", comment="#", header=None,
                    names=["comp", "m_lo", "m_hi", "n"] + ["u", "g", "r", "i", "z", "y", "F146"])
    for comp in range(4):
        s = t[t.comp == comp]
        print(f"comp {comp}: {len(s)} bins (n per bin {s.n.min()}-{s.n.max()}, {int((s.n == 0).sum())} empty); "
              "M_r at 0.1/0.3/0.5/0.7/0.9 =",
              [round(float(np.interp(m, (s.m_lo + s.m_hi) / 2, s.r)), 2) for m in (0.1, 0.3, 0.5, 0.7, 0.9)],
              " M_F146 =", [round(float(np.interp(m, (s.m_lo + s.m_hi) / 2, s.F146)), 2) for m in (0.1, 0.3, 0.5, 0.7, 0.9)])


def main(argv=None):
    ap = argparse.ArgumentParser(description="Lens-light mass -> absolute magnitude table (see the module docstring).",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--input", default=BC.DEFAULT_INPUT, help="Besancon catalogue (damaged rows are skipped and counted)")
    ap.add_argument("--dwarfs", choices=["besancon", "empirical"], default="empirical",
                    help="empirical: thin disc and bulge on the Pecaut & Mamajek sequence below 0.7 Msun; "
                         "besancon: every component's own medians")
    ap.add_argument("--out", default=None,
                    help=f"output table; default {BC.STAGING_ROOT}/<input>_<dwarfs>_fillnone/lens_ml.dat. "
                         "Writing components/lens_ml.dat needs it explicitly.")
    ap.add_argument("--dmax", type=float, default=None,
                    help="kpc; only for --input Besancon/bos9.dat (default 1.5): nearby dwarfs only")
    a = ap.parse_args(argv)
    a.input = resolve(a.input)
    legacy = os.path.basename(a.input).startswith("bos9")
    if a.out is None:
        a.out = os.path.join(here, BC.STAGING_ROOT, BC.variant_name(a.input, a.dwarfs, "none"), "lens_ml.dat")
    else:
        a.out = resolve(a.out)
    if legacy:
        if a.dwarfs != "besancon":
            raise SystemExit("--input bos9 is the legacy table (nearby, pooled dwarfs): pass --dwarfs besancon")
        a.dmax = 1.5 if a.dmax is None else a.dmax
        table_legacy(a)
    else:
        if a.dmax is not None:
            raise SystemExit("--dmax applies to bos9.dat only: bos10 has no magnitude limit")
        rows, comments, rep = table_modern(a)
        write_modern(a, rows, comments, rep)
    summarise(a.out)


if __name__ == "__main__":
    main()
