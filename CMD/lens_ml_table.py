"""Main-sequence mass -> absolute magnitude table for LUMINOUS LENSES (Deviation 74, step M2).

WHY. The simulator draws a lens's mass from its population's mass function, with no star attached,
and until Deviation 74 never added the lens's own light. For the ordinary-star (bulge) population a
lens below the main-sequence turnoff is a living star; its light dilutes both the magnification and
the astrometric centroid shift. The per-component CMD files in components/ cannot supply its
magnitudes: they were built through apply_visibility_filter, which keeps only stars detectable at
their Besancon distance, so at low mass only anomalously bright stars survive (their median M_r at
0.09 Msun is ~9.9, against ~16 for real M dwarfs). This script goes back to the RAW Besancon
catalogue (no visibility filter), keeps dwarfs (CL = 5), and computes their absolute magnitudes with
the SAME MIST bolometric corrections (AB) that made the source magnitudes, so lens and source light
are on one system.

OUTPUT  components/lens_ml.dat: one row per (component, mass bin):
    comp  m_lo  m_hi  n  Mab_u Mab_g Mab_r Mab_i Mab_z Mab_y Mab_F146      (medians, AB)
comp: 0 thin disk (Besancon Pop 1-7), 1 bulge (Pop 10), 2 thick disk (Pop 8, 11), 3 halo (Pop 9) --
the GalacticComponent order of Bulge.h. Mass bins 0.08-1.00 Msun in 0.02 steps; below 0.08 Msun
(brown dwarfs) the lens is treated as dark.

    cd CMD && ../.roman/bin/python lens_ml_table.py [--input Besancon/bos9.dat]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.argv_saved = list(sys.argv)
here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)
# Import the BC class without running BolometricCorrection.py's module-level pipeline.
src = open("BolometricCorrection.py").read()
cut = src.index("# ---------------------------------------------------------------------------\n# Main")
ns = {}
exec(compile(src[:cut], "BolometricCorrection.py", "exec"), ns)
MIST = ns["MISTBolometricCorrection"]

COLS = ["Teff", "logg", "Pop", "Mass", "Mbol", "[M/H]", "[a/Fe]", "CL"]
COMP_OF_POP = {**{p: 0 for p in range(1, 8)}, 10: 1, 8: 2, 11: 2, 9: 3}   # as save_components()
BINS = np.round(np.arange(0.08, 1.0001, 0.02), 4)
BANDS = ["LSST_u", "LSST_g", "LSST_r", "LSST_i", "LSST_z", "LSST_y", "Roman_F146"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="Besancon/bos9.dat")
    ap.add_argument("--dmax", type=float, default=1.5, help="kpc; only dwarfs nearer than this")
    a = ap.parse_args()

    # WHY ONLY NEARBY DWARFS. The raw catalogue is itself magnitude-limited (V <= 29), so at bulge
    # distances only anomalously bright low-mass dwarfs survive: median Mv at 0.5 Msun is 9.06 for
    # stars beyond 6 kpc against 10.38 within 1 kpc (0.7 Msun: 7.34 vs 7.96). Within ~1.5 kpc the
    # limit does not bite, so those stars trace the model's own mass-luminosity relation. They are
    # disk stars; each component's magnitudes are then made at that component's median metallicity.
    with open(a.input) as f:
        header = f.readline().lstrip("#").split()
    near, meta = [], []
    for chunk in pd.read_csv(a.input, sep=r"\s+", skiprows=1, header=None, names=header,
                             usecols=COLS + ["Dist"], chunksize=3_000_000):
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

    # Same BC machinery as BolometricCorrection.py (AV = 0: extinction is applied in the simulator).
    rubin = MIST("Rubin"); rubin.input_data = df.reset_index(drop=True); rubin.interp()
    roman = MIST("F146");  roman.input_data = rubin.input_data; roman.interp()
    out = roman.input_data
    path = os.path.join("components", "lens_ml.dat")
    nrow = 0
    with open(path, "w") as f:
        f.write(f"# lens mass -> absolute magnitude (AB): median Teff, logg, Mbol of dwarfs (CL=5) within "
                f"{a.dmax} kpc in {a.input}, MIST BCs at each component's median [M/H]; "
                "lens_ml_table.py, Deviation 74\n")
        f.write("# comp m_lo m_hi n Mab_u Mab_g Mab_r Mab_i Mab_z Mab_y Mab_F146\n")
        for r in out.sort_values(["comp", "bin"]).itertuples(index=False):
            mags = [getattr(r, x) for x in BANDS]
            if not np.all(np.isfinite(mags)):
                continue
            f.write(f"{int(r.comp)} {BINS[r.bin]:.2f} {BINS[r.bin + 1]:.2f} {int(r.n)} "
                    + " ".join(f"{v:.4f}" for v in mags) + "\n")
            nrow += 1
    print(f"wrote {path}: {nrow} rows")
    t = pd.read_csv(path, sep=r"\s+", comment="#", header=None,
                    names=["comp", "m_lo", "m_hi", "n"] + ["u", "g", "r", "i", "z", "y", "F146"])
    for comp in range(4):
        s = t[t.comp == comp]
        print(f"comp {comp}: {len(s)} bins (n per bin {s.n.min()}-{s.n.max()}); M_r at 0.1/0.3/0.5/0.7/0.9 =",
              [round(float(np.interp(m, (s.m_lo + s.m_hi) / 2, s.r)), 2) for m in (0.1, 0.3, 0.5, 0.7, 0.9)])


if __name__ == "__main__":
    main()
