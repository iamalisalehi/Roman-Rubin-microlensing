"""Builds the source/neighbour catalogue CMD/components/{thin_disk,bulge,thick_disk,halo}.dat.

WHAT IT IS (Deviation 81, 2026-10-02). Each file is a sample of EVERY star of one Galactic component
in the Besancon catalogue (Besancon/bos9.dat: 0.1 deg^2 at (l, b) = (0.5, -1.4), 0-8 kpc, Av = 0),
with absolute AB magnitudes from MIST bolometric corrections. The simulator draws a source AND each
of its blend neighbours uniformly from these lists, with a distance from its own density model, and
counts them with Nstart = rho / <m>. A draw therefore stands for a random member of the population
Nstart counts, and the lists must be that population -- complete, not "the stars someone could see".

There is NO visibility filter. Until Deviation 81 one kept a star only if it was visible in F146 AND
>= 1 LSST band at its Besancon distance, against thre/satu parsed from Bulge.h. That excluded
Roman-only stars, stars visible only when magnified, stars saturated at their catalogue distance --
40% of the bulge -- and, because the neighbours come from the same lists, made every blend too
bright. Visibility is decided where it belongs: per event in the simulator, with magnification,
dust and each survey's own limits.

DARK ENTRIES. A star MIST cannot place on its grid (brown dwarfs, white dwarfs; Besancon Typ
9.0-9.2) is kept with every magnitude = DARK_MAG: it is part of what Nstart counts, so dropping it
would bias the list bright, but it contributes no light.

MEAN MASSES. <m> per component over the whole list (dark entries included) is printed and written
to components/provenance.txt; Bulge.h's MEANMASS_* must equal it, so that Nstart counts the very
population the draws come from.

LOW-MASS DWARFS OF THE BULGE, THICK DISC AND HALO (Deviation 82). Besancon's dwarfs in these
components are 2-3 mag brighter at 0.1 Msun than its thin-disc dwarfs (and than real M dwarfs). Below
MS_FIX_HI each such dwarf (CL = 5) is moved, in every band, by the difference between
components/lens_ml.dat -- the relation lens light already uses (lens_ml_table.py: nearby thin-disc
dwarfs, MIST BCs at the component's metallicity) -- and its component's own running median at that
mass; the star's scatter about the median is kept. The shift tapers linearly to zero between
MS_FIX_LO and MS_FIX_HI, above which stars may be evolving and are left alone. The thin disc is the
reference and is not touched. Run lens_ml_table.py first.

SUBSAMPLE. A component with more than MAX_ROWS stars is reduced to a uniform random subsample
(fixed seed): a luminosity function needs no more, and the simulator holds the lists in RAM.

    cd CMD && ../.roman/bin/python BolometricCorrection.py
"""
import glob
import numpy as np
import pandas as pd

from scipy.interpolate import RegularGridInterpolator


def mh_to_feh(mh, afe):
    return mh - np.log10(0.638 * 10**afe + 0.362)


INPUT_COLUMNS = ["Teff", "logg", "Pop", "Age", "Mass", "Mbol", "[M/H]", "[a/Fe]", "CL", "Typ", "Av"]
FILTER_ORDER = ["LSST_u", "LSST_g", "LSST_r", "LSST_i", "LSST_z", "LSST_y", "Roman_F146"]
DARK_MAG = 99.0                       # "no light": 10^(-0.4*99) is zero in every sum
DARK_TYP = (9.0, 9.2)                 # Besancon white dwarfs: no MIST track, faint -> dark
MAX_ROWS = 3_500_000                  # per component; only the bulge (5.9M) exceeds it
SUBSAMPLE_SEED = 20261002
# Besancon Pop codes -> component, as in Bulge.h's GalacticComponent order of the files.
COMPONENTS = {"thin_disk": list(range(1, 8)), "bulge": [10], "thick_disk": [8, 11], "halo": [9]}
MS_FIX_LO, MS_FIX_HI = 0.6, 0.7       # Msun: full shift below LO, none above HI (Deviation 82)
MS_FIX_COMP = {"bulge": 1, "thick_disk": 2, "halo": 3}   # lens_ml.dat comp codes
MS_FIX_BIN = 0.02                     # Msun, running-median bin
# Upper age bounds read_cmd() CHECKs (helper.cpp); a violation stops the build, it does not drop.
AGE_MAX = {"thin_disk": 10, "bulge": 10, "thick_disk": 13, "halo": 14}


class MISTBolometricCorrection:
    def __init__(self, phot):
        if phot.lower() in ["lsst", "rubin"]:
            self.phot = "lsst"
            directory = "./Rubin"
        elif phot.lower() in ["roman", "wfirst"]:
            self.phot = "roman"
            directory = "./Roman"
        elif phot.lower() in ["f_146", "f146"]:
            self.phot = "f146"
            directory = "./Roman"
        else:
            raise ValueError("Unknown photometric system")

        tables = []
        for fname in glob.glob(directory + "/*"):
            with open(fname, "r") as f:
                for i, line in enumerate(f):
                    if i == 3:
                        self.header_table = line.lstrip("#").split()
                        break

            tables.append(np.loadtxt(fname))

        self.table = pd.DataFrame(np.concatenate(tables, axis=0), columns=self.header_table)

        print("Loaded BC table")
        print(self.table.shape)

    def read_input(self, path):
        """Reads only the columns this pipeline uses, in chunks (bos9.dat is 2.7 GB)."""
        with open(path, "r") as f:
            self.header_input = f.readline().lstrip("#").split()
        parts = [c for c in pd.read_csv(path, sep=r"\s+", skiprows=1, header=None,
                                         names=self.header_input, usecols=INPUT_COLUMNS,
                                         chunksize=2_000_000)]
        self.input_data = pd.concat(parts, ignore_index=True)
        print("Loaded catalog", self.input_data.shape)
        if (self.input_data["Av"] != 0).any():
            raise RuntimeError("catalogue has Av != 0: it must be generated without extinction "
                               "(the simulator applies its own dust once, at the drawn distance)")

    def _prepare(self):
        feh = mh_to_feh(self.input_data["[M/H]"].values, self.input_data["[a/Fe]"].values)

        self.input_columns = np.column_stack([
            np.log10(self.input_data["Teff"].values),
            self.input_data["logg"].values,
            np.zeros(len(self.input_data)),  # AV = 0, always -- extinction is applied
                                               # exactly once, downstream, by
                                               # interpExtinctionAlongSightline in
                                               # Lensing.cpp, at each star's *simulated*
                                               # distance -- not here, at its Besancon
                                               # distance.
            feh,
            self.input_data["[a/Fe]"].values,
        ])

        self.table = self.table.sort_values(["lgTef", "logg", "Av", "Fe_H", "a_Fe"])

        self.grid = (
            np.sort(np.unique(self.table["lgTef"])),
            np.sort(np.unique(self.table["logg"])),
            np.sort(np.unique(self.table["Av"])),
            np.sort(np.unique(self.table["Fe_H"])),
            np.sort(np.unique(self.table["a_Fe"]))
        )

        print("\nGrid dimensions")
        for i, g in enumerate(self.grid):
            print(i, len(g), g.min(), g.max())

        print("\nInput ranges")
        for i in range(5):
            print(i, self.input_columns[:, i].min(), self.input_columns[:, i].max())

        n_expected = np.prod([len(g) for g in self.grid])

        print("\nGrid check")
        print("Rows:", len(self.table))
        print("Expected:", n_expected)

        if len(self.table) != n_expected:
            raise RuntimeError("Table is not a complete regular grid.\n RegularGridInterpolator cannot be used.")

        self.shape = tuple(len(g) for g in self.grid)

    def interp(self):
        self._prepare()

        if self.phot == "lsst":
            self.filter_names = ["LSST_u", "LSST_g", "LSST_r", "LSST_i", "LSST_z", "LSST_y"]

        elif self.phot == "roman":
            self.filter_names = ["Roman_F062", "Roman_F087", "Roman_F106", "Roman_F129", "Roman_F146", "Roman_F158", "Roman_F184", "Roman_F213", "Roman_Grism", "Roman_Prism"]

        else:
            self.filter_names = ["Roman_F146"]

        for name in self.filter_names:
            print(f"Interpolating {name}")

            values = (self.table[name].to_numpy().reshape(self.shape))

            interp = RegularGridInterpolator(self.grid, values, bounds_error=False, fill_value=np.nan)

            self.input_data[name] = self.input_data["Mbol"] - interp(self.input_columns)

        print("NaN values:\n")
        for filt in self.filter_names:
            print(filt, self.input_data[filt].isna().sum())

        print(self.input_data.head())


def fix_low_mass_dwarfs(sub, name, lens_ml):
    """Deviation 82: put the component's unevolved dwarfs on the lens-light M-L relation."""
    if name not in MS_FIX_COMP:
        return sub, 0
    t = lens_ml[lens_ml["comp"] == MS_FIX_COMP[name]]
    centres = ((t["m_lo"] + t["m_hi"]) / 2).to_numpy()
    sel = ((sub["CL"] == 5) & (sub["mass"] < MS_FIX_HI) & (sub["LSST_r"] < DARK_MAG)).to_numpy()
    m = sub["mass"].to_numpy()[sel]
    taper = np.clip((MS_FIX_HI - m) / (MS_FIX_HI - MS_FIX_LO), 0.0, 1.0)
    bins = np.floor(m / MS_FIX_BIN).astype(int)
    bands = {"Roman_F146": "F146", "LSST_u": "u", "LSST_g": "g", "LSST_r": "r", "LSST_i": "i",
             "LSST_z": "z", "LSST_y": "y"}
    for col, short in bands.items():
        v = sub[col].to_numpy()[sel]
        med = pd.Series(v).groupby(bins).median()
        target = np.interp((med.index.to_numpy() + 0.5) * MS_FIX_BIN, centres, t[short].to_numpy())
        shift = pd.Series(target - med.to_numpy(), index=med.index)
        new = v + taper * shift.reindex(bins).to_numpy()
        sub.loc[sub.index[sel], col] = new
    return sub, int(sel.sum())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build():
    rubin = MISTBolometricCorrection("Rubin")
    rubin.read_input("./Besancon/bos9.dat")
    rubin.interp()
    roman = MISTBolometricCorrection("F146")
    roman.input_data = rubin.input_data
    roman.interp()
    df = roman.input_data

    mags = df[FILTER_ORDER].to_numpy()
    no_bc = ~np.isfinite(mags).all(axis=1)
    wd = df["Typ"].between(*DARK_TYP).to_numpy()
    dark = no_bc | wd
    df.loc[dark, FILTER_ORDER] = DARK_MAG
    print(f"dark entries: {int(dark.sum())} of {len(df)} "
          f"(no MIST BC {int(no_bc.sum())}, white dwarfs {int(wd.sum())})")

    out = pd.DataFrame({"mass": df["Mass"], "logT": np.log10(df["Teff"]), "Mbol": df["Mbol"],
                        "Age": df["Age"], "Pop": df["Pop"]})
    for f in ["Roman_F146"] + FILTER_ORDER[:6]:
        out[f] = df[f]
    out["CL"], out["Typ"] = df["CL"], df["Typ"]
    if not np.isfinite(out.to_numpy()).all():
        raise RuntimeError("non-finite value left after the dark-entry substitution")

    lens_ml = pd.read_csv("./components/lens_ml.dat", sep=r"\s+", comment="#", header=None,
                          names=["comp", "m_lo", "m_hi", "n", "u", "g", "r", "i", "z", "y", "F146"])
    rng = np.random.default_rng(SUBSAMPLE_SEED)
    prov = ["# CMD/components provenance (BolometricCorrection.py, Deviation 81)",
            "# component  n_catalogue  n_written  n_dark  mean_mass  median_Mr_lum"]
    for name, pops in COMPONENTS.items():
        sub = out[out["Pop"].isin(pops)].copy()
        sub, n_fix = fix_low_mass_dwarfs(sub, name, lens_ml)
        print(f"{name}: {n_fix} low-mass dwarfs moved onto the lens_ml.dat relation")
        bad = (sub["Age"] > AGE_MAX[name]).sum()
        if bad:
            raise RuntimeError(f"{name}: {bad} stars older than read_cmd's bound {AGE_MAX[name]}")
        n_cat, mmean = len(sub), float(sub["mass"].mean())
        if n_cat > MAX_ROWS:
            sub = sub.iloc[np.sort(rng.choice(n_cat, MAX_ROWS, replace=False))]
        lum = sub["LSST_r"] < DARK_MAG
        sub.to_csv(f"./components/{name}.dat", sep=" ", index=False, float_format="%.4f")
        prov.append(f"{name} {n_cat} {len(sub)} {int((~lum).sum())} {mmean:.4f} "
                    f"{sub.loc[lum, 'LSST_r'].median():.3f}")
        print(f"{name:10s}: catalogue {n_cat:9d}  written {len(sub):9d}  dark {int((~lum).sum()):7d}  "
              f"<m> = {mmean:.4f}")
    open("./components/provenance.txt", "w").write("\n".join(prov) + "\n")
    print("\n".join(prov))
    print("\nBulge.h: N1..N4 = written counts (thin, bulge, thick, halo); MEANMASS_* = mean_mass.")


if __name__ == "__main__":
    build()
