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

    rng = np.random.default_rng(SUBSAMPLE_SEED)
    prov = ["# CMD/components provenance (BolometricCorrection.py, Deviation 81)",
            "# component  n_catalogue  n_written  n_dark  mean_mass  median_Mr_lum"]
    for name, pops in COMPONENTS.items():
        sub = out[out["Pop"].isin(pops)]
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
