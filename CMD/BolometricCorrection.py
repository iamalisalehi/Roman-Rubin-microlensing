import glob
import numpy as np
import pandas as pd

from scipy.interpolate import RegularGridInterpolator


def mh_to_feh(mh, afe):
    return mh - np.log10(0.638 * 10**afe + 0.362)


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
        with open(path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    self.header_input = line.lstrip("#").split()
                else:
                    break

        self.input_data = pd.DataFrame(np.loadtxt(path), columns=self.header_input)

        print("Loaded catalog")
        print(self.input_data.shape)
        feh = mh_to_feh(self.input_data["[M/H]"], self.input_data["[a/Fe]"])

        # Note: these Av/FeH/aFe diagnostics describe Besancon's own catalog values.
        # Since Av is now forced to 0 before interpolation (see _prepare below),
        # "Av > 6" no longer predicts which rows get dropped as NaN -- it's
        # informational only, showing how much extinction Besancon itself assigned.
        print("Av > 6      :", np.sum(self.input_data["Av"] > 6))
        print("FeH > 0.5   :", np.sum(feh > 0.5))
        print("FeH < -3.0  :", np.sum(feh < -3.0))
        print("aFe > 0.6   :", np.sum(self.input_data["[a/Fe]"] > 0.6))
        print("aFe < -0.2  :", np.sum(self.input_data["[a/Fe]"] < -0.2))

    def _prepare(self):
        feh = mh_to_feh(self.input_data["[M/H]"].values, self.input_data["[a/Fe]"].values)

        self.input_columns = np.column_stack([
            np.log10(self.input_data["Teff"].values),
            self.input_data["logg"].values,
            np.zeros(len(self.input_data)),  # AV = 0, always. Extinction is applied
                                               # exactly once, downstream, by
                                               # interpExtinctionAlongSightline in
                                               # Lensing.cpp, at each star's *simulated*
                                               # distance -- not here, at its Besancon
                                               # distance. Feeding Besancon's own Av
                                               # into this grid was the double-extinction
                                               # / NaN root cause.
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
            print(filt,self.input_data[filt].isna().sum())
        
        print(self.input_data.head())

    def save_components(self, breakdown_components=False):
        df = pd.DataFrame()
        df['mass'] = self.input_data['Mass']
        df['logT'] = np.log10(self.input_data['Teff'])
        df['Mbol'] = self.input_data['Mbol']
        df['Age'] = self.input_data['Age']
        df['Pop'] = self.input_data['Pop']
        df[self.filter_names] = self.input_data[self.filter_names]
        df['CL'] = self.input_data['CL']
        df['Typ'] = self.input_data['Typ']
        df = df.dropna()

        if not breakdown_components:
            df.to_csv("regular_output.dat", sep=" ", index=False, float_format="%.4f")
        else:
            thin_disk = df[df["Pop"].between(1, 7)].copy()
            thick_disk = df[df["Pop"].isin([8, 11])].copy()
            halo = df[df["Pop"] == 9]
            bulge = df[df["Pop"] == 10]

            thin_disk.to_csv("./components/thin_disk.dat", sep=" ", index=False, float_format="%.4f")
            thick_disk.to_csv("./components/thick_disk.dat", sep=" ", index=False, float_format="%.4f")
            halo.to_csv("./components/halo.dat", sep=" ", index=False, float_format="%.4f")
            bulge.to_csv("./components/bulge.dat", sep=" ", index=False, float_format="%.4f")


rubin = MISTBolometricCorrection("Rubin")
rubin.read_input("./Besancon/bos9.dat")
rubin.interp()

roman = MISTBolometricCorrection("F146")
roman.input_data = rubin.input_data
roman.interp()
for name in rubin.filter_names:
    roman.filter_names.append(name)

roman.save_components(True)
