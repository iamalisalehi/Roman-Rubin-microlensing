import glob
import re
import numpy as np
import pandas as pd

from scipy.interpolate import RegularGridInterpolator


def mh_to_feh(mh, afe):
    return mh - np.log10(0.638 * 10**afe + 0.362)


# ---------------------------------------------------------------------------
# Detection thresholds -- auto-parsed from Bulge.h so this can't silently
# drift out of sync with the actual C++ simulation's thre[]/satu[] constants.
# Falls back to MANUAL_THRE/MANUAL_SATU below if Bulge.h can't be found or
# parsed. Order matches Bulge.h and read_cmd exactly: indices 0-5 = LSST
# u,g,r,i,z,y; index 6 = Roman F146.
# ---------------------------------------------------------------------------

BULGE_H_PATH = "../Bulge.h"

FILTER_ORDER = ["LSST_u", "LSST_g", "LSST_r", "LSST_i", "LSST_z", "LSST_y", "Roman_F146"]

# Fallback values -- only used if Bulge.h can't be read/parsed. Fill these in
# with your actual current thre[]/satu[] arrays from Bulge.h if you'd rather
# not rely on auto-parsing, e.g.:
# MANUAL_THRE = [23.4, 24.6, 24.3, 23.6, 22.9, 21.7, 23.0]
# MANUAL_SATU = [15.2, 16.3, 16.0, 15.3, 14.6, 13.4, 12.0]
MANUAL_THRE = None
MANUAL_SATU = None

# Per-component age ceilings, mirroring read_cmd()'s CHECK(age_X[j] <= ...)
# bounds exactly. Keep these in sync with helper.cpp if those ever change.
AGE_MAX = {"thin_disk": 10, "bulge": 10, "thick_disk": 13, "halo": 14}


def parse_array_from_bulge_h(path, name):
    """Pulls a `constexpr std::array<double, M> NAME = {...};`-style literal
    directly out of Bulge.h. Returns None (triggering the MANUAL_* fallback)
    if the file or the named array can't be found."""
    try:
        text = open(path, "r").read()
    except OSError as e:
        print(f"[warn] could not read {path} ({e})")
        return None

    m = re.search(rf"\b{name}\s*=\s*\{{([^}}]*)\}}", text)
    if not m:
        print(f"[warn] could not find '{name}' array in {path}")
        return None

    return [float(v) for v in m.group(1).split(",") if v.strip()]


def get_thresholds():
    thre_list = parse_array_from_bulge_h(BULGE_H_PATH, "thre") or MANUAL_THRE
    satu_list = parse_array_from_bulge_h(BULGE_H_PATH, "satu") or MANUAL_SATU

    if thre_list is None or satu_list is None:
        raise RuntimeError(
            "Could not obtain thre[]/satu[] -- either point BULGE_H_PATH at your "
            "real Bulge.h, or fill in MANUAL_THRE/MANUAL_SATU by hand at the top "
            "of this file."
        )
    if len(thre_list) != 7 or len(satu_list) != 7:
        raise RuntimeError(
            f"Expected 7 values each (6 LSST + F146), got thre={len(thre_list)}, "
            f"satu={len(satu_list)} -- check the arrays in Bulge.h match the "
            f"current M=7 filter set."
        )

    thre = dict(zip(FILTER_ORDER, thre_list))
    satu = dict(zip(FILTER_ORDER, satu_list))
    print(f"Using detection thresholds:\n  thre = {thre}\n  satu = {satu}")
    return thre, satu


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

        if "Dist" not in self.input_data.columns:
            raise RuntimeError(
                "No 'Dist' column found in the Besancon catalog header -- the "
                "visibility filter needs each star's own Besancon-assigned distance "
                "to estimate an apparent magnitude. Check the exact column name in "
                "your catalog's header line."
            )

        feh = mh_to_feh(self.input_data["[M/H]"], self.input_data["[a/Fe]"])

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


def apply_visibility_filter(df, thre, satu):
    """Keeps only stars detectable in at least F146 AND one LSST filter.

    Apparent magnitude is estimated from each star's own Besancon-assigned
    Dist (kpc) via the same distance-modulus formula Lensing.cpp uses
    (Map = Mab + 5*log10(Ds*100)) -- just using Besancon's own Dist in place
    of the C++ simulation's later, independently re-drawn Ds, since that's
    the only distance information available at this stage of the pipeline.

    Extinction is intentionally left out here, for the same reason it's
    zeroed in _prepare(): it's applied exactly once, downstream, by the C++
    code's 3D dust map. That makes this filter deliberately optimistic --
    it keeps a star if it could plausibly be seen unreddened -- rather than
    risking throwing away something that would only look invisible because
    of a cut made before real extinction is even known. Anything that's
    only marginally detectable will still have to pass the simulation's own
    per-event detectability check later; this filter only exists to drop
    stars that could never be seen under any circumstance.
    """
    dist_kpc = df["Dist"].to_numpy()
    n_bad_dist = int((dist_kpc <= 0).sum())
    if n_bad_dist:
        print(f"[warn] {n_bad_dist} rows have Dist <= 0 kpc -- treating as not visible in any filter.")

    dist_mod = np.where(dist_kpc > 0, 5.0 * np.log10(np.clip(dist_kpc, 1e-12, None) * 100.0), np.inf)

    visible = {}
    for filt in thre:
        app_mag = df[filt].to_numpy() + dist_mod
        visible[filt] = (app_mag > satu[filt]) & (app_mag <= thre[filt])

    f146_ok = visible["Roman_F146"]
    lsst_ok = np.zeros(len(df), dtype=bool)
    for filt in ["LSST_u", "LSST_g", "LSST_r", "LSST_i", "LSST_z", "LSST_y"]:
        lsst_ok |= visible[filt]

    keep = f146_ok & lsst_ok
    n_before = len(df)
    n_after = int(keep.sum())
    print(f"\nVisibility filter (F146 AND >=1 LSST filter): "
          f"kept {n_after}/{n_before} ({100*n_after/max(n_before,1):.2f}%), "
          f"dropped {n_before - n_after}.")

    return df[keep].copy()

# ---------------------------------------------------------------------------
# Stellar type/class filtering.
#
# Besancon's CL (luminosity class) and Typ (spectral type) codes are
# specific to your query -- I don't have a verified, universal numeric
# mapping to hand you, and I'd rather you check your own catalog's
# documentation than have me guess at what excludes what. Check either
# the header/notes for your bos9.dat query, or the model's own docs at
# https://model.obs-besancon.fr/, to find which CL/Typ values correspond
# to the star types you don't want, then fill these in:
# ---------------------------------------------------------------------------

# CL values to drop entirely, e.g. [7] to drop one specific luminosity
# class code. Leave as [] to skip CL-based filtering.
EXCLUDE_CL = []

# Typ ranges to drop, as a list of (min, max) tuples (inclusive), e.g.
# [(7.0, 9.0)] to drop everything with Typ between 7 and 9.
# Leave as [] to skip Typ-based filtering.
EXCLUDE_TYP_RANGES = [(9.0, 9.2)]


def apply_type_filter(df, exclude_cl=None, exclude_typ_ranges=None):
    """Drops stars matching unwanted CL (luminosity class) and/or Typ
    (spectral type) codes. Prints the actual CL/Typ distribution present
    in the catalog first, so you can sanity-check EXCLUDE_CL/
    EXCLUDE_TYP_RANGES against what's really there, then reports exactly
    how many rows each exclusion rule drops."""
    exclude_cl = exclude_cl or []
    exclude_typ_ranges = exclude_typ_ranges or []

    print(f"\n{'=' * 60}\nType/class filter\n{'=' * 60}")
    print("CL value counts in catalog:")
    print(df["CL"].value_counts().sort_index().to_string())
    print(f"\nTyp distribution: min={df['Typ'].min():.2f} max={df['Typ'].max():.2f}")

    keep = np.ones(len(df), dtype=bool)

    if exclude_cl:
        cl_bad = df["CL"].isin(exclude_cl)
        print(f"\nCL in {exclude_cl}: dropping {int(cl_bad.sum())} rows")
        keep &= ~cl_bad

    for lo, hi in exclude_typ_ranges:
        typ_bad = df["Typ"].between(lo, hi)
        print(f"Typ in [{lo}, {hi}]: dropping {int(typ_bad.sum())} rows")
        keep &= ~typ_bad

    n_before = len(df)
    n_after = int(keep.sum())
    print(f"\nType/class filter total: kept {n_after}/{n_before} "
          f"({100 * n_after / max(n_before, 1):.2f}%), dropped {n_before - n_after}.")

    return df[keep].copy()

def save_components(df, filter_names):
    """Splits by population, applies each component's own age ceiling
    (mirroring read_cmd()'s CHECK(age_X[j] <= ...) exactly), and writes the
    four CMD files. Reports before/after row counts at each step so you can
    see exactly what N1..N4 need to become in Bulge.h."""
    out = pd.DataFrame()
    out['mass'] = df['Mass']
    out['logT'] = np.log10(df['Teff'])
    out['Mbol'] = df['Mbol']
    out['Age'] = df['Age']
    out['Pop'] = df['Pop']
    out[filter_names] = df[filter_names]
    out['CL'] = df['CL']
    out['Typ'] = df['Typ']

    n_before_dropna = len(out)
    out = out.dropna()
    print(f"\ndropna() (catches any remaining interpolation failures): "
          f"kept {len(out)}/{n_before_dropna}")

    components = {
        "thin_disk":  out[out["Pop"].between(1, 7)].copy(),
        "thick_disk": out[out["Pop"].isin([8, 11])].copy(),
        "halo":       out[out["Pop"] == 9].copy(),
        "bulge":      out[out["Pop"] == 10].copy(),
    }

    print(f"\n{'=' * 60}\nPer-component age filter\n{'=' * 60}")
    new_N = {}
    for name, sub in components.items():
        age_max = AGE_MAX[name]
        n_before = len(sub)
        sub_filtered = sub[sub["Age"] <= age_max].copy()
        n_after = len(sub_filtered)
        print(f"  {name:12s}: Age <= {age_max:2d}  ->  kept {n_after}/{n_before} "
              f"(dropped {n_before - n_after})")
        components[name] = sub_filtered
        new_N[name] = n_after

    components["thin_disk"].to_csv("./components/thin_disk.dat", sep=" ", index=False, float_format="%.4f")
    components["thick_disk"].to_csv("./components/thick_disk.dat", sep=" ", index=False, float_format="%.4f")
    components["halo"].to_csv("./components/halo.dat", sep=" ", index=False, float_format="%.4f")
    components["bulge"].to_csv("./components/bulge.dat", sep=" ", index=False, float_format="%.4f")

    print(f"\n{'=' * 60}\nFinal row counts -- update Bulge.h to match\n{'=' * 60}")
    print(f"  N1 (thin_disk)  = {new_N['thin_disk']}")
    print(f"  N2 (bulge)      = {new_N['bulge']}")
    print(f"  N3 (thick_disk) = {new_N['thick_disk']}")
    print(f"  N4 (halo)       = {new_N['halo']}")

    return new_N


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

thre, satu = get_thresholds()

rubin = MISTBolometricCorrection("Rubin")
rubin.read_input("./Besancon/bos9.dat")
rubin.interp()

roman = MISTBolometricCorrection("F146")
roman.input_data = rubin.input_data
roman.interp()
for name in rubin.filter_names:
    roman.filter_names.append(name)

n_total = len(roman.input_data)
print(f"\nTotal stars before any visibility/age filtering: {n_total}")

visible_data = apply_visibility_filter(roman.input_data, thre, satu)
typed_data = apply_type_filter(visible_data, EXCLUDE_CL, EXCLUDE_TYP_RANGES)
new_N = save_components(typed_data, roman.filter_names)

