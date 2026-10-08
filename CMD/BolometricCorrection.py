"""Builds the source/neighbour catalogues CMD/components/{thin_disk,bulge,thick_disk,halo}.dat.

Each file is a sample of every star of one Galactic component in a Besancon catalogue (one
0.05-0.1 deg^2 field at (l, b) = (0.5, -1.4), Av = 0), with absolute AB magnitudes from MIST
bolometric corrections. The simulator draws a source and each of its blend neighbours uniformly from
these lists, with a distance from its own density model, and counts them with Nstart = rho / <m>, so
the lists must be the complete population Nstart counts. There is no visibility filter: visibility
is decided per event in the simulator, with magnification, dust and each survey's own limits.

Dark entries: a star MIST cannot place on its grid (white dwarfs, Besancon Typ 9.0-9.2, plus a few
NaNs) is kept with every magnitude = DARK_MAG. It is part of what Nstart counts but contributes no
light.

Mean masses: <m> per component over the whole list (dark entries included) is written to
<out-dir>/provenance.txt, over the written list and over the full pre-subsample population.
config/data_products.h's MEANMASS_* (generated from provenance.txt by tools/sync_data_products.py)
must equal the population value, so that Nstart counts the population the draws come from.

Usage (paths are relative to CMD/; the script changes into it):

    cd CMD && ../.roman/bin/python BolometricCorrection.py [options]

  --input PATH       the Besancon catalogue (default Besancon/bos10). Rows whose field count differs
                     from the header's are skipped and counted (the file is streamed through
                     `awk 'NF == <header fields>'`); the count and line numbers go to provenance.txt.
  --out-dir DIR      where the four lists, their *_aux.npz files and provenance.txt go. Default
                     components_staging/<variant>, <variant> = <input>_<dwarfs>_fill<floor>.
                     Writing the production lists needs the explicit `--out-dir components`.
  --dwarfs MODE      what is done to the low-mass dwarfs (CL = 5, Typ < 9):
        besancon     nothing: every star keeps Besancon's Teff, Mbol and logg.
        empirical    (default) the thin-disc and bulge dwarfs below 0.7 Msun are moved onto the
                     Pecaut & Mamajek (2013) dwarf sequence (CMD/empirical/, Mamajek's table version
                     2022.04.16; the 41 rows with 0.075 <= Msun <= 1.0). See below.
        dev82        the bulge, thick-disc and halo dwarfs are moved, in every band, onto the
                     lens_ml.dat relation. Needs a lens_ml.dat (components/, or --lens-ml).
  --fill-floor MODE  none (default), kroupa (alpha = 1.3) or koshimoto (alpha = 1.16): add synthetic
                     dwarfs below each component's mass floor. See below.
  --no-grid-clamp    switch the BC-grid clamp off (stars outside the MIST grid become dark).
  --lens-ml PATH     the lens_ml.dat that --dwarfs dev82 reads (default components/lens_ml.dat).

Empirical dwarf sequence. Besancon's thin-disc M dwarfs are 0.1-0.3 mag too faint at a given mass,
and its bulge/thick-disc/halo dwarfs 2-3 mag too bright, against real M dwarfs (Mann et al. 2019
M_K-mass, Benedict et al. 2016 M_V-mass; thick disc and halo pass, thin disc and bulge do not).
bos10's masses, Teff, Mbol, logg, [M/H], [a/Fe] and ages carry no injected noise, so the correction is
made in physical space at the star's true mass m. For the thin disc (Pop 1-7) and the bulge (Pop 10),
each dwarf with m < 0.7 Msun is shifted in Mbol and in log10 Teff by

    taper(m) * [ X_E(m) - median_c X(m) ],      X = Mbol, log10 Teff,
    taper(m) = clip((0.7 - m) / 0.1, 0, 1),

where median_c X(m) is the component's own median of X in 0.01-Msun true-mass bins from 0.07 to 0.80
Msun (all of the component's dwarfs in the input, before subsampling; bins with >= 20 stars; empty
bins interpolated over, flat beyond the populated range) and X_E(m) is the empirical sequence at that
mass (Mbol_E = 4.74 - 2.5 logL_E, log10 Teff_E, interpolated linearly in log10 m). The star keeps its
scatter about the median, which is Besancon's age/metallicity spread. logg of the shifted stars is
recomputed: logg = 4.438 + log10 m + 4 log10(Teff / 5772) + 0.4 (Mbol - 4.74). The bolometric
corrections are taken afterwards. The thick disc (Pop 8, 11) and halo (Pop 9) are unchanged.

Filling the mass floor. Besancon's bulge, thick disc and halo have no stars below a floor m_f (the
component's minimum dwarf mass, ~0.154-0.159 Msun), though the thin disc reaches 0.073. With
--fill-floor, N_add synthetic dwarfs per component are drawn between 0.08 Msun and m_f from a power law
dN/dm = n0 (m / m_mid)^-alpha, continuous with the catalogue at the floor: n0 is the number of the
component's dwarfs with m in [m_f, m_f + 0.03) divided by 0.03, at m_mid = m_f + 0.015; N_add =
round(integral). Windows are half-open with a 1e-6 tolerance because catalogue masses are quantised to
0.001 Msun. Each synthetic star copies Pop, Age, [M/H], [a/Fe] from a random real dwarf of the same
component with m in [m_f, m_f + 0.05), has CL = 5, takes Typ from a random thin-disc dwarf within
+-0.005 Msun of its mass, and gets the empirical Mbol, Teff, logg. Draws use seed FILL_SEED, in the
order bulge, thick disc, halo; per component: masses, donors, Typ.

BC-grid clamp (always on unless --no-grid-clamp; BC lookup only, Besancon's [M/H], [a/Fe] are what is
stored). The MIST grid spans [Fe/H] = -3.0..+0.5 and [a/Fe] = -0.2..+0.6; outside it
RegularGridInterpolator returns NaN. [a/Fe] is clipped to the grid; [M/H] is then shifted so that
[Fe/H] = [M/H] - log10(0.638 10^[a/Fe] + 0.362) lies inside the grid (1e-6 from the edges). The counts
of what moved are in provenance.txt.

Subsample. A component with more than MAX_ROWS stars (synthetic ones included) is reduced to a uniform
random subsample (seed SUBSAMPLE_SEED), since the simulator holds the lists in RAM. It is drawn before
the BCs are computed (a star's BC does not depend on the others), except with --dwarfs dev82, whose
running medians need every star.

Output in <out-dir>: thin_disk.dat bulge.dat thick_disk.dat halo.dat (header
`mass logT Mbol Age Pop Roman_F146 LSST_u LSST_g LSST_r LSST_i LSST_z LSST_y CL Typ`, %.4f), one
<name>_aux.npz per component, row-aligned to the written list (mass, Teff_used, Mbol_used, logg_used,
M_H and a_Fe as Besancon has them, Pop, Age, CL, Typ, synthetic, shifted, dark), and provenance.txt.

dev82 shift. Below MS_FIX_HI each bulge, thick-disc and halo dwarf (CL = 5) is moved, in every band,
by the difference between components/lens_ml.dat and its component's own running median at that mass;
the star's scatter about the median is kept. The shift tapers linearly to zero between MS_FIX_LO and
MS_FIX_HI. The thin disc is the reference and is not touched. Run lens_ml_table.py first.

Other code (analysis/b4, b6, lens_ml_table.py) imports the part of this file above the "# Main" marker
without building anything: keep that part free of module-level work and of __file__.
"""
import argparse
import contextlib
import glob
import io
import os
import resource
import subprocess
import tempfile
import time

import numpy as np
import pandas as pd

from scipy.interpolate import RegularGridInterpolator


def mh_to_feh(mh, afe):
    return mh - np.log10(0.638 * 10**afe + 0.362)


INPUT_COLUMNS = ["Teff", "logg", "Pop", "Age", "Mass", "Mbol", "[M/H]", "[a/Fe]", "CL", "Typ", "Av"]
FILTER_ORDER = ["LSST_u", "LSST_g", "LSST_r", "LSST_i", "LSST_z", "LSST_y", "Roman_F146"]
DARK_MAG = 99.0                       # "no light": 10^(-0.4*99) is zero in every sum
DARK_TYP = (9.0, 9.2)                 # Besancon white dwarfs: no MIST track, faint -> dark
MAX_ROWS = 3_500_000                  # per component; only the bulge exceeds it
SUBSAMPLE_SEED = 20261002
FILL_SEED = 20261005
# Besancon Pop codes -> component, in the order of include/common.h's GalacticComponent.
COMPONENTS = {"thin_disk": list(range(1, 8)), "bulge": [10], "thick_disk": [8, 11], "halo": [9]}
MS_FIX_LO, MS_FIX_HI = 0.6, 0.7       # Msun: full shift below LO, none above HI
MS_FIX_COMP = {"bulge": 1, "thick_disk": 2, "halo": 3}   # lens_ml.dat comp codes
MS_FIX_BIN = 0.02                     # Msun, running-median bin
# Upper age bounds read_cmd() CHECKs (src/galaxy/catalogue.cpp); a violation stops the build, it does not drop.
AGE_MAX = {"thin_disk": 10, "bulge": 10, "thick_disk": 13, "halo": 14}

STAGING_ROOT = "components_staging"
DEFAULT_INPUT = "Besancon/bos10"
EMPIRICAL_TABLE = "empirical/EEM_dwarf_UBVIJHK_colors_Teff.txt"
EMPIRICAL_ROWS = 41                   # rows of the table with finite Msun, Teff, logL and 0.075 <= Msun <= 1.0
EMPIRICAL_MASS = (0.075, 1.0)
EMPIRICAL_COMPONENTS = ("thin_disk", "bulge")      # option E moves these two only
FILL_COMPONENTS = ("bulge", "thick_disk", "halo")  # the ones with a mass floor
FILL_ALPHA = {"kroupa": 1.3, "koshimoto": 1.16}
FILL_MASS_LO = 0.08                   # Msun: synthetic dwarfs start at the hydrogen-burning limit
FILL_NORM_WIDTH = 0.03                # Msun: window above the floor that sets n0
FILL_DONOR_WIDTH = 0.05               # Msun: window above the floor the attributes are copied from
FILL_TYP_WIDTH = 0.005                # Msun: Typ is copied from a thin-disc dwarf this close in mass
MASS_EPS = 1.0e-6                     # catalogue masses are quantised to 0.001 Msun; window-edge tolerance
OFFSET_EDGES = np.round(np.arange(0.07, 0.80 + 1e-9, 0.01), 4)   # 0.01-Msun bins of the median curve
OFFSET_NMIN = 20                      # stars a bin needs for its medians to count
OFFSET_REPORT = (0.16, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65)            # masses quoted in provenance
LOG_TEFF_SUN = np.log10(5772.0)
# The MIST grid and the clamp of (Besancon's [M/H], [a/Fe]) onto it.
AFE_LO, AFE_HI, FEH_LO, FEH_HI, GRID_EPS = -0.2, 0.6, -3.0, 0.5, 1.0e-6

# What the catalogue reader keeps, and how (Mass stays float64: window edges are decided on it).
CATALOGUE_DTYPES = {"Teff": np.float32, "logg": np.float32, "Pop": np.int8, "Age": np.float32,
                    "Mass": np.float64, "Mbol": np.float32, "[M/H]": np.float32, "[a/Fe]": np.float32,
                    "CL": np.int8, "Typ": np.float32, "Dist": np.float32}
CATALOGUE_KEYS = {"[M/H]": "M_H", "[a/Fe]": "a_Fe"}      # column name -> key in the dict read_catalogue returns
DEFAULT_COLUMNS = ("Teff", "logg", "Pop", "Age", "Mass", "Mbol", "[M/H]", "[a/Fe]", "CL", "Typ")


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
        """Reads only the columns this pipeline uses (the catalogue is multi-GB), skipping damaged rows."""
        self.input_data = pd.concat(list(iter_catalogue_chunks(path, INPUT_COLUMNS, 2_000_000)),
                                    ignore_index=True)
        print("Loaded catalog", self.input_data.shape)
        if (self.input_data["Av"] != 0).any():
            raise RuntimeError("catalogue has Av != 0: it must be generated without extinction "
                               "(the simulator applies its own dust once, at the drawn distance)")

    def _ensure_grid(self):
        """Sort the table onto its regular grid once, check it is complete, remember the axes."""
        if getattr(self, "_grid_ready", False):
            return
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

        n_expected = np.prod([len(g) for g in self.grid])

        print("\nGrid check")
        print("Rows:", len(self.table))
        print("Expected:", n_expected)

        if len(self.table) != n_expected:
            raise RuntimeError("Table is not a complete regular grid.\n RegularGridInterpolator cannot be used.")

        self.shape = tuple(len(g) for g in self.grid)
        self._interps = {}
        self._grid_ready = True

    def _prepare(self):
        feh = mh_to_feh(self.input_data["[M/H]"].values, self.input_data["[a/Fe]"].values)

        self.input_columns = np.column_stack([
            np.log10(self.input_data["Teff"].values),
            self.input_data["logg"].values,
            np.zeros(len(self.input_data)),  # AV = 0: the simulator applies extinction once, at the
                                               # star's simulated distance (src/galaxy/extinction.cpp)
            feh,
            self.input_data["[a/Fe]"].values,
        ])

        self._ensure_grid()

        print("\nInput ranges")
        for i in range(5):
            print(i, self.input_columns[:, i].min(), self.input_columns[:, i].max())

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

    def bc(self, name, lgTef, logg, feh, afe):
        """Bolometric correction of band `name` (Mbol - M_band) at AV = 0; NaN off the grid.

        The interpolator of each band is built once and kept, so a build can feed it in chunks."""
        self._ensure_grid()
        f = self._interps.get(name)
        if f is None:
            f = RegularGridInterpolator(self.grid, self.table[name].to_numpy().reshape(self.shape),
                                        bounds_error=False, fill_value=np.nan)
            self._interps[name] = f
        return f(np.column_stack([lgTef, logg, np.zeros(len(lgTef)), feh, afe]))


# ---------------------------------------------------------------------------
# Reading the catalogue (the damaged rows, the typed columns)
# ---------------------------------------------------------------------------

def variant_name(input_path, dwarfs, fill_floor):
    """<input stem>_<dwarfs>_fill<floor>, e.g. bos10_empirical_fillnone: names a staging directory."""
    stem = os.path.splitext(os.path.basename(os.path.normpath(input_path)))[0]
    return f"{stem}_{dwarfs}_fill{fill_floor}"


def iter_catalogue_chunks(path, usecols=None, chunksize=1_000_000, report=None):
    """Yields DataFrame chunks of the catalogue's GOOD rows; skipped ones are counted in `report`.

    A good row has as many fields as the header. The file goes through awk, which drops the others
    (pandas would shift their columns or pad them with NaN) and notes their file line numbers (the
    header is line 1). report["n_skipped"] and report["bad_lines"] are filled once the generator is
    exhausted; an awk failure raises."""
    with open(path, "r") as f:
        header = f.readline().lstrip("#").split()
    nf = len(header)
    bad = tempfile.TemporaryFile(mode="w+")
    prog = 'NR > 1 { if (NF == %d) print; else print NR > "/dev/stderr" }' % nf
    proc = subprocess.Popen(["awk", prog, path], stdout=subprocess.PIPE, stderr=bad,
                            env={**os.environ, "LC_ALL": "C"})
    done = False
    try:
        for chunk in pd.read_csv(proc.stdout, sep=r"\s+", header=None, names=header, usecols=usecols,
                                 chunksize=chunksize):
            yield chunk
        done = True
    finally:
        proc.stdout.close()
        if not done:
            proc.kill()
        rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"awk failed on {path} (exit {rc})")
    bad.seek(0)
    lines = [int(x) for x in bad.read().split()]
    bad.close()
    if report is not None:
        report["n_skipped"], report["bad_lines"], report["n_fields"] = len(lines), lines, nf
    print(f"{path}: skipped {len(lines)} rows with a field count != {nf}"
          + (f" (file lines {', '.join(map(str, lines))})" if lines else ""))


def read_catalogue(path, columns=DEFAULT_COLUMNS):
    """The catalogue's good rows as typed NumPy arrays: {key: array}, plus a report.

    Keys are the column names, except [M/H] -> M_H and [a/Fe] -> a_Fe. Float32 where that loses nothing
    (Teff, logg, Mbol, Age, Typ, the metallicities are printed with <= 4 decimals), Mass in float64
    (window edges are decided on it), Pop and CL in int8. Av must be 0 everywhere (the simulator applies
    its own dust, once); a non-finite value, or a non-integral Pop or CL, raises."""
    columns = list(columns)
    parts = {c: [] for c in columns}
    report = {}
    av_max, n = 0.0, 0
    for chunk in iter_catalogue_chunks(path, columns + ["Av"], 1_000_000, report):
        av_max = max(av_max, float(np.abs(chunk["Av"].to_numpy()).max()))
        for c in columns:
            v = chunk[c].to_numpy()
            if not np.isfinite(v).all():
                raise RuntimeError(f"{path}: non-finite value in column {c}")
            dt = CATALOGUE_DTYPES[c]
            if np.issubdtype(dt, np.integer) and not (v == np.round(v)).all():
                raise RuntimeError(f"{path}: non-integral value in column {c}")
            parts[c].append(v.astype(dt))
        n += len(chunk)
    if av_max != 0.0:
        raise RuntimeError("catalogue has Av != 0: it must be generated without extinction "
                           "(the simulator applies its own dust once, at the drawn distance)")
    cat = {}
    for c in columns:
        cat[CATALOGUE_KEYS.get(c, c)] = np.concatenate(parts[c])
        parts[c] = None
    report["n_rows"] = n
    print(f"read {n} rows of {path}")
    return cat, report


def component_codes(pop):
    """int8 array: index into COMPONENTS (thin, bulge, thick, halo); raises on a Pop code it does not know."""
    lut = np.full(256, -1, np.int8)
    for k, pops in enumerate(COMPONENTS.values()):
        lut[pops] = k
    comp = lut[pop.astype(np.int64)]
    if (comp < 0).any():
        raise RuntimeError(f"{int((comp < 0).sum())} stars with a Pop code outside {sorted(sum(COMPONENTS.values(), []))}: "
                           f"{np.unique(pop[comp < 0])}")
    return comp


# ---------------------------------------------------------------------------
# The empirical dwarf sequence and the component's own median curve
# ---------------------------------------------------------------------------

def qmass(m):
    """Mass as float64 rounded to 1e-6 Msun, so that 0.186 compares equal to 0.186 whatever its dtype."""
    return np.round(np.asarray(m, dtype=np.float64), 6)


def taper(m):
    """1 below MS_FIX_LO, 0 above MS_FIX_HI, linear between."""
    return np.clip((MS_FIX_HI - m) / (MS_FIX_HI - MS_FIX_LO), 0.0, 1.0)


def logg_formula(mass, logT, mbol):
    """logg from mass, log10 Teff and Mbol: 4.438 + log10 m + 4 log10(Teff/5772) + 0.4 (Mbol - 4.74)."""
    return 4.438 + np.log10(mass) + 4.0 * (logT - LOG_TEFF_SUN) + 0.4 * (mbol - 4.74)


class EmpiricalDwarfs:
    """Pecaut & Mamajek (2013) mean dwarf sequence: Mbol_E(m) and log10 Teff_E(m).

    Parsed from Mamajek's table (version 2022.04.16): whitespace columns, dots for missing values,
    comment lines starting with '#', and free-text notes after the table, where the parsing stops.
    The rows kept have finite Msun, Teff and logL with 0.075 <= Msun <= 1.0 (41 of them); both
    quantities are interpolated linearly in log10(mass) and held flat beyond the first/last row."""

    def __init__(self, path=EMPIRICAL_TABLE):
        names, rows = None, []
        with open(path, "r") as f:
            for line in f:
                s = line.strip()
                if names is None:
                    if s.startswith("#SpT"):
                        names = s.lstrip("#").split()
                    continue
                if not s or s.startswith("#"):
                    break                       # the table ends; what follows is notes
                rows.append(s.split())
        if names is None or not rows:
            raise RuntimeError(f"{path}: no table found")
        for r in rows:
            if len(r) != len(names):
                raise RuntimeError(f"{path}: row {r[0]} has {len(r)} fields, header has {len(names)}")

        def col(name):
            out = np.full(len(rows), np.nan)
            j = names.index(name)
            for i, r in enumerate(rows):
                try:
                    out[i] = float(r[j])
                except ValueError:
                    pass                        # '...', '.....', '19.25:' -> missing
            return out

        mass, teff, logl = col("Msun"), col("Teff"), col("logL")
        keep = (np.isfinite(mass) & np.isfinite(teff) & np.isfinite(logl)
                & (mass >= EMPIRICAL_MASS[0]) & (mass <= EMPIRICAL_MASS[1]))
        if int(keep.sum()) != EMPIRICAL_ROWS:
            raise RuntimeError(f"{path}: {int(keep.sum())} usable rows, expected {EMPIRICAL_ROWS}")
        order = np.argsort(mass[keep], kind="stable")
        self.mass, self.teff, self.logL = mass[keep][order], teff[keep][order], logl[keep][order]
        self.spt = [rows[i][0] for i in np.flatnonzero(keep)[order]]
        if not (np.diff(self.mass) > 0).all():
            raise RuntimeError(f"{path}: the usable masses are not strictly increasing")
        self._x, self._lt = np.log10(self.mass), np.log10(self.teff)
        self.path, self.version = path, "Mamajek 2022.04.16 (Pecaut & Mamajek 2013)"

    def logT(self, m):
        return np.interp(np.log10(m), self._x, self._lt)

    def mbol(self, m):
        return 4.74 - 2.5 * np.interp(np.log10(m), self._x, self.logL)


def median_curve(mass, mbol, logT):
    """A component's own median Mbol and log10 Teff per 0.01-Msun true-mass bin, 0.07-0.80 Msun.

    Bins with fewer than OFFSET_NMIN stars are interpolated over, and the curve is flat beyond the
    populated range. Returns {centres, n, ok, mbol, logT}."""
    edges = OFFSET_EDGES
    nb = len(edges) - 1
    ib = np.digitize(qmass(mass), edges) - 1
    inr = (ib >= 0) & (ib < nb)
    g = pd.DataFrame({"b": ib[inr], "Mbol": np.asarray(mbol, np.float64)[inr],
                      "lT": np.asarray(logT, np.float64)[inr]}).groupby("b").agg(
        n=("Mbol", "size"), Mbol=("Mbol", "median"), lT=("lT", "median")).reindex(range(nb))
    n = g["n"].fillna(0).to_numpy().astype(int)
    ok = n >= OFFSET_NMIN
    if int(ok.sum()) < 2:
        raise RuntimeError("median curve: fewer than 2 populated mass bins")
    centres = 0.5 * (edges[:-1] + edges[1:])
    return {"centres": centres, "n": n, "ok": ok,
            "mbol": np.interp(centres, centres[ok], g["Mbol"].to_numpy()[ok]),
            "logT": np.interp(centres, centres[ok], g["lT"].to_numpy()[ok])}


def curve_at(curve, key, m):
    return np.interp(m, curve["centres"], curve[key])


# ---------------------------------------------------------------------------
# The MIST grid clamp and the seven simulator bands
# ---------------------------------------------------------------------------

def clamp_to_grid(mh, afe):
    """(mh, afe) put onto the MIST grid for the BC lookup: [a/Fe] clipped, [M/H] shifted until [Fe/H] fits."""
    afe_c = np.clip(afe, AFE_LO, AFE_HI)
    feh = mh_to_feh(mh, afe_c)
    feh_c = np.clip(feh, FEH_LO + GRID_EPS, FEH_HI - GRID_EPS)
    return mh + (feh_c - feh), afe_c


def clamp_counts(mh, afe, wd=None, chunk=2_000_000):
    """How many stars the clamp moves: [a/Fe] clipped, [Fe/H] moved, either (and either, not a white dwarf)."""
    out = dict(n=len(mh), afe=0, feh=0, any=0, any_not_wd=0)
    for s in range(0, len(mh), chunk):
        m, a = np.asarray(mh[s:s + chunk], np.float64), np.asarray(afe[s:s + chunk], np.float64)
        afe_c = np.clip(a, AFE_LO, AFE_HI)
        feh = mh_to_feh(m, afe_c)
        f_afe = np.abs(a - afe_c) > 1e-4
        f_feh = np.abs(np.clip(feh, FEH_LO + GRID_EPS, FEH_HI - GRID_EPS) - feh) > 1e-4
        both = f_afe | f_feh
        out["afe"] += int(f_afe.sum())
        out["feh"] += int(f_feh.sum())
        out["any"] += int(both.sum())
        if wd is not None:
            out["any_not_wd"] += int((both & ~wd[s:s + chunk]).sum())
    return out


class BCSet:
    """The seven simulator bands (FILTER_ORDER: LSST ugrizy, Roman F146) from the MIST tables, in chunks."""

    def __init__(self, quiet=True):
        with contextlib.redirect_stdout(io.StringIO()) if quiet else contextlib.nullcontext():
            self.rubin = MISTBolometricCorrection("Rubin")
            self.roman = MISTBolometricCorrection("F146")
            self.rubin._ensure_grid()
            self.roman._ensure_grid()
        for t in (self.rubin, self.roman):
            if (abs(t.grid[3].min() - FEH_LO) > 1e-9 or abs(t.grid[3].max() - FEH_HI) > 1e-9
                    or abs(t.grid[4].min() - AFE_LO) > 1e-9 or abs(t.grid[4].max() - AFE_HI) > 1e-9):
                raise RuntimeError("the MIST grid is not [Fe/H] -3..+0.5, [a/Fe] -0.2..+0.6: "
                                   "update AFE_LO/AFE_HI/FEH_LO/FEH_HI")

    def mags(self, logT, logg, mbol, mh, afe, clamp=True, chunk=400_000, bands=None):
        """(N, len(bands)) absolute AB magnitudes (default FILTER_ORDER); NaN where the grid has no BC."""
        bands = list(bands or FILTER_ORDER)
        n = len(logT)
        out = np.full((n, len(bands)), np.nan)
        for s in range(0, n, chunk):
            sl = slice(s, min(s + chunk, n))
            m, a = np.asarray(mh[sl], np.float64), np.asarray(afe[sl], np.float64)
            if clamp:
                m, a = clamp_to_grid(m, a)
            feh = mh_to_feh(m, a)
            for j, name in enumerate(bands):
                src = self.rubin if name.startswith("LSST") else self.roman
                out[sl, j] = np.asarray(mbol[sl], np.float64) - src.bc(name, logT[sl], logg[sl], feh, a)
        return out


def fix_low_mass_dwarfs(sub, name, lens_ml, return_mask=False):
    """dev82 mode: put the component's unevolved dwarfs on the lens-light M-L relation."""
    if name not in MS_FIX_COMP:
        return (sub, 0, np.zeros(len(sub), bool)) if return_mask else (sub, 0)
    t = lens_ml[lens_ml["comp"] == MS_FIX_COMP[name]]
    centres = ((t["m_lo"] + t["m_hi"]) / 2).to_numpy()
    sel = ((sub["CL"] == 5) & (sub["mass"] < MS_FIX_HI) & (sub["LSST_r"] < DARK_MAG)).to_numpy()
    m = sub["mass"].to_numpy()[sel]
    tpr = np.clip((MS_FIX_HI - m) / (MS_FIX_HI - MS_FIX_LO), 0.0, 1.0)
    bins = np.floor(m / MS_FIX_BIN).astype(int)
    bands = {"Roman_F146": "F146", "LSST_u": "u", "LSST_g": "g", "LSST_r": "r", "LSST_i": "i",
             "LSST_z": "z", "LSST_y": "y"}
    for col, short in bands.items():
        v = sub[col].to_numpy()[sel]
        med = pd.Series(v).groupby(bins).median()
        target = np.interp((med.index.to_numpy() + 0.5) * MS_FIX_BIN, centres, t[short].to_numpy())
        shift = pd.Series(target - med.to_numpy(), index=med.index)
        new = v + tpr * shift.reindex(bins).to_numpy()
        sub.loc[sub.index[sel], col] = new
    return (sub, int(sel.sum()), sel) if return_mask else (sub, int(sel.sum()))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def git_state():
    """`git rev-parse --short HEAD`, plus a note if the two builder scripts differ from it (read-only)."""
    try:
        head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True,
                              check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--", "BolometricCorrection.py",
                                "lens_ml_table.py"], capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown (git unavailable)"
    return head + (" + uncommitted changes in CMD/BolometricCorrection.py or lens_ml_table.py" if dirty else "")


def make_floor_fill(name, cat, idx, dwarf, thin_ref, emp, alpha, rng):
    """Synthetic dwarfs between FILL_MASS_LO and the component's mass floor (see the module docstring).

    Returns (arrays keyed like the working set, info). thin_ref = (sorted thin-disc dwarf masses, their Typ)."""
    d = idx[dwarf[idx]]
    mq = qmass(cat["Mass"][d])
    m_f = float(mq.min())
    n_win = int(((mq >= m_f - MASS_EPS) & (mq < m_f + FILL_NORM_WIDTH - MASS_EPS)).sum())
    n0 = n_win / FILL_NORM_WIDTH                       # stars per Msun at m_mid
    m_mid = m_f + 0.5 * FILL_NORM_WIDTH
    p = 1.0 - alpha
    x0, x1 = FILL_MASS_LO ** p, m_f ** p
    n_add = int(round(n0 * m_mid ** alpha * (x1 - x0) / p))
    u = rng.random(n_add)
    mass = (x0 + u * (x1 - x0)) ** (1.0 / p)           # inverse CDF of the truncated power law
    pool = d[(mq >= m_f - MASS_EPS) & (mq < m_f + FILL_DONOR_WIDTH - MASS_EPS)]
    donor = pool[rng.integers(0, len(pool), n_add)]
    thin_m, thin_typ = thin_ref
    lo = np.searchsorted(thin_m, mass - FILL_TYP_WIDTH - 1e-9, side="left")
    hi = np.searchsorted(thin_m, mass + FILL_TYP_WIDTH + 1e-9, side="right")
    if (hi <= lo).any():
        raise RuntimeError(f"{name}: no thin-disc dwarf within +-{FILL_TYP_WIDTH} Msun of a synthetic star")
    typ = thin_typ[lo + np.floor(rng.random(n_add) * (hi - lo)).astype(np.int64)]
    logT = emp.logT(mass)
    mbol = emp.mbol(mass)
    syn = {"mass": mass, "logT": logT, "Mbol": mbol, "logg": logg_formula(mass, logT, mbol),
           "M_H": cat["M_H"][donor].astype(np.float64), "a_Fe": cat["a_Fe"][donor].astype(np.float64),
           "Age": cat["Age"][donor].astype(np.float64), "Pop": cat["Pop"][donor],
           "CL": np.full(n_add, 5, np.int8), "Typ": typ.astype(np.float64)}
    info = {"alpha": alpha, "m_f": m_f, "n_window": n_win, "n0_per_Msun": n0, "m_mid": m_mid,
            "m_lo": FILL_MASS_LO, "n_add": n_add, "n_donors": int(len(pool)),
            "mass_min": float(mass.min()) if n_add else float("nan"),
            "mass_max": float(mass.max()) if n_add else float("nan")}
    return syn, info


def materialise(cat, idx, syn, rows, n_cat):
    """The working arrays (float64) of the component's rows `rows` (indices into real stars, then synthetic)."""
    real = idx[rows[rows < n_cat]]
    s = rows[rows >= n_cat] - n_cat
    W = {"mass": cat["Mass"][real], "logT": np.log10(cat["Teff"][real].astype(np.float64)),
         "Mbol": cat["Mbol"][real].astype(np.float64), "logg": cat["logg"][real].astype(np.float64),
         "M_H": cat["M_H"][real].astype(np.float64), "a_Fe": cat["a_Fe"][real].astype(np.float64),
         "Age": cat["Age"][real].astype(np.float64), "Pop": cat["Pop"][real], "CL": cat["CL"][real],
         "Typ": cat["Typ"][real].astype(np.float64)}
    nr = len(real)
    if syn is not None:
        for k in W:
            W[k] = np.concatenate([W[k], syn[k][s]])
    W["synthetic"] = np.concatenate([np.zeros(nr, bool), np.ones(len(W["mass"]) - nr, bool)])
    W["shifted"] = np.zeros(len(W["mass"]), bool)
    return W


def apply_empirical_shift(W, curve, emp):
    """Empirical-dwarf shift on the working arrays: shift the real dwarfs below MS_FIX_HI onto the empirical sequence."""
    sel = (W["CL"] == 5) & (W["Typ"] < DARK_TYP[0]) & (W["mass"] < MS_FIX_HI) & ~W["synthetic"]
    m = W["mass"][sel]
    t = taper(m)
    W["Mbol"][sel] += t * (emp.mbol(m) - curve_at(curve, "mbol", m))
    W["logT"][sel] += t * (emp.logT(m) - curve_at(curve, "logT", m))
    W["logg"][sel] = logg_formula(m, W["logT"][sel], W["Mbol"][sel])
    W["shifted"] |= sel
    return int(sel.sum())


def is_white_dwarf(typ):
    return (typ >= DARK_TYP[0] - 1e-4) & (typ <= DARK_TYP[1] + 1e-4)


def write_component(path, W, mags):
    cols = {"mass": W["mass"], "logT": W["logT"], "Mbol": W["Mbol"], "Age": W["Age"], "Pop": W["Pop"],
            "Roman_F146": mags[:, 6]}
    for j, f in enumerate(FILTER_ORDER[:6]):
        cols[f] = mags[:, j]
    cols["CL"], cols["Typ"] = W["CL"], W["Typ"]
    step = 1_000_000
    for s in range(0, len(W["mass"]), step):
        df = pd.DataFrame({k: v[s:s + step] for k, v in cols.items()})
        df.to_csv(path, sep=" ", index=False, float_format="%.4f", mode="w" if s == 0 else "a", header=(s == 0))


def build(a):
    t_start = time.time()
    out_dir = a.out_dir
    explicit = out_dir is not None
    if out_dir is None:
        out_dir = os.path.join(STAGING_ROOT, variant_name(a.input, a.dwarfs, a.fill_floor))
    if a.dwarfs == "dev82" and a.fill_floor != "none":
        raise SystemExit("--dwarfs dev82 is the legacy path: it cannot be combined with --fill-floor")
    if os.path.abspath(out_dir) == os.path.abspath("components"):
        if not explicit:
            raise SystemExit("refusing to write components/ without an explicit --out-dir components")
        print("NOTE: writing the PRODUCTION lists in components/")
    os.makedirs(out_dir, exist_ok=True)
    clamp = not a.no_grid_clamp
    print(f"input {shown(a.input)}; dwarfs {a.dwarfs}; fill-floor {a.fill_floor}; clamp {clamp}; out {shown(out_dir)}")

    cat, rep = read_catalogue(a.input)
    comp = component_codes(cat["Pop"])
    dwarf = (cat["CL"] == 5) & (cat["Typ"] < DARK_TYP[0])
    wd_all = is_white_dwarf(cat["Typ"])
    names = list(COMPONENTS)
    bad = {nm: int((cat["Age"][comp == k] > AGE_MAX[nm]).sum()) for k, nm in enumerate(names)}
    if any(bad.values()):
        raise RuntimeError(f"stars older than read_cmd's bound {AGE_MAX}: {bad}")

    emp = EmpiricalDwarfs() if (a.dwarfs == "empirical" or a.fill_floor != "none") else None
    curves, offsets = {}, {}
    if a.dwarfs == "empirical":
        for nm in EMPIRICAL_COMPONENTS:
            m = (comp == names.index(nm)) & dwarf
            curves[nm] = median_curve(cat["Mass"][m], cat["Mbol"][m], np.log10(cat["Teff"][m].astype(np.float64)))
            at = np.array(OFFSET_REPORT)
            offsets[nm] = {"n_ok": int(curves[nm]["ok"].sum()), "n_bins": len(curves[nm]["ok"]),
                           "dMbol": emp.mbol(at) - curve_at(curves[nm], "mbol", at),
                           "dlogT": emp.logT(at) - curve_at(curves[nm], "logT", at)}
    thin_ref = None
    if a.fill_floor != "none":
        m = (comp == 0) & dwarf
        order = np.argsort(qmass(cat["Mass"][m]), kind="stable")
        thin_ref = (qmass(cat["Mass"][m])[order], cat["Typ"][m][order].astype(np.float64))
    lens_ml = None
    if a.dwarfs == "dev82":
        lens_ml = pd.read_csv(a.lens_ml, sep=r"\s+", comment="#", header=None,
                              names=["comp", "m_lo", "m_hi", "n", "u", "g", "r", "i", "z", "y", "F146"])

    bcs = BCSet()
    rng_sub = np.random.default_rng(SUBSAMPLE_SEED)
    rng_fill = np.random.default_rng(FILL_SEED)
    info = {}
    # the synthetic stars are drawn in a fixed order (bulge, thick disc, halo) from one generator
    syn_all, fill_info = {}, {}
    if a.fill_floor != "none":
        for nm in FILL_COMPONENTS:
            k = names.index(nm)
            syn_all[nm], fill_info[nm] = make_floor_fill(nm, cat, np.flatnonzero(comp == k), dwarf, thin_ref, emp,
                                                         FILL_ALPHA[a.fill_floor], rng_fill)

    for k, nm in enumerate(names):
        t0 = time.time()
        idx = np.flatnonzero(comp == k)
        syn = syn_all.get(nm)
        n_cat = len(idx)
        n_syn = 0 if syn is None else len(syn["mass"])
        total = n_cat + n_syn
        mass_cat = float(cat["Mass"][idx].mean())
        mass_pop = (float(cat["Mass"][idx].sum()) + (float(syn["mass"].sum()) if n_syn else 0.0)) / total
        pop_wd = int(wd_all[idx].sum())
        # clamp counts over the whole population (real stars + synthetic)
        cc_pop = clamp_counts(cat["M_H"][idx], cat["a_Fe"][idx], wd_all[idx])
        if n_syn:
            cs = clamp_counts(syn["M_H"], syn["a_Fe"])
            for key in ("n", "afe", "feh", "any"):
                cc_pop[key] += cs[key]
            cc_pop["any_not_wd"] += cs["any"]
        n_shift_pop = 0
        if a.dwarfs == "empirical" and nm in EMPIRICAL_COMPONENTS:
            n_shift_pop = int((dwarf[idx] & (cat["Mass"][idx] < MS_FIX_HI)).sum())

        if a.dwarfs == "dev82" or total <= MAX_ROWS:
            rows_all = np.arange(total)
            rows = rows_all
        else:
            rows = np.sort(rng_sub.choice(total, MAX_ROWS, replace=False))
        W = materialise(cat, idx, syn, rows, n_cat)
        if a.dwarfs == "empirical" and nm in EMPIRICAL_COMPONENTS:
            apply_empirical_shift(W, curves[nm], emp)

        mags = bcs.mags(W["logT"], W["logg"], W["Mbol"], W["M_H"], W["a_Fe"], clamp=clamp)
        no_bc = ~np.isfinite(mags).all(axis=1)
        wd = is_white_dwarf(W["Typ"])
        dark = no_bc | wd
        mags[dark] = DARK_MAG
        print(f"{nm}: dark entries {int(dark.sum())} of {len(dark)} (no MIST BC {int(no_bc.sum())}, "
              f"white dwarfs {int(wd.sum())}, of which no BC and not a white dwarf {int((no_bc & ~wd).sum())})")
        examples = [(float(W["mass"][i]), float(10 ** W["logT"][i]), float(W["logg"][i]),
                     float(W["M_H"][i]), float(W["a_Fe"][i])) for i in np.flatnonzero(no_bc & ~wd)[:10]]

        n_fix = 0
        if a.dwarfs == "dev82":
            sub = pd.DataFrame({"mass": W["mass"], "CL": W["CL"]})
            for j, f in enumerate(FILTER_ORDER):
                sub[f] = mags[:, j]
            sub, n_fix, sel = fix_low_mass_dwarfs(sub, nm, lens_ml, return_mask=True)
            print(f"{nm}: {n_fix} low-mass dwarfs moved onto the lens_ml.dat relation")
            for j, f in enumerate(FILTER_ORDER):
                mags[:, j] = sub[f].to_numpy()
            W["shifted"] |= sel
            del sub
            if total > MAX_ROWS:
                keep = np.sort(rng_sub.choice(total, MAX_ROWS, replace=False))
                W = {key: v[keep] for key, v in W.items()}
                mags, dark = mags[keep], dark[keep]
                wd, no_bc = wd[keep], no_bc[keep]
        if not (np.isfinite(mags).all() and all(np.isfinite(W[key]).all() for key in
                ("mass", "logT", "Mbol", "Age", "logg", "Typ"))):
            raise RuntimeError(f"{nm}: non-finite value left after the dark-entry substitution")
        n_w = len(W["mass"])
        lum = ~dark
        cc_w = clamp_counts(W["M_H"], W["a_Fe"], wd)
        mr = mags[lum, 2]
        info[nm] = dict(
            n_catalogue=n_cat, n_synthetic=n_syn, n_population=total, n_written=n_w,
            n_dark=int(dark.sum()), n_dark_wd=int((dark & wd).sum()), n_dark_nobc=int((dark & ~wd).sum()),
            n_pop_wd=pop_wd, mean_mass_catalogue=mass_cat, mean_mass_population=mass_pop,
            mean_mass_written=float(W["mass"].mean()), median_Mr_lum=float(np.median(mr)),
            max_Mr_lum=float(mr.max()), n_Mr_gt_20=int((mr > 20.0).sum()),
            n_synthetic_written=int(W["synthetic"].sum()), n_shifted_pop=n_shift_pop, n_shifted_written=int(W["shifted"].sum()),
            clamp_pop=cc_pop, clamp_written=cc_w, nobc_examples=examples,
            mass_min_written=float(W["mass"].min()))

        write_component(os.path.join(out_dir, nm + ".dat"), W, mags)
        f32 = lambda v: np.asarray(v, np.float32)                       # noqa: E731
        np.savez_compressed(os.path.join(out_dir, nm + "_aux.npz"),
                            mass=f32(W["mass"]), Teff_used=f32(10.0 ** W["logT"]), Mbol_used=f32(W["Mbol"]),
                            logg_used=f32(W["logg"]), M_H=f32(W["M_H"]), a_Fe=f32(W["a_Fe"]), Pop=W["Pop"],
                            Age=f32(W["Age"]), CL=W["CL"], Typ=f32(W["Typ"]), synthetic=W["synthetic"],
                            shifted=W["shifted"], dark=dark)
        print(f"{nm:10s}: catalogue {n_cat:9d}  synthetic {n_syn:8d}  written {n_w:9d}  dark {int(dark.sum()):7d}  "
              f"<m>_pop = {mass_pop:.4f}  <m>_written = {info[nm]['mean_mass_written']:.4f}  "
              f"({time.time() - t0:.0f} s)")
        del W, mags, dark, wd, no_bc, rows, idx

    elapsed = time.time() - t_start
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    prov = provenance_text(a, out_dir, clamp, rep, info, offsets, fill_info, emp, elapsed, peak)
    with open(os.path.join(out_dir, "provenance.txt"), "w") as f:
        f.write(prov)
    print("\n" + prov)
    print("To adopt: run `python3 tools/sync_data_products.py` from the repo root -- it sets N1..N4 (written counts: thin, bulge, thick, halo) and MEANMASS_* (mean_mass_population) in config/data_products.h.")


def provenance_text(a, out_dir, clamp, rep, info, offsets, fill_info, emp, elapsed, peak_mb):
    L = ["# CMD star lists: provenance (BolometricCorrection.py)",
         f"# built {time.strftime('%Y-%m-%d %H:%M:%S')}; git {git_state()}; {elapsed / 60:.1f} min; peak RSS {peak_mb:.0f} MB",
         f"input: {shown(a.input)}  rows read {rep['n_rows']}  rows skipped (field count != {rep['n_fields']}): "
         f"{rep['n_skipped']}" + (f" at file lines {rep['bad_lines']}" if rep["n_skipped"] else ""),
         f"options: dwarfs={a.dwarfs} fill_floor={a.fill_floor} grid_clamp={'on' if clamp else 'off'} "
         f"max_rows={MAX_ROWS} subsample_seed={SUBSAMPLE_SEED} fill_seed={FILL_SEED}",
         f"output: {shown(out_dir)}",
         ("# the subsample is drawn after the BCs and the dev82 shift (their running medians need every star)"
          if a.dwarfs == "dev82" else "# the subsample is drawn before the BCs (a star's BC does not depend on the others)")
         + "; counts below are over the WRITTEN list unless they say 'population'",
         "# component n_catalogue n_synthetic n_population n_written n_dark n_dark_wd n_dark_nobc "
         "mean_mass_population mean_mass_catalogue mean_mass_written median_Mr_lum"]
    for nm, i in info.items():
        L.append(f"{nm} {i['n_catalogue']} {i['n_synthetic']} {i['n_population']} {i['n_written']} "
                 f"{i['n_dark']} {i['n_dark_wd']} {i['n_dark_nobc']} {i['mean_mass_population']:.4f} "
                 f"{i['mean_mass_catalogue']:.4f} {i['mean_mass_written']:.4f} {i['median_Mr_lum']:.3f}")
    L.append("# mean_mass_population -> config/data_products.h MEANMASS_* via tools/sync_data_products.py (real catalogue stars + synthetic, before the subsample);"
             " n_written -> N1..N4")
    L.append("# stars with no MIST BC that are not white dwarfs (mass, Teff, logg, [M/H], [a/Fe]):")
    for nm, i in info.items():
        L.append(f"{nm}: {i['n_dark_nobc']}" + "".join(
            f"  ({e[0]:.3f}, {e[1]:.0f}, {e[2]:.2f}, {e[3]:.3f}, {e[4]:.3f})" for e in i["nobc_examples"]))
    L.append("# BC-grid clamp ([a/Fe] clipped to -0.2..0.6, [M/H] shifted until [Fe/H] in -3..+0.5): "
             + ("ON" if clamp else "OFF (stars outside the grid are dark)"))
    L.append("# component population: n  afe_clipped  feh_moved  either  either_not_wd  |  written: n  either")
    for nm, i in info.items():
        p, w = i["clamp_pop"], i["clamp_written"]
        L.append(f"{nm} {p['n']} {p['afe']} {p['feh']} {p['any']} {p['any_not_wd']} | {w['n']} {w['any']}")
    L.append("# luminosity: faintest luminous M_r and how many luminous entries are fainter than M_r = 20 "
             "(read_cmd CHECKs M_r <= 20):")
    for nm, i in info.items():
        L.append(f"{nm}: max M_r {i['max_Mr_lum']:.3f}  n(M_r > 20) {i['n_Mr_gt_20']}  lowest mass written {i['mass_min_written']:.4f}")
    if a.dwarfs == "empirical":
        L.append(f"# empirical dwarf sequence: {emp.version}, {EMPIRICAL_TABLE}: {len(emp.mass)} rows "
                 f"({emp.spt[0]}..{emp.spt[-1]}), Msun {emp.mass[0]}-{emp.mass[-1]}; shifted dwarfs: CL = 5, Typ < 9, "
                 f"m < {MS_FIX_HI}, taper {MS_FIX_HI} -> {MS_FIX_LO}; thick disc and halo untouched")
        for nm in EMPIRICAL_COMPONENTS:
            i, o = info[nm], offsets[nm]
            L.append(f"{nm}: dwarfs shifted: {i['n_shifted_pop']} in the population, {i['n_shifted_written']} in the "
                     f"written list; median curve: {o['n_ok']} of {o['n_bins']} 0.01-Msun bins with >= {OFFSET_NMIN} stars")
            L.append(f"  offset X_E(m) - median_c X(m) before the taper, at m = "
                     + " / ".join(f"{m:.2f}" for m in OFFSET_REPORT))
            L.append("    dMbol " + " / ".join(f"{x:+.3f}" for x in o["dMbol"]))
            L.append("    dlogT " + " / ".join(f"{x:+.4f}" for x in o["dlogT"]))
    elif a.dwarfs == "dev82":
        L.append(f"# dev82 shift: lens_ml = {shown(a.lens_ml)}; dwarfs moved (written list): "
                 + ", ".join(f"{nm} {i['n_shifted_written']}" for nm, i in info.items()))
    if a.fill_floor != "none":
        L.append(f"# floor fill ({a.fill_floor}, alpha = {FILL_ALPHA[a.fill_floor]}): dN/dm = n0 (m/m_mid)^-alpha from "
                 f"{FILL_MASS_LO} Msun to the floor m_f; n0 = dwarfs in [m_f, m_f+{FILL_NORM_WIDTH}) / {FILL_NORM_WIDTH}; "
                 f"attributes from real dwarfs in [m_f, m_f+{FILL_DONOR_WIDTH}); Typ from the thin disc within "
                 f"+-{FILL_TYP_WIDTH} Msun")
        for nm, fi in fill_info.items():
            i = info[nm]
            L.append(f"{nm}: m_f = {fi['m_f']:.4f}  dwarfs in window {fi['n_window']}  n0 = {fi['n0_per_Msun']:.1f} per Msun "
                     f"at m_mid = {fi['m_mid']:.4f}  N_add = {fi['n_add']}  N_add/N_c = {fi['n_add'] / i['n_catalogue']:.4f}  "
                     f"mass range of the draws {fi['mass_min']:.4f}-{fi['mass_max']:.4f}  donors {fi['n_donors']}  "
                     f"synthetic in the written list {i['n_synthetic_written']}")
    return "\n".join(L) + "\n"


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Build the per-component star lists from a Besancon catalogue (see the module docstring).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--input", default=DEFAULT_INPUT, help="Besancon catalogue (damaged rows are skipped and counted)")
    ap.add_argument("--out-dir", default=None,
                    help=f"output directory; default {STAGING_ROOT}/<input>_<dwarfs>_fill<floor>. "
                         "Writing components/ needs it explicitly.")
    ap.add_argument("--dwarfs", choices=["besancon", "empirical", "dev82"], default="empirical",
                    help="besancon: no change; empirical: thin disc and bulge dwarfs below 0.7 Msun onto the "
                         "Pecaut & Mamajek sequence; dev82: the legacy shift onto lens_ml.dat")
    ap.add_argument("--fill-floor", choices=["none", "kroupa", "koshimoto"], default="none",
                    help="add synthetic dwarfs between 0.08 Msun and the bulge/thick-disc/halo mass floor, "
                         "with alpha = 1.3 (kroupa) or 1.16 (koshimoto)")
    ap.add_argument("--no-grid-clamp", action="store_true",
                    help="no BC-grid clamp: stars off the MIST grid go dark (reproduces the legacy lists)")
    ap.add_argument("--lens-ml", default="components/lens_ml.dat", help="lens_ml.dat for --dwarfs dev82 (read only)")
    return ap.parse_args(argv)


def shown(path):
    """A path as it is printed: relative to CMD/ when it lies inside it, absolute otherwise."""
    rel = os.path.relpath(path)
    return path if rel.startswith("..") else rel


def main(argv=None):
    a = parse_args(argv)
    # a relative user path is taken from where the command was typed, falling back to CMD/ (the defaults
    # are CMD-relative); everything is made absolute before the script changes into CMD/
    here = os.path.dirname(os.path.abspath(__file__))
    a.input = os.path.abspath(a.input) if os.path.exists(a.input) else os.path.join(here, a.input)
    a.lens_ml = os.path.abspath(a.lens_ml) if os.path.exists(a.lens_ml) else os.path.join(here, a.lens_ml)
    if a.out_dir is not None:
        a.out_dir = os.path.abspath(a.out_dir)
    os.chdir(here)
    build(a)


if __name__ == "__main__":
    main()
