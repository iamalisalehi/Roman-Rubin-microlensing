"""
qc_bolometric_correction.py

Pre-flight validation for the four CMD component files produced by
BolometricCorrection.py (thin_disk.dat, bulge.dat, thick_disk.dat, halo.dat),
meant to run BEFORE those files are ever handed to the C++ simulation's
read_cmd().

This mirrors every CHECK() read_cmd performs, using the exact same bounds,
but reports ALL violations found per file (with counts and example row
indices) rather than throwing on the first one -- read_cmd's own CHECK()s
should never actually fire if this script reports a clean pass. It also
runs a few extra sanity checks the C++ side doesn't perform at all
(column order/name matching, duplicate rows, a broad physical-plausibility
band on magnitudes, and an int-vs-double parsing hazard check on Age).

Usage:
    python qc_bolometric_correction.py

Edit the CONFIGURATION block below to match your file layout.
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Configuration -- edit these to match your setup
# ---------------------------------------------------------------------------

COMPONENT_FILES = {
    "thin_disk":  "../CMD/components/thin_disk.dat",
    "bulge":      "../CMD/components/bulge.dat",
    "thick_disk": "../CMD/components/thick_disk.dat",
    "halo":       "../CMD/components/halo.dat",
}

# Set to None to skip auto-detection and rely purely on MANUAL_N below.
BULGE_H_PATH = "../Bulge.h"

# Fallback expected row counts (N1..N4 in Bulge.h), used only for any
# component where BULGE_H_PATH couldn't be read or didn't contain the
# matching constant. Fill these in as a backup if you don't have Bulge.h
# handy, or leave as None to skip the row-count check for that component.
MANUAL_N = {
    "thin_disk":  None,
    "bulge":      None,
    "thick_disk": None,
    "halo":       None,
}

# The exact column order read_cmd expects. read_cmd reads POSITIONALLY
# (operator>> in sequence), not by column name, so order matters as much
# as presence -- this is deliberately checked as an exact list, not a set.
EXPECTED_COLUMNS = [
    "mass", "logT", "Mbol", "Age", "Pop",
    "Roman_F146", "LSST_u", "LSST_g", "LSST_r", "LSST_i", "LSST_z", "LSST_y",
    "CL", "Typ",
]

# Mirrors read_cmd()'s per-component CHECK() bounds exactly. Keep this in
# sync with helper.cpp if those bounds ever change there.
BOUNDS = {
    "thin_disk":  {"mab_r_max": 20.0, "age_max": 10, "pop_valid": set(range(1, 8))},
    "bulge":      {"mab_r_max": 18.0, "age_max": 10, "pop_valid": {10}},
    "thick_disk": {"mab_r_max": 20.0, "age_max": 8,  "pop_valid": {8, 11}},
    "halo":       {"mab_r_max": 20.0, "age_max": 9,  "pop_valid": {9}},
}

CL_MAX = 7
TYP_MAX = 9.0

# Broad plausibility band for any magnitude column -- not asserted anywhere
# in the C++ code, but a value outside this range almost certainly means a
# computation bug (e.g. a NaN that dropna() somehow missed, or an
# interpolator returning garbage) rather than a genuine faint/bright star.
MAG_PLAUSIBLE_MIN, MAG_PLAUSIBLE_MAX = -15.0, 40.0

MAG_COLUMNS = ["Roman_F146", "LSST_u", "LSST_g", "LSST_r", "LSST_i", "LSST_z", "LSST_y"]


# ---------------------------------------------------------------------------
# Bulge.h constant extraction
# ---------------------------------------------------------------------------

def extract_N_from_bulge_h(path):
    """Pulls N1..N4 directly out of Bulge.h so this script can't silently
    drift out of sync with the compiled C++ constants. Falls back to
    MANUAL_N (component-by-component) if the file or a given constant
    can't be found."""
    const_to_component = {"N1": "thin_disk", "N2": "bulge", "N3": "thick_disk", "N4": "halo"}
    result = dict(MANUAL_N)

    if path is None:
        return result

    try:
        text = Path(path).read_text()
    except OSError as e:
        print(f"[warn] could not read {path} ({e}); using MANUAL_N as fallback.")
        return result

    for const_name, component in const_to_component.items():
        m = re.search(rf"\b{const_name}\s*=\s*(\d+)", text)
        if m:
            result[component] = int(m.group(1))
        else:
            print(f"[warn] could not find '{const_name}' in {path}; "
                  f"using MANUAL_N['{component}'] instead.")

    return result


# ---------------------------------------------------------------------------
# Per-file checks
# ---------------------------------------------------------------------------

def check_component(name, path, expected_n):
    print(f"\n{'=' * 72}\n{name}  :  {path}\n{'=' * 72}")
    issues = []
    warnings = []

    try:
        df = pd.read_csv(path, sep=r"\s+")
    except Exception as e:
        print(f"  FATAL: could not read file at all: {e}")
        return False

    # --- column presence AND order (read_cmd reads positionally) ---
    actual_columns = list(df.columns)
    if actual_columns != EXPECTED_COLUMNS:
        issues.append(
            "Column order/name mismatch -- read_cmd reads columns positionally, "
            "so this WILL misalign every field if run as-is.\n"
            f"      expected: {EXPECTED_COLUMNS}\n"
            f"      actual:   {actual_columns}"
        )
        # If columns don't match at all, the rest of the checks below would
        # just raise KeyError noise -- report what we found and stop here.
        print(f"\n  {len(issues)} ISSUE(S) FOUND:")
        for i, msg in enumerate(issues, 1):
            print(f"   {i}. {msg}")
        return False

    # --- row count vs. Bulge.h's compiled-in constant ---
    n_rows = len(df)
    if expected_n is not None:
        if n_rows != expected_n:
            issues.append(
                f"Row count mismatch: file has {n_rows} rows, Bulge.h's constant "
                f"expects {expected_n}. read_cmd's CHECK(j == N...) WILL fire on this "
                f"file as-is -- either update the constant in Bulge.h to {n_rows}, "
                f"or find out why the row count changed (e.g. dropna() removing more "
                f"rows than expected)."
            )
    else:
        warnings.append(
            "No expected row count available for this component (Bulge.h constant "
            "not found/not provided) -- row-count mismatch against read_cmd's "
            "CHECK(j == N...) cannot be checked here."
        )

    # --- NaN / inf anywhere (dropna() should have already removed these) ---
    n_nan = int(df.isna().sum().sum())
    numeric_df = df.select_dtypes(include=[np.number])
    n_inf = int(np.isinf(numeric_df.to_numpy()).sum())
    if n_nan:
        issues.append(f"{n_nan} NaN values found across the file (dropna() should have removed all of these)")
    if n_inf:
        issues.append(f"{n_inf} inf values found across the file")

    # --- exact duplicate rows ---
    n_dup = int(df.duplicated().sum())
    if n_dup:
        warnings.append(f"{n_dup} exact duplicate rows (not necessarily wrong, but worth a glance)")

    b = BOUNDS[name]

    def report(cond_series, label):
        bad = ~cond_series
        n_bad = int(bad.sum())
        if n_bad:
            bad_idx = df.index[bad][:5].tolist()
            issues.append(f"{label}: {n_bad} row(s) fail (first offending row indices: {bad_idx})")

    # --- exact mirrors of read_cmd's per-row CHECK()s ---
    report(df["mass"] >= 0.0, "mass >= 0.0")
    report(df["logT"] >= 0.0, "logT >= 0.0")
    report(df["LSST_r"] <= b["mab_r_max"], f"LSST_r <= {b['mab_r_max']}")
    report(df["Age"] <= b["age_max"], f"Age <= {b['age_max']}")
    report(df["CL"] <= CL_MAX, f"CL <= {CL_MAX}")
    report(df["Typ"] <= TYP_MAX, f"Typ <= {TYP_MAX}")

    # --- population purity: did the Pop-based split in save_components() work? ---
    bad_pop = ~df["Pop"].isin(b["pop_valid"])
    if bad_pop.any():
        found = sorted(df.loc[bad_pop, "Pop"].unique().tolist())
        issues.append(
            f"Pop contamination: {int(bad_pop.sum())} row(s) have Pop values outside "
            f"{sorted(b['pop_valid'])} for this component (found: {found}) -- "
            f"these rows shouldn't be in this file at all."
        )

    # --- broad magnitude plausibility band (not checked in C++ at all) ---
    for col in MAG_COLUMNS:
        out_of_band = ~df[col].between(MAG_PLAUSIBLE_MIN, MAG_PLAUSIBLE_MAX)
        if out_of_band.any():
            warnings.append(
                f"{col}: {int(out_of_band.sum())} row(s) fall outside the broad "
                f"[{MAG_PLAUSIBLE_MIN}, {MAG_PLAUSIBLE_MAX}] plausibility band -- "
                f"almost certainly a computation artifact, not a real faint/bright star."
            )

    # --- Age fractional-value hazard: int-vs-double mismatch risk in the C++ struct ---
    frac_ages = int((df["Age"] % 1 != 0).sum())
    if frac_ages:
        warnings.append(
            f"{frac_ages} row(s) have non-integer Age values (e.g. 0.15 Gyr, from the "
            f"youngest thin-disk sub-population). Confirm the corresponding age_* "
            f"field in the CMD struct (Bulge.h) is declared as a floating-point type, "
            f"not int -- reading '0.15' into an int with operator>> only consumes "
            f"the leading '0' and desynchronizes every field read after it for that row."
        )

    # --- did every filter actually get computed independently, or all identical? ---
    if df[MAG_COLUMNS].nunique(axis=1).eq(1).any():
        n_flat = int(df[MAG_COLUMNS].nunique(axis=1).eq(1).sum())
        warnings.append(
            f"{n_flat} row(s) have IDENTICAL values across all {len(MAG_COLUMNS)} "
            f"filter columns -- worth a look, since distinct filters should almost "
            f"never produce exactly equal magnitudes for a real star."
        )

    # --- distribution summary, for eyeballing ---
    print(df[["mass", "logT", "Mbol", "Age", "Roman_F146", "LSST_r"]].describe().to_string())

    if warnings:
        print(f"\n  {len(warnings)} WARNING(S) (not blocking, worth a look):")
        for i, msg in enumerate(warnings, 1):
            print(f"   - {msg}")

    if issues:
        print(f"\n  {len(issues)} ISSUE(S) FOUND (these WILL trip read_cmd's CHECK()s):")
        for i, msg in enumerate(issues, 1):
            print(f"   {i}. {msg}")
        return False

    print("\n  All checks passed.")
    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    n_map = extract_N_from_bulge_h(BULGE_H_PATH)

    results = {}
    for name, path in COMPONENT_FILES.items():
        if not Path(path).exists():
            print(f"\n{'=' * 72}\n{name}  :  {path}\n{'=' * 72}\n  FATAL: file does not exist.")
            results[name] = False
            continue
        results[name] = check_component(name, path, n_map.get(name))

    print(f"\n{'=' * 72}\nSUMMARY\n{'=' * 72}")
    all_ok = True
    for name, ok in results.items():
        print(f"  {name:12s} : {'PASS' if ok else 'FAIL'}")
        all_ok = all_ok and ok

    if all_ok:
        print("\nAll four component files are ready for read_cmd(). "
              "Safe to run the C++ simulation -- its CHECK()s should now be a pure failsafe.")
        sys.exit(0)
    else:
        print("\nFix the issue(s) above before running the C++ simulation -- "
              "read_cmd()'s CHECK()s WILL fire otherwise.")
        sys.exit(1)


if __name__ == "__main__":
    main()
