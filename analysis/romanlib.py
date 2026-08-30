"""Shared reader for the Roman+Rubin forecast outputs.

Every analysis script goes through this module. That is deliberate: the C++ side reports
"not measured" as an explicit -1.0 sentinel rather than NaN, and three separate bugs during
development came from a -1 being summed as though it were a measurement. Encoding the
sentinel rules once, here, is what stops that class of error reappearing in Python.

The two rules that matter
-------------------------
1. **A sentinel is never a measurement.** -1.0 in any sigma, condition number or relative
   error means "this could not be determined". It is not a small error, not a large one, and
   must never enter a sum, a mean, a ratio or a histogram.

2. **Gate on okA/okB, never on flagi.** `flagi` is set inside FisherM and is NOT reset per
   event, so it reads as the previous characterised event's value on rows where nothing was
   characterised (OPEN_ITEMS.md). `okA_J` and `okB_J` are reset every event and are the
   flags to trust.

Two indexing systems that look alike
------------------------------------
- `magb_*` / `blend_*` are per FILTER: u g r i z y are Rubin's, F146 is Roman's.
- `mbs0` / `fb0` are Rubin (r-band); `mbs1` / `fb1` are Roman (F146). These are per
  TELESCOPE and are what the Fisher matrix actually fits.
"""

from __future__ import annotations

import os
import re

import numpy as np
import pandas as pd

SENTINEL = -1.0

# Detection taxonomy -- see DetClass in Bulge.h.
DET_CLASS = {
    0: "none",
    1: "joint-only",       # neither telescope alone would have found it
    2: "Rubin+joint",
    3: "Roman+joint",
    4: "both+joint",
    5: "ANOMALY",          # a telescope detected it but the joint test did not
}

# Characterisability taxonomy -- see SynergyClass in Bulge.h.
SYN_CLASS = {
    0: "none",
    1: "both-alone",
    2: "Rubin-only",       # only Rubin characterises alone; Roman still sharpens the joint fit
    3: "Roman-only",
    4: "joint-only",       # NEITHER alone, but the joint fit works -- pure rescue
}

# Where t0 fell relative to Roman's observing seasons -- see T0Zone in Bulge.h.
T0_ZONE = {
    0: "in-season",
    1: "in-gap",           # bracketed by Roman data: the gap-filling regime
    2: "off-mission",      # before Roman's first epoch or after its last: Rubin-only by construction
}

SURVEYS = {"joint": "J", "rubin": "L", "roman": "R"}

# Column layout of MapLMC2.dat, one row per AGGREGATED sightline, in the order the
# `fil3 <<` block in Bulge_LSST.cpp writes them. Each of the first 22 quantities is
# written as a pair: [0] over all recorded events, [1] over detected events only.
_MAP_PAIRS = ["tE", "RE", "piE", "tetE", "Vt", "u0", "Ml", "opd", "Dl", "Ds", "vl",
              "vs", "mbs", "fb", "fwhm", "vsn", "DelT", "Struc", "murel", "Map",
              "nbl", "Ext"]
MAP_COLS = ([f"{n}_{i}" for n in _MAP_PAIRS for i in (0, 1)]
            + ["EffiD", "EffiL", "log10_EFF", "log10_Gamma", "log10_Neven",
               "Eru0", "ErtE", "Erfb", "ErpiE", "ErtetE", "Erml", "Erdl", "Ermul", "Ermus",
               "nsim", "numd0", "numd1", "nerr", "nri", "nde",
               "log10_Rostart", "log10_Nstart", "log10_nstart"])


def load_events(path):
    """Read the per-event table (test2.dat) written by the `filg_in <<` block.

    Column names come from the file's own `#` header, not from a list hardcoded here, so a
    schema change surfaces as a loud mismatch rather than a silent misalignment.
    """
    with open(path) as fh:
        header = fh.readline()
    if not header.startswith("#"):
        raise ValueError(
            f"{path}: no '#' header line. Files written before Step D1 have no header and a "
            f"different column set; they cannot be read with this loader."
        )
    cols = header.lstrip("#").split()

    df = pd.read_csv(path, sep=r"\s+", comment="#", header=None, names=cols)
    if df.shape[1] != len(cols):
        raise ValueError(f"{path}: header names {len(cols)} columns, data has {df.shape[1]}")

    df["detClsName"] = df["detCls"].map(DET_CLASS)
    df["synClassName"] = df["synClass"].map(SYN_CLASS)
    df["t0zoneName"] = df["t0zone"].map(T0_ZONE)
    return df


def load_sightlines(path):
    """Read MapLMC2.dat, one row per aggregated sightline.

    Note this file is opened in APPEND mode by the simulator, so a re-run without clearing
    it first silently concatenates two runs (OPEN_ITEMS.md). `nri`/`nde` restarting from
    zero part-way down the file is the signature.
    """
    df = pd.read_csv(path, sep=r"\s+", header=None, names=MAP_COLS)
    return df


def load_provenance(path):
    """Parse run_provenance.txt into a dict of strings.

    The sightline-outcome block is appended at the END of a run, so its absence means the
    run did not finish -- which is exactly when you must not quote a density.
    """
    prov = {}
    with open(path) as fh:
        for line in fh:
            m = re.match(r"^#\s+(\w+)\s+(.*?)\s*(?:#.*)?$", line)
            if m:
                prov[m.group(1)] = m.group(2).strip()
    return prov


# ---------------------------------------------------------------------------------------
# Sentinel-aware accessors. Use these instead of touching the columns directly.
# ---------------------------------------------------------------------------------------

def sigma(df, param, survey):
    """1-sigma forecast for `param` from `survey`, NaN where it was not measured.

    param  : "tE" | "piE"    (photometric, gated on okA)
             "tetE"          (astrometric, gated on okB)
             "Ml"            (derived from tetE and piE, gated on its own positivity)
    survey : "joint" | "rubin" | "roman"

    Gating is on the ok flag AND on positivity. Both are needed: the flag can be set while
    an individual parameter is still a sentinel, because each survey partition fits its own
    active parameter subset (activePhotParams in Bulge.h). An event Roman detects with no
    Rubin epochs has a valid joint fit in which the Rubin blend fraction was never a free
    parameter -- exactly the case that aborted the 2026-08-29 run.
    """
    q = SURVEYS[survey]
    if param in ("tE", "piE"):
        col, gate = f"sig{param}_{q}", f"okA_{q}"
    elif param == "tetE":
        col, gate = f"sigtetE_{q}", f"okB_{q}"
    elif param == "Ml":
        col, gate = f"relMl_{q}", None
    else:
        raise KeyError(f"unknown param {param!r}")

    v = df[col].astype(float)
    ok = (v > 0.0)
    if gate is not None:
        ok &= (df[gate] == 1)
    return v.where(ok)


def characterized(df, survey):
    """Abrams et al. 2025 characterization criterion: tE > 2*sigma_tE AND piE > 2*sigma_piE.

    Deliberately their criterion, not a stricter single-parameter one, so our Rubin-alone
    numbers are directly comparable to their published values. They note it is appropriately
    looser than sigma_tE/tE < 0.1 because two parameters are constrained at once.

    Returns a boolean Series; events where either sigma is unmeasured are False, never NaN.
    """
    s_tE = sigma(df, "tE", survey)
    s_piE = sigma(df, "piE", survey)
    return ((df["tE"] > 2.0 * s_tE) & (df["piE"] > 2.0 * s_piE)).fillna(False)


def ratio_joint_over(df, param, survey):
    """Per-event sigma_joint / sigma_<survey>, NaN unless BOTH were measured.

    Always per event, then aggregate -- never a ratio of two separately-averaged sigmas.
    Adding data cannot worsen a Fisher forecast, so this is bounded above by 1 in exact
    arithmetic. In practice a handful of events exceed 1 by ~1e-3; every one of them has a
    photometric condition number above 1e9, where double precision has already lost most of
    its digits (see OPEN_ITEMS.md). A violation on a WELL-conditioned event would be a
    partitioning bug; on an ill-conditioned one it is round-off on a forecast that was
    meaningless anyway.
    """
    return sigma(df, param, "joint") / sigma(df, param, survey)


def detected(df, survey):
    """Boolean: did this survey's own detection test fire?"""
    return df[{"joint": "detJ", "rubin": "detL", "roman": "detR"}[survey]] == 1


def check_monotonicity(df, params=("tE", "piE", "tetE", "Ml"), tol=1e-9):
    """Assert sigma_joint <= sigma_single wherever both exist. Returns a list of violations.

    This is a physics invariant, not a preference: the joint Fisher matrix is the sum of the
    per-survey ones, so adding data cannot increase a forecast error. A violation means a
    bug in the partitioning, not a marginal case.
    """
    bad = []
    for p in params:
        for s in ("rubin", "roman"):
            r = ratio_joint_over(df, p, s).dropna()
            if len(r) and r.max() > 1.0 + tol:
                bad.append((p, s, float(r.max()), int((r > 1.0 + tol).sum())))
    return bad


# Where the simulator writes run_provenance.txt, and where a copy is sometimes kept
# beside an archived run. Searched in order. Hardcoding only the repo-root name meant every
# figure silently carried NO provenance stamp, because the file the simulator actually
# writes lives under files/MONTLMC/files/.
PROVENANCE_SEARCH = ("run_provenance.txt",
                     "files/MONTLMC/files/run_provenance.txt")


def find_provenance(explicit=None, near=None):
    """Locate run_provenance.txt: an explicit path, then beside the events file, then the
    standard locations. Returns None if there is none -- callers must say so on the figure
    rather than print an unlabelled plot."""
    cands = []
    if explicit:
        cands.append(explicit)
    if near:
        cands.append(os.path.join(os.path.dirname(os.path.abspath(near)),
                                  "run_provenance.txt"))
    cands.extend(PROVENANCE_SEARCH)
    for c in cands:
        if c and os.path.exists(c):
            return c
    return None


def describe(path_events, path_prov=None):
    """One-line provenance summary to print at the top of every figure-producing script."""
    parts = [f"events={os.path.basename(path_events)}"]
    path_prov = find_provenance(path_prov, near=path_events)
    if path_prov and os.path.exists(path_prov):
        prov = load_provenance(path_prov)
        for k in ("git_commit", "stride", "events_target", "sightlines_aggregated"):
            if k in prov:
                parts.append(f"{k}={prov[k]}")
        if "sightlines_aggregated" not in prov:
            parts.append("INCOMPLETE-RUN")
    else:
        parts.append("provenance=NOT FOUND")
    return "  ".join(parts)
