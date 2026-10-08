"""Shared reader for the Roman+Rubin forecast outputs.

Every analysis script goes through this module. The C++ side reports "not measured" as an explicit
-1.0 sentinel rather than NaN, and the sentinel rules are encoded once here.

Two rules:
1. A sentinel is never a measurement. -1.0 in any sigma, condition number or relative error means
   "could not be determined" and must never enter a sum, a mean, a ratio or a histogram.
2. Gate on okA/okB, never on flagi. `flagi` is set inside FisherM and is not reset per event, so it
   carries the previous characterised event's value on rows where nothing was characterised.
   `okA_J` and `okB_J` are reset every event.

Two indexing systems that look alike:
- `magb_*` / `blend_*` are per FILTER: u g r i z y are Rubin's, F146 is Roman's.
- `mbs0` / `fb0` are Rubin (r-band); `mbs1` / `fb1` are Roman (F146). These are per TELESCOPE and
  are what the Fisher matrix fits.
"""

from __future__ import annotations

import io
import os
import re
import warnings

import numpy as np
import pandas as pd

from cparams import P            # the C++ headers' own numbers (analysis/cparams.py)

SENTINEL = -1.0

# Detection taxonomy -- see DetClass in include/fisher/fisher.h.
DET_CLASS = {
    0: "none",
    1: "joint-only",       # neither telescope alone would have found it
    2: "Rubin+joint",
    3: "Roman+joint",
    4: "both+joint",
    5: "ANOMALY",          # a telescope detected it but the joint test did not
}

# Characterisability taxonomy -- see SynergyClass in include/fisher/fisher.h.
SYN_CLASS = {
    0: "none",
    1: "both-alone",
    2: "Rubin-only",       # only Rubin characterises alone; Roman still sharpens the joint fit
    3: "Roman-only",
    4: "joint-only",       # NEITHER alone, but the joint fit works -- pure rescue
}

# Where t0 fell relative to Roman's observing seasons -- see T0Zone in include/surveys/schedule.h.
T0_ZONE = {
    0: "in-season",
    1: "in-gap",           # bracketed by Roman data: the gap-filling regime
    2: "off-mission",      # before Roman's first epoch or after its last: Rubin-only by construction
}

SURVEYS = {"joint": "J", "rubin": "L", "roman": "R"}

# Column layout of MapLMC2.dat, one row per aggregated sightline, in the order the `fil3 <<` block in
# src/sim/sightline.cpp writes them. Each of the first 22 quantities is a pair: [0] over all recorded
# events, [1] over detected events only.
_MAP_PAIRS = ["tE", "RE", "piE", "tetE", "Vt", "u0", "Ml", "opd", "Dl", "Ds", "vl",
              "vs", "mbs", "fb", "fwhm", "vsn", "DelT", "Struc", "murel", "Map",
              "nbl", "Ext"]
MAP_COLS = ([f"{n}_{i}" for n in _MAP_PAIRS for i in (0, 1)]
            + ["EffiD", "EffiL", "log10_EFF", "log10_Gamma", "log10_Neven",
               "Eru0", "ErtE", "Erfb", "ErpiE", "ErtetE", "Erml", "Erdl", "Ermul", "Ermus",
               "nsim", "numd0", "numd1", "nerr", "nri", "nde",
               "log10_Rostart", "log10_Nstart", "log10_nstart",
               # w_area is the deg^2 of sky this sightline stands for (constant unless
               # --stride-roman is used); lon/lat tie a map row to the events it produced.
               # Older files lack these three; load_sightlines() detects that by width.
               "w_area", "lon", "lat",
               # Median 5-sigma depth per LSST band (st.rubinDepthMed; -inf where the band has no
               # visit), which Rubin's acceptance in preselectEvent compares each draw's peak with.
               # Older files lack these six (NaN here).
               "depth5_u", "depth5_g", "depth5_r", "depth5_i", "depth5_z", "depth5_y"])
_N_DEPTH, _N_E1 = 6, 3


def _narrow(df, keep64):
    """float64 -> float32 and int64 -> int32, except the columns named in keep64."""
    out = {}
    for c in df.columns:
        if c in keep64:
            out[c] = df[c]
        elif df[c].dtype == np.float64:
            out[c] = df[c].astype(np.float32)
        elif df[c].dtype == np.int64:
            out[c] = df[c].astype(np.int32)
        else:
            out[c] = df[c]
    return pd.DataFrame(out, index=df.index)


def sightline_index(df):
    """(codes, keys): an int32 sightline code per row, and keys[code] = (lon, lat) to 3 dp.

    The keys are the rounded values of each sightline's first row, so dict lookups keyed on them
    (nsim, density_profile) match the (lon, lat) tuples used elsewhere, without a per-row tuple list.
    """
    lon = df["lon"].to_numpy(np.float64).round(3)
    lat = df["lat"].to_numpy(np.float64).round(3)
    # rint(x*1000) is the integer numpy's round(3) went through, so equal packs <=> equal pairs.
    packed = (np.rint(lon * 1000).astype(np.int64) + 1_000_000) * 10_000_000 \
        + (np.rint(lat * 1000).astype(np.int64) + 1_000_000)
    _, first, codes = np.unique(packed, return_index=True, return_inverse=True)
    keys = [(float(lon[i]), float(lat[i])) for i in first]
    return codes.astype(np.int32).ravel(), keys


def sightline_groups(codes):
    """Yield (code, row positions) per sightline, from sightline_index()'s codes."""
    order = np.argsort(codes, kind="stable")
    bounds = np.flatnonzero(np.diff(codes[order])) + 1
    for pos in np.split(order, bounds):
        if len(pos):
            yield int(codes[pos[0]]), pos


def load_events(path, keep=None, chunksize=None, usecols=None, narrow=False,
                keep64=("lon", "lat")):
    """Read the per-event table (test5.dat) written by the `filg_in <<` block.

    Column names come from the file's own `#` header, so a schema change surfaces as a loud mismatch.

    keep, chunksize -- read in chunks and keep only the rows `keep(chunk)` selects. The production
        table is millions of rows x ~90 float64 columns (~4 GB resident), enough to exhaust memory
        on a small machine; filtering per chunk holds the peak at one chunk plus the survivors.
        `keep` takes a chunk and returns a boolean mask over it, e.g.

            df = R.load_events(path, keep=lambda c: c["detJ"] == 1, chunksize=500_000)

        The default (keep=None) reads the whole file in one pass.

    usecols -- keep only these columns. For statistics over all draws, which cannot discard rows.
        The `detCls`/`synClass`/`t0zone` label columns are added only if their source column survives.

    narrow -- store floats as float32 and integers as int32, halving memory, for statistics that
        need every draw. float32 keeps 7 significant digits; callers must do arithmetic in float64.
        Columns in `keep64` stay float64 because they are matched, not computed with: lon/lat key
        the sightline lookup, and float64(float32(-0.319)) differs from the map file's -0.319.
    """
    with open(path) as fh:
        header = fh.readline()
    if not header.startswith("#"):
        raise ValueError(
            f"{path}: no '#' header line. Files from older versions of the simulator have no header and a "
            f"different column set; they cannot be read with this loader."
        )
    cols = header.lstrip("#").split()

    reader_kw = dict(sep=r"\s+", comment="#", header=None, names=cols)
    if usecols is not None:
        missing = [c for c in usecols if c not in cols]
        if missing:
            raise ValueError(f"{path}: no such column(s) {missing}")
        reader_kw["usecols"] = list(usecols)
    if keep is None and chunksize is None and not narrow:
        df = pd.read_csv(path, **reader_kw)
    else:
        parts = []
        with pd.read_csv(path, chunksize=chunksize or 500_000, **reader_kw) as it:
            for chunk in it:
                if keep is not None:
                    chunk = chunk[keep(chunk)]
                if narrow:
                    chunk = _narrow(chunk, keep64)
                parts.append(chunk)
        df = (pd.concat(parts, ignore_index=True) if parts
              else pd.DataFrame(columns=cols))

    expected = len(reader_kw.get("usecols", cols))
    if df.shape[1] != expected:
        raise ValueError(f"{path}: expected {expected} columns, data has {df.shape[1]}")

    for src, dst, mapping in (("detCls", "detClsName", DET_CLASS),
                              ("synClass", "synClassName", SYN_CLASS),
                              ("t0zone", "t0zoneName", T0_ZONE)):
        if src in df.columns:
            df[dst] = df[src].map(mapping)
    return df


def load_sightlines(path):
    """Read MapLMC2.dat, one row per aggregated sightline.

    The simulator opens this file in append mode, so a re-run without clearing it concatenates two
    runs (`nri`/`nde` restarting from zero part-way down is the signature).

    Malformed lines are dropped with a warning rather than killing the read. The stream is never
    flushed, so a killed run loses its buffered tail and the next run's first row lands on the
    fragment, leaving one line of the wrong width.
    """
    with open(path) as fh:
        lines = fh.readlines()
    widths = {}
    for i, line in enumerate(lines, start=1):
        if line.strip():
            widths.setdefault(len(line.split()), []).append(i)
    if len(widths) > 1:
        expected = max(widths, key=lambda w: len(widths[w]))
        bad = sorted(i for w, ii in widths.items() if w != expected for i in ii)
        warnings.warn(
            f"{path}: dropping {len(bad)} malformed line(s) (line numbers {bad[:5]}"
            f"{'...' if len(bad) > 5 else ''}); their sightlines are absent from the result. "
            f"A killed run loses this file's buffered tail.")
        kept = [ln for ln in lines if ln.strip() and len(ln.split()) == expected]
        path = io.StringIO("".join(kept))
        ncol = expected
    else:
        ncol = len(pd.read_csv(path, sep=r"\s+", header=None, nrows=1).columns)
    # Positional file with no header: the width identifies the vintage (no w_area/lon/lat; no
    # depths; current). Columns a file lacks come back NaN.
    full = len(MAP_COLS)
    if ncol not in (full - _N_DEPTH - _N_E1, full - _N_DEPTH, full):
        raise ValueError(f"{path}: {ncol} columns; expected {full - _N_DEPTH - _N_E1}, "
                         f"{full - _N_DEPTH} or {full}")
    df = pd.read_csv(path, sep=r"\s+", header=None, names=MAP_COLS[:ncol])
    for c in MAP_COLS[ncol:]:
        df[c] = np.nan
    return df


def load_provenance(path):
    """Parse run_provenance.txt into a dict of strings.

    The sightline-outcome block is appended at the end of a run, so its absence means the run
    did not finish.
    """
    prov = {}
    with open(path) as fh:
        for line in fh:
            m = re.match(r"^#\s+(\w+)\s+(.*?)\s*(?:#.*)?$", line)
            if m:
                prov[m.group(1)] = m.group(2).strip()
    return prov


# ---- Sentinel-aware accessors. Use these instead of touching the columns directly. ----

def sigma(df, param, survey, noise="W"):
    """1-sigma forecast for `param` from `survey`, NaN where it was not measured.

    param  : "tE" | "piE"    (photometric, gated on okA)
             "tetE"          (astrometric, gated on okB)
             "Ml"            (derived from tetE and piE, gated on its own positivity)
    survey : "joint" | "rubin" | "roman"
    noise  : astrometric noise variant for "tetE" and "Ml": "W" white (the main
             columns), "N" nominal, "P" pessimistic (columns sigtetE_N<q>, relMl_P<q>, gated on
             okB_N<q> ...). They exist for the joint and Roman partitions; Rubin's is the same
             in all three, so noise is ignored for "rubin". Photometric params ignore it.

    Gating is on the ok flag and on positivity: the flag can be set while an individual parameter
    is still a sentinel, because each partition fits its own active parameter subset
    (activePhotParams in include/fisher/fisher.h). E.g. an event Roman detects with no Rubin epochs
    has a valid joint fit in which the Rubin blend fraction was never free.
    """
    q = SURVEYS[survey]
    if noise not in ("W", "N", "P"):
        raise KeyError(f"unknown noise variant {noise!r}")
    v_ = "" if (noise == "W" or q == "L") else noise      # suffix: "", "N" or "P"
    if param in ("tE", "piE"):
        col, gate = f"sig{param}_{q}", f"okA_{q}"
    elif param == "tetE":
        col, gate = f"sigtetE_{v_}{q}", f"okB_{v_}{q}"
    elif param == "Ml":
        col, gate = f"relMl_{v_}{q}", None
    else:
        raise KeyError(f"unknown param {param!r}")

    v = df[col].astype(float)
    ok = (v > 0.0)
    if gate is not None:
        ok &= (df[gate] == 1)
    return v.where(ok)


def characterized(df, survey):
    """Abrams et al. 2025 characterization criterion: tE > 2*sigma_tE AND piE > 2*sigma_piE.

    Their criterion is used so Rubin-alone numbers are directly comparable to their published
    values; it is looser than sigma_tE/tE < 0.1 because two parameters are constrained at once.

    Returns a boolean Series; events where either sigma is unmeasured are False, never NaN.
    """
    s_tE = sigma(df, "tE", survey)
    s_piE = sigma(df, "piE", survey)
    return ((df["tE"] > 2.0 * s_tE) & (df["piE"] > 2.0 * s_piE)).fillna(False)


def ratio_joint_over(df, param, survey):
    """Per-event sigma_joint / sigma_<survey>, NaN unless BOTH were measured.

    Per event, then aggregate -- never a ratio of separately averaged sigmas. Adding data cannot
    worsen a Fisher forecast, so this is bounded above by 1 in exact arithmetic. A handful of events
    exceed 1 by ~1e-3, all with photometric condition number above 1e9, where double precision has
    lost most of its digits. A violation on a well-conditioned event would be a partitioning bug.
    """
    return sigma(df, param, "joint") / sigma(df, param, survey)


def detected(df, survey):
    """Boolean: did this survey's own detection test fire?"""
    return df[{"joint": "detJ", "rubin": "detL", "roman": "detR"}[survey]] == 1


def area_weight(df, prov=None):
    """Per-event sky-area weight in deg^2: the area of sky each row stands for.

    Sightlines inside Roman's footprint (~2.6% of the scanned region) are visited on a finer grid
    than those outside, so the sample is not proportional to sky area. Per-sightline quantities,
    statistics conditional on a selection (in-footprint events, per-field tables) and per-event
    ratios are unaffected. Anything pooled across the whole scan (survey-wide yield, histogram over
    all detections) must be weighted by this.

    Returns a float Series aligned to df. Falls back to the run's constant `area_per_sightline`
    for tables without a w_area column, and to 1.0 (plain counts) with no provenance.
    """
    if "w_area" in df.columns:
        return df["w_area"].astype(float)
    if prov and "area_per_sightline" in prov:
        return pd.Series(float(prov["area_per_sightline"]), index=df.index)
    return pd.Series(1.0, index=df.index)


def event_weight(df, sightlines, nsim_override=None):
    """Per-event importance weight for a pooled statistic.

    A pooled fraction or median is a statement about the sky only if the events are distributed
    like real events. The simulator draws

        Dl  with density proportional to rho(Dl) sqrt(Ds x(1-x)),  x = Dl/Ds
        Ml  from the Kroupa IMF and remnant map (number-weighted)
        v   from the component Gaussians (unweighted)
        u0, t0 uniform

    while the event rate carries an extra factor `R_E * v_t`:

        Gamma = int dDl dM d^2v  n(Dl) phi(M) f(v) * 2 u0m * R_E(M, Dl) * v_t

    Dividing the rate by the sampling density leaves, per drawn event,

        W = [w_area * Nstart / nsim] * sqrt(Ml) * Vt * Z(Ds)

    The bracket converts one draw into sky events at that sightline (area, stars per deg^2, draws
    taken); `sqrt(Ml) * Vt` is the rate weighting the sampler omits; `Z(Ds)` is the normaliser of
    the lens-distance sampler (galaxy_model.lens_distance_norm). Constants common to all events
    (2 u0m, Tobs, the Einstein-radius coefficient, the mean lens mass) cancel in any weighted
    fraction, median or ratio and are omitted: this is not an absolute yield.

    Limitation: the exact source-star weight carries 1/<m> for the source's Galactic component, but
    the event table stores the lens's component only, so the draw-average `Nstart` is used. That is
    unbiased between sightlines and drops a within-sightline factor of at most ~1.5.

    Per-event quantities (sigma_joint/sigma_single, paired comparisons) are already weight-free;
    weight only when pooling across events. Quote the Kish effective sample size
    (sum w)^2 / sum(w^2) beside any weighted number.

    df          : event table from load_events()
    sightlines  : map table from load_sightlines(), for `nsim` and the (lon, lat) join
    nsim_override : optional {(lon, lat) rounded to 3 dp: nsim}, for sightlines missing from a
                    damaged map file. The run log prints `nsim` for every sightline.

    Returns a float Series aligned to df. Raises if any event's sightline has no `nsim`.
    """
    import galaxy_model as G

    if "w_area" not in df.columns:
        raise ValueError("event table has no w_area column (older simulator version), so no pooled weight")

    codes, keys = sightline_index(df)
    nsim = {}
    if sightlines is not None and "lon" in sightlines:
        for lon, lat, n in zip(sightlines["lon"], sightlines["lat"], sightlines["nsim"]):
            if np.isfinite(lon):
                nsim[(round(float(lon), 3), round(float(lat), 3))] = float(n)
    if nsim_override:
        nsim.update({(round(k[0], 3), round(k[1], 3)): float(v)
                     for k, v in nsim_override.items()})

    # Sightlines with no draw count have two causes:
    #   barren: a sightline that drew stars but ended with no characterised event skips both the
    #     map-row write and the `nsim:` print (src/sim/sightline.cpp). Its rows carry no detection
    #     and no characterisation, so weight 0 changes no weighted statistic or count.
    #   killed run: the map stream is unflushed, so an interrupted run loses its buffered tail and
    #     real sightlines go missing. Those rows carry detections; weighting them 0 would delete
    #     part of the sky.
    # So weight 0 only rows that demonstrably contribute nothing; raise otherwise.
    missing = sorted(k for k in keys if k not in nsim)
    zero_weight = np.zeros(len(df), dtype=bool)
    if missing:
        miss = set(missing)
        zero_weight = np.isin(codes, [c for c, k in enumerate(keys) if k in miss])
        countable = [c for c in ("detL", "detR", "detJ", "okA_J", "okB_J") if c in df.columns]
        if countable:
            carried = int((df.loc[zero_weight, countable] == 1).any(axis=1).sum())
        else:
            # No detection columns (a paired-satellite file, every row a detection): nothing
            # can be shown to be empty, so nothing may be dropped.
            carried = int(zero_weight.sum())
        if carried:
            raise ValueError(
                f"{len(missing)} sightline(s) in the event table have no nsim, e.g. "
                f"{missing[:3]}, and {carried:,} of their rows carry a detection or a "
                f"characterisation. That is the signature of a killed run whose map file lost "
                f"its buffered tail, NOT of barren sightlines; weighting them "
                f"at zero would delete part of the sky. Recover nsim from the run log and pass "
                f"nsim_override.")
        print(f"  note: {zero_weight.sum():,} rows from {len(missing)} barren sightline(s) "
              f"carry no detection and get weight 0 (they have no nsim to normalise by)")
        for k in missing:
            nsim[k] = 1.0      # placeholder; these rows are zeroed at the end regardless

    # One density profile per sightline, not per event.
    n_draws = np.empty(len(df))
    nstart = np.empty(len(df))
    Z = np.empty(len(df))
    Ds = df["Ds"].to_numpy()
    for c, pos in sightline_groups(codes):
        k = keys[c]
        prof = G.density_profile(*k)
        n_draws[pos] = nsim[k]
        nstart[pos] = prof.Nstart
        Z[pos] = G.lens_distance_norm(prof, Ds[pos])

    w = (df["w_area"].to_numpy() * nstart / n_draws
         * np.sqrt(df["Ml"].to_numpy()) * df["Vt"].to_numpy() * Z)
    # Barren rows are zeroed after the formula, the only place a weight becomes 0.
    w[zero_weight] = 0.0
    return pd.Series(w, index=df.index)


def keep_weightable(map_path=None, log_paths=()):
    """A `keep` predicate for load_events that drops rows no statistic can use.

    A barren sightline (drew stars, characterised nothing) writes its rows but no map row and no
    `nsim:` line, so those rows have no draw count, get weight 0 (see event_weight) and carry no
    detection. At the scan's western edge such a sightline runs to the full --maxdraws cap, so
    they can be most of the table; dropping them is a memory saving, not a cut.

    Returns None when there is nothing to filter against, so a caller can pass the result
    straight through to load_events(keep=...) unconditionally.
    """
    known = set()
    if map_path and os.path.exists(map_path):
        sl = load_sightlines(map_path)
        if "lon" in sl:
            known |= {(round(float(a), 3), round(float(b), 3))
                      for a, b in zip(sl["lon"], sl["lat"]) if np.isfinite(a)}
    known |= set(nsim_from_logs(log_paths).keys())
    if not known:
        return None

    def keep(chunk):
        k = zip(chunk["lon"].round(3), chunk["lat"].round(3))
        return pd.Series([p in known for p in k], index=chunk.index)

    return keep


def nsim_from_logs(paths):
    """{(lon, lat): nsim} from the run log's per-sightline report.

    The map file is written unflushed, so a killed run loses its buffered tail and the
    sightlines it finished last have no map row. The log still has them.
    """
    out, lon, lat = {}, None, None
    for p in paths or ():
        with open(p, errors="replace") as fh:
            for line in fh:
                m = re.match(r"^longtitude:\s*(\S+)\s+latitude:\s*(\S+)", line)
                if m:
                    lon, lat = round(float(m.group(1)), 3), round(float(m.group(2)), 3)
                elif line.startswith("nsim:") and lon is not None:
                    out[(lon, lat)] = float(line.split()[1])
    return out


def attach_weight(df, map_path=None, log_paths=(), unweighted=False):
    """Return (weights, label) for a pooled statistic over `df`.

    The plumbing every figure script needs: read the sightline table, patch in `nsim` for
    sightlines a killed run's map file lost, and build the per-event weight.

    `unweighted=True` returns ones and says so in the label. A script must never fall back to
    unweighted silently: the unweighted sample over-represents long-tE events roughly tenfold.
    """
    if unweighted:
        return pd.Series(1.0, index=df.index), "unweighted"
    if not map_path:
        raise ValueError(
            "pooled statistics need the event-rate weight. Pass the map file "
            "(--map), adding --log for a run whose map file lost rows, or pass --unweighted "
            "to say deliberately that this figure is of the raw sample.")
    sl = load_sightlines(map_path)
    w = event_weight(df, sl, nsim_override=nsim_from_logs(log_paths))
    return w, f"event-rate weighted, N_eff = {kish_neff(w):,.0f}"


# ---- Absolute yields ----
#
# event_weight() is the event rate with every factor common to all events stripped out, because a
# fraction does not need them. An absolute yield does:
#
#   one draw's rate per source star  gamma_i = (F / <M>) * 2 u0m * kappa * Z(Ds) * sqrt(M) * v_t
#   expected detected events         N = T * sum_k Omega_k Nstar_k / nsim_k * sum_i gamma_i det_i
#                                      = T * 2 u0m * kappa * F * sum_i W_i det_i / <M>
#
# with kappa = sqrt(4 G Msun / c^2) and the unit conversions that make Z (Msun/pc^3 kpc^1.5) times
# sqrt(M) times v_t (km/s) a rate. <M> is the mean lens mass of the population as drawn (the rate
# is per lens, and F rho / <M> is the lens number density), taken per Galactic component because
# the `bulge` population's mass function differs between them.
#
# u0m = 3.0 (u0 is drawn uniform on [0.001, u0m]) and t0 is uniform on [2 d, Tobs - 2 d], read from
# config/parameters.h.
U0M = P.u0m
T0_MARGIN_DAYS = P.T0_MARGIN_DAYS
_G, _C, _MSUN = 6.67430e-11, 2.99792458e8, 1.98847e30
_PC = 3.0856775814913673e16                           # m
_KPC = 1.0e3 * _PC
RATE_UNIT = (2.0 * U0M * np.sqrt(4.0 * _G * _MSUN / _C**2)   # m^0.5
             * _KPC**1.5 / _PC**3                            # Z's units -> m^-1.5
             * 1.0e3)                                        # v_t km/s -> m/s
# RATE_UNIT * Z * sqrt(M) * Vt / <M> is gamma_i in s^-1 per source star (F = 1).

# Per-filter single-visit depth and saturation (ugrizy, F146) from config/parameters.h `thre` / `satu`.
# Roman's F146 entries are what the C++ uses. For Rubin the C++ uses per-sightline medians of the
# matched visits' 5-sigma depths (st.rubinDepthMed) with saturation = depth - RUBIN_SATU_BELOW_M5;
# the map file carries them (depth5_*) and acceptance_probability reads them. The fixed ugrizy
# entries here are SRD-style values, used only with fixed_rubin_depths=True (an approximation).
THRE = P.thre
SATU = P.satu
FILTERS = ["u", "g", "r", "i", "z", "y", "F146"]


def draw_rate(df):
    """gamma_i: the rate (s^-1 per source star, F = 1) that draw i stands for.

    Independent of event_weight() (Z is recomputed here per sightline), so the optical-depth
    check in y1_absolute_yield.py tests the constants rather than re-deriving them.
    """
    import galaxy_model as G
    Z = np.empty(len(df))
    Ds = df["Ds"].to_numpy()
    codes, keys = sightline_index(df)
    for c, pos in sightline_groups(codes):
        Z[pos] = G.lens_distance_norm(G.density_profile(*keys[c]), Ds[pos])
    return (RATE_UNIT * Z * np.sqrt(df["Ml"].to_numpy()) * df["Vt"].to_numpy()
            / mean_lens_mass(df))


def mean_lens_mass(df):
    """Per-row <M>: the mean of `Ml` over ALL draws of the same lens component (`struc`).

    Every draw is written to the table and the mass is drawn before any detection test, so a plain
    mean over rows is the sampler's own mean. Needs a table that still holds the undetected draws.
    """
    if not (df["detJ"] == 0).any():
        raise ValueError("mean_lens_mass needs the undetected draws too; this table has "
                         "only detections, whose masses are biased toward detectable ones")
    m = df.groupby("struc")["Ml"].transform("mean")
    return m.to_numpy()


def acceptance_probability(df, sightlines, fixed_rubin_depths=False):
    """P that the simulator kept this draw, rebuilt from the table (preselectEvent in src/sim/draw.cpp).

    A drawn star is kept for light-curve generation if Rubin could see its peak (Mpeak below the
    depth and baseline above saturation in >= 2 of ugrizy) AND a uniform draw falls below its r-band
    blend fraction, or likewise for Roman in F146. That thinning follows Sajadian & Makler
    (criterion ii): the blend fraction is the probability of "realising" one star of an unresolved
    blend, which counts events per resolved object. Both light curves are generated whenever either
    survey accepts, so detection is independent of which acceptance fired and 1/P undoes the
    thinning exactly.

    Rubin's depth is the sightline's median 5-sigma depth per band, from the map file's depth5_*
    columns (`sightlines` = load_sightlines()), matched to each row on (lon, lat) to 3 dp; saturation
    is that depth - RUBIN_SATU_BELOW_M5. A band with no visit has depth -inf and cannot be seen. A row
    whose sightline has no map row, or whose map has no depths (NaN), raises ValueError.
    fixed_rubin_depths=True instead uses the fixed THRE/SATU for Rubin, an approximation for runs
    whose map file has no depths; it can disagree with the C++ for any draw whose peak lies between
    the fixed and the per-sightline depth.
    """
    u0 = df["u0"].to_numpy()
    A0 = (u0**2 + 2.0) / np.sqrt(u0**2 * (u0**2 + 4.0))
    if fixed_rubin_depths:
        depth = np.broadcast_to(np.asarray(THRE[:6], dtype=float), (len(df), 6))
        satu = np.broadcast_to(np.asarray(SATU[:6], dtype=float), (len(df), 6))
    else:
        codes, keys = sightline_index(df)
        cols = [f"depth5_{f}" for f in FILTERS[:6]]
        table = {(round(a, 3), round(b, 3)): d for a, b, d in
                 zip(sightlines["lon"], sightlines["lat"], sightlines[cols].to_numpy())}
        missing = [k for k in keys if k not in table]
        if missing:
            raise ValueError(f"{len(missing)} sightline(s) of the table have no row in the map "
                             f"file (first: {missing[0]}); their Rubin depths are unknown")
        per_sl = np.array([table[k] for k in keys])
        if np.isnan(per_sl).any():
            raise ValueError("the map file has no per-sightline Rubin depths (an older "
                             "simulator version); rerun, or pass fixed_rubin_depths=True for the "
                             "approximate fixed-depth acceptance")
        depth = per_sl[codes]
        satu = depth - P.RUBIN_SATU_BELOW_M5
    seen = np.empty((len(df), 7), dtype=bool)
    for i, f in enumerate(FILTERS):
        mb, fb = df[f"magb_{f}"].to_numpy(), df[f"blend_{f}"].to_numpy()
        mpeak = mb - 2.5 * np.log10(A0 * fb + 1.0 - fb)
        if i < 6:
            seen[:, i] = (mpeak <= depth[:, i]) & (mb > satu[:, i])
        else:
            seen[:, i] = (mpeak <= THRE[i]) & (mb > SATU[i])
    pL = np.where(seen[:, :6].sum(axis=1) > 1, df["blend_r"].to_numpy(), 0.0)
    pR = np.where(seen[:, 6], df["blend_F146"].to_numpy(), 0.0)
    return 1.0 - (1.0 - pL) * (1.0 - pR)


def yield_weight(df, sightlines, tobs_days, nsim_override=None):
    """Per-row expected number of events in the window, for F = 1, per resolved object.

    y_i = T * [w_area Nstar / nsim] * gamma_i. Summing y_i over detected rows gives the
    detected yield with every lens drawn from the population and F = 1; multiply by F.
    Divide each y_i by acceptance_probability() to count events on every star instead.
    """
    w = event_weight(df, sightlines, nsim_override=nsim_override).to_numpy()
    # event_weight = w_area Nstart/nsim sqrt(M) Vt Z; the rest of gamma is RATE_UNIT / <M>.
    T = (tobs_days - 2.0 * T0_MARGIN_DAYS) * 86400.0
    y = T * RATE_UNIT * w / mean_lens_mass(df)
    return y


def weighted_quantile(values, weights, q):
    """Quantile `q` of `values` under `weights`, by cumulative weight. NaN if empty."""
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    ok = np.isfinite(v) & np.isfinite(w) & (w > 0)
    v, w = v[ok], w[ok]
    if v.size == 0:
        return np.nan
    o = np.argsort(v)
    c = np.cumsum(w[o]) / w[o].sum()
    return float(v[o][np.searchsorted(c, q)])


def weighted_median(values, weights):
    return weighted_quantile(values, weights, 0.5)


def weighted_fraction(mask, weights):
    """Weighted fraction of `mask` being true, over the whole weighted sample."""
    m = np.asarray(mask.fillna(False) if hasattr(mask, "fillna") else mask, dtype=bool)
    w = np.asarray(weights, dtype=float)
    tot = w.sum()
    return float(w[m].sum() / tot) if tot > 0 else np.nan


def kish_neff(w):
    """Effective sample size of a weighted sample: (sum w)^2 / sum(w^2).

    A weighted fraction over N events with a spread of weights carries the Poisson precision of
    this many unweighted ones. Report it whenever you report a weighted number.
    """
    w = np.asarray(w, dtype=float)
    return float(w.sum()**2 / np.square(w).sum()) if w.size else 0.0


def is_stratified(prov):
    """True if the run used --stride-roman, i.e. the sightlines stand for unequal sky areas.

    Scripts that quote an absolute yield from such a run must apply area_weight().
    """
    return bool(prov) and str(prov.get("stratified", "0")).strip() == "1"


def check_monotonicity(df, params=("tE", "piE", "tetE", "Ml"), tol=1e-9):
    """Assert sigma_joint <= sigma_single wherever both exist. Returns a list of violations.

    The joint Fisher matrix is the sum of the per-survey ones, so adding data cannot increase a
    forecast error; a violation indicates a partitioning bug.
    """
    bad = []
    for p in params:
        for s in ("rubin", "roman"):
            r = ratio_joint_over(df, p, s).dropna()
            if len(r) and r.max() > 1.0 + tol:
                bad.append((p, s, float(r.max()), int((r > 1.0 + tol).sum())))
    return bad


# Where the simulator writes run_provenance.txt, and where a copy is sometimes kept beside an
# archived run. Searched in order.
PROVENANCE_SEARCH = ("run_provenance.txt",
                     "files/MONTLMC/files/run_provenance.txt")


# Roman's per-exposure astrometric error, mirrored from src/surveys/noise.cpp errRomanA. The anchors
# are F146 Vega magnitudes (Lam et al. 2026); the simulator's magnitudes are AB, converted by
# m_Vega = m_AB - 1.0324. Runs that convert say so in their provenance ("# roman_noise"); runs
# without that line used the AB magnitude as if Vega, and analyses of them must do the same.
F146_AB_MINUS_VEGA = 1.0324
ROMAN_AST = dict(floor=1.1, mflr=20.62, mbkg=23.5, sbkg=10.0, slope_src=0.33285, slope_bkg=0.4)


def roman_ast_vega_offset(prov_path):
    """AB - Vega offset the run that wrote `prov_path` applied before errRomanA (0 if it did not convert)."""
    if prov_path and os.path.exists(prov_path):
        if "# roman_noise" in open(prov_path).read():
            return F146_AB_MINUS_VEGA
    return 0.0


def roman_ast_error(mag_ab, ab_minus_vega=F146_AB_MINUS_VEGA):
    """errRomanA, vectorised: per-exposure, per-coordinate astrometric error [mas] at F146 AB."""
    a = ROMAN_AST
    m = np.asarray(mag_ab, dtype=float) - ab_minus_vega
    out = np.where(m > a["mbkg"], a["sbkg"] * 10.0 ** (a["slope_bkg"] * (m - a["mbkg"])),
                   np.where(m > a["mflr"], a["floor"] * 10.0 ** (a["slope_src"] * (m - a["mflr"])),
                            a["floor"]))
    return np.maximum(out, a["floor"])


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


def population(prov):
    """Which lens population produced a run, from its provenance.

    Provenance without a `population` line is the Kroupa-plus-remnants bulge population (the only
    one before the line was added), reported as "bulge (implied)".
    """
    if not prov:
        return "unknown"
    return prov.get("population", "bulge (implied)")


def assert_same_population(provs, what="this figure"):
    """Refuse to pool runs from different lens populations.

    The event-rate weight carries a sqrt(Ml) factor that is only correct for the mass function
    actually sampled, so two populations in one pooled statistic give a wrong number.
    Comparing them side by side (one curve per population) is not pooling and does not come here.
    """
    names = {population(p) for p in provs if p}
    if len(names) > 1:
        raise ValueError(
            f"{what} mixes lens populations {sorted(names)}. Each population has its own mass "
            f"function, and the pooled weight is only valid for the one that was sampled. "
            f"Plot them as separate series instead of pooling them.")
    return names.pop() if names else "unknown"


def describe(path_events, path_prov=None):
    """One-line provenance summary to print at the top of every figure-producing script."""
    parts = [f"events={os.path.basename(path_events)}"]
    path_prov = find_provenance(path_prov, near=path_events)
    if path_prov and os.path.exists(path_prov):
        prov = load_provenance(path_prov)
        parts.append(f"population={population(prov)}")
        for k in ("git_commit", "stride", "events_target", "sightlines_aggregated"):
            if k in prov:
                parts.append(f"{k}={prov[k]}")
        if is_stratified(prov):
            # Loud: absolute yields from a stratified run are wrong by the oversampling factor
            # unless weighted.
            parts.append(f"STRATIFIED(stride_roman={prov.get('stride_roman', '?')};"
                         f" weight by w_area)")
        if "sightlines_aggregated" not in prov:
            parts.append("INCOMPLETE-RUN")
    else:
        parts.append("provenance=NOT FOUND")
    return "  ".join(parts)
