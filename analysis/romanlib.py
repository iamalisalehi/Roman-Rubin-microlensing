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

import io
import os
import re
import warnings

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
               "log10_Rostart", "log10_Nstart", "log10_nstart",
               # Step E1, appended in this order. w_area is the deg^2 of sky this
               # sightline stands for -- constant across an unstratified run, NOT constant
               # once --stride-roman is used. lon/lat are the sightline's position, which
               # this file did not previously record at all: without them a map row could
               # not be tied to the events it produced, and `nsim` -- the draw count any
               # pooled yield needs as its denominator -- was unreachable from the event
               # table. Files written before Step E1 have none of the three;
               # load_sightlines() detects that by width.
               "w_area", "lon", "lat"])


def load_events(path, keep=None, chunksize=None, usecols=None):
    """Read the per-event table (test5.dat) written by the `filg_in <<` block.

    Column names come from the file's own `#` header, not from a list hardcoded here, so a
    schema change surfaces as a loud mismatch rather than a silent misalignment.

    keep, chunksize -- read in chunks and keep only the rows `keep(chunk)` selects.
        The production table is 5.57M rows x 90 float64 columns, which is ~4 GB resident
        before pandas' parse buffers are counted. On a machine with less than about 12 GB
        that is an OOM kill, not a slow read -- and the kill is silent, exit status 0 with
        an empty stdout, which looks exactly like a script that did nothing (OPEN_ITEMS.md).
        Filtering per chunk holds the peak at one chunk plus the surviving rows.

        `keep` is called with each chunk and must return a boolean mask over it. Every
        analysis here begins by discarding undetected events -- 98.7% of the table -- so:

            df = R.load_events(path, keep=lambda c: c["detJ"] == 1, chunksize=500_000)

        The default (keep=None) reads the whole file in one pass, unchanged, so existing
        callers behave exactly as before.

    usecols -- keep only these columns. The other lever on the same problem: a statistic over
        ALL draws (the intrinsic tE distribution, Deviation 41) cannot discard 98.7% of the
        rows, so it has to discard columns instead -- eight of ninety is ~380 MB rather than
        ~4 GB. The `detCls`/`synClass`/`t0zone` label columns are only added if their source
        column survives the selection.
    """
    with open(path) as fh:
        header = fh.readline()
    if not header.startswith("#"):
        raise ValueError(
            f"{path}: no '#' header line. Files written before Step D1 have no header and a "
            f"different column set; they cannot be read with this loader."
        )
    cols = header.lstrip("#").split()

    reader_kw = dict(sep=r"\s+", comment="#", header=None, names=cols)
    if usecols is not None:
        missing = [c for c in usecols if c not in cols]
        if missing:
            raise ValueError(f"{path}: no such column(s) {missing}")
        reader_kw["usecols"] = list(usecols)
    if keep is None and chunksize is None:
        df = pd.read_csv(path, **reader_kw)
    else:
        parts = []
        with pd.read_csv(path, chunksize=chunksize or 500_000, **reader_kw) as it:
            for chunk in it:
                parts.append(chunk if keep is None else chunk[keep(chunk)])
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

    Note this file is opened in APPEND mode by the simulator, so a re-run without clearing
    it first silently concatenates two runs (OPEN_ITEMS.md). `nri`/`nde` restarting from
    zero part-way down the file is the signature.

    Malformed lines are dropped with a warning rather than killing the read. The simulator
    never flushes this stream, so a killed run loses its buffered tail and the next run's
    first row lands on the fragment, leaving one line of the wrong width (OPEN_ITEMS.md, and
    the v3 file has exactly that at line 687). Silently mis-parsing it would shift every
    column; failing outright would make the whole file unreadable for one bad line.
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
            f"A killed run loses this file's buffered tail -- see OPEN_ITEMS.md.")
        kept = [ln for ln in lines if ln.strip() and len(ln.split()) == expected]
        path = io.StringIO("".join(kept))
        ncol = expected
    else:
        ncol = len(pd.read_csv(path, sep=r"\s+", header=None, nrows=1).columns)
    if ncol == len(MAP_COLS) - 3:
        # Written before Step E1 added the per-sightline area weight. Positional file with
        # no header, so the only way to tell is the width -- and guessing wrong shifts every
        # column by one, which produces a plausible plot of the wrong quantity.
        df = pd.read_csv(path, sep=r"\s+", header=None, names=MAP_COLS[:-3])
        for c in ("w_area", "lon", "lat"):
            df[c] = np.nan
        return df
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


def area_weight(df, prov=None):
    """Per-event sky-area weight in deg^2 -- the area of sky each row stands for (Step E1).

    WHEN YOU NEED THIS. Roman's GBTDS footprint is ~2.6% of the scanned region, so a uniform
    scan spends 97% of its draws where the joint Fisher matrix is simply Rubin's. Step E1
    stratifies the scan: sightlines inside the footprint are visited on a finer grid than
    those outside, and the sample is then deliberately NOT proportional to sky area.

    That changes nothing computed at a single sightline, and nothing conditional on a
    selection you make yourself. A statistic over in-footprint events only (F3 panel (a), the
    F4 footprint panels), or a per-field table (F1), or a per-event ratio -- none of them
    move, because the sightline you drew from is not an input to any of them.

    What DOES move is anything pooled across the whole scan: a survey-wide yield, a
    histogram over all joint detections, the "all detections" check panel of F4. Weight those
    by this, or they will describe a sky in which Roman covers whatever fraction of the
    SAMPLE the stratification bought rather than the 2.6% of the SKY it actually covers.

    Returns a float Series aligned to df. Falls back to the run's constant
    `area_per_sightline` for tables written before Step E1, and to 1.0 (with everything
    equally weighted, i.e. plain counts) if there is no provenance to fall back to.
    """
    if "w_area" in df.columns:
        return df["w_area"].astype(float)
    if prov and "area_per_sightline" in prov:
        return pd.Series(float(prov["area_per_sightline"]), index=df.index)
    return pd.Series(1.0, index=df.index)


def event_weight(df, sightlines, nsim_override=None):
    """Per-event importance weight for a POOLED statistic (Step W1, Deviation 41).

    WHY A WEIGHT IS NEEDED AT ALL. A pooled fraction or median is only a statement about the sky
    if the events it averages are distributed like real events. They are not. The simulator draws

        Dl  with density proportional to rho(Dl) sqrt(Ds x(1-x)),  x = Dl/Ds
        Ml  from the Kroupa IMF and remnant map -- NUMBER-weighted
        v   from the component Gaussians -- unweighted
        u0, t0 uniform

    while the event rate carries a factor `R_E * v_t` on top of the population:

        Gamma = int dDl dM d^2v  n(Dl) phi(M) f(v) * 2 u0m * R_E(M, Dl) * v_t

    Dividing the rate by the sampling density cancels rho, phi, f and the distance factor and
    leaves, per drawn event,

        W = [w_area * Nstart / nsim] * sqrt(Ml) * Vt * Z(Ds)

    The bracket converts one draw into sky events at that sightline (area, stars per deg^2, draws
    taken); `sqrt(Ml) * Vt` is the rate weighting the sampler omits; `Z(Ds)` is the normaliser of
    the lens-distance sampler (galaxy_model.lens_distance_norm). Constants -- 2 u0m, Tobs, the
    Einstein-radius coefficient, the mean lens mass -- are identical for every event and cancel in
    any weighted fraction, median or ratio, so they are not included: THIS IS NOT AN ABSOLUTE YIELD.

    WHAT IT DOES NOT FIX. The exact source-star weight carries 1/<m> for the SOURCE's Galactic
    component; the event table stores the lens's component only, so `Nstart` (the draw-average of
    the same thing) is used instead. That is unbiased between sightlines and drops a within-
    sightline factor of at most ~1.5.

    WHEN NOT TO USE IT. Anything computed per event -- sigma_joint/sigma_single for one event,
    H3's paired comparison -- is already weight-free. Weight only when pooling ACROSS events.
    Quote the Kish effective sample size (sum w)^2 / sum(w^2) beside any weighted number: on the
    v3 table the weight takes 8,894 footprint events to an effective 2,564.

    df          : event table from load_events()
    sightlines  : map table from load_sightlines(), for `nsim` and the (lon, lat) join
    nsim_override : optional {(lon, lat) rounded to 3 dp: nsim}, for sightlines missing from a
                    damaged map file. The run log prints `nsim` for every sightline.

    Returns a float Series aligned to df. Raises if any event's sightline has no `nsim`.
    """
    import galaxy_model as G

    if "w_area" not in df.columns:
        raise ValueError("event table predates Step E1: no w_area column, so no pooled weight")

    key = list(zip(df["lon"].round(3), df["lat"].round(3)))
    nsim = {}
    if sightlines is not None and "lon" in sightlines:
        for lon, lat, n in zip(sightlines["lon"], sightlines["lat"], sightlines["nsim"]):
            if np.isfinite(lon):
                nsim[(round(float(lon), 3), round(float(lat), 3))] = float(n)
    if nsim_override:
        nsim.update({(round(k[0], 3), round(k[1], 3)): float(v)
                     for k, v in nsim_override.items()})

    missing = sorted({k for k in set(key) if k not in nsim})
    if missing:
        raise ValueError(
            f"{len(missing)} sightline(s) in the event table have no nsim, e.g. {missing[:3]}. "
            f"The map file is written without flushing and loses its tail when a run is killed "
            f"(OPEN_ITEMS.md); recover nsim from the run log and pass nsim_override.")

    # One density profile per sightline, not per event: ~9,500 grid points each.
    n_draws = np.empty(len(df))
    nstart = np.empty(len(df))
    Z = np.empty(len(df))
    Ds = df["Ds"].to_numpy()
    order = pd.Series(np.arange(len(df))).groupby(pd.Series(key)).groups
    for k, pos in order.items():
        pos = np.asarray(pos)
        prof = G.density_profile(*k)
        n_draws[pos] = nsim[k]
        nstart[pos] = prof.Nstart
        Z[pos] = G.lens_distance_norm(prof, Ds[pos])

    w = (df["w_area"].to_numpy() * nstart / n_draws
         * np.sqrt(df["Ml"].to_numpy()) * df["Vt"].to_numpy() * Z)
    return pd.Series(w, index=df.index)


def nsim_from_logs(paths):
    """{(lon, lat): nsim} from the run log's per-sightline report.

    The map file is written unflushed, so a killed run loses its buffered tail and the
    sightlines it finished last have no map row (Deviation 42). The log still has them.
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

    `unweighted=True` returns ones and says so in the label -- an explicit choice, which is
    the point. A script must never fall back to unweighted silently: the unweighted sample
    over-represents long-tE events roughly tenfold (Deviation 41), so a figure that quietly
    dropped the weight would look finished and be wrong.
    """
    if unweighted:
        return pd.Series(1.0, index=df.index), "unweighted"
    if not map_path:
        raise ValueError(
            "pooled statistics need the event-rate weight (Deviation 41). Pass the map file "
            "(--map), adding --log for a run whose map file lost rows, or pass --unweighted "
            "to say deliberately that this figure is of the raw sample.")
    sl = load_sightlines(map_path)
    w = event_weight(df, sl, nsim_override=nsim_from_logs(log_paths))
    return w, f"event-rate weighted, N_eff = {kish_neff(w):,.0f}"


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

    A run where this is True and a downstream script ignores area_weight() is a run whose
    absolute yields are wrong -- so scripts that quote one should check it and say so.
    """
    return bool(prov) and str(prov.get("stratified", "0")).strip() == "1"


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


def population(prov):
    """Which lens population produced a run, from its provenance.

    Runs made before 2026-09-17 have no `population` line because there was only one: the
    Kroupa-plus-remnants bulge population, then selected at compile time. Reporting them as
    "bulge (implied)" is accurate -- that is what they are -- while still distinguishing them
    from a run that says so itself.
    """
    if not prov:
        return "unknown"
    return prov.get("population", "bulge (implied)")


def assert_same_population(provs, what="this figure"):
    """Refuse to pool runs from different lens populations.

    The event-rate weight carries a sqrt(Ml) factor that is only correct for the mass function
    actually sampled (Deviation 45). Two populations in one pooled statistic is therefore not a
    presentation choice but a wrong number, and it is an easy mistake to make once the files
    are named testbh.dat and testns.dat and differ by two characters.

    Comparing them side by side -- one curve per population -- is fine and is the point of the
    population figures; that is not pooling, and does not come through here.
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
        # First, because it is the thing that makes two otherwise identical tables mean
        # different things.
        parts.append(f"population={population(prov)}")
        for k in ("git_commit", "stride", "events_target", "sightlines_aggregated"):
            if k in prov:
                parts.append(f"{k}={prov[k]}")
        if is_stratified(prov):
            # Loud, because a stratified run whose absolute yields are quoted unweighted is
            # wrong by the oversampling factor and looks entirely normal.
            parts.append(f"STRATIFIED(stride_roman={prov.get('stride_roman', '?')};"
                         f" weight by w_area)")
        if "sightlines_aggregated" not in prov:
            parts.append("INCOMPLETE-RUN")
    else:
        parts.append("provenance=NOT FOUND")
    return "  ".join(parts)
