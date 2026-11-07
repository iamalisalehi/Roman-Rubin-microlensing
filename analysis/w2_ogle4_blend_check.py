#!/usr/bin/env python3
"""Observational check of the simulated blend fractions against OGLE-IV (design section 3.6).

    .roman/bin/python analysis/w2_ogle4_blend_check.py EVENTS --map MAP --ogle TABLE3 [--log LOG ...] [-o out.png]

EVENTS is a bulge-population event table (test5.dat), MAP its sightline table (MapLMC5.dat), TABLE3
Table 3 of Mroz et al. (2019, ApJS 244, 29): 5,790 OGLE-IV events in the low-cadence fields, read
by the byte-by-byte description in the file's own header.

Simulated sample: Rubin-detected events (detL == 1) with source I-band magnitude
I_s = magb_i - 2.5 log10(blend_i) - (i_AB - I_Vega) <= 21, u0 <= 1, and Galactic latitude within the
OGLE sample's range. The simulated i magnitude is AB (LSST i) and OGLE's Is_best is Vega, so the
simulation is converted with romanlib.I_AB_MINUS_VEGA (0.4 mag) before the cut and the magnitude
bins. OGLE sample: I_s <= 21, u0 <= 1 (best fit), raw (unweighted) counts, f_s > 1 in the top bin.
The f_S distribution (bins < 0.2, 0.2-0.6, 0.6-0.9, >= 0.9) is printed overall and in source
magnitude bins 17-19, 19-20, 20-21, for the simulation in both counting conventions (every-star
first, the headline; N and Kish N_eff) beside OGLE.

Caveats (the check is of shape and trend, and of which convention lies closer):
  - the simulated blend disc uses Rubin's i seeing (FWHM 0.95") not OGLE's (~1.2-1.4"), so the
    simulation should come out less blended at fixed source magnitude;
  - the two surveys use different detection pipelines;
  - OGLE's fitted f_s carries fit noise and "negative blending" (f_s > 1).
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R          # noqa: E402

CONVS = ("every-star", "per-object")
FS_EDGES = (0.2, 0.6, 0.9)
FS_LABELS = ("< 0.2", "0.2 - 0.6", "0.6 - 0.9", ">= 0.9")
MAG_BINS = ((17.0, 19.0), (19.0, 20.0), (20.0, 21.0))
I_S_MAX = 21.0          # published cuts, inclusive: I_s <= 21, u0 <= 1
U0_MAX = 1.0


def read_ogle_table3(path):
    """Mroz et al. 2019 Table 3 as a DataFrame, columns placed by the file's byte-by-byte description."""
    spec = []
    with open(path) as fh:
        for line in fh:
            m = re.match(r"#\s+(\d+)-\s*(\d+)\s+\S+\s+(\S+)\s", line)
            if m:
                spec.append((m.group(3), int(m.group(1)) - 1, int(m.group(2))))
    need = ["fs_best", "Is_best", "u0_best", "glat"]
    have = {s[0] for s in spec}
    if not set(need) <= have:
        raise ValueError(f"{path}: byte description lacks {sorted(set(need) - have)}")
    keep = [s for s in spec if s[0] in need]
    df = pd.read_fwf(path, colspecs=[(a, b) for _, a, b in keep], names=[n for n, _, _ in keep],
                     comment="#", header=None)
    if df[need].isna().any().any():
        raise ValueError(f"{path}: unparsed values in {need}")
    return df


def fs_fractions(fs, w):
    """Weighted fraction of fs in each bin (the top bin is open: fs >= 0.9, including fs > 1)."""
    code = np.digitize(np.asarray(fs, dtype=float), FS_EDGES)
    w = np.asarray(w, dtype=float)
    tot = w.sum()
    return np.array([w[code == k].sum() / tot if tot > 0 else np.nan for k in range(4)])


def weighted_cdf(fs, w, grid):
    fs, w = np.asarray(fs, dtype=float), np.asarray(w, dtype=float)
    o = np.argsort(fs)
    c = np.cumsum(w[o]) / w.sum()
    return np.concatenate([[0.0], c])[np.searchsorted(fs[o], grid, side="right")]


def subsets(mag):
    yield "overall (I_s <= 21)", np.ones(len(mag), dtype=bool)
    for lo, hi in MAG_BINS:
        yield f"I_s {lo:g} - {hi:g}", (mag >= lo) & (mag < hi)


def render(sim, ogle, weights):
    smag, sfs = sim["Is"].to_numpy(), sim["fs"].to_numpy()
    omag, ofs = ogle["Is_best"].to_numpy(), ogle["fs_best"].to_numpy()
    L = []
    for (name, sm), (_, om) in zip(subsets(smag), subsets(omag)):
        cols, head, ns = [], ["f_S bin"], []
        for conv in CONVS:
            w = weights[conv][sm]
            cols.append(fs_fractions(sfs[sm], w) if w.sum() > 0 else np.full(4, np.nan))
            head.append(f"sim {conv}")
            ns.append((int(sm.sum()), R.kish_neff(w[w > 0]) if sm.any() else 0.0))
        cols.append(fs_fractions(ofs[om], np.ones(int(om.sum()))) if om.any() else np.full(4, np.nan))
        head.append("OGLE-IV")
        L += [f"### {name}", "",
              "| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
        for k, lab in enumerate(FS_LABELS):
            L.append(f"| {lab} | " + " | ".join("n/a" if np.isnan(c[k]) else f"{c[k]:.3f}" for c in cols) + " |")
        L.append("| N | " + " | ".join(f"{n[0]:,}" for n in ns) + f" | {int(om.sum()):,} |")
        L.append("| N_eff | " + " | ".join(f"{n[1]:,.0f}" for n in ns) + f" | {int(om.sum()):,} |")
        L.append("")
        for c in cols:
            if not np.isnan(c).all() and abs(c.sum() - 1.0) > 1e-9:
                raise RuntimeError(f"{name}: f_S bin fractions sum to {c.sum()}, not 1")
    return L


def make_figure(sim, ogle, weights, path):
    import plotstyle as ps
    ps.use_paper_style()
    import matplotlib.pyplot as plt
    grid = np.linspace(0.0, 1.0, 401)
    fig, axes = plt.subplots(1, 4, figsize=(ps.WIDTH["double"], 2.4), sharey=True)
    colors = {"every-star": ps.SURVEY["joint"], "per-object": ps.SURVEY["roman"]}
    sets = list(zip(subsets(sim["Is"].to_numpy()), subsets(ogle["Is_best"].to_numpy())))
    for ax, ((name, sm), (_, om)) in zip(axes, sets):
        for conv in CONVS:
            w = weights[conv][sm]
            if w.sum() > 0:
                ax.plot(grid, weighted_cdf(sim["fs"].to_numpy()[sm], w, grid), color=colors[conv],
                        label=f"simulation, {conv}")
        if om.any():
            ax.plot(grid, weighted_cdf(ogle["fs_best"].to_numpy()[om], np.ones(int(om.sum())), grid),
                    color="k", ls="--", label="OGLE-IV")
        ax.set_xlabel("$f_S$ (LSST $i$ / OGLE $I$)")
        ax.set_title(name, fontsize=7)
    axes[0].set_ylabel("cumulative fraction")
    axes[0].legend(fontsize=6)
    fig.tight_layout()
    ps.save_figure(fig, os.path.splitext(path)[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("events")
    ap.add_argument("--map", default=None, help="the run's sightline table (event-rate weights)")
    ap.add_argument("--log", action="append", default=[], help="run log(s), when the map file lost rows")
    ap.add_argument("--unweighted", action="store_true", help="deliberately use the raw simulated sample")
    ap.add_argument("--ogle", required=True, help="Mroz et al. 2019 Table 3 (ogle4_table3.dat)")
    ap.add_argument("--provenance", default=None, help="the run's run_provenance.txt (default: found beside MAP or EVENTS)")
    ap.add_argument("--chunksize", type=int, default=500_000)
    ap.add_argument("-o", "--out", default=None, help="write the CDF figure here (PNG; a PDF beside it)")
    a = ap.parse_args()

    pv = (R.find_provenance(a.provenance, near=a.map or a.events)
          or R.find_provenance(None, near=a.events))
    if pv is None:
        sys.exit(f"{a.events}: no run_provenance.txt found, so the lens population cannot be checked; "
                 "pass --provenance PATH")
    pop = R.population(R.load_provenance(pv))
    if not pop.startswith("bulge"):
        sys.exit(f"{a.events}: population is {pop!r}; the OGLE-IV events are bulge lenses, so this check "
                 "needs a bulge-population run")

    ogle_all = read_ogle_table3(a.ogle)
    lat_lo, lat_hi = float(ogle_all["glat"].min()), float(ogle_all["glat"].max())
    ogle = ogle_all[(ogle_all["Is_best"] <= I_S_MAX) & (ogle_all["u0_best"] <= U0_MAX)].reset_index(drop=True)

    df = R.load_events(a.events, keep=lambda c: c["detL"] == 1, chunksize=a.chunksize).reset_index(drop=True)
    if df.empty:
        sys.exit(f"{a.events}: no Rubin-detected events (detL == 1)")
    weights, labels = {}, []
    for conv in CONVS:
        w, lab = R.attach_weight(df, a.map, a.log, a.unweighted, conv)
        weights[conv] = np.asarray(w, dtype=float)
        labels.append(lab)
    n_det = len(df)
    sim_all = pd.DataFrame({"Is": df["magb_i"] - 2.5 * np.log10(df["blend_i"]) - R.I_AB_MINUS_VEGA, "fs": df["blend_i"],
                            "u0": df["u0"], "lat": df["lat"]})
    ok = ((sim_all["Is"] <= I_S_MAX) & (sim_all["u0"] <= U0_MAX)
          & (sim_all["lat"] >= lat_lo) & (sim_all["lat"] <= lat_hi)).to_numpy()
    sim = sim_all[ok].reset_index(drop=True)
    weights = {c: w[ok] for c, w in weights.items()}
    if sim.empty:
        sys.exit(f"{a.events}: no simulated events pass the OGLE-like cuts")

    head = ["# OGLE-IV blend-fraction check", "",
            f"Simulation: {a.events}. Weights: " + "; ".join(f"{c}: {l}" for c, l in zip(CONVS, labels)) + ".",
            f"Simulated Rubin-detected events {n_det:,}; after cuts (I_s <= {I_S_MAX:g}, u0 <= {U0_MAX:g}, "
            f"{lat_lo:.2f} <= b <= {lat_hi:.2f} deg) {len(sim):,}.",
            f"Population: {pop}. Simulated I_s is Vega: LSST i (AB) minus {R.I_AB_MINUS_VEGA:g} mag "
            "(romanlib.I_AB_MINUS_VEGA); OGLE's Is_best is Vega.",
            f"OGLE-IV Table 3: {len(ogle_all):,} events; after cuts (I_s <= {I_S_MAX:g}, u0 <= {U0_MAX:g}) "
            f"{len(ogle):,}; raw (unweighted) counts; f_s > 1 counted in the top bin.", "",
            "Caveats: the simulated blend disc uses Rubin's i seeing (FWHM 0.95 arcsec), not OGLE's (~1.2-1.4 arcsec), so",
            "the simulation should be less blended at fixed source magnitude; the detection pipelines differ;",
            "OGLE's fitted f_s carries fit noise and negative blending (f_s > 1). The check is of shape and trend,",
            "and of which counting convention lies closer.", ""]
    print("\n".join(head + render(sim, ogle, weights)))
    if a.out:
        make_figure(sim, ogle, weights, a.out)


if __name__ == "__main__":
    main()
