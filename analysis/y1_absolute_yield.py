#!/usr/bin/env python3
"""Step Y: absolute event yields, and what the raw Monte Carlo counts mean (Deviation 54).

WHY THIS EXISTS. Every number the analysis layer has produced so far is a fraction, median or
ratio, in which the rate's overall constants cancel (Deviation 41). The raw counts printed in
the report -- "110,144 detected events" for the black-hole run -- are Monte Carlo SAMPLE SIZES,
set by the detections each sightline was told to collect (--events/--lenses) and capped at
--maxdraws; with every lens drawn from the chosen population, they say nothing about the sky.
This script restores the constants and turns the same draws into an expected number of events.

THE FORMULA (full derivation in DEVIATIONS.md 54 and the report's yield section). One draw i,
with source at Ds, lens mass M and transverse velocity v_t, stands for a rate per source star

    gamma_i = (F / <M>) * 2 u0m * kappa * Z(Ds) * sqrt(M) * v_t,     kappa = sqrt(4 G / c^2)

(F: the population's fraction of the Galaxy's stellar mass; <M>: its mean mass as drawn;
Z: the lens-distance sampler's normaliser). Summing over sightlines k,

    N_det = T * sum_k Omega_k Nstar_k / nsim_k * sum_{i in k} gamma_i det_i      (1)

with T the window t0 is drawn over. Replace det_i by any selection -- Roman detected it, the
mass is measured to 10% -- to get that subset's yield. N is linear in F, so it is computed for
F = 1 and scaled.

TWO CONVENTIONS, both reported. The simulator keeps a drawn star only with probability equal
to its blend fraction (the legacy method: Sajadian & Makler, criterion ii), which counts events
per RESOLVED OBJECT. Weighting each draw by 1/P(kept) instead counts events on every star,
faint blended sources included. (1) as written is the first; the Sajadian papers' yields are
in that convention, so comparisons with them use it.

CHECKS, printed and written before any yield:
  tau   (pi / 2u0m) * <gamma * tE> over draws must reproduce the optical depth the C++ computes
        independently (optical_depth(), column opt_1e6), sightline by sightline. This tests the
        constants, the units and Z: a factor error anywhere shows here.
  OGLE  for the `bulge` population only (F = 1 by definition: every lens is a star or remnant),
        the rate per source star with I < 21, u0 < 1, against Mroz et al. (2019).
  legacy  the map file's own Neven, summed over the scan, which has known defects (it averages
        eps/tE over DETECTED events, has no sqrt(M) v weight and no F) -- reported, not trusted.

    .roman/bin/python analysis/y1_absolute_yield.py \\
        --run bh=runs/prod_bh_20260917 --run ns=runs/prod_ns_20260917 \\
        --run bulge=../roman_runs/2026-09-06_v3_h7 -o figures/yield_20260922
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R          # noqa: E402

COLS = ["tE", "Ml", "Vt", "Ds", "u0", "opt_1e6", "struc", "lon", "lat", "w_area",
        "detL", "detR", "detJ", "relMl_J", "relMl_R", "ndw_R"] + \
       [f"magb_{f}" for f in R.FILTERS] + [f"blend_{f}" for f in R.FILTERS]

# Fractions of the Galaxy's stellar mass in the population. Literature, checked 2026-09-22:
#   BH 0.004-0.005  Sajadian & Makler (arXiv:2608.16448), their F
#   BH ~0.01        Sweeney et al. 2022: NS+BH together ~1% (with natal kicks), so BH below it
#   BH ~0.03        Gould 2000 (bulge census 69:22:6:3 MS:WD:NS:BH by mass); Olejak et al. 2020
#                   (1.2e8 BHs x 14 Msun over ~6e10 Msun); Lam et al. 2020 (2e8 BHs, 5-16 Msun)
#   NS ~0.06        Gould 2000, no natal kicks; Sweeney et al. 2022 find 40% of NSs escape
DEFAULT_F = {"bh": [0.005, 0.01, 0.03], "ns": [0.005, 0.01, 0.03, 0.06], "bulge": [1.0]}

YEAR_S = 365.25 * 86400.0

# Mroz et al. 2019 (ApJS 244, 29), OGLE-IV, |l| < 3 deg, southern fields, sources I < 21,
# u0 < 1: Gamma = G0 exp(cG (3 - |b|)), tau300 = T0 exp(cT (3 - |b|)).
MROZ_SOUTH = dict(G0=14.3e-6, cG=0.52, T0=1.39e-6, cT=0.41)
# i_AB - I_Vega for a typical reddened bulge giant; the CMD magnitudes are AB. The cut is
# I_Vega < 21, i.e. i_AB < 21 + this. Approximate; the rate per star depends on it weakly.
I_AB_MINUS_VEGA = 0.4


class Run:
    def __init__(self, name, directory, chunksize):
        self.name, self.dir = name, directory
        prov_path = os.path.join(directory, "files/MONTLMC/files/run_provenance.txt")
        self.prov = R.load_provenance(prov_path) if os.path.exists(prov_path) else {}
        self.population = R.population(self.prov)
        self.pop_key = self.population.split()[0]
        tag = self.prov.get("population_tag", "5")
        self.table = os.path.join(directory, f"test{tag}.dat")
        self.map = os.path.join(directory, f"files/MONTLMC/files/MapLMC{tag}.dat")
        self.logs = [os.path.join(directory, f) for f in ("run.log", "run2.log")
                     if os.path.exists(os.path.join(directory, f))]
        self.tobs_days = float(self.prov.get("Tobs_days", 3652.43))
        print(f"[{name}] reading {self.table} ({self.population})", flush=True)
        self.df = R.load_events(self.table, usecols=COLS, chunksize=chunksize,
                                keep=R.keep_weightable(self.map, self.logs))
        self.sl = R.load_sightlines(self.map)
        override = R.nsim_from_logs(self.logs)
        df = self.df
        df["y"] = R.yield_weight(df, self.sl, self.tobs_days, nsim_override=override)
        df["gamma"] = R.draw_rate(df)
        df["P"] = R.acceptance_probability(df)
        df["key"] = list(zip(df["lon"].round(3), df["lat"].round(3)))
        roman_keys = set(df.loc[df["ndw_R"] > 0, "key"])
        df["foot"] = df["key"].isin(roman_keys)
        self.nsim = {**{(round(a, 3), round(b, 3)): n for a, b, n in
                        zip(self.sl.get("lon", []), self.sl.get("lat", []),
                            self.sl.get("nsim", [])) if np.isfinite(a)}, **override}
        print(f"[{name}] {len(df):,} weightable draws, {int(df.detJ.sum()):,} joint "
              f"detections, {df.key.nunique()} sightlines", flush=True)


# ---------------------------------------------------------------------------------------------
def raw_numbers(run):
    df = run.df
    per = df.groupby("key").agg(rows=("detJ", "size"), w=("w_area", "first"))
    nsim = np.array([run.nsim.get(k, np.nan) for k in per.index])
    return {
        "sightlines": len(per),
        "area_deg2": per["w"].sum(),
        "area_footprint_deg2": per.loc[[k in set(df.key[df.foot]) for k in per.index], "w"].sum(),
        "draws (table rows)": len(df),
        "draws (sum of nsim)": np.nansum(nsim),
        "rows == nsim on every sightline": bool(np.allclose(per["rows"].to_numpy(), nsim)),
        "joint detections (MC count)": int(df.detJ.sum()),
        "Rubin detections (MC count)": int(df.detL.sum()),
        "Roman detections (MC count)": int(df.detR.sum()),
        "mean lens mass <M> [Msun], per component": df.groupby("struc")["Ml"].mean().round(4).to_dict(),
    }


def tau_check(run):
    """(pi / 2u0m) <gamma tE> per sightline vs the C++'s own optical depth, same draws."""
    df = run.df
    est = np.pi / (2.0 * R.U0M) * df["gamma"] * df["tE"] * 86400.0     # tE days -> s
    g = pd.DataFrame({"key": df["key"], "est": est, "cpp": df["opt_1e6"] * 1e-6})
    per = g.groupby("key").agg(est=("est", "mean"), cpp=("cpp", "mean"), n=("est", "size"))
    per = per[per.n >= 200]
    r = per["est"] / per["cpp"]
    return {"sightlines (>=200 draws)": len(per),
            "ratio median": r.median(), "ratio 16-84%": (r.quantile(.16), r.quantile(.84)),
            "ratio pooled (sum est / sum cpp)": (per.est * per.n).sum() / (per.cpp * per.n).sum(),
            "tau C++ median": per.cpp.median()}


def ogle_check(run):
    """Rate per source star, I < 21, u0 < 1, |l| < 3, b < 0, against Mroz et al. 2019."""
    df = run.df
    isrc = df["magb_i"] - 2.5 * np.log10(df["blend_i"].clip(lower=1e-12))   # source's own i
    sel = (isrc < 21.0 + I_AB_MINUS_VEGA) & (df["lon"].abs() < 3.0) & (df["lat"] < 0)
    g = pd.DataFrame({"b": df["lat"][sel].round(1),
                      "gam": df["gamma"][sel] / R.U0M * YEAR_S,
                      "tau": df["opt_1e6"][sel] * 1e-6})
    per = g.groupby("b").agg(gam=("gam", "mean"), tau=("tau", "mean"), n=("gam", "size"))
    per = per[per.n >= 500]
    m = MROZ_SOUTH
    per["gam_ogle"] = m["G0"] * np.exp(m["cG"] * (3.0 - per.index.to_numpy().__abs__()))
    per["tau_ogle"] = m["T0"] * np.exp(m["cT"] * (3.0 - per.index.to_numpy().__abs__()))
    per["gam_ratio"] = per.gam / per.gam_ogle
    per["tau_ratio"] = per.tau / per.tau_ogle
    return per


def legacy_neven(run):
    sl = run.sl
    return float(np.nansum(10.0 ** sl["log10_Neven"] * sl["w_area"]))


# ---------------------------------------------------------------------------------------------
SELECTIONS = [
    ("joint (either survey)",     lambda d: d.detJ == 1),
    ("Rubin detects",             lambda d: d.detL == 1),
    ("Roman detects",             lambda d: d.detR == 1),
    ("both detect",               lambda d: (d.detL == 1) & (d.detR == 1)),
    ("Rubin only",                lambda d: (d.detL == 1) & (d.detR == 0)),
    ("Roman only",                lambda d: (d.detR == 1) & (d.detL == 0)),
    ("Roman, sigma(M)/M < 1%",    lambda d: (d.detR == 1) & (d.relMl_R > 0) & (d.relMl_R < .01)),
    ("Roman, sigma(M)/M < 5%",    lambda d: (d.detR == 1) & (d.relMl_R > 0) & (d.relMl_R < .05)),
    ("Roman, sigma(M)/M < 10%",   lambda d: (d.detR == 1) & (d.relMl_R > 0) & (d.relMl_R < .10)),
    ("joint, sigma(M)/M < 10%",   lambda d: (d.detJ == 1) & (d.relMl_J > 0) & (d.relMl_J < .10)),
]


def yields(df, scale=1.0):
    """Rows of (selection, N per object, err, N all stars, err, N_mc) at F = 1."""
    out = []
    y_obj = df["y"].to_numpy() * scale
    P = df["P"].to_numpy()
    for name, f in SELECTIONS:
        s = f(df).to_numpy()
        bad = int((s & (P <= 0)).sum())
        y_all = np.where(P > 0, y_obj / np.where(P > 0, P, 1.0), 0.0)
        out.append(dict(selection=name, n_mc=int(s.sum()), p_zero=bad,
                        N_obj=y_obj[s].sum(), err_obj=np.sqrt((y_obj[s] ** 2).sum()),
                        N_all=y_all[s].sum(), err_all=np.sqrt((y_all[s] ** 2).sum())))
    return pd.DataFrame(out)


def mass_subset(df, lo, hi, full_lo, full_hi):
    """Yield for a log-flat mass function over [lo, hi] from draws of one over [full_lo, full_hi].

    Draws with M in [lo, hi] are log-flat within it, so they are a sample of the narrower
    population; nsim counted every draw, so the kept ones are up-weighted by 1/P(M in [lo,hi]),
    and <M> is that of the narrower range (the rate is per lens: F rho / <M>).
    """
    keep = (df["Ml"] >= lo) & (df["Ml"] <= hi)
    frac = np.log(hi / lo) / np.log(full_hi / full_lo)
    mean_full = (full_hi - full_lo) / np.log(full_hi / full_lo)
    mean_sub = (hi - lo) / np.log(hi / lo)
    return yields(df[keep], scale=(mean_full / mean_sub) / frac)


def fmt(n, e):
    if n == 0:
        return "0"
    return f"{n:.3g} ± {e:.2g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", action="append", required=True, help="name=directory")
    ap.add_argument("--F", action="append", default=[],
                    help="population=f1,f2,... (default: the literature grid in DEFAULT_F)")
    ap.add_argument("-o", "--out", required=True, help="output directory")
    ap.add_argument("--chunksize", type=int, default=500_000)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    Fgrid = dict(DEFAULT_F)
    for spec in a.F:
        k, v = spec.split("=")
        Fgrid[k] = [float(x) for x in v.split(",")]

    lines = ["# Absolute yields (Step Y, Deviation 54)", "",
             "Generated by `analysis/y1_absolute_yield.py`. N = expected events whose peak falls "
             "in the 10-yr window, over the scanned sky; ± is the Monte Carlo error only. "
             "'per object' = the legacy blend-fraction convention (events per resolved object, "
             "as in Sajadian & Sahu 2023 / Sajadian & Makler); 'all stars' = every star counted.",
             ""]
    csv_rows = []
    for spec in a.run:
        name, directory = spec.split("=", 1)
        run = Run(name, directory, a.chunksize)
        commit = run.prov.get("git_commit", "?")
        lines += [f"## `{name}` -- {run.population}, commit `{commit}`", ""]

        raw = raw_numbers(run)
        lines += ["### Raw Monte Carlo numbers (sample sizes, NOT yields)", ""]
        lines += [f"- {k}: {v:,.4g}" if isinstance(v, float) else f"- {k}: {v}"
                  for k, v in raw.items()] + [""]

        tc = tau_check(run)
        lines += ["### Check 1: optical depth from the rate constants vs the C++", ""]
        lines += [f"- {k}: " + (f"{v[0]:.4f}-{v[1]:.4f}" if isinstance(v, tuple) else f"{v:.4g}")
                  for k, v in tc.items()] + [""]
        print(f"[{name}] tau check: {tc}", flush=True)

        if run.pop_key == "bulge":
            og = ogle_check(run)
            lines += ["### Check 2: rate per star (I < 21, u0 < 1) vs OGLE-IV (Mroz et al. 2019)",
                      "", "| b | draws | Gamma model [1e-6/yr] | Gamma OGLE | ratio | tau model [1e-6] "
                      "| tau OGLE | ratio |", "|---|---|---|---|---|---|---|---|"]
            for b, r in og.iterrows():
                lines.append(f"| {b:.1f} | {int(r.n):,} | {r.gam*1e6:.2f} | {r.gam_ogle*1e6:.2f} "
                             f"| {r.gam_ratio:.2f} | {r.tau*1e6:.2f} | {r.tau_ogle*1e6:.2f} "
                             f"| {r.tau_ratio:.2f} |")
            lines.append("")
            print(og.to_string(), flush=True)

        leg = legacy_neven(run)
        lines += ["### Check 3: the legacy map-file Neven, summed (F = 1, per object)", "",
                  f"- sum Neven x w_area = {leg:,.4g} events", ""]

        for scope, sub in (("whole scan", run.df), ("Roman footprint only", run.df[run.df.foot])):
            t = yields(sub)
            if t.p_zero.sum():
                print(f"[{name}] WARNING: {t.p_zero.sum()} selected draws with rebuilt "
                      f"P(kept) = 0 -- the acceptance rebuild disagrees with the C++", flush=True)
            Fs = Fgrid.get(run.pop_key, [1.0])
            lines += [f"### Yields, {scope}", "",
                      "| selection | MC draws | " + " | ".join(
                          f"F={F:g} per object | F={F:g} all stars" for F in Fs) + " |",
                      "|---|---|" + "---|---|" * len(Fs)]
            for _, r in t.iterrows():
                cells = []
                for F in Fs:
                    cells += [fmt(F * r.N_obj, F * r.err_obj), fmt(F * r.N_all, F * r.err_all)]
                    csv_rows.append(dict(run=name, population=run.population, commit=commit,
                                         scope=scope, selection=r.selection, F=F,
                                         n_mc=r.n_mc, N_per_object=F * r.N_obj,
                                         err_per_object=F * r.err_obj,
                                         N_all_stars=F * r.N_all, err_all_stars=F * r.err_all))
                lines.append(f"| {r.selection} | {r.n_mc:,} | " + " | ".join(cells) + " |")
            lines.append("")

            if run.pop_key == "bh":
                t2 = mass_subset(sub, 3.0, 50.0, 3.0, 1000.0)
                lines += [f"#### Comparison row: the same draws re-weighted to log-flat 3-50 Msun "
                          f"({scope}) -- NOT the population used, for comparing with "
                          f"Sajadian & Sahu 2023", "",
                          "| selection | MC draws | " + " | ".join(
                              f"F={F:g} per object | F={F:g} all stars" for F in Fs) + " |",
                          "|---|---|" + "---|---|" * len(Fs)]
                for _, r in t2.iterrows():
                    cells = []
                    for F in Fs:
                        cells += [fmt(F * r.N_obj, F * r.err_obj), fmt(F * r.N_all, F * r.err_all)]
                    lines.append(f"| {r.selection} | {r.n_mc:,} | " + " | ".join(cells) + " |")
                lines.append("")
        del run

    with open(os.path.join(a.out, "y1_yields.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    pd.DataFrame(csv_rows).to_csv(os.path.join(a.out, "y1_yields.csv"), index=False)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
