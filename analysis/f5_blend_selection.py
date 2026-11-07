#!/usr/bin/env python3
"""Blending, follow-up labels and the dark-lens test, per run (design sections 3.4 and 3.5).

    .roman/bin/python analysis/f5_blend_selection.py EVENTS --map MAP [--log LOG ...] --ms MS [-o out.md]

EVENTS is the run's event table (test5.dat or testbh.dat), MAP its sightline table (MapLMC5.dat; the
run_provenance.txt beside it names the run's selection), MS the main-sequence table
(CMD/components/ms_hypothesis.dat). Four sets of events are tabulated: detected by Rubin (detL), by
Roman (detR), by the joint fit (detJ), and characterised (joint, Abrams et al. criterion). For each:

  - the share with f_S >= 0.6 in Rubin and in Roman;
  - the shares of the follow-up classes (romanlib.blend_class);
  - for the characterised events only (a mass error exists only for them; n/a elsewhere): the share
    with P(dark) > 0.75 in each partition (joint, Rubin alone, Roman alone), and the share for which
    no verdict exists (no mass error, or no fitted blend fraction);
  - when the table has lensLum, again for characterised events: the joint P(dark) > 0.75 share among
    luminous lenses (the false-dark rate) and among dark lenses, each over the lenses of that kind
    that have a verdict.

Every number is given in both counting conventions, every-star first (the headline), with the
Kish N_eff of the weights; the dark-lens logic lives in romanlib.dark_lens_probability.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R          # noqa: E402

CONVS = ("every-star", "per-object")
PARTS = ("joint", "rubin", "roman")


def share(mask, w):
    """(weighted share of mask over w, N_eff of w) for the rows w is given on."""
    return R.weighted_fraction(mask, w), R.kish_neff(w[w > 0])


def build(df, ms, weights):
    """Rows of (label, {convention: (value, N_eff)}) for each set; df holds detected rows only."""
    sets = {"detected, Rubin (detL)": R.detected(df, "rubin"),
            "detected, Roman (detR)": R.detected(df, "roman"),
            "detected, joint (detJ)": R.detected(df, "joint"),
            "characterised (joint)": R.detected(df, "joint") & R.characterized(df, "joint")}
    char = sets["characterised (joint)"]
    pd_ = {p: R.dark_lens_probability(df[char], p, ms).reindex(df.index) for p in PARTS}
    cls = R.blend_class(df)
    out = {}
    for name, sel in sets.items():
        rows = []
        for conv in CONVS:
            w = weights[conv].to_numpy()[sel.to_numpy()]
            s = df[sel]
            c = cls[sel]
            m = lambda mask: share(mask.reset_index(drop=True), w)
            r = [("n events (unweighted)", None), ("f_S(Rubin) >= 0.6", m(s["fb0"] >= R.BLEND_CLEAN_FS)),
                 ("f_S(Roman) >= 0.6", m(s["fb1"] >= R.BLEND_CLEAN_FS))]
            r += [(f"class {k}", m(c == k)) for k in R.BLEND_CLASSES]
            is_char = name.startswith("characterised")
            for p in PARTS:
                v = pd_[p][sel]
                r.append((f"P(dark) > {R.DARK_P_MIN}, {p}", m(v > R.DARK_P_MIN) if is_char else None))
                r.append((f"no verdict, {p}", m(v.isna()) if is_char else None))
            if "lensLum" in s.columns:
                v = pd_["joint"][sel]
                for lab, k in (("luminous lenses", 1), ("dark lenses", 0)):
                    sub = (s["lensLum"] == k) & v.notna()
                    tag = "characterised " + lab + " with a verdict"
                    if not is_char:
                        r.append((f"P(dark) > {R.DARK_P_MIN}, joint, among {tag}"
                                  + (" (false-dark rate)" if k == 1 else ""), None))
                    elif sub.any():
                        r.append((f"P(dark) > {R.DARK_P_MIN}, joint, among {tag}"
                                  + (" (false-dark rate)" if k == 1 else ""),
                                  share((v[sub] > R.DARK_P_MIN).reset_index(drop=True), w[sub.to_numpy()])))
                    if k == 1 and is_char:
                        good = sub & (s["relMl_J"] < 0.3)
                        if good.any():
                            r.append((f"P(dark) > {R.DARK_P_MIN}, joint, among {tag} and relMl_J < 0.3"
                                      " (false-dark rate)",
                                      share((v[good] > R.DARK_P_MIN).reset_index(drop=True), w[good.to_numpy()])))
            rows.append((conv, r, int(sel.sum())))
        out[name] = rows
    return out


def render(res, label):
    L = [f"## {label}", ""]
    for name, rows in res.items():
        (_, ra, n), (_, rb, _) = rows
        L += [f"### {name}  (n = {n})", "", "| quantity | every-star | N_eff | per-object | N_eff |", "|---|---|---|---|---|"]
        for (lab, a), (_, b) in zip(ra[1:], rb[1:]):
            if a is None or b is None:
                L.append(f"| {lab} | n/a | n/a | n/a | n/a |")
            else:
                L.append(f"| {lab} | {a[0]:.4f} | {a[1]:,.0f} | {b[0]:.4f} | {b[1]:,.0f} |")
        L.append("")
    return L


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("events")
    ap.add_argument("--map", required=True, help="the run's sightline table (event-rate weights)")
    ap.add_argument("--log", action="append", default=[], help="run log(s), when the map file lost rows")
    ap.add_argument("--ms", required=True, help="ms_hypothesis.dat")
    ap.add_argument("--chunksize", type=int, default=500_000)
    ap.add_argument("-o", "--out", default=None, help="write the markdown here too")
    a = ap.parse_args()

    ms = R.load_ms_hypothesis(a.ms)
    det = lambda c: (c[["detJ", "detL", "detR"]] == 1).any(axis=1)
    df = R.load_events(a.events, keep=det, chunksize=a.chunksize).reset_index(drop=True)
    if df.empty:
        sys.exit(f"{a.events}: no detected events")
    weights, labels = {}, []
    for conv in CONVS:
        weights[conv], lab = R.attach_weight(df, a.map, a.log, False, conv)
        labels.append(lab)
    res = build(df, ms, weights)
    text = "\n".join([f"# Blending and dark lenses: {a.events}", "",
                      "Weights: " + "; ".join(labels) + ".", ""]
                     + render(res, os.path.basename(a.events)))
    print(text)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(text + "\n")


if __name__ == "__main__":
    main()
