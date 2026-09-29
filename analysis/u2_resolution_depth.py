#!/usr/bin/env python3
"""Step U2: how much Roman's image-resolution probability depends on its assumed depth.

WHY THIS EXISTS. Step R1's per-epoch count (Deviation 48) calls both lensed images detectable
when each is between saturation and the filter's single-visit depth, `thre` in Bulge.h. For F146
that is 29.0 mag, commented "(value needs to change)": a placeholder. Roman's own photometric
error table, files/sigma_roman.txt -- the noise model every Roman epoch in the simulation uses --
reaches 5 sigma at F146 = 25.5 and ends at 27.0. So the resolution count accepts minor images down
to 29 mag, where a single exposure has SNR ~ 1, and the faint minor image is exactly what decides
resolvability: the images separate only as it fades (A_- -> 1/u^4). The counts cannot be redone
from the table -- the per-epoch quantities never leave the light-curve loop -- so this script
rebuilds them semi-analytically, per detected event, on Roman's real epoch list:

    u(t) = sqrt(u0^2 + ((t - t0)/tE)^2)                       (no parallax)
    A_+- = A/2 +- 1/2,   m_+- = m_src - 2.5 log10 A_+-,   m_src = m_base - 2.5 log10 f_b
    m_tot = m_base - 2.5 log10(f_b A + 1 - f_b)               (the recorded, blended magnitude)
    epoch counts if  12 <= m_tot <= 29  (the simulator's recording gate),
                     12 <= m_+- <= depth,
                     theta_E sqrt(u^2+4) >= bar,  bar = D sigma_a(m_tot)  or  the PSF FWHM
    resolvable if >= 3 epochs count  (Sajadian & Makler's N >= 3)

with sigma_a = errRomanA() mirrored from helper.cpp. VALIDATION FIRST: at depth = 29.0 the model
must reproduce the simulator's own nres{5,20,PSF}_R >= 3, event by event and in the weighted
fraction; only then is the depth-25.5 number a statement about the depth rather than about the
model's approximations (it has no parallax and no per-field timing).

    .roman/bin/python analysis/u2_resolution_depth.py \\
        --run bulge=runs/prod_bulge_20260924 --run bh=runs/prod_bh_20260924 \\
        --run ns=runs/prod_ns_20260924 --mean-mass figures/yield_20260925 -o figures/u1_20260929
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import romanlib as R                 # noqa: E402
import u1_report_numbers as U        # noqa: E402

# Mirrored from Bulge.h.
SATU_F146 = 12.0
THRE_F146_CODE = 29.0                # the placeholder the production runs used
PSF_FWHM_MAS = 105.0
ROMAN_AST_FLOOR = 1.1
ROMAN_AST_MFLR, ROMAN_AST_MBKG, ROMAN_AST_SBKG = 20.62, 23.5, 10.0
ROMAN_AST_SLOPE_SRC, ROMAN_AST_SLOPE_BKG = 0.33285, 0.4
FOV_ROMAN_DEG = 0.3003
SIGMA_ROMAN = "files/sigma_roman.txt"


def err_roman_a(m):
    """helper.cpp errRomanA(), vectorised."""
    e = np.where(m > ROMAN_AST_MBKG,
                 ROMAN_AST_SBKG * 10.0 ** (ROMAN_AST_SLOPE_BKG * (m - ROMAN_AST_MBKG)),
                 np.where(m > ROMAN_AST_MFLR,
                          ROMAN_AST_FLOOR * 10.0 ** (ROMAN_AST_SLOPE_SRC * (m - ROMAN_AST_MFLR)),
                          ROMAN_AST_FLOOR))
    return np.maximum(e, ROMAN_AST_FLOOR)


def depth_at_snr(snr):
    """F146 magnitude at which sigma_roman.txt's error reaches 1.0857/snr."""
    t = np.loadtxt(SIGMA_ROMAN)
    return float(t[np.argmax(t[:, 1] >= 1.0857 / snr), 0])


def roman_times(l0=0.4, b0=-1.2):
    """Roman's unique epoch times at one footprint sightline (every footprint sightline has the
    same 50,401: one mission-wide visit list)."""
    rom = np.loadtxt("Baseline/RomanBaseline.dat", comments="#", usecols=(3, 4, 5))
    m = np.hypot(rom[:, 0] - l0, rom[:, 1] - b0) <= FOV_ROMAN_DEG
    return np.unique(rom[m, 2])


def count_epochs(ev, times, depths):
    """{(depth, bar): per-event count of resolvable epochs}, bars 5, 20, PSF."""
    out = {(d, b): np.zeros(len(ev), np.int32) for d in depths for b in ("5", "20", "PSF")}
    te, u0, tE, t0 = (ev[c].to_numpy(float) for c in ("tetE", "u0", "tE", "t0"))
    mb = ev["magb_F146"].to_numpy(float)
    fb = np.clip(ev["blend_F146"].to_numpy(float), 1e-12, 1.0)
    msrc = mb - 2.5 * np.log10(fb)
    for i in range(len(ev)):
        tau = (times - t0[i]) / tE[i]
        u = np.sqrt(u0[i] ** 2 + tau * tau)
        root = np.sqrt(u * u + 4.0)
        A = (u * u + 2.0) / (u * root)
        am = A / 2.0 - 0.5
        ok = am > 0
        mtot = mb[i] - 2.5 * np.log10(fb[i] * A + 1.0 - fb[i])
        rec = ok & (mtot >= SATU_F146) & (mtot <= THRE_F146_CODE)
        mp = msrc[i] - 2.5 * np.log10(A / 2.0 + 0.5)
        mm = np.where(ok, msrc[i] - 2.5 * np.log10(np.where(ok, am, 1.0)), 99.0)
        sep = te[i] * root
        sa = err_roman_a(mtot)
        for d in depths:
            both = rec & (mp >= SATU_F146) & (mp <= d) & (mm >= SATU_F146) & (mm <= d)
            out[(d, "5")][i] = int((both & (sep >= 5.0 * sa)).sum())
            out[(d, "20")][i] = int((both & (sep >= 20.0 * sa)).sum())
            out[(d, "PSF")][i] = int((both & (sep >= PSF_FWHM_MAS)).sum())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", action="append", required=True, metavar="NAME=DIR")
    ap.add_argument("--mean-mass", required=True)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()

    d5, d3 = depth_at_snr(5.0), depth_at_snr(3.0)
    depths = (THRE_F146_CODE, d3, d5)
    print(f"depths: code {THRE_F146_CODE}, SNR3 {d3:.2f}, SNR5 {d5:.2f} (from {SIGMA_ROMAN})")
    times = roman_times()
    print(f"Roman epochs: {len(times):,}")
    if "t0" not in U.COLS:
        U.COLS.append("t0")

    rows = []
    for spec in a.run:
        name, directory = spec.split("=", 1)
        run = U.Run(name, directory, a.mean_mass)
        d = run.df
        det = (d.detR == 1).to_numpy()
        ev = d[det].reset_index(drop=True)
        W, y = ev["W"].to_numpy(), ev["y"].to_numpy()
        b = U.Boot(len(ev), seed=U.SEED + 10)
        allm = np.ones(len(ev), bool)
        print(f"[{name}] counting epochs for {len(ev):,} Roman detections", flush=True)
        cnt = count_epochs(ev, times, depths)
        # --- validation against the simulator's own counts, at the code's depth ---
        for bar, col in (("5", "nres5_R"), ("20", "nres20_R"), ("PSF", "nresPSF_R")):
            sim = ev[col].to_numpy() >= U.RESOLVE_MIN_EPOCHS
            mod = cnt[(THRE_F146_CODE, bar)] >= U.RESOLVE_MIN_EPOCHS
            agree = float(np.mean(sim == mod))
            fs, fm = b.frac(sim, allm, W), b.frac(mod, allm, W)
            print(f"[{name}] D={bar:3s} depth 29: simulator {100*fs[0]:.3f}%  model "
                  f"{100*fm[0]:.3f}%  event agreement {100*agree:.2f}%", flush=True)
            rows.append(dict(population=name, bar=bar, depth="29.0 (simulator)",
                             P=100 * fs[0], err=100 * fs[1], N1=y[sim].sum(),
                             N1_err=np.sqrt((y[sim] ** 2).sum()), n=int(sim.sum()),
                             agreement=np.nan))
            for dep in depths:
                mod = cnt[(dep, bar)] >= U.RESOLVE_MIN_EPOCHS
                f = b.frac(mod, allm, W)
                rows.append(dict(population=name, bar=bar, depth=f"{dep:.2f} (model)",
                                 P=100 * f[0], err=100 * f[1], N1=y[mod].sum(),
                                 N1_err=np.sqrt((y[mod] ** 2).sum()), n=int(mod.sum()),
                                 agreement=100 * agree if dep == THRE_F146_CODE else np.nan))
        # --- how faint are Roman's detected sources? (the same placeholder, for detection) ---
        mb = ev["magb_F146"].to_numpy(float)
        src = mb - 2.5 * np.log10(np.clip(ev["blend_F146"].to_numpy(float), 1e-12, 1.0))
        for lab, v in (("baseline", mb), ("source", src)):
            for cut in (d5, 27.0):
                f = b.frac(v > cut, allm, W)
                rows.append(dict(population=name, bar=f"{lab} F146 > {cut:.2f}",
                                 depth="Roman detections", P=100 * f[0], err=100 * f[1],
                                 N1=y[v > cut].sum(), n=int((v > cut).sum())))
        f4 = (d["ndw_R"] > 0).to_numpy()
        q = run.boot.quantile(d["magb_F146"].to_numpy(float), f4, d["W"].to_numpy())
        rows.append(dict(population=name, bar="median baseline F146, all footprint detections",
                         depth="", P=q[0], err=q[1], n=q[2]))
        q = run.boot.quantile(err_roman_a(d["magb_F146"].to_numpy(float)), f4,
                              d["W"].to_numpy())
        rows.append(dict(population=name, bar="median sigma_a(baseline) [mas], footprint dets",
                         depth="", P=q[0], err=q[1], n=q[2]))
        del run

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(a.out, "u2_resolution_depth.csv"), index=False)
    with pd.option_context("display.width", 200, "display.max_rows", 500):
        print(out.to_string(index=False))


if __name__ == "__main__":
    main()
