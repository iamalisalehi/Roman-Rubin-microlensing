#!/usr/bin/env python3
"""M9: production run-time estimate from a timed sample of sightlines (2026-10-02).

Since Deviation 77 every sightline has its own random stream, so the sightlines timed by
runs/m9_timing.sh (and runs/m9_timing_extra.sh) are exactly the ones a production run would simulate.
The scan (analysis/gbtds_geometry, --stride 10 --stride-roman 5) splits into three strata whose costs
differ by orders of magnitude:
    C  footprint sightlines ON a Roman detector (either roll)   -- ~50,000 Roman epochs each
    U  footprint-stratum sightlines on no detector (chip gaps, outline margin) -- Rubin only
    O  outside the footprint (coarse grid)                      -- Rubin only
Per stratum: mean CPU per sightline (user time minus the run's start-up, measured on the sample),
bootstrap error; total = sum N_stratum x mean + one start-up. Wall time for k concurrent processes on
this laptop (2 physical cores + hyper-threading) uses the measured parallel wall/CPU ratio.

    .roman/bin/python analysis/m9_time_estimate.py
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gbtds_geometry as G                         # noqa: E402

ROOT = G.ROOT
STARTUP_S = 21.0          # under 3-process load: the smallest total of any timed sightline (bh, outside)


def main():
    g = G.scan_sightlines(10, 5)
    od = G.on_detector(g["lon"], g["lat"])
    cls = np.where(g["fine"], np.where(od.any(axis=1), "C", "U"), "O")
    n = {s: int((cls == s).sum()) for s in "CUO"}
    rows = []
    for f in glob.glob(os.path.join(ROOT, "runs", "m9_timing*_*", "times.csv")):
        t = pd.read_csv(f, comment=None)
        t = t[t["pop"].isin(["bulge", "bh", "ns"])]
        rows.append(t)
    t = pd.concat(rows, ignore_index=True)
    t["index"] = t["index"].astype(int)
    t["stratum"] = [cls[i] for i in t["index"]]
    t["cpu"] = (t.user_s.astype(float) + t.sys_s.astype(float) - STARTUP_S).clip(lower=0.0)
    t["ratio"] = t.wall_s.astype(float) / (t.user_s.astype(float) + t.sys_s.astype(float))
    rng = np.random.default_rng(1)
    print(f"scan: {len(cls)} sightlines -- C {n['C']} (on a Roman detector), U {n['U']}, O {n['O']}")
    print(f"wall/CPU while 3 processes share the laptop: median {t.ratio.median():.2f}")
    out = []
    for pop, gp in t.groupby("pop"):
        tot, var = STARTUP_S, 0.0
        parts = []
        for s in "CUO":
            x = gp[gp.stratum == s].cpu.to_numpy()
            if len(x) == 0:
                parts.append(f"{s}: no sample"); continue
            m = x.mean()
            boot = [rng.choice(x, len(x)).mean() for _ in range(2000)]
            e = np.std(boot)
            tot += n[s] * m
            var += (n[s] * e) ** 2
            parts.append(f"{s}: {len(x)} timed, mean {m/60:.1f} +- {e/60:.1f} min CPU")
        h, eh = tot / 3600, np.sqrt(var) / 3600
        out.append(dict(pop=pop, cpu_h=h, err_h=eh))
        print(f"{pop:5s}: total CPU {h:.1f} +- {eh:.1f} h   [{'; '.join(parts)}]")
    o = pd.DataFrame(out)
    r = t.ratio.median()
    print(f"\nAll three in parallel on this laptop: wall ~ max(CPU) x {r:.2f} = "
          f"{o.cpu_h.max() * r:.1f} h (+- {o.err_h.max() * r:.1f}); summed CPU {o.cpu_h.sum():.1f} h")
    o.to_csv(os.path.join(ROOT, "runs", "m9_time_estimate.csv"), index=False)


if __name__ == "__main__":
    main()
