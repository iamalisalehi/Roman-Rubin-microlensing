# Notes for the next overview report (after the 2026-10 production runs)

**What this is.** A running notebook, filled in AS the pre-production fixes and the new runs happen,
so that writing `Report/overview_v2/overview_v2_report.tex` later is assembly, not archaeology. The
new report must cover everything `Report/overview/overview_report.tex` covers (user, 2026-10-01).
It starts from `Report/TEMPLATE_report.tex` like every report (`Report/README.md`).

Plan these notes follow: `/home/ali/.claude/plans/i-want-to-do-resilient-hummingbird.md`
(Step 1 footprint, Step 2 dust, Step 3 astrometric floor; then the minor items; then the runs).

## Checklist for writing (from Report/README.md and the user's standing requests)

- Audience is a scientific advisor; abstract gives yields and physics, never MC sample sizes; no
  "for a newcomer" framing; flowcharts, not code history.
- Every pooled number event-rate weighted, quoted with N_eff; yields quoted with their F; sample
  sizes kept apart from yields.
- Every number from a script's output file, cited in a `%` comment beside it.
- State known-bug status of the runs used (which fixes they include; which open items remain).
- Sajadian & Sahu (2023): compare with the published paper only, collegial framing ("this
  simulation / the paper"); Dr. Sajadian is supervisor and collaborator.
- Penny et al. 2019: both 27,000 and 54,000 are detections (no "bracket").
- Build clean: `grep -c '^TODO:' *.log` = 0, no undefined refs, figures read at page size.

## Model changes since the old report (fill in as each lands)

| Step | Change | Deviation | Commit | Moves which results |
|---|---|---|---|---|
| 1 | Adopted GBTDS layout (spring/autumn, 18-SCA detectors); scan region by distance rule; Roman start day 730 -> 306 (2027-02-11); Rubin visit list 3,686 -> 12,308; footprint area post-stratified | 69 | 078756b | Roman footprint yields, Rubin whole-scan totals, gap geometry (timeline vs Rubin seasons) |
| 2 | Extinction tables rebuilt from DECaPS + Marshall (dustmaps), regular grid, one file; hardened reader | 70 | (pending) | everything, most near the plane (U5 estimated x0.4-0.9) |
| 3 | Astrometric reference position audit; correlated-floor bracket | | | every theta_E / mass number |

## Section by section

Old section -> what the new one needs. "Old source" is where the old report's content came from.

### 1. The science question
- Old: 1.1 microlensing and dark-lens mass; 1.2 complementary surveys. Physics text; mostly reusable.
- New: update the Roman description to the adopted layout (six fields, 1.7 deg^2, two roll angles).

### 2. How the simulation works
- 2.1 Inputs. Old source: `make_timeline.py` -> `timeline.tex`. **Regenerate**: Roman now starts
  2027-02-11 (sim day 306), so the timeline and the Rubin-only stretches change (~0.8 yr before Roman,
  more after). Cite the vendored layout (`Baseline/gbtds_layout/`, Penny's gbtds_optimizer, 2026.4.3).
- 2.2 One event; 2.3 three Fisher matrices (TikZ). Update 2.3 if Step 3 changes the astrometric
  parameter set (reference position) or adds the noise bracket.
- 2.4 Covering the sky. Old: `make_footprints.py` -> `figures/footprint_20260929/`. **Regenerate**
  both figures on the new geometry; the "simulated vs real" figure should now agree by construction
  -- say so and show it.
- 2.5 Events -> sky (weights, w_area). Unchanged in method; restate the new scan area.

### 3. The analyses, and what each tells us
- Table of analyses; add any new ones (noise-bracket columns, u6 on the new tables).

### 4. Results (all regenerate on the new runs)
Old figure -> script: f2_gap_filling (f2), f3_characterization_map (f3), f4_fisher_precision (f4),
h5_astrometric_shift / h5_astrometry_summary (h5), h3a_where / h3c_honest (h3), p6_synergy,
p7_precision, p6_resolution (u1/u2), y3_yield_vs_F (y1->y3), samples (S3 + s2 plotter,
make_sample_table.py). Pooled numbers: `analysis/u1_report_numbers.py`.
- 4.1 sample; 4.2 who detects what; 4.3 gap filling (re-measure: U5 could only bound it);
  4.4 Roman for Rubin; 4.5 combination; 4.6 precision; 4.7 astrometric shift and theta_E (now
  with the noise bracket); 4.8 satellite parallax; 4.9 two baselines; 4.10 bh/ns; 4.11 image
  resolution (Rubin's fractions were never dust-corrected); 4.12 absolute yields; 4.13 samples
  (new S3-style sample runs needed on the new geometry).
- Old-vs-new to show: U5's dust-corrected predictions (`figures/u1_20260929/u5_corrected_numbers.md`)
  against the new runs (they should agree within U5's errors; if not, explain); U3's footprint
  estimate (+1..+5% Roman) against the actual change.

### 5. How certain are these results?
- 5.1 MC uncertainty (u1 bootstrap) -- regenerate.
- 5.2 Modelling uncertainties -- re-rank: footprint and dust items should leave the list or shrink to
  residuals (e.g. LSSTCam outline vs circle, CCM89 for F146); the astrometric floor becomes the
  bracket width.
- 5.3 "Simulated footprint against the real one" -> becomes a short validation (overlay agrees).
- 5.4 "Dust near the plane" -> becomes a validation: new tables vs VVV (`u6_vvv_check.py`).
- 5.5 Published comparisons: Penny et al. (dust-corrected estimate was 49.4 +- 2.9 per deg^2 per
  day, (23 +- 5)% below their 63.8 -- re-measure); SS23 (rate-weighted vs unweighted, u7).
- 5.6 Summary: what to believe.

### 6. Summary
### 7. Next steps (whole-footprint cluster run; open items in plain language)
### Appendices: all sample events; gallery

## Log of what each step produced (append as you go)

<!-- Step N (date, Deviation, commit): outputs, headline numbers, figures, commands. -->

### Step 1 -- footprint (2026-10-01, Deviation 69)
- **What the new report should say (Sec. 2.1 / 2.4).** Roman: the adopted GBTDS layout (Penny's
  gbtds_optimizer 2026.4.3, centres confirmed by STScI: spring l 0.5, autumn l 0.35, b -1.4; GC at
  b -0.22), each field 18 detectors of 0.125 deg; the telescope rolls 180 deg between spring and
  autumn, so a star in a chip gap in one roll is usually seen in the other. Exact areas: 1.684 deg^2
  per roll, 1.878 seen in at least one roll, 1.487 in both; 0.982 deg^2 of the footprint stratum is
  in gaps in both. Mission 2027-02-11 -> 2031-10-25 (sim days 306-2024): 6 high-cadence + 4
  low-cadence seasons on the published windows; Rubin-only baseline ~0.8 yr before, ~4.5 yr after.
- Rubin: 12,308 visits from baseline_v5.1.0 reach the scan (was 3,686: the old box undercounted
  sightlines near its edge). Scan region: within 3.94 deg of any Roman field centre, 69.3 deg^2,
  2,013 sightlines at the production flags (286 footprint at 0.1 deg).
- **Method point worth a sentence:** point-sampling the footprint aliases against the chip gaps
  (0.1-deg grid: 13% low area); the sky-area weights are post-stratified by coverage class, so
  Roman's sky area is exact at any grid step. Flowchart "Covering the sky" needs a box for this.
- Figures: `figures/footprint_20261001/footprint_simulated` (scan, Rubin visits, both rolls'
  detectors) and `footprint_vs_gbtds` (overlay on the Aladin tiles: spring agrees to ~0.02 deg;
  the screenshot's autumn tiles are at the spring centres, 0.15 deg from STScI's -- say so in the
  caption or drop the autumn overlay). Commands:
  `.roman/bin/python Report/overview/make_footprints.py [--log runs/<new>/run.log]`,
  `.roman/bin/python Report/overview/make_timeline.py` -> `Report/overview_v2/timeline.tex`.
- Supersedes: old Sec. 5.3 "simulated footprint against the real one" (Deviation 59/60, U3) --
  becomes a short validation paragraph; Table 10 (U3 estimate +1..+5% Roman) can be compared
  with the actual change once the runs exist.

### Step 2 -- extinction (2026-10-01, Deviation 70)
- **What the new report should say (Sec. 2.1 inputs, and Sec. 5's dust subsection becomes a
  validation):** A_V(d) from DECaPS (Zucker+2025) where it is reliable and unsaturated, and from
  Marshall+2006's near-infrared 3D map, calibrated onto DECaPS's scale (k = A_Ks/A_V = 0.0830,
  median over the five-field block at 8 kpc), where DECaPS cannot see -- within ~1 deg of the plane,
  from ~4.5 kpc. Both through dustmaps 1.0.14. On a 0.05/0.1-deg grid (15,965 positions).
- **Validation numbers to quote** (`figures/ext_20261001/v1_ext_vvv.md`): A_V(8 kpc) / VVV per |b|
  bin 0.99 / 0.98 / 0.99 / 0.91 / 0.88 (|b| 0-0.5 / 0.5-1 / 1-1.5 / 1.5-2.5 / 2.5-6); GC/five-field
  contrast 4.74 vs VVV 5.52 -> the GC field ~10% thin (known Marshall shortfall). Old tables were
  0.17 near the plane: one sentence of history, no more (advisor audience).
- **A figure worth making:** A_V(8 kpc) map of the new tables over the scan, with Roman's
  detectors, and/or tables vs VVV scatter by |b| (data: `v1_ext_vvv.csv`).
- **Old-vs-new to show once the runs exist:** the new runs against U5's dust-corrected predictions
  (`figures/u1_20260929/u5_corrected_numbers.md`), which used k 0.0805 on the old fields -- expect
  agreement within U5's errors, and explain any difference (new footprint too).
- Supersedes old Sec. 5.5 ("The dust near the Galactic plane: what is wrong and how to correct it")
  and the dust row of Sec. 5.2's uncertainty list; the GC residual becomes a small row.

