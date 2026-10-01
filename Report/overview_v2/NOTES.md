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

#### Draft text: why the extinction model was changed

*(Written for the report's Sec. 2.1 / Sec. 5; reword freely, but keep the chain of evidence.)*

Every simulated source is dimmed by the dust between it and us, A_lambda = A_V(d) x (A_lambda/A_V),
with A_V(d) the extinction to distance d along its line of sight. Toward the inner bulge that
column is large, A_V ~ 5-30 mag, and it decides which sources each survey can see: through the
same dust, Rubin's r band loses 4.3 times as many magnitudes as Roman's F146. So the extinction
model sets the absolute yields of both surveys, and it sets their ratio.

The earlier runs took A_V(d) from two optical 3D dust maps, Bayestar19 and DECaPS, split by
declination as the dustmaps documentation suggests (Bayestar north of -30 deg). That put almost
all of Roman's footprint on Bayestar. Optical maps work by measuring the reddening of individual
stars at known distances; behind the dense dust lanes within ~1 deg of the plane the optical stars
are too faint to be seen, so the maps stop increasing and the extinction to bulge distances is a
lower limit. Bayestar's own reliability flag says exactly that: it was false beyond ~4 kpc on all
139 of the footprint's Bayestar sightlines, i.e. for every bulge source. DECaPS, which goes deeper,
is reliable over most of the scan but saturates (its stated limit is A_V ~ 12) on the
Galactic-centre field.

We checked this against an independent measurement that does not share the problem: the VVV
reddening map (Surot et al. 2020), built from the near-infrared colours of bulge red-clump stars,
which see through A_V ~ 30. Comparing the column to 8 kpc, the old tables had 17% of VVV's
extinction within 0.5 deg of the plane and 31% at 0.5-1 deg; and the ratio of the
Galactic-centre field's column to the five-field block's -- a test that does not depend on the
extinction law -- was 0.57 in the old tables against 4.89 in VVV. The dust was not merely
miscalibrated: it was missing where it matters most.

The new tables use DECaPS where DECaPS can see (it agrees with VVV to ~10% away from the plane),
and switch to the near-infrared 3D map of Marshall et al. (2006) where it cannot. They now track
VVV to 0.88-0.99 in every latitude bin, with no fall toward the plane. One residual remains: the
Galactic-centre field is ~10% thin (Marshall's own known shortfall there), worth ~0.5 mag in F146;
it is stated as a limitation rather than tuned away.

#### Draft text: what the Marshall et al. (2006) 3D map is

*(Many readers know 2D extinction maps -- SFD/Planck dust emission, or the VVV and Gonzalez et al.
2012 bulge maps -- which give one number per line of sight: the TOTAL column, or the column to the
red clump. A 3D map gives the extinction as a function of DISTANCE along each line of sight, which
a simulation needs because its sources sit at different distances.)*

Marshall, Robin, Reylé, Schultheis & Picaud (2006, A&A 453, 635) built a 3D extinction map of the
inner Galaxy from 2MASS near-infrared photometry. The idea: in each direction, the Besançon model
of the Galaxy predicts how many stars of each type lie at each distance and what their intrinsic
J-Ks colours are. The observed stars are redder than predicted by their colour excess E(J-Ks),
which grows with the dust in front of them. Matching the observed colour distribution to the
model's, distance bin by distance bin, gives the extinction A_Ks as a function of distance. It
covers |l| <= 100 deg, |b| <= 10 deg on a 15-arcmin grid (over 64,000 lines of sight), and reaches
the bulge because the near infrared is ~10 times less extinguished than the optical (A_Ks ~ 0.08
A_V): stars behind A_V ~ 30 are still detected in Ks. Its limits: 15-arcmin resolution, coarser
than DECaPS; and its distances come from a Galaxy model rather than from the stars themselves.
The authors show the result is insensitive to moderate changes in the model. (The Besançon model
is also the stellar-population model behind this simulation's synthetic colour-magnitude
diagrams.)

How it is used here: Marshall tabulates A_Ks; we put it on DECaPS's A_V scale with
k = A_Ks/A_V = 0.0830, measured where both maps are valid (the median ratio at 8 kpc over Roman's
five-field block, 580 positions). It is used only where DECaPS stops seeing -- beyond DECaPS's own
reliable distance, and outright once its calibrated column passes DECaPS's sensitivity limit
(A_V = 12), which happens at 53% of table positions, typically from ~4.5 kpc.

