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
| M1 | Roman: AB->Vega before the astrometric error curve (x0.45 at the median source); photometric curve anchored to STScI 66-s depth 25.45 AB, interpolated; depth/saturation sourced | 72 | 5b10cdb |
| M1b | Rubin: per-visit depth/saturation gate (visit fiveSigmaDepth, -8.3); FWHM = median seeingFwhmGeom of the bulge visits; gamma from Ivezic 2019 | 73 | 0e22ab6 |
| M2 | Luminous lenses (bulge population, MS stars) add their light to the blend; astrometric centroid light-weighted (source, lens, blend) | 74 | 7a572b3 |
| M3 | Astrometric derivative steps swept; central stencil | 75 | 0a518f0 |
| M4 | Telescopes with < 3 epochs left out of the photometric fit (fixes the x8.3 mass event); gap geometry (t0zone, dt_edge) from the observed peak t0obs | 76 | (next) |
| 3 | Astrometric reference position freed; sqrt(2) per-coordinate fix; correlated-floor bracket | 71 | 39b50c0 | every theta_E / mass number |

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

### Step 3 -- the astrometric floor (2026-10-01, Deviation 71, commit 39b50c0)
- **Two corrections to state plainly (they move every theta_E and mass number):** (A) the fit now
  solves for the source's reference position, as every real astrometric fit does (Lam+2026 Eq. 1-2);
  before, theta_E was partly "measured" from the source's absolute position. Fixture: sigma(theta_E)
  x1.4-10 larger. (B) per-coordinate errors: the sources quote 1D (x or y) precisions; the old code
  doubled each coordinate's variance. sigma x1/sqrt(2).
- **The floor, in one paragraph for Sec. 5:** Roman's per-exposure centroiding floor of 1% of a pixel
  (1.1 mas per coordinate; Sanderson+2019, Lam+2026, McKinnon & van der Marel 2026) is everywhere
  treated as white noise that averages down as sqrt(N). The GBTDS's several-pixel and sub-pixel
  dithers sample the pixel phase, which is what makes that plausible (Lam+2026). Correlated
  systematics -- distortion residuals (few x 0.1% pixel; Bellini 2024), crowding biases (~mas, fixed
  per roll), frame alignment -- are acknowledged in all three papers and quantified in none
  ("it is not clear how best to incorporate these effects", McKinnon & van der Marel). We therefore
  report a bracket (Step 3c). Note also: a time-CONSTANT error is absorbed by the reference position
  and costs nothing; only correlations on timescales shorter than the event matter.
- Literature PDFs/text used: arXiv 1712.05420, 2608.24998, 2602.00310, 2601.10789 (refs to add to
  Report/refs.bib: McKinnon & van der Marel 2026 PASP 138; Kaczmarek+2026 A&A; Lam+2026 already cited).
- **Implemented (3c) and piloted.** Report it as: theta_E and masses are quoted under the white-noise
  assumption the literature uses (W), with the nominal (N) and pessimistic (P) bracket beside every
  number. Pilot (bulge lenses, stub, unweighted): sigma(theta_E) N/W ~1.5, P/W ~12; theta_E to 10%:
  7.9% / 1.3% / 0% (W/N/P). The old report's astrometric headlines (H5, Sec. 4.7; bh/ns masses in
  4.10; yields of masses to 10% in 4.12; image-resolution bars use sigma_a too -- check R1's D*sigma_a
  bars, which use erra per coordinate and are unaffected by fix A) must all be redone.
- **Figure worth making:** sigma(theta_E)/theta_E CDF for W/N/P, per population (bulge/ns/bh), from
  the production tables; and the fraction to 10% vs lens mass for the three variants.
- **Method box for Sec. 2.3 (three Fisher matrices):** the astrometric matrix now has a free
  reference position per frame (telescope; per roll in N; per season in P) and day-correlated noise;
  a one-line flowchart addition.

#### Draft text: why three astrometric noise models (W, N, P), and why these three

*(For the report's Sec. 2.3 and Sec. 5. The argument is the point; the numbers are placeholders
until the production runs.)*

**What has to be decided.** Roman measures the Einstein radius theta_E from the tiny shift of the
source's centroid during the event (~0.1-1 mas for ordinary lenses, a few mas for black holes),
using ~50,000 exposures, each with a per-coordinate precision of 1.1 mas for bright stars and
several mas for typical bulge sources. The forecast therefore rests on how the per-exposure errors
COMBINE. If they are independent, N exposures beat the error down by sqrt(N) ~ 220; if part of
the error is shared between exposures, that part does not average down.

**What is known.** The 1.1 mas is a centroiding floor of 1% of a 0.11-arcsec pixel, the level
demonstrated with HST (Sanderson et al. 2019; Lam et al. 2026; McKinnon & van der Marel 2026). All
published Roman forecasts, including the ones we compare with, add it as independent ("white")
noise, and there is a physical reason to expect most of it to behave that way: the GBTDS will
dither by several pixels and by sub-pixel steps, so a star falls on a different part of a pixel
and of the detector from one exposure to the next, and the pixel-level errors that set the floor
are re-drawn each time (Lam et al. 2026). But the same papers name systematics that would NOT be
re-drawn -- residuals of the geometric-distortion solution (expected at "a few x 0.1% of a pixel",
Bellini 2024), biases from crowding (a neighbour can shift a faint star by ~2 mas in HST bulge data),
time-dependent instrument effects (readout, detector, plate-scale drifts), and frame alignment --
and none of them gives the size or the timescale of the shared part. McKinnon & van der Marel say
it outright: "it is not clear how best to incorporate these effects without specific knowledge
about the true Roman performance." No choice of a single number can be defended, so we carry a
bracket, and choose its members by what is physically distinct rather than by arbitrary scaling.

**The key observation that makes a bracket possible.** A real astrometric fit always solves for the
source's reference position (its position at some epoch). Any error that is the SAME in every
exposure of a frame -- a fixed bias in that frame's zero point -- is absorbed by that reference
position and costs nothing. So the question is not "is the floor correlated?" but "over what
TIMESCALE is it correlated, compared with the event?" Two timescales matter physically, and the
three models are three answers to them:

1. *Frame timescale (offsets).* How often is the astrometric frame effectively re-set? Each
   re-set adds a free offset, and each free offset throws away the information in the absolute
   position between frames.
2. *Short timescale (within a day).* Do exposures taken close together share an error that
   exposures on different days do not? The day is the natural unit: it holds ~120 GBTDS exposures
   cycling through one dither sequence under one thermal and guiding state, and Lam et al. bin
   Roman astrometry to one day for the same reason. Such a shared error is modelled as a
   per-coordinate sigma_c common to all of a day's exposures (it averages down across days, not
   within one).

| model | frame offsets (free) | shared error per day | what it represents |
|---|---|---|---|
| **W** white | one per telescope | 0 | the literature's assumption: dithering makes the floor white; the frame is stable for the whole mission (Roman at L2). Used for every published comparison. |
| **N** nominal | one per Roman **roll** (spring / autumn) | 0.3 mas | the two systematics with a physical anchor: the telescope turns 180 deg between spring and autumn, so the PSF and the detector layout flip and any crowding or PSF-model bias changes sign -- a separate frame per roll; and distortion residuals at the level Bellini (2024) expects (0.3 mas = 0.27% of a pixel, "a few x 0.1%"), shared by a day's exposures. |
| **P** pessimistic | one per Roman **season** (10) | 1.1 mas (the whole floor) | every season is calibrated independently (plate-scale and distortion drifts between seasons of ~6 months, which Lam et al. flag as time-dependent systematics), and the floor does not average down within a day at all -- a day of ~120 exposures is no better than one. This contradicts the HST experience that 100 stacked exposures reach ~0.1 mas (Sanderson et al.), which is why it is a pessimistic bound rather than an expectation. |

The models are nested: N adds nuisance parameters and correlated noise to W, and P adds more of
both to N, so sigma_W <= sigma_N <= sigma_P for every event; the code asserts this. Rubin's
astrometry is the same in all three (one frame, white noise): Rubin contributes little to theta_E
and its own systematics are a separate question.

**What the bracket does not cover.** (i) Extra white scatter from crowding (Whitaker et al. find
a factor 2-5 for faint stars near bright ones) would scale every model's per-exposure error, not
change their ratios; (ii) an error correlated on timescales between a day and a season (e.g. a
slow drift within a season) lies between N and P; (iii) the true split can only come from Roman
data or from the Roman project's own astrometry simulations, and a measured split would replace the
bracket.

**How to quote results.** Main numbers in W, so they compare directly with published forecasts
(Sajadian & Sahu 2023; Lam et al. 2026; Kaczmarek et al. 2026), with N and P given beside every
astrometric number (theta_E precision, lens masses, yields of masses to 10%). The pilot already
shows the spread matters: on bulge lenses theta_E to 10% falls from 7.9% (W) to 1.3% (N) to 0% (P).

**And the two corrections underneath all three.** Independently of the bracket, two errors in the
earlier forecasts were fixed: the reference position was not a free parameter (the fit could use
the star's absolute position, making sigma(theta_E) a median 6.8x too small in the pilot), and the
per-coordinate error was counted twice in variance (sigma sqrt(2) too large). The first dominates;
the earlier reports' theta_E and mass precisions were optimistic for this reason, not because of
the floor.

### M1 -- Roman's magnitude system and noise (2026-10-02, Deviation 72)
- **Say in Sec. 2.1 (noise):** all magnitudes are AB (MIST bolometric corrections). Roman's
  photometric error: Penny et al. 2019's single-epoch F146 curve (1 mmag floor), anchored to STScI's
  current 5-sigma point-source depth for a 66-s exposure, 25.45 AB (25.37 in 57 s); saturation 14.8 AB
  (Penny et al. 2019). Roman's astrometric error: Lam et al. 2026's per-exposure curve, whose anchors
  are Vega magnitudes, evaluated at m_Vega = m_AB - 1.03 (synphot, STScI throughput, CALSPEC Vega).
  Per-exposure precision at the median Roman detection (F146 ~ 21.8 AB): ~1.2 mas per coordinate.
- **Correction to quote:** the earlier forecasts used AB magnitudes in a Vega-calibrated curve,
  overstating Roman's astrometric errors ~2.2x -- in the opposite direction to the reference-position
  correction (Step 3). Net effect on theta_E: measure on the new runs.
- **Decision to flag to the reader (and the user):** Penny's curve was rescaled to the current STScI
  depth rather than by exposure time alone (which would be 0.26 mag deeper).

### M1b -- Rubin's survey constants (2026-10-02, Deviation 73)
- **Say in Sec. 2.1 (noise and detection):** each Rubin epoch is recorded if the source is between
  that visit's saturation and its own 5-sigma depth (OpSim fiveSigmaDepth, baseline_v5.1.0), and its
  error follows Ivezic et al. 2019 eq. 5 at that depth with a 5 mmag floor; blending and the PSF bar use
  the median geometric seeing of the bulge visits (0.90-1.11" by band). Table of per-band medians
  (depth 23.3 / 24.3 / 23.9 / 23.4 / 22.9 / 22.0, seeing) is worth including -- they are what Rubin
  delivers toward the bulge, below its all-sky design numbers.

### M2 -- luminous lenses and the blended centroid (2026-10-02, Deviation 74)
- **Say in Sec. 2.2 (one event):** in the ordinary-star population, lenses below the turnoff are
  main-sequence stars whose light (magnitudes from the Besancon model's mass-luminosity relation with
  MIST bolometric corrections, dimmed by the dust in front of the lens) is blended with the source;
  black holes and neutron stars are dark. Measured astrometry is the light-weighted centroid of
  source, lens and blends. Pilot: a luminous lens supplies a median 24% of the F146 baseline light
  (0.1% in r).
- **Physics worth a sentence:** lens light both dilutes the shift and pulls the centroid toward the
  lens; it is why theta_E for ordinary lenses is harder than for black holes beyond the sqrt(M) factor.

### M3/M4 (2026-10-02, Deviations 75-76)
- **Gap filling (Sec. 4.3):** "in a gap" is now judged by the OBSERVED peak time (parallax included),
  not the model's reference t0 -- 10% of detections change zone in the pilot. Say so in the method.
- Numerical: astrometric derivatives verified step-independent (sweep); no text needed beyond a line.

