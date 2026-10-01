# OPEN_ITEMS.md

Running list of known-but-deferred discrepancies surfaced while working through
`JOINT_FIT_REFACTOR_PLAN.md`. Items are removed once resolved; unresolved items are carried
into the whitepaper's open-items discussion at Phase G, Step G3.

---


## AUDIT 2026-10-01: every open entry re-checked against the code, the reports and the sources

The user warned that entries here can be outdated or wrong (several were, in the overview-report
work). Each open entry below was re-checked on 2026-10-01 against the CURRENT code (commit 39b50c0),
the overview report, and the primary sources. Verdicts:

| entry | verdict | evidence |
|---|---|---|
| Astrometric Fisher CHECK aborts | **RESOLVED** (outdated) | the CHECK no longer exists; null derivatives are accepted (FisherM, Deviation 71 rewrite) |
| Joint flag_det is Rubin-only | true, diagnostic only | `flag_det` set only in the Rubin branch (l. ~2002), written to stdout only |
| §5.1 visit-list model (whitepaper) | true, document only | whitepaper still describes one merged list; Phase G |
| Delta2[] unswept, tetE/piE biased stencil | **true, now important** | no sweep in DEVIATIONS; tetE/piE still use sig2; theta_E is now the headline |
| fb stencil depends on fb bin | true, negligible | co.bb[] bins unchanged; C3 showed < 0.3% effect |
| table still test2.dat in repo root | **partly outdated** | now `./test<tag>.dat` in the run directory; still truncated by a fresh run there |
| stale TODO(Ali) (RomanBaseline read) | true | l. ~885; also a stale "no Roman astrometric error model exists yet" TODO at l. ~988 (wrong since H4) |
| flagi stale on uncharacterised events | true | not in the per-event reset block |
| two output files in append mode | **partly outdated** | EfLMC/EfLMC_B/magC0/datC0 now truncate on a fresh run; MapLMC and LpLMC still always append |
| startup check on an unread file | **RESOLVED** (outdated) | the `fil0` input stream is gone; LpLMC is created if missing |
| numd[0] counts drawn stars | true | records.push_back (l. 2489) sits outside the visibility gate (closes l. 2486) |
| BH remnant single slope | true, model choice | helper.cpp remnantMass, BH_MI_SLOPE |
| five events sigma_joint > sigma_Rubin, cond > 1e9 | historical; re-check on new runs | kMaxCondition 1e12 unchanged |
| binary git stamp stale | true | Makefile captures GIT_COMMIT at make time; no forced rebuild |
| F3 rests on 1,950 events / F2 long-tE on 124 | superseded by the stratified runs | re-check sample sizes on the new runs |
| F1/F2/F3 OOM | **probably RESOLVED** | romanlib.load_events streams in chunks with `keep`; F1 uses it |
| two Phase F panels unweighted | **RESOLVED** (outdated) | f3 and f4 call romanlib.attach_weight |
| E1 tE stratification not implemented | true, deferred by design | |
| Roman halo orbit not modelled | true, negligible | |
| no astrometric-shift analysis product | **RESOLVED** (outdated) | analysis/h5_astrometric_shift.py, h5_astrometry_summary.py, h5_crosscheck.py |
| astrometric shift not diluted by blending | **true, now important** | lightcurve(): pos1c adds def1c undiluted (l. ~4000) |
| provenance stamp stale (discipline) | true | same as the git-stamp entry |
| --dry-run destroys previous output | **true** | outputs are opened/truncated (l. ~1080-1141) before the dry-run exit (l. ~1620); a NOTE is now printed |
| dchiP/dchiA absolute | true, diagnostic | l. ~2282 |
| startup ~170 s parsing extinction files | **RESOLVED** (Deviation 70) | single table file; dry run 26 s |
| H7 cost mechanism | superseded (marked so) | |
| G2 has no published number | true, deferred | |
| pooled stats / per-mass efficiencies | RESOLVED (already marked inside) | |
| H3 satellite numbers unweighted | open until a new paired run | the new runs use --pair-satellite |
| map stream never flushed | code FIXED (marked inside); data note only for v3 | v3 data deleted 2026-10-01 -> entry can close |
| Roman photometric error placeholder | **true, and its source is now identified** | sigma_roman.txt = Penny+2019 Fig. 4 (W149 = F146, AB, 46.8-s Cycle-7 exposure, 1 mmag floor); nearest-neighbour lookup; the adopted GBTDS exposure is 66 s |
| magC0/datC0 dump dead | true | gate `save < 0` with save = 0 can never fire |
| t0 is not the observed peak | true | lightcurve() references parallax to Earth at t = 0 |
| CRITICAL AlAv inverted | FIXED IN CODE (Dev. 53), closes with the re-runs | |
| event rate too flat in latitude | open, re-measure on new runs | measured through the old dust |
| Neven ~300x | open, analysis | |
| one DET_ANOMALY | open | code path exists (nDetClass[DET_ANOMALY]); not reproducible until per-sightline seeding |
| Roman detects off-mission peaks | open: a reporting convention | measured 2.9-4.4% |
| x8.3 sigma_joint > sigma_Roman on Ml | **open; the entry's proposed cause is incomplete** | with okA/okB set on both partitions, min-of-two-routes is monotone (each route is); a violation needs a partition whose matrix was rejected (condition cut) -- investigate the event, do not presume a fix |
| Penny comparison | open, after the runs | |
| no --seed / --end-index | true | seed = 42 compile-time; --start-index skips without consuming RNG |
| F146 depth thre[6] = 29 placeholder | true | Bulge.h; the curve's own 5 sigma is 25.52 (AB, 46.8 s) |
| stale statements in older documents | true, documents | |
| Rubin FoV a circle (new) | true; LSSTCam facts to be re-verified when worked on | |
| GC dust 10% thin (new) | true (measured, Deviation 70) | |

**The hard-coded survey constants in Bulge.h (user's request, same day)**, checked against Ivezić
et al. 2019 (ApJ 873, 111; arXiv:0805.2366, Tables 1-2), the OpSim database the visit list comes
from (baseline_v5.1.0, the 12,308 bulge visits), Lam et al. 2026 and STScI:

| constant | value in code | verdict |
|---|---|---|
| `thre[0-5]` (u..y depth gate) | 23.4 24.6 24.3 23.6 22.9 21.7 | **= the SRD MINIMUM single-visit 5-sigma depths** (Ivezić Table 1, "min."). Not what these visits reach: median per-visit m5 of the bulge visits 23.31 24.33 23.93 23.36 22.86 22.00 (0.1-0.4 mag shallower in u-z, 0.3 deeper in y). The noise model (errlsstM) already uses each visit's own m5, so the recording gate and the error bars disagree. |
| `satu[0-5]` (saturation) | 15.2 16.3 16.0 15.3 14.6 13.4 | **= thre - 8.3 in every band**, i.e. derived, not sourced. Ivezić quotes only "the LSST saturation limit at r ~ 16" (r 16.0 agrees). |
| `gama` | 0.037 0.038 0.039 0.039 0.040 0.040 | an older version of Ivezić Table 2's gamma (0.038 0.039 0.039 0.039 0.039 0.039); < 3% in the error. 6 values for M = 7 (gama[6] = 0, unused). |
| `delta2` = 0.005 mag | systematic photometric floor | **verified**: Ivezić requirement 3, "photometric repeatability should achieve 5 mmag precision at the bright end". |
| `FWHM[0-5]` (Rubin PSF; also sets the blending disc) | 1.221 1.101 0.993 0.967 0.952 0.937 | close to the GEOMETRIC PSF width of the actual bulge visits (median seeingFwhmGeom 1.114 1.042 0.982 0.949 0.932 0.899): an older OpSim's value, 1-10% wide (u worst). Not Ivezić's zenith theta_eff. Per-visit seeing spans 0.75-1.5" (16-84%), but blending uses one value per band. |
| `FWHM[6]` (F146) | 0.105" | STScI's F146 PSF FWHM; **Lam et al. 2026 use 115 mas** -- check which (centre vs field average) before relying on the PSF resolution bar. |
| `sigma[]` (per-band extinction scatter) | 0.022 ... 0.04 | unsourced; F146 value commented "PLACEHOLDER: K-band value". Small (0.02-0.04 mag) next to the dust maps' own errors. |
| `seeing`, `msky`, `Cm`, `Dci`, `km`, `cade1` | | **unused** anywhere; older than Ivezić 2019. Dead constants. |
| `LSST_AST_FLOOR` = 10 mas | per visit per coordinate | **verified** (Ivezić: "10 mas per observation per coordinate"). |
| `RUBIN_REF_BANDS` = {2} | one Rubin blend fraction and baseline (r) for all six bands in the Fisher fit | a modelling simplification, documented in Bulge.h; a real fit has per-band fb and mbs. Not sourced, not wrong; state it in the report. |
| LSSTCam facts in the new Rubin-FoV entry | 21 rafts, 189 CCDs, 9.6 deg^2 | **verified** (Ivezić Sec. 2.6.2, Table 1); the "fill factor ~0.9" is NOT in Ivezić -- re-verify in M8. |

**Three NEW problems found by this audit** (entries below): the AB/Vega mismatch in Roman's
astrometric error, the unsourced F146 saturation limit, and the 46.8-s exposure of the photometric
table. The overview report repeats one of them: its "2.75 mas at the median F146 = 21.82" evaluates
the Vega-calibrated curve at an AB magnitude.

## Astrometric Fisher CHECK aborts when a perturbation is unresolvable (from the fixture)

`FisherM`'s astrometric branch asserts

    CHECK(!(s.pos1c == l.soux[i] && s.pos2c == l.souy[i]));

which aborts the whole program (uncaught `std::runtime_error`) whenever a parameter perturbation
leaves *both* modelled source coordinates bit-identical to the stored ones. Two ways that happens:

1. **An epoch lands exactly on `t0`.** There the source-motion terms `mus1*(timh - t0)` and
   `mus2*(timh - t0)` vanish identically, so perturbing `mus1` or `mus2` changes nothing at all.
   Found immediately by `tests/fisher_fixture.cpp`, whose regular cadence grid hits it easily.
   In production `t0` is drawn continuously so exact coincidence is measure-zero -- but not
   impossible, and it would be an unexplained hard crash if it ever happened.
2. **An epoch enormously far from peak.** The astrometric deflection falls off as 1/u^2, so at
   large `|t - t0|/tE` a `piE` perturbation can change the modelled position by less than a
   double can represent.

Neither case is a modelling error -- they are legitimate points where a particular derivative is
zero or unresolvable. Aborting the program is the wrong response.

**What's needed:** replace the assert with per-parameter handling -- detect a null derivative for
that (epoch, parameter) pair and skip its contribution, rather than testing the AND of both
coordinates and crashing. This belongs with Step C4's "explicit not-characterizable outcome
instead of crash-or-garbage" work, which is the same class of problem.

**Workaround meanwhile:** the fixture's `t0` values are deliberately offset off the cadence grid,
and it carries an explicit guard that reports this cause rather than leaving an opaque abort.

## Joint `flag_det` is Rubin-only (from Step C0)

`flag_det` (`Bulge_LSST.cpp`, in the Rubin epoch branch) is described as the joint run-test flag,
but it is only ever set inside the **Rubin** branch and is gated on `ndw_L`. There is no equivalent
update in the Roman branch, so an event with a strong persistent Roman signal and no Rubin signal
leaves `flag_det` at 0.

**Why this is not currently a live bug:** `flag_det` does not feed the detection decision. `detL`,
`detR`, and `detJ` are built from `flag_det_L` and `flag_det_R` (the explicit per-survey flags added
in Step B1), and `FFG[0]` is their union. `flag_det` is only written to the diagnostic output stream.

**What's needed:** either make it a genuine joint flag (set from either branch) or retire it and stop
writing it to output, since as it stands the column is mislabeled — a reader would reasonably take it
for a joint quantity. Decide when Step D1 reworks `EventRecord` and the output table.

## §5.1 visit-list model (from Step A2)

`whitepaper.tex` §5.1 still describes a single merged Roman+Rubin visit list with `FoV = 1.75°`.
The code now uses two separate visit lists and cadences: `FoV` (1.75°) for Rubin via
`BulgeBaseline.dat`, `FoVRoman` (0.28°, currently a placeholder — see the `TODO(Ali)` next to its
declaration in `Bulge.h`) for Roman via `RomanBaseline.dat`, matched independently by two calls to
`matchVisibleEpochs()`. Fix in Phase G, Step G3.

## Astrometric finite-difference steps (`Delta2[]`) are unswept, and two use a biased stencil (from Step C3)

Step C3 swept and retuned the **photometric** steps (`Delta1[]`) and fixed the stencil for the
photometric `tE` and `piE` (see `DEVIATIONS.md` entries 10 and 11). The astrometric matrix was
left untouched, and has both of the same problems:

- **`Delta2[]` has never been swept.** Its steps came from the same legacy codebase as the
  photometric ones, which turned out to be ~4 orders of magnitude too large. There is no reason
  to assume the astrometric ones are better placed, and every reason to expect they are not.
- **`tetE` and `piE` still use `sig2`**, the first-order biased stencil (`Bulge_LSST.cpp`, the
  `Delta2[j] * sig2[h]` lines). `mus1`/`mus2` correctly use the central `sig`.

**Why it was not fixed at the same time:** the effect could not be verified. The astrometric
branch also depends on `errlsstA()` standing in for a real Roman F146 astrometric error model
(first item in this file), so its sigmas are not yet trustworthy in absolute terms regardless of
the stencil. Fixing one of two unverifiable inputs in isolation buys nothing checkable.

**What's needed:** extend `runSweep()` in `tests/fisher_fixture.cpp` to sweep `Delta2[]` (the
harness already carries `deltaScale` and the CSV/plot pipeline; `covarian` would need the
astrometric equivalent), then retune and switch the two `sig2` uses to `sig`. Expect the same
qualitative outcome: absolute astrometric sigmas shrink toward their true, larger values, while
the paired joint/single ratios hold up.

**Where this should land:** before Step G1, which reruns the astrometric matrix to separate
satellite parallax from temporal-baseline parallax, and before any absolute `tetE` precision
number goes in the whitepaper.

## `fb` derivative stencil depends on which blend-fraction bin the event lands in (from Step C3)

The blend-fraction step `co.bb[]` (`Bulge_LSST.cpp`, inside `FisherM`'s data loop) is chosen by
binning on `s.fb[tt]`: the middle bin gives `{-0.07, +0.07}`, a proper central difference, but the
outer bins give `{+0.07, +0.15}` and `{-0.07, -0.15}` — two forward (or two backward) differences,
carrying the same first-order bias as `sig2` in `DEVIATIONS.md` entry 11. So `fb`'s derivative
accuracy depends on where `fb` happens to sit, which is not a property anyone chose deliberately.

**Why it is not urgent:** the C3 sweep showed `fb0`/`fb1` flat to <0.3% even at the *unscaled*
steps, and `kFDStepScale` now shrinks them by 1e-4, making the residual bias negligible. The bin
structure and the bound-safety clamp are also now far from binding, since a step of ~7e-6 cannot
push `fb` out of [0,1].

**What's needed:** replace the binned step with a symmetric `{-h, +h}` pair at a fixed small `h`,
keeping a clamp only as a guard. This would also let the bin table and the clamp block be deleted,
simplifying `FisherM`'s data loop. Do it when `FisherM` is next opened for other reasons.

---

## RESOLVED 2026-10-01 (Deviation 69) — Sky coverage is Penny et al.'s, not the current GBTDS footprint — and Rubin's FoV overlap is unmodelled (raised by Ali, 2026-08-24)

**Resolution.** The simulator now uses the adopted layout (`Baseline/gbtds_layout/`, Penny's
`gbtds_{spring,autumn}_2026.4.3`, centres confirmed by STScI's GBTDS page): 18 detectors per field,
spring and autumn rolls, coverage decided per visit by "on a detector". The scan region is a
distance rule around the 12 placements and the Rubin visit list was rebuilt for it (12,308 visits).
**The Rubin partial-overlap concern below does not apply to this code**: coverage is decided per
sightline, which is a point -- a Rubin pointing either images it or not, however little of the
pointing's field overlaps a Roman field -- so partial overlap is handled exactly at the grid's
resolution; blending (`s.blend[]`, `s.fb[]`) is computed per source, not per pointing. What
remains is that Rubin's field is a circle, not LSSTCam's real outline (new entry below). Kept for
the record:

**MEASURED 2026-09-29 (Deviation 59):** against the real tile layout in
`Whitepaper/roman_967_both_aladinX.png` (spring + autumn), the modelled five-field block sits at
b = -1.2 where the real one is centred at b = -1.40 (spans -1.80..-0.97), and the modelled GC field
at (0, -0.125) against the real (0.06, -0.22). 79% of the simulated footprint lies on real
detectors; it covers 59% of the real tiles (1.97 deg^2 over both seasons, 1.71 per season), missing
mostly b < -1.5. **Effect estimated (Deviation 60):** Roman footprint yields change by only +1 to +5%
(area +16.5%, lower density farther from the plane -10 to -13%); Rubin's footprint detections -1 to
-10%; per-event fractions ~unchanged. A run on the adopted layout is still the remedy. Figure:
`figures/footprint_20260929/footprint_vs_gbtds.pdf` (`Report/overview/make_footprints.py`).

**What is wrong:** the region scanned by the Monte Carlo (`l1`/`l2`/`b1`/`b2`/`wid` in `Bulge.h`)
and the Roman field centres (`FIELDS_L_B` in `Baseline/generateRomanBaseline.py`) both descend
from the Penny et al. survey design. That design is now superseded — STScI's published GBTDS
footprint has changed since. The authoritative source is the Roman documentation already recorded
in `DEVIATIONS.md` entry 16:
<https://roman-docs.stsci.edu/roman-community-defined-surveys/galactic-bulge-time-domain-survey>

Note this is a *different* correction from entry 16.4. That fixed `FoVRoman`'s units (an area used
as a radius); the field **centres and the total footprint** were left as they were.

**The subtlety that makes this more than a coordinate update:** Rubin's field of view is very much
larger than Roman's (`FoV` vs `FoVRoman` in `Bulge.h`). A Rubin pointing whose centre lies well
outside a GBTDS field can still cover part of that field. The current `matchVisibleEpochs` treats
both instruments identically — a sightline matches an epoch when the sightline is within the
instrument's radius of the pointing centre — so this partial overlap is either counted as full
coverage or missed entirely, depending only on centre separation. Two consequences:

1. **Joint coverage is mis-stated.** Sightlines inside a GBTDS field that Rubin observes only via
   the edge of a large pointing are exactly the ones the joint-fit science case depends on.
2. **Blending is wrong there too.** Blend fractions (`s.blend[]`, `s.fb[]`) depend on how much
   Rubin actually sees at that position, so a coverage error propagates into the photometry and
   from there into `fb0`/`mbs0` and the Fisher forecast — not just into the event counts.

**Also required:** the Rubin baseline (`Baseline/readbaselineBulge.py` → `BulgeBaseline.dat`, and
`Nl` in `Bulge.h`) must be regenerated for whatever region the corrected footprint implies. The
two baselines share the simulation clock and the scan region, so they cannot be updated
independently.

**Why it is not being done now:** it changes what sky the forecast describes, so every detection
and precision number would have to be re-taken afterwards. It should be done as its own step, with
the new footprint read off the STScI documentation rather than inferred, and it interacts with the
run-scaling grid (Step 4) — the stride guard checks the grid against `FIELDS_L_B`, so the guard
protects the change but the field list itself is what needs replacing.

**Deliberately deferred at Ali's request, 2026-08-24.** Wanted, not optional.

---

## The per-event table is still called `test2.dat` and still lives in the repo root

**What.** The 88-column analysis table that Step D1 built is written to `./test2.dat`, a name left
over from when it was a debug dump. It sits in the repo root next to the binary, is truncated by
every run, and is `.gitignore`d, so a result can be overwritten by the next smoke test with
nothing to say it happened.

**Where it should go.** `files/MONTLMC/files/`, alongside the other outputs and `run_provenance.txt`
— ideally named for the run, so a production run and a stub run cannot clobber each other.

**Why not now.** `tests/c3_live_compare.py` reads it by path and by position, and is being actively
edited. Renaming mid-flight would break a comparison in progress for no scientific gain. It is a
five-minute change whenever that script settles.

**Note for whoever does it:** that script's `read_csv(..., header=None)` will need
`comment="#"` added, because D1 gave the file a header line. Its column-count assertion fails
loudly first, which is the intended behaviour.

---

## `TODO(Ali)` at `Bulge_LSST.cpp` (RomanBaseline.dat read) is stale

It asks for a generator sourced from Roman's own season structure rather than an LSST OpSim.
Commit `74e5e18` delivered exactly that (`Baseline/generateRomanBaseline.py`, against the ROTAC
2025 / STScI GBTDS design). The comment should be deleted; left in place only because removing it
is cosmetic and belongs with whatever step next touches that read.

---

## `co->flagi` is stale on uncharacterized events

`flagi` is set inside `FisherM`, which only runs for detected events, and it is **not** in the
per-event reset block that clears `okA`/`okB`/`condA`/`condB`/`Era`/`Erb`/`relMl`. So the `flagi`
column of the per-event table reads `1` on rows where nothing was characterized — it is the
previous characterized event's value.

**Harmless today**, because the only consumer (the per-field precision average) pairs it with
`co->okA[SJOINT]`, which *is* reset. But it is a trap for anything reading the table: **use
`okA_J`, not `flagi`, as the characterizability flag.**

Not fixed in Step D1 because `flagi` is also read inside `FisherM` itself and adding it to the
reset changes behaviour rather than only bookkeeping — it needs its own look at what `flagi` is
actually supposed to mean.

---

## Two output files open in append mode, so a re-run silently doubles them

**What.** `LpLMC2.dat` and `MapLMC2.dat` are opened with `std::ios::app`, not truncated:

```
Bulge_LSST.cpp:496    std::ofstream fil3(fnGam, std::ios::app);      // MapLMC2.dat
Bulge_LSST.cpp:1149   std::ofstream fil0_append(fnLDt, std::ios::app); // LpLMC2.dat
```

Every other output (`EfLMC2.dat`, `EfLMC2B.dat`, `magC0.dat`, `datC0.dat`, `test2.dat`) is
truncated at startup. These two are not, and nothing in the run says so: a second run appends its
rows to the first run's and the file ends up holding two runs' worth of events with no marker
between them. A stub run followed by a production run produces a `MapLMC2.dat` that is 36
sightlines of diagnostics glued to the front of the real dataset.

**How it was found.** Launching the Step D1 production run, 2026-08-29. The files had to be
cleared by hand first, and there is nothing in the code, the provenance block or the output that
would have caught it if they hadn't been.

**What it should be.** Either truncate them like every other output, or — if the append is
deliberate, e.g. accumulating across `IMnum` values the way `BHLSSTMONTS.dat` is cleared only when
`IMnum == 1` — say so in a comment and record the pre-existing row count in `run_provenance.txt`
so a downstream reader can tell where this run's rows begin.

**Why not now.** It needs a decision about whether the append was ever intentional, which means
looking at what `IMnum` is for. Changing it blind risks breaking a workflow that relies on the
accumulation.

---

## The startup file check makes an unread file's existence a precondition, and names no file when it fails

**What.** `LpLMC2.dat` is opened twice. Line 490 opens it for *reading*; nothing ever reads from
that stream. Its only purpose is to be tested at line 523:

```
Bulge_LSST.cpp:490    std::ifstream fil0(fnLDt);
Bulge_LSST.cpp:523    if (!fil0 || !fil2 || !fil2b || !fil3 || !fil4 || !fil5) {
                          std::cerr << "Cannot open one or more files!" << std::endl;
```

The actual writing is done by a separate append stream at line 1149. So the file must **exist** for
the run to start, even though its contents are never used — and if it does not, the run dies after
reading the baselines, the extinction maps and the whole CMD set, several minutes in, with a
message that names none of the six files it tested.

**How it was found.** The first attempt at the Step D1 production run, 2026-08-29, died here after
~4 minutes because `LpLMC2.dat` had been deleted to clear the append (see the item above). The log
showed the schedule guard and `read_cmd` succeeding, then `Cannot open one or more files!` with no
indication of which.

**What it should be.** Drop `fil0` if the existence of that file is genuinely not a precondition,
and check each stream separately with the filename in the message — the same discipline the Roman
schedule guard already follows, which prints the margins that failed rather than a bare verdict.

**Why not now.** Cosmetic on its own, and it belongs with whatever step next touches that block —
most naturally the `test2.dat` renaming item above, which moves output paths anyway.

---

## `numd[0]` counts drawn stars, not observable ones, so `EffiD` is 100% by construction

**What.** `icon` is incremented only inside the visibility gate:

```
Bulge_LSST.cpp:1149   if (flagf > 0 and ndw > 2) { //if star is visible
Bulge_LSST.cpp:1151       icon +=1;
```

but `records.push_back(...)` sits **outside** that block, so every drawn star gets a record, and
`numd[0]` — which counts records — is the number of stars *drawn*, not the number observable. The
two are equal only where every draw is observable.

**Consequences.**

- `EffiD = numd[0] * 100 / nsim` is described as "probability of detecting stars" but is 100% by
  construction, since a record is pushed on every draw.
- The per-sightline `[0]` means (`tE[0]`, `u0[0]`, `Ml[0]`, `fb[0]`, `mbs[0]`, `Map[0]`, `Ext[0]`,
  …) are averages over **drawn** stars, not observed ones.
- `EffiL = numd[1] * 100 / numd[0]` is detected-per-drawn, not detected-per-observable.

**How it was found.** `CHECK(icon == numd[0])` aborted the 2026-08-29 full-region run at
l = -3.499, b = -1.98, where 5 of 15 drawn stars were observable. It had never fired before because
every run prior to `81a6b04` used the dense 0.1x0.1 deg stub patch, where every draw *is*
observable and the equality holds by accident of the field.

**What was done.** The assertion was loosened to `CHECK(icon <= numd[0])`, which is the invariant
that actually holds. **No computed value changed** — only the assertion. The denominators are
untouched.

**What is still open.** Which denominator each quantity *should* use. If `EffiD` is meant to be the
fraction of drawn stars that are observable, it should be `icon / nsim`. If the `[0]` means are
meant to describe the observable population, the sums need the visibility gate. Both change
published numbers, so neither belongs in a bookkeeping fix — it needs a decision about what each
quantity is for, and a re-take of anything already quoted from them.

---

## Black-hole remnant masses use a single proportional slope

`remnantMass()` maps an initial mass above 20 Msun to a black hole by `Ml = 0.24 * Mi`,
giving ~4.8 Msun at Mi = 20 and ~28.8 at Mi = 120. That spans the observed stellar-mass
black hole range and is monotone, but it is the crudest link in the mass chain.

Real remnant masses depend on metallicity and on mass loss during the progenitor's life in
ways no single slope reproduces: at low metallicity weaker winds leave heavier black holes,
and the relation is not monotone across the pair-instability region. There is also no
attempt at a natal-kick velocity distribution, so black hole lenses share the kinematics of
their progenitors.

**Why it matters here.** Black holes are the long-tE tail, and the long-tE regime is exactly
where the joint fit earns its keep -- Roman's thetaE and Rubin's piE combining into a lens
mass. The *shape* of the black hole mass distribution therefore feeds directly into the
headline precision numbers, not just into a tail nobody looks at.

**Why not now.** Choosing an initial-final mass relation for black holes is a literature
decision with real spread between prescriptions, and it should be made deliberately rather
than folded into a bug fix. The white dwarf branch is on firmer ground (Kalirai et al. 2008,
calibrated on open clusters) and the neutron star branch uses the canonical 1.4 Msun, so
this is the one segment carrying an arbitrary choice.

---

## Five events have sigma_joint > sigma_Rubin, on matrices with condition number > 1e9

Adding data cannot worsen a Fisher forecast: the joint information matrix is the sum of the
per-survey ones, so `sigma_joint <= sigma_single` is an identity, not an approximation.
`romanlib.check_monotonicity` asserts it, and on the 2026-08-30 Kroupa run it reports five
violations out of 74,812 joint-detected events, in both `tE` and `piE`, with a maximum
excess of **1.0011** -- one part in a thousand.

All five have the same signature:

| tE (d) | piE | ndw_L | ndw_R | condA_J | condA_L | sigtE_J | sigtE_L |
|---|---|---|---|---|---|---|---|
| 10.16 | 0.111 | 2361 | 49972 | 8.55e10 | 8.53e10 | 284.05 | 283.75 |
| 10.39 | 0.180 | 1075 | 49634 | 1.64e11 | 1.64e11 | 54.32 | 54.32 |
| 12.21 | 0.172 | 893 | 46908 | 5.33e10 | 5.33e10 | 147.28 | 147.27 |
| 5.95 | 0.180 | 2360 | 49972 | 7.01e10 | 7.01e10 | 59.12 | 59.10 |
| 9.30 | 0.055 | 2371 | 49984 | 1.28e09 | 1.28e09 | 71.58 | 71.58 |

**This is round-off, not a partitioning bug.** A double carries about 16 significant digits;
inverting a matrix with condition number 1e11 loses about 11 of them, leaving ~5 -- so a
relative discrepancy of 1e-3 between two nearly identical inversions is the expected size,
not a surprise. Note also that every one of these forecasts is meaningless on its own terms:
`sigma_tE = 284 d` on a `tE = 10 d` event is not a measurement, and none of them pass the
characterization criterion.

**What is still worth doing.** The pipeline currently reports a sigma for any matrix GSL
manages to invert, with no conditioning floor. A cut -- refuse to report when
`condA > 1e8`, say, and set the sentinel instead -- would remove this class of artifact
outright and would also stop absurd sigmas propagating into medians. Choosing the threshold
needs a look at the condition-number distribution across the run, which is why it is
recorded here rather than applied.

Until then, `check_monotonicity` will keep printing this warning on every figure script. It
should stay noisy: the day it fires on a *well*-conditioned event, that IS a bug.

---

## The binary's git stamp goes stale whenever `make` has nothing to do

`Makefile:18` captures `GIT_COMMIT` at compile time and bakes it into the binary, which
writes it to `run_provenance.txt`. But the stamp only refreshes when something actually
recompiles. Editing a source file, building, and then committing leaves the binary carrying
the *pre-commit* label -- and a later `make` reports "Nothing to be done" and keeps it.

That is exactly what happened to the 2026-08-30 production run. `run_provenance.txt` and
every figure derived from it read `git_commit=d6dd293-dirty`, while the sources that were
actually compiled are the content of `5c74fbd` ("Lens the bulge with bulge stars, not
MACHOs"). **The science is unaffected -- the right code ran -- but the label on 2.5 GB of
output names the wrong commit.**

**Fix options, none yet chosen:** make the stamp a `.PHONY` prerequisite so every build
refreshes it; or have the binary refuse to run when `git diff --quiet HEAD` fails; or emit
the stamp at *run* time by shelling out, which costs a subprocess but cannot go stale.

**Interim rule for anyone running a production job:** `make clean && make` *after* the last
commit, never before, and check the `git_commit` line in `run_provenance.txt` against
`git rev-parse --short HEAD` before starting.

---

## F3's Roman-alone comparison rests on 1,950 events

The (`tE`, `piE`) characterization map's panel (a) -- "what Rubin adds to Roman" -- can only
use events inside Roman's footprint, because outside it Roman characterizes nothing for a
reason that is geometric, not physical. On the 2026-08-30 run that is **1,950 of 74,812**
joint-detected events, because only ~39 of 1,706 scanned sightlines fall inside the GBTDS
footprint at the current stride.

At 0.5 dex cells that leaves 16 coloured cells out of 56, most of the plane grey. The signal
in those cells is real (ratios up to 2.0) but the error bars are not small, and no cell
carries enough events to quote a trend within it.

**What would fix it:** sampling that concentrates draws inside the Roman footprint rather
than spreading them uniformly over the scan region -- which is Step E1's stratification
question wearing a different hat. A footprint-weighted run would buy roughly a factor of 20
more Roman-observed events for the same wall-clock, at the cost of needing explicit weights
to recover survey-wide totals.

Do not quote panel (a)'s cell values as precise until this is addressed. Panel (b), on all
74,812 events, does not have this problem.

---

## F2's long-`tE` bins rest on 124 events, and that is where the null result lives

The `sigma_piE` version of F2 (`DEVIATIONS.md` 24) reports that the joint fit adds nothing to
Roman's parallax precision for events longer than 300 days -- the ratio is flat at ~0.98
across the whole season gap. That is the negative answer to the plan's headline long-`tE`
prediction, so it deserves more support than it currently has.

It rests on **124 events**, spread over six `dt_edge` bins of 16-26 events each. The
supporting diagnostic (94% of those events have `piE` measured at better than 2 sigma by
Roman alone, median `sigma_piE` = 0.0098) is a strong and internally consistent explanation,
but the medians themselves are drawn from small samples, and the `> 300 d` bin is the one the
Kroupa mass function populates most sparsely -- long events need heavy lenses, which are
rare, and the remnant tail is exactly where the mass function is least certain.

**Why it matters scientifically:** "Rubin does not improve Roman's parallax measurement of
long events" is a claim about black-hole and neutron-star lens characterization, which is the
science case the long-`tE` regime exists for in this forecast. A null result at n = 124 is
suggestive; it is not yet quotable in a paper.

**Why it is deferred:** the fix is not a change to the analysis, it is more events in that
corner of parameter space, which is Step E1's stratified sampling -- the same fix the F3
footprint item needs. Doing it now would mean a second production run before the sampling
question is settled.

**What the fix involves:** stratify the draws in `tE` (Step E1) so the long-`tE` bins are
populated to comparable depth as the 30-100 d bin, and carry explicit weights to recover
survey-wide totals. Re-run F2 with `--param piE` on the stratified table and check whether
the flat ~0.98 line survives.

## F1/F2/F3 still load the whole 5.57M-row table and will be OOM-killed on a small machine (found during Step F4)

`analysis/romanlib.load_events()` read the entire per-event table into one DataFrame. For
`test5.dat` that is 5,571,168 rows x 90 float64 columns -- about 4 GB resident before
pandas' parse buffers are counted. On a machine with ~7 GB of RAM the kernel kills the
process.

**The failure mode is what makes this dangerous, not the failure itself.** The kill is
silent: the shell reports **exit status 0** with an empty stdout. It looks exactly like a
script that ran and chose to print nothing, and it produces no figure, no CSV and no error.
It was mistaken for a no-op for several minutes during Step F4 before the memory arithmetic
was checked.

**What was done:** `load_events()` gained optional `keep=` and `chunksize=` arguments that
filter per chunk while reading, holding the peak at one chunk plus the surviving rows. The
default path is byte-for-byte unchanged, so no existing caller's behaviour moved.
`analysis/f4_fisher_precision.py` uses it (`keep=lambda c: c["detJ"] == 1`), and the
production figure was verified identical to the reference computed the old way.

**What is deliberately NOT done:** `f1_results_table.py`, `f2_gap_filling.py` and
`f3_characterization_map.py` still call `load_events(path)` with no filter. They were run
successfully on 2026-08-30/31 when more memory happened to be free, so their outputs stand
-- but **they may not be reproducible on this machine as they are.** They are not being
changed here because each one selects a different subset (F2 gates on mission scope and the
Roman footprint, F3 on joint detection, F1 aggregates per field), so each needs its own
`keep` predicate written and its output re-verified against the existing CSV. That is a
step of its own, not a side-effect of F4.

**The fix, when it is taken:** give each of the three a `keep=` predicate that discards
undetected events during the read -- they are 98.66% of the table and every one of these
scripts throws them away immediately anyway -- then diff the regenerated CSV against the
committed one to prove nothing moved.

**Still not taken, and worked around instead during Step H6 (2026-09-12).** With ~2 GB free the
v3 table (6,075,044 rows) cannot be loaded at all, so the F-series was run against a streaming
`awk` extract of the rows with `detL | detR | detJ` set -- 82,888 rows, every column kept, the
`#` header preserved so `load_events()` still reads names from the file:

    head -1 test5.dat > det.dat
    awk 'NR>1 && ($38==1 || $39==1 || $40==1)' test5.dat >> det.dat   # detL detR detJ

F2, F3 and F4 are unaffected, because each discards undetected rows itself. **F1 is affected in
three fields** -- `N_events`, `N_neither` and `frac_gap_seen_by_rubin` become conditional on
detection -- so an F1 run on such an extract must not be used to quote a detection efficiency.
The workaround is a shell one-liner and does not close this item; the predicates still want
writing.

---

## Two Phase F panels pool events across sightlines and do not yet apply the area weight (raised during Step E1)

**What is wrong:** Step E1 makes the scan stratified — sightlines inside Roman's footprint can
be visited on a finer grid than those outside — so the sightlines no longer stand for equal
pieces of sky. Every event row now carries `w_area`, the deg² its sightline represents, and
`romanlib.area_weight()` returns it. **Two existing panels pool events across sightlines and
still count them one-for-one:**

- `analysis/f3_characterization_map.py` **panel (b)** — joint over Rubin-alone over all
  joint-detected events.
- `analysis/f4_fisher_precision.py` in its all-detections mode (`f4_fisher_all_kroupa.png`).

Everything else in Phase F is safe by construction and was checked one at a time: F1 tabulates
per Roman field (all in-footprint, all at the same fine step), F2 is restricted to the
footprint by construction (Deviation 23.1), F3 panel (a) and the F4 footprint panels are
in-footprint selections. Per-event ratios and per-sightline efficiencies are unaffected by
construction — which sightlines were visited is not an input to either.

**Why it matters scientifically:** unweighted, those two panels would describe a sky in which
Roman covers whatever fraction of the *sample* the stratification bought, instead of the ~2.6%
of the *sky* it actually covers. On the existing unstratified table the weight is a constant
and the panels are correct as they stand; the error appears only once a stratified run exists,
and it will not look like an error — F4's all-detections mass ratio would simply drift off
1.000, which is exactly the number Deviation 25.3 uses as a *bug detector* for the survey
partitioning. A weighting mistake would be read as a partitioning mistake.

**Why it is deferred:** applying the weight changes what those two panels plot — a weighted
CDF and a weighted 2D histogram rather than counts — and therefore changes committed figures
that currently stand as results. That is a science decision about how the all-sky panels should
be normalized, not a mechanical fix, and there is no stratified run to validate against yet.
Doing it now would mean editing figures to match a table that does not exist.

**The interim guard, already in place:** `romanlib.is_stratified()` reads the `stratified` flag
that the simulator now writes into `run_provenance.txt`, and `describe()` — which every figure
script already prints into its footer and its stdout — appends
`STRATIFIED(stride_roman=N; weight by w_area)` when it is set. A stratified run therefore
cannot be plotted by these scripts without saying so on the figure itself.

**What the fix involves:** pass `weights=R.area_weight(df, prov)` into the histogram and CDF
calls in those two code paths, re-run both against the unstratified table first (the weights
are constant there, so every number must come out unchanged — that is the regression), then
re-run against the stratified table.

---

## Step E1's `tE` stratification is not implemented; only the sky half of E1 is (raised during Step E1)

**What is wrong:** `JOINT_FIT_REFACTOR_PLAN.md` Step E1 asks for a fixed number of events **per
`tE` bin**, reweighted afterwards by the `tE` distribution to recover absolute yields, with bin
edges 1–5, 5–10, 10–20, 20–30, 30–60, 60–90, 100–200, 200–500, 500–1000 days. What was built is
the sky-position half: stratify the *scan* toward Roman's footprint, weight by `w_area`. The
`tE` axis is untouched, and `FunctE`'s 100 linear bins from 0 to 50 years remain what they were.

**Why it matters scientifically:** the plan wants the `tE > 200` d tail populated because that is
where the black-hole lens result lives, and a population-weighted draw produces few of them.

**Why it is deferred:** Deviation 26.1 argues the plan misdiagnoses which axis is starving those
bins. `tE` is not a drawn quantity — it falls out of the lens mass, the two distances and the
relative proper motion — so stratifying in it means acceptance sampling on a derived value,
carrying a second independent weight through every pooled statistic in Phase F, on top of the
sky weight just introduced. And the bins the plan wants filled are *footprint* bins: the
long-`tE` null rests on 124 events because it is a subset of the 1,950 in-footprint events, not
because long events are rare in the draw. The sky stratification multiplies that count by ~23
at `--stride-roman 2` with no acceptance sampling and no new weight to reason about.

**What the fix involves, if it is still wanted after a stratified run:** an acceptance test on
`l->tE` placed immediately after `func_lens()` and before the light-curve loop — which is where
it is nearly free, since the expensive per-draw work is the loop over ~50,000 epochs and
`FisherM`, all of it downstream of the point where `tE` is known. Accept with probability
`p[bin]`, carry weight `1/p[bin]`, and keep `nstE`/`ndtE` counting only accepted draws so the
per-bin efficiency stays unbiased (acceptance depends on the bin alone, so nothing *within* a
bin is distorted). The trap is the same one as above and worse: two weights now multiply, and
every absolute yield needs both.

---

## Roman's halo orbit around L2 is not modelled (raised 2026-09-04, planning Phase H)

**What is wrong:** Step H1 will place Roman at the mean L2 point. Roman actually flies a halo
orbit about L2 with an amplitude of order 10⁵–10⁶ km, so its true offset from Earth varies
through the mission.

**Why it matters scientifically:** it modulates an effect that is itself only ~10⁻³ in `u`, so
it is a fractional correction to a small term. It could matter for the narrow
high-magnification niche where satellite parallax does anything at all, because there the
sensitivity to `Δu` is highest.

**Why it is deferred:** the leading term has to exist before its correction is worth having,
and Step H3's experiment will say whether the effect is large enough for the correction to
change any conclusion. Modelling the halo orbit also needs an actual ephemeris or orbit
specification, which is not in the repository.

**What the fix involves:** replace the constant `L2_OFFSET_AU` displacement with a
time-dependent one, sourced from a published Roman orbit specification, and re-run H3's
comparison to see whether anything moves.

---

## No astrometric-shift analysis product exists (raised 2026-09-04, planning Phase H)

**What is wrong:** the astrometric microlensing signal — the centroid deflection
`δθ = θE · u / (u² + 2)` — is computed per epoch in `lightcurve()` (`s.def1c`/`s.def2c`, stored
in `l.soux`/`l.souy`) and feeds the astrometric Fisher matrix, but **no analysis script reads
it**. F4 plots `sigma_tetE`, which is the forecast *precision* on the angular Einstein radius,
not the shift itself, so nothing in the project says how large the astrometric signal is or how
often it exceeds the per-epoch astrometric precision.

**Why it matters scientifically:** whether the astrometric signal is *detectable* is a
different question from whether `θE` is *forecastable*, and the lens mass needs `θE` and `piE`
together. The maximum shift is `θE/√8` at `u = √2` — i.e. it peaks *outside* the photometric
peak, which interacts with Roman's seasonal gaps in a way nothing has yet looked at.

**Why it is deferred:** it is Step H5 of `PHASE_H_PLAN.md`. Its blocker, Step H4, is now
done — `errRomanA()` exists and Roman no longer borrows Rubin's astrometric error — so H5 is
unblocked and is simply not written yet.

**What the fix involves:** `analysis/h5_astrometric_shift.py` per `PHASE_H_PLAN.md` H5, then
re-run F4 on a post-H4 table and record how far the `tetE` and `Ml` panels moved. That
number is itself a result: it says how much the `errlsstA` placeholder was distorting them.

---

## RESOLVED 2026-09-05 — Rubin's astrometric error was a mission average used per epoch

**Status: closed by commit `d393dc3`.**

The suspicion recorded here was right and the cause was worse than "unsourced".
`files/sigmaA_LSST.txt` is a **mission-averaged** curve and `errlsstA()` fed it straight into
`l.erra[]`, which is a **per-epoch** error that `FisherM` divides by once per epoch and then
sums over ~2,300 epochs — applying the sqrt(N) averaging twice. Rubin's per-epoch astrometry
was **26.7x better than reality** and its astrometric Fisher information ~715x too large.

Two independent checks fix the factor. The table's bright-star floor is 0.3739576 mas, and
10.0/0.3739576 = 26.74 — exactly Ivezic et al. 2019's "assumed astrometric accuracy of 10 mas
per observation per coordinate" over sqrt(715), and ~715 visits is the right order for a
ten-year all-band LSST count. Scaled by the same 26.74 the faint end reads 132 mas at
r = 24.44 against an independent seeing-limited estimate FWHM/SNR ~ 700/5 ~ 140 mas. The
table's shape was right; only its normalisation was wrong, so the fix renormalises rather than
replaces, and does it in code so the delivered input file is untouched.

Roman is now the better astrometer, which is the physically obvious ordering and was not true
before. 10 mas is the conservative reading: Rubin commissioning (SITCOMTN-159) reports a 3-7
mas single-visit systematic, so the delivered floor may be better than assumed here.

**Every `tetE` and lens-mass number involving Rubin astrometry moves.** Nothing has been
re-run against it yet; the 2026-09-05 v2 production run is the first that will be.


## The astrometric shift is not diluted by blending (found in Step H5)

**What is wrong:** the simulator's modelled centroid, `s.pos1c`/`s.pos2c` in `lightcurve()`,
adds the lensing deflection `s.def1c`/`s.def2c` undiluted — it is the SOURCE's centroid. A real
measurement sees the centroid of every star in the aperture, so an unlensed blend of fraction
`(1 - fb)` drags the measured shift down by roughly `fb`. The per-epoch astrometric *error* is
evaluated at the blended magnitude, so the photon-noise half of blending is modelled; the
centroid-dilution half is not.

**Why it matters scientifically:** it makes every astrometric signal optimistic by a factor of
about the source-flux fraction, and the bulge is a crowded field where `fb` is often well below
1. `theta_E`, and hence every lens mass that uses it, inherits the optimism. **For Rubin this is
large** — the median `blend_r` is 0.15, so a real Rubin centroid would move ~7x less than the
simulator's. For Roman it is currently moot, but only because `blend_F146` is pinned at exactly
1.0 by the entry below; fix that first and this one grows teeth.

**Why it is deferred:** it is a change to the forward model, not to an analysis script, and it
would move every astrometric number in the project at the same time as Step H4 already has.
Doing both in one step would make it impossible to attribute the change to either.

**What the fix involves:** multiply the deflection term by the per-telescope source-flux
fraction where `pos1c`/`pos2c` are formed, and confirm against the fixture that `sigma_tetE`
degrades by roughly `1/fb`. Meanwhile `analysis/h5_astrometric_shift.py` panel (a) plots the
undiluted shift and the `x fb` version side by side, so the size of the simplification is
visible on the figure rather than hidden.

---

## RESOLVED 2026-09-05 — Roman was unblended because a Poisson count was drawn as a Gaussian

**Status: closed by commit `58d2863`.** The diagnosis recorded here was wrong; checking it is
what found the real cause.

**This entry guessed the input population was at fault** — that `Nstart` counts only stars
above some completeness limit, so the faint unresolved population that actually blends a Roman
pixel was missing from the density. **It is not.** `Nstart` is built from the mass density
divided by the mean mass of each population, so it is complete down the whole IMF and stops at
no survey's detection limit. The geometry was right too. The statistics were wrong.

The mean number of field stars in a seeing disc is a surface density times a disc area: ~12.7
for Rubin's 0.993" r-band disc at a bulge `Nstart` of 2.13e8 per deg^2, and **~0.14 for
Roman's 0.105" F146 disc**. The old code added `RandN(sqrt(mean), 2.0)` and clamped up to 1. At
12.7 that is a passable Gaussian approximation to a Poisson count. At 0.14 it is qualitatively
wrong: the truncated Gaussian can add at most 2*sqrt(0.14) = 0.75, and the clamp then rounds
**every** draw to exactly 1. A mean of 0.14 neighbours was rendered as "no neighbours, always".

The fix is a Poisson draw, `1 + Poisson(mean)` — the `1 +` because we are looking AT a source,
so the disc is conditioned to contain it, and Slivnyak's theorem says the rest of a Poisson
field is unchanged by that conditioning. `max(1, Poisson)` would absorb the first neighbour
into the source and leave Roman unblended 99% of the time, the same bug one step on.

**Measured, same stub patch, before and after:** `blend_F146` fraction below 1 goes from
**0.000 to 0.146**, against the predicted 1 - exp(-0.142) = 0.133; median stays 1.0, because
Roman really is mostly unblended at 0.105" — that part of the old answer was right. `fb1` is
off its boundary for the same 14.6% of events instead of pinned for all of them. Rubin moves as
predicted and only slightly, `nsbl_r` 14.3 -> 15.0, the same Palm-conditioning correction.

**Consequences to re-examine once the v2 run lands:** Roman's detection efficiency was
`testR <= blend[6]` with `blend[6] == 1`, so Roman accepted every detectable event while Rubin
accepted ~15%. That asymmetry fed every Roman-against-Rubin yield comparison, including the
gap-filling claim. It should now be smaller, and the size of the change is a result in itself.


## The provenance stamp goes stale silently, and discipline has now failed twice

**Status:** open. Found on 2026-09-05, while checking the freshly launched production run.

**What happens.** `Makefile` captures the git description into `GIT_COMMIT` at the moment
`make` runs, and compiles it into the binary so every run's provenance block can name the
source that produced it. But `make`'s dependency graph is over *files*, and a commit changes no
file. So after `git commit`, a subsequent `make` says "Nothing to be done for 'all'" and the
binary silently keeps the stamp it was built with — which names the previous commit, with a
`-dirty` suffix, because the tree was dirty when it was built.

The Makefile already knows this and says so in a comment: *"the value is captured when make
runs, so `make clean && make` is what refreshes it after a commit -- an incremental build keeps
the stamp its objects were built with."* The mitigation is therefore a rule a human has to
remember.

**It has now failed twice.** Commit `585432b` ("Rebuild so the next run is labelled with the
commit that produced it") exists solely because the stamp read `d6dd293-dirty` when the source
was `7a1b591`. On 2026-09-05 the identical thing happened again: the run launched with a binary
stamping `6c97375-dirty` when the source was `e8f4135`. A rule that has been broken every time
it has been tested is not a working control.

**Why it matters more than it looks.** The stamp is the only link between a multi-GB table and
the code that produced it. Every figure, every quoted number and every caveat in `PROGRESS.md`
section 4 is ultimately justified by "this run contains change X", and the stamp is the evidence.
A stamp naming a commit that predates the change it is being used to vouch for inverts that: a
future session would look up `6c97375`, find no H2 columns in it, and conclude either that the
table is corrupt or that the columns came from somewhere unrecorded.

**Why it was not fixed on the spot.** The run it was found on was already in flight and the
finding is a build-system change, not a physics change; folding it into a launched run would
break the rule that a step does one thing. The specific run was rescued instead, by proof
rather than by restart — see `PROGRESS.md` section 5c.

**What the fix involves.** Make the stamp a file so `make` can reason about it. Generate a
`git_commit.h` from a phony rule that rewrites it only when the description actually changes,
and have the sources include it; then a commit invalidates exactly the objects that embed it and
nothing else, and no one has to remember anything:

```make
.PHONY: force
git_commit.h: force
	@echo '#define GIT_COMMIT "$(GIT_COMMIT)"' > $@.tmp
	@cmp -s $@.tmp $@ || mv $@.tmp $@   # only touch it when it really changed
	@rm -f $@.tmp
```

Two smaller things belong with it. The binary should be able to print its own stamp
(`./roman --version`) so it can be checked without launching a run or reading a provenance
block. And `--dry-run` should print the stamp it would write, since `--dry-run` is already the
documented way to inspect a run's cost before committing hours to it, and the stamp is part of
what you would want to inspect.

---

## RESOLVED 2026-09-05 — DET_ANOMALY was never 0; the run total was never counted

**Status: closed by commit `a5bb600`.** Kept because the wrong hypothesis in it is instructive.

**What was suspected.** The 2026-09-05 partial run logged 433 DET_ANOMALY events in 23 of its
24 Roman-covered sightlines where the 2026-08-30 run reported 0 in 5.57M events, and this entry
named **H1 (satellite parallax)** as the leading suspect, on the reasoning that it is the only
one of E1a/H1/H2/H4 that touches the photometric model.

**That was wrong, and the experiment this entry specified is what showed it.** A paired A/B on
`--no-satellite-parallax` — same seed, same stub patch inside the footprint, `--events` equal to
`--maxdraws` and `--lenses` set out of reach so both variants execute exactly 150 draws per
sightline whatever they detect, hence identical events — came out **indistinguishable**:

| | none | Rubin+joint | Roman+joint | both+joint | DET_ANOMALY |
|---|---|---|---|---|---|
| satellite parallax on | 1197 | 36 | 81 | 35 | **15** |
| satellite parallax off | 1197 | 36 | 81 | 35 | **15** |

**The actual cause.** `nDetClass[DET_ANOMALY]`, the per-sightline counter, was incremented.
`NDetClassTot[DET_ANOMALY]`, the run total, was **not** — even though it is read at the end of
the run to decide whether to print the note explaining what DET_ANOMALY means. Every run since
the counter was introduced reported exactly zero anomalies in its summary and the explanatory
note never fired once. **The 2026-08-30 run's "0 ANOMALY" is that bug, not a measurement.**
The anomalies were in its log all along; nobody read that log, because it was unavailable.

**So nothing regressed.** The rate is what Deviation C-D.2 predicted: thresholding `dchi`
against `2 * ndw` is a bar on the MEAN per-epoch improvement, and where Roman contributes
50,401 epochs against Rubin's ~2,300, Rubin's epochs lift the joint bar by ~4,728 while adding
little signal. Any event clearing Roman's own bar by less than that fails the joint one — about
1% of footprint draws. `detJ` is still forced monotone, so no table was ever wrong.

**Still open, and now with evidence behind it:** the `2 * ndw` threshold form itself, which
Deviation C-D.2 deliberately left alone. A chi-squared statistic should be thresholded on its
total, not its mean. With Roman's epoch count 22x Rubin's, that asymmetry is no longer
hypothetical, and the anomaly rate is the measurement of how hard it now bites.

---

## `--dry-run` destroys the previous run's output, and it has already cost one table

**Status: open. Found 2026-09-05 the hard way — it destroyed the v1 partial production table.**

**What happens.** `Bulge_LSST.cpp` opens the event table with `std::ofstream head(testf);` —
truncating mode, the default — and writes the column header, at around line 632. The
`--dry-run` early exit is at around line 916. So **every `--dry-run` truncates `test5.dat` and
zeroes the append-mode files in `files/MONTLMC/files/` before it decides not to simulate
anything.** A flag whose entire purpose is "build the scan and tell me the cost without
running" is destructive to the previous run's output.

**What it cost.** At 18:07 on 2026-09-05, `./roman --dry-run --stub --stride 2` was run from
the worktree to check the freshly rebuilt provenance stamp. The worktree's `test5.dat` was at
that moment still symlinked to the v1 partial run directory, so the dry run truncated it from
773 MB and 1,643,024 rows to a 620-byte header, and zeroed its six MONTLMC outputs.

The scientific loss is small — `run.log` (97 MB) and `README.md` survived, so the cost model,
the detection-class totals and the 433 DET_ANOMALY events were all already extracted into
`PROGRESS.md` and `DEVIATIONS.md`, and v1 was superseded by v2 in any case. But `PROGRESS.md`
§5c said that table was preserved, and it was not.

**Why this is dangerous rather than annoying.** The per-event write is
`filg_in.open(testf, std::ios::app)` **inside the event loop**, so the path is re-resolved on
every event. That means a live run cannot be protected by repointing the symlink out of harm's
way mid-run — the running process would simply follow it. While a production run is in flight
through the worktree symlink, **any** `./roman` invocation from the worktree, `--dry-run` and
`--help`-adjacent flag checks included, destroys that run's table.

**What the fix involves.** Move the header write after the `--dry-run` exit, which is where it
belongs — a dry run should touch no output file at all. Two things belong with it:

1. **Refuse to truncate a non-empty table without being told to.** The run already knows the
   file name; if `test5.dat` exists and is larger than its header, require an explicit
   `--overwrite` (or write to `test5.dat.N`). Losing hours of simulation to a mistyped command
   should not be possible.
2. **Say what is about to be destroyed.** The startup banner already prints the configuration;
   printing "truncating an existing 773 MB table" would have made this visible immediately.

**Until then, the working rule:** never invoke `./roman` from the worktree while a production
run is writing through the worktree symlinks, not even `--dry-run`. Use a copy of the binary in
an isolated directory (the pattern in `setup_variant.sh`) for any check made while a run is in
flight.

---

## `dchiP` and `dchiA` are still absolute values while `dchiL` is signed

**Status: open. Created deliberately by Step H7 on 2026-09-06, not inherited.**

Step H7 made the lensing detection statistic signed:

```cpp
dchiL   = chi3  - chi1;             // signed: chi3 is the flat-baseline residual
dchiP   = std::fabs(chi2  - chi1);  // still an absolute value
dchiA   = std::fabs(chi2a - chi1a); // still an absolute value
```

**Why `dchiL` had to change.** Only `chi3 - chi1 > 0` — the lensing model fitting *better* than
a flat baseline — is evidence of lensing. Under `fabs`, a lensing model fitting far *worse*
than the baseline would have cleared the detection bar. More importantly it breaks the
monotonicity H7 exists to establish: an instrument whose epochs all fall outside the event
contributes a small negative difference from noise, and `|dchi_L + dchi_R|` can then fall below
`max(|dchi_L|, |dchi_R|)`, which is exactly the `DET_ANOMALY` condition. Signed, the sum is
exact and the joint test cannot contradict a single-survey detection.

**Why the other two were left alone.** `dchiP` (parallax) and `dchiA` (astrometric deflection)
are not detection tests — nothing thresholds them. They are reported as the *size* of a
perturbation, for which an absolute value is defensible. Changing them would alter two more
output columns for no acceptance criterion, and Step H7's remit is the detection decision.

**Why it is still an open item.** Three columns that look like the same kind of quantity now
have two different sign conventions, which is exactly the sort of thing that produces a wrong
plot two months from now. Either give `dchiP`/`dchiA` the same signed convention — noting that
"the parallax model fits worse than the no-parallax model" is a meaningful negative that
`fabs` currently hides — or rename the two absolute ones so the asymmetry is visible at the
call site. Do it as its own step with its own before/after, not as a drive-by.

---

## Step H7's cost regression: RESOLVED by profiling, and its size was overstated

**Status: RESOLVED 2026-09-07 by `perf` profiling. The mechanism is identified and the
magnitude was overstated by roughly a factor of 3. The original entry is kept verbatim below,
because the two hypotheses it refutes are still correctly refuted; only its conclusion is
withdrawn.**

### The answer

**H7 costs more because it detects more, and for no other reason.** `FisherM` -- and, inside
it, `lightcurve` and the GSL matrix work -- runs only under `if (detL or detR or detJ)`.
Lowering the bar to a fixed 500 roughly doubles how often that branch is taken, and the cost
follows the branch. Nothing got slower per call.

Profiled with `perf record -F 999 -e cycles:u` on scan index 716 (a footprint sightline,
`ndd (Roman) = 50401`) under both binaries -- v2's actual `a5bb600` build and v3's `f959c8c`
build, verified against each run's `run_provenance.txt` -- with a fixed draw budget
(`--maxdraws 100`, `--events`/`--lenses` unreachable) so both simulated the same population.
The pairing held: 104 events simulated against 101.

| | pre-H7 | post-H7 | ratio |
|---|---|---|---|
| characterised detections (`LpLMC5.dat` rows) | 7 | 18 | **2.57x** |
| detection-gated path (`FisherM` + `lightcurve` + `gsl_*`) | 21.3 s | 48.6 s | **2.28x** |
| everything else | 228.3 s | 254.0 s | **1.11x** |

The gated path grows *slightly slower* than the detection count, so the cost **per Fisher call
is unchanged or marginally lower**. The "roughly 3x per characterised event" conclusion in the
entry below is withdrawn. Everything outside that branch is unchanged to within the
measurement's error.

Note this is not a regression to regret: the extra cost buys the detections H7 exists to
recover. Per *detected* event, H7 is slightly cheaper than what it replaced.

### The magnitude: ~1.4-2.0x per footprint sightline, not 2.5-6.6x

Two independent estimates, both far below the old claim:

- **From the profile.** Startup is a large one-time cost and has to be removed first (see the
  next item). After removing it the sightline itself costs **80.0 s -> 137.1 s, a factor
  1.72**; composing the two rows of the table above instead gives 1.35. The spread between
  those is the honest precision of a one-sightline profile sitting on a ~170 s startup.
- **From the corrected A/B.** Same fixed budget, 900 s each, startup removed: pre-H7 does
  6 footprint sightlines, post-H7 does 4, giving **~2.0x**.

**Why the old number was ~3x too large.** The claim "the pre-H7 binary completed 18 footprint
sightlines to the post-H7 binary's 3" counted `NEW STEP` lines, which count *every* sightline
entered, not footprint ones. Split by stratum, those same two logs say:

| in 900 s, from index 716 | footprint | outside |
|---|---|---|
| pre-H7 | **6** | 13 |
| post-H7 | **4** | 0 |

Six against four, not eighteen against three. The pre-H7 side spent most of its run on cheap
outside sightlines and was credited with footprint throughput it never had. This is the third
time in this investigation that a per-stratum cost was quoted from a number that mixed strata.
**Any statement about sightline cost must name the stratum and be counted with
`ndd (Roman): 50401` versus `ndd (Roman): 0`, never with `NEW STEP`.**

### Method note: why the cycle counts could not simply be subtracted

`perf` sampled in frequency mode and was throttled to 808-924 Hz rather than the requested
999; its "Event count (approx.)" is the sum of adapted sample periods. Across runs under
different contention those totals are not comparable -- subtracting them made a footprint
sightline look like 4% of its own run, which is impossible. Within a single profile the
per-symbol share of samples is reliable, and so is the run's wall span; every number above
uses only those two.

Startup was separated using float parsing as a tracer: `std::num_get::_M_extract_float` runs
only at startup, so its share of a paired run divided by its share of a pure-startup run
(`--start-index 5000`, which enters no sightline) is the fraction of that run that was
startup. The check that this works: it recovers **171.0 s and 167.6 s** from two independent
profiles for a quantity that must be identical.

Scripts: `perf_h7.sh`, `perf_startup.sh`, `perf_cmp.py`, `h7_final.py`, `h7_attrib.py`.

## Startup costs ~170 s, ~40% of it parsing ~1000 extinction files, and every chunk pays it

**Status: open, measured 2026-09-07. Not a correctness issue.**

Running with `--start-index 5000`, so that no sightline is ever entered, still takes **116 s of
wall time** and, under load, ~170 s of CPU. Its profile is dominated by text parsing:
`_M_extract_float` alone is **27.7%** and the whole iostream float-parsing group is ~43%. That
is `readBayestar` reading the ~1000 `files/ext/bayestar_*.txt` files with `>>` on an
`ifstream`.

Across a full 1829-sightline run this is negligible, ~0.5% of the total. It matters only
because the runs are now **chunked**: every resume pays it again, so a 10-chunk run spends
~30 minutes re-reading extinction files it has already read.

Worth doing only if chunking gets much finer, and the fix is routine (`std::from_chars`, or
cache the parsed grid to one binary file). Filed rather than fixed because it is not this
step's business.

---

<!-- Superseded 2026-09-07 by the entry above. Kept verbatim: its two refuted
     hypotheses remain correct and useful; only its concluding magnitude is withdrawn. -->

## Step H7 made footprint sightlines ~6x more expensive, and the mechanism is unknown

**Status: SUPERSEDED 2026-09-07 — see the resolved entry above. Left as written for the
record. Its two refuted hypotheses stand; its concluding magnitude (2.5-6.6x, and "roughly 3x
per Fisher call") is withdrawn, having rested on a footprint-sightline count that mixed
strata.**

**What was predicted.** H7 lowers the detection bar from `2*ndw` (100,802 inside the footprint)
to a fixed 500, so detections arrive about twice as fast. The per-sightline loop stops on its
`--events`/`--lenses` targets, so reaching those targets in fewer draws should make every
sightline *cheaper*. That is what was written in `PROGRESS.md` §5e and in Deviations 33.

**What happened.** Outside the footprint, exactly that: **14.3 s -> 6.8 s** per sightline.
Inside it, the opposite and larger: **~346 s** per sightline, against somewhere between 52 s
and 136 s for v2. The footprint term dominates the total, so the run got *slower*: projected
~8.8-13 h for v2 against **~16 h** for v3.

**The v3 number is solid and the v2 number is not.** Both v3 runs agree independently
(299 and 329 s/sightline at n=11, 339 and 346 at n=41), and v3's outside cost rests on 618
samples. But v2's per-stratum split came from a two-run solve that was already flagged as
weakly determined: its footprint cost is 52 s if its outside sightlines cost 14.3 s and 136 s
if they cost 6 s, and v2's outside cost was never measured cleanly — the same mistake this file
records elsewhere, made a second time. **So the regression is a factor of 2.5-6.6, direction
certain, magnitude not.** Quoting "six times" without that range overstates what was measured.

**Why it is strange.** The comparison is at *identical* sky positions — the v2 and v3 scans are
the same deterministic grid, and the first twelve footprint sightlines were verified to match
on `nri`, lon and lat. At those same positions v3 does **less** work:

| per footprint sightline | v2 (pre-H7) | v3 (post-H7) |
|---|---|---|
| scored events (lightcurves built) | 333-917 | 300-435 |
| detections (= Fisher calls) | 56-87 | 50-81 |

Fewer lightcurves built, no more Fisher matrices, the same 50,401 Roman epochs and the same
`ndw` — and six times the wall clock. Nothing in the H7 diff obviously accounts for that.

**The first hypothesis was tested and is REFUTED (2026-09-06).** It proposed that a bar of 500
admits marginal events whose Fisher matrices are near-singular, so that `invert_matrix` does
expensive work that is then thrown away. Characterisation rates on footprint events say
otherwise:

| detections characterised (`okA == 1`) | v3 post-H7 | v2 pre-H7 |
|---|---|---|
| joint | 2405/2419 = **99.4%** | 3276/3301 = 99.2% |
| Rubin | 1308/1320 = 99.1% | 1928/1954 = 98.7% |
| Roman | 1693/1709 = 99.1% | 2201/2225 = 98.9% |

Marginal detections invert *at least as reliably* as the overwhelming ones did. Nothing is being
computed and discarded. Whatever costs the time, it is not failed inversions.

**What the same measurement makes stranger.** Per footprint sightline, v3 builds **307**
lightcurves against v2's **510**, and makes **59** Fisher calls against v2's **51**. Fewer
lightcurves, ~16% more Fisher calls, and several times the wall clock. For the arithmetic to
work, the cost *per Fisher call* would have to have risen roughly fivefold — but `FisherM`
evaluates a fixed finite-difference stencil whose cost goes as `ndw * Nx`, and neither changed.
**A second hypothesis was tested and is also REFUTED (2026-09-06).** It proposed that the old
`2*ndw` bar systematically rejected high-epoch events — an event with more epochs faced a
proportionally higher bar — so that H7 now admits exactly the events whose Fisher matrices cost
the most, since `FisherM` scales with `ndw`. The table says no:

| mean `ndw` over joint detections, footprint | v3 post-H7 | v2 pre-H7 |
|---|---|---|
| all footprint events | 51,299 | 51,217 |
| joint detected | **51,086** | **50,914** |
| detected/all ratio | 0.996 | 0.994 |

A 0.3% difference. Inside the footprint `ndw` is dominated by Roman's 50,401 epochs and is
near-constant across events, so `2*ndw` was effectively a *constant* bar as well — just a very
high one. Its scaling with epoch count did its damage *between partitions* (the joint bar
sitting above Roman's own), which is what H7 fixed, not between individual events.

**What is established, and what is not.** The controlled A/B settles that the regression is
real and not an artefact of comparing two production runs: with a *fixed* draw budget, so both
sides build identical lightcurves, the pre-H7 binary completed 18 footprint sightlines to the
post-H7 binary's 3 in the same 900 s, while making only 1.8x as many Fisher calls. That is
roughly 3x the cost per Fisher call, at identical `ndw`, with `FisherM` itself unchanged.

**No further hypotheses should be invented here.** Two have been proposed and both refuted by
measurement; a third would be guessing. Settle it by profiling — `perf record` on one footprint
sightline for each binary at the same `--start-index` and the same fixed draw budget, which is
about 20 minutes and will name the function. Until that is done, the honest statement is that
H7 costs roughly 3x per characterised event for reasons not yet identified.

**How to settle it cheaply.** Run one footprint sightline under `perf record` (or simply time
`FisherM` and `invert_matrix` with a counter) against the pre-H7 binary at the same
`--start-index`, with a fixed draw budget so both sides see identical events — the pattern in
`h7_ab.sh`. About 20 minutes, and it competes with the production runs for cores, so it belongs
between chunks rather than during one.

**Why it might matter beyond runtime.** If marginal detections are producing ill-conditioned
Fisher matrices, the question is not only what they cost but whether their `sigma` values are
trustworthy. `okA`/`okB` gate on successful inversion, so nothing invalid should reach the
table — but the *fraction* of detections that fail characterisation is now a number worth
reporting, because H7 changed which events are offered to the Fisher step at all. Check
`sigmaas`/`okA` rates on the v3 table against v2 before quoting any yield.



## Step H3's paired Fisher comparison fails its own physical checks -- RESOLVED

**Status: RESOLVED 2026-09-11. Cause found: `FisherM` differenced every derivative against the
model value cached at light-curve generation time, which belongs to the observer the run
actually used. Re-evaluating with Roman moved to Earth made that reference inconsistent with
the perturbed model, injecting a constant `dm_sat/Delta` into every derivative that the
symmetric stencil happened to cancel for `u0`/`tE`/`piE`/`xi`/`t0` and did not cancel for the
cross-telescope rows, `fb`'s outer bins, or anything astrometric. Fixed by recomputing the
reference under the observer in force; `./fishertest` byte-identical across the change, so no
production result moves. Full account in DEVIATIONS.md 36.**

**Neither of the two candidates guessed below was the cause.** Recording that, because the
diagnostic that would have separated them -- the condition numbers -- would not have found
this, whereas asking "what is each derivative actually differenced against?" did.

**The conditioning guess is now refuted by data, not merely unused.** `h3_pair.dat` records both
matrices' condition numbers for both observer positions since `a3c323d`. Measured on the fixed
run, over Roman-covered events:

| matrix | satellite median | no-satellite median | ratio no-sat / sat |
|---|---|---|---|
| photometric `condA` | 1.094e+06 | 1.214e+06 | **1.0094** |
| astrometric `condB` | 4.545 | 4.548 | **1.0001** |

The two observer geometries are equally conditioned to within 1%. There is no near-degeneracy
at `satScale = 0`, so that hypothesis is closed. The other guess -- that the per-epoch
photometric weights were not recomputed on the flip -- was also wrong, and in a way worth
keeping: those weights are deliberately held fixed, which is what makes the comparison
controlled. Changing them would have confounded the geometry change with a noise-model change.

The original entry is kept verbatim below. SUPERSEDED.

---

**Status: open, measured 2026-09-11 on the full paired run. Blocking the H3 result; not
blocking anything else.**

`--pair-satellite` characterises each detected event twice, at L2 and at Earth, in one run. The
plumbing is verified -- events with no Roman epochs near the peak return a ratio of exactly
1.000000 on n = 29,171 -- but for the 2,264 Roman-covered events the comparison fails two
checks that follow from the physics rather than from statistics:

- the improvement **shrinks** as `du_sat` grows (log-log correlation +0.227), where the whole
  mechanism requires it to grow;
- `sigma_tE` is worse for 85.5% of events, median ratio 4.88, which cannot be an information
  statement about the same event at the same epochs.

The diagnostic that points at the cause: it is the `satScale = 0` forecast that swings wildly
between the improving and worsening groups (`sigma_piE` 0.904 vs 0.165), not the `satScale = 1`
one.

**Two untested candidates, recorded rather than asserted.** (1) The per-epoch photometric
weights entering `FisherM` are computed from the `satScale = 1` magnitudes and are not
recomputed when the offset is flipped. (2) With `satScale = 0` the two observatories see model
curves identical in shape, which may leave a near-degeneracy that clears the `okA`
condition-number gate while leaving marginalised errors unstable. Recording the condition
numbers of both matrices in `h3_pair.dat` would separate these cheaply and was not done.

Until this is settled H3 reports only `du_sat` (median 0.0022, max 0.039 theta_E) and states
that the precision consequence is unmeasured. `analysis/h3_satellite_parallax.py` enforces this
in code: it withholds the ratios when the checks fail.


## ADDRESSED 2026-10-01 BY A BRACKET (Deviation 71) -- Roman's 1.1 mas astrometric floor is treated as independent per exposure, and it may not be

**STATUS 2026-10-01.** The literature (Sanderson+2019, Lam+2026, McKinnon & van der Marel 2026,
Kaczmarek+2026; read in full, Deviation 71) treats the floor as white, with the GBTDS's sub-pixel
dithers as the physical justification, and quantifies no correlated part. So every event now carries
three forecasts: W (white; the main columns), N (offset per Roman roll + 0.3 mas per coordinate
shared within each day) and P (offset per Roman season + the full 1.1 mas shared within each day);
columns `sigtetE_{N,P}{J,R}`, `relMl_{N,P}{J,R}`, `okB_{N,P}{J,R}`; `romanlib.sigma(..., noise=)`.
**Two corrections made on the way, which move every theta_E number:** the source's reference
position is now free (it was not: theta_E was partly measured from the absolute position), and the
per-coordinate error is no longer doubled in variance. **The text below is kept for the record but
is partly wrong:** "in the fully-correlated limit ... theta_E would not be measurable at all" -- with
the reference position free, an error constant over the mission is absorbed and costs nothing; only
correlation on timescales shorter than the event matters, which is what N and P model. Remove this
entry once the N/P bracket has been measured on production runs and reported. The remaining open
question is the TRUE split, which only Roman data (or the Roman project's astrometry simulations)
can give.

**Status: open, identified 2026-09-11 while producing the H5 astrometric results. This is the
dominant caveat on every theta_E and lens-mass number in the project.**

Step H4 gives Roman a per-exposure astrometric error with a floor of 1.1 mas, and documents that
floor as **1% of the 110 mas pixel -- a centroiding systematic, not photon noise**, which is why
it does not improve for brighter stars. That reasoning is about brightness. It leaves the
separate question untouched: does the floor average down over *exposures*?

`FisherM` assumes it does, and completely. Each epoch contributes
`(dtheta/dparam)^2 / erra[i]^2` as an independent term, so N epochs buy a factor `sqrt(N)`. The
numbers that follow from that assumption are large:

- ~50,000 Roman exposures per detected event (median `ndwR` = 49,977; a GBTDS field gets 50,401
  12.1-minute exposures and `--stride-roman` strides sightlines, not epochs);
- the per-exposure precision for these sources is **6.69 mas** at the median, not the 1.1 mas
  bright-source floor -- `errRomanA` at each source's own F146 magnitude, median F146 ~ 22.98;
- so an independent error averages to **0.0300 mas**;
- against a median peak centroid shift of **0.1106 mas**, i.e. a signal at **3.68 sigma** of the
  averaged precision, even though only **0.11%** of events have a shift exceeding a SINGLE
  exposure's precision.

**The whole astrometric result rests on that factor of ~224 in averaging.** If any component of the 1.1 mas is
correlated between exposures -- and a centroiding systematic is exactly the kind of error that
usually is, through the PSF model, the distortion solution, or the reference frame -- then the
effective precision is worse, approaching the full 6.69 mas per-exposure value in the
fully-correlated limit, and every `sigma_tetE` and `relMl` in the project is optimistic. In that
limit the median signal sits nearly two orders of magnitude *below* the noise and theta_E would
not be measurable at all.

So the honest statement of the H5 result is conditional: 33.2% of Roman-observed detections
measure theta_E to better than 10% **if the per-exposure astrometric error is independent
between exposures**. That conditional is not currently stated in the whitepaper.

**What would settle it.** The error model would need a two-component form -- a white part that
averages and a correlated floor that does not,
`sigma_eff^2 = sigma_white^2/N + sigma_corr^2` -- with the split taken from the Roman
astrometry literature rather than assumed. Sanderson et al. 2019's claim of ~0.1 mas from
stacking ~100 exposures is itself evidence that *some* averaging is real (1.1/sqrt(100) = 0.11,
which matches their number almost exactly, implying the floor is close to fully white over 100
exposures); whether that continues for another 500x in N is the open question, and nothing in
the two sources cited by H4 addresses it.

Recording rather than fixing, because inventing a split would be worse than naming the
assumption.

---

## Step G2 as the plan writes it has no published number to compare against; the replacement is designed but deferred

**Status: open, raised 2026-09-15. The user chose to skip G2 for now and do E2 first.** Full
reasoning: `DEVIATIONS.md` entry 39 and `PROGRESS.md` §5i.

**What is wrong.** `JOINT_FIT_REFACTOR_PLAN.md` Step G2 asks for the Rubin-alone characterised
fraction (`tE > 2 sigma_tE` and `piE > 2 sigma_piE`) at Abrams et al. 2025's bulge field,
l = 0.33°, b = 2.82°, compared against "their published value". Abrams et al.
(arXiv:2309.15310) **publish no absolute value for that quantity** -- §3.7 shows only
OpSim-to-OpSim ratio maps (Figs. 11-14) and an unlabelled histogram (Fig. 10). Their only
absolute efficiencies are Table 6's `sigma_tE/tE < 0.1` fractions, which are sky-wide,
parallax-free, single-mean-star and per-band-flux, on OpSims up to `baseline_v3.0`.

**Why it matters scientifically.** G2 is the one external check on the Rubin branch -- the answer
to "why trust a new simulator". The advisor's paper (arXiv:2608.16448) cannot substitute: it runs
the same parent code, so a shared defect would pass silently.

**Why deferred.** User decision, 2026-09-15. It also needs a new dependency and two large
downloads (below).

**What the fix involves -- Option A, a paired comparison on identical inputs:**

1. **Anchor.** Install `rubin_sim` in `.roman/` (not installed as of 2026-09-15; needs network,
   plus its support data, possibly several GB -- check disk and the ~2 GB free RAM first).
   Download `baseline_v3.0_10yrs.db` (~800 MB) and reproduce one Table 6 row with
   `rubin_sim`'s `MicrolensingMetric`: Fisher metric 0.05 / 0.12 / 0.22 / 0.80 in `tE` bins
   10-20 / 20-30 / 30-60 / 200-500 d. This proves their metric is being run as they ran it
   before anything is compared against it.
2. **Cadence.** Query `Baseline/baseline_v5.1.0_10yrs.db` (on disk) for a 3.5° square around
   RA = 263.89°, Dec = -27.16°. **Do not reuse `Baseline/BulgeBaseline.dat`**: it spans
   b = -3.59 to +1.10 and has zero visits within 1.75° of that field.
3. **One event list, their sampling.** `tE` uniform within each bin, `u0` uniform on [0, 1],
   `t0` uniform over the survey, source magnitudes u 25.2, g 25.0, r 24.5, i 23.4, z 22.8,
   y 22.5, blend fraction 0.5 (Abrams et al. §2.3).
4. **Both codes on every event.** `rubin_sim`'s Fisher metric, and this pipeline's `FisherM` in
   a new standalone driver in the style of `tests/fisher_fixture.cpp` (production code
   untouched), with `piE` and `xi` held fixed by subset inversion of `inputA[SRUBIN]` so both
   fit the same geometric parameters (method and its `Era` cross-check: Deviation 39).
5. **Compare per event**, not only per bin: the `sigma_tE` ratio distribution, then the
   `sigma_tE/tE < 0.1` fraction per bin.
6. **Fix the matching rules before running**, so agreement cannot be tuned into existence. Two
   differences are expected and are not bugs; measure each separately if the codes disagree:
   (a) Abrams et al. fit a source/blend flux pair **per band**, this pipeline one r-band pair
   (`RUBIN_REF_BANDS = {2}`, unfinished Step C2, Deviation 3), so theirs should be the larger
   sigma; (b) `calc_mag_error_m5` in `rubin_sim` against `errlsstM`, same Ivezic et al. 2019
   family, constants to be compared.

**Option B, cheap and complementary:** ask Martin Makler (co-author of both papers) for the
numbers behind Abrams et al.'s bulge-field parallax characterisation (Fig. 10). It is the only
route to testing the parallax quantity the plan actually named.

**Rejected, Option C:** comparing directly against Table 6. A sky-wide v3.0 number against a
single-field v5.1 number cannot distinguish a bug from a setup difference.

---

## Pooled per-event statistics give every sightline the same number of events, whatever its event rate

**Status: open, found 2026-09-15 while preparing Step E2. DERIVED AND MEASURED 2026-09-15 --
Deviation 41. The framing below is superseded:** the per-sightline terms move pooled fractions
by under a point; the dominant term is *within* each sightline. Draws are not rate-weighted in
lens mass or velocity, so the weight is
`W = w_area * Nstart / nsim * sqrt(Ml) * Vt * Z(Ds)`. It halves-to-sixths the F4 "better than
10%" fractions and moves the median `tE` of joint detections from 73 d to 23 d. The F2
short-`tE` headline survives (0.25 -> 0.26). **RESOLVED 2026-09-16 (Steps W1 and W2,
Deviations 41 and 43):** the weight is `romanlib.event_weight()` on top of
`analysis/galaxy_model.py`, validated against OGLE-IV, and `f1`-`f4` now apply it, report `N_eff`
and refuse to run unweighted unless told to. The whitepaper's pooled numbers were regenerated.
**What is left is H3/H5, tracked as its own item below.** Original text kept below.

**What is wrong.** The per-sightline loop stops on a *count* (Deviation 40): outside Roman's
footprint almost every sightline contributes exactly 50 detections, inside it ~57-84. So a
sightline's share of the per-event table is set by the stopping rule, not by how many events the
sky there actually produces -- which varies with stellar density, extinction and optical depth
across the scan. Any statistic pooled over events from several sightlines therefore weights
sightlines roughly equally per event rather than by expected yield. `w_area` corrects only the
sky-area part of that (Step E1), and **no analysis script applies even that**: `f1`-`f4` never
call `romanlib.area_weight()`, and nothing in `analysis/` reads `nsim`.

The code already knew. The map-file writer's comment (`Bulge_LSST.cpp`, Step E1 columns) says
the columns were added so that "the correct pooled weight, `w_area/nsim`, is computable" -- the
join exists (`lon`, `lat` in both files), the weight was never applied. Whether `w_area/nsim` is
itself sufficient, or also needs the per-sightline star count (`log10 Nstart` is in the map file)
and event rate (`log10 Gamma`), is **not yet derived** -- derive it before applying anything.

**Why it matters scientifically.** It does not touch anything computed per event (a ratio of two
forecasts for the same event, H3's paired comparison) or per sightline. It does touch every
pooled fraction and median: the whitepaper's "fraction measured better than 10%" numbers, F2's
per-bin medians, F3's map, F4's distributions. Inside the footprint the sightlines are fairly
alike and the effect may be small; pooled across footprint and outside, where detection rates
differ most, it may not be.

**Why deferred.** It is an analysis-layer fix, not a simulation one, and it needs the weight
derived first. It is also coupled to E2: any new stopping rule changes the per-sightline counts,
and so changes how much this matters.

**What the fix involves.** (1) Derive the per-event weight for a pooled statistic from how a draw
maps to sky events (`nsim`, `Nstart`, `Gamma`, `w_area`). (2) Add it to `romanlib.py` beside
`area_weight()`, joined on (`lon`, `lat`). (3) Re-run one headline pooled number both ways --
the six-field "better than 10%" fractions are the cheapest -- and record the shift. Only then
decide whether F1-F4 need regenerating.

---

## The per-mass and per-parallax detection efficiencies have never been computed

**Status: RESOLVED 2026-09-17 the same day it was found -- Deviation 46.** All six missing
curves are now computed (measured: their `EfLMC` columns went from summing to exactly 0 to being
populated, at no measurable cost), and the `Nhalo`/`Nself` columns no longer print `-nan`. The
original text is kept below because it explains why the efficiency-versus-mass curve could not be
drawn from any run before this date -- including the v3 production run, whose `EfLMC5.dat` still
has those columns as zeros and cannot be repaired without re-running.

Original text:

**What is wrong.** `FuncMl()` and `FuncPi()` bin an event into a lens-mass or parallax bin so the
detection efficiency can be reported against those axes, exactly as `FunctE()` does for `tE`.
**Both call sites are commented out** (`Bulge_LSST.cpp:1796-97`), so `NsMl`/`NdMl` and
`Nspi`/`Ndpi` stay at zero and the corresponding columns of `EfLMC<tag>.dat` are columns of
zeros. `FunctE` is called and its `tE` efficiency is real; the other two are not.

**Why it matters scientifically.** Detection efficiency versus lens mass is *the* plot for a
black-hole population study -- it is how a yield becomes a statement about which masses a survey
can find. Right now that curve cannot be made from the simulator's own output; it would have to be
reconstructed in the analysis layer from the per-event table, which is possible but is a different
estimator with different sampling noise.

**Why deferred.** Re-enabling them changes what a run computes and costs CPU per draw, so it is a
behaviour change that belongs in its own step with its own before/after numbers -- not a silent
fix folded into the population work.

**What the fix involves.** Uncomment the two calls, confirm the histogram accumulators are
incremented on both the simulated and detected sides (the `tE` pair is the model), and check the
cost per draw. `FuncMl`'s stale `CHECK(l.Ml >= 3.0)` was already corrected to the population's own
lower bound, so the assertion no longer blocks it.

---

## H3's satellite-parallax numbers are not event-rate weighted, and its existing data cannot be

**Status: H5 is DONE (2026-09-17, Deviation 44). H3 stays open until a new paired run.**

**What is wrong.** `h3_satellite_parallax.py` now weights everything it pools, but the only
paired file that exists (`2026-09-11_h3_footprint/h3_pair.dat`, 2,673 events) was written before
the simulator recorded `Ml`, `Dl`, `Ds` and `Vt` in that row. The weight needs
`sqrt(Ml) * Vt * Z(Ds)` and **none of it is recoverable**: `tetE` and `piE` give `Ml`, and
`tetE/tE` gives `murel`, but `pirel = 1/Dl - 1/Ds` is one equation in two unknowns, so neither
`Vt` nor `Z(Ds)` can be had. Joining to the v3 event table does not work either -- the paired run
used `--start-index 518` and consumed a different stretch of the RNG stream, so **0 of its 2,673
events match any v3 event** (measured, on `lon, lat, tE, u0`); only the 134 sightlines coincide.

So the satellite-parallax section still quotes raw-sample statistics: the 0.9924 median
`sigma_piE` ratio over 528 Roman-covered events, and the 36.7% improving by more than 1%.

**Why it matters scientifically.** Each of those is an aggregate over events, and the raw sample
over-represents long, slow, massive lenses about tenfold (Deviation 41), so each will move. The
direction is not obvious: H3's gain concentrates in SHORT, high-magnification events, which the
weight favours, so the satellite-parallax median could well improve rather than shrink.

**What does NOT move:** the pairing itself. Every H3 and H5 comparison is the same event
evaluated twice, and that per-event ratio carries no weight. The sign test and the bit-exact
control are unaffected.

**Why deferred.** It needs a production run, not an analysis change. The simulator now writes the
four columns (Deviation 44), so the work is already done on the code side.

**What the fix involves.** Re-run `--pair-satellite --events 60 --lenses 20 --nerr 2
--stride-roman 5 --start-index 518` on a binary built from 2026-09-17 or later, then run
`h3_satellite_parallax.py` with `--map` and `--log` pointing at that run. The script refuses a
legacy file unless `--unweighted` is passed, so it cannot silently report raw medians. The
bootstrap interval and the sign test still need weighted versions -- both are currently computed
on the raw sample and are labelled as such in the script's output.

---

## The map file's stream is never flushed, so a killed run loses its last sightline rows

**Status: the code is FIXED (2026-09-16, Deviation 42) -- `fil3`, `fil2` and `fil2b` now flush
per sightline, verified by killing a stub run under both binaries. This entry stays open for the
DATA: the v3 files are still short six map rows and six efficiency blocks and still carry the
merged line 687. Any v3 analysis needing per-sightline quantities must take `nsim` from the run
logs, and read the map through `romanlib.load_sightlines()`, which skips the bad line.**

Original text, describing the code before the fix:

**What is wrong.** `fil3` (`MapLMC5.dat`) is an `std::ofstream` opened `ios::app` and written
with `"\n"`, never flushed per row. Every production pause so far has been a kill, and a kill
discards the buffer. The v3 map file lost six chunk-1 rows (l = 0.281, b = -0.94 .. -0.14) and
carries one 122-field line: the lost buffer's 52-field fragment with the next run's first row
appended. The 2026-09-10 re-run was stopped the same way and lost the same tail, so comparing
the two copies could not catch it.

**Why it matters.** The map file is the only on-disk source of per-sightline `nsim`, `Nstart` and
`Gamma`, which the pooled weight needs. `test5.dat` is unaffected because it is opened, written
and closed per event.

**Why deferred.** It is a C++ change and outside the current step. The v3 data are recoverable:
`run.log`/`run2.log` print `nsim` for every sightline, and `Nstart` is a deterministic function of
(l, b).

**What the fix involves.** Flush `fil3` after each row, one line, negligible cost at ~1 row per
minutes of CPU. Any reader should also reject lines whose field count is not 70. For v3, either
patch the six rows and split line 687 from the logs, or have the analysis read `nsim` from the
logs.

## Roman's photometric error model is a flagged placeholder, and every Roman error bar inherits it

**What is wrong.** `errRomanM()` (`Bulge_LSST.cpp`, just above `main`) carries
`TODO(Ali): PLACEHOLDER`. It is a **nearest-neighbour** lookup of magnitude in
`files/sigma_roman.txt` -- no interpolation -- and the comment above it says outright that what
that file encodes has not been confirmed: whether a flat mag-vs-error curve with no per-visit
depth term is right for Roman, and whether the file's magnitude sampling is fine enough for
nearest-neighbour to be acceptable. Rubin's counterpart, `errlsstM()`, uses each visit's own
5-sigma depth and is not affected.

**Why it matters scientifically.** It sets `errgR`, the F146 photometric error on every Roman
epoch. That feeds Roman's detection chi-square, every Roman photometric Fisher matrix, and --
from Step S1 on -- the Roman error bars drawn in every sample light-curve figure. A coarse file
sampled nearest-neighbour makes the error a step function of magnitude, so a figure of a
brightening source would show the error bars jumping in discrete steps rather than shrinking
smoothly.

**Why deferred.** Found while designing Step S1, which is a dump of existing values and must not
change them: fixing the model changes every Roman forecast, which is a result-moving change to be
made deliberately and re-verified with `fishertest`, not slipped in beside a plotting step.

**Provenance (2026-09-29):** the user identifies the curve as Penny et al. (2019)'s and cites it so
in `Report/overview/`; its 5-sigma point is F146 = 25.52 and it ends at 27.0 mag. See also the entry
"Roman's F146 single-visit depth is a 29 mag placeholder" (end of file): `thre[6]` = 29.0 does not match it.

**What the fix involves.** Establish what `sigma_roman.txt` is (its provenance, and whether it is
per exposure -- which is what one row of `RomanBaseline.dat` is), then replace the lookup with
linear interpolation in magnitude, and confirm with `fishertest` and a stub run how far Roman's
sigmas move. Until then, a figure caption quoting Roman error bars should say they come from a
tabulated per-exposure model.

## The legacy magC0.dat / datC0.dat demo dump is dead code

**What is wrong.** Inherited from the LMC code. `fil4` (`magC0.dat`, a dense model curve) and
`fil5` (`datC0.dat`, noisy sampled epochs) are written only when `flagm > 0`. `flagm` is set only
inside `if (test < 1.0 && save < 0 && gPop->legacyId == 1)`, `save` is initialised to 0, and its
only increment is inside that block. The condition can never be true, so both files are always
0 bytes (confirmed in the repo and in `runs/prod_bh_20260917/`). The `BHLSSTMONTS.dat` write,
gated on the same `flagm`, is dead with it. `fil4`/`fil5` are still opened, truncated and
checked on every run.

**Why it matters.** Mostly it does not -- except as a trap. The obvious "fix" of making the gate
reachable would silently change the run, because the `fil5` block calls `RandN()` inline and
consumes the RNG stream, which changes every subsequent event. It is also superseded: Step S1's
`--dump-samples` does the same job for both telescopes, in both frames, for any population,
without touching the RNG (Deviation 50).

**Why deferred.** Removing it is a cleanup with no scientific effect, and it touches the file
setup that `runctl.sh`'s resume logic depends on (DEVIATIONS.md's output-file table lists
`magC0.dat` and `datC0.dat` as truncated on start and appended on resume).

**What the fix involves.** Delete the `flagm` gate's three write blocks, the `fil4`/`fil5`
streams, `save`, and `initial` (which only exists to widen that dead block's time range); update
the output-file table in DEVIATIONS.md; confirm a stub run is byte-identical before and after.
Do **not** instead make the gate reachable.

## The table's t0 is not the observed peak, and t0zone / dt_edge / nep_pk are all computed from it

**What is wrong.** `lightcurve()` measures the observer's displacement from Earth's position at
**t = 0** (the start of the simulation), so `u0` and `t0` are the closest approach and its time
for the straight line in THAT gauge -- in effect, as seen by an observer parked at Earth's t = 0
position. The event an Earth observer actually sees is the parallax-bent trajectory, whose
closest approach is at a different time and a different separation. Found while drawing the
Step S2 sample figures (Deviation 51). On the ten S1 acceptance-test events (bulge, unweighted,
so indicative only):

| event | tE [d] | observed peak - t0 [d] | as a fraction of tE | u0 (table) | u_min (observed) |
|---|---:|---:|---:|---:|---:|
| both_003 | 82.4 | +14.0 | 0.17 | 0.514 | 0.249 |
| both_006 | 27.0 | +5.9 | 0.22 | 0.738 | 1.969 |
| any_004 | 134 | +3.6 | 0.03 | 0.857 | 0.229 |
| gap_filler_010 | 6.37 | -2.7 | 0.42 | 0.318 | 0.123 |
| astrometric_001 | 1142 | +351 | 0.31 | 0.355 | 0.122 |

Median |shift| over the ten is ~0.12 tE. `both_003` is the sharp case: the table says its `t0`
fell **in a mid-mission gap**; its observed peak is **inside a high-cadence Roman season**, and
Roman recorded it.

**Why it matters scientifically.** Three table columns are computed from the parameter `t0`
rather than from the observed peak: `t0zone` and `dt_edge` (the gap-filling axes -- F2's
headline plots against `dt_edge` and splits on `t0zone`), and `nepL_pk`/`nepR_pk` (the
+-2 tE coverage window, which the S1 sample selectors also use). An event near a season edge
can be put on the wrong side of it. For short events the shift is a sizeable fraction of tE,
which is exactly the regime where the gap-filling result is claimed. The direction of any bias
is not obvious -- the shift is symmetric-looking in the sample -- so it needs measuring, not
guessing.

**What it does NOT affect.** The Fisher forecasts: the marginalised sigma(pi_E) is invariant
under the reparametrisation (it redefines u0, t0, tE as functions of pi_E, not pi_E itself),
and the photometric model and its derivatives are self-consistent in the t = 0 gauge. The
simulation's no-parallax chi-square (`dchiP`) is in the same gauge and so overstates the
parallax signal, but it feeds only the dead `BHLSSTMONTS.dat` write (see the magC0/datC0 item).
Nor the detection test, which compares the true light curve to a flat baseline.

**Why deferred.** It is a question about F2's classification, discovered in a plotting step.
Changing what `t0zone` means changes a headline number and must be done deliberately.

**What the fix involves.** Measure first: compute the observed peak time per event (the minimum
of `u(t)` in the geocentric frame -- the S2 plotter's `Geocentric.tref` does this from a dense
curve; the C++ could do it cheaply in the time loop), re-derive `t0zone`/`dt_edge` from it on a
stub run, and count how many events change zone, split by tE. If it is material, add
`t_peak`, `u_min`, `t0zone_pk`, `dt_edge_pk` as appended columns, leave the old ones, and rerun
F2 against both.

## **CRITICAL.** The extinction law is inverted: `AlAv()` takes 1/lambda twice

**Status 2026-09-22: FIXED IN CODE (Step E1, Deviation 53), NOT YET IN THE DATA.** Every
existing production table, figure, the report and the whitepaper still carry the inverted law.
This item closes when `bulge`, `bh` and `ns` are re-run and their products regenerated.

**Found 2026-09-22 while choosing sightlines for Step S3. Not fixed -- awaiting the user's
decision, because the fix moves every result in the project.**

**What is wrong.** `helper.cpp`:

```cpp
double AlAv(double lambda_um, double Rv) {
    double x = 1.0 / lambda_um;             // wavelength -> wavenumber
    return CCM89_a(x) + CCM89_b(x) / Rv;    // but CCM89_a/b take a WAVELENGTH and do 1/x again
}
```

so CCM89 is evaluated at wavenumber = lambda. The law comes out **inverted**: extinction rises
from blue to red instead of falling. `A_lambda / A_V` at the bulge's R_V = 2.5:

| band | code | correct CCM89 | code / correct |
|---|---:|---:|---:|
| u | 0.072 | 1.694 | 0.04 |
| g | 0.112 | 1.216 | 0.09 |
| r | 0.169 | 0.854 | 0.20 |
| i | 0.231 | 0.627 | 0.37 |
| z | 0.290 | 0.460 | 0.63 |
| y | 0.346 | 0.381 | 0.91 |
| F146 | 0.751 | 0.197 | 3.8 |

"Correct" was computed two independent ways: the code's own `CCM89_a/b` called with the
wavelength, and a from-scratch Python transcription of Cardelli, Clayton & Mathis (1989), which
also returns A_V/A_V = 0.999 at 0.55 um. The maps really are A_V (`maps.py`: `av = rv * ebv`),
so the error is entirely here. Introduced in `9919917` (2026-07-26, "fixed some roman values"),
before the refactor began, so **every production run carries it**: v3 (`test5.dat`), `bh`,
`ns`, and every number in `Report/populations/populations_report.tex` and the whitepaper.

**How it showed up.** Profiling the best `rubin_only` sightline (l = -0.619, b = -1.04) in the
`bh` table: source baselines averaged r ~ 19.5 and F146 ~ 28.1-29.3. A reddened bulge star
must be BRIGHTER in the near-infrared than in r; the CMD files themselves are right (first
bulge row, a 0.61 Msun K dwarf: absolute r 7.83, F146 6.79).

**Why it matters scientifically -- in the direction of every headline.** Toward the low-latitude
bulge A_V is ~5-20 mag. With the inverted law, the optical bands see a small fraction of the
true extinction and F146 sees several times too much: Rubin's sources are several magnitudes
too bright and Roman's too faint. Rubin's detections, blending and precision are overstated;
Roman's understated. The "Rubin-only ~89% of detections" split, the joint-vs-single gains, the
gap-filling result and the astrometric numbers all move, and the direction favours Rubin.

**Measured size (2026-09-22).**
- *Survey-wide*, from the `bh` table's applied r-band extinction `Ai_r` (A_V = Ai_r / 0.169,
  approximate for disk sources, which have R_V = 3.1): over 110,144 detections, A_V has median
  **3.8** (16-84%: 2.3-8.9, max 21.0); 38% have A_V > 5, 12% > 10. The mean correction is
  **r +3.7 mag fainter, F146 3.0 mag brighter** -- a 6.7 mag swing in the Rubin-Roman comparison
  for the typical detected event. (Indicative only: the detected set was itself selected under
  the bug.)
- *At the best "Rubin-only" sightline* (l = -0.619, b = -1.04): mean A_V **12.6**, so F146 was
  made 7.0 mag too faint and r 8.6 mag too bright. Its Rubin-only events -- the densest in the
  survey -- are **largely an artifact of the bug**: Roman "missed" sources the law hid from it.
- *Stub (the S1 acceptance flags, scratch build with the one-line fix)*: a lightly extincted
  patch (A_V ~ 2.3). Blended r baseline 17.75 -> 19.35, F146 21.39 -> 20.53; the source's own
  r - F146 goes from -0.7 to +1.7 (IR-bright, as it must be). Roman-only detections 31 -> 41,
  both 18 -> 15, Rubin-only 42 -> 42, over ~90 detections -- small numbers, and the RNG stream
  diverges once any magnitude changes a detection, so this is not a like-for-like comparison.

**Why deferred.** It is a one-line fix, but it invalidates every production data product and
means re-running `bh`, `ns` and `bulge` (~13.5 h CPU each) and regenerating every figure, the
report and the whitepaper numbers. That is the user's call, not a side effect of a plotting step.
It also blocks Step S3: sample figures drawn now would illustrate the wrong photometry.

**What the fix involves.** `return CCM89_a(lambda_um) + CCM89_b(lambda_um) / Rv;` -- and a unit
test pinning A_V/A_V = 1 at 0.55 um and the table above, since nothing caught this for two
months. Then: `fishertest` (unaffected -- it uses no photometry), a stub comparison, the three
production runs, and every downstream product.

## RESOLVED (Step Y, Deviation 54, 2026-09-23) -- No absolute yield exists: the "detected events" counts are Monte Carlo sample sizes, and the rate has no compact-object abundance

**Resolved by** `analysis/y1_absolute_yield.py` and the `romanlib` rate functions. The diagnosis
below was right: the counts were sample sizes at F = 1. The speculation below about the CAUSE of
the legacy shortfall was not tested and is superseded -- see the new item on the legacy `Neven`.
The yield normalisation is validated to 0.1% against the C++ optical depth, and bracketed by
Penny et al. (2019). Kept for the record; do not act on the "candidate reasons" list.

**Raised 2026-09-22 by the user:** the `bh`/`ns` detection counts (110,144 / 93,685 in the report's
table) looked ~4 orders of magnitude above the literature -- Sajadian & Sahu 2023 (AJ 165, 96):
Roman detects 56-77 isolated stellar-mass BHs (2-50 Msun); Sajadian & Makler (arXiv:2608.16448):
Rubin detects ~2 (LMC) and ~0.3 (SMC) IBH events for a BH mass fraction F = 5e-3.

**What those counts are.** The number of Monte Carlo draws that passed detection. It is set by the
compute budget (`--events 300 --lenses 50`, capped at 50,000 draws per sightline), not by the sky,
and every draw's lens IS a BH (or NS). It is not a yield and must not be read as one. The report
says so in its caveats ("No absolute yield is quoted") but its Table 1 row is labelled "detected
events", which invites exactly that reading.

**What the code's own legacy rate implies** (`Neven = nstart * Gamma * 10`, deg^-2 per 10 yr, in
MapLMC column 49 as log10; summed as Neven x w_area over the 1612 aggregated sightlines, 59.3
deg^2, split per sightline by the unweighted detL/detR fractions):

| | joint | Rubin | Roman | both |
|---|---|---|---|---|
| `bh` | 7,111 | 6,747 | 738 | 379 |
| `ns` | 20,660 | 19,620 | 1,814 | 808 |

**These assume ALL Galactic mass is in the chosen population**: `s.opt` (the optical depth that
feeds `Gamma`) is computed from the full stellar mass density in `optical_depth()`, and nothing
multiplies it by the population's mass fraction F. With F = 5e-3 (Sajadian & Makler's value) the
`bh` numbers become Roman ~3.7, Rubin ~34 over 10 yr -- i.e. Roman is ~15-20x BELOW Sajadian &
Sahu, not 10^4 above. Candidate reasons, none yet measured: the inverted extinction law (makes
F146 ~3 mag too faint, so it suppresses Roman specifically), the 3-1000 Msun range (<M^-1/2> is
1.9x smaller than for 2-50 Msun at fixed tau), `--stride-roman 5`, a different F, and the legacy
formula's own unweighted averages of tau and 1/tE (the Deviation 41 bias, never checked here).
Rubin toward the bulge vs the Magellanic Clouds is not like for like (source density and tau are
each ~10-100x higher toward the bulge).

**Why deferred.** Every published number in the report is a fraction, median or ratio, which a
global F cancels out of. An absolute yield is a new deliverable, not a fix; it needs the user's
choice of F per population and a validation of the legacy `Neven` against a known rate (OGLE-IV /
KMTNet for the `bulge` population) before any number is quoted.

**What doing it would involve.** Multiply the rate by F_pop (a per-population constant in
`LensPopulation`), recompute Neven with the event-rate weight rather than unweighted draw
averages, validate on `bulge` against OGLE-IV/KMTNet rates and Roman's Penny et al. 2019 yield,
then compare `bh` against Sajadian & Sahu on Roman's footprint with their mass range. Meanwhile,
relabel the report's "detected events" row as the Monte Carlo sample size.


## The model's event rate is too flat in Galactic latitude

**NOTE 2026-09-30 (Deviation 64):** the comparison counts I < 21 sources through the model's
extinction tables, which are 3-6x too thin within 1 deg of the plane (Deviation 63). Its conclusion
"footprint yields are more likely under- than over-stated" therefore does not carry over to the
dust-corrected yields; the report no longer says it. Re-measure on the re-run with rebuilt tables
(or with the U5 reference dust) before chasing the density model.

**RE-MEASURED 2026-09-26 on the post-fix bulge run** (`figures/yield_20260925/bulge/y1_yields.md`):
model/observed rate 1.38 at b = -5.1 falling to 0.48 at b = -0.1, crossing 1 near b = -2.3 (was
1.47 -> 0.60, crossing -1.7). Still too shallow; at the footprint (b ~ -1.4) the model is 25-40% low.

Measured in Step Y against OGLE-IV (Mroz et al. 2019, southern fields, sources with I < 21,
u0 < 1, 0.1-deg latitude bins, >= 500 draws each). The ratio of modelled to observed event rate
per star runs monotonically from 1.47 at b = -5.1 deg to 0.60 at b = -0.1 deg, crossing 1 near
b = -1.7 deg. The optical depth ratio behaves the same way but more weakly, 1.41 down to 0.90.

**Why it matters.** The model overproduces events in the outer fields and underproduces them in
the innermost ones, so any yield summed over the full 59 deg^2 scan carries a systematic of a few
tens of per cent, and the SHAPE of the sky maps (figure `p7_sky`) is wrong in the same way.
Yields confined to Roman's footprint, which sits at small |b|, are much less affected -- but they
sit where the model is if anything LOW, so footprint yields are more likely under- than
over-stated.

**Why deferred.** It is a property of the Besancon density normalisation and/or the bulge bar
model, not of anything Step Y added, and diagnosing it means going back into the density profile
and the star-count normalisation `Nstar_k`. That is a separate investigation from the yield, and
the yield is validated independently (tau against the C++ to 0.1%, Penny et al. 2019 for Roman).

**What the fix would involve.** Compare `G.density_profile` and the CMD star counts against the
OGLE-IV star-count and optical-depth maps bin by bin in (l, b), and decide whether the discrepancy
lives in the lens density, the source counts, or the extinction that sets which sources are
bright enough to enter the I < 21 sample. The extinction law was inverted when these data were
made (Deviation 53), which affects the I < 21 cut directly, so this must be re-measured on
post-fix data before being chased.

## The legacy map-file `Neven` is ~300x below the rate-based yield, and nobody knows why

Step Y computed both. Summed over the scan with `w_area`, the map file's own `Neven` gives
4.78e4 (bulge), 7,111 (bh) and 2.07e4 (ns) at F = 1, against 3.94e5, and the new estimate is
~300x higher on the stub. Two differences are known -- `Neven` averages `eps/tE` over DETECTED
events only and carries no `sqrt(M) v_t` weight (the Deviation 41 bias), and its `nstart`
(4e5/deg^2) is 0.002 of the new `Nstar` on the stub and 0.04-0.06 in production -- but neither
has been shown to account for the size of the gap, and `nstart` does not scale with `nsim`.

**Why it matters.** `Neven` is the number the legacy code has always printed, and it feeds
`figures/*_sky`. If it is wrong by a large factor, every legacy yield statement is too, including
any inherited from the advisor's LMC+ELT work.

**Why deferred.** Nothing in the current analysis uses it -- Step Y reports it for continuity and
explicitly does not trust it. Diagnosing it means reading the legacy `Neven` accumulation in
`Bulge_LSST.cpp` line by line against the new formula, which is a step of its own.

**What the fix would involve.** Either derive `Neven` from the same weighted sum Step Y uses and
delete the legacy accumulator, or find the missing normalisation. Do not "fix" it in passing.

## RESOLVED 2026-09-30: the Sajadian & Sahu (2023) characterisation gap was the sample definition (Deviation 65)

Not a forecast difference. SS23's published code counts u0 <= 1 events peaking inside Roman's mission,
UNWEIGHTED by the event rate; counted that way ours are 1.07 (tE), 1.40 (piE), 1.28 (Ml) times theirs.
Left: tetE 0.81x, from their 1.6-8x smaller per-exposure astrometric errors (a modelling choice,
covered by the astrometric-floor items). Full record, numbers and code references: DEVIATIONS 65.

## One DET_ANOMALY event in the post-extinction-fix bulge run (2026-09-24)

**What is wrong.** `runs/prod_bulge_20260924` (commit a5028fe) reports
`ANOMALY(single-not-joint) 1` in its RUN TOTALS: one event where a single-survey detection
test passed and the joint test did not. Since Step H7 all three tests share one fixed bar and
chi2 accumulates over both instruments, so the code's own warning says this count "should be
ZERO by construction" and a non-zero value means a sign convention has crept into the signed
lensing statistic. The 09-17 bh and ns runs had 0.

**Where.** Scan position 597 (0-based; the 598th `NEW STEP`), l = -0.319 deg, b = -0.94 deg,
inside Roman's footprint; run.log line ~7,216,041.

**Why it matters, and why it is deferred.** detJ is still made monotone downstream, so no output
table is wrong; one event in 85,061 detections changes no number. But it is evidence that the
"by construction" argument has a hole, possibly one the extinction fix exposed (fainter Rubin
sources, brighter Roman ones). Deferred because the production runs are in flight and it does
not bias them. **Fix would involve:** finding the row in test5.dat (detL/detR set, raw joint
test failed), recomputing its three delta-chi2 values, and checking the sign of each survey's
contribution.

## Roman "detects" events whose peak falls years outside its mission (2026-09-24)

**MEASURED 2026-09-26 (Deviation 55):** events with no Roman epoch within +-2 tE of t0 are 2.9% (bh),
3.9% (ns) and 4.4% (bulge) of the Roman-detected yield (`figures/yield_20260925/y5_peak_coverage.md`).
Peaks outside the mission are 46% of the bh Roman yield, but long events keep Roman on the magnified
part. A peak-coverage requirement would lower Roman yields by at most 4.4%. Kept open only for the
report-convention decision (quote Roman yields with or without it).

**What was seen.** Sample astrometric_b003 (bulge, batch b): detR = 1 with **0 Roman epochs
within +-2 tE of the peak** -- the peak is ~4 yr after Roman's last season (t0zone = 2), so Roman
saw only the far photometric wing. All three batch-b astrometric events have t0zone = 2.

**Why it matters.** Plausible physics -- 50,000 Roman exposures can accumulate delta-chi2 >= 500
from a small, slow wing brightening -- but it changes how "Roman detections" in the yields should
be read: some are wing detections of events Rubin sees peak, with no Roman coverage of the peak.
It also bears on the detection test's sensitivity to baseline systematics that the simulation
does not model (a real survey fits the baseline; a slow wing is degenerate with it).

**Deferred because** the production runs are in flight and this is a question of
interpretation, not a code defect. **What resolving it would involve:** from the post-fix tables,
the fraction of detR = 1 events with t0zone = 2 (or with nep_pk_R = 0), weighted; and deciding
whether the report should quote Roman yields with a peak-coverage requirement as well.


## One bulge event violates sigma_joint <= sigma_Roman on the lens mass by a factor 8.3 (2026-09-26)

**What is wrong.** `f3_characterization_map.py` on the post-fix bulge detections
(`figures/wp_20260926/f3.log`) reports `('Ml', 'roman', 8.34, 1)`: one event whose joint lens-mass
error is 8.3x its Roman-alone error. Every other violation in the same check is round-off sized
(<= 1.0007 on tE/piE, 10-20 events). The joint information matrix is the sum of the per-survey
ones, so a real violation of this size on sigma(Ml) points at the mass propagation rather than the
matrices -- e.g. `ErrorCal` taking piE from the astrometric matrix for one partition and the
photometric for the other ("keeps whichever is tighter", per partition).

**Likely benign:** `f4_fisher_precision.py` reports the same single point above the line and
attributes it to a photometric matrix with condition number > 1e9 (`figures/wp_20260926/f4.log`),
i.e. double precision has lost the answer and neither forecast is meaningful. Keep open until the
row is inspected.

**Why deferred.** One event in 85,061; no pooled number moves. **Fix would involve:** finding the
row (relMl_J > relMl_R, both valid), printing its per-partition sigpiE (photometric and astrometric)
and sigtetE, and checking which piE each partition's mass used.

## Roman's detection rate differs from Penny et al.'s |u0| < 3 figure by ~25% (2026-09-26; re-read 2026-09-30)

**STATUS 2026-09-30, later (Deviation 66): the "bracket" below is a misreading -- still open.** Penny
et al. (2019, Sec. 2 and Table 2) give DETECTIONS at both |u0| < 1 (~27,000) and |u0| < 3 (~54,000);
three times as many events *occur* at |u0| < 3. Our detections (u0 to 3) therefore compare with 63.8
directly, not "below" it. As simulated 80.4 +- 1.4 = (26 +- 2)% above; dust-corrected 49.4 +- 2.9 =
(23 +- 5)% below. The dust correction flips the sign; it does not resolve the comparison. The overview
report now says this. The populations report, postfix report and whitepaper were corrected the
same day (Deviation 67); Deviation 54's verification bullet is left as the historical record. The
"fix would involve" below still stands and is now the only way to settle the comparison: restrict
our count to Penny's six high-cadence seasons, their |u0| range and t0 window, and redo it with the
corrected dust (or on the re-run).

**Superseded STATUS 2026-09-30 (Deviation 64), kept for the record:** "accounted for by the dust,
pending the re-run ... (23 +- 5)% below Penny's |u0| < 3 figure and above the |u0| < 1 one -- where
this entry says it should sit."

**What is wrong.** Per deg^2 of footprint per day of Roman coverage (693 d), Roman detections
(u0 < 3, delta-chi2 >= 500) are 80.4 post-fix, against Penny et al. (2019) 31.9 (|u0| < 1) and
63.8 (|u0| < 3). Pre-fix the figure was 40.4, inside that bracket, and the whitepaper called the
bracket "predicted, not fitted". Our count should sit BELOW the |u0| < 3 figure, since the
delta-chi2 cut rejects most low-magnification events.

**Why it matters.** It is the one external check on Roman's absolute normalisation; the internal
tau check (1.001) validates the rate constants, not the source counts or the detection efficiency.

**Why deferred / candidates.** Known differences push the right way -- ten seasons including four
low-cadence ones (Penny: six high-cadence), a different detection criterion, a different field
layout and source-count model -- but none has been quantified. The OGLE-IV comparison says our
footprint rate per star is 25-40% LOW, which makes an over-count from the rate unlikely; the source
counts or the detection efficiency are the places to look. **Fix would involve:** restricting our
count to Penny's six high-cadence seasons and |u0| < 1, and comparing per-star rather than per-area.

## The simulator cannot yet be split safely across cluster jobs: no seed option, no end index (2026-09-26)

**What is wrong.** The RNG is `mt19937_64` seeded with the compile-time constant `seed = 42`
(`Bulge.h:32`), and `--start-index` skips sightlines *without* consuming random numbers
(`Bulge_LSST.cpp`, `if (iScan < cfg.startIndex) continue;`). A scan split into chunks by start
index therefore starts every chunk from the same random stream: chunk k's first sightline sees
exactly the deviates chunk 0's first sightline saw. There is also no `--end-index`, so a chunk
runs to the end of the scan, and the analysis layer assumes one run directory per population.

**Why it matters.** Draws in different chunks would be correlated (different sightlines, same
uniform deviates), so pooled Monte Carlo errors would be underestimated. Harmless for the laptop
runs (one chunk each, or resumes that redo a sightline on purpose), but it is exactly what a
cluster production run would do.

**Why deferred.** No cluster run is scheduled yet; it is named, without the technical detail, as a
challenge in `Report/overview/` "Next steps" (the report's earlier cluster section, with the cost
table below, was removed 2026-09-28 at the user's request). **Fix would involve:** a `--seed` option written to `run_provenance.txt`, an
`--end-index`, a documented per-chunk seed rule, a merge step for tables and map files that
checks every sightline appears exactly once, per-sightline timestamps in the log (cost is
currently unmeasurable by stratum), and a validation that a chunked stub run reproduces a
single-process run's pooled numbers within Monte Carlo error.

**Related fact found while costing it.** Every footprint sightline stops on the `--events 300`
target, not `--lenses 50`: the bh and ns tables hold exactly 44,100 = 147 x 300 footprint rows
(bulge 44,119), with 28-55% of footprint draws detected. The whitepaper's "Stopping criteria"
paragraph ("the second floor binds") is true of the outside sightlines only. So a finer
`--stride-roman` grows CPU time but barely grows the tables, and `--events` is the lever on
footprint statistics.

**Cost projection (moved here 2026-09-28 from the old onboarding report's cluster section, so it is
not lost).** Units are the laptop CPU-hours of the measured runs (i7-6500U, two cores, under
hyper-threading contention); a cluster core is probably faster, which is not assumed. The CPU split
between footprint and outside sightlines is unmeasured post-fix (no timestamps in `run.log`): pre-fix
a profiled footprint sightline cost ~340 s against ~7 s outside (~80% of a run in the footprint),
and the fix grew the outside share, so the footprint share is bracketed at f = 40-80%. Footprint cost
scales with the footprint sightline count; outside cost and table size stay as they are.

| `--stride-roman` | footprint step | footprint sightlines | CPU-h (bulge / ns / bh) | footprint N_eff | table size |
|---|---|---|---|---|---|
| 5 (current) | 0.10 deg | 147 | 16 / 23 / 27 (measured) | 1x | 3-6 GB |
| 2 | 0.04 deg | 907 | 50-82 / 72-120 / 82-136 | ~6x | +1-2% |
| 1 | 0.02 deg | 3656 | 170-320 / 250-470 / 280-530 | ~25x | +6-12% |

Raising `--events` adds draws mainly in the footprint (where it binds) and so grows footprint
statistics without changing the sky grid; raising `--lenses` adds draws mainly outside and grows the
table roughly in proportion. Which target binds where should be re-checked in a pilot. The three
populations are independent runs and parallelise trivially. Operational traps for such a run have
their own items: stale git stamp (`make clean && make` after the last commit), `--dry-run`
truncating outputs, append-mode outputs, ~170 s start-up per chunk, silent OOM kills in analysis
(run under `analysis/memrun.py`); resume with `runs/runctl.sh stop/continue`.

## Roman's F146 single-visit depth is a 29 mag placeholder, and the image-resolution count used it (2026-09-29)

**What is wrong.** `thre` in `Bulge.h` gives F146 a single-visit depth of 29.0 mag, commented
"(value needs to change)". Roman's own photometric error table, `files/sigma_roman.txt` (the noise
model every simulated Roman epoch uses), reaches 5 sigma at 25.52 mag and 3 sigma at 26.11, and
ends at 27.0 -- beyond which `errRomanM()`'s nearest-neighbour lookup returns 0.83 mag, an
underestimate for a 28-29 mag point. `RomanBaseline.dat`'s own `sig5` column says 24.0 and is unused.

**Why it matters scientifically.** Three places read `thre[6]`: (1) the Step R1 image-resolution
count (Deviation 48), which calls the faint minor image detectable down to 29 mag -- and resolution
happens exactly as the minor image fades; (2) the "could Roman see it magnified" acceptance gate;
(3) the per-epoch recording gate. Measured by `analysis/u2_resolution_depth.py` (Deviation 58):
at a 25.52 mag depth Roman's P(resolvable) falls from 65.5 to 45.4% (bh, D=5), 14.9 to 2.9% (ns),
5.1 to 1.4% (bulge); the overview report now quotes the corrected values. For detection it is
small: 0.3-1% of Roman's detections have a baseline fainter than 25.5 mag, 0-0.02% fainter than 27.

**Why deferred.** Changing `thre[6]` moves every Roman output (acceptance, recorded epochs,
resolution) and needs a re-run of all three populations; the resolution numbers are corrected
analytically meanwhile, validated to 99.5-100% per event against the simulator at 29 mag.

**What the fix involves.** Set `thre[6]` to the error table's 5-sigma point (25.5) -- or better,
derive it from the same table at run time so the two cannot drift -- extend or clamp
`sigma_roman.txt` beyond 27 mag, re-run `fishertest`, and fold the change into the next production
runs (it belongs with the Roman photometric-model item above, since both are the same curve).

## RESOLVED 2026-10-01 (Deviation 69) — The scan's corner cut removes the wrong corner (raised by the user, 2026-09-29)

**Resolution.** The box and its corner cut are gone; the scan region is every point within
2 x 1.75 deg + the field reach (3.944 deg) of any Roman field centre, either roll -- the
"distance test" proposed below. Kept for the record:

**What is wrong.** `Bulge_LSST.cpp` drops every grid point with `lon < lx and lat > bx`
(Bulge.h: lx = 1.0053 - 0.2 - 1.75 = -0.9447, bx = -1.64 + 0.2 + 1.75 = 0.31), i.e. the upper corner
on the NEGATIVE-longitude side (upper right in a map with l increasing leftward, as in the overview
report's Fig. 6; the figure reproduces the run's 1,829 sightlines exactly, so it is the code, not the
plot). The user's intent was to drop the corner NOT needed to cover Roman's Galactic-centre field.
Measured with Rubin's real pointings (`analysis/u3_footprint_offset.py`, Deviation 60): the cut
corner is imaged in Rubin exposures that also contain the GC field over 4.85 deg^2 (52,866
cell-visits); the opposite upper corner over 2.65 deg^2 (17,374). The cut therefore removes 4.9 of
the 57.5 deg^2 that Rubin images together with a Roman field (the scan keeps 91%).

**Why it matters.** Only Rubin-only events outside Roman's footprint are affected: estimated from
simulated sightlines with the same Rubin visit counts, the corner would add ~1.0% (bulge), 2.4% (bh),
1.2% (ns) to the whole-scan Rubin yields. No footprint result changes.

**Why deferred.** A code change to the scan region needs new production runs; the effect is small
and confined to whole-scan Rubin totals.

**What the fix involves.** Decide which corner (if any) should be cut -- the opposite corner is still
imaged with the five-field block over 3.56 deg^2, so neither upper corner is empty of Roman-overlapping
Rubin exposures -- then flip the l inequality (or replace the rectangle cut by a distance test: keep
a point if it is within 1.75 + 1.75 + 0.30 deg of any Roman field centre, which is the actual
design intent), update the provenance line, and re-run with the next production set. The scan box
itself (l1/l2/b1/b2) was built from an older field layout (l -0.219..1.413, b -1.64..-0.85) that
excludes the current GC field; with the distance test that stops mattering.

## FIXED IN CODE 2026-10-01 (Deviation 70); results pre-fix until the re-runs -- CRITICAL: the dust is three to four times too thin within 1 deg of the Galactic plane (2026-09-29)

**STATUS 2026-10-01.** Plan steps 1-3 below are done (Deviation 70): `maps.py` builds
`files/ext/ext_tables.dat` from DECaPS + Marshall (k re-measured 0.0830) on a regular grid;
`readExtinction` replaces `readBayestar` and refuses bad tables; validated against VVV (A_V(8 kpc)
tables/VVV 0.88-0.99 in every |b| bin, was 0.17 near the plane) and on a pilot (every draw's A_r
matches the tables). A ~10% residual remains on the Galactic-centre field (new entry below).
**Still open: steps 4-5** -- the production re-runs, and checking them against U5's corrected
numbers. Every number in the existing reports remains pre-fix until then. Remove this entry when
the re-runs are checked.

**CURRENT STATE 2026-09-30 (Deviation 63; supersedes the Deviation 62 update below).** Verified,
not assumed: maps.py follows the dustmaps-documented rule (Bayestar19 north of dec -30, DECaPS
south); 139 of Roman's 147 sightlines took Bayestar, whose OWN reliable_dist flag is false at
>= 4 kpc on all 139 -- every bulge source sits beyond the range Bayestar vouches for. DECaPS has
data at all 1,829 scan sightlines, but it is sensitive to A_V ~ 12 only and is SATURATED on the
Galactic-centre field even where its flag says reliable. Independent check, the VVV E(J-Ks) map
(Surot+2020, `analysis/u6_vvv_check.py`): GC-field/five-field contrast VVV 4.89, Marshall 4.46,
DECaPS 1.07, simulator 0.57; over the whole scan (A_V at 8 kpc / VVV's, all on DECaPS's scale)
the simulator is 0.17 (|b| < 0.5), 0.31 (0.5-1), 0.88 (1-1.5), 1.02-1.11 farther out; DECaPS
0.57 / 0.86 / 0.95 / 0.88-0.89. So neither map alone is right: the fix is DECaPS where it can see,
the near-infrared map where it cannot. The dust-corrected numbers in the report
(`figures/u1_20260929/u5_corrected_numbers.md`) use exactly that reference, which tracks VVV to
0.88-0.99 in every |b| bin. Size of the effect with it: Roman footprint detections x0.62-0.67
(five-field block x0.83-0.88, GC field x0.07-0.13 -- u4_hybrid.log 0.068/0.124/0.126), Roman 10% masses x0.43-0.54, Rubin footprint
x0.29-0.38, Rubin whole scan x0.44-0.46; per-event fractions move by < 4 points. The numbers under
"Why it matters" below are the Deviation 61 estimate (Marshall at 0.11) and are superseded. The
concrete fix plan is under "What the fix involves" below.

**UPDATE 2026-09-30 (Deviation 62; partly superseded -- DECaPS is NOT close to right on the GC
field): the fault is Bayestar plus maps.py's dec = -30 rule, not the optical maps in general.** DECaPS covers all 1,829 sightlines and agrees with Marshall to ~10-25%
except within ~0.5 deg of the plane (GC field 7.9 vs 24); where the model used Bayestar it has
0.20x (|b|<0.5) and 0.35x (0.5-1) Marshall's dust. Fix: drop the rule, DECaPS everywhere, near-IR
within ~0.5 deg; check against Marshall; re-run. Dust-corrected headline numbers with errors:
`analysis/u5_corrected_numbers.py` -> `figures/u1_20260929/u5_corrected_numbers.md` (the report now
quotes these).

**What is wrong.** `files/ext/` (maps.py) takes A_V(d) from Bayestar19 (dec > -30 deg) and DECaPS
(dec < -30 deg), both built from OPTICAL photometry. Toward the inner bulge their stars cannot be
seen through the dust lanes near the plane, the profiles saturate after a few kpc, and the
extinction to bulge sources is a lower limit. Against the Marshall et al. (2006) near-infrared 3D
map (2MASS; dustmaps' MarshallQuery, local copy in `dustmaps/marshall/`), at 8 kpc with
A_Ks/A_V = 0.11, the model/Marshall A_V ratio (median over sightlines) is 0.23 for |b| < 0.5,
0.38 for 0.5-1.0, and 1.1-1.3 farther out. Examples: Roman's Galactic-centre field A_V 2.9 vs ~24;
b = -0.94: 2.1 vs 7.9. Found while checking the footprint-mismatch estimate (the model's event
density rose toward the plane because its dust thinned toward the plane). Measured by
`analysis/u4_dust_check.py` -> `figures/u1_20260929/u4_dust_check.csv` (Deviation 61).

**Why it matters scientifically.** 45-47% of Roman's and 51-53% of Rubin's detected yield come from
sightlines where the model has less than half Marshall's dust. Estimated with the dust corrected
(each draw dimmed by the extra extinction, detection efficiency read off at the new magnitude):
Roman footprint yields x0.76 (bulge) / 0.81 (bh) / 0.78 (ns) -- the GC field x0.22-0.34, the
five-field block x0.95-0.97; Rubin footprint x0.57-0.61 (GC field -> 0); Rubin whole scan (bulge)
x0.56 (1.75e5 -> 9.8e4); share of detections outside the footprint 64% -> ~57%. The extra dust
costs 4.3x more in r than in F146, so every Roman-Rubin comparison moves in Roman's favour.
Per-event results (Roman-for-Rubin, gap filling, precision fractions) were not recomputed; ~a
quarter of the footprint events behind them are in the GC field.

**Why deferred.** Fixing it means new extinction files and new production runs of all three
populations; the report states the estimate and its limits (Section 5.5).

**What the fix involves (plan, 2026-09-30, Deviation 63).** Does it need code changes and new
runs? maps.py: yes (the fix itself). helper.cpp `readBayestar`: yes (robustness, not physics).
Lensing.cpp / Bulge_LSST.cpp: no -- the way extinction is applied (A_V(Ds) x CCM89(R_V) + scatter)
is not the fault. Bulge.h: only if the number of tables changes (NFILES = 2518). New production
runs: yes -- the simulator reads the tables and dims every source at run time, so no table change
reaches the results without re-running; U5 is the estimate until then.
1. **maps.py -- build every table from the reference U5 uses.** `Dust.reference_profile(...,
   "hybrid", 0.0805)` in `analysis/u5_corrected_numbers.py` is the reference implementation:
   (a) query DECaPS (`DECaPSQueryLite`, A_V = 3.32 E(B-V)) for EVERY sightline -- drop the dec rule;
   Bayestar only as a fallback where DECaPS has no data (none in this scan), and then only within
   its own reliable_dist; (b) keep DECaPS out to its last `reliable_dist` distance, beyond it add
   Marshall's further increase; (c) from the first distance at which Marshall's A_Ks/0.0805
   reaches 12 (DECaPS's stated limit), use Marshall's A_Ks/0.0805 itself; (d) 0.0805 = median
   A_Ks(Marshall)/A_V(DECaPS) at 8 kpc over Roman's five-field block, re-measure it with
   `u6_vvv_check.py` if anything changes; (e) exit non-zero if any finished table has a non-finite
   value (no more "TOTAL DROPOUT" notes); (f) write the provenance (maps, calibration, date) to a
   file that is NOT `.txt` (`readBayestar` reads every `.txt` in files/ext/ as a table).
2. **helper.cpp `readBayestar` -- refuse bad tables.** Check the stream after every row; exit
   with the file name on a failed parse, a non-finite value, a row count != NROWS, or l, b that
   change within a file (the check exists, commented out). Optional, for reproducibility: sort the
   paths before reading -- `directory_iterator` order is filesystem-dependent, and it decides ties
   in `nearestSightline()` (0.08% of draws in the current runs).
3. **Validate before any run.** (a) Re-run the u6 comparison on the NEW tables: A_V(8 kpc)/VVV
   should be ~0.9-1.0 in every |b| bin (currently 0.17 at |b| < 0.5), with no fall toward the
   plane, and GC/five-field contrast ~4.5-4.9 (currently 0.57); (b) `make extinctiontest &&
   ./extinctiontest`; (c) a short pilot, then the U4 reconstruction check that per-draw A_r matches
   the new tables.
4. **Re-run bulge, bh and ns production** together with the other fixes that need new runs (F146
   depth, corner cut, adopted GBTDS layout), so the runs are repeated once.
5. **Check the new runs against U5's corrected numbers** (they should agree within U5's stated
   errors; if not, the efficiency re-weighting misses something -- understand it before use).
   Then re-measure what U5 could NOT correct (Deviation 64): the gap-filling table (U5 re-weights
   medians by survival only; Rubin's sources dim 4.3x as much as Roman's, so the reported gains are
   upper limits), Rubin's image-resolution fractions, the joint/Roman 10%-mass ratio (re-weighted
   with F146 alone, so an upper limit), the OGLE-IV latitude comparison (it selects I < 21 sources
   through the dust), and the Penny et al. rate (dust-corrected estimate 49.4 +- 2.9 per deg^2 per
   day, inside their bracket).
Not part of this fix, stated in the report as a limit: the V-to-F146/r conversion stays CCM89 at
R_V 2.5 for the bulge; a near-infrared law for F146 from A_Ks directly would be the next refinement.

## FIXED IN CODE 2026-10-01 (Deviation 70); existing runs affected -- 78 extinction tables are empty (all NaN) and the simulator silently reads them as zero dust (2026-09-30)

**STATUS 2026-10-01.** The rebuilt tables come from DECaPS, which covers every position (no
dropouts), `maps.py` exits rather than write a non-finite value, and `readExtinction` refuses one.
The 2026-09-24 runs still carry the zero-dust patch; it disappears with the re-runs.

**What is wrong.** 78 of the 2,518 `files/ext/bayestar_*.txt` tables hold `nan` at every distance
("TOTAL DROPOUT -- needs neighbor fallback" in maps.py; the fallback was never written). All 78 lie
at dec -30.0..-29.1, just north of maps.py's dec = -30 switch, so they were sent to Bayestar, whose
southern edge is ragged there and has no data. `readBayestar()` (helper.cpp) reads each table with
`>>` and never checks the stream: the first `nan` fails the read, the rest of the file is never read,
and sources whose nearest table is one of these get essentially no dust. Measured on the production
detections: median A_r 0.0002 mag (max 0.017; half exactly 0) against ~3 mag elsewhere.

**Why it matters scientifically.** 20 scan sightlines use these tables -- a contiguous patch at
l ~ 0.1-0.7+, b ~ -2.1..-2.9, just south of Roman's footprint -- carrying 1.3-1.7% of all simulated
detections (none in Roman's footprint). Unextincted, their Rubin detections are strongly
over-produced, so the whole-scan Rubin totals are inflated. The dust-corrected estimate
(`analysis/u5_corrected_numbers.py`) now treats these sightlines as the zero-dust ones they were.

**Why deferred.** Same fix and same re-run as the CRITICAL dust entry above.

**What the fix involves.** Part of the CRITICAL dust entry's plan above (Deviation 63). (1) Rebuild
the tables from DECaPS, which has data at all 1,829 scan sightlines (removes the dropouts). (2) In `readBayestar()`, check `fin` after each row and refuse a
table that fails to parse or contains a non-finite value, instead of continuing on a failed stream.
(3) Have maps.py exit on a total dropout rather than print a note.

## Stale or wrong statements left in the older documents (2026-09-30, found in Deviations 66-67)

**What is wrong.** Three things the overview report now gets right are still wrong elsewhere:
1. The Sajadian & Sahu *characterisation* deficit (~1/3 of theirs, "concentrated on the parallax")
   is presented as unexplained in `Report/populations/populations_report.tex` (abstract, SS23
   section), `Report/postfix/postfix_report.tex` (Characterisation paragraph) and
   `Whitepaper/whitepaper.tex` (~l. 1200: "the parallax deficit has another cause"). Deviation 65
   traced it to the sample definition (their u0 <= 1, on-mission, unweighted counting).
2. The whitepaper's OpSim table caption names `baseline_v5.3.0_10yrs.db`; the Rubin visit list was
   built from `baseline_v5.1.0_10yrs.db` (Deviation 67, 23/23 rows verified).
3. `Whitepaper/refs.bib` `Biswas2019OpSim` is garbled: DOI 10.3847/1538-4365/ab4b44 does not exist;
   arXiv:1905.02887 is Biswas et al. 2020, ApJS, doi 10.3847/1538-4365/ab72f2, with a different
   title. The rest of the whitepaper bib has not been audited against Crossref.
4. (2026-10-01, Deviation 68) The SS23 paragraphs of the populations report, postfix report and
   whitepaper still use "ours vs. theirs" / "Against SS23" framing and cite the released code
   (u0 draw, unweighted counting, astrometric precision file) as if it were the paper. The overview
   report now compares with the published paper only, in "this simulation / the paper" terms.

**Why it matters.** The populations and postfix reports may already have been read; the whitepaper
is the eventual paper. Each states something a reader would take as a finding.

**Why deferred.** Not asked for: the request was the overview report plus the Penny misreading.
The whitepaper is due a full reconciliation anyway (JOINT_FIT_REFACTOR_PLAN Phase G).

**Fix would involve.** For (1), a note in each older report pointing to the resolution, and a
rewrite of the whitepaper's SS23 paragraph from the overview's (which also covers (4)). For (2), one word. For (3), a
Crossref pass over `Whitepaper/refs.bib` as was done for `Report/refs.bib` in Deviation 66.

## Rubin's field of view is a 1.75-deg circle, not LSSTCam's outline (2026-10-01, from Deviation 69)

**What is wrong.** `matchVisibleEpochs` gives a Rubin visit to every sightline within `FoV` = 1.75
deg of the pointing centre. LSSTCam's focal plane is a square-ish mosaic of 21 rafts with the
corners cut (9.6 deg^2, the same area as the circle), with gaps between sensors (fill factor ~0.9),
and it rotates with `rotSkyPos` from visit to visit (the column is in the OpSim database but not
extracted). Also: the 40 earliest OpSim visits that reach the scan (MJD 60981-61141, before the
simulation's day 0) are dropped, which shortens the Rubin baseline before day 0 for events peaking
in the first months.

**Why it matters scientifically.** Per sightline, the visit count is right on average but wrong at
the field edges (the square reaches 2.1 deg along its diagonals and 1.6 deg along its sides), and
chip gaps remove ~10% of visits at random. Rubin's per-event cadence, hence its detection
efficiency and Fisher precision, inherit that. The effect on pooled yields should be at the
few-percent level because rotation and dithering average the outline; it has not been measured.

**Why deferred.** Not one of the three pre-production fixes; the circle is the standard
approximation, and the user asked for the minor items after those three.

**What the fix involves.** Extract `rotSkyPos` in `readbaselineBulge.py`; add a Rubin coverage
predicate built from the LSSTCam raft/sensor layout (rubin_sim/`lsst.obs.lsst` geometry) rotated per
visit; the predicate interface in `matchVisibleEpochs` already allows it. Then the readbaselineBulge
reach must use the outline's maximum radius (~2.1 deg) instead of 1.75.

## The rebuilt dust is ~10% thin on the Galactic-centre field against VVV (2026-10-01, from Deviation 70)

**What is wrong.** On the adopted GC field the rebuilt tables give A_V(8 kpc) ~ 26, against VVV's
~29 on the same scale: tables/VVV 0.90 at 8 kpc (0.87 at 7.5, 0.92 at 8.5, 0.94 at 9). The
GC-field/five-field contrast is 4.74 in the tables and 5.52 in VVV (law-free). Everywhere else the
tables track VVV to 0.88-0.99 (`figures/ext_20261001/v1_ext_vvv.md`). The tables use Marshall's
near-infrared map there (DECaPS saturates), and Marshall's own contrast is known to fall ~10% short of
VVV's (U6: 4.46 vs 4.89 on the notional fields).

**Why it matters scientifically.** ~2.6 mag of A_V is ~0.5 mag in F146 (A_F146/A_V = 0.197) on the
GC field: Roman's GC-field yields are somewhat overstated. Rubin is unaffected in practice (A_r ~ 22
already). The five-field block, which carries most of Roman's yield, is unaffected (calibrated
there).

**Why deferred.** The remaining reference, VVV, is a 2D column map not in the dustmaps library (the
user asked for maps built with dustmaps), and calibrating Marshall to it on one field would add a
free parameter fitted to 37 sightlines. The residual is stated instead.

**What the fix involves.** Either scale the near-infrared part of the profile by the local VVV/
Marshall column ratio at the bulge distance (a 2D correction applied to a 3D shape: needs VVV at
every table position, ~16,000 VizieR boxes, and a choice of how to distribute the extra dust in
distance), or bracket it: U5's `aks_lo` variant (k 0.0734, +10% Marshall) gives the size of the
effect on the yields without rebuilding.

## Roman's astrometric error is a Vega-calibrated curve evaluated at AB magnitudes (2026-10-01, audit)

**What is wrong.** `errRomanA(magF146)` (helper.cpp) uses Lam et al. (2026)'s anchors, which are in
F146 VEGA magnitudes ("F146_Vega < 20.62" floor, "F146_Vega < 23.5" background; their Fig. 5). The
simulator's magnitudes are AB: the MIST bolometric-correction tables in CMD/Roman and CMD/Rubin are
headed "Roman (AB)" and "LSST (AB)". F146's AB - Vega offset is about +1.04 mag (EXOZIPPy issue #313,
computed for W149; consistent with STScI's Roman-STScI-000825 Table 4, 2MASS J 0.913 / H 1.391,
which bracket F146). So every Roman exposure's astrometric error is read ~1 mag too faint.

**Why it matters scientifically.** On the curve's slopes (0.333 dex/mag source-dominated, 0.4
background-dominated) that overstates sigma_ast by ~2.2-2.6x for typical bulge sources, and puts
the 1.1 mas floor at AB 20.62 instead of AB ~21.66. Every astrometric number is pessimistic by
this, in the opposite direction to Deviation 71's fixes. Photometry is NOT affected: Rubin's model
and sigma_roman.txt (Penny+2019, AB) are both AB.

**What the fix involves.** Convert before the lookup, m_Vega = m_AB - dAB, with dAB computed from
the F146 throughput and a Vega spectrum (synphot) rather than taken from a code issue; record it in
Bulge.h with its source; check u2_resolution_depth.py and any Python mirror of errRomanA.

## Roman's F146 saturation limit satu[6] = 12.0 is an unsourced placeholder (2026-10-01, audit)

`satu` in Bulge.h gives F146 a single-visit saturation of 12.0 mag, commented "(value needs to
change)". It gates which epochs are recorded and the image-resolution count. Penny+2019 Fig. 4
marks the single-read saturation for their exposure; STScI's WFI pages give saturation times. Set
it from a source in the right system (AB) and exposure (66 s, with up-the-ramp reads), together with
thre[6].

## Roman's photometric error table is for a 46.8-s exposure; the GBTDS uses 66 s (2026-10-01, audit)

`files/sigma_roman.txt` reproduces Penny+2019 Fig. 4 (AB; "the Cycle 7 design's assumed exposure
time (46.8 s)"; 1 mmag floor). The adopted GBTDS exposure is 66 s (Lam et al. 2026; STScI), so the
table's faint end is ~0.1-0.2 mag pessimistic, and the error at fixed magnitude ~10-20% too large.
Options: rescale (SNR ~ sqrt(t) where background-limited) or replace with a 66-s model (Wilson et
al. 2023, STScI/Pandeia). Part of the existing "photometric error placeholder" entry's fix.

## Rubin's depth gate is the SRD minimum, not each visit's own depth (2026-10-01, audit)

**What is wrong.** Rubin epochs are recorded only if the source is between `satu[fi]` and `thre[fi]`
(Bulge_LSST.cpp, the Rubin branch, and the pre-selection `Mpeak <= thre[i]`), with `thre` = the SRD
minimum single-visit depths. But each visit's own 5-sigma depth (`sig5`, OpSim's fiveSigmaDepth) is
in the visit list and already sets that epoch's error (errlsstM). For the bulge visits the real
depths are 0.1-0.4 mag shallower than `thre` in u-z and 0.3 mag deeper in y (medians above), and
vary by ~0.4 mag (16-84%) from visit to visit. So epochs with SNR < 5 are recorded in good-seeing
bands, and y-band epochs between 21.7 and the visit's depth are dropped.

**Why it matters.** It sets which Rubin epochs count, near the detection limit, where most bulge
sources are. Size not measured.

**What the fix involves.** Gate each Rubin epoch on its own `sig5` (and a saturation limit tied to
it, e.g. sig5 - 8.3 as the current constants imply, or a sourced bright limit); keep a per-sightline
median for the pre-selection. Re-derive `FWHM[0-5]` from the visit list's seeingFwhmGeom medians and
update `gama` to Ivezić 2019 Table 2 in the same step; delete the unused constants.

