# PROGRESS.md — where this project stands

**Last updated:** 2026-09-15 (later), when the pooled weight was derived and measured —
**Deviation 41**. The premise was wrong in a way that matters more than expected.
- **The weight.** Draws are rate-weighted in lens distance only, not in mass or velocity, so
  every pooled statistic needs `W = w_area * Nstart / nsim * sqrt(Ml) * Vt * Z(Ds)`, not just a
  per-sightline weight.
- **Effect on v3.** F4's "better than 10%" fractions fall 2-6x (`tE` joint 20.4% -> 10.1%, `piE`
  15.9% -> 2.6%). The median `tE` of joint detections falls from 73 d to 23 d.
- **What survives.** F2's short-`tE` gap-filling headline, 0.25 -> 0.26.
- **Step W1 is done and committed:** `analysis/galaxy_model.py` (the `Disk_model` port),
  `romanlib.event_weight()` / `kish_neff()`, and `analysis/w1_pooled_weight_check.py`, which
  regenerates the numbers above. F1 on the v3 extract is byte-identical to
  `figures/f1_results_table_v3.csv`, so nothing existing moved.
- **Step W2 is the next decision and has not been done:** apply the weight inside F1-F4, report
  `N_eff` beside each weighted number, regenerate the figures and update the whitepaper. **Until
  then every pooled distribution, fraction or median in the whitepaper is unweighted.** Per-event
  ratios and H3's pairs are unaffected either way.
- **Validated against OGLE-IV** (Deviation 41). Weighting all 6.07M draws gives an intrinsic mean
  `tE` of **24.0 d**, against Mroz et al. 2019's efficiency-corrected **22 d** in the central
  bins; unweighted it is 56.2 d. The weight is what brings the model onto the published value.
- **Map file damaged, and the cause is now fixed** (Deviation 42). `MapLMC5.dat` lost six chunk-1
  rows (l 0.281, b -0.94..-0.14) and has one merged line; `EfLMC5`/`EfLMC5B` lost the same six
  blocks. Cause: the per-sightline streams were never flushed and every pause has been a kill.
  `Bulge_LSST.cpp` now flushes them, verified by killing a stub run under both binaries (0 rows
  on disk before, 7 of 7 after). **The v3 files stay damaged**: take `nsim` from the run logs.
  §"Chunk 2" below is wrong to call the map file undamaged.

**Last updated:** 2026-09-15, when **Step G2 was started** and the advisor's paper on the parent
code (Sajadian & Makler, arXiv:2608.16448) was read against it — **§5i**. **G2 is deferred by
the user**; the designed replacement ("Option A") is in `OPEN_ITEMS.md`. **Step E2 is done**
(Deviation 40): the plan's premise was wrong — the "well-conditioned" stopping floor had counted
detections since the first commit. The floor now counts `okA[SJOINT]`; no per-bin floor was
added, because pooled precision is set by `--stride-roman`. Output byte-identical on the stub
regression; at production `--nerr 2` it changes no stopping decision. **Next, in order of
consequence:** derive and apply the pooled per-sightline weight (`OPEN_ITEMS.md`, "Pooled
per-event statistics give every sightline the same number of events") — it touches every pooled
number in the whitepaper — then E1b's decision, then G2 Option A when the user wants it. G1 needs nothing: it was
done as H3, astrometric half included. **§1 and the top of §5 below are stale** (they still
describe the v2 run as in flight and E1a as unrun); §5h and §5i are current.

**Last updated:** 2026-09-12, when Step H6 (the collaborator whitepaper) was written against
the post-H7 v3 table — see §5h. **Phase H is complete.**

**Last updated:** 2026-09-06, when Step H7 changed the detection threshold and the v2
production runs were killed because of it — see §5e.

**Last updated:** 2026-09-05, when the stratified production run was launched. The steps it
carries are E1a (footprint-stratified sampling), H1 (satellite parallax), H2 (the satellite
observable in the table) and H4 (Roman's astrometric error model). Step H5 (the astrometric
shift as a product) is also done and is analysis-only.
**The v2 runs were killed on 2026-09-06 at 734 and 618 sightlines, deliberately, because
Step H7 changed which events count as detected — see §5e.** Their tables are kept as the only
pre-H7 measurement of the Rubin-astrometry and blending fixes, and are labelled provisional.
Nothing taken from them is a result.

**A production run was in flight from 2026-09-05 18:31 — see §5d.** It is the v2 run, and it
carries four fixes made that afternoon on top of E1a/H1/H2/H4: Rubin's astrometric error
renormalised from a mission average to a per-visit one, Roman's blending made Poisson so it can
be blended at all, `--start-index` so an interrupted scan can be continued, and the repaired
`DET_ANOMALY` run total. The regression fixture, which had not compiled since H1, was repaired
first so those changes could be checked at all.

**Every number in §4 is now superseded twice over** — once by E1a/H1/H2/H4 and again by the
Rubin astrometry and blending fixes. Do not quote §4 without saying so.
**Branch:** `joint-fisher-refactor` (never commit to `main`).
**Head at last update:** `738b54d` — "Report what each survey measures, not just how much the
joint fit adds" (Step F4).
**Step E1a is committed on top of that**; see §2, Phase E.

**Nothing about the existing results has moved.** E1a changes how the *next* run samples the
sky; with its new flag absent, the simulator behaves exactly as it did. `test5.dat` and every
figure made from it still stand.

**Two caveats on the numbers in §4, both now fixed in code but not yet in any run.**
1. *Satellite parallax.* Every `piE` number in this file comes from a run in which
   `lightcurve()` put both telescopes at the centre of the Earth, so it contains only the
   annual Earth-orbit signal and is a **lower bound** on the real pair. Step H1 fixes it
   (Deviations 27, 28).
2. *Roman's astrometric error.* Every `tetE` and lens-mass number comes from a run in which
   Roman's epochs were weighted by a **stale Rubin astrometric error from a different
   timestep** — worse than the documented `errlsstA` placeholder. Step H4 fixes it
   (Deviation 29.2).

Both land in the next production run. Until then, quote §4 with these caveats.

---

## 0. Read this first if you are a new session

This project has deviated far enough from `JOINT_FIT_REFACTOR_PLAN.md` that the plan alone
no longer tells you where things are. **This file is the entry point.** Read it, then read
whichever of the four documents below your task touches.

| Document | What it is | When to read it |
|---|---|---|
| `PROGRESS.md` (this file) | Where we are, what is done, what is next | **Always, first** |
| `JOINT_FIT_REFACTOR_PLAN.md` | The original roadmap and the binding working style (§0) | Before any code step. §0 is not optional |
| `DEVIATIONS.md` | Every place the implementation departed from the plan, and why | Before touching anything the plan describes — the plan is often out of date, this file is not |
| `OPEN_ITEMS.md` | Known problems deliberately **not** fixed | Before "fixing" anything you notice. It may already be a recorded decision |
| `PHASE_H_PLAN.md` | The roadmap for satellite parallax, the astrometric shift, and the collaborator draft of the whitepaper | Before starting any of those three. It supersedes the original plan's Step G1 |
| `ORIENTATION.md` | Struct-by-struct tour of the C++ (Step A1's output) | When you need to know what a legacy variable means |

**The plan's original text is deliberately never edited.** It records what was believed when
it was written, which is worth preserving. Corrections live in `DEVIATIONS.md`.

**The working agreement in `JOINT_FIT_REFACTOR_PLAN.md` §0 governs everything** and
supersedes generic development habits: one step at a time, teaching brief before any edit,
wait for approval, explain the physics and not just the code, define every acronym, flag
uncertainty rather than inventing an explanation, and never silently fix something noticed
in a later step. One `git commit` per step, named for the physical change.

---

## 1. Where we are in one paragraph

Phases A through F of the plan are done. The pipeline now reports Roman-only, Rubin-only and
joint detection **and** precision separately for every simulated event, with three
independent Fisher matrices per event and a per-event output table. A full production run
finished on 2026-08-30 against a corrected bulge lens population: 5,571,168 simulated
events, 74,812 joint detections, over 1,489 aggregated sightlines. All three Phase F
analysis products (F1 table, F2 gap-filling figure, F3 characterization map) have been built
and run against it, and **the gap-filling signal — the thesis's second novelty claim — is
visible in the data for the first time.** A fourth product, F2 repeated for `sigma_piE`, has
since settled the plan's long-`tE` annual-parallax prediction: it is not supported, because
Roman alone already measures the parallax of long events (Deviation 24). A fifth, F4, reports
the absolute forecast precision per survey rather than a ratio, and reaches Deviation 24's
conclusion independently: inside Roman's footprint the joint fit improves a typical lens mass
by 4% over Roman alone and by a factor of five over Rubin alone (Deviation 25). Phase E's sampling work is now half done: **Step E1a** stratifies the scan
toward Roman's footprint and gives every event the sky area its sightline stands for, which is
what turns all four sample-limited results into quotable ones — they are limited by the same
1,950 in-footprint events, not by four different things (Deviation 26). It is implemented and
verified but **not yet run**: choosing `--stride-roman` is a wall-clock decision and is the
next thing that needs a human. What remains after that is Step E1b (the `tE` stratification the
plan asks for, which may prove unnecessary — `OPEN_ITEMS.md`), Step E2, Phase G (physics
separation, validation against published Rubin-only numbers, and reconciling the whitepaper),
and the deferred GBTDS-footprint item.

---

## 2. What has been done, in order

Commits are on `joint-fisher-refactor`. Where a step deviated from the plan, the
`DEVIATIONS.md` entry number is given; read it before assuming the plan describes the code.

### Phases A–B — orientation and per-survey bookkeeping
- `ORIENTATION.md` written (Step A1); `OPEN_ITEMS.md` opened (Step A2).
- Per-survey detection bookkeeping added: every light-curve datum is tagged in `lens::tele[]`
  with which observatory produced it (0 = Rubin, 1 = Roman), and the detection chi-squared
  is accumulated separately for Rubin-only, Roman-only and joint streams.
- Roman astrometric errors still use `errlsstA()` as an explicit placeholder — **`OPEN_ITEMS.md`**.

### Phase C — the Fisher matrices
- `a27a6ce` — `fb` perturbation was binned off Rubin's blend fraction and applied to
  whichever telescope produced the epoch, pushing Roman's past 1.0 and crashing the program.
  Pre-existing. **Deviation 1.**
- `5f39836` — `t0` added to the photometric parameter vector, **appended at index 5** rather
  than inserted second as the plan said, to avoid reflowing six hard-coded index sites.
  **Deviation 2.** The code's order is authoritative:
  `0 u0, 1 tE, 2 fb0, 3 piE, 4 xi, 5 t0, 6 mbs0, 7 fb1, 8 mbs1`.
- `385a6c1` — each telescope gets its own free source flux and blend fraction (Step C2b).
- `d8411f3` — the Fisher matrix is normalized before inversion and its condition number
  recorded (`condA_*`, `condB_*` in the output).
- `1683480`, `bd1cb6b` — derivative steps placed on their convergence plateau (Step C3),
  verified by `./fishertest --sweep` through `tests/c3_step_sweep.py`.
- `eb192e8` — events labelled by joint-fit synergy (`SynergyClass`).
- `7419e17` — `tests/fisher_fixture.cpp` / `make fishertest`: a data-free Fisher regression
  harness that runs in ~3 s. **Run it before and after any change to `FisherM` and diff the
  output.**

### Phase C–D (a section inserted between C and D) — taxonomy and run statistics
- `e47390a` — `DetClass` replaced the old five-way label; the joint detection test was made
  monotone in the data. **Deviations 14.1–14.5.**
- `2492559` — blend flux was accumulating across every star ever drawn. Pre-existing.

### Phase D — the survey data, the schedule, and the output table
- `df873ad` — **the survey data had never actually reached the simulation.** Epoch counts
  were capped at 1000, and `BulgeBaseline.dat` was doubled with a silent short read.
  **Deviation 15.**
- `74e5e18` — Roman's mission placed on the real GBTDS schedule from STScI's published
  design: `MISSION_START_DAY = 730` as a runtime parameter, real season windows (10 seasons
  over days 730.000–2447.965), and `FoVRoman` corrected — it had been confusing an area with
  a radius. **Deviation 16.**
- `d4151fe` — the simulator had never been compiled with optimization. **Deviation 17.**
- `81a6b04` — the scan covers the whole survey region instead of one corner. **Deviation 18.**
- `30c9add`, `9fa4f88` — every event records where it peaked relative to Roman's seasons:
  `dt_edge` (days to the nearest season boundary, negative inside a season) and `t0zone`
  (in-season / in-gap / off-mission). A guard refuses a schedule whose "seasons" hold a
  single epoch. **Deviation 19.** `dt_edge` is the independent variable of the F2 figure.
- `aad8d55` — only the light-curve slots an event actually used are reset.

### The production runs and what went wrong on the way
- `682c978` — the widened scan stalled forever on sightlines neither telescope observes.
  Empty-sightline skip plus a `maxDraws` cap. **Deviation 20.**
- `d6dd293` — a `-1` sentinel was dragging a per-field mean precision negative.
  **Deviation 21.**
- `5c74fbd` — **the lens mass function was the LMC simulation's MACHO range** (3–5000 solar
  masses), giving a median lens of 386.9 solar masses and a median `tE` of 953 days. Replaced
  with a Kroupa (2001) IMF plus stellar remnants. **Deviation 22 — read this one.** Every
  number produced before this commit describes a MACHO population, not the bulge.

### Phase F — the analysis layer
- `df6f39b` — `analysis/romanlib.py` (the shared sentinel-aware loader), `analysis/f1_results_table.py`,
  `analysis/f2_gap_filling.py`. **Deviation 23.**
- `2526bd0` — `analysis/f3_characterization_map.py`.
- `4f3bb0d` — `f2_gap_filling.py` gained `--param {tE,piE}`, and the `sigma_piE`
  version of F2 was made. It answers the question Deviation 23.2 left open, and the answer
  is negative: the plan's long-`tE` annual-parallax prediction fails on `piE` too.
  **Deviation 24.**
- Step F4 — `analysis/f4_fisher_precision.py`, a fourth Phase F product the plan does not
  contain: the absolute forecast precision per survey partition, read straight off the three
  Fisher matrices instead of as a ratio. `romanlib.load_events()` gained chunked filtered
  reading in the same step, because the unfiltered read is OOM-killed on this machine and the
  kill is silent (exit 0, no output). **Deviation 25**, and an `OPEN_ITEMS.md` entry for
  F1/F2/F3, which still read unfiltered.

---

### Phase E — sampling strategy
- Step E1a — the scan is stratified in **sky position**. `--stride-roman N` visits sightlines
  inside Roman's GBTDS footprint on a finer grid than the rest; every sightline carries the
  deg² of sky it stands for, written as `w_area` into every event row and into the map file
  (which also, for the first time, records `lon`/`lat`, so an event can be tied back to the
  sightline that produced it). `--dry-run` builds the grid and reports the strata without
  drawing a star, so the cost of a choice can be read before a multi-hour run rather than
  during one. **Deviation 26.** The `tE` half of the plan's Step E1 is deliberately not done
  and is argued against in 26.1 and `OPEN_ITEMS.md`.

  **What it buys, at `--stride 10` (`--dry-run`, full region):**

  | `--stride-roman` | footprint step | footprint sightlines | total sightlines | footprint sample vs now |
  |---|---|---|---|---|
  | absent | 0.20° | 39 | 1,706 | 1× (this is the current run) |
  | 5 | 0.10° | 147 | 1,829 | ~3.8× |
  | 2 | 0.04° | 907 | 2,591 | ~23× |
  | 1 | 0.02° | 3,656 | 5,357 | ~94× |

  Footprint sightlines are the expensive ones (~50,000 Roman epochs against ~2,400 Rubin), so
  wall clock grows faster than the sightline count. **This choice is not made; it is the next
  decision.**

### Phase H — satellite parallax and the astrometric shift
- Step H1 — **Roman is at L2.** `lightcurve()` takes a `tele` argument and returns the
  trajectory that observatory actually sees; Roman's heliocentric position is Earth's scaled
  by `(1 + L2_OFFSET_AU)`, `L2_OFFSET_AU = 0.01003`. Threaded through all four call sites,
  including the three inside `FisherM` — a derivative evaluated with a different observer than
  its datum makes the matrix inconsistent. `--no-satellite-parallax` puts Roman back on Earth,
  which is the "off" run of Step H3's experiment and the step's own regression.
  **Deviation 28.** The `t = 0` gauge trap that would have deleted the whole signal is
  described there and in `PHASE_H_PLAN.md` H1.

  Verified: `|Δu| = piE · D_perp` to machine precision; the observers differ by 8.708e-03 in
  projected position at `t = 0` (non-zero is the whole point); exact coincidence when switched
  off. **`D_perp` is the PROJECTED separation and ran 0.87–0.99 of the full L2 offset across a
  year** — `L2_OFFSET_AU · piE` is the ceiling on the effect, not its value (Deviation 28.2).

- Step H2 — **the satellite observable is in the table.** `du_sat` (the observer separation in
  Einstein radii at `t0`, `piE · D_perp/AU`) plus `nepL_pk`/`nepR_pk` (epochs from each survey
  within ±2 `tE` of `t0`, i.e. while the event is magnified — satellite parallax needs
  *contemporaneous* coverage, which `ndw_L`/`ndw_R` cannot express). **Deviation 30.**

- Step H5 — **the astrometric shift is a product.** `analysis/h5_astrometric_shift.py`, six
  panels. The headline is a reconciliation: the centroid shift is far below a single exposure
  and only clears the noise after averaging, so *detectable* and *forecastable* are different
  questions with opposite answers. **Deviation 31.** The numbers first recorded here were
  pre-H7 and are superseded by §5f, measured on the v3 table: median shift **0.1106 mas**
  against **6.69 mas** per exposure, **0.11%** of events above a single exposure, **3.68 sigma**
  against the stack, **42.5%** with an astrometric peak across a season edge.

- Step H4 — **Roman has its own astrometric error.** `errRomanA()` in `helper.cpp`: a 1.1 mas
  centroiding floor (1% of the 110 mas pixel) for `F146 <= 20.62`, rising to 10 mas at 23.5 and
  at 0.4/mag beyond, from Sanderson et al. 2019 (arXiv:1712.05420) and arXiv:2608.24998.
  **Per exposure** — one row of `RomanBaseline.dat` is one 12.1-minute exposure (measured), so
  the 0.1 mas daily-binned figure would have overstated Roman's astrometry tenfold.
  **Deviation 29.** It also fixed a pre-existing bug: the Roman branch computed `errsR` and
  then stored `errs`, Rubin's error from a different timestep, into `l->erra[]`
  (**Deviation 29.2**).

## 3. The current data, and what is wrong with its label

| Item | Path |
|---|---|
| Per-event table (2.5 GB, 5,571,168 rows, 90 columns, `#` header) | `test5.dat` |
| Per-sightline aggregates | `files/MONTLMC/files/MapLMC5.dat` |
| Run provenance | `files/MONTLMC/files/run_provenance.txt` |
| F1 results table | `f1_kroupa.csv` |
| F2 gap-filling figure and its data, `sigma_tE` | `f2_kroupa.png`, `f2_kroupa.csv` |
| F2 gap-filling figure and its data, `sigma_piE` | `f2_piE_kroupa.png`, `f2_piE_kroupa.csv` |
| F3 characterization map and its data | `f3_kroupa.png`, `f3_kroupa.csv` |
| F4 Fisher precision figure, Roman footprint | `f4_fisher_kroupa.png`, `f4_fisher_kroupa.csv` |
| F4 Fisher precision figure, all joint detections | `f4_fisher_all_kroupa.png`, `f4_fisher_all_kroupa.csv` |
| Previous (MACHO-population) run, kept for comparison only | `runs/macho_final_20260830/` |

Run configuration: stride 10, `maxdraws` 5e4, `IMnum = 5` (Kroupa + remnants), seed 42.

**Sightline accounting:** 1,489 aggregated + 140 no-coverage + 77 barren = 1,706 exactly,
with 77 capped. 68.24 deg² scanned, 59.56 deg² produced events. Zero assertion failures.

**Detections:** 5,496,356 none (98.66%) · 73,489 Rubin+joint (1.32%) · 788 Roman+joint ·
535 both+joint · 0 joint-only · 0 ANOMALY.

**The provenance label is wrong and the data are fine.** `run_provenance.txt` says
`git_commit=d6dd293-dirty`; the sources that were actually compiled are the content of
`5c74fbd`. The binary was built after patching `Lensing.cpp` but before committing it, and
the later `make` said "Nothing to be done" and kept the stale stamp. Every figure carries
the wrong commit in its footer. **`OPEN_ITEMS.md` has the fix options and the interim rule:
`make clean && make` after the last commit, never before.**

The binary has since been rebuilt clean at `7a1b591` and stamps correctly, and `./fishertest`
passes with all assertions held. The **next** run will be labelled properly; `test5.dat` and the
figures already made from it keep the wrong label and should be cited with that caveat.

---

## 4. The results, as they currently stand

### F2 — the gap-filling figure. This is the headline.

1,363 events in scope (joint-detected, peaking within Roman's mission, and **inside Roman's
footprint** — that last restriction is load-bearing, see Deviation 23.1). Ratio defined on
1,341; 84 rescued; characterized joint 709 against Roman alone 625.

Median `sigma_joint / sigma_Roman` for `tE`, against days from the nearest season edge:

| tE bin | mid-season | deep in gap |
|---|---|---|
| 10–30 d | 0.99 | **0.014** |
| 30–100 d | 0.97 | 0.49 |
| 100–300 d | 0.92 | 0.90 |
| 300+ d | ~1.0 | ~1.0 |

The yield panel shows the rescue fraction rising through the gap, reaching ~0.30 for the
30–100 d bin.

**The ordering is the reverse of what the plan predicted.** See Deviation 23.2.

### F2 in `sigma_piE` — the plan's long-`tE` claim, tested and not supported

`f2_piE_kroupa.png`, same 1,363 events, same binning, ratio in `piE` instead of `tE`:

| tE bin | mid-season | in gap (+22.5 d) | deep in gap (+52.5 d) | n |
|---|---|---|---|---|
| 10–30 d | 0.99 | 0.29 | **0.056** | 329 |
| 30–100 d | 0.97 | 0.82 | 0.55 | 476 |
| 100–300 d | 0.94 | 0.93 | 0.90 | 267 |
| 300+ d | 0.98 | 0.97 | **0.99** | 124 |

The ordering is the same as for `sigma_tE` — the drop deepens with **short** `tE`, and the
long-`tE` bin is flat at ~1 straight through the gap. The plan expected the opposite, on the
grounds that Roman's gaps would prevent the annual parallax from being sampled for long
events. They do not: a 300+ day event spans several Roman seasons, and **94% of long
in-scope events have `piE` measured to better than 2 sigma by Roman alone** (median
`sigma_piE` = 0.0098 on a median `piE` of 0.21). There is nothing left for Rubin to add.

**The governing variable is whether Roman saw the event at all, not which parameter is being
forecast.** Full argument, diagnostics and the second-order result (`piE` gains ~2.4× less
from gap-filling than `tE` does) in **Deviation 24**. The long-`tE` null rests on 124 events
— see `OPEN_ITEMS.md` before quoting it.

### F1 — the per-field results table

24 rows: all 6 Roman field centres × 4 `tE` bins. Totals across the table: 962 events
characterized by the joint fit against 717 by Roman alone, **+245 characterized events from
combining the two.** The gain is present in every field and grows with `tE` bin as a
fraction, e.g. field F3 in the 30–100 d bin: joint 69, Roman 38, ΔN = +31.

The fraction of gap-peaking events seen by Rubin rises with `tE` from ~0.04 (10–30 d) to
0.21–0.56 (300+ d), which is the yield statistic behaving as expected.

Rubin-alone Fisher singularity is essentially zero (≤1.4% in two cells, 0 elsewhere).

### F4 — what the Fisher matrices forecast, before any ratio is taken

`f4_fisher_kroupa.png`. F1–F3 all report how much the joint fit *adds*; this one reports
what each survey partition actually *measures*, which is the scale those ratios are ratios
of. Six panels: cumulative distributions of the fractional 1σ forecast on `tE`, `piE`
(photometric matrix) and `tetE` (astrometric matrix), the derived lens mass
`Ml = tetE / (kappa * piE)`, the per-event joint-against-single scatter that exhibits the
`sigma_joint <= sigma_single` invariant, and the condition-number distribution.

Fraction of the 1,950 in-footprint joint-detected events measured to better than 10%:

| Parameter | joint | Roman alone | Rubin alone |
|---|---|---|---|
| `tE` | **34.8%** | 24.3% | 11.8% |
| `piE` | **26.9%** | 20.4% | 9.3% |
| `tetE` | **90.1%** | 88.5% | 55.6% |
| `Ml` (derived) | **32.2%** | 27.2% | 7.5% |

Median per-event `sigma_joint / sigma_single` on the lens mass: **0.96 against Roman alone,
0.19 against Rubin alone.** Zero points above the 1:1 line.

**This is Deviation 24's conclusion reached a second way, from precision alone and with no
reference to season geometry.** Inside Roman's footprint, adding Rubin to Roman buys ~4% on
a typical mass; adding Roman to Rubin buys a factor of five. Rubin's contribution is
concentrated in events Roman never saw, not spread over the ones it did — gap-filling is a
*yield* effect, not a precision effect.

`f4_fisher_all_kroupa.png` is the same figure over all 74,812 joint detections and is a
**check, not a science figure**: the joint and Rubin curves lie exactly on top of each other
and the mass ratio is exactly 1.000, because Roman observed only 2.6% of the sample and on
the other 97.4% the joint matrix *is* Rubin's matrix. Any departure from 1.000 there would
be a partitioning bug.

**Caveat, inherited not introduced:** the `tetE` and `Ml` panels rest on the astrometric
matrix, and Roman's per-epoch astrometric error is still `errlsstA()` as a placeholder
(`OPEN_ITEMS.md`). The ordering of the curves is robust; the absolute fractions in those two
panels are only as good as that placeholder.

### F3 — the (`tE`, `piE`) characterization map

Two panels on one shared colour scale, in the format of Abrams et al. 2025 Figures 11–14.
Panel (a), joint over Roman-alone inside the footprint: 1,950 events, 16 coloured cells,
ratios up to 2.0, concentrated at `tE` of 10–100 d. Panel (b), joint over Rubin-alone over
all 74,812 joint-detected events: 25 coloured cells, max ratio 1.48, plus two hatched cells
where Rubin alone characterizes nothing and the joint fit does.

**Panel (a) rests on only 1,950 events** because so few sightlines fall inside the GBTDS
footprint. Do not quote its cell values as precise — `OPEN_ITEMS.md`.

---

## 5. What is next

**The user has asked for three things, in this order: satellite parallax, the astrometric
shift, then the whitepaper brought up to date for potential collaborators. `PHASE_H_PLAN.md` is
the roadmap for all three** — read it before starting any of them; it also explains why the
original plan's Step G1 cannot be run as written.

0. **The production run — v2 IN FLIGHT since 2026-09-05 18:31 (§5d). The first attempt that
   day was stopped at 526/1,829 sightlines (§5c).** E1a (stratified sampling), H1
   (satellite parallax), H2 (the satellite observable) and H4 (Roman astrometric errors) are
   all in the code and all change what a run produces; nothing downstream could move until one
   run existed with all four. The user chose `--stride-roman 5`. The run was stopped by the
   user after the laptop hibernated mid-run and the remaining cost was measured at ~2.5 h.
   See §5c for what it produced and what it cost. **Every result in §4 still predates all four
   changes and stays caveated.**

   Both blockers raised after the first attempt are now cleared. The `DET_ANOMALY` jump was
   **not** a regression and **not** H1 — the run-total counter was simply never incremented, so
   every run had reported zero (`OPEN_ITEMS.md`, and Deviation 32). And the run is now resumable
   via `--start-index`, so an interruption no longer costs the whole scan.

   **The `--no-satellite-parallax` twin that Step H3 needs is scheduled, not forgotten.** H3 is
   the difference between the two runs. A watcher checks memory at the primary's midpoint and
   either starts the twin alongside it or holds it until the primary finishes; §5d has the rule
   it applies and where the twin's output goes.

Then, in roughly the order that makes sense:

1. **Run the stratified scan.** Step E1a is built, verified and committed; what it needs is
   a `--stride-roman` and a machine. Use `./roman --dry-run --stride-roman N` to see the
   sightline counts first (§2, Phase E). `--stride-roman 2` is the one that makes all four
   sample-limited results quotable (~23× the in-footprint sample), and it is also the
   expensive one; `5` is the cheap version at ~3.8×. Re-run F1–F4 against the new table
   afterwards, and read the two `OPEN_ITEMS.md` entries about area weighting **before**
   quoting anything pooled across the whole sky from it.

1b. **Step E1b — stratify in `tE`.** The other half of the plan's Step E1, deliberately not
   built. Deviation 26.1 argues it may be unnecessary once E1a has run, because the thin bins
   are footprint bins rather than rare-`tE` bins. Decide after seeing the stratified table.
2. **Step E2 — revisit the per-sightline stopping criteria.** The current third floor (2
   well-conditioned events) was written when there was one Fisher matrix; with three
   matrices and stratified bins it needs replacing.
3. **Phase H is complete.** H1, H2, H4, H5, H7, H3 (§5g) and H6 (§5h) are all done. **H3
   replaced the original plan's Step G1**, which cannot be run as written (Deviation 27), and
   **G3 was absorbed into H6**. What Phase H leaves open is in `OPEN_ITEMS.md`, and the item
   that matters most is the exposure-independence assumption behind every `theta_E` number.

4. **The rest of Phase G.** G2: validate the Rubin-alone branch against Abrams et al. 2025 at
   l = 0.33°, b = 2.82°, **reweighting to their sampling first** or the comparison will look
   like a bug. G3 is absorbed into H6.
5. **The GBTDS footprint item** — deferred by the user, but wanted. The sky coverage is
   Penny et al.'s and the footprint has since changed; Rubin's large field of view can see
   corners of the GBTDS region whose centre it is not pointed at, which matters for blending
   too, and `BulgeBaseline.dat` must change with it. **`OPEN_ITEMS.md` has the full text.**

---

## 5b. Where the code currently lives, and how to re-verify it

**Step E1a is committed on a branch that has NOT been merged.** It is
`worktree-e1-stratified-sampling` (commit `5e805d1`, on top of `738b54d`), pushed to `origin`.
It was developed in a git worktree so that the user's checkout was never touched. To review:

```bash
git diff joint-fisher-refactor..worktree-e1-stratified-sampling
```

**The regression recipe, which every step that touches the simulator should repeat.** Step E1a
established it and Step H1's acceptance criteria reuse it: build a binary from the previous
commit, run both on the same tiny configuration, and diff every output file.

```bash
# baseline binary from the previous commit, built somewhere outside the tree
git show <prev-commit>:Bulge_LSST.cpp > /tmp/base/Bulge_LSST.cpp   # and Bulge.h
g++ -O2 -std=c++17 -DGIT_COMMIT='"base"' -I/tmp/base -o /tmp/base/roman_base \
    /tmp/base/Bulge_LSST.cpp Lensing.cpp helper.cpp -lgsl -lgslcblas -lm

# same tiny run for each; ~5 min per binary, 9 sightlines, ~48 event rows
./roman --stub --stride 2 --events 2 --lenses 1 --nerr 0 --maxdraws 500
```

Clear `test5.dat` and the append-mode files in `files/MONTLMC/files/` between runs or the two
runs concatenate (`OPEN_ITEMS.md`). Compare `test5.dat` on the shared column prefix, and
`MapLMC5.dat` / `LpLMC5.dat` / `EfLMC5*.dat` byte for byte.

**Working in a worktree:** the data files are gitignored and therefore absent, so symlink the
inputs in and keep the outputs private — never symlink `files/MONTLMC/files/`, because those
are opened in append mode and a test run would concatenate itself onto the production outputs.
`Baseline/BulgeBaseline.dat`, `Baseline/RomanBaseline.dat`, `CMD/components/*.dat`,
`files/density/*`, `files/ext/*`, `files/sigmaA_LSST.txt`, `files/sigma_roman.txt` are the
inputs; `files/MONTLMC/files/{LpLMC,EfLMC,EfLMC*B,MapLMC}<IMnum>.dat` must exist (they can be
empty) or the run exits with "Cannot open one or more files!".

## 5c. The run that was attempted on 2026-09-05, and what it cost

**Outcome: stopped at 526 of 1,829 sightlines, no usable production table.** Launched 01:54,
stopped 09:29 on the user's instruction. The laptop hibernated for most of that window, so the
run consumed only **3,080 s (51 min) of CPU** across 7 h 35 min of wall clock. With ~2.5 h of
CPU still to go and the expensive footprint band immediately ahead, it was not worth resuming
in place.

**The partial table no longer exists.** It was truncated to its header on 2026-09-05 at 18:07
by a `./roman --dry-run` run from the worktree while the worktree symlinks still pointed here:
`--dry-run` opens the event table in truncating mode and writes the header *before* it decides
not to simulate (`OPEN_ITEMS.md`). `run.log` and `README.md` survive, so everything below that
was extracted from the run is intact, and v2 supersedes it in any case. What follows describes
what it held.

It was **not** a production table: the sky coverage stops around `lon = -0.2` and is systematically
one-sided, and it reached only **22 of 147 footprint sightlines**, fewer than the finished
2026-08-30 run's 39. Its use is as a pipeline artefact — it carries the E1a, H1, H2 and H4
columns, so H3 and the re-run F-series can be developed against it without waiting for a full
run. Its final `MapLMC5.dat` line is a partially flushed buffer, so drop the last line when
parsing.

**What it produced that is worth keeping.**

1. **A measured cost model.** Outside the footprint a sightline costs **~2.2 s**; inside it,
   **~18-40 s** — roughly **8-18x**. The reason is in the log: a footprint sightline gives every
   event **50,401 Roman epochs** against Rubin's ~2,300, and H1 calls `lightcurve()` twice per
   Roman epoch. **A full `--stride-roman 5` run costs about 3.5 h of uninterrupted CPU.**
   Note that both naive extrapolations fail: the scan opens on the barren western edge where
   every sightline runs to the 50,000-draw cap, so the first minutes are the slowest and most
   row-producing and overestimate the total (the first 24 sightlines took 15.1 s each); the
   middle of the scan is all cheap outside sightlines and underestimates it (2.2 s each).
2. **The `DET_ANOMALY` jump**, which is the more important of the two and has its own
   `OPEN_ITEMS.md` entry. Detection classes over the partial run: 393,331 none · 26,683
   Rubin+joint · 337 Roman+joint · 281 both+joint · 0 joint-only · **433 DET_ANOMALY**.

**It is not resumable, and that is now the binding constraint.** There is no checkpoint and no
`--start-index`; the sightline vector is rebuilt from scratch every launch and the outputs are
opened in append mode. On a laptop that hibernates, a 3.5 h single-shot run is a poor bet. The
cheap fix is a `--start-index N` that skips the first N entries of the already-deterministic
`scan` vector; the append-mode outputs already concatenate correctly. **Built as `91d4200`.**
Note the resume index is the number of sightlines the scan **entered**,
`grep -c 'NEW STEP'` on the interrupted run's log — **not** the `MapLMC5.dat` line count, which
this paragraph originally proposed and which is wrong. A no-coverage or barren sightline is
entered and counted but writes no map row, so the map file undercounts: on the v2 run it read
624 against 710 entered. Resuming at the map count would re-simulate those 86 sightlines and
append duplicate rows to `test5.dat`.

Launched 2026-09-05 01:54 local, from the worktree, at commit `e8f4135`:

```bash
./roman --events 300 --lenses 50 --stride-roman 5
```

**`--events 300 --lenses 50` are not the defaults and are not decoration.** The built-in
defaults are 850 and 150. The 2026-08-30 production run passed 300/50 explicitly, so taking
the defaults here would have changed the per-sightline statistical budget by ~2.8x at the same
moment as the stratification and the new physics, and the two tables would not have been
comparable. The only intended differences from the 2026-08-30 run are the ones that were
chosen: stratified sky sampling, satellite parallax, Roman's own astrometric errors, and the
three H2 columns.

**Grid, from `--dry-run`:** 147 footprint sightlines at 0.1 deg (1.47 deg^2) and 1,682 outside
sightlines at 0.2 deg (66.47 deg^2); 1,829 total over 67.94 deg^2. The footprint is 8.0% of
the sightlines and 2.2% of the area — the imbalance that `w_area` exists to undo. Against the
2026-08-30 run's 39 footprint sightlines this is **3.77x the in-footprint sample**, which is
the sample that limits all four of the results in §4.

**The output does not live in the worktree.** It is multi-GB and gitignored, and a worktree can
be deleted with the session that made it, so `test5.dat` and the six append-mode files in
`files/MONTLMC/files/` are symlinks into

```
/home/ali/Documents/PhD/Offline_project/roman_runs/2026-09-05_stride-roman-5/
```

which also holds `run.log`, `run.pid` and `run_provenance.txt`. This is a stronger guard than
clearing the files would have been: the run physically cannot append onto the 2026-08-30
outputs, which stay untouched at the repo root. **If you symlink outputs for a run of your own,
give it its own directory** — the standing warning against symlinking `files/MONTLMC/files/`
is about pointing it at outputs that already matter, not about the technique.

**Checking on it:** `wc -l` on `MONTLMC/MapLMC5.dat` in that directory counts finished
sightlines out of 1,829. `run.log` is block-buffered, so an empty log is not a failed run —
that was diagnosed once already and is not worth diagnosing twice.

**This run's provenance stamp reads `6c97375-dirty`. Do not believe it, and do not distrust the
run because of it.** The binary was built at 00:22 on 2026-09-05, while H2 was written but not
yet committed; H2 was then committed as `1c4f04e` and H5 as `e8f4135`, and neither commit
touched a file `make` watches, so `make` reported "Nothing to be done" and the binary kept the
stamp from its dirty build. The `-dirty` suffix is doing its job — it truthfully says "this
commit plus edits" — it just cannot say which commit those edits became.

**The gap was closed by proof rather than by restarting.** `Bulge_LSST.cpp` was recompiled from
the clean `e8f4135` tree with identical flags and the running binary's stamp, and the two object
files were compared after `objcopy --strip-debug`: the disassembly of `.text` is identical apart
from the filename in objdump's own header line. The running binary therefore contains exactly
`e8f4135`'s source, and this run's table may be attributed to `e8f4135`. The output confirms it
independently — the header carries `du_sat nepL_pk nepR_pk` (H2), and the provenance block
carries `stratified 1`, `stride_roman 5` (E1a) and `satellite_parallax 1`,
`L2_offset_AU 0.0100267` (H1).

The reason this keeps happening, and the build-system change that would stop it happening a
third time, is an `OPEN_ITEMS.md` entry of its own. **Until that lands, `make clean && make`
after the last commit and before any run is not optional** — the interim rule that commit
`585432b` recorded the first time this went wrong.

**Early cost, for whoever plans the next one.** The scan starts at `lon = -3.719 deg`, the
western edge, well outside Roman's footprint. Those sightlines are barren, run to the
`maxdraws` cap of 50,000 draws each, and are therefore both the slowest and the most
row-producing in the whole scan: the first 24 took 15.1 s each and wrote 317 MB. Do not
extrapolate a total from the first few minutes of a scan in this geometry — you will
overestimate. The 2026-08-30 run capped on 77 of 1,706 sightlines and finished in ~3.5 h at
2.5 GB.

---

## 5d. The v2 run — KILLED 2026-09-06 for Step H7, superseded by §5e

Launched 2026-09-05 18:31 local at commit `a5bb600`, built with `make clean && make` so the
stamp is current and clean (no `-dirty`) — the rule §5c had to be written about twice:

```bash
./roman --events 300 --lenses 50 --stride-roman 5
```

Same targets as the 2026-08-30 run and the first 2026-09-05 attempt, for the same reason (§5c):
the defaults are 850/150 and changing the statistical budget at the same moment as the physics
would leave two tables that cannot be compared.

**Output:** `/home/ali/Documents/PhD/Offline_project/roman_runs/2026-09-05_v2_stride-roman-5/`,
outside the repository, symlinked in — same arrangement and same reasons as §5c.
**Twin:** `.../2026-09-05_v2_twin_nosat/`, a fully isolated directory with its own input
symlinks so it can run concurrently without sharing the worktree's output paths. The twin is
`--no-satellite-parallax` with everything else identical, so the two runs draw identical events
and Step H3's comparison is paired. The scheduling rule: at the primary's midpoint, start the
twin alongside if at least 1.8 GB is available (an ~800 MB run plus 1 GB of headroom), and
otherwise hold it until the primary exits rather than risk both.

**What this run carries that no previous run did:**

| | |
|---|---|
| E1a | stratified sky sampling, `w_area` on every row |
| H1 | satellite parallax — Roman at L2, not at the centre of the Earth |
| H2 | `du_sat`, `nepL_pk`, `nepR_pk` columns |
| H4 | Roman's own astrometric error model, and the stale-`errs` fix |
| `d393dc3` | Rubin astrometry renormalised per-visit — **26.7x worse, and correct** |
| `58d2863` | Roman blending Poisson — **14.6% of Roman events now blended, was 0%** |
| `a5bb600` | `DET_ANOMALY` run total actually counted |

**The cost model, measured three times and wrong the first two.** Runtime per sightline is
not one number and not even two. Measured on the primary at 65/147 footprint and 655/1682
outside sightlines, 3.55 h of CPU in:

| | measured | how |
|---|---|---|
| outside sightline | ~13 s | joint solve of the primary and twin totals |
| footprint sightline, run average | ~50-70 s | same solve |
| footprint sightline, **in the plane** | **~360 s** | direct 20-min window, `nri` 40, lat -0.14 |

Three traps, in the order they were fallen into. **The opening of the scan is not
representative**: the first columns are empty sky that skips instantly, which is where an early
"~2.2 s per outside sightline" came from; it does not survive contact with Rubin's coverage,
where an outside sightline still carries ~2,300 epochs. **The overall average is not
representative either**, because it is dominated by the ~92% of sightlines that are outside.
And **the footprint stratum is itself strongly inhomogeneous**: cost rises by a factor of ~5-7
as the scan crosses the Galactic plane, where the source density and optical depth peak, so
more draws and more Fisher matrices are needed per sightline. The footprint columns run `nri`
33-40 so far (lon -0.619 to +0.081), widening from 4 to 12 sightlines per column as they
approach lat 0.

The practical consequence: **any ETA built from a single seconds-per-sightline number is wrong,
and wrong by a factor of several in whichever direction the current stratum differs from the
average.** Progress must be split by stratum and re-measured as the scan moves. `watch_v3.sh`
does this; the earlier watchers did not, and additionally counted `MapLMC5.dat` rows, which
stall through skipped sightlines.

**First result out of the run, and it is about H7 rather than H3.** The restored
`DET_ANOMALY` counter reports for the first time, and the pathology is not a corner case.
Summed over the primary's first 637 scored sightlines:

| stratum | sightlines | detections | DET_ANOMALY |
|---|---|---|---|
| outside the footprint | 574 | 30,144 (all Rubin+joint) | **0** |
| inside the footprint | 63 | 3,154 | **854 (21.3%)** |

Outside, `ANOMALY` is zero everywhere and every sightline reports Rubin+joint 50/50, which is
the correct degenerate case: with no Roman epochs the joint test *is* Rubin's test. Inside,
**about one detection in five is found by a single telescope and then rejected by the joint
test** — the `2*ndw` threshold pathology of Step H7, now measured on production data instead of
argued from the source. Roman's 50,401 epochs put its own bar at 100,802; adding Rubin's ~2,300
near-flat epochs lifts the joint bar to 105,530, so the joint test is *harder to pass than
either telescope alone*, which is backwards.

No table is wrong because of it — `detJ` is forced monotone (`if (detL or detR) detJ = 1;`) —
but this is the number that justifies H7 blocking H3, H5 and H6, and it is large enough that
H3/H5 results taken off these runs must be labelled provisional.

**Expect the numbers to move a long way, and mostly against Rubin.** Rubin's astrometric
information falls by ~715x, so any `tetE` or lens mass that leaned on Rubin astrometry gets
worse. Roman's detection efficiency should fall too, because `testR <= blend[6]` was accepting
everything while `blend[6]` was pinned at 1. Both changes make the forecast more conservative
and both were bugs, so a smaller joint gain in the next F-series is the expected outcome, not a
regression to hunt.

---

## 5e. The v3 runs, post-H7 — in flight

Launched 2026-09-06 00:55 local at commit `f959c8c`, `make clean && make`, stamp clean.

```bash
./roman --events 300 --lenses 50 --stride-roman 5                            # primary
./roman --events 300 --lenses 50 --stride-roman 5 --no-satellite-parallax    # twin
```

Same budget and same grid as v2, deliberately, so the only difference between the v2 and v3
tables is Step H7. 1,829 sightlines: 147 inside Roman's footprint at 0.1 deg, 1,682 outside at
0.2 deg.

**Both runs are in fully isolated directories**, each with its own input symlinks and **its own
copy of the binary**:

- `/home/ali/Documents/PhD/Offline_project/roman_runs/2026-09-06_v3_h7/`
- `/home/ali/Documents/PhD/Offline_project/roman_runs/2026-09-06_v3_h7_twin_nosat/`

This is a change from v2, where the primary was reached through the worktree's output symlinks.
Nothing in the worktree points at either run, so a `./roman` invoked there — `--dry-run`
included — cannot truncate them. That is how the v1 table was destroyed.

**What v3 carries that v2 did not:** Step H7 (`f959c8c`) — the detection bar is a fixed
`dchi > 500` on the signed statistic instead of `2*ndw`, applied identically to the Rubin-only,
Roman-only and joint tests. `DET_ANOMALY` is zero by construction; if the run reports any, the
construction is broken. Deviations 33.

**Expect detections to rise and every pre-H7 yield to be superseded.** On the paired test the
detection rate went 10.1% -> 18.7% inside the footprint and 8.5% -> 13.4% outside.

**The runtime prediction made here was wrong, and is corrected in place.** It said runtime per
sightline should *fall*, because the loop stops on its targets and detections now arrive about
twice as fast. That reasoning fails because the three targets are ANDed: the loop runs until
`icon >= --events` **and** `nlens >= --lenses` **and** `nerr >= --nerr`. `icon` counts events
simulated and H7 does not touch it, so `--events 300` stays the binding constraint, the number
of draws is unchanged, and the faster-arriving detections buy nothing — they only add Fisher
calls. Outside the footprint the cost still fell, **14.3 s -> 6.8 s** per sightline over 556 of
them, because there are few detections either way. Inside it the cost rose.

**How much it rose was also overstated, and is corrected here (2026-09-07).** This section
previously said ~52 s -> ~317 s, "six times the wall clock", and that the mechanism was
unknown. Profiling settled both. The mechanism is that `FisherM` runs only under
`if (detL or detR or detJ)`, so doubling the detection rate roughly doubles the Fisher work and
nothing else changes: on a paired profile of the same footprint sightline the detection-gated
path went 21.3 s -> 48.6 s (2.28x) against 2.57x more characterised detections, while
everything else went 1.11x. **Cost per Fisher call is unchanged.** The per-sightline factor is
**~1.4-2.0x, not 6x**; the old figure came from counting `NEW STEP` lines, which mix footprint
and outside sightlines. `OPEN_ITEMS.md` has the measurement and the arithmetic error in full.

So v3 is modestly slower overall, not dramatically: with the corrected footprint cost, v3's
savings outside the footprint offset most of what it spends inside, and the projected totals
are roughly **14 h for v2 against ~16 h for v3** — not 8.8 h against 16 h. The extra time buys
the detections H7 exists to recover, and per *detected* event H7 is slightly cheaper.

The comparison is at identical sky positions (same deterministic grid; the first twelve
footprint sightlines were checked to match on `nri`, lon and lat).

**Run in chunks.** At the user's instruction these run until told to pause, then are checkpointed
and resumed. The resume index is `grep -c 'NEW STEP' run<k>.log` for the chunk, plus the
`--start-index` that chunk began at — **not** the `MapLMC5.dat` row count.

**Each chunk writes its own log** (`run.log`, `run2.log`, ...). This is not tidiness: if every
chunk appended to one log, `grep -c 'NEW STEP'` would count the whole history and the resume
index would be wrong by the sum of all previous chunks. `resume_v3.sh` picks the next free name.

### Chunk 1: 2026-09-06 00:55 -> 09:50, paused at scan index 774

| | primary | twin |
|---|---|---|
| sightlines entered | 775 | 775 |
| footprint done | 78/147 | 78/147 |
| outside done | 689/1682 | 689/1682 |
| CPU | 8.88 h | 8.88 h |
| table rows (after trim) | 1,936,653 | 1,938,264 |
| `DET_ANOMALY` | **0** | **0** |

**Resume both with `--start-index 774`.**

**The last sightline was incomplete and its rows were removed.** Both runs were stopped
mid-sightline at scan index 774 (lon 0.281, lat -0.04), which had written 314 rows to the
primary's table and 118 to the twin's, plus 58 and 29 characterised events to `LpLMC5.dat`.
Those events are real but the sightline never reached its `--events`/`--lenses` budget, so its
event count does not represent the sky area it stands for and its `w_area`-weighted yield would
be an undercount. Leaving them and resuming at 775 would orphan them in the table; leaving them
and resuming at 774 would append a second, complete copy on top. Both files were therefore
truncated back to the last complete sightline (lon 0.281, lat -0.14) and the run resumes AT 774
so that sightline is done properly.

A trap for whoever automates this: `test5.dat` writes lon/lat as `0.281 -0.04` and
`LpLMC5.dat` writes the same values as `0.2810000 -0.0400000`. A string comparison finds the
rows in one file and silently misses them in the other. Compare numerically.

### Between chunks 1 and 2: the H7 cost regression was profiled and resolved

`perf` was used on one footprint sightline under each binary with a fixed draw budget, while no
production run was in flight. Result above and in `OPEN_ITEMS.md`: the extra cost is entirely
`FisherM` and its callees running more often, cost per call unchanged, and the regression is
~1.4-2.0x rather than the 2.5-6.6x previously recorded.

Two things worth keeping from how it was measured. **`perf`'s extrapolated cycle counts are not
comparable across runs** that experienced different machine load, because frequency-mode
sampling adapts the period and was throttled to 808-924 Hz; only within-run sample shares and
wall spans were used. And **startup is ~170 s of CPU** (~43% of it parsing the ~1000 extinction
files), which a one-sightline profile is mostly made of — it has to be subtracted before any
per-sightline number means anything. It was measured directly by running `--start-index 5000`,
past the end of the grid, so no sightline is entered.

### Chunk 2: resumed 2026-09-07 00:35 at scan index 774 -- FINISHED, and it exposed a data-loss bug

Chunk 2 entered 1055 sightlines and completed the grid (774 + 1055 = 1829). Both runs ended
cleanly with full summaries, `ANOMALY(single-not-joint) = 0`, 920 aggregated, 76 no-coverage,
59 barren, 59 capped on `--maxdraws`. Footprint accounting closes exactly: 82 + 66 - 1 (index
774 redone) = 147/147, outside 693 + 989 = 1682/1682.

**But the resume had silently truncated the event table.** `test5.dat` came out holding only
chunk 2 -- 4,138,393 rows, first data row at lon 0.281, lat -0.04, which is index 774. Chunk 1's
1,936,653 rows were destroyed at the moment chunk 2 started, by an `ofstream` opened in the
default (truncating) mode purely to write the column header. `EfLMC5.dat` and `EfLMC5B.dat` went
the same way. `MapLMC5.dat` (1606 rows) and `LpLMC5.dat` (82,888 rows) survived because they
were already opened `ios::app`. Fixed, with the rule and the verification, in Deviations 34.

**How it got missed:** the chunk-1 pause was verified thoroughly -- row counts, column counts,
truncation offsets -- and nothing was re-checked after the *resume*, which is where the loss
happened. A checkpoint is not verified until the run that follows it has been shown to still
contain the run that preceded it.

### Chunk 1 REDONE (part one), 2026-09-10 -- RECOVERED AND VERIFIED

Indices 0-773 re-run for both primary and twin into fresh directories (`2026-09-10_v3_part1`,
`..._twin_nosat`), then spliced in front of the surviving chunk-2 table. Fresh directories
deliberately: running with `--start-index 0` inside the v3 directories would now correctly, and
loudly, truncate them and destroy chunk 2.

**The recovery is exact, and that was proved rather than assumed.** The re-run's `MapLMC5.dat`
rows came out **686/686 byte-identical** to the chunk-1 rows that survived in the v3 directory,
for both runs. Independently, the computed trim offsets came out at **922,213,653** (primary)
and **905,619,093** (twin) bytes -- exactly the offsets the 2026-09-06 checkpoint truncated to.
The re-run reproduced chunk 1 down to the byte, as a fixed-seed `mt19937_64` advanced as a
single stream must.

| | primary | twin |
|---|---|---|
| combined `test5.dat` rows | **6,075,045** | **6,072,298** |
| header lines | 1 | 1 |
| first data row | lon -3.719, lat -4.14 | same |
| last data row | lon 5.081, lat 1.26 | same |
| `EfLMC5(B).dat` | 69,286 + 92,920 = 162,206 = 101 x 1606 | same |

Subtracting chunk 2 leaves 1,936,652 and 1,938,263 recovered data rows, matching the originals.
`MapLMC5.dat` (1606 rows) and `LpLMC5.dat` were never damaged and were left untouched. The
chunk-2-only table is kept beside the combined one as `test5_chunk2only.dat`.

**The v3 tables are complete: all 1829 sightlines, 147/147 footprint, 1682/1682 outside.**

### Step H3's paired run, and its provenance label

`roman_runs/2026-09-11_h3_paired/`, launched 2026-09-11 with `--pair-satellite --events 60
--lenses 20 --nerr 2 --stride-roman 5`. A reduced per-sightline budget on purpose: H3 needs
per-event ratios, not area-weighted yields, so how completely each sightline is sampled does
not enter the statistic.

**Its provenance stamp says `4570e53-dirty`, and the sources are the content of `85aa501`.**
The binary was built after patching in `--pair-satellite` and before committing it -- exactly
the trap recorded in section 3 of this file, walked into a second time. The `-dirty` suffix
means the stamp is at least honest: it says the sources were not the named commit. The mapping
is verified rather than asserted -- `git diff 85aa501 -- Bulge_LSST.cpp Lensing.cpp helper.cpp
Bulge.h` is empty, so the binary's sources are exactly that commit's.

The run was NOT restarted for this, because the data are unaffected and 520 sightlines were
already done. The worktree binary has since been rebuilt clean (`make clean && make` after the
last commit) and now stamps without a suffix. Any H3 figure must carry `85aa501`, not the
string in the file.

### A throttled laptop, not a cost model, explains the runtime numbers

Part one redid **provably identical** work -- byte-identical outputs -- in **3 h 08 m** of wall
clock, against the 8.88 h of CPU that chunk 1 was recorded as taking for the same 775
sightlines. Identical work cannot take 2.8x the CPU-seconds unless the clock itself was slower,
so the overnight run was throttled (thermal or power-save).

This retires a puzzle recorded earlier in this file: chunk 1 and chunk 2 could not be reconciled
into a single two-stratum cost model -- solving them together gave a *negative* footprint cost.
The reason was not that footprint sightlines vary by 3x with position. It is that the two chunks
ran at different clock speeds. **Per-sightline costs measured in different sessions on this
laptop are not comparable**, and any ETA built from them has to be bracketed accordingly.

### Chunk 2 (original launch record): resumed 2026-09-07 00:35 at scan index 774

Both runs relaunched with `resume_v3.sh 774`, logging to `run2.log`. Next resume index is
`774 + $(grep -c 'NEW STEP' <dir>/run2.log)`. About 8.6 h of work remained at the start of this
chunk (69 footprint and 993 outside sightlines).

### The v2 tables, kept and labelled

`2026-09-05_v2_stride-roman-5/` (734 sightlines, 806 MB) and `2026-09-05_v2_twin_nosat/` (618
sightlines) are kept. They are the only measurement of the Rubin-astrometry and blending fixes
uncontaminated by H7, and they carry the 21.3% anomaly measurement that justified H7. **They are
not results**: every detection in them was decided by the threshold H7 replaced.

---

## 5g. Step H3's satellite-parallax result, measured (2026-09-11)

`roman_runs/2026-09-11_h3_footprint/`, `--pair-satellite --events 60 --lenses 20 --nerr 2
--stride-roman 5 --start-index 518`, binary `a3c323d`.

**Why `--start-index 518`.** Roman-covered events -- the only ones H3 can speak about -- exist
only in the footprint stratum, which runs from scan ordinal 518 to 1070 of 1829. Starting at 0
on a machine then carrying ~45 competing processes would have taken ~9 hours just to reach
sightline 518. The cost is that the run consumes a different stretch of the RNG stream from the
superseded one, so the two are not comparable row by row. That comparison was a convenience:
`./fishertest` is byte-identical across the DEVIATIONS 36 fix, which is what establishes that
the `satScale = 1` path is untouched.

### The result

Final sample: the run was stopped after 137 sightlines (scan ordinals 518-654), giving 2,673
paired rows.

| | median `sigma(piE)` ratio | improved | n |
|---|---|---|---|
| control, no Roman epochs at peak | **1.000000** (90.6% bit-exactly 1) | 7.3% | 2,114 |
| Roman covers the peak, **joint** | **0.992353** | 84.5% | 528 |
| Roman covers the peak, **Roman alone** | **0.991085** | 87.3% | 526 |

Sign test on the Roman-covered events: 446 of 527 non-tied improve, one-sided
**`p = 2.1e-62`**. Bootstrap 95% CI on the joint median: **[0.99151, 0.99347]**.

**The estimate did not move as the sample grew.** The median ratio was 0.9927 at n = 71,
0.9924 at n = 111 and 0.99235 at n = 528 -- stable to four decimals across a sevenfold increase,
with the confidence interval narrowing from [0.9913, 0.9959] to [0.99151, 0.99347] as expected.
That is what makes stopping early defensible rather than convenient.

**Coverage caveat.** Those 137 sightlines contain **36 of the scan's 147 footprint sightlines**
(24%), spanning `lon -0.62..-0.12`, `lat -1.44..+0.06` -- a contiguous piece of Roman's
footprint, not all of it. The H3 statistic is a per-event ratio of two forecasts for the same
event, so it is far less sensitive to which sightlines were sampled than an area-weighted yield
would be; but stellar density and extinction do vary across the footprint, so the event mix is
not guaranteed representative. Anyone quoting these numbers as a footprint-wide average should
finish the stratum first (`--start-index 518`, run to ordinal 1070).

**Satellite parallax improves `sigma(piE)` by 0.4-0.9%, and the effect is overwhelmingly
significant.** Highly significant and very small is the physically expected combination: the
two observers are separated by a median of 0.00214 Einstein radii, so the perturbation is tiny,
while the paired design removes essentially all the noise that would otherwise hide it.

Elsewhere it does nothing, correctly:

- `sigma(theta_E)`: median ratio **0.999857**. The astrometric Einstein radius comes from the
  deflection amplitude, not from a parallax baseline.
- `sigma(tE)`: median **0.998**, 7.0% of events worse. (Before the fix: 4.881 and 85.5%.)
- lens mass `relMl`: a sub-percent gain, inherited from `piE`.
- condition numbers, no-satellite / satellite: **1.0070** photometric, **1.0006** astrometric.

### Two things the analysis got wrong on the way, both now recorded

**The validation gate tested a confounded variable** (DEVIATIONS 37). It required the gain to
grow with `du_sat`. But `du_sat = piE * D_perp/AU` and `D_perp` is one observatory's orbit,
identical for every event: `du_sat/piE` spans a factor of **1.14** while `piE` spans **35.6**,
and `corr(log piE, log du_sat) = +0.9985`. The check was testing a dependence on `piE`. It is
replaced by five checks that do follow from the physics -- the control being bit-exact,
`sigma_tE` intact, `sigma(theta_E)` untouched, the two geometries equally conditioned, and a
sign test against the control's exact null. **The `sigma_tE` check is kept unchanged, so the new
gate would still have refused the corrupted run.**

**`nepR_pk` was tried as the replacement monotone axis and also shows no trend** -- and the
reason is physics, not noise. The lowest tercile, median **43** Roman epochs near the peak,
already shows the full gain against the highest tercile with 16,410 epochs. A
simultaneous baseline is a geometric constraint: once a few epochs see the source from both
positions at once, the offset is constrained, and further epochs only reduce photon noise on a
term that is already small. **The gain saturates almost immediately**, so there is no monotone
axis to test.

This supersedes the statement in DEVIATIONS 35.2 that the precision consequence is unmeasured.
It is measured, and it is under one percent.

---

## 5f. The astrometric deflection, measured on the post-H7 v3 table (2026-09-11)

Step H5's script existed but its numbers were labelled **provisional** throughout this file,
because they came off pre-H7 runs and H7 changed which events reach the Fisher step at all. Run
against the v3 table -- which is post-H4 (Roman's own astrometric error model, `6c97375`) and
post-H7 (`f959c8c`) -- they are provisional no longer.

Sample: **82,888 detected events**, of which **8,894 (10.7%)** were observed by Roman at all.
Roman is the only astrometric instrument here worth the name, so every number below is over
those 8,894.

### The signal is far below a single exposure

The centroid of the source's unresolved images is displaced by
`delta(u) = theta_E * u/(u^2+2)`, which peaks at `u = sqrt(2)`, **not** at closest approach. So
the largest deflection an event ever reaches is `theta_E/sqrt(8)` for the **74.1%** whose
trajectory crosses `u = sqrt(2)`, and `theta_E*u0/(u0^2+2)` for the rest.

| max centroid shift | [mas] |
|---|---|
| 5th percentile | 0.0258 |
| median | **0.1106** |
| 95th percentile | 0.4601 |

Roman's per-exposure astrometric precision for **these** sources is **6.69 mas** (median). The
1.1 mas figure quoted for Roman is the bright-source centroiding floor and does not apply here:
`errRomanA` evaluated at each source's own F146 magnitude gives 6.69 mas at the median, which
inverts to F146 ~ 22.98 -- faint bulge main-sequence stars, exactly what a GBTDS microlensing
source is.

So the median event's peak astrometric signal is about **one sixtieth of what a single exposure
can measure**, and only **0.11%** of events have a peak shift exceeding single-exposure
precision at all. The measurement exists entirely through averaging: ~50,000 exposures take the
precision to **0.0300 mas** (median), against which the median 0.1106 mas signal stands at
**3.68 sigma**.

That factor of ~224 in averaging is the whole astrometric result, and that it is available is an
assumption -- see `OPEN_ITEMS.md`. It is the dominant caveat on everything in this section.

### The astrometric peak is not the photometric peak, and seasons cut it

Because the deflection peaks at `u = sqrt(2)` rather than at closest approach, an event's
astrometric maximum falls at `|t - t0| = tE*sqrt(2 - u0^2)` -- twice, symmetrically, up to
~1.41 tE either side of the photometric peak. Measured on this sample, **42.5%** of events have
an astrometric peak that crosses a season edge, and of the events whose photometric peak Roman
does catch inside a season, **52.8%** have an astrometric peak outside one. Catching the
brightening is not the same as catching the wobble.

### theta_E is measured, and it is Roman's alone

Fractional precision `sigma(theta_E)/theta_E`, over events where the astrometric matrix inverted:

| partition | median | better than 10% | better than 1% |
|---|---|---|---|
| joint | **0.2523** | **33.2%** | 4.3% |
| Rubin only | 1.9282 | 0.7% | 0.0% |
| Roman only | 0.2686 | 32.8% | 4.2% |

Rubin's ground-based astrometry is not a capability here -- a median fractional error of 193%.
Its contribution to `theta_E` is a **1.7%** tightening of Roman's number: the paired,
event-by-event median of `sigma_joint/sigma_Roman` is **0.9829**. (Comparing the two medians
instead gives 0.2523/0.2686 = 0.939, a 6% gain -- the ratio of medians is not the median of
ratios, and for a paired quantity the per-event statistic is the right one.) Against Rubin alone
the same paired median is **0.1293**. **Any claim that the joint fit measures the Einstein radius should
attribute it to Roman.**

Note against the pre-H4 record: §4 of this file quotes F4's "90.1% of events measure `tetE` to
better than 10%". It is now **33.2%**. H4 replacing Rubin's error curve with Roman's real one
made the astrometry substantially worse, which is correct -- the placeholder was flattering
Roman -- and it is the size of the distortion H4's verification step asked for.

### The lens mass is the quantity that needs both telescopes

`Ml = theta_E/(kappa piE)` takes one observable from each matrix: `theta_E` from Roman's
astrometry, `piE` from the photometric baseline. It is the only place a change in either shows
up as a change in the science.

| Ml measured better than | joint | Rubin only | Roman only | joint gets it, **neither alone** |
|---|---|---|---|---|
| 100% | 36.0% | 6.6% | 28.1% | **608** |
| 30% | 17.8% | 0.8% | 13.9% | **338** |
| 10% | 7.5% | 0.0% | 5.9% | **145** |

That last column is the synergy result, and it is a real one: 338 events get a lens mass to
better than 30% from the joint fit that neither survey delivers by itself.

**Counted with a precision threshold, deliberately.** A first pass counted events where
`relMl` was merely not the `-1` sentinel and found **zero** masses exclusive to the joint fit --
a meaningless number, because `relMl_L` exists for 99.9% of these events with a *median
fractional error of 15* (1500%). An inverted matrix is not a measured mass. Any future count of
"how many events measure X" in this project needs a threshold, not a sentinel check.

### Two validation checks, both passed

1. **`sigma_joint <= sigma_single`, event by event.** `F[SJOINT] == F[SRUBIN] + F[SROMAN]`
   exactly, so the joint matrix sees strictly more information and its marginalised error cannot
   exceed either single-survey error. **0 violations** over 8,889 joint-vs-Rubin and 8,893
   joint-vs-Roman comparisons.
2. **Precision improves with the size of the deflection.** If it did not, the astrometric matrix
   would not be responding to the deflection at all:

   | max-shift quintile | median shift [mas] | median `sigma(theta_E)/theta_E` |
   |---|---|---|
   | 1 | 0.0357 | 0.9670 |
   | 2 | 0.0702 | 0.4153 |
   | 3 | 0.1106 | 0.2801 |
   | 4 | 0.1729 | 0.1413 |
   | 5 | 0.3403 | 0.0701 |

   Monotone across all five, `corr(log delta_max, log sigma/theta_E) = -0.481`.

Both were computed by a second script (`analysis/h5_crosscheck.py`) written from the physics
rather than from `analysis/h5_astrometric_shift.py`, so the two are an independent cross-check of
each other and not one code path agreeing with itself. Where they overlap they agree to every
digit either prints:

| quantity | `h5_astrometric_shift.py` | `h5_crosscheck.py` |
|---|---|---|
| median max centroid shift | 0.110635 mas | 0.1106 mas |
| fraction reaching `u = sqrt(2)` | 0.741286 | 74.1% |
| per-exposure precision | 6.69197 mas | 6.6938 mas |
| fraction of shifts above one exposure | 0.112435% | 0.1124% |
| sqrt(N)-averaged precision | 0.0300273 mas | 0.0299 mas |
| median signal significance | 3.68449 sigma | 3.69 sigma |
| `sigma_tetE` < 10%, joint / Rubin / Roman | 0.331572 / 0.0073083 / 0.327524 | 33.2% / 0.7% / 32.8% |
| violations of `sigma_joint <= sigma_single` | 0 and 0 | 0 and 0 |

Two independent implementations of `errRomanA` reading the same 2.9 GB table, agreeing on the
per-exposure precision to four significant figures, is what makes the 6.69 mas number safe to
build on -- and it is the number I had wrong before, so it was worth checking twice.

---

## 5h. Step H6 — the whitepaper, written as a blueprint (2026-09-12)

`Whitepaper/whitepaper.tex`, **45 pages**, compiles clean: no undefined citations, no undefined
references. Build: `pdflatex; bibtex whitepaper; pdflatex; pdflatex` from `Whitepaper/`.

**The framing is the user's, and it changed the step.** "The purpose of H6 is for a reader to be
able to use it as a blueprint to do what we are trying to do here, both scientifically and
computationally." So project archaeology is out — no "an earlier version of this code…", no "no
number predating this survives" — and anything instructive inside it was converted into design
guidance a reimplementer can act on. **Deviation 38** has the conversion table, entry by entry.

**What is new in the document**

- **§7 Implementation** (new section): code layout; the three traps (the per-filter vs
  per-telescope indexing, the telescope tag, the `-1.0` sentinel); the output schema;
  determinism and the provenance block; what a full scan costs and where the time goes.
- **§9 Results**: rewritten on the post-H7 v3 table, with seven figures.
- **§10 Open items**: reordered by how much each would move the numbers, led by the
  exposure-independence assumption. Two stale entries retired (H4 landed; H3 is measured).
- **§4.2** rewritten: `errRomanA`'s three regimes with their sources, and the two factors of ten
  (per-exposure vs daily-binned, bright-source floor vs these sources' 6.69 mas).

**Two false claims were removed, not restyled.** §8.4 said the satellite-parallax effect was
unmeasured — true when written on the morning of 2026-09-11, false by that evening. §4.2 said
Roman's astrometry came from an ELT-scale `FWHM` entry and a zero `gama[6]`; `errRomanA()` reads
neither array.

### The F-series was re-run on the v3 table, and the headline numbers moved

Everything in the old §8.1 came from `f1_kroupa.csv`/`f2_kroupa.csv` (2026-08-30), predating H1,
E1a, H4 **and H7**. Re-run against
`/home/ali/Documents/PhD/Offline_project/roman_runs/2026-09-06_v3_h7/test5.dat`:

| | pre-H7 | post-H7 v3 |
|---|---|---|
| characterised jointly, by neither survey alone (in the six fields) | 245 | **351** |
| in-gap median `sigma_piE(joint)/sigma_piE(Roman)` by `tE` bin | 0.31 / 0.80 / 0.95 / 0.984 | **0.433 / 0.936 / 0.978 / 0.991** |
| in-gap median `sigma_tE` ratio, short bin | 0.13 | **0.250** |

**The in-season control is new and is the reason the in-gap row means anything:** the same
statistic for events peaking *inside* a season is **0.98 in every `tE` bin**. Rubin's in-season
contribution is ~2%; the in-gap short-`tE` 0.250 is against that floor, not against 1.0.

Other numbers the document now carries: 82,888 detections of 6,075,044 draws; detection classes
76,146 Rubin+joint / 3,935 Roman+joint / 2,647 both / **160 joint-only**, anomaly 0; inside the
six fields 8,902 detections, characterised joint 2,795, Roman 2,042, Rubin 1,113; fraction
measured better than 10% (joint / Roman / Rubin) `tE` 20.4/14.3/6.3, `piE` 15.9/11.9/5.3,
`tetE` 33.2/32.8/0.7, `Ml` 7.5/5.9/0.0.

### Where the products are

| | |
|---|---|
| detected-event extract (82,888 rows, all columns) | scratch only — regenerate with the `awk` one-liner in `OPEN_ITEMS.md` |
| F1 table | `figures/f1_results_table_v3.csv` |
| F2 / F3 / F4 figures + CSVs | `figures/f2_gap_filling_v3.*`, `f2_gap_filling_piE_v3.*`, `f3_characterization_map_v3.*`, `f4_fisher_precision_v3.*` |
| H3 figures | `figures/h3a_where.png`, `h3b_corner.png`, `h3c_honest.png` (regenerated; `TEMPORAL_GAIN` in `analysis/h3_satellite_parallax.py` is now the v3 in-gap `piE` medians) |
| H5 figures | `figures/h5_astrometric_shift_v3.png`, `h5_astrometry_summary.png` |

All seven figures in the paper carry `git_commit=f959c8c` in their own footers. `whitepaper.tex`
reaches them through `\graphicspath{{./}{../figures/}}` rather than by copying.

**A trap for the next session:** `f1`/`f2`/`f3` still cannot load the 2.7 GB table on this
machine (~2 GB free). The F-series above ran against the `awk` extract. F2/F3/F4 are unaffected
by that — they filter to detections themselves — but **F1's `N_events`, `N_neither` and
`frac_gap_seen_by_rubin` become conditional on detection**, so do not quote a detection
efficiency off that run. `OPEN_ITEMS.md` carries the one-liner and the caveat.

---

## 5i. Step G2 started, and the advisor's paper read against the parent code (2026-09-15)

**Nothing in the code has changed.** Two documentation products and one decision pending.

**1. Deviation 39** — Sajadian & Makler 2026 (arXiv:2608.16448) and its code
(`LSSTHealpix.cpp`, github.com/SSajadian54/MapsMLRubin) are the parent of this pipeline. The
entry tables which inherited choices were bugs and which were right for the LMC black-hole target.
Two facts a later session would otherwise get wrong: **the Deviation 4 overwrite bug is NOT in the
advisor's code** (it came in with this repo's GSL port, `e40716a`), and **the advisor's code
expects a per-visit `sigmaA_LSST.txt`** (10-74 mas), independently confirming Deviation 32.3.

**The scratch measurement** in that entry (fixing `t0` and `mbs0` vs marginalising them, on the
fixture's events, Rubin partition): sigmas shrink by ×1.0-4.2 typically, ×133 for a sparse 5-day
event, and by ×1.9-2.5 on the 900-day event, where the baseline magnitude, not `t0`, is the
degeneracy. The program lived in the session scratchpad and is gone; Deviation 39 gives the exact
method (four subset inversions of `inputA[SRUBIN]`, checked against `Era` to 0.0).

**2. G2's premise, checked against Abrams et al. 2025 (arXiv:2309.15310) — it does not hold as
written** (plan §0.1 rule 6):

- The plan asks to compare the **parallax** characterised fraction (`tE > 2σ`, `piE > 2σ`) at
  l = 0.33°, b = 2.82°. **Abrams et al. publish no absolute value for it** — §3.7 gives only
  OpSim-to-OpSim ratio maps (Figs. 11-14) and a histogram without numbers (Fig. 10).
- The only absolute efficiencies are **Table 6's Fisher metric**: fraction with `σtE/tE < 0.1` in
  bins 10-20 / 20-30 / 30-60 / 200-500 d, e.g. `baseline_v3.0_10yrs` = 0.05 / 0.12 / 0.22 / 0.80.
  But that is **sky-wide** (TRILEGAL N²-weighted HEALPix), **no parallax**, one mean star
  (r = 24.5), fixed 50% blend, `u0` in [0, 1], analytic Fisher over `tE, t0, u0` and a
  **per-band** `F_S, F_B` pair, and OpSims only up to v3.0.
- **Our Rubin visit file does not reach that field**: `Baseline/BulgeBaseline.dat` spans
  b = -3.59 to +1.10; zero visits within 1.75° of (0.33, 2.82). `Baseline/baseline_v5.1.0_10yrs.db`
  (773 MB) is on disk, so a field-specific visit list can be queried. **`rubin_sim` is not
  installed** in `.roman/`.
- Our Rubin light curve uses **one** reference-band flux pair (`RUBIN_REF_BANDS = {2}`, the
  unfinished half of Step C2, Deviation 3); Abrams et al. fit one per band.

**Pending: the user chooses the comparison** (options in the session's G2 brief; recommendation
there is to run Abrams et al.'s own public metric, `rubin_sim`'s `MicrolensingMetric`, and this
pipeline's `FisherM` on the same OpSim, field and parameter-space sample, paired event by event,
and to anchor the `rubin_sim` run by reproducing a Table 6 row first).

---

## 6. Traps a new session will otherwise fall into

- **The two indexing systems.** `s.blend[i]` / `s.magb[i]` are indexed by **filter**
  (0–5 = Rubin `ugrizy`, 6 = Roman `F146`). `s.fb[tt]` / `s.mbs[tt]` are indexed by
  **telescope** (0 = Rubin r-band, 1 = Roman F146). They look alike and are not.
- **`-1.0` is a sentinel, never a measurement.** Three separate bugs have come from summing
  one. On the Python side, always go through `analysis/romanlib.py`.
- **`okA` means the matrix inverted, not that every parameter was fitted.** Which parameters
  were free is decided per event by `activePhotParams()`.
- **`flagi` is stale on uncharacterized events.** Gate on `okA`/`okB`.
- **`test5.dat` takes ~6 minutes just to parse.** When iterating on a figure, cache a subset
  first rather than re-reading 5.57M rows each time.
- **`IMnum` doubles as the mass-function selector and the output-file suffix.** Changing it
  changes which files get written.
- **Data files are gitignored and absent after a fresh clone.** `./roman` must be run from
  the repo root; every data path is hardcoded and relative.
- **Two output files open in append mode**, so a re-run without clearing them silently
  concatenates two runs — `OPEN_ITEMS.md`.
- **`--dry-run` is destructive.** It truncates `test5.dat` and zeroes the `files/MONTLMC/files/`
  outputs before it exits, because the header write comes first. It has already destroyed one
  partial production table. **Never invoke `./roman` from the worktree while a run is writing
  through the worktree symlinks** — the per-event write re-opens the path every time, so a live
  run cannot be protected by moving the symlink. `OPEN_ITEMS.md`.
- **The `analysis/` scripts need the project venv, `.roman/bin/python`.** This trap entry
  previously said `/usr/bin/python3`, which is wrong and would cost a session an hour:
  `/usr/bin/python3` has numpy but **no pandas**, and every `analysis/` script uses pandas.
  Measured 2026-09-11: `/usr/bin/python3 -c "import pandas"` raises `ModuleNotFoundError`,
  while `/home/ali/Documents/PhD/Offline_project/Roman/.roman/bin/python` has pandas,
  matplotlib and scipy. The `/usr/local/bin/python3` that `which python3` finds has neither.
  Corrected in place rather than left standing, because a trap list that is itself wrong is
  worse than no trap list.
- **`lightcurve()` takes a telescope argument** since Step H1: Roman observes from L2 and
  Rubin from the ground, so the two see different observer displacements. Passing the wrong one
  silently removes satellite parallax. Any table written before H1 has none — do not describe
  its `piE` forecasts as including it. Deviations 27, 28.
- **The `t = 0` parallax gauge.** `lightcurve()` subtracts the observer displacement at
  `t = 0`, which is what makes `u0` and `t0` mean what they mean. When two observers exist they
  must share one origin — referencing each to its own `t = 0` cancels exactly the offset that
  *is* the satellite parallax, silently. `PHASE_H_PLAN.md` Step H1.
