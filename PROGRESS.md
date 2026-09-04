# PROGRESS.md — where this project stands

**Last updated:** 2026-09-04, after Steps E1a (footprint-stratified sampling), H1
(satellite parallax) and H4 (Roman's astrometric error model).
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

0. **The production run — this is the bottleneck now.** E1a (stratified sampling), H1
   (satellite parallax) and H4 (Roman astrometric errors) are all in the code and all change
   what a run produces. Nothing downstream can move until one run exists with all three.
   Two decisions, both the user's: `--stride-roman` (see item 1) and whether to also do the
   `--no-satellite-parallax` twin run that Step H3 needs. **Every result in §4 predates all
   three changes.**

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
3. **The rest of Phase H.** H1 and H4 are done. Remaining: **H2** (record `du_sat` and
   contemporaneous-coverage per event), **H3** (the satellite-parallax experiment and its three
   figures — needs the twin runs), **H5** (the astrometric-shift analysis product, now
   unblocked by H4), **H6** (the whitepaper for collaborators). Full text and the dependency
   graph in `PHASE_H_PLAN.md`. **H3 replaces the original plan's Step G1**, which cannot be run
   as written (Deviation 27).

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
- **`lightcurve()` has no telescope argument**, so both observatories sit at the centre of the
  Earth and there is no satellite parallax. Do not describe any current `piE` forecast as
  including it. Deviation 27; Step H1 fixes it.
- **The `t = 0` parallax gauge.** `lightcurve()` subtracts the observer displacement at
  `t = 0`, which is what makes `u0` and `t0` mean what they mean. When two observers exist they
  must share one origin — referencing each to its own `t = 0` cancels exactly the offset that
  *is* the satellite parallax, silently. `PHASE_H_PLAN.md` Step H1.
