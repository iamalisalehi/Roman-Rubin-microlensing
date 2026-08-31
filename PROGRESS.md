# PROGRESS.md — where this project stands

**Last updated:** 2026-08-31, after the `sigma_piE` version of F2.
**Branch:** `joint-fisher-refactor` (never commit to `main`).
**Head at last update:** `4f3bb0d` — "Test the plan's parallax claim on the parallax itself".

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
Roman alone already measures the parallax of long events (Deviation 24). What remains is
Phase E's sampling work (which would sharpen the Roman-footprint statistics and the
sample-limited long-`tE` bins), Phase G (physics separation, validation against published
Rubin-only numbers, and reconciling the whitepaper), and the deferred GBTDS-footprint item.

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

---

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

In roughly the order that makes sense:

1. **Step E1 — stratify sampling in `tE`, and weight the Roman footprint.** This is what
   turns panel (a) of F3, the whole short-`tE` corner, **and the long-`tE` null of the new
   `sigma_piE` figure** (124 events) from suggestive into quotable. A footprint-weighted run
   buys roughly 20× more Roman-observed events for the same wall-clock, at the cost of
   explicit weights to recover survey-wide totals. This is now the single highest-value next
   step: three separate results are sample-size-limited and all three are limited by the
   same thing.
2. **Step E2 — revisit the per-sightline stopping criteria.** The current third floor (2
   well-conditioned events) was written when there was one Fisher matrix; with three
   matrices and stratified bins it needs replacing.
3. **Phase G.** G1: separate satellite parallax from temporal-baseline parallax by rerunning
   with Rubin's observer position forced to Roman's. G2: validate the Rubin-alone branch
   against Abrams et al. 2025 at l = 0.33°, b = 2.82°, **reweighting to their sampling
   first** or the comparison will look like a bug. G3: fold `OPEN_ITEMS.md` into the
   whitepaper's open-items discussion.
4. **The GBTDS footprint item** — deferred by the user, but wanted. The sky coverage is
   Penny et al.'s and the footprint has since changed; Rubin's large field of view can see
   corners of the GBTDS region whose centre it is not pointed at, which matters for blending
   too, and `BulgeBaseline.dat` must change with it. **`OPEN_ITEMS.md` has the full text.**

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
