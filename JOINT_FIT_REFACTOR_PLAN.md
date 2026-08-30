# Roman + Rubin Joint Microlensing Forecast — Guided Refactor Plan

**Repository:** `github.com/iamalisalehi/Roman-Rubin-microlensing`
**Primary files:** `Bulge.h`, `Bulge_LSST.cpp`, `Lensing.cpp`, `helper.cpp`
**Python pipeline:** `interpolate2.py`, `maps2.py`, `BolometricCorrection.py`, `generateRomanBaseline.py`
**Document:** `whitepaper.tex`, `refs.bib`

---

## 0. How to work through this document

Read this whole section before touching any code.

### 0.1 The working agreement

I (the user) am a PhD student in astronomy. I did not write most of this C++ from scratch — it was adapted from my advisor's legacy simulation for a different science case (LSST + ELT toward the LMC), and I have refactored parts of it. **There are functions and variables in here whose purpose I do not fully understand.** I want to end this process understanding the code, not just having a working version of it.

So this is a *teaching* refactor. The rules:

1. **One step at a time.** Do exactly one numbered step, then stop and wait for me to say "continue." Do not batch steps. Do not work ahead. Do not silently fix something you notice in a later step — tell me about it and let me decide when to address it.
2. **Explain before you edit.** For every step, first produce the teaching brief (Section 0.2), then show me the proposed diff, then wait for approval, then apply it.
3. **Explain the physics and the astronomy, not just the code.** When you touch a variable, tell me what physical quantity it represents, what units it is in, and what would go wrong scientifically if it were wrong. "This is a `double` that gets passed to `FisherM`" is not an acceptable explanation.
4. **Never assume I know an acronym or a variable name.** If you mention `piE`, `tetE`, `murel`, `fb`, `nsbl`, `dchiL`, `ndw`, `coun`, `minc`, `cade` — define it, in that message, in words.
5. **Flag uncertainty honestly.** If you cannot tell what a legacy variable does, say so explicitly and propose a way to find out (grep for all uses, add a temporary print, check against the whitepaper). Do not invent a plausible-sounding explanation. A wrong confident explanation is worse than "I don't know yet."
6. **If a step reveals that my premise was wrong, stop and say so.** This plan was written from a partial reading of the code. If you open the file and find that something I asserted is not actually true, tell me before proceeding.

### 0.2 Required format for every step

For each numbered step, produce, in this order:

**A. Teaching brief** (before any edit)
- *What this part of the code currently does* — narrative prose, not a line-by-line paraphrase.
- *Every variable and function you are about to touch*, in a table: name, type, physical meaning, units, where it is set, where it is read.
- *Why the current behaviour is wrong or insufficient* for the science goal (Section 0.3).
- *What the change will do*, and what observable difference it makes to the output.
- *What could break* — including numerically, not just at compile time.

**B. Proposed diff** — show it, do not apply it yet. Wait for my approval.

**C. Verification** (after I approve and you apply)
- `g++ -std=c++17 -fsyntax-only <file>` must pass. This is my standard check; do it every time.
- State exactly what I should look for in the output to confirm the change worked.
- If the step has a numerical acceptance criterion (several do), state whether it passed and show the numbers.

**D. Commit** — one `git commit` per step, with a message that names the physical change, not the mechanical one. Good: `Separate Roman and Rubin detection bookkeeping so per-survey yields are recoverable`. Bad: `add ndw_R variable`.

### 0.3 The science goal this refactor serves

Read this carefully; every step below exists to serve it.

The project forecasts the combined microlensing yield and parameter-recovery precision from the **Nancy Grace Roman Space Telescope's Galactic Bulge Time Domain Survey (GBTDS)** and **Rubin Observatory's LSST**, toward the Galactic bulge.

Roman observes the bulge in short, intense bursts: six high-cadence seasons of ~72 days, imaging each field roughly every 12.1 minutes in the wide F146 filter, with long gaps in between (including a multi-year gap in the middle of the mission). Rubin observes the bulge far less intensely (~3-day cadence, six optical bands `ugrizy`) but *continuously* across the whole decade.

**The scientific question:** which microlensing events can be detected and characterized by *both* facilities, and how much does combining them help? The expected answer is asymmetric:
- **Short-duration events** (small Einstein crossing time `tE`, i.e. low-mass lenses and free-floating planets): Roman dominates during its seasons and Rubin adds nearly nothing. But events peaking *in a Roman gap* are invisible to Roman entirely — Rubin either catches them or nobody does.
- **Long-duration events** (large `tE`, i.e. black holes and neutron stars): these span multiple Roman seasons. Rubin's coverage of the gaps is what lets the *annual parallax* signal be sampled, which is what turns a `tE` measurement into a lens *mass* measurement.

**I want quantitative numbers for both regimes**, including the small-gain short-`tE` regime. A near-unity result there is a publishable result, not a null.

The two headline novelty claims of the thesis are:
1. **Joint stitched Fisher-matrix fits** across both telescopes (this does not currently exist in the literature).
2. **Quantifying how Rubin fills the gaps in Roman's observing seasons.**

**The single most important structural fact:** the code currently produces *one* detection boolean and *one* Fisher matrix per simulated event, computed over the merged Roman+Rubin data stream. That architecture **cannot answer the question above**, because it cannot separate "Roman found it" from "Rubin found it" from "the combination found it." Most of this plan is about fixing that.

### 0.4 What is already right — do not "fix" these

- **`tele[]` tagging.** Each light-curve datum is tagged with which observatory produced it, and `FisherM` already does `tt = int(l.tele[i])` to index `s.fb[tt]` and `s.mbs[tt]`. This is the key architectural asset that makes the whole refactor cheap. Preserve it.
- **Detection-then-Fisher ordering.** The Monte Carlo detection test runs first; `FisherM` is only invoked on events that pass. This is correct — forecasting precision for an event no survey would ever flag is meaningless. Keep this ordering.
- **Separate photometric and astrometric Fisher matrices.** Two matrices, not one. This is deliberate: breaking the lens-mass degeneracy needs both the microlensing parallax `piE` (from photometry) and the angular Einstein radius `tetE` (from Roman's astrometry). Keep both.
- **Real cadence on both sides.** Real Rubin OpSim visits and a real ROTAC-2025-derived Roman visit list, not idealized cadences.
- **In-memory `std::vector<EventRecord>` buffering** instead of writing to disk and re-reading (the legacy code lost precision on a text round-trip).

### 0.5 Setup before Step 1

Create a branch: `git checkout -b joint-fisher-refactor`. Every step commits to this branch. Do not touch `main`.

Establish a build baseline: run the syntax check on all `.cpp` files and record the result, so we know whether any later failure is ours.


### 0.6 Keeping the record — added 2026-08-31, and binding from here on

*This subsection was added after Phase F. Everything above it is the plan as originally
written and is left untouched by design.*

Sessions end — they run out of context, the laptop hibernates, the power goes out. All three
have happened during this work. **Knowledge that lives only in a session transcript is lost
when that session dies.** So writing it down is part of doing a step, not cleanup afterwards,
and it does not wait to be asked for.

Three documents, maintained by every session:

1. **`DEVIATIONS.md`** — one numbered entry per departure from this plan, per pre-existing
   bug found, and per result that contradicts what this plan expected. Give what the plan
   said, what was done instead, why, the commit, and the verification numbers.
2. **`OPEN_ITEMS.md`** — one entry per problem noticed and deliberately not fixed. Give what
   is wrong, why it matters scientifically, why it is deferred, and what the fix would
   involve. This is the mechanism that makes rule 1 of §0.1 — never silently fix something
   noticed in a later step — actually work.
3. **`PROGRESS.md`** — the entry point for a new session: head commit, what is done, where
   the current data live, the headline numbers, what is next, and the traps. Update it at the
   end of any session that changed the state of the work.

**This plan's own text is never edited** (this subsection and any future rules being the
exception). It records what was believed when it was written, which is worth preserving —
including where it turned out to be wrong. Corrections go in `DEVIATIONS.md`.

Write as the work happens, not in a batch at the end: an interrupted session must still
leave a usable record. Before starting a long-running job, record what it is, where its
output will land, and how to tell whether it finished.

---

## PHASE A — Orientation

### Step A1. Map the code and confirm my assumptions

**Do not edit anything in this step.**

Read `Bulge.h`, `Bulge_LSST.cpp`, `Lensing.cpp`, `helper.cpp` and produce a written orientation document (`ORIENTATION.md`, committed) covering:

1. **The struct inventory.** For `source`, `lens`, `astromet`, `covarian`, and any others: every member, its physical meaning, its units, and whether it is per-filter, per-epoch, or scalar. Flag any member you cannot identify.

2. **The main loop structure.** Narrate, in prose, the nesting: sky position → field visit-matching → star generation → light-curve time loop → detection test → Fisher call → per-field aggregation. Give approximate line numbers.

3. **The filter convention.** `M` is the number of filters. Confirm the index order (I believe `0-5` are LSST `ugrizy` and `6` is Roman `F146`, based on `s.fb[0] = s.blend[2]  // r-LSST` and `s.fb[1] = s.blend[6]  // F146`). Confirm whether this convention is consistent between the C++ and the Python pipeline (`BolometricCorrection.py` has a `FILTER_ORDER`) — a column-order mismatch between them has bitten this project before.

4. **The two indexing systems that look alike but aren't.** `s.blend[i]` and `s.magb[i]` are indexed by *filter* (0..M-1). `s.fb[tt]` and `s.mbs[tt]` are indexed by *telescope* (0 = Rubin r-band, 1 = Roman F146). Confirm this and state it plainly — I want to be sure I'm not confusing them.

5. **The epoch-matching machinery.** Explain `matchVisibleEpochs()`, `ls->ct[]`, the cursors (`gi`/`giR`, `sq`/`sqR`), `FoV` vs `FoVRoman`, `minc`, `cade`, and the adaptive timestep `dt`. Specifically: how does a simulated timestep get identified as landing on a real visit, and how does the code know which observatory that visit belongs to?

6. **A list of every variable you could not confidently identify.** This list is an output of the step, not a failure.

**Acceptance:** I read `ORIENTATION.md` and confirm it matches my understanding, or correct you. Do not proceed until I do.

### Step A2. Known-blocker sweep

Address the three sizing/config issues already identified, as one small step, because later steps will crash without them:

- **`coun`** is currently `500`. This is the allocated length of the per-event light-curve arrays (`timn`, `magn`, `errm`, `soux`, `souy`, `erra`, `tele`). The Roman baseline is **309,084 visits** across 6 fields and 10 seasons; a star in a high-cadence season will produce far more than 500 usable epochs. Explain what `coun` bounds, what happens on overflow (does `CHECK(ndw <= ndd)` catch it, or is there a silent buffer overrun?), and raise it appropriately. Tell me the memory cost of the new value.
- **`N3` and `N4`** in `Bulge.h` need updating after the Besançon age-filter cutoffs were raised (thick disk to ≤13 Gyr, halo to ≤14 Gyr). Explain what N3 and N4 count and where they are used before changing them.
- **Whitepaper drift (note only, no edit yet):** `whitepaper.tex` §5.1 still describes a *single merged visit list* with `FoV = 1.75°`. That contradicts the separate-cadence refactor already in the code. Log this in a running `OPEN_ITEMS.md` file; we fix the document in Phase G.

---

## PHASE B — Per-survey bookkeeping

This phase is what makes the science question answerable at all. Nothing downstream works without it.

### Step B1. Split detection bookkeeping by observatory

**Current behaviour:** the detection machinery accumulates over the merged stream. There is one `ndw` (number of usable data points), one set of running χ² statistics (`chi1`, `chi2`, `chi3` for the lensed model, the no-parallax model, and the unlensed baseline; plus `chi1a`, `chi2a`, `chi3a` for the astrometric equivalents), one three-consecutive-3σ-point run test (`flag0`, `flag1`, `flag2` → `flag_det`), and one final detection decision using `dchiL` (the lensing Δχ²), `dchiP` (the parallax Δχ²) and `dchiA` (the astrometric deflection Δχ²).

**Why that's wrong:** it yields a single boolean per event. We cannot recover "Roman detected it," "Rubin detected it," or "only the combination detected it." Both novelty claims depend on that distinction.

**The change:** duplicate the accumulators per observatory. Keep the merged/joint versions too — we want all three.

- `ndw` → `ndw` (joint), `ndw_L` (Rubin), `ndw_R` (Roman)
- χ² accumulators → `_L` and `_R` variants alongside the existing joint ones
- the run test (`flag0/1/2`, `flag_det`) → independent per observatory, since three consecutive Rubin points and three consecutive Roman points are different statements
- final decision → produce a **detection label** per event in `{neither, Rubin-only, Roman-only, both}`

The `tele[]` tag already tells you which observatory each datum belongs to, so this is bookkeeping, not redesign.

**Teach me:** what Δχ² actually measures here and why the threshold is `> 2.0 * ndw` (a significance floor that scales with how much data the event has); why three consecutive 3σ points is the criterion (it's in the same family as OGLE's original Early Warning System, which used five); and why the run test must be per-observatory rather than allowed to mix Roman and Rubin points.

**Acceptance:** for a small test run, print the counts in each of the four detection categories for one field. Sanity check with me: Roman-only should dominate in-season, Rubin-only should appear for gap-peaking events.

### Step B2. Make pre-selection per-survey

**Current behaviour:** before any light curve is generated, a cheap pre-selection requires the source to be detectable (between saturation and threshold) in **at least two filters** (`fdet > 1.0`), and draws a random number against `s->blend[2]` — the **LSST r-band** blend fraction — to decide whether this star is the one actually being monitored.

**Why that's wrong:** a heavily reddened bulge source that is bright in Roman's near-infrared F146 but invisible in Rubin's optical `ugrizy` gets killed *before a light curve is ever generated*. That is exactly the population where Roman is uniquely strong, and exactly the population that makes the "Roman-only" and "both" categories interesting. The current cut systematically biases those categories low in the reddest fields — and the bulge is the reddest sightline there is.

**The change:** pre-selection becomes a per-survey OR — `detectable_by_Rubin` (≥2 of the six optical bands within saturation/threshold limits) **OR** `detectable_by_Roman` (F146 within its limits) — and the blend draw uses each survey's own band rather than r for everything.

**Teach me:** what "blend fraction" (`blend[i]`, and the derived `fb`) physically means — the fraction of the total flux in the photometric aperture that comes from the source star itself, as opposed to the lens and unresolved neighbours — and why it is degenerate with `tE` and `u0` in a light-curve fit. Also explain the saturation/threshold arrays (`satu[]`, `thre[]`) and where they come from.

### Step B3. Per-instrument blending

**Current concern:** Roman's PSF is roughly 0.1 arcsec; Rubin's seeing toward the bulge is roughly 0.7 arcsec. The number of stars blended into one photometric aperture therefore differs enormously between them.

**Investigate first, then fix.** `s.blend[i]` is already per-filter, which is good. But check whether `s.nsbl[i]` — the number of blended stars — is computed with a single crowding radius applied across all `M` filters, or with a per-instrument PSF. If it's a single radius, F146 blending is badly overestimated.

**Why this matters for the science, not just for realism:** part of Roman's advantage in the joint fit is that it *deblends* what Rubin cannot. If we model Roman's crowding as if it had Rubin's PSF, we **understate the joint gain** — which biases the headline result in the conservative direction, but wrongly.

Report what you find before changing anything. If it already handles this correctly, say so and move on.

---

## PHASE C — Fisher matrix correctness

Everything in this phase is about making the *precision* numbers trustworthy before we start comparing them.

### Step C1. Add `t0` to the photometric parameter set

**Current behaviour:** the photometric Fisher matrix is built over five parameters: `u0` (impact parameter — the minimum lens–source separation in units of the Einstein radius), `tE` (Einstein crossing time, days), `fb` (blend fraction), `piE` (microlensing parallax), and `xi` (trajectory angle, radians). The matrix dimension is `Nx`. Derivatives are finite-difference, with step sizes in `co.Delta1[]`.

**Why that's wrong:** `t0` — the time of closest approach, i.e. the epoch of peak magnification — is a free parameter of any real microlensing fit and is **missing**. Omitting a genuinely free parameter makes the forecast uncertainties on the remaining parameters too optimistic, because you are implicitly asserting perfect knowledge of when the peak occurred.

**Why this specifically threatens *this* project:** the bias from omitting `t0` is worse for sparsely sampled data than for densely sampled data. Rubin is the sparsely sampled one. So the omission systematically **flatters the Rubin-alone column** — precisely the comparison the entire paper rests on. This is not a cosmetic correction.

**The change:** extend the photometric parameter vector to `{u0, t0, tE, fb, piE, xi}`, increment `Nx`, add the corresponding `Delta1[]` step, and add the perturbation/restore branches in both the `j` and `k` loops of `FisherM`. Note the reference implementation in the literature (Abrams et al. 2025, ApJS 276, 10) fits `tE, t0, u0` plus blend and baseline flux per passband — so this brings us in line with the standard.

**Teach me:** what a Fisher matrix *is* in this context — the local curvature of the χ² surface around the true parameters, whose inverse gives the Cramér–Rao lower bound on achievable parameter uncertainties — and why we use it instead of actually fitting hundreds of thousands of simulated light curves. Then explain concretely why adding a parameter can only ever *increase* the forecast uncertainties on the others.

**Acceptance:** on a handful of test events, σ(`tE`) must increase (or stay equal) after adding `t0`. If it decreases anywhere, something is wrong — stop and investigate.

### Step C2. Per-band source and blend flux

**Current behaviour:** `l->magn[ndw] = magni[2];` with the comment `// scale to r-LSST band, but the errors are different`. Every light-curve point is stored as an r-band-equivalent magnitude. The Fisher matrix carries a **single** `fb`, selected per-telescope via `s.fb[tt]`.

**Why that's wrong for a joint fit:** aligning all data to one reference band is standard practice for a *single-observatory* relative metric (Abrams do exactly this, using r as the reference baseline). But for a **joint** fit it destroys information. F146 and r have genuinely different source fluxes and different blend fluxes for the same star, and **the flux ratio between bands is one of the things that breaks the `fb`–`tE`–`u0` degeneracy.** Collapsing to one band throws away part of the very signal the paper is trying to measure, and can bias the recovered blend fraction.

**The change:** give each band its own source flux and blend flux as free Fisher parameters, rather than a single shared `fb`. Store the model magnitude in its native band. Be explicit with me about how much this grows `Nx` and what that does to runtime and conditioning — this is the most invasive change in the plan and it interacts with Step C4.

**Teach me:** why microlensing magnification is *achromatic* (the same `A(t)` multiplies the source flux in every band) but the *observed* light curve is not (because `fb` differs per band), and why that chromatic difference is information rather than noise.

### Step C3. Finite-difference step-size convergence test

**Current behaviour:** the derivative steps are hard-coded constants — `co.Delta1[0] = 0.1507586576` for `u0`, `co.Delta1[1] = l.tE * 0.254674` for `tE`, `co.Delta1[3] = 0.2509463534656 * l.piE` for `piE`, and so on. These came from the legacy LMC codebase and were tuned for a different problem.

**Why this matters more than it looks:** I intend to claim precision improvements that may be at the few-percent level for short-`tE` events. A finite-difference derivative with a badly chosen step produces errors at exactly that scale. If I cannot show that σ is stable against the step size, my "gain" could be a differencing artifact.

**The change:** add a diagnostic mode (compile flag or runtime switch — your call, tell me which and why) that, for a handful of representative events spanning the `tE` range, sweeps each `Delta1[]` entry over about a decade and records the recovered σ. Produce a plot. σ should show a plateau; we pick steps in the middle of it.

**Acceptance:** σ must be stable to well under 1% across the plateau for every parameter. If a parameter has no plateau, stop — that parameter's derivative is not trustworthy and we need to discuss it. This plot goes in the whitepaper appendix; a referee will ask for it.

### Step C4. Conditioning and safe inversion

**Current behaviour:** the Fisher matrix is assembled in raw physical units and inverted with GSL.

**Why that's a problem:** the parameters span wildly different magnitudes — `tE` ~ 30 (days), `u0` ~ 0.3 (dimensionless), `piE` ~ 0.1, `fb` ~ 0.5, `xi` in radians. A matrix built from such heterogeneous scales is badly conditioned, and inverting a badly conditioned matrix introduces errors at the few-percent level. Again: exactly the level of the effect I am trying to measure.

**The change:**
- Normalize before inverting: divide row and column *i* by `sqrt(F_ii)`, invert the resulting correlation-like matrix, then rescale the result. Explain to me why this is mathematically equivalent but numerically far better behaved.
- Compute and **store the condition number per event**.
- Replace any crash-or-garbage behaviour on singular matrices with an explicit "not characterizable" outcome.

**Teach me:** what a condition number is and what it means for a matrix to be ill-conditioned, in plain terms — and why an ill-conditioned Fisher matrix is often telling you something physically true (the data genuinely cannot constrain some parameter combination) rather than merely being a numerical nuisance.

### Step C5. Three Fisher matrices per event

**This is the step the entire thesis novelty claim rests on.**

**Current behaviour:** `FisherM` loops `for (int i = 0; i < ndw; ++i)` over **all** data points, producing one matrix.

> **Correction (Step C0, 2026-08-18).** This step's original text said `FisherM` *accumulates* over
> all data points. It did not. Both accumulation sites used `gsl_matrix_set` (overwrite) rather than
> a running sum, and the zeroing loops run once *before* the data loop, so every epoch clobbered the
> previous one and the final matrix held a single data point's contribution -- rank 1, hence singular,
> hence `invert_matrix` nudged the diagonal by `1e-10` and returned an inverse of order `1e10`. Every
> sigma reported before Step C0 was meaningless (a real detected event showed relative sigma of 1.8e7
> on `tE`). Fixed in commit `d5c8867`; the partitioning described below is now built on a matrix that
> genuinely sums information over epochs. This also explains the degenerate synthetic-event result
> during Step C1, which was misdiagnosed at the time as a bad test fixture.

**The change:** in the same pass over the data, accumulate **three** matrices:
- `F_joint` — all points
- `F_roman` — points where `tele[i] == 1`
- `F_rubin` — points where `tele[i] == 0`

Invert each independently. The derivatives are identical; only the accumulation is partitioned. This should be roughly thirty lines.

**Critical detail — the paired-comparison property.** Because all three matrices come from the *same* simulated event with the *same* random draws, the comparison between them is **paired**. Everything that dominates the scatter between events — source magnitude, `u0`, blending, sightline extinction — cancels in the per-event ratio. This is what makes a 2% median gain measurable with a few thousand events instead of millions.

**Therefore, an absolute rule for all downstream analysis:** never compute `mean(sigma_joint) / mean(sigma_roman)`. Always compute the ratio **per event** and then take the median and quartiles of the distribution of ratios. Write this rule as a comment in the code where the ratios are formed, so future-me does not undo it.

**Handling the singular cases.** For short-`tE` events, `F_rubin` will frequently be singular or near-singular — Rubin alone genuinely cannot constrain the event. Do **not** let those events silently drop out of the averages: that would select only the Rubin-favourable cases and bias the short-`tE` gain upward. Record them as an explicit category and report both the fraction that are singular and the median ratio among the well-conditioned subset.

**Acceptance:** for events peaking deep inside a Roman high-cadence season, `sigma_joint / sigma_roman` should be very close to 1. For long-`tE` events peaking in a Roman gap, it should be meaningfully below 1. If you see the opposite, stop — adding data can never *increase* Fisher information, so a ratio above 1 is a bug, not a result.

---

## PHASE D — The output table

### Step D1. Extend `EventRecord` to the full per-event row

Every figure in the paper is a cut on one flat table. Build it once, properly.

Each simulated event contributes one row:

**True (input) quantities:** `tE`, `t0`, `u0`, `piE`, `xi`, `tetE` (angular Einstein radius), lens mass `Ml`, lens distance `Dl`, source distance `Ds`, relative proper motion `murel`, per-band baseline magnitudes and blend fractions, sightline `(lon, lat)`, and the optical depth for that sightline.

**Detection:** `det_roman`, `det_rubin` (booleans), the derived four-way label, `ndw_R`, `ndw_L`.

**Precision — three columns for each quantity:** σ(`tE`), σ(`piE`), σ(`tetE`), σ(`Ml`), each for Roman-alone, Rubin-alone, and joint. Nine to twelve columns.

**Diagnostics:** condition number for each of the three matrices; a singular/well-conditioned flag for each.

**Gap geometry — this one is essential and easy to forget:** `dt_to_season_edge`, the signed time from `t0` to the nearest Roman observing-season boundary, and a flag for whether `t0` falls inside a Roman season, inside a gap, or outside the mission entirely. **The gap-filling result is a plot against this variable**, so it must be computed and stored at simulation time.

Keep the in-memory `std::vector<EventRecord>` approach. Write out as a single well-headed file per field.

### Step D2. Check `t0` sampling covers the gaps

**Investigate:** confirm that `t0` is drawn uniformly across the **full ten-year window** including the Roman gaps, not just within observing windows.

**Why:** if `t0` sampling is tied to the observation window, we undersample exactly the gap-peaking events the entire paper is about, and the gap-filling result will be biased or simply absent. Report what you find. This is a small check with a large consequence.

---

## PHASE E — Sampling strategy

### Step E1. Stratify in `tE`, then reweight

**The problem:** a population-weighted Monte Carlo produces very few long-`tE` events, because they are intrinsically rare. But the `tE > 200` day tail is where the headline black-hole result lives.

**The change:** run a fixed number of events **per `tE` bin**, then reweight by the Besançon-derived `tE` distribution to recover absolute yields per square degree. This decouples "how precisely do I know the gain in this bin" from "how many such events exist," and lets us push N up only where the effect is interesting.

**Use these bin edges**, which match the reference literature so our Rubin-alone column stays directly checkable against published values: 1–5, 5–10, 10–20, 20–30, 30–60, 60–90, 100–200, 200–500, 500–1000 days.

**Teach me:** importance sampling and reweighting, concretely, with the actual weight expression for this case. Be explicit about which reported quantities are per-bin (unweighted, conditional on `tE` being in that bin) and which are absolute yields (reweighted). Mixing these up is an easy way to publish a wrong number.

### Step E2. Revisit the stopping criteria

**Current behaviour:** each field runs until three floors are met — at least 850 detected stars, at least 150 detected lensing events, and at least 2 events with a well-conditioned Fisher matrix.

**The issue:** the third floor is very weak, and it was written when there was only one Fisher matrix. With three matrices and stratified `tE` bins, the criteria need rethinking — we now need enough well-conditioned events *per bin per matrix*. Propose a replacement and explain the run-time trade-off.

---

## PHASE F — Analysis and results

### Step F1. The results table generator

Write a Python analysis script that reads the `EventRecord` output and produces, per (field, `tE` bin):

| Quantity | Regime where it is the headline |
|---|---|
| N detected: Roman-only / Rubin-only / both / neither | both regimes |
| Fraction of gap-peaking events recovered by Rubin | **short `tE`** |
| Median σ_joint/σ_tE,Roman (well-conditioned only), with quartiles | **long `tE`** |
| Median σ_joint/σ_piE,Roman, with quartiles | **long `tE`** (watch the short-`tE` FFP corner) |
| ΔN meeting the characterization criterion, joint minus Roman-alone | **long `tE`** |
| Fraction where Rubin-alone Fisher is singular | both regimes |

**Characterization criterion:** use `tE > 2*sigma_tE` **and** `piE > 2*sigma_piE`. This is deliberately the same criterion used by Abrams et al. 2025 for their two-parameter parallax characterization (they note it is appropriately looser than the single-parameter `sigma_tE/tE < 0.1` because two parameters are being constrained at once). Using their criterion makes our Rubin-alone numbers directly comparable to published values.

**Remember the three currencies.** Short-`tE` gain is mostly a *yield* statistic (Roman literally cannot see gap-peaking events, so σ_Roman is infinite and the ratio is undefined). Long-`tE` gain is a *precision* statistic. Do not force one metric across both regimes — that produces a misleadingly flat answer.

### Step F2. The gap-filling figure

Plot median `sigma_joint / sigma_Roman` against `dt_to_season_edge`, split by `tE` bin.

Expected shape: flat at ~1 for events peaking mid-season; dropping as `t0` moves into a gap; the drop deepening with `tE`. **This curve *is* the gap-filling result** — the thing that has been advocated in white papers for years and never simulated.

### Step F3. The (`tE`, `piE`) characterization map

A 2D histogram in `log(tE)` vs `log(piE)`, colored by the ratio of joint-characterized to Roman-alone-characterized fractions.

**Deliberately mirror the figure format of Abrams et al. 2025 (their Figures 11–14)**, which show exactly this plane with an OpSim-ratio colour scale. The community already knows how to read that plot; matching the format is worth real effort. Note that in that plane, more massive lenses (black holes) sit toward long `tE` and small `piE`, while low-mass lenses sit at short `tE` and large `piE`.

---

## PHASE G — Physics separation, validation, documentation

### Step G1. Separate the two parallax mechanisms

This is a cheap experiment with a genuinely novel result.

**The physics I need you to get right, because the literature framing is easy to over-claim:** Roman sits at L2, roughly 0.01 AU from Earth. The induced difference in the impact parameter seen by the two observatories is of order `0.01 * piE` — around 10⁻³ for typical bulge events. So **simultaneous "satellite parallax" between Roman and Rubin is a narrow niche** (high magnification, small projected Einstein radius, free-floating-planet-like lenses), not a headline. Abrams et al. say as much: they note that for short events a single telescope usually cannot measure `piE`, and that Roman–Rubin satellite parallax is likely viable only for short events because the baseline between the telescopes is short.

**The dominant joint gain is temporal, not spatial:** Rubin extends the time baseline across Roman's gaps so that the *annual* (Earth-orbit) parallax signal is actually sampled.

**The experiment:** run the joint fit twice — once with real geometry, and once with Rubin's observer position forced to Roman's (which kills the spatial baseline while preserving the timing). The difference isolates the satellite-parallax contribution from the temporal-baseline contribution.

Do the same for the astrometric matrix.

**Deliverable:** two separate figures with two separate captions. Do not merge them into one "parallax improvement" number — that would conflate two different physical effects and invite exactly the criticism the honest framing avoids.

### Step G2. Validate against published Rubin-only results

Run the **Rubin-alone** branch of the pipeline on the `baseline_v3.0_10yrs` OpSim in the bulge field used by Abrams et al. — RA = 263.89°, Dec = −27.16° (l = 0.33°, b = 2.82°) — and compare the characterized fraction against their published value.

**Essential caveat, or the comparison will look like a bug in our code:** their simulated sample is a *parameter-space* survey, not a population. They draw `u0` uniform in [−1, 1] and `log(tE)` uniform from about 5 to 600 days, and for the MAF metrics they use a single mean TRILEGAL star with a fixed ~50% blend fraction. Ours is population-weighted, per-star Besançon + MIST photometry with real 3D dust and real per-band blending. **We must reweight to their sampling before comparing, or state the difference explicitly.**

This validation figure is worth a lot: it disarms the most obvious referee objection, which is "why should we trust a new simulator?"

### Step G3. Sync the whitepaper

Update `whitepaper.tex` to match the code as it now stands:

- **§5.1** currently describes a single merged visit list with `FoV = 1.75°`. Rewrite for the separate Roman/Rubin cadence handling, `FoVRoman`, and the duplicate-timestamp deduplication in the overlap zone between adjacent Roman fields (Roman's field of view exceeds half the field spacing, so overlap produces tied timestamps that must be skipped, not treated as errors).
- **§5.3–5.4** rewrite for the per-survey detection labels and the revised stopping criteria.
- **§6** rewrite for the six-parameter photometric vector, per-band fluxes, and the three-matrix structure.
- Add the step-size convergence plot and the Abrams validation to an appendix.

**Style requirements for the whitepaper prose:** methodology at the level of physical reasoning, not code description — a reader should understand *why* the simulation does something, not just *what function is called*. Flag open items explicitly with "pending" or "known issue" language rather than silently omitting them. Keep cross-references between sections tight and bidirectional.

Then reconcile `OPEN_ITEMS.md` — anything resolved gets removed, anything still open gets carried into the whitepaper's open-items discussion.

---

## Deferred — do not address unless I ask

These are known open items being resolved in a deliberate order. Note them if you encounter them; do not fix them.

- Tile geometry details for the whitepaper Locations section
- The Roman astrometric error model decision (F146 noise constants: the ~100 mas FWHM replacing the current 20 mas, and the γ value for the F146 ≈ 22 transition, are known from the literature but not yet transcribed into `Bulge.h`)
- Whether `BulgeBaseline.dat` is a merged visit list
- Remaining TODO-flagged variables from the cadence refactor
- Cross-model benchmarking against genulens / SynthPop / MaBμlS-2

---

## A closing note on expectations

I expect the short-`tE` joint gain to be small and the long-`tE` gain to be large. **That is the hypothesis, not the desired outcome.** If the numbers come out otherwise, report them as they are. And if the short-`tE` ratio comes out at essentially 1.0, that is a result worth stating precisely — nobody has ever put a number on how much Rubin helps Roman for a ten-day event, and "essentially nothing during seasons, but it recovers X% of them from the gaps" is a clean, citable statement that makes the long-`tE` claim more credible, not less.

Never present a ratio above 1.0 as a finding. Adding data cannot decrease Fisher information; such a result is a bug, and we stop and find it.
