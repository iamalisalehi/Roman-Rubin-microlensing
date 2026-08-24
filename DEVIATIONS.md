# DEVIATIONS.md

Where the implementation departed from `JOINT_FIT_REFACTOR_PLAN.md`, and why.

The plan's original text is deliberately **left intact** — it records what was believed when it was
written, which is worth preserving. This file is the single place to look for what actually
happened instead. Each entry gives: what the plan said, what was done, why, and the commit.

Entries are listed in the order the work happened.

---

## 1. Unplanned bug fix: `fb` perturbation step binned off the wrong telescope

**Commit:** `c838c96` (inserted during Step C1)

**Not in the plan.** `FisherM` chose the finite-difference step for the blend fraction by binning
against `s.fb[0]` (Rubin) unconditionally, then applied that step to `s.fb[tt]` — whichever
telescope produced the epoch. Roman's F146 blend fraction is routinely much higher (its PSF is
~10x smaller, per Step B3), so a step sized for a low Rubin `fb` pushed Roman's past the physical
bound of 1.0 and tripped `CHECK(s.fb[tt] <= 1.0)`, an uncaught `std::runtime_error` that
hard-crashed the program.

**Why it was fixed immediately rather than deferred:** it made *any* numerical verification of
*any* Fisher-matrix step impossible. Confirmed pre-existing by reproducing it on the unmodified
branch tip. Surfaced to the user, who chose to fix it as its own step rather than fold it into C1.

---

## 2. Step C1: `t0` appended at index 5 rather than inserted second

**Commit:** `5f39836`

**Plan said:** extend the photometric parameter vector to `{u0, t0, tE, fb, piE, xi}` — i.e. `t0`
inserted in second position.

**Done instead:** `{u0, tE, fb, piE, xi, t0}` — `t0` appended at index 5.

**Why:** the plan's ordering was illustrative, but adopting it literally would have shifted the
index of every parameter after `u0`. Those indices are hard-coded in roughly six places through
`co.resu[]` — the output printers, the per-field aggregation, `corr1`, and `f1`/`f2` — where the
photometric entries `resu[0..4]` sit contiguously against astrometric entries `resu[5..8]` and
derived quantities `resu[9..14]`. Reflowing all of that was an invasive, error-prone change with
no physical benefit: a Fisher matrix does not care what order its parameters are listed in.

**Consequence to be aware of:** the plan text and the code disagree on parameter layout. The code
is authoritative. Anything reading `Delta1[]`, `Era[]` or matrix rows by index must use the code's
order.

---

## 3. Step C2 split into C2a and C2b, with C2b deferred

**Commit:** `639309a` (C2a). C2b **not yet implemented.**

**Plan said:** Step C2 gives each band its own source flux and blend flux as free Fisher
parameters, and stores the model magnitude in its native band.

**Done instead:** the step was split. C2a — an unplanned preparatory refactor — replaced three
separate hard-coded uses of LSST's r-band with a configurable `RUBIN_REF_BANDS` list, defaulting
to `{2}` so it is a verified no-op. C2b, the actual plan content, was deferred.

**Why the split:** the user asked for the Rubin reference band to be selectable (any single LSST
band, or a custom combination) rather than fixed to r. That is a self-contained change with an
exact acceptance test — byte-identical output under the default — so it was worth landing on its
own before the much larger parameter-set change.

**C2b since implemented.** `Nx` 6 -> 9, layout
`{u0, tE, fb0, piE, xi, t0, mbs0, fb1, mbs1}`. Each telescope's source-flux fraction and baseline
magnitude are now free parameters affecting only that telescope's epochs. In these coordinates
`(mbs, fb)` is a bijective reparametrization of the plan's `(source flux, blend flux)`:
`F_src = fb * 10^(-0.4 mbs)`, `F_bl = (1-fb) * 10^(-0.4 mbs)`.

Because a single-survey partition carries no information about the other telescope's flux
parameters, each partition inverts only the submatrix it can constrain (`activePhotParams`). That
set is also epoch-aware: an event with no Roman epochs drops `fb1`/`mbs1` from the *joint* set too,
without which the joint matrix went singular on exactly the gap-peaking events the project is about.

The cost the plan asked to be explicit about: derivative evaluations per epoch go 54 -> 108
(exactly 2x), measured fixture runtime 1.81 s -> 4.47 s (2.5x), and conditioning worsens by roughly
an order of magnitude in the fixture, far more on some live events (joint condition number 612 ->
4.5e6 on a 936-day event). Step C4 landing first is what makes that safe to absorb.

**The scientific consequence was larger than expected and is worth flagging loudly.** Marginalizing
over the flux parameters inflates sigma(tE) by 12-680% depending on event, and substantially
*reduces* the measured joint gain: joint/Roman-alone for a long in-season event went 0.723 -> 0.993,
and joint/Rubin-alone for a long gap-peaking event went 0.212 -> 0.979. The earlier, larger gains
were inflated by holding baseline flux and blend fraction fixed at their true values, which
over-credits data that only measures baseline — for a gap-peaking event Roman's epochs see a flat
light curve, which constrains `mbs1` but says almost nothing about `fb1`. The gap-filling claim
survives strongly (joint/Roman = 0.0065 for `long_ingap`); the symmetric "Roman helps Rubin" and
in-season joint gains largely do not.

**Follow-up fix (same step).** The first C2b commit (`385a6c1`) left `fb0` perturbing `s.fb[tt]`
rather than `s.fb[0]` — correct before the split, when index 2 was a single telescope-selected
blend fraction, but wrong once index 7 became Roman's own. On Roman epochs parameters 2 and 7 both
moved `s.fb[1]`, so the joint matrix carried a duplicated direction and reported
`sigma(fb1) = 1.14e-2` against Roman-alone's `6.25e-3` — a joint fit worse than a single survey,
which is impossible.

That commit's message and an earlier version of the fixture comment both explained the violation as
a legitimate consequence of the partitions having different parameter sets. **That explanation was
wrong.** By the Schur complement the inequality is a theorem: partition the joint parameters into A
(the single survey's active set) and B (the other survey's flux parameters); the other survey
contributes nothing to B, so the joint fit's effective information on A is
`F_thisSurvey[A,A] + Schur(F_other)`, and a Schur complement of a positive semi-definite matrix is
positive semi-definite. So `sigma_joint <= sigma_single` for every parameter the survey constrains,
including its own flux parameters. A violation is always a bug.

Fixed, and the fixture assertion restored to full strength rather than restricted to the shared
geometric parameters. The corruption was confined to the flux parameters — every `sigma(tE)` ratio
is byte-identical before and after — so the joint-gain findings above are unaffected.

**Still outstanding from Step C2:** model magnitudes are not stored in their native band. Every
Rubin epoch, whatever filter it was actually taken in, still collapses to the single representative
band chosen by `RUBIN_REF_BANDS`, so the chromatic information across `ugrizy` never reaches the
Fisher matrix. This is the one part of the plan's Step C2 not done.

---

## 4. Unplanned Step C0: Fisher matrices were never accumulated over epochs

**Commit:** `d5c8867`. Plan text corrected inline at Step C5 by `1efe919`.

**Not in the plan** — and it contradicted the plan's stated premise. Step C5's original text said
"`FisherM` accumulates `for (int i = 0; i < ndw; ++i)` over **all** data points." It did not. Both
accumulation sites used `gsl_matrix_set` (overwrite) rather than a running sum, and the zeroing
loops run once *before* the data loop, so every epoch clobbered the previous one and the final
matrix held a single data point's contribution — rank 1, hence singular, hence `invert_matrix`
nudged the diagonal by `1e-10` and returned an inverse of order `1e10`.

**Every sigma the code produced before this fix was meaningless.** A real detected event reported a
relative sigma of 1.8e7 on `tE`. After the fix the same event reports 0.0030.

**Why it went ahead of everything else:** Phase C is entirely about making precision numbers
trustworthy, and every later step builds on this one. In particular Step C5 partitions this
accumulation — partitioning a broken sum would have meant nothing.

**Related correction:** the degenerate synthetic-event result seen during Step C1 was recorded at
the time as a badly conditioned test fixture. That diagnosis was wrong; it was this bug. Corrected
in `OPEN_ITEMS.md`.

---

## 5. Unplanned: a synthetic Fisher-matrix fixture, built before Step C5

**Commit:** `7419e17` (`tests/fisher_fixture.cpp`, `make fishertest`)

**Not in the plan.** The plan assumes numerical acceptance criteria can be checked against the live
Monte Carlo. In practice detection efficiency in the test field measured 1 event per 886, 2600 and
10735 star draws on three separate runs, each taking tens of minutes. Step C1's acceptance
criterion ("on a handful of test events") ended up resting on a single event.

**Why it was built before C5 rather than after:** C3's step-size sweep and C5's per-event ratio
distributions both need many events, not one. The user chose fixture-first when offered the
alternative of widening the sky field.

**Design constraint worth remembering:** the fixture reads **no data files**, because `Baseline/`,
`CMD/` and the extinction maps are gitignored and absent after a fresh clone. It is a regression
harness, not a validation of absolute precision — events are hand-picked rather than
population-weighted and the error model is a flat stand-in for `errlsstM`/`errRomanM`.

---

## 6. Step C5 moved ahead of Step C2b

**Commit:** `1410cec`

**Plan order:** C2 → C3 → C4 → C5.

**Done instead:** C5 before C2b.

**Why:** the user objected to the C2b design on scientific grounds — they want the measurement of
how much each survey *helps* the other, and had no interest in Fisher results for events only one
telescope saw. Investigating that showed C5 subsumes most of C2b: once accumulation is partitioned
by `tele[i]`, `F_rubin` perturbs only `s.fb[0]` and `F_roman` only `s.fb[1]`, so both already have
unambiguous per-survey blend-fraction semantics. Only `F_joint` genuinely pools the two. C5 was
therefore both the higher-value step and the one that makes C2b smaller and verifiable.

**Decision recorded at the same time:** the joint/single gain ratio is computed only for events
where both surveys have data (`ndw_L > 0 && ndw_R > 0`); zero-epoch cases are counted as a
category, per the plan's Step F1 "three currencies" rule.

---

## 7. Part of Step C4 pulled forward into Step C5

**Commit:** `1410cec`

**Plan said:** Step C4 replaces "any crash-or-garbage behaviour on singular matrices with an
explicit *not characterizable* outcome."

**Done early:** `invert_matrix` now returns success/failure instead of nudging a singular matrix's
diagonal by `1e-10`, and partitions with fewer epochs than free parameters are rejected up front.
Invalid partitions report sigma `-1` as an explicit sentinel.

**Why:** C5 creates single-survey partitions that are *routinely* empty or rank-deficient — a short
event peaking in a Roman gap genuinely has no Roman data. Without this, C5's own output would have
been garbage on exactly the events that matter most.

**Since completed.** Step C4 was taken next and finished the remaining work: symmetric
normalization before inversion (`F~ = D F D` with `D = diag(1/sqrt(F_ii))`, inverted, then
rescaled — algebraically identical, numerically far better behaved), per-event condition numbers
computed from the normalized matrix and stored in the output, and rejection above a documented
threshold. All sigmas came out bit-identical to the pre-C4 values, confirming the normalization is
the no-op it should be.

**Caveat worth carrying forward:** on the current event set the `kMaxCondition` threshold has never
actually fired — every rejection so far came from the zero-epoch or exactly-singular checks that
already existed. Its behaviour on a genuinely near-singular matrix is therefore still unverified.
A fixture event deliberately constructed to sit near the threshold would close that gap.

**Related fix in the same commit:** `invert_matrix` was destroying its input, because
`gsl_linalg_LU_decomp` overwrites in place. It now works on a scratch copy. This was found by the
fixture's `F[joint] == F[rubin] + F[roman]` assertion and would otherwise have silently broken
C4's condition-number work, which needs the undecomposed matrix.

---

## 8. Unplanned: explicit synergy classification

**Not in the plan.** Added at the user's request after Step C5.

Events where one survey contributes real information but cannot characterize the event alone are
the strongest evidence for the joint fit. The live 936-day event is the motivating case: Roman
*detected* it and contributed 760 epochs, cannot characterize it alone (`okA[SROMAN] = 0`), yet
adding its data cuts sigma(tE) from Rubin's 4.41 d to 2.79 d.

A naive analysis computing `sigma_joint / sigma_roman` hits the not-characterizable sentinel on
exactly these events, and dropping the row would discard the best synergy cases — the mirror image
of the selection bias the plan warns about at Step C5. `SynergyClass` (`Bulge.h`) labels every
event `none` / `both-alone` / `rubin-only-alone` / `roman-only-alone` / `joint-only`, the last
being a pure joint-fit rescue where neither survey alone suffices. Recorded per event in
`EventRecord` and the output.

**Rule for downstream analysis:** never drop a row because a single-survey sigma is missing.
Classify it.

---

## 9. Step C1's acceptance criterion, closed late

**Commits:** `5f39836` (change), then closed properly once the fixture existed.

The plan required sigma(`tE`) to increase or stay equal on "a handful of test events." At the time
C1 was committed this could only be checked on **one** event, for the detection-efficiency reasons
in entry 5. It was closed later across all fixture events and every survey partition.

The result quantified the specific bias the plan predicted — that omitting `t0` "systematically
flatters the Rubin-alone column." Rubin's sigma(`tE`) was understated by up to a factor of 3.2,
while Roman's densely sampled partition barely moved (ratios ~1.00-1.02).

---

## 10. Step C3: the sweep window in the plan was in the wrong place, and the result was not a confirmation

**Commit:** this step. **Plot:** `c3_step_sweep.png`, from `c3_sweep.csv` via `tests/c3_step_sweep.py`.

**Plan said:** sweep each `Delta1[]` entry "over about a decade", expect a plateau, pick steps in
the middle of it, and confirm sigma is stable to well under 1%. The tone of the step is that this
is a check the code is expected to pass.

**What happened instead:** over a decade either side of the legacy steps there is no plateau for
any nonlinear parameter — `sigma(u0)` varies by four orders of magnitude across that window. The
window was in the wrong place. The legacy steps are ~25% of the parameter value (`u0` perturbed by
0.15 on a `u0` of ~0.3; `tE` and `t0` by `0.25*tE`), which is enormous for a derivative step, so
the entire assumed window sits up the truncation-error branch.

Widening the sweep to nine decades (`1e-8` to `10` times the legacy steps) exposes both walls of
the expected U and a very wide plateau: sigma is flat to **<0.2%** over `1e-6..1e-3` of the legacy
values, with round-off taking over below and truncation above. So the plan's acceptance criterion
is met — but only at steps four orders of magnitude smaller than the ones in use.

**What this cost.** The legacy steps were not merely imprecise, they were biased, and unevenly:

| parameter | error in sigma at the legacy step (median / worst over the fixture) |
|---|---|
| `u0`  | 57% / 110% |
| `t0`  | 71% / 92% |
| `tE`  | 26% / 55% |
| `piE` | 9.6% / 93% |
| `xi`, `fb0`, `fb1` | 0.05–0.3% (worst 3.5%) |
| `mbs0`, `mbs1` | 0.00% — linear in the model, exactly flat at every step, the sweep's control |

**Where the error lived, and why it matters more than the percentages suggest.** Diagonalizing the
normalized Fisher matrix (`./fishertest --eigen`) shows the six best-constrained eigenvalues are
essentially unchanged by the step correction. The entire discrepancy sits in the flattest
directions: for `short_inseason` the smallest normalized eigenvalue was inflated by a factor of
~4e7, and the condition number fell from a healthy-looking 3.4e2 to its true 1.5e10. Since the
normalized matrix has unit diagonal its trace is fixed at the dimension, so this was information
*redistributed* into the degenerate directions, not created. In physical terms, truncation error
was papering over real parameter degeneracies — reporting a 0.2% parallax measurement for a 5-day
event, which cannot measure annual parallax at all.

**Consequence to be aware of: every absolute sigma printed before this step is too small**, badly
so for short events. This step invalidates absolute-precision numbers, it does not merely refine
them. Condition numbers rise correspondingly and are now genuine; `kMaxCondition = 1e12` is
consequently much closer to firing than DEVIATIONS entry 7 recorded, and short events now land
within about two decades of it.

**What survives, and this is the load-bearing part:** the joint-vs-single-survey *ratios* are
essentially unchanged, because Step C5's paired comparison cancels the common inflation.
`sigma(tE)` joint / best-single, legacy steps → corrected steps:

| event | legacy | corrected | | event | legacy | corrected |
|---|---|---|---|---|---|---|
| short_inseason | 0.9995 | 0.9981 | | long_inseason | 0.9935 | 0.9717 |
| short_ingap | 1.0000 | 1.0000 | | long_ingap | 0.9788 | 0.9827 |
| mid_inseason | 0.9980 | 0.9902 | | verylong | 0.7588 | 0.7951 |
| mid_ingap | 1.0000 | 1.0000 | | | | |

Same ordering, same physics, gains slightly larger in most cases. The two events where the gain
*shrank* (`long_ingap`, `verylong`) both have the correction hitting the *helper* survey harder
than the leader: in `long_ingap` Roman is a weak helper (peaks in a Roman gap) and its `sigma(tE)`
inflated 2.14x against Rubin's 1.28x; in `verylong` Rubin's contribution is the parallax channel
and its `sigma(piE)` inflated 2.18x against Roman's 1.10x. Read that mechanism as a hypothesis
from seven fixture events, not an established property — it wants re-testing on live events.

**Implementation choice:** a single `kFDStepScale = 1.0e-4` multiplying `Delta1[]` (and the
telescope-keyed `bb[]` blend-fraction steps), rather than five rewritten constants. This keeps the
change auditable against the legacy values and keeps the sweep — which is defined in units of the
production step — directly comparable. `1e-5` reproduces every sigma to 0.3%, so the choice is not
delicate.

## 11. Unplanned: `tE` and `piE` were on a first-order, biased finite-difference stencil

**Commit:** same as entry 10. **Not in the plan** — found while doing C3's sweep.

`Bulge.h` defines two stencils. `sig = {+1,-1}` averages a forward and a backward difference,
giving a **central difference, accurate to O(h^2)**. `sig2 = {+0.5,+1.0}` averages two *forward*
differences (at `h/2` and `h`); nothing cancels, and the result is `f'(x) + (3/8) f''(x) h + O(h^2)`
— **first order and biased**. Verified numerically on a smooth test function: `sig` gains 2.00
decades of accuracy per decade of step reduction, `sig2` gains 1.00, and `sig2` carries ~22x the
error at equal step.

`sig2` was applied to `tE` and `piE` in the photometric matrix — two of the three parameters the
thesis novelty claims rest on — and to `tetE` and `piE` in the astrometric one. The pattern is
exactly the strictly-positive parameters, which suggests the legacy intent was to avoid stepping
them through zero. That is not a real constraint: the steps are a small fraction of the parameter,
and after entry 10 they are four orders of magnitude smaller still.

**Done:** the photometric `tE` and `piE` now use `sig`. The astrometric uses are left alone pending
a `Delta2[]` sweep (logged in `OPEN_ITEMS.md`) — the astrometric matrix has its own unswept steps
and a placeholder Roman error model, so fixing its stencil in isolation would be a change whose
effect could not be verified.

**Numerical effect at the corrected step size: none measurable** — every fixture sigma and every
joint/single ratio is unchanged to four significant figures. This is expected, since the bias term
scales with `h` and `h` is now tiny. It is a correctness fix that removes a trap: it means the
`tE` and `piE` derivatives no longer degrade faster than the others if steps are ever revisited.

---

## 12. Unplanned bug fix: `Fluxb` accumulated blend flux across every star ever drawn

**Not in the plan.** Found while running the Step C3 live verification; present since the initial
commit (`e40716a`), inherited from the legacy LMC code.

`func_source` (`Lensing.cpp`) builds a source star plus the unresolved neighbours blended into the
same seeing disc, summing their fluxes into `s.Fluxb[i]` with `+=`. The neighbour count `s.nsbl[i]`
is recomputed at the top of every call. **`s.Fluxb[i]` was never reset.** It was zero-initialised
once at construction and thereafter only ever incremented, so every star inherited the blend flux
of every star drawn before it.

**What that does.** Two derived quantities are built from it:

    s.magb[i]  = -2.5 * log10(s.Fluxb[i]);                 // blended baseline magnitude
    s.blend[i] = pow(10, -0.4 * s.Map[i]) / s.Fluxb[i];     // source's share of the blend

so `magb` brightened without limit as `magb(N) ~ magb(1) - 2.5*log10(N)`, and the blend fraction
decayed as `1/N`. Since a microlensing event magnifies only the source, `A_obs = A*fb + (1 - fb)`,
a collapsing `fb` washes the event out entirely.

**Measured on the 300k-draw run that exposed it** (preserved at
`runs/test2_c3run_STALLED_20260821.dat`): `fb` fell from 0.130 on the first star to 0.0017 on the
*second*; mean blended magnitude drifted 17.6 -> 3.5 across the run, brighter than Sirius. Back-
extrapolating `magb + 2.5*log10(N)` from any block recovers `magb(1) ~ 17.2-17.5` against a measured
17.635 on row 1 -- constant across a 29x range in N, which is proof of clean unbroken accumulation.

By ~3000 stars `magb` passes the saturation guard `magb[i] > satu[i]` and `blend -> 0` drives the
acceptance draw `testL <= blend[2]` to zero. **Detection then becomes impossible, permanently.**
That run produced 9 detectable events in its first 2965 draws and then none at all in the following
290,000, while still consuming 100% CPU.

**Consequences, which reach further than the stall:**

- Only the **first** star of any run ever had correct blending.
- **No run could ever satisfy the 20-detection stopping criterion**, so `./roman` never terminated
  on its own. Every previous "completed" run was a manual kill. The pre-C3 baseline stopping at
  2297 rows with `icon = 8` is exactly this.
- Every detection-efficiency figure this project has quoted (1 per 886 / 2600 / 10735 draws) was
  measuring the bug, not the astronomy.
- `fb` and `mbs` were wrong on every live Fisher-scored event past the first.

**The fix** is one line -- zero `Fluxb` in the loop that already resets `nsbl`, where its absence
was easy to miss.

**Verified.** `./fishertest` is bit-identical before and after, confirming the Fisher path is
untouched (the fixture hand-builds its events and never calls `func_source`, which is also why
every Step C3 result stands). After the fix `magb` scatters around 18.1 with no trend and `fb`
varies 0.0005-0.42 with crowding instead of decaying. A short run reached the per-field stopping
criterion in **24 draws** -- 24 usable light curves, 5 detections -- the first time any field has
ever completed.

**Not affected:** Step C3, the fixture, and the whitepaper appendix, all of which are fixture-based.
**Affected and needing re-measurement:** every live yield, efficiency and per-field aggregate.

**Immediately exposed a second, previously unreachable bug** -- see entry 13.

## 13. `nstE` / `ndtE` are never filled, so the efficiency and event-rate outputs are dead

**Not in the plan.** Surfaced the moment entry 12's fix let a field complete for the first time.

`FunctE` (`helper.cpp`) returns the `tE`-histogram bin `gg` for each event, and `gg` is stored in
`EventRecord`. But **nothing ever increments `l->nstE[gg]` or `l->ndtE[gg]`** -- the counters for
"simulated in this `tE` bin" and "detected in this `tE` bin". They are zeroed per field
(`Bulge_LSST.cpp:338`) and then only read, at the per-field aggregation (`849-851`) and at

    EFF += l->ndtE[gg] / (l->tE / year);   // 1/years

So `EFF` sums zeros. Everything downstream of it is identically zero too: the detection efficiency
`EFF`, the event rate `Gamma = 2/pi * opd * EFF / u0m`, and the expected event count `Neven`. The
run then aborts on `CHECK(EFF > 0.0)`. The neighbouring histogram calls (`FuncMl`, `FuncPi`,
`Funcu0`, ...) are commented out at the same site, which suggests a block of histogram-filling code
was disabled or lost and the consumers were never updated.

**Why it was never seen:** this code lives after the per-field `do/while`, which entry 12's bug made
unreachable -- no field had ever finished.

**What's needed:** decide what each counter means and fill it. `nstE[gg]` should plausibly increment
for every simulated event in the bin and `ndtE[gg]` for every detected one -- but "detected" now has
three meanings (`detL`, `detR`, `detJ`), so this needs a decision rather than a guess, and it is
really the per-survey efficiency question of Phase F. Until then the per-field efficiency, event
rate and yield outputs cannot be trusted and `./roman` cannot complete a field.

---

## 14. New plan section inserted between Phases C and D: detection taxonomy and run statistics

**Added at the user's request.** The plan goes straight from Phase C (the Fisher matrices) to Phase
D (the output table). Fixing entry 12 let a field complete for the first time, which made the
per-field aggregation code reachable and exposed a cluster of problems that belong together and
are not Phase C work. They are grouped here as **Phase C-D**.

### C-D.1 A detection taxonomy, replacing the old five-way label

The user's framing: an event is only meaningfully detected if the **joint** fit detects it, because
the joint stream contains strictly more data than either survey alone. That leaves four ways an
event can be detected, distinguished by which telescopes *also* detect it unaided:

| class | meaning |
|---|---|
| `DET_JOINT_ONLY`  | neither telescope alone, but the combination does — **the most interesting class, and expected to be rare** |
| `DET_RUBIN_JOINT` | Rubin alone, and the joint fit |
| `DET_ROMAN_JOINT` | Roman alone, and the joint fit |
| `DET_BOTH_JOINT`  | both telescopes alone, and the joint fit — more interesting than either single class |

The old labelling (`nDetBoth` / `nDetRubinOnly` / `nDetRomanOnly` / `nDetJointOnly` / `nDetNeither`)
ignored `detJ` except in its last branch, so "Rubin-only" meant "Rubin detected", regardless of
whether the joint fit did. Replaced by `DetClass` (`Bulge.h`) and `nDetClass[]`.

### C-D.2 The joint detection test was not monotone in the data

`DET_ANOMALY` exists to catch the case that should be impossible: a single telescope detecting an
event the joint test misses. **It occurred immediately** — 1 of 5 detections in the first field.

The cause is a threshold form, not a coding slip. All three tests compare `dchi` against
`2 * ndw`, i.e. they threshold the **mean per-epoch** chi-squared improvement. Since
`dchi_joint ~ dchi_L + dchi_R` and `ndw = ndw_L + ndw_R`, pooling a survey with many low-signal
epochs raises the joint bar without contributing signal. Observed live: a `tE = 631` d event
cleared Roman's bar on its 760 epochs, then failed the joint test because Rubin's 989 near-flat
epochs lifted the joint threshold by ~1978 while adding almost no `dchi`.

This contradicts a basic principle — conditioning on more data cannot destroy information — so
`detJ` is now made monotone: if either survey alone detects, the joint detects. `detJ_raw`
preserves the unmodified test so the rate of the inconsistency stays measurable as `DET_ANOMALY`.

**Left open, deliberately:** the `2 * ndw` threshold form itself. A chi-squared detection statistic
should be thresholded on its total, not its mean, so that more data helps rather than hurts.
Changing it moves `detL` and `detR` too, which makes it a science decision about detection
criteria rather than a bug fix. It belongs with Phase F's yield work.

### C-D.3 Filling `nstE` / `ndtE` (closes entry 13)

`gg = FunctE(*l)` is now evaluated for **every** simulated draw rather than only inside the
detection branch, and `nstE[gg]` incremented there, so the efficiency denominator counts
everything simulated. `ndtE[gg]` is incremented on `detJ`, the taxonomy's definition of detected.
`EFF` went from identically zero to `0.2016` on the first field, and `Gamma` and `Neven` with it.

### C-D.4 The not-characterizable sentinel was poisoning the precision averages

The per-field means `Eru0`, `ErtE`, `ErpiE` and the rest summed `co.resu[]` over every event that
reached `FisherM`. Step C4 reports sigma `-1` for a partition it cannot characterize, and
`ErrorCal` divides that by the parameter value, so one such event contributed `resu[3] = -697` and
dragged the mean fractional `piE` error negative, tripping `CHECK(ErpiE > 0.0)`.

The means are now taken only over events whose joint partition is characterizable, with `nErAvg`
tracking the true denominator, and the `CHECK`s apply only when `nErAvg > 0` — a field with no
characterizable event has no mean precision, which is a legitimate outcome rather than a crash.

**This is not the selection bias DEVIATIONS entry 8 warns about.** That rule governs the
joint-vs-single *ratio* statistics, where a missing single-survey sigma is itself the result and
dropping the row would discard the best synergy cases. Here we are forming a *mean precision*, and
an event with no measurement has no precision to contribute; including the sentinel would be
averaging a placeholder.

A related replay bug was fixed in the same place: the per-field aggregation loop restores each
event's `resu[]` from its `EventRecord`, but `co->okA[SJOINT]` was not restored, so it held
whatever the last `FisherM` call had left. `okJoint` is now replayed with the rest.

### C-D.5 Run-wide statistics

Per-field counts are too small to quote for a rare class, so `DetClass` counts are also accumulated
run-wide and broken down by `tE` bin, printed at the end of the run with `sqrt(N)` Poisson
uncertainties and as a percentage of all simulated events. The class is additionally persisted per
event (`EventRecord::detCls`, and a new output column) so any cut can be made offline without
rerunning. The `tE` breakdown matters because the entire science case is that the joint gain is
`tE`-dependent.

**Acceptance for this section:** a run completes without aborting, `DET_ANOMALY` is zero, and the
run totals give each class with its Poisson error — including an honest number, possibly zero, for
`DET_JOINT_ONLY`.

---

## 15. Step D2 found the survey data never reached the simulation

**Step D2's own question is answered and clean.** `t0` is drawn uniformly on `[2, Tobs-2]`
(`Lensing.cpp:312`), the full ten-year window, with no reference to any visit list. Both baselines
share its clock (day 0 = first Rubin bulge visit, MJD 61141.312 = 2026-04-11), the light-curve loop
integrates continuously over `[-100, 3752]` d rather than hopping between observing windows, and
the pre-selection gate is peak-magnitude-only with no time dependence. No gap-peaking event is
discarded for *when* it peaks.

The investigation instead surfaced two defects that made that correctness moot.

### 15.1 `ct` was capped at 1000, silently truncating both surveys

`matchVisibleEpochs` records a sightline's matching visit indices into `ct`, allocated flat at
`1000` for both instruments, with the cap hardcoded twice (`for (i < 1000)` reset, `if (ndd >= 999)
break`). The break was silent -- no warning, no `CHECK`. The schedule was therefore truncated to
its first 999 entries *in time order*, i.e. **the survey stopped early**:

| survey | visits on sightline | kept | stream ended | fraction of mission |
|---|---:|---:|---|---:|
| Rubin | 2,412 | 999 | day 1,572.9 | 43% |
| Roman | 51,514 | 999 | day 8.4 | **0.5%** |

Roman's design is six 72-day high-cadence seasons at 12.1-minute sampling separated by ~111-day
gaps. 999 epochs at 12.1 minutes buys 8.4 days -- not one season, and not one gap. So the
simulation contained no Roman season structure at all: no in-season vs in-gap distinction, no
multi-season baseline for the annual parallax that converts `tE` into a lens mass, and no route to
an event Roman alone detects beyond day 8.

**This invalidates the detection numbers in entry 14 / commit `e47390a`.** `DET_JOINT_ONLY = 0` was
recorded there as "consistent with Poisson limits." It is not a Poisson statement: joint-only means
"neither survey alone, but the combination succeeds," and with Roman reduced to 8.4 days that class
cannot be populated by construction. The zero measured this bug. The C-D *code* is unaffected --
`fishertest` is bit-identical across this fix -- but every C-D *number* must be re-taken.

`ct` is now sized to each instrument's own visit count (`Nl` / `NlRoman`; 14.7 KB and 1.24 MB, once)
and the cap derives from `ct.size()`. It is now unreachable by construction, and kept only as a
guard that aborts loudly, since silent truncation is exactly what is being removed.

**Cost, measured not estimated:** 2,454 -> 40,048 integrator steps per event (16.3x), because
`dt = min(cade, cadeR)` must step at 12.1 minutes once Roman's real epochs exist. Wall clock rose
only 6.9x (7.2 -> 50 s/field), because the per-event buffer reset at `Bulge_LSST.cpp:384-387` clears
all `coun = 312,770` slots across 7 arrays -- 2.2M writes per event regardless of the ~2,000 epochs
actually used -- and that fixed cost currently dominates. **Deferred, not fixed:** clearing only the
`ndw` entries actually written is a large and easy speedup, and belongs with the run-scaling work.

### 15.2 `BulgeBaseline.dat` was doubled, and a short read zero-filled silently

`readbaselineBulge.py` opened the output in append mode and had been run twice, so the file was
`header + 3,686 visits + header + 3,686 visits`. `Nl = 7373` read through the embedded second
header; `operator>>` failed on `#ID`, set the stream's failbit, and every subsequent extraction
became a no-op leaving the row at its zero-initialised value -- 3,687 phantom visits at
`(l,b) = (0,0)`, `tim = 0`, `filter = u`, `sig5 = 0`.

Nothing caught it: zeros pass every `CHECK` in the read loop, since `(0,0)` lies inside the bulge
region, `tim = 0` inside `[0, Tobs]`, and `filter = 0` is a valid `u` index. The phantoms were
harmless only by accident -- `(0,0)` is 1.03 deg from the nearest simulated sightline, well inside
Rubin's 1.75 deg FoV, so they *did* match, and survived only because `matchVisibleEpochs` rejects
epochs not strictly later than the last kept one, and because the 999-cap fired first. Fixing 15.1
would have exposed them.

Fixed at the root (append -> truncate), file regenerated, `Nl` corrected to 3,686, and both baseline
reads now abort with a diagnostic on a failed extraction instead of zero-filling. The regenerated
file is byte-identical to the old file's first copy apart from one header word, so day 0 did not
move. Incidental: `readbaselineBulge.py` loaded `layout_7f_3.centers` from the wrong directory
(aborting the script) and never used the result; path corrected. Plot outputs moved from `./jpg/`
to `../pics/`, where they were already being kept by hand.

### 15.3 Decisions taken, for the steps that follow

- **`MISSION_START_DAY = 730` d, as a variable, not a constant.** Roman's mission currently sits
  flush at day 0 of Rubin's window. Rubin has already started and Roman has not, so Roman's 4.7-year
  baseline moves to days 730-2446, leaving Rubin-only stretches of ~2 yr before and ~3.3 yr after.
  Those are not wasted draws: they are the control arm for "how much does adding Roman help", and
  the long-`tE` events that span them are where the parallax baseline pays. The value is expected to
  change, so it stays a named parameter.
- **Run scaling: stride in space, never in time.** The ten-year window is kept intact -- shortening
  it would destroy the gap geometry Phase D exists to measure, and `t0` is a per-event draw so a
  shorter window buys no proportional speedup anyway. Sightlines are instead strided across the full
  production region rather than taken as one contiguous corner (which samples a single extinction
  column and one bulge density regime). The stride must hit Roman's six field centres: at
  `FoVRoman = 0.28` deg only 13 of the current 36 sightlines see Roman at all, so a careless stride
  yields a Rubin-only run. The event budget (`icon`/`nlens`/`nerr`) is the honest cost lever.
- **Rubin and Roman already observe on independent cadences** -- separate `matchVisibleEpochs` calls,
  FoVs, epoch lists, cursors and measured cadences. The single shared `dt` is the *integrator* step,
  not a cadence: it must be small enough not to step over either survey's next visit, which is what
  protects Roman's 12.1-minute sampling from being erased by Rubin's ~3-day one. No change needed.

**Also noted, not acted on:** plan Section 0.3 describes a "multi-year gap in the middle" of Roman's
mission. `RomanBaseline.dat` has no such gap -- the middle four of ten seasons are low-cadence
(3-day), per ROTAC 2025. The code matches ROTAC; the plan text describes an earlier design. Logged
for Phase G.

---

## 16. Roman survey model reconciled against STScI's published GBTDS design

Follows entry 15. `MISSION_START_DAY` became a real parameter, and three defects in the
Roman survey model were corrected against two STScI sources the user designated as citable
for the whitepaper:

- **Design page** — https://roman-docs.stsci.edu/roman-community-defined-surveys/galactic-bulge-time-domain-survey
- **First-two-years schedule** — https://roman-docs.stsci.edu/roman-community-defined-surveys/roman-observations-in-the-first-two-years-of-science-operations

### 16.1 `MISSION_START_DAY = 730`, and it is now a runtime parameter

Roman's mission previously sat flush at day 0 of Rubin's window, asserting that both surveys
start together. Rubin has started and Roman has not. Roman now occupies days **730 -> 2448**
(2028-04-10 -> 2032-12-23), leaving ~2.0 yr of Rubin-only baseline before and ~3.3 yr after.

Those Rubin-only stretches are **wanted, not waste**: they are the control arm for "how much
does adding Roman help", and the span over which long-`tE` events accumulate parallax
baseline. Since the value is expected to change it is overridable (`--mission-start DAYS`),
and the generator now refuses a value that would push the mission past `Tobs` rather than
letting the C++ read abort on `CHECK(ro->tim[i] <= Tobs)` with no explanation.

Day 0 of that clock is the earliest Rubin bulge visit, **MJD 61141.312 = 2026-04-11**, set by
`readbaselineBulge.py` subtracting it. The **real** GBTDS start (2027-02-11) is sim day
**306**, available as `--mission-start 306`; 730 is the user's deliberate choice, not an
approximation of it.

### 16.2 The season gaps are not a half-year split

`SEASON_PERIOD_DAYS = year/2` produced a uniform 110.6-day gap between all ten seasons. It
was an inferred number, flagged as ASSUMPTION 2 in the generator's own header. STScI's
published windows show the gaps **alternate**, because Roman's sun-angle constraint sets the
visibility windows rather than arithmetic:

| window | dates | span | gap before |
|---|---|---:|---:|
| high | 2027-02-11 -> 04-20 | 69 d | -- |
| high | 2027-08-15 -> 10-25 | 72 d | 117 d |
| high | 2028-02-11 -> 04-21 | 71 d | 109 d |
| low  | 2028-08-16 -> 10-24 | 70 d | 117 d |

Replaced by `SEASON_PATTERN = [(0.0, 69.0), (185.0, 72.0)]`, a spring/fall pattern repeated
annually: gaps now alternate 116.0 d and 108.2 d. **This matters for Phase D specifically** --
`dt_to_season_edge` is the variable the gap-filling result is plotted against, so a wrong gap
distribution biases the headline figure directly.

Two independent cross-checks from the design page confirm the pattern rather than merely
permitting it: it quotes "~70.5 days" allocated per high-cadence season, and (69+72)/2 = 70.5
exactly; and six such seasons total 423 d = 96.6% of its quoted "438 observing days", matching
its "~97% devoted to high cadence" figure. The schedule page independently confirms
`HIGH_CADENCE_SEASONS = {0,1,2,7,8,9}` -- its first three published windows are high-cadence
and the fourth is low.

### 16.3 Low cadence 3 d -> 5 d (conflict recorded, not resolved)

The design page states low-cadence seasons use "~1.5 hour observing units that are repeated
every five days". `LOW_CADENCE_DAYS = 5.0`.

**The conflict flagged as ASSUMPTION 1 is not resolved by this.** The ROTAC 2025 overguide
(arXiv:2505.10574) cited in the literature review says 3-day; the schedule page gives no
cadence at all. This is a choice between sources. The whitepaper must cite the STScI design
page for the 5-day figure, not ROTAC.

### 16.4 `FoVRoman` confused an area with a radius

`matchVisibleEpochs` tests `sqrt(dl^2 + db^2) <= FoVRoman`, i.e. treats it as a **radius in
degrees**, while its comment cited "~0.28 deg^2 total" -- an **area**. The design page quotes
**1.7 deg^2 over six WFI fields** = 0.2833 deg^2 per field, so the equal-area radius is
`sqrt(0.2833/pi) = 0.3003` deg. The old 0.28 implied `pi*0.28^2 = 0.2463` deg^2, ~13% small.

Small number, large effect: sightlines with any Roman coverage went from **13 of 36 (36%) to
21 of 36 (58%)**. This parameter decides whether a sightline sees Roman at all, so it gates
every joint-detection statistic.

### Net effect

`NlRoman` 309,084 -> **302,406** (high-cadence seasons now total 423 d rather than 432 d, and
low-cadence visits drop with the 5-day step). Per-sightline `ndd (Roman)` 51,514 -> 50,401.
`fishertest` bit-identical, so the Fisher path is untouched. Entry 15's warning still stands:
no detection statistic should be quoted until a run is re-taken on this corrected model.

---

## 17. The simulator had never been compiled with optimization

Not a deviation from the plan so much as a correction to every runtime number recorded in
it. Raised while starting the run-scaling work (the step that decides how much of the sky a
run can afford), because that step is entirely a question of cost per sightline.

### 17.1 What was found

`Makefile:3` read `CXXFLAGS = -g -Wall -Wextra -std=c++17`. There is no `-O` flag, and GCC's
default is `-O0`. The Monte Carlo, the light-curve integrator and `FisherM` have therefore
all been running unoptimized for the entire life of this project, including in CI.

Measured on the same seeded run, counting events completed in a fixed 300 s wall-clock window:

| build | events / 300 s | sightlines |
|---|---:|---:|
| `-O0` | 257 | 7 |
| `-O2` | **1280** | **25** |
| `-O3` | 1280 | 25 |

**5.0x.** `-O3` gives nothing further, so `-O2` is what was adopted.

### 17.2 Why this is safe, and why `-ffast-math` is not

GCC does not reassociate floating-point arithmetic unless explicitly told to, because FP
addition is not associative -- `(a+b)+c` and `a+(b+c)` can differ in the last bits. So `-O2`
is required to leave results unchanged, and the measurement agrees: `-O2` output is
bit-identical to `-O0` across `EfLMC2.dat` (562 lines), `EfLMC2B.dat` (568), `LpLMC2.dat`,
`test2.dat`, stdout, and the whole `fishertest` table.

`-ffast-math` would break exactly this guarantee. The chi-squared accumulators and the Fisher
matrix element sums are long reductions over thousands of epochs; permitting reassociation
there would change the forecast sigmas silently and unreproducibly. It must not be added.

`-g` is kept so `CHECK()` aborts remain debuggable, and `-O2` introduced no new compiler
warnings (the same six pre-existing `set but not used` lines).

### 17.3 Consequence for numbers already recorded

Every wall-clock figure in this document is an `-O0` figure and is now ~5x pessimistic. In
particular entry 15's cost measurements -- "7.2 s/field -> 50 s/field", "36-field stub run
~4.3 min -> ~30 min" -- should be read as roughly 10 s/field and ~6 min at `-O2`. The
*relative* 6.9x slowdown from the `ct` fix stands; only the absolute times move.

This materially changes what the run-scaling step is choosing between: a 5x larger sky
sample, or a 5x larger event budget, is available for the same wall clock that was assumed
when that step was scoped.

### 17.4 A wrong hypothesis, recorded

The per-event buffer clear (`Bulge_LSST.cpp:411`) was flagged in entry 15 as "currently
dominating runtime" and was fixed first on that basis. That claim was wrong. The arithmetic
looked convincing -- 306,092 slots x 7 arrays = 2.1M writes per draw to reset the ~2,000
actually used -- but ~17 MB of sequential memory traffic is a few ms against an event that
costs ~1.2 s, i.e. about 0.1%. Measured back-to-back the fix gave 253 vs 257 events per
300 s, inside run-to-run scatter.

The fix was kept anyway (commit `aad8d55`): it is provably equivalent, it removes real if
minor waste, and the prefix form preserves the clean-slate invariant that would expose a
future partial-write bug. But it is a hygiene change, not a speedup, and the initial
estimate should have been measured before it was asserted.

---

## 18. Step 4: the run now scans the survey region instead of one corner

Phase D's run-scaling step, sized against the cost model measured in entry 17.

### 18.1 What the run was

`Bulge_LSST.cpp` hardcoded `lon 0.5 -> 0.6`, `lat -1.0 -> -0.9` -- a 0.1 x 0.1 deg patch --
with the real production bounds commented out on the lines directly above. Every detection
number this pipeline has produced describes that one patch. The production region is ~66.4
deg^2, and stellar density, extinction and Roman coverage all vary strongly across it.

The event budget was likewise hardcoded (`icon < 20 or nlens < 5 or nerr < 1.0`) with the
production `850/150/2.0` commented out beside it. Both are cost levers, and needing a
recompile to move them is how a run ends up with no record of what produced it.

### 18.2 The cost model this was sized against

Measured at `-O2` on a timestamped run, per sightline at the stub budget:

| sightline type | n | mean | min | max |
|---|---:|---:|---:|---:|
| Roman-covered | 19 | **10.6 s** | 7.1 | 17.1 |
| Rubin-only | 5 | **0.5 s** | 0.4 | 0.6 |

A Roman-covered sightline costs **22x** a Rubin-only one -- ~50,400 Roman epochs against
Rubin's ~2,400, and the light-curve integrator scales with that.

The design consequence is the opposite of what was expected. Roman's six fields cover
`6*pi*0.3003^2 = 1.700 deg^2`, only **2.56%** of the scan region, which looked like an
argument for stratified sampling -- a coarse grid overall plus a dense one inside the fields.
It is not: the 97.4% of the region Roman never sees costs ~2% of the run, because being 22x
cheaper and ~40x more numerous very nearly cancel. Total cost tracks total sightline count, so
a **uniform stride** is the right design and no stratification is needed.

Full-region cost by stride, at the production budget:

| stride | step | sightlines | Roman-covered | production run |
|---:|---:|---:|---:|---:|
| 5 | 0.10 deg | 6,712 | 147 | 2.4 d |
| 8 | 0.16 deg | 2,660 | 52 | 21.9 h |
| **10** | **0.20 deg** | **1,706** | **39** | **14.7 h** |
| 15 | 0.30 deg | 778 | 15 | 6.4 h |
| 20 | 0.40 deg | 435 | 11 | 3.9 h |

**Stride 10 at the full budget was chosen** (Ali, 2026-08-24): the finest grid that still
finishes unattended overnight, at ~6-7 Roman-covered sightlines per GBTDS field.

### 18.3 What changed

`main()` takes `argc/argv` and a `RunConfig`: `--stride`, `--events`, `--lenses`, `--nerr`,
`--stub`, `--help`. Defaults are the production values above. Nothing about a run's size
requires editing source any more.

Two guards, because a stride that steps over Roman's fields is invisible after the fact:

1. A square grid of spacing `h` only reaches every point of a disk of radius `r` when
   `h <= r*sqrt(2)`. Strides giving more than `FoVRoman*sqrt(2) = 0.4247` deg are refused.
2. Stronger, and the one that would survive a change to the field list: the six distinct
   pointings are read back out of `RomanBaseline.dat` and the grid is checked against each.
   Missing any is fatal. This fires correctly -- it caught `--stub`, which reaches only 2 of
   6, so it is advisory rather than fatal in stub mode.

Why this matters: a grid that misses a field produces zero Roman epochs there and reports
**joint-labelled columns built from Rubin data alone**, with no crash and no warning. Same
failure mode as the `ct` truncation in entry 15.

### 18.4 A silent off-by-one in the old loop

`for (lon = 0.5; lon <= 0.6; lon += dd)` lost its last row **and** last column. `dd = 0.02` is
not representable in binary, so after five additions the accumulator sits a few ulp above the
bound and the comparison fails. The stub grid advertised 36 sightlines and ran **25**. The
scan is now index-driven (`lonMin + i*gridStep`), which is exact and lets the sightline count
be computed before the run -- which the provenance block and the coverage guard both need.

This means `--stub` cannot byte-reproduce the old output: it now runs 36 sightlines, and the
11 extra ones shift the shared RNG stream, so everything after the first longitude column
differs. That is the correction landing, not a regression, and it was proven rather than
asserted: rebuilt with the stub bounds trimmed to the 25 points the old loop actually visited,
the new code is **bit-identical over all 1650 lines** of the pre-Step-4 run.

### 18.5 Provenance

Every run writes `files/MONTLMC/files/run_provenance.txt` (and echoes it to stdout) before any
science output: git commit with a `-dirty` marker, build time, stride, grid step, lon/lat
range, corner cut, sightline counts, Roman field coverage, all three budget targets, `Tobs`,
`Nl`, `NlRoman`, both fields of view, and the RNG seed.

The entry that is easiest to lose is `area_per_sightline`. `Neven` is a **surface density in
deg^-2** and nothing in the C++ sums sightlines into a survey-wide yield -- that happens
downstream, where each sightline must be weighted by the area it stands for. At stride N that
area is `(N*dd)^2`, not `dd^2`, so **a strided run aggregated as if unstrided is wrong by
N^2** -- a factor of 100 at the default. It is written out explicitly rather than left to be
inferred.

### 18.6 `srand` removed

`srand(time(0))` was called at the top of `main()`. Nothing in any of the four sources or the
fixture ever calls `rand()`; the RNG is the seeded `mt19937_64` in `Bulge.h`. Its only effect
was to make runs look clock-seeded, which is exactly the opposite of what 18.5 is for. Removed
at Ali's request as part of this step.

### 18.7 Verification

`fishertest` bit-identical (all of this lives inside `main()`, which the fixture excludes).
Six pre-existing warnings, none new. Default config reports 1,706 sightlines / 39
Roman-covered / 6 of 6 fields, matching the independent Python cost model exactly.

Entry 15's warning still stands, and now has a run worth taking: no detection statistic should
be quoted until a full-region run is taken on this model.

---

## 19. Step D1: the per-event row becomes the analysis table

Every figure in the paper is a cut on one flat table. Step D1 is where that table stops being a
debug dump and becomes the thing the analysis reads. Steps C4 and C5 had already delivered most
of it — per-survey epoch counts, three detection booleans, three sigmas each for `tE`/`piE`/`tetE`,
the `DetClass` and `SynergyClass` labels, three photometric condition numbers. Four things were
missing, one of which the headline result cannot be produced without.

### 19.1 Gap geometry: the independent variable of the gap-filling result

Roman can only observe the bulge when the Sun angle permits, so its visit list is a comb of
~70-day observing seasons separated by ~110-day gaps. The joint-fit science claim is a statement
about what happens to events peaking **in those gaps**, so where `t0` falls relative to the
seasons is the x-axis of the headline plot — and it has to be computed at simulation time,
because nothing downstream can reconstruct it from the row.

Two columns now carry it:

- **`dt_edge`** — signed days from `t0` to the nearest season boundary, **negative inside a
  season**, positive outside. Signed that way so the plot reads left to right: `x < 0` is "Roman
  was watching", `x > 0` is "Roman was not", and the joint-over-Roman gain should grow with `x`.
- **`t0zone`** — 0 in-season, 1 mid-mission gap, 2 outside the mission.

**Why `t0zone` is three-valued and not a boolean.** An event peaking at day 300 and one peaking at
day 1050 both have "no Roman data at `t0`", and they are not the same object at all. Day 1050 sits
between seasons 1 and 2: Roman brackets it with dense photometry 63 days before and 45 days after,
so a long-`tE` event's wings are measured even though its peak was missed. That is the case the
joint fit is meant to rescue. Day 300 is before Roman launches — Rubin-only by construction, with
nothing to rescue. Pooling them would dilute the measured gain with events that were never
candidates, which is the easiest available way to wash the effect out. On the current schedule the
split of simulated peaks is roughly **19% in-season, 28% mid-mission gap, 53% off-mission**, so
this is not a small correction to a rare category.

**The windows are derived, not restated.** `buildRomanSchedule` (`helper.cpp`) clusters the epoch
times in `RomanBaseline.dat`: a spacing larger than `SEASON_GAP_MIN_DAYS = 20` d starts a new
season. The schedule already lives in `Baseline/generateRomanBaseline.py`; a second copy in the
C++ would drift silently, and the schedule is expected to change (see the deferred GBTDS-footprint
item in `OPEN_ITEMS.md`). Deriving means the C++ can never disagree with the visit list it is
actually integrating.

The 20-day threshold has a wide margin — the largest spacing *inside* a season is 5.0 d (the
low-cadence seasons' five-day sampling) and the smallest gap *between* seasons is 108.2 d — but
the margin is **checked at runtime, not assumed**. If a future schedule ever samples a season more
sparsely than the threshold, or packs seasons closer together than it, the clustering fails while
`dt_edge` and `t0zone` still look perfectly reasonable. `main()` refuses to run in that case
rather than emit gap geometry that is quietly fiction, and the recovered season count and both
measured margins go into `run_provenance.txt`.

### 19.2 sigma(lens mass), per survey

`Ml` is never fitted. It follows from the two Einstein-radius observables,

    Ml = tetE / (kappa * piE),   kappa = 8.144 mas/Msun

and because that is a pure ratio the fractional errors add in quadrature. It is now computed for
all three Fisher partitions (`relMl_J/L/R`) rather than the joint only.

The reason it belongs in the table three times over is that its two ingredients come from
different instruments' different strengths. `tetE` comes from the **astrometric** matrix — the
sub-milliarcsecond centroid wobble, which is Roman's regime. `piE` comes from the **photometric**
matrix over a long time baseline — the annual parallax distortion of the light curve, which is
Rubin's. A mass measurement can therefore exist in the joint fit that exists in neither survey
alone. That is the black-hole result, and it is invisible unless the mass precision is stored per
survey.

### 19.3 A silent bug: sigma(piE) could come out negative

`ErrorCal` took the better of the two independent parallax routes with
`resu[3] = MIN(resu[3], resu[8])` — photometric (`Era[3]`) against astrometric (`Erb[3]`). But
`Erb[3]` is `-1.0` when the astrometric matrix is singular, and that `-1.0` is an explicit
"not characterizable" **sentinel**, not a measurement. `MIN` picked it in preference to a real
sigma every time.

This is the source of the `resu[3] = -697` event already recorded in entry 8, which at the time
was worked around by *excluding* such events from the field average rather than fixed. It
propagated into `resu[9]` (mass) and `resu[10]` (distance), so both were wrong on those events
too. Now fixed with the same sentinel-aware rule used for `relMl[]`: take the minimum only among
routes that actually produced a measurement, and report `-1` if neither did.

Fixed here rather than left open because D1 adds `relMl_J`, which is the same physical quantity as
`resu[9]`. Leaving the old one unguarded would have put two columns in the same row that claim to
be the fractional mass error and disagree with each other on exactly the events where the
astrometric fit failed.

### 19.4 A silent bug: the first event of every run was dropped from the table

`filg_in` was constructed **already open** (`std::ofstream filg_in(testf)`), and the per-event
write then called `.open()` on it. Calling `open()` on an already-open `ofstream` sets `failbit`
and does nothing, so that write was silently discarded. The matching `close()` cleared the way, so
every *subsequent* event wrote fine — one row lost per run, always the first.

Verified by reproduction rather than inferred:

    iter 0 fail=1     file contains:  row1
    iter 1 fail=0                     row2
    iter 2 fail=0

The in-memory `records` vector was never affected, so no aggregate statistic was ever wrong; only
the flat table lost a row. `filg_in` is now left closed and the file is truncated (and the header
written) in a scope of its own.

The open/append/close **per event** is deliberate and stays: it flushes each row to disk as it is
produced, so a 15-hour run that is interrupted keeps everything it had computed. At ~1.2 s of
physics per event the syscalls are not measurable.

### 19.5 The row, and a deviation from the plan's wording

The row goes from 57 to **88 columns**, all appended — never interleaved, since the positional
aggregate initialiser at the `push_back` site and every column index downstream depend on the
order. New: `t0`, `xi`, `lon`, `lat`, Roman's `mbs1`/`fb1`, the per-filter `magb_*`/`blend_*`
arrays, `relMl_{J,L,R}`, `okB_{J,L,R}`, `condB_{J,L,R}`, `dt_edge`, `t0zone`.

`okB` is not redundant with `okA`: a row can have `okRubin = 1` (photometry fine) while
`sigtetE_L = -1` (astrometry singular), and before this nothing in the row explained why. An
analysis script would read that `-1` as a measurement.

`mbs1`/`fb1` today duplicate `magb[6]`/`blend[6]`, and `mbs0`/`fb0` duplicate `magb[2]`/`blend[2]`,
because `RUBIN_BANDS` is `{r}`. They are kept separate because they are the quantities the Fisher
matrix actually fits (photometric parameters 6–8), and the duplication ends the moment a band is
added.

The file now carries a `#` header line naming every column, from `eventTableHeader()` in
`Bulge_LSST.cpp`. An unlabelled 88-column matrix is unusable six months later, and mis-numbering a
column by one produces a plausible plot of the wrong quantity.

**Deviation:** the plan (Step D1) says "write out as a single well-headed file per field." After
Step 4 that would mean **1,706 files**. Instead there is one run-level file with `lon`/`lat`
columns, so "per field" is a one-line filter in pandas.

**Not done:** the file is still `./test2.dat` in the repo root. The name is meaningless and the
location is wrong, but `tests/c3_live_compare.py` reads it and is being actively edited. Recorded
in `OPEN_ITEMS.md` instead of renamed mid-flight.

### 19.6 A stale-value bug the verification caught — in the new code

The first D1 run produced `relMl_J` values on rows where `okA_J = 0`, `okB_J = 0` and
`sigpiE_J = sigtetE_J = -1`. A mass precision on an event that has no parallax and no angular
Einstein radius is impossible, and the repeats made the cause obvious: rows 6 and 7 both carried
`0.023221`, rows 9–11 all carried `0.104189`. It was the previous *characterized* event's value
leaking forward.

`ErrorCal` is only reached for events that pass detection and produce an invertible matrix.
Everything else in `covarian` that survives across events is therefore explicitly reset to its
sentinel at the top of each draw — `okA`, `okB`, `nepochA`, `condA`, `condB`, `Era`, `Erb`. The new
`relMl` array was not added to that list, so it alone kept the last value written.

687 of 905 rows were affected. Fixed by adding `relMl` to the per-event reset. Worth recording
because the failure mode is silent and plausible-looking — every leaked value was a perfectly
reasonable mass precision, just belonging to a different event — and because the same discipline
must be applied to anything added to `covarian` in future: **if it is written inside `FisherM` or
`ErrorCal`, it must be reset in the per-event loop, or it becomes the previous event's answer.**

Two related things noticed and deliberately *not* changed here:

- `co->flagi` is not in that reset list either, so it reads `1` on uncharacterized rows. It is
  harmless today because the only consumer pairs it with `okA[SJOINT]`, which *is* reset — but the
  `flagi` column in the table is stale for those rows and must not be used as a characterizability
  flag. Use `okA_J`. Recorded in `OPEN_ITEMS.md`.
- `relMl_J` and `rel_Ml` (`resu[9]`) now agree on every row where the joint fit is characterizable,
  which is the check that confirms the sentinel fix in 19.3 and the new per-survey computation are
  the same calculation.

### 19.7 Verification

Paired stub run (`--stub --events 20 --lenses 5 --nerr 1`, 36 sightlines) against a binary rebuilt
from `81a6b04`:

- **Columns 1–57 bit-identical across all 909 paired rows.** Nothing that existed before D1 moved,
  including `rel_piE` and `rel_Ml` — the sentinel fix in 19.3 changes those columns only on events
  where the astrometric matrix is singular *and* the photometric parallax is measurable, which did
  not occur in this stub. The fix is therefore verified as no-regression here, not as exercised;
  the negative-`resu[3]` case is rare and known from a previous full run.
- **The post-D1 file has exactly one more row than the pre-D1 file**, and the pre-D1 file's first
  row is `icon = 2`. That is 19.4 directly: the dropped first event is back.
- **Header 88 names, data 88 columns.**
- **Gap geometry**: 21.8% in-season / 26.3% mid-mission gap / 52.0% off-mission, against the
  19/28/53 predicted from the season windows and a uniform `t0`. `dt_edge < 0` on exactly the
  `t0zone == 0` rows and nowhere else. In-season range [−35.7, −0.1] d, elsewhere [+0.04, +1200.6] d.
- **Season derivation** checked independently of the run, by a standalone harness linking
  `buildRomanSchedule` against `RomanBaseline.dat`: 10 seasons over days 730.000–2447.965, matching
  the Python analysis of the same file window-for-window, with `maxInSeasonSpacing = 5.0` and
  `minSeasonGap = 108.247` against the 20 d threshold.
- **Photometry columns**: all `blend_*` in (0, 1]; `fb1 == blend_F146`, `fb0 == blend_r`,
  `mbs1 == magb_F146`, confirming the filter/telescope index mapping is the one intended.
- **No stray negatives**: no value of `rel_piE` or `rel_Ml` is negative other than the `-1` sentinel.
- `fishertest` bit-identical. Six pre-existing warnings, none new.

### 19.8 Closing the two verification gaps — and a hole in the guard

Two things 19.7 could confirm only indirectly are now covered by data-free checks in
`tests/fisher_fixture.cpp`, so they stay covered:

**`checkSentinelDiscipline`** exercises the `ErrorCal` fix of 19.3 directly, on all four
combinations of photometric/astrometric characterizability, rather than waiting for a live event
in that state to turn up. The discriminating case is *photometry good, astrometry singular*:
`resu[3]` must come out **positive**. Under the old unguarded `MIN` it was negative, every time.
The check also asserts that with both routes available the *smaller* is taken (not merely one of
them), that a mass is reported only when both of its ingredients exist, and that `relMl[SJOINT]`
agrees with `resu[9]` exactly — the two are the same physical quantity computed independently.

Writing it exposed a trap worth recording: **`ErrorCal` does not read `Era`/`Erb`, it recomputes
them** from the diagonals of the inverse covariance matrices. Setting them by hand in a test is
silently overwritten. The inverses are what a test has to construct.

**`checkSeasonClustering`** builds synthetic schedules — no data files — and checks both halves of
the runtime guard: a healthy schedule (3 seasons, 70 d each, 5 d cadence, 110 d gaps) must cluster
correctly and pass, with the `dt_edge` sign convention and all three `t0zone` states verified on
probe points; a pathological one must be caught.

**That second case found a real hole in the guard.** A schedule whose in-season cadence *exceeds*
`SEASON_GAP_MIN_DAYS` makes every epoch look like a season boundary, so each "season" ends up
holding exactly one epoch. Both margins the guard checked looked healthy: `maxInSeasonSpacing`
stayed **0** — no spacing was ever classified as in-season, so it was never updated — and
`minSeasonGap` was 25 d, comfortably above the threshold. Nine single-epoch "seasons" sailed
through, and `dt_edge`/`t0zone` would have been computed off them without complaint.

Fixed by adding `minSeasonLength` (the shortest season's duration) to `RomanSchedule` and
requiring it to be strictly positive. A season holding one epoch has zero length, which is the
one signature that case cannot hide. The real schedule's shortest season is 65 d, so the margin
is wide. Both the guard and the fixture's copy of the predicate were updated together, and the
value is reported in `run_provenance.txt`.

The pre-existing fixture event table is bit-identical — the new checks are purely additive and
print above it.
