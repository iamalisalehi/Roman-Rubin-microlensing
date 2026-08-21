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
