# OPEN_ITEMS.md

Running list of known-but-deferred discrepancies surfaced while working through
`JOINT_FIT_REFACTOR_PLAN.md`. Items are removed once resolved; unresolved items are carried
into the whitepaper's open-items discussion at Phase G, Step G3.

---

## Roman/F146 astrometric error model (from Step B1)

Step B1 needed a per-epoch Roman astrometric error to build `chi1a_R/chi2a_R/chi3a_R`
(Roman-only astrometric χ² bookkeeping) and `dchiA_R`. No real Roman astrometric error
function or dataset exists yet, so `Bulge_LSST.cpp`'s Roman branch currently calls
`errlsstA(*ls, magni[fiR])` — LSST's astrometric-error curve, evaluated at Roman's own
F146 magnitude/epoch — as an explicit placeholder (search `errsR` for the `TODO(Ali)`
comment at the call site). This is the same placeholder value the pre-B1 code was
already implicitly using for `l->erra[]` on Roman epochs (just made deterministic and
epoch-correct instead of a stale leftover from LSST's last-executed epoch).

**What's needed:** a real `errRomanA()`-style function (mirroring `errlsstA()`) reading
a real F146 astrometric-error-vs-magnitude dataset (mirroring `sigmaA_LSST.txt` /
`lsst::mag,err`), using the noise constants named in this plan's Deferred section (the
~100 mas FWHM replacing the current 20 mas, and the γ value for the F146≈22 transition).

**Where this should land:** before Phase C's astrometric Fisher-matrix work (`inputB`/
`inverB`, params `tetE, mus1, mus2, piE`) and Step G1 (which explicitly reruns the
astrometric matrix to isolate satellite-parallax from temporal-baseline parallax) —
both currently consume `l->erra[]` values for Roman epochs that trace back to this same
placeholder. Detection itself (Step B1's `FFG[0]` gate) is unaffected, since it's driven
by `dchiL` (photometric), not `dchiA`.

## Controlled Fisher-matrix test harness (from Step C1)

Step C1 needed to verify a numerical claim (adding `t0` to the photometric Fisher parameter
set must not *decrease* `sigma(tE)` for any event). Two attempts to check this concretely both
failed:

- The live Monte Carlo's detection efficiency in the current test field (`s->lon` 0.5-0.6,
  `s->lat` -1.0 to -0.9) is far too low (about 1 detected-and-Fisher-scored event per 2600 star
  draws) to gather "a handful" of real events inside a reasonable wall-clock budget.
- A hand-built synthetic single event (temporary code, not committed) produced a numerically
  degenerate Fisher matrix — absurd sigma values across *every* parameter, not just `t0`.

**CORRECTION (Step C0).** The second failure was originally recorded here as "the test
configuration itself was ill-conditioned rather than the code being wrong." That diagnosis was
wrong. The real cause was the Fisher-accumulation defect fixed in Step C0: `FisherM` used
`gsl_matrix_set` (overwrite) instead of accumulating the per-epoch information sum, so the matrix
held only the last data point's contribution and was rank-1 by construction. The synthetic test
was reporting a genuine bug, and the same pathology was present in the live Monte Carlo (a real
detected event reported relative sigma of 1.8e7 on `tE`). Every Fisher number produced before
Step C0 is meaningless.

**What's still needed:** a small, deliberately well-conditioned synthetic lens/source/light-curve
fixture (or a short list of a few), checked in as an actual reusable test rather than throwaway
code, that a future Fisher-matrix change (Steps C2-C5 all touch `FisherM`) can run quickly to
confirm sigma values move in the expected direction — without depending on the live Monte Carlo's
detection efficiency. The "why did it degenerate" question is now answered, so such a fixture
should be straightforward to build and trust.

**Where this should land:** before or alongside Step C3 (the finite-difference step-size
convergence sweep), since both need the same kind of controlled, repeatable single-event setup.

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
