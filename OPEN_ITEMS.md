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
