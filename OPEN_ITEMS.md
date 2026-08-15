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

## §5.1 visit-list model (from Step A2)

`whitepaper.tex` §5.1 still describes a single merged Roman+Rubin visit list with `FoV = 1.75°`.
The code now uses two separate visit lists and cadences: `FoV` (1.75°) for Rubin via
`BulgeBaseline.dat`, `FoVRoman` (0.28°, currently a placeholder — see the `TODO(Ali)` next to its
declaration in `Bulge.h`) for Roman via `RomanBaseline.dat`, matched independently by two calls to
`matchVisibleEpochs()`. Fix in Phase G, Step G3.
