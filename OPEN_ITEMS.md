# OPEN_ITEMS.md

Running list of known-but-deferred discrepancies surfaced while working through
`JOINT_FIT_REFACTOR_PLAN.md`. Items are removed once resolved; unresolved items are carried
into the whitepaper's open-items discussion at Phase G, Step G3.

---

## §5.1 visit-list model (from Step A2)

`whitepaper.tex` §5.1 still describes a single merged Roman+Rubin visit list with `FoV = 1.75°`.
The code now uses two separate visit lists and cadences: `FoV` (1.75°) for Rubin via
`BulgeBaseline.dat`, `FoVRoman` (0.28°, currently a placeholder — see the `TODO(Ali)` next to its
declaration in `Bulge.h`) for Roman via `RomanBaseline.dat`, matched independently by two calls to
`matchVisibleEpochs()`. Fix in Phase G, Step G3.
