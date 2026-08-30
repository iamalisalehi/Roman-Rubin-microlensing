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

## Astrometric finite-difference steps (`Delta2[]`) are unswept, and two use a biased stencil (from Step C3)

Step C3 swept and retuned the **photometric** steps (`Delta1[]`) and fixed the stencil for the
photometric `tE` and `piE` (see `DEVIATIONS.md` entries 10 and 11). The astrometric matrix was
left untouched, and has both of the same problems:

- **`Delta2[]` has never been swept.** Its steps came from the same legacy codebase as the
  photometric ones, which turned out to be ~4 orders of magnitude too large. There is no reason
  to assume the astrometric ones are better placed, and every reason to expect they are not.
- **`tetE` and `piE` still use `sig2`**, the first-order biased stencil (`Bulge_LSST.cpp`, the
  `Delta2[j] * sig2[h]` lines). `mus1`/`mus2` correctly use the central `sig`.

**Why it was not fixed at the same time:** the effect could not be verified. The astrometric
branch also depends on `errlsstA()` standing in for a real Roman F146 astrometric error model
(first item in this file), so its sigmas are not yet trustworthy in absolute terms regardless of
the stencil. Fixing one of two unverifiable inputs in isolation buys nothing checkable.

**What's needed:** extend `runSweep()` in `tests/fisher_fixture.cpp` to sweep `Delta2[]` (the
harness already carries `deltaScale` and the CSV/plot pipeline; `covarian` would need the
astrometric equivalent), then retune and switch the two `sig2` uses to `sig`. Expect the same
qualitative outcome: absolute astrometric sigmas shrink toward their true, larger values, while
the paired joint/single ratios hold up.

**Where this should land:** before Step G1, which reruns the astrometric matrix to separate
satellite parallax from temporal-baseline parallax, and before any absolute `tetE` precision
number goes in the whitepaper.

## `fb` derivative stencil depends on which blend-fraction bin the event lands in (from Step C3)

The blend-fraction step `co.bb[]` (`Bulge_LSST.cpp`, inside `FisherM`'s data loop) is chosen by
binning on `s.fb[tt]`: the middle bin gives `{-0.07, +0.07}`, a proper central difference, but the
outer bins give `{+0.07, +0.15}` and `{-0.07, -0.15}` — two forward (or two backward) differences,
carrying the same first-order bias as `sig2` in `DEVIATIONS.md` entry 11. So `fb`'s derivative
accuracy depends on where `fb` happens to sit, which is not a property anyone chose deliberately.

**Why it is not urgent:** the C3 sweep showed `fb0`/`fb1` flat to <0.3% even at the *unscaled*
steps, and `kFDStepScale` now shrinks them by 1e-4, making the residual bias negligible. The bin
structure and the bound-safety clamp are also now far from binding, since a step of ~7e-6 cannot
push `fb` out of [0,1].

**What's needed:** replace the binned step with a symmetric `{-h, +h}` pair at a fixed small `h`,
keeping a clamp only as a guard. This would also let the bin table and the clamp block be deleted,
simplifying `FisherM`'s data loop. Do it when `FisherM` is next opened for other reasons.

---

## Sky coverage is Penny et al.'s, not the current GBTDS footprint — and Rubin's FoV overlap is unmodelled (raised by Ali, 2026-08-24)

**What is wrong:** the region scanned by the Monte Carlo (`l1`/`l2`/`b1`/`b2`/`wid` in `Bulge.h`)
and the Roman field centres (`FIELDS_L_B` in `Baseline/generateRomanBaseline.py`) both descend
from the Penny et al. survey design. That design is now superseded — STScI's published GBTDS
footprint has changed since. The authoritative source is the Roman documentation already recorded
in `DEVIATIONS.md` entry 16:
<https://roman-docs.stsci.edu/roman-community-defined-surveys/galactic-bulge-time-domain-survey>

Note this is a *different* correction from entry 16.4. That fixed `FoVRoman`'s units (an area used
as a radius); the field **centres and the total footprint** were left as they were.

**The subtlety that makes this more than a coordinate update:** Rubin's field of view is very much
larger than Roman's (`FoV` vs `FoVRoman` in `Bulge.h`). A Rubin pointing whose centre lies well
outside a GBTDS field can still cover part of that field. The current `matchVisibleEpochs` treats
both instruments identically — a sightline matches an epoch when the sightline is within the
instrument's radius of the pointing centre — so this partial overlap is either counted as full
coverage or missed entirely, depending only on centre separation. Two consequences:

1. **Joint coverage is mis-stated.** Sightlines inside a GBTDS field that Rubin observes only via
   the edge of a large pointing are exactly the ones the joint-fit science case depends on.
2. **Blending is wrong there too.** Blend fractions (`s.blend[]`, `s.fb[]`) depend on how much
   Rubin actually sees at that position, so a coverage error propagates into the photometry and
   from there into `fb0`/`mbs0` and the Fisher forecast — not just into the event counts.

**Also required:** the Rubin baseline (`Baseline/readbaselineBulge.py` → `BulgeBaseline.dat`, and
`Nl` in `Bulge.h`) must be regenerated for whatever region the corrected footprint implies. The
two baselines share the simulation clock and the scan region, so they cannot be updated
independently.

**Why it is not being done now:** it changes what sky the forecast describes, so every detection
and precision number would have to be re-taken afterwards. It should be done as its own step, with
the new footprint read off the STScI documentation rather than inferred, and it interacts with the
run-scaling grid (Step 4) — the stride guard checks the grid against `FIELDS_L_B`, so the guard
protects the change but the field list itself is what needs replacing.

**Deliberately deferred at Ali's request, 2026-08-24.** Wanted, not optional.

---

## The per-event table is still called `test2.dat` and still lives in the repo root

**What.** The 88-column analysis table that Step D1 built is written to `./test2.dat`, a name left
over from when it was a debug dump. It sits in the repo root next to the binary, is truncated by
every run, and is `.gitignore`d, so a result can be overwritten by the next smoke test with
nothing to say it happened.

**Where it should go.** `files/MONTLMC/files/`, alongside the other outputs and `run_provenance.txt`
— ideally named for the run, so a production run and a stub run cannot clobber each other.

**Why not now.** `tests/c3_live_compare.py` reads it by path and by position, and is being actively
edited. Renaming mid-flight would break a comparison in progress for no scientific gain. It is a
five-minute change whenever that script settles.

**Note for whoever does it:** that script's `read_csv(..., header=None)` will need
`comment="#"` added, because D1 gave the file a header line. Its column-count assertion fails
loudly first, which is the intended behaviour.

---

## `TODO(Ali)` at `Bulge_LSST.cpp` (RomanBaseline.dat read) is stale

It asks for a generator sourced from Roman's own season structure rather than an LSST OpSim.
Commit `74e5e18` delivered exactly that (`Baseline/generateRomanBaseline.py`, against the ROTAC
2025 / STScI GBTDS design). The comment should be deleted; left in place only because removing it
is cosmetic and belongs with whatever step next touches that read.

---

## `co->flagi` is stale on uncharacterized events

`flagi` is set inside `FisherM`, which only runs for detected events, and it is **not** in the
per-event reset block that clears `okA`/`okB`/`condA`/`condB`/`Era`/`Erb`/`relMl`. So the `flagi`
column of the per-event table reads `1` on rows where nothing was characterized — it is the
previous characterized event's value.

**Harmless today**, because the only consumer (the per-field precision average) pairs it with
`co->okA[SJOINT]`, which *is* reset. But it is a trap for anything reading the table: **use
`okA_J`, not `flagi`, as the characterizability flag.**

Not fixed in Step D1 because `flagi` is also read inside `FisherM` itself and adding it to the
reset changes behaviour rather than only bookkeeping — it needs its own look at what `flagi` is
actually supposed to mean.

---

## Two output files open in append mode, so a re-run silently doubles them

**What.** `LpLMC2.dat` and `MapLMC2.dat` are opened with `std::ios::app`, not truncated:

```
Bulge_LSST.cpp:496    std::ofstream fil3(fnGam, std::ios::app);      // MapLMC2.dat
Bulge_LSST.cpp:1149   std::ofstream fil0_append(fnLDt, std::ios::app); // LpLMC2.dat
```

Every other output (`EfLMC2.dat`, `EfLMC2B.dat`, `magC0.dat`, `datC0.dat`, `test2.dat`) is
truncated at startup. These two are not, and nothing in the run says so: a second run appends its
rows to the first run's and the file ends up holding two runs' worth of events with no marker
between them. A stub run followed by a production run produces a `MapLMC2.dat` that is 36
sightlines of diagnostics glued to the front of the real dataset.

**How it was found.** Launching the Step D1 production run, 2026-08-29. The files had to be
cleared by hand first, and there is nothing in the code, the provenance block or the output that
would have caught it if they hadn't been.

**What it should be.** Either truncate them like every other output, or — if the append is
deliberate, e.g. accumulating across `IMnum` values the way `BHLSSTMONTS.dat` is cleared only when
`IMnum == 1` — say so in a comment and record the pre-existing row count in `run_provenance.txt`
so a downstream reader can tell where this run's rows begin.

**Why not now.** It needs a decision about whether the append was ever intentional, which means
looking at what `IMnum` is for. Changing it blind risks breaking a workflow that relies on the
accumulation.

---

## The startup file check makes an unread file's existence a precondition, and names no file when it fails

**What.** `LpLMC2.dat` is opened twice. Line 490 opens it for *reading*; nothing ever reads from
that stream. Its only purpose is to be tested at line 523:

```
Bulge_LSST.cpp:490    std::ifstream fil0(fnLDt);
Bulge_LSST.cpp:523    if (!fil0 || !fil2 || !fil2b || !fil3 || !fil4 || !fil5) {
                          std::cerr << "Cannot open one or more files!" << std::endl;
```

The actual writing is done by a separate append stream at line 1149. So the file must **exist** for
the run to start, even though its contents are never used — and if it does not, the run dies after
reading the baselines, the extinction maps and the whole CMD set, several minutes in, with a
message that names none of the six files it tested.

**How it was found.** The first attempt at the Step D1 production run, 2026-08-29, died here after
~4 minutes because `LpLMC2.dat` had been deleted to clear the append (see the item above). The log
showed the schedule guard and `read_cmd` succeeding, then `Cannot open one or more files!` with no
indication of which.

**What it should be.** Drop `fil0` if the existence of that file is genuinely not a precondition,
and check each stream separately with the filename in the message — the same discipline the Roman
schedule guard already follows, which prints the margins that failed rather than a bare verdict.

**Why not now.** Cosmetic on its own, and it belongs with whatever step next touches that block —
most naturally the `test2.dat` renaming item above, which moves output paths anyway.

---

## `numd[0]` counts drawn stars, not observable ones, so `EffiD` is 100% by construction

**What.** `icon` is incremented only inside the visibility gate:

```
Bulge_LSST.cpp:1149   if (flagf > 0 and ndw > 2) { //if star is visible
Bulge_LSST.cpp:1151       icon +=1;
```

but `records.push_back(...)` sits **outside** that block, so every drawn star gets a record, and
`numd[0]` — which counts records — is the number of stars *drawn*, not the number observable. The
two are equal only where every draw is observable.

**Consequences.**

- `EffiD = numd[0] * 100 / nsim` is described as "probability of detecting stars" but is 100% by
  construction, since a record is pushed on every draw.
- The per-sightline `[0]` means (`tE[0]`, `u0[0]`, `Ml[0]`, `fb[0]`, `mbs[0]`, `Map[0]`, `Ext[0]`,
  …) are averages over **drawn** stars, not observed ones.
- `EffiL = numd[1] * 100 / numd[0]` is detected-per-drawn, not detected-per-observable.

**How it was found.** `CHECK(icon == numd[0])` aborted the 2026-08-29 full-region run at
l = -3.499, b = -1.98, where 5 of 15 drawn stars were observable. It had never fired before because
every run prior to `81a6b04` used the dense 0.1x0.1 deg stub patch, where every draw *is*
observable and the equality holds by accident of the field.

**What was done.** The assertion was loosened to `CHECK(icon <= numd[0])`, which is the invariant
that actually holds. **No computed value changed** — only the assertion. The denominators are
untouched.

**What is still open.** Which denominator each quantity *should* use. If `EffiD` is meant to be the
fraction of drawn stars that are observable, it should be `icon / nsim`. If the `[0]` means are
meant to describe the observable population, the sums need the visibility gate. Both change
published numbers, so neither belongs in a bookkeeping fix — it needs a decision about what each
quantity is for, and a re-take of anything already quoted from them.
