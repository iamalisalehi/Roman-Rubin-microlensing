# OPEN_ITEMS.md

Running list of known-but-deferred discrepancies surfaced while working through
`JOINT_FIT_REFACTOR_PLAN.md`. Items are removed once resolved; unresolved items are carried
into the whitepaper's open-items discussion at Phase G, Step G3.

---

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

---

## Black-hole remnant masses use a single proportional slope

`remnantMass()` maps an initial mass above 20 Msun to a black hole by `Ml = 0.24 * Mi`,
giving ~4.8 Msun at Mi = 20 and ~28.8 at Mi = 120. That spans the observed stellar-mass
black hole range and is monotone, but it is the crudest link in the mass chain.

Real remnant masses depend on metallicity and on mass loss during the progenitor's life in
ways no single slope reproduces: at low metallicity weaker winds leave heavier black holes,
and the relation is not monotone across the pair-instability region. There is also no
attempt at a natal-kick velocity distribution, so black hole lenses share the kinematics of
their progenitors.

**Why it matters here.** Black holes are the long-tE tail, and the long-tE regime is exactly
where the joint fit earns its keep -- Roman's thetaE and Rubin's piE combining into a lens
mass. The *shape* of the black hole mass distribution therefore feeds directly into the
headline precision numbers, not just into a tail nobody looks at.

**Why not now.** Choosing an initial-final mass relation for black holes is a literature
decision with real spread between prescriptions, and it should be made deliberately rather
than folded into a bug fix. The white dwarf branch is on firmer ground (Kalirai et al. 2008,
calibrated on open clusters) and the neutron star branch uses the canonical 1.4 Msun, so
this is the one segment carrying an arbitrary choice.

---

## Five events have sigma_joint > sigma_Rubin, on matrices with condition number > 1e9

Adding data cannot worsen a Fisher forecast: the joint information matrix is the sum of the
per-survey ones, so `sigma_joint <= sigma_single` is an identity, not an approximation.
`romanlib.check_monotonicity` asserts it, and on the 2026-08-30 Kroupa run it reports five
violations out of 74,812 joint-detected events, in both `tE` and `piE`, with a maximum
excess of **1.0011** -- one part in a thousand.

All five have the same signature:

| tE (d) | piE | ndw_L | ndw_R | condA_J | condA_L | sigtE_J | sigtE_L |
|---|---|---|---|---|---|---|---|
| 10.16 | 0.111 | 2361 | 49972 | 8.55e10 | 8.53e10 | 284.05 | 283.75 |
| 10.39 | 0.180 | 1075 | 49634 | 1.64e11 | 1.64e11 | 54.32 | 54.32 |
| 12.21 | 0.172 | 893 | 46908 | 5.33e10 | 5.33e10 | 147.28 | 147.27 |
| 5.95 | 0.180 | 2360 | 49972 | 7.01e10 | 7.01e10 | 59.12 | 59.10 |
| 9.30 | 0.055 | 2371 | 49984 | 1.28e09 | 1.28e09 | 71.58 | 71.58 |

**This is round-off, not a partitioning bug.** A double carries about 16 significant digits;
inverting a matrix with condition number 1e11 loses about 11 of them, leaving ~5 -- so a
relative discrepancy of 1e-3 between two nearly identical inversions is the expected size,
not a surprise. Note also that every one of these forecasts is meaningless on its own terms:
`sigma_tE = 284 d` on a `tE = 10 d` event is not a measurement, and none of them pass the
characterization criterion.

**What is still worth doing.** The pipeline currently reports a sigma for any matrix GSL
manages to invert, with no conditioning floor. A cut -- refuse to report when
`condA > 1e8`, say, and set the sentinel instead -- would remove this class of artifact
outright and would also stop absurd sigmas propagating into medians. Choosing the threshold
needs a look at the condition-number distribution across the run, which is why it is
recorded here rather than applied.

Until then, `check_monotonicity` will keep printing this warning on every figure script. It
should stay noisy: the day it fires on a *well*-conditioned event, that IS a bug.

---

## The binary's git stamp goes stale whenever `make` has nothing to do

`Makefile:18` captures `GIT_COMMIT` at compile time and bakes it into the binary, which
writes it to `run_provenance.txt`. But the stamp only refreshes when something actually
recompiles. Editing a source file, building, and then committing leaves the binary carrying
the *pre-commit* label -- and a later `make` reports "Nothing to be done" and keeps it.

That is exactly what happened to the 2026-08-30 production run. `run_provenance.txt` and
every figure derived from it read `git_commit=d6dd293-dirty`, while the sources that were
actually compiled are the content of `5c74fbd` ("Lens the bulge with bulge stars, not
MACHOs"). **The science is unaffected -- the right code ran -- but the label on 2.5 GB of
output names the wrong commit.**

**Fix options, none yet chosen:** make the stamp a `.PHONY` prerequisite so every build
refreshes it; or have the binary refuse to run when `git diff --quiet HEAD` fails; or emit
the stamp at *run* time by shelling out, which costs a subprocess but cannot go stale.

**Interim rule for anyone running a production job:** `make clean && make` *after* the last
commit, never before, and check the `git_commit` line in `run_provenance.txt` against
`git rev-parse --short HEAD` before starting.

---

## F3's Roman-alone comparison rests on 1,950 events

The (`tE`, `piE`) characterization map's panel (a) -- "what Rubin adds to Roman" -- can only
use events inside Roman's footprint, because outside it Roman characterizes nothing for a
reason that is geometric, not physical. On the 2026-08-30 run that is **1,950 of 74,812**
joint-detected events, because only ~39 of 1,706 scanned sightlines fall inside the GBTDS
footprint at the current stride.

At 0.5 dex cells that leaves 16 coloured cells out of 56, most of the plane grey. The signal
in those cells is real (ratios up to 2.0) but the error bars are not small, and no cell
carries enough events to quote a trend within it.

**What would fix it:** sampling that concentrates draws inside the Roman footprint rather
than spreading them uniformly over the scan region -- which is Step E1's stratification
question wearing a different hat. A footprint-weighted run would buy roughly a factor of 20
more Roman-observed events for the same wall-clock, at the cost of needing explicit weights
to recover survey-wide totals.

Do not quote panel (a)'s cell values as precise until this is addressed. Panel (b), on all
74,812 events, does not have this problem.

---

## F2's long-`tE` bins rest on 124 events, and that is where the null result lives

The `sigma_piE` version of F2 (`DEVIATIONS.md` 24) reports that the joint fit adds nothing to
Roman's parallax precision for events longer than 300 days -- the ratio is flat at ~0.98
across the whole season gap. That is the negative answer to the plan's headline long-`tE`
prediction, so it deserves more support than it currently has.

It rests on **124 events**, spread over six `dt_edge` bins of 16-26 events each. The
supporting diagnostic (94% of those events have `piE` measured at better than 2 sigma by
Roman alone, median `sigma_piE` = 0.0098) is a strong and internally consistent explanation,
but the medians themselves are drawn from small samples, and the `> 300 d` bin is the one the
Kroupa mass function populates most sparsely -- long events need heavy lenses, which are
rare, and the remnant tail is exactly where the mass function is least certain.

**Why it matters scientifically:** "Rubin does not improve Roman's parallax measurement of
long events" is a claim about black-hole and neutron-star lens characterization, which is the
science case the long-`tE` regime exists for in this forecast. A null result at n = 124 is
suggestive; it is not yet quotable in a paper.

**Why it is deferred:** the fix is not a change to the analysis, it is more events in that
corner of parameter space, which is Step E1's stratified sampling -- the same fix the F3
footprint item needs. Doing it now would mean a second production run before the sampling
question is settled.

**What the fix involves:** stratify the draws in `tE` (Step E1) so the long-`tE` bins are
populated to comparable depth as the 30-100 d bin, and carry explicit weights to recover
survey-wide totals. Re-run F2 with `--param piE` on the stratified table and check whether
the flat ~0.98 line survives.

## F1/F2/F3 still load the whole 5.57M-row table and will be OOM-killed on a small machine (found during Step F4)

`analysis/romanlib.load_events()` read the entire per-event table into one DataFrame. For
`test5.dat` that is 5,571,168 rows x 90 float64 columns -- about 4 GB resident before
pandas' parse buffers are counted. On a machine with ~7 GB of RAM the kernel kills the
process.

**The failure mode is what makes this dangerous, not the failure itself.** The kill is
silent: the shell reports **exit status 0** with an empty stdout. It looks exactly like a
script that ran and chose to print nothing, and it produces no figure, no CSV and no error.
It was mistaken for a no-op for several minutes during Step F4 before the memory arithmetic
was checked.

**What was done:** `load_events()` gained optional `keep=` and `chunksize=` arguments that
filter per chunk while reading, holding the peak at one chunk plus the surviving rows. The
default path is byte-for-byte unchanged, so no existing caller's behaviour moved.
`analysis/f4_fisher_precision.py` uses it (`keep=lambda c: c["detJ"] == 1`), and the
production figure was verified identical to the reference computed the old way.

**What is deliberately NOT done:** `f1_results_table.py`, `f2_gap_filling.py` and
`f3_characterization_map.py` still call `load_events(path)` with no filter. They were run
successfully on 2026-08-30/31 when more memory happened to be free, so their outputs stand
-- but **they may not be reproducible on this machine as they are.** They are not being
changed here because each one selects a different subset (F2 gates on mission scope and the
Roman footprint, F3 on joint detection, F1 aggregates per field), so each needs its own
`keep` predicate written and its output re-verified against the existing CSV. That is a
step of its own, not a side-effect of F4.

**The fix, when it is taken:** give each of the three a `keep=` predicate that discards
undetected events during the read -- they are 98.66% of the table and every one of these
scripts throws them away immediately anyway -- then diff the regenerated CSV against the
committed one to prove nothing moved.

---

## Two Phase F panels pool events across sightlines and do not yet apply the area weight (raised during Step E1)

**What is wrong:** Step E1 makes the scan stratified — sightlines inside Roman's footprint can
be visited on a finer grid than those outside — so the sightlines no longer stand for equal
pieces of sky. Every event row now carries `w_area`, the deg² its sightline represents, and
`romanlib.area_weight()` returns it. **Two existing panels pool events across sightlines and
still count them one-for-one:**

- `analysis/f3_characterization_map.py` **panel (b)** — joint over Rubin-alone over all
  joint-detected events.
- `analysis/f4_fisher_precision.py` in its all-detections mode (`f4_fisher_all_kroupa.png`).

Everything else in Phase F is safe by construction and was checked one at a time: F1 tabulates
per Roman field (all in-footprint, all at the same fine step), F2 is restricted to the
footprint by construction (Deviation 23.1), F3 panel (a) and the F4 footprint panels are
in-footprint selections. Per-event ratios and per-sightline efficiencies are unaffected by
construction — which sightlines were visited is not an input to either.

**Why it matters scientifically:** unweighted, those two panels would describe a sky in which
Roman covers whatever fraction of the *sample* the stratification bought, instead of the ~2.6%
of the *sky* it actually covers. On the existing unstratified table the weight is a constant
and the panels are correct as they stand; the error appears only once a stratified run exists,
and it will not look like an error — F4's all-detections mass ratio would simply drift off
1.000, which is exactly the number Deviation 25.3 uses as a *bug detector* for the survey
partitioning. A weighting mistake would be read as a partitioning mistake.

**Why it is deferred:** applying the weight changes what those two panels plot — a weighted
CDF and a weighted 2D histogram rather than counts — and therefore changes committed figures
that currently stand as results. That is a science decision about how the all-sky panels should
be normalized, not a mechanical fix, and there is no stratified run to validate against yet.
Doing it now would mean editing figures to match a table that does not exist.

**The interim guard, already in place:** `romanlib.is_stratified()` reads the `stratified` flag
that the simulator now writes into `run_provenance.txt`, and `describe()` — which every figure
script already prints into its footer and its stdout — appends
`STRATIFIED(stride_roman=N; weight by w_area)` when it is set. A stratified run therefore
cannot be plotted by these scripts without saying so on the figure itself.

**What the fix involves:** pass `weights=R.area_weight(df, prov)` into the histogram and CDF
calls in those two code paths, re-run both against the unstratified table first (the weights
are constant there, so every number must come out unchanged — that is the regression), then
re-run against the stratified table.

---

## Step E1's `tE` stratification is not implemented; only the sky half of E1 is (raised during Step E1)

**What is wrong:** `JOINT_FIT_REFACTOR_PLAN.md` Step E1 asks for a fixed number of events **per
`tE` bin**, reweighted afterwards by the `tE` distribution to recover absolute yields, with bin
edges 1–5, 5–10, 10–20, 20–30, 30–60, 60–90, 100–200, 200–500, 500–1000 days. What was built is
the sky-position half: stratify the *scan* toward Roman's footprint, weight by `w_area`. The
`tE` axis is untouched, and `FunctE`'s 100 linear bins from 0 to 50 years remain what they were.

**Why it matters scientifically:** the plan wants the `tE > 200` d tail populated because that is
where the black-hole lens result lives, and a population-weighted draw produces few of them.

**Why it is deferred:** Deviation 26.1 argues the plan misdiagnoses which axis is starving those
bins. `tE` is not a drawn quantity — it falls out of the lens mass, the two distances and the
relative proper motion — so stratifying in it means acceptance sampling on a derived value,
carrying a second independent weight through every pooled statistic in Phase F, on top of the
sky weight just introduced. And the bins the plan wants filled are *footprint* bins: the
long-`tE` null rests on 124 events because it is a subset of the 1,950 in-footprint events, not
because long events are rare in the draw. The sky stratification multiplies that count by ~23
at `--stride-roman 2` with no acceptance sampling and no new weight to reason about.

**What the fix involves, if it is still wanted after a stratified run:** an acceptance test on
`l->tE` placed immediately after `func_lens()` and before the light-curve loop — which is where
it is nearly free, since the expensive per-draw work is the loop over ~50,000 epochs and
`FisherM`, all of it downstream of the point where `tE` is known. Accept with probability
`p[bin]`, carry weight `1/p[bin]`, and keep `nstE`/`ndtE` counting only accepted draws so the
per-bin efficiency stays unbiased (acceptance depends on the bin alone, so nothing *within* a
bin is distorted). The trap is the same one as above and worse: two weights now multiply, and
every absolute yield needs both.

---

## Roman's halo orbit around L2 is not modelled (raised 2026-09-04, planning Phase H)

**What is wrong:** Step H1 will place Roman at the mean L2 point. Roman actually flies a halo
orbit about L2 with an amplitude of order 10⁵–10⁶ km, so its true offset from Earth varies
through the mission.

**Why it matters scientifically:** it modulates an effect that is itself only ~10⁻³ in `u`, so
it is a fractional correction to a small term. It could matter for the narrow
high-magnification niche where satellite parallax does anything at all, because there the
sensitivity to `Δu` is highest.

**Why it is deferred:** the leading term has to exist before its correction is worth having,
and Step H3's experiment will say whether the effect is large enough for the correction to
change any conclusion. Modelling the halo orbit also needs an actual ephemeris or orbit
specification, which is not in the repository.

**What the fix involves:** replace the constant `L2_OFFSET_AU` displacement with a
time-dependent one, sourced from a published Roman orbit specification, and re-run H3's
comparison to see whether anything moves.

---

## No astrometric-shift analysis product exists (raised 2026-09-04, planning Phase H)

**What is wrong:** the astrometric microlensing signal — the centroid deflection
`δθ = θE · u / (u² + 2)` — is computed per epoch in `lightcurve()` (`s.def1c`/`s.def2c`, stored
in `l.soux`/`l.souy`) and feeds the astrometric Fisher matrix, but **no analysis script reads
it**. F4 plots `sigma_tetE`, which is the forecast *precision* on the angular Einstein radius,
not the shift itself, so nothing in the project says how large the astrometric signal is or how
often it exceeds the per-epoch astrometric precision.

**Why it matters scientifically:** whether the astrometric signal is *detectable* is a
different question from whether `θE` is *forecastable*, and the lens mass needs `θE` and `piE`
together. The maximum shift is `θE/√8` at `u = √2` — i.e. it peaks *outside* the photometric
peak, which interacts with Roman's seasonal gaps in a way nothing has yet looked at.

**Why it is deferred:** it is Step H5 of `PHASE_H_PLAN.md`. Its blocker, Step H4, is now
done — `errRomanA()` exists and Roman no longer borrows Rubin's astrometric error — so H5 is
unblocked and is simply not written yet.

**What the fix involves:** `analysis/h5_astrometric_shift.py` per `PHASE_H_PLAN.md` H5, then
re-run F4 on a post-H4 table and record how far the `tetE` and `Ml` panels moved. That
number is itself a result: it says how much the `errlsstA` placeholder was distorting them.

---

## RESOLVED 2026-09-05 — Rubin's astrometric error was a mission average used per epoch

**Status: closed by commit `d393dc3`.**

The suspicion recorded here was right and the cause was worse than "unsourced".
`files/sigmaA_LSST.txt` is a **mission-averaged** curve and `errlsstA()` fed it straight into
`l.erra[]`, which is a **per-epoch** error that `FisherM` divides by once per epoch and then
sums over ~2,300 epochs — applying the sqrt(N) averaging twice. Rubin's per-epoch astrometry
was **26.7x better than reality** and its astrometric Fisher information ~715x too large.

Two independent checks fix the factor. The table's bright-star floor is 0.3739576 mas, and
10.0/0.3739576 = 26.74 — exactly Ivezic et al. 2019's "assumed astrometric accuracy of 10 mas
per observation per coordinate" over sqrt(715), and ~715 visits is the right order for a
ten-year all-band LSST count. Scaled by the same 26.74 the faint end reads 132 mas at
r = 24.44 against an independent seeing-limited estimate FWHM/SNR ~ 700/5 ~ 140 mas. The
table's shape was right; only its normalisation was wrong, so the fix renormalises rather than
replaces, and does it in code so the delivered input file is untouched.

Roman is now the better astrometer, which is the physically obvious ordering and was not true
before. 10 mas is the conservative reading: Rubin commissioning (SITCOMTN-159) reports a 3-7
mas single-visit systematic, so the delivered floor may be better than assumed here.

**Every `tetE` and lens-mass number involving Rubin astrometry moves.** Nothing has been
re-run against it yet; the 2026-09-05 v2 production run is the first that will be.


## The astrometric shift is not diluted by blending (found in Step H5)

**What is wrong:** the simulator's modelled centroid, `s.pos1c`/`s.pos2c` in `lightcurve()`,
adds the lensing deflection `s.def1c`/`s.def2c` undiluted — it is the SOURCE's centroid. A real
measurement sees the centroid of every star in the aperture, so an unlensed blend of fraction
`(1 - fb)` drags the measured shift down by roughly `fb`. The per-epoch astrometric *error* is
evaluated at the blended magnitude, so the photon-noise half of blending is modelled; the
centroid-dilution half is not.

**Why it matters scientifically:** it makes every astrometric signal optimistic by a factor of
about the source-flux fraction, and the bulge is a crowded field where `fb` is often well below
1. `theta_E`, and hence every lens mass that uses it, inherits the optimism. **For Rubin this is
large** — the median `blend_r` is 0.15, so a real Rubin centroid would move ~7x less than the
simulator's. For Roman it is currently moot, but only because `blend_F146` is pinned at exactly
1.0 by the entry below; fix that first and this one grows teeth.

**Why it is deferred:** it is a change to the forward model, not to an analysis script, and it
would move every astrometric number in the project at the same time as Step H4 already has.
Doing both in one step would make it impossible to attribute the change to either.

**What the fix involves:** multiply the deflection term by the per-telescope source-flux
fraction where `pos1c`/`pos2c` are formed, and confirm against the fixture that `sigma_tetE`
degrades by roughly `1/fb`. Meanwhile `analysis/h5_astrometric_shift.py` panel (a) plots the
undiluted shift and the `x fb` version side by side, so the size of the simplification is
visible on the figure rather than hidden.

---

## RESOLVED 2026-09-05 — Roman was unblended because a Poisson count was drawn as a Gaussian

**Status: closed by commit `58d2863`.** The diagnosis recorded here was wrong; checking it is
what found the real cause.

**This entry guessed the input population was at fault** — that `Nstart` counts only stars
above some completeness limit, so the faint unresolved population that actually blends a Roman
pixel was missing from the density. **It is not.** `Nstart` is built from the mass density
divided by the mean mass of each population, so it is complete down the whole IMF and stops at
no survey's detection limit. The geometry was right too. The statistics were wrong.

The mean number of field stars in a seeing disc is a surface density times a disc area: ~12.7
for Rubin's 0.993" r-band disc at a bulge `Nstart` of 2.13e8 per deg^2, and **~0.14 for
Roman's 0.105" F146 disc**. The old code added `RandN(sqrt(mean), 2.0)` and clamped up to 1. At
12.7 that is a passable Gaussian approximation to a Poisson count. At 0.14 it is qualitatively
wrong: the truncated Gaussian can add at most 2*sqrt(0.14) = 0.75, and the clamp then rounds
**every** draw to exactly 1. A mean of 0.14 neighbours was rendered as "no neighbours, always".

The fix is a Poisson draw, `1 + Poisson(mean)` — the `1 +` because we are looking AT a source,
so the disc is conditioned to contain it, and Slivnyak's theorem says the rest of a Poisson
field is unchanged by that conditioning. `max(1, Poisson)` would absorb the first neighbour
into the source and leave Roman unblended 99% of the time, the same bug one step on.

**Measured, same stub patch, before and after:** `blend_F146` fraction below 1 goes from
**0.000 to 0.146**, against the predicted 1 - exp(-0.142) = 0.133; median stays 1.0, because
Roman really is mostly unblended at 0.105" — that part of the old answer was right. `fb1` is
off its boundary for the same 14.6% of events instead of pinned for all of them. Rubin moves as
predicted and only slightly, `nsbl_r` 14.3 -> 15.0, the same Palm-conditioning correction.

**Consequences to re-examine once the v2 run lands:** Roman's detection efficiency was
`testR <= blend[6]` with `blend[6] == 1`, so Roman accepted every detectable event while Rubin
accepted ~15%. That asymmetry fed every Roman-against-Rubin yield comparison, including the
gap-filling claim. It should now be smaller, and the size of the change is a result in itself.


## The provenance stamp goes stale silently, and discipline has now failed twice

**Status:** open. Found on 2026-09-05, while checking the freshly launched production run.

**What happens.** `Makefile` captures the git description into `GIT_COMMIT` at the moment
`make` runs, and compiles it into the binary so every run's provenance block can name the
source that produced it. But `make`'s dependency graph is over *files*, and a commit changes no
file. So after `git commit`, a subsequent `make` says "Nothing to be done for 'all'" and the
binary silently keeps the stamp it was built with — which names the previous commit, with a
`-dirty` suffix, because the tree was dirty when it was built.

The Makefile already knows this and says so in a comment: *"the value is captured when make
runs, so `make clean && make` is what refreshes it after a commit -- an incremental build keeps
the stamp its objects were built with."* The mitigation is therefore a rule a human has to
remember.

**It has now failed twice.** Commit `585432b` ("Rebuild so the next run is labelled with the
commit that produced it") exists solely because the stamp read `d6dd293-dirty` when the source
was `7a1b591`. On 2026-09-05 the identical thing happened again: the run launched with a binary
stamping `6c97375-dirty` when the source was `e8f4135`. A rule that has been broken every time
it has been tested is not a working control.

**Why it matters more than it looks.** The stamp is the only link between a multi-GB table and
the code that produced it. Every figure, every quoted number and every caveat in `PROGRESS.md`
section 4 is ultimately justified by "this run contains change X", and the stamp is the evidence.
A stamp naming a commit that predates the change it is being used to vouch for inverts that: a
future session would look up `6c97375`, find no H2 columns in it, and conclude either that the
table is corrupt or that the columns came from somewhere unrecorded.

**Why it was not fixed on the spot.** The run it was found on was already in flight and the
finding is a build-system change, not a physics change; folding it into a launched run would
break the rule that a step does one thing. The specific run was rescued instead, by proof
rather than by restart — see `PROGRESS.md` section 5c.

**What the fix involves.** Make the stamp a file so `make` can reason about it. Generate a
`git_commit.h` from a phony rule that rewrites it only when the description actually changes,
and have the sources include it; then a commit invalidates exactly the objects that embed it and
nothing else, and no one has to remember anything:

```make
.PHONY: force
git_commit.h: force
	@echo '#define GIT_COMMIT "$(GIT_COMMIT)"' > $@.tmp
	@cmp -s $@.tmp $@ || mv $@.tmp $@   # only touch it when it really changed
	@rm -f $@.tmp
```

Two smaller things belong with it. The binary should be able to print its own stamp
(`./roman --version`) so it can be checked without launching a run or reading a provenance
block. And `--dry-run` should print the stamp it would write, since `--dry-run` is already the
documented way to inspect a run's cost before committing hours to it, and the stamp is part of
what you would want to inspect.

---

## RESOLVED 2026-09-05 — DET_ANOMALY was never 0; the run total was never counted

**Status: closed by commit `a5bb600`.** Kept because the wrong hypothesis in it is instructive.

**What was suspected.** The 2026-09-05 partial run logged 433 DET_ANOMALY events in 23 of its
24 Roman-covered sightlines where the 2026-08-30 run reported 0 in 5.57M events, and this entry
named **H1 (satellite parallax)** as the leading suspect, on the reasoning that it is the only
one of E1a/H1/H2/H4 that touches the photometric model.

**That was wrong, and the experiment this entry specified is what showed it.** A paired A/B on
`--no-satellite-parallax` — same seed, same stub patch inside the footprint, `--events` equal to
`--maxdraws` and `--lenses` set out of reach so both variants execute exactly 150 draws per
sightline whatever they detect, hence identical events — came out **indistinguishable**:

| | none | Rubin+joint | Roman+joint | both+joint | DET_ANOMALY |
|---|---|---|---|---|---|
| satellite parallax on | 1197 | 36 | 81 | 35 | **15** |
| satellite parallax off | 1197 | 36 | 81 | 35 | **15** |

**The actual cause.** `nDetClass[DET_ANOMALY]`, the per-sightline counter, was incremented.
`NDetClassTot[DET_ANOMALY]`, the run total, was **not** — even though it is read at the end of
the run to decide whether to print the note explaining what DET_ANOMALY means. Every run since
the counter was introduced reported exactly zero anomalies in its summary and the explanatory
note never fired once. **The 2026-08-30 run's "0 ANOMALY" is that bug, not a measurement.**
The anomalies were in its log all along; nobody read that log, because it was unavailable.

**So nothing regressed.** The rate is what Deviation C-D.2 predicted: thresholding `dchi`
against `2 * ndw` is a bar on the MEAN per-epoch improvement, and where Roman contributes
50,401 epochs against Rubin's ~2,300, Rubin's epochs lift the joint bar by ~4,728 while adding
little signal. Any event clearing Roman's own bar by less than that fails the joint one — about
1% of footprint draws. `detJ` is still forced monotone, so no table was ever wrong.

**Still open, and now with evidence behind it:** the `2 * ndw` threshold form itself, which
Deviation C-D.2 deliberately left alone. A chi-squared statistic should be thresholded on its
total, not its mean. With Roman's epoch count 22x Rubin's, that asymmetry is no longer
hypothetical, and the anomaly rate is the measurement of how hard it now bites.

---

## `--dry-run` destroys the previous run's output, and it has already cost one table

**Status: open. Found 2026-09-05 the hard way — it destroyed the v1 partial production table.**

**What happens.** `Bulge_LSST.cpp` opens the event table with `std::ofstream head(testf);` —
truncating mode, the default — and writes the column header, at around line 632. The
`--dry-run` early exit is at around line 916. So **every `--dry-run` truncates `test5.dat` and
zeroes the append-mode files in `files/MONTLMC/files/` before it decides not to simulate
anything.** A flag whose entire purpose is "build the scan and tell me the cost without
running" is destructive to the previous run's output.

**What it cost.** At 18:07 on 2026-09-05, `./roman --dry-run --stub --stride 2` was run from
the worktree to check the freshly rebuilt provenance stamp. The worktree's `test5.dat` was at
that moment still symlinked to the v1 partial run directory, so the dry run truncated it from
773 MB and 1,643,024 rows to a 620-byte header, and zeroed its six MONTLMC outputs.

The scientific loss is small — `run.log` (97 MB) and `README.md` survived, so the cost model,
the detection-class totals and the 433 DET_ANOMALY events were all already extracted into
`PROGRESS.md` and `DEVIATIONS.md`, and v1 was superseded by v2 in any case. But `PROGRESS.md`
§5c said that table was preserved, and it was not.

**Why this is dangerous rather than annoying.** The per-event write is
`filg_in.open(testf, std::ios::app)` **inside the event loop**, so the path is re-resolved on
every event. That means a live run cannot be protected by repointing the symlink out of harm's
way mid-run — the running process would simply follow it. While a production run is in flight
through the worktree symlink, **any** `./roman` invocation from the worktree, `--dry-run` and
`--help`-adjacent flag checks included, destroys that run's table.

**What the fix involves.** Move the header write after the `--dry-run` exit, which is where it
belongs — a dry run should touch no output file at all. Two things belong with it:

1. **Refuse to truncate a non-empty table without being told to.** The run already knows the
   file name; if `test5.dat` exists and is larger than its header, require an explicit
   `--overwrite` (or write to `test5.dat.N`). Losing hours of simulation to a mistyped command
   should not be possible.
2. **Say what is about to be destroyed.** The startup banner already prints the configuration;
   printing "truncating an existing 773 MB table" would have made this visible immediately.

**Until then, the working rule:** never invoke `./roman` from the worktree while a production
run is writing through the worktree symlinks, not even `--dry-run`. Use a copy of the binary in
an isolated directory (the pattern in `setup_variant.sh`) for any check made while a run is in
flight.

---

## `dchiP` and `dchiA` are still absolute values while `dchiL` is signed

**Status: open. Created deliberately by Step H7 on 2026-09-06, not inherited.**

Step H7 made the lensing detection statistic signed:

```cpp
dchiL   = chi3  - chi1;             // signed: chi3 is the flat-baseline residual
dchiP   = std::fabs(chi2  - chi1);  // still an absolute value
dchiA   = std::fabs(chi2a - chi1a); // still an absolute value
```

**Why `dchiL` had to change.** Only `chi3 - chi1 > 0` — the lensing model fitting *better* than
a flat baseline — is evidence of lensing. Under `fabs`, a lensing model fitting far *worse*
than the baseline would have cleared the detection bar. More importantly it breaks the
monotonicity H7 exists to establish: an instrument whose epochs all fall outside the event
contributes a small negative difference from noise, and `|dchi_L + dchi_R|` can then fall below
`max(|dchi_L|, |dchi_R|)`, which is exactly the `DET_ANOMALY` condition. Signed, the sum is
exact and the joint test cannot contradict a single-survey detection.

**Why the other two were left alone.** `dchiP` (parallax) and `dchiA` (astrometric deflection)
are not detection tests — nothing thresholds them. They are reported as the *size* of a
perturbation, for which an absolute value is defensible. Changing them would alter two more
output columns for no acceptance criterion, and Step H7's remit is the detection decision.

**Why it is still an open item.** Three columns that look like the same kind of quantity now
have two different sign conventions, which is exactly the sort of thing that produces a wrong
plot two months from now. Either give `dchiP`/`dchiA` the same signed convention — noting that
"the parallax model fits worse than the no-parallax model" is a meaningful negative that
`fabs` currently hides — or rename the two absolute ones so the asymmetry is visible at the
call site. Do it as its own step with its own before/after, not as a drive-by.

---

## Step H7 made footprint sightlines ~6x more expensive, and the mechanism is unknown

**Status: open, measured 2026-09-06 on the v3 runs. Not blocking — the runs are correct and
proceeding — but unexplained, and the explanation may matter.**

**What was predicted.** H7 lowers the detection bar from `2*ndw` (100,802 inside the footprint)
to a fixed 500, so detections arrive about twice as fast. The per-sightline loop stops on its
`--events`/`--lenses` targets, so reaching those targets in fewer draws should make every
sightline *cheaper*. That is what was written in `PROGRESS.md` §5e and in Deviations 33.

**What happened.** Outside the footprint, exactly that: **14.3 s -> 6.8 s** per sightline.
Inside it, the opposite and larger: **~346 s** per sightline, against somewhere between 52 s
and 136 s for v2. The footprint term dominates the total, so the run got *slower*: projected
~8.8-13 h for v2 against **~16 h** for v3.

**The v3 number is solid and the v2 number is not.** Both v3 runs agree independently
(299 and 329 s/sightline at n=11, 339 and 346 at n=41), and v3's outside cost rests on 618
samples. But v2's per-stratum split came from a two-run solve that was already flagged as
weakly determined: its footprint cost is 52 s if its outside sightlines cost 14.3 s and 136 s
if they cost 6 s, and v2's outside cost was never measured cleanly — the same mistake this file
records elsewhere, made a second time. **So the regression is a factor of 2.5-6.6, direction
certain, magnitude not.** Quoting "six times" without that range overstates what was measured.

**Why it is strange.** The comparison is at *identical* sky positions — the v2 and v3 scans are
the same deterministic grid, and the first twelve footprint sightlines were verified to match
on `nri`, lon and lat. At those same positions v3 does **less** work:

| per footprint sightline | v2 (pre-H7) | v3 (post-H7) |
|---|---|---|
| scored events (lightcurves built) | 333-917 | 300-435 |
| detections (= Fisher calls) | 56-87 | 50-81 |

Fewer lightcurves built, no more Fisher matrices, the same 50,401 Roman epochs and the same
`ndw` — and six times the wall clock. Nothing in the H7 diff obviously accounts for that.

**The first hypothesis was tested and is REFUTED (2026-09-06).** It proposed that a bar of 500
admits marginal events whose Fisher matrices are near-singular, so that `invert_matrix` does
expensive work that is then thrown away. Characterisation rates on footprint events say
otherwise:

| detections characterised (`okA == 1`) | v3 post-H7 | v2 pre-H7 |
|---|---|---|
| joint | 2405/2419 = **99.4%** | 3276/3301 = 99.2% |
| Rubin | 1308/1320 = 99.1% | 1928/1954 = 98.7% |
| Roman | 1693/1709 = 99.1% | 2201/2225 = 98.9% |

Marginal detections invert *at least as reliably* as the overwhelming ones did. Nothing is being
computed and discarded. Whatever costs the time, it is not failed inversions.

**What the same measurement makes stranger.** Per footprint sightline, v3 builds **307**
lightcurves against v2's **510**, and makes **59** Fisher calls against v2's **51**. Fewer
lightcurves, ~16% more Fisher calls, and several times the wall clock. For the arithmetic to
work, the cost *per Fisher call* would have to have risen roughly fivefold — but `FisherM`
evaluates a fixed finite-difference stencil whose cost goes as `ndw * Nx`, and neither changed.
**A second hypothesis was tested and is also REFUTED (2026-09-06).** It proposed that the old
`2*ndw` bar systematically rejected high-epoch events — an event with more epochs faced a
proportionally higher bar — so that H7 now admits exactly the events whose Fisher matrices cost
the most, since `FisherM` scales with `ndw`. The table says no:

| mean `ndw` over joint detections, footprint | v3 post-H7 | v2 pre-H7 |
|---|---|---|
| all footprint events | 51,299 | 51,217 |
| joint detected | **51,086** | **50,914** |
| detected/all ratio | 0.996 | 0.994 |

A 0.3% difference. Inside the footprint `ndw` is dominated by Roman's 50,401 epochs and is
near-constant across events, so `2*ndw` was effectively a *constant* bar as well — just a very
high one. Its scaling with epoch count did its damage *between partitions* (the joint bar
sitting above Roman's own), which is what H7 fixed, not between individual events.

**What is established, and what is not.** The controlled A/B settles that the regression is
real and not an artefact of comparing two production runs: with a *fixed* draw budget, so both
sides build identical lightcurves, the pre-H7 binary completed 18 footprint sightlines to the
post-H7 binary's 3 in the same 900 s, while making only 1.8x as many Fisher calls. That is
roughly 3x the cost per Fisher call, at identical `ndw`, with `FisherM` itself unchanged.

**No further hypotheses should be invented here.** Two have been proposed and both refuted by
measurement; a third would be guessing. Settle it by profiling — `perf record` on one footprint
sightline for each binary at the same `--start-index` and the same fixed draw budget, which is
about 20 minutes and will name the function. Until that is done, the honest statement is that
H7 costs roughly 3x per characterised event for reasons not yet identified.

**How to settle it cheaply.** Run one footprint sightline under `perf record` (or simply time
`FisherM` and `invert_matrix` with a counter) against the pre-H7 binary at the same
`--start-index`, with a fixed draw budget so both sides see identical events — the pattern in
`h7_ab.sh`. About 20 minutes, and it competes with the production runs for cores, so it belongs
between chunks rather than during one.

**Why it might matter beyond runtime.** If marginal detections are producing ill-conditioned
Fisher matrices, the question is not only what they cost but whether their `sigma` values are
trustworthy. `okA`/`okB` gate on successful inversion, so nothing invalid should reach the
table — but the *fraction* of detections that fail characterisation is now a number worth
reporting, because H7 changed which events are offered to the Fisher step at all. Check
`sigmaas`/`okA` rates on the v3 table against v2 before quoting any yield.

