# PHASE H — Satellite parallax, astrometric shift, and the collaborator draft

**Written:** 2026-09-04, at the end of the Step E1a session.
**Status:** H1 and H4 **done** (2026-09-04, Deviations 28 and 29); H2, H3, H5, H6 planned.
**Supersedes:** `JOINT_FIT_REFACTOR_PLAN.md` Step G1, whose premise turns out to be false —
see §0.2 and `DEVIATIONS.md` entry 27.

This document is a *plan*, in the same sense and the same working style as
`JOINT_FIT_REFACTOR_PLAN.md`: **its §0 working agreement governs every step here** — one step
at a time, teaching brief before any edit, explain the physics and not just the code, define
every acronym, flag uncertainty rather than inventing an explanation, one `git commit` per step
named for the physical change.

Read `PROGRESS.md` first, then this.

---

## 0. Why this phase exists

The user asked for three things, in this order:

1. **Satellite parallax.** Rubin is on the ground; Roman is at the Sun–Earth L2 point. The two
   observatories therefore see a microlensing event from different places, and the difference
   is information about the lens. Calculate it, put it in the analysis, plot it.
2. **The astrometric shift.** Same treatment: calculate, analyse, plot.
3. **The whitepaper**, brought fully up to date with everything the code now does, good enough
   to hand to potential collaborators. (To be revised once more when the project is finished;
   this pass is for the collaborator conversation, not the journal.)

### 0.1 What the code actually does today — check this before believing anything else

**There is no satellite parallax in the simulation, and there never has been. Both telescopes
are treated as if they sat at the centre of the Earth.**

The observer's projected displacement is built in `lightcurve()` (`Bulge_LSST.cpp`, the
"Light curves and Astrometry" block) as `as.ue_n1` / `as.ue_n2`, from Earth's orbit alone:

```cpp
void lightcurve(source & s, lens & l, astromet & as, double timh)   // <-- no telescope argument
{
    ...
    dvex = +vearth * std::sin(omegae * tt + M_PI / 2.0) / omegae;
    dvey = -vearth * std::cos(omegae * tt + M_PI / 2.0) / omegae;
    as.Ve_n1 =  std::cos(tetp) * dvex * std::sin(l.deltao) - dvey * std::cos(l.deltao);
    ...
}
```

The function has **no telescope parameter**, and all four of its call sites pass the same
geometry regardless of which observatory produced the epoch — even though `l.tele[i]` (0 =
Rubin, 1 = Roman) is in scope at three of them:

| Call site | What it is |
|---|---|
| the light-curve fill loop | generating the data |
| `FisherM`, first-derivative loop | ∂(model)/∂θ |
| `FisherM`, second-derivative loop | the mixed terms |
| `FisherM`, the chi-squared/reference loop | the baseline |

So every `piE` forecast the project has produced — including F4's — contains **only the annual
(Earth-orbit) parallax**, sampled at whatever times the two surveys happened to observe. It
contains no contribution from the ~0.01 AU baseline between Earth and L2. Those forecasts are
therefore a *lower bound* on what the real pair of observatories can do. This is recorded in
`OPEN_ITEMS.md`.

### 0.2 Why this kills Step G1 as written

`JOINT_FIT_REFACTOR_PLAN.md` Step G1 proposes: *"run the joint fit twice — once with real
geometry, and once with Rubin's observer position forced to Roman's ... The difference isolates
the satellite-parallax contribution from the temporal-baseline contribution."*

**That experiment cannot be run, because there is no spatial baseline to switch off.** Forcing
Rubin's observer position to Roman's would change nothing: they are already the same position.
Run today, G1 would return exactly zero and the zero would be read as a physics result rather
than as a missing term.

Phase H therefore does G1's job in the opposite order: **first put the satellite geometry in
(H1), then the two-run experiment becomes meaningful (H3)** — and it is the same experiment,
with the "off" run being today's code rather than a special mode.

### 0.3 The physics, stated honestly up front

The plan's own framing (G1) is right and should survive into the whitepaper:

- Roman at L2 is ~1.5 × 10⁶ km from Earth ≈ **0.01 AU**.
- The difference in the impact parameter seen by the two observatories is
  `Δu ≈ (D_perp / AU) · piE`, where `D_perp` is the observer separation projected
  perpendicular to the line of sight and `piE` is the microlensing parallax amplitude.
- For a typical bulge event with `piE ≈ 0.1–0.3`, `Δu ≈ 10⁻³`. That is a **small** effect.
- It matters where the light curve is *steep* in `u` — high magnification (small `u0`), short
  `tE`, small projected Einstein radius. Free-floating-planet-like lenses are the niche.

**So the expected result is a small gain concentrated in a narrow corner of parameter space,
and a near-null over the bulk of the sample. That is a publishable result, not a failure** —
it is the same posture the plan takes about the short-`tE` gain, and it is the honest framing
that keeps the headline claim (temporal gap-filling, Deviations 23–25) credible. Do not merge a
satellite-parallax number and a temporal-baseline number into one "parallax improvement"
figure; they are different physical effects and the plan is explicit that conflating them
invites exactly the criticism the honest framing avoids.

### 0.4 Scheduling — read this before launching any production run

Step E1a (footprint-stratified sampling) is built and waiting for a `--stride-roman` choice and
a machine. Step H1 changes the *physics* of every light curve that has a Roman epoch.

**Recommendation: do H1 first, then launch ONE production run that has both.** A stratified run
launched before H1 would have to be thrown away and repeated once satellite parallax exists,
and these runs cost 15 h and up. The alternative — run E1a now for the sampling-limited results
(F3 panel (a), the long-`tE` `piE` null) and re-run after H1 — is defensible if those results
are needed for a specific deadline, but it costs a second full run. **This is a scheduling
decision for the user, not a technical one.**

---

## PHASE H STEPS

### Step H1. Give the observer a position

**Goal:** `lightcurve()` learns which telescope is asking, and returns the trajectory that
telescope actually sees.

**Teach me, in the brief, before any edit:**
- What frame `ue_n1` / `ue_n2` live in. They are built from `dvex`/`dvey` (Earth's orbital
  position components) rotated by `l.deltao`, `s.FI` and `tetp` into a basis tied to the
  event's position on the sky. Say what that basis is, in words, and which way `n1` and `n2`
  point. Note that `vearth` is defined as `omegae` in `Bulge.h`, so `dvex = sin(...)` — these
  are dimensionless position components in units of AU, not velocities, despite the names.
- The derivation of the L2 offset in that same basis: L2 is anti-sunward of Earth on the
  Sun–Earth line, so to leading order Roman's position is Earth's scaled by `(1 + 0.01)`
  radially. Write down the projected separation `D_perp` and confirm the `Δu ≈ 0.01 · piE`
  scaling falls out.
- Which of `u0`, `t0` and `piE` are *definitions relative to an observer* and which are not.

**THE TRAP, and it is the whole step.** `lightcurve()` currently runs its `ig` loop twice and
**subtracts the value at `t = 0`**:

```cpp
if (ig == 0) { int1  = as.Ve_n1; int2  = as.Ve_n2; }
if (ig == 1) { int1 -= as.Ve_n1; int2 -= as.Ve_n2; }
```

That is a gauge choice: the parallax offset is defined to vanish at `t = 0`, so that `u0` and
`t0` mean what they mean. There is already a comment at the FisherM reference loop noting that
"an epoch at `t == 0` has zero parallax offset BY CONSTRUCTION".

**If each observer subtracts its own `t = 0` position, the constant offset between the two
observers is erased — and that constant offset IS the satellite parallax.** The step would then
compile, run, produce plausible-looking output, and measure nothing.

#### How to overcome it — the derivation, and the two lines it comes down to

Write `P(·)` for the map that takes a heliocentric position and returns its projection into the
`(n1, n2)` sky basis of this event. In the code that is

```cpp
Ve_n1 =  cos(tetp)*dvex*sin(deltao) - dvey*cos(deltao);
Ve_x  = -cos(tetp)*dvex*cos(deltao) - dvey*sin(deltao);
Ve_n2 = -sin(FI)*Ve_x + cos(FI)*sin(tetp)*dvex;
```

**`P` is linear in `(dvex, dvey)`** — every term is a constant times one component. That single
fact is what makes this easy, and it should be stated in the brief because the whole solution
rests on it.

What the code computes today is `ue(t) = P(X_E(t)) - P(X_E(0))`, where `X_E(t)` is Earth's
heliocentric position: a displacement measured **from Earth's position at `t = 0`**. That
choice of origin is the gauge, and it is what fixes the meaning of `u0` and `t0`.

To leading order Roman sits on the Sun–Earth line, beyond Earth, at
`X_R(t) = (1 + f) · X_E(t)` with `f = L2_OFFSET_AU ≈ 0.01`. The correct trajectory for each
telescope is its own position measured **from the same single origin** — keep Earth at `t = 0`,
so that every existing `u0` and `t0` keeps its meaning and Rubin is untouched:

```
ue_Rubin(t) = P(X_E(t))       - P(X_E(0))                     <- exactly what the code does now
ue_Roman(t) = P((1+f)·X_E(t)) - P(X_E(0))
            = [P(X_E(t)) - P(X_E(0))]  +  f · P(X_E(t))        <- by linearity of P
            = ue_Rubin(t)              +  f · P(X_E(t))
```

So the fix is: **keep the existing differenced term exactly as it is, and add `f` times the
UN-differenced projection evaluated at `t`.** The `ig` loop and its `t = 0` subtraction do not
change at all; the satellite term is simply never passed through them.

In code, inside the existing `ig` loop, capture the undifferenced value at `ig == 0`:

```cpp
if (ig == 0) { int1 = as.Ve_n1; int2 = as.Ve_n2;  abs1 = as.Ve_n1; abs2 = as.Ve_n2; }
if (ig == 1) { int1 -= as.Ve_n1; int2 -= as.Ve_n2; }
...
const double f = (tele == 1) ? L2_OFFSET_AU : 0.0;   // 0 for Rubin => today's expression
as.ue_n1 = int1 + f * abs1;
as.ue_n2 = int2 + f * abs2;
```

**Why this is right and the naive version is wrong, in one sentence each.** The naive version
computes `P(X_R(t)) - P(X_R(0))`, which by linearity is `(1+f)·[P(X_E(t)) - P(X_E(0))]` — a
0.01 *rescaling* of Earth's annual parallax ellipse, with the constant inter-observer offset
gone. This version computes `P(X_R(t)) - P(X_E(0))`, which keeps it.

**The sanity check that distinguishes them, and it must be in the verification.** At `t = 0`
the two observers must NOT agree: `ue_Roman(0) - ue_Rubin(0) = f · P(X_E(0)) ≠ 0`. Print it for
one event. If it comes out zero, the gauge has eaten the signal. The magnitude should be
`|Δu| ≈ f · piE ≈ 0.01 · piE`, i.e. ~10⁻³ for a typical bulge event — check that too, because
it is the number §0.3 says the physics requires and it falls out of the code independently.

**The change:**
- `Bulge.h`: `constexpr double L2_OFFSET_AU = 1.5e6 / 1.496e8;` — with the numbers and their
  source in the comment, not a bare `0.01`.
- `lightcurve(source&, lens&, astromet&, double timh, int tele)`; `tele == 1` (Roman) adds the
  L2 displacement, `tele == 0` (Rubin) is exactly today's expression.
- Thread `tele` through **all four** call sites. The three inside `FisherM` are the ones that
  matter most: **if a derivative is evaluated with a different observer than the datum it is
  differentiating, the Fisher matrix is inconsistent** and the error it produces is not a
  forecast of anything. Use `int(l.tele[i])` at each, exactly as `FisherM` already does for
  `s.fb[tt]` / `s.mbs[tt]`.

**Acceptance criteria** (both required):
1. With `L2_OFFSET_AU` forced to 0, a `--stub` run is **byte-identical** to the previous
   commit's, on every column of `test5.dat`, `MapLMC5.dat`, `LpLMC5.dat`, `EfLMC5*.dat`. This
   is the same regression pattern Step E1a used; the recipe is in `PROGRESS.md`.
2. With it at its real value, **only rows with Roman epochs move**. An event with
   `ndw_R == 0` must be bit-identical. Check this explicitly — it is the cheapest possible test
   that `tele` was threaded correctly and not, say, applied to every epoch.
3. `./fishertest` still passes. Extend the fixture with one event whose epochs are split
   across both telescopes and whose `piE` is large, so the satellite term has something to bite
   on, and record its sigmas as the new regression baseline.

**What could break, checked against the code:**
- `s.ut0` and the "without parallax" deflections `def1a`/`def2a` are formed by subtracting
  `l.piE * as.ue_n*`, so folding the satellite term into `ue_n*` means "without parallax" now
  also means "without the satellite offset". That is the consistent reading — it is one
  observer-displacement effect, not two — but say so in the brief rather than letting it happen
  silently, because `dchiP` (the parallax contribution to chi-squared) is derived from it and
  its meaning widens accordingly.
- **The detection thresholds do not move.** `detL`, `detR` and `detJ` are thresholded on
  `dchiL` (the *lensing* effect), not on `dchiP`; `dchiP` is recorded in the output table and
  never gates anything. Verified by reading the detection block. So H1 cannot change which
  events are detected except through the light curve itself, which is the intended effect.
- The astrometric expressions `s.pos1b`/`pos2b` and `l.pos1`/`pos2` carry `- ue_n* · pis` and
  `- ue_n* · pil` (source and lens parallax). They pick up the satellite term automatically
  once `ue_n*` carries it, which is correct and is what Step H5 needs — **do not special-case
  them out.**

---

### Step H2. Record the satellite-parallax observable per event

**Goal:** make the effect *visible* in the output table, not only implicit in the Fisher matrix.

Add to `EventRecord` and the per-event table:
- `du_sat` — the observer separation in Einstein radii, `(D_perp / AU) · piE`, evaluated at
  `t0`. This is the amplitude of the effect for that event, and it is the natural x-axis of
  every plot in H3.
- `nep_both` — whether the event has epochs from both telescopes *while it is magnified*
  (say, within ±2 `tE` of `t0`). Satellite parallax needs contemporaneous coverage; an event
  Roman saw in season 3 and Rubin saw in season 7 has none, however many epochs each has.

Both are cheap and both are things the analysis will otherwise have to reconstruct badly from
other columns. Append them to `eventTableHeader()` and the row writer together, and verify with
`head -1 | wc -w` against a data line, as Step D1 established.

---

### Step H3. The satellite-parallax experiment and its figures

**This is Step G1 done properly**, now that there is something to switch off.

**The experiment:** two production runs, identical seed and configuration, differing only in
`L2_OFFSET_AU` (real vs 0). Match events by row; every event appears in both.

**Why this and not an analytic estimate:** the question is not "how big is `Δu`" — that is
`0.01 · piE` and needs no simulation. The question is how much the *forecast precision on
`piE`* improves once the two observatories no longer see the identical trajectory, which
depends on the cadence, the photometric errors, and the correlation structure of the Fisher
matrix. Only the pipeline answers that.

**Figures (separate, separately captioned — do not merge):**
- **H3a — where the effect lives.** `sigma_piE(with L2) / sigma_piE(without)` for the joint
  fit, against `du_sat`, coloured by `tE`. Expect ~1 everywhere except a corner.
- **H3b — the corner.** The same ratio in the (`tE`, `u0`) plane, restricted to events with
  contemporaneous coverage (`nep_both`). This is where the high-magnification, short-`tE`
  niche should appear, if it appears at all.
- **H3c — the honest comparison.** Satellite-parallax gain and temporal-baseline gain side by
  side, as two distributions on one axis, with the caption saying in words that they are
  different physical effects. The temporal number already exists (Deviations 23–25); this is
  what puts the satellite number next to it at the right scale.

**Report the null loudly if it is a null.** State the fraction of events with any improvement
better than 1%, and the properties of the ones that have it. "Roman–Rubin satellite parallax
adds nothing measurable outside N events with `u0 < x` and `tE < y`" is a citable sentence and
disarms an obvious referee question.

**Everything goes through `analysis/romanlib.py`** — sentinel rules, `okA`/`okB` gating, and
now `area_weight()` if the run is stratified.

---

### Step H4. A real astrometric error model for Roman  *(blocks H5)*

**This is a dependency, not a nicety.** `Bulge_LSST.cpp` currently does:

```cpp
errsR = errlsstA(*ls, magni[fiR]);   // Roman astrometric error -- PLACEHOLDER
```

Roman's per-epoch astrometric uncertainty is Rubin's model evaluated at Roman's magnitude.
Every `tetE` number in the project — including F4's "90.1% of events measure `tetE` to better
than 10%" — rests on it, which `PROGRESS.md` §4 and `OPEN_ITEMS.md` already flag. **Any
astrometric-shift figure produced before this is fixed is a figure about `errlsstA`.**

#### The model, researched 2026-09-04

Roman's per-exposure astrometric precision is anchored by two independent sources:

- **Sanderson et al. 2019**, *Astrometry with the Wide-Field Infrared Survey Telescope*
  (arXiv:1712.05420), §1.1: *"single-exposure precision for well-exposed point sources is
  0.01 pixel, or about 1.1 mas"*, and a factor ~10 improvement to 0.1 mas by stacking ~100
  exposures.
- **Black hole astrometric binaries in the Roman GBTDS** (arXiv:2608.24998), Figure 5 and
  surrounding text, which is GBTDS-specific and gives the magnitude dependence: *"Centroiding
  of 1% of a pixel is a commonly assumed astrometric precision, corresponding to a floor of
  1.1 mas for Roman"*; the floor *"impacts bright sources F146_Vega < 20.62 mag"*; sources
  become background-dominated around *"F146_Vega < 23.5 mag, which corresponds to
  σ_ast ≈ 10 mas"*; pixels are *"0.11 arcsec"*; each GBTDS exposure is *"66 seconds"* at a
  *"12.1 minute"* cadence, ~100 exposures per day. Their underlying curve comes from Pandeia
  and the Roman astrometry simulation tool of Bellini et al. 2024.

**PER EXPOSURE, and this is a factor of ten waiting to be got wrong.** The 0.1 mas figure that
appears in both sources is the *daily-binned* precision — ~100 exposures stacked. Our
`l.erra[i]` is a per-epoch error and **one row of `RomanBaseline.dat` is one 12.1-minute
exposure**: measured directly, the median inter-epoch gap for field (l, b) = (0.4, −1.2) is
0.008403 d = 12.1 min, with 50,401 epochs per field × 6 fields = 302,406 = `NlRoman`. So the
per-exposure number, 1.1 mas, is the one to use. Using 0.1 mas would overstate Roman's
astrometry tenfold and would flatter every `tetE` and lens-mass forecast in the project.

The model to implement, `errRomanA(mag)` returning mas:

```
                | 1.1                                    m <= 20.62      centroiding floor
sigma_ast(m) =  | 1.1 * 10^(0.3329 * (m - 20.62))        20.62 < m <= 23.5
                | 10.0 * 10^(0.4   * (m - 23.5))         m > 23.5        background dominated
```

- The **floor**, 1.1 mas, is 1% of the 110 mas pixel. It is a centroiding systematic, not
  photon noise, so it does not improve for brighter stars.
- The **middle slope**, 0.3329 per magnitude, is not a free choice and not a physical constant:
  it is `log10(10 / 1.1) / (23.5 - 20.62)`, the slope that connects the two anchor points the
  GBTDS paper quotes. Pure source-dominated photon noise (`SNR ∝ sqrt(counts)`) would give
  0.2/mag and pure background domination (`SNR ∝ counts`) gives 0.4/mag; 0.333 sits between
  them because the transition is already under way across this range. **Say this in the code
  comment.** It is an interpolation between two published points, not a derived SNR curve, and
  a later step could replace it with a Pandeia-derived table exactly as `errlsstA` reads
  `sigmaA_LSST.txt`.
- Beyond 23.5 the slope is the background-dominated 0.4/mag, continuous with the branch below.

Sanity values: 1.1 mas at m ≤ 20.62; 2.0 mas at 22; 10 mas at 23.5; 25 mas at 24.5.

**Shape and placement.** `errlsstA` interpolates a 96-row table (`files/sigmaA_LSST.txt`);
`errRomanA` needs no table and no state, so it is a closed-form function of magnitude and lives
beside `errlsstA` in `helper.cpp` (note that its photometric sibling `errRomanM` lives in
`Bulge_LSST.cpp` instead, because that one needs the `roman` struct). Constants and their
citations go in `Bulge.h`, not inline.

**Wire it in at exactly one site:** `Bulge_LSST.cpp`, where the Roman branch fills
`l.erra[ndw]`. The Rubin branch and the per-event scalar `s.errA` (Rubin r-band) must not
change.

**Verification:**
- An event with `ndw_R == 0` must be bit-identical to the previous commit.
- Roman epochs must move, and in the right direction: `errlsstA` at a typical bulge F146
  magnitude against `errRomanA` at the same magnitude — state both numbers and the ratio, so
  the size of the placeholder's error is on the record.
- Re-run F4 and report how far the `tetE` and `Ml` panels moved. **That number is a result**:
  it says how much the placeholder was distorting the two panels `PROGRESS.md` §4 already
  flags as resting on it.

---

### Step H5. The astrometric shift as a first-class product

**What already exists.** The astrometric microlensing signal is computed — `s.def1c`/`s.def2c`
in `lightcurve()`, the centroid deflection `δθ = θE · u / (u² + 2)` — it is stored per epoch in
`l.soux`/`l.souy`, and it already feeds the astrometric Fisher matrix (`inputB`, `Ny = 4`,
parameters `0 tetE, 1 mus1, 2 mus2, 3 piE`) which is computed per survey partition like the
photometric one. **Nothing needs deriving; what is missing is that no analysis product reads
it directly.** F4 plots `sigma_tetE`, which is the *forecast* on the Einstein radius, not the
shift itself.

**Teach me, in the brief:** why the maximum centroid shift is `θE / √8` at `u = √2`, why that
makes the astrometric signal peak *outside* the photometric peak, and why that matters for a
survey with seasonal gaps — an event whose photometric peak falls in a Roman gap may still have
its astrometric maximum inside a season, or the reverse.

**New product, `analysis/h5_astrometric_shift.py`:**
- distribution of the maximum centroid shift over the sample, against `tE` and `θE`;
- the fraction of events whose maximum shift exceeds Roman's per-epoch astrometric precision
  (this is the number that says whether astrometric microlensing is *detectable*, as distinct
  from whether `θE` is *forecastable*);
- per-survey `sigma_tetE` and `sigma_mus`, and the joint gain, in the F4 CDF style;
- the mass–distance solution: `Ml = θE / (κ · piE)` needs *both* matrices, so show the joint
  constraint in the (`piE`, `θE`) plane for a few representative events — this is the figure
  that shows why the project needs both telescopes and both matrices at once.

**And the satellite term reaches astrometry too.** `s.pos1b`/`pos2b` carry `- ue_n1 * pis`, the
*source* parallax, which becomes observer-dependent the moment H1 lands. H1 must therefore be
threaded through the astrometric expressions as well, not only the photometric ones — do not
let H1 stop at `s.ux`/`s.uy`.

---

### Step H6. The whitepaper, for collaborators

**Audience:** potential collaborators, not referees. It has to be *correct* and *complete about
what is not yet done*, and it does not have to be finished science. Flag open items explicitly
with "pending" / "known issue" language rather than omitting them — that is the plan's stated
style requirement and it is exactly right for this audience.

Carry over `JOINT_FIT_REFACTOR_PLAN.md` Step G3's list, which is still accurate:

- **§5.1** — currently describes a single merged visit list with `FoV = 1.75°`. Rewrite for
  separate Roman/Rubin cadence handling, `FoVRoman`, and the duplicate-timestamp deduplication
  in the overlap between adjacent Roman fields.
- **§5.3–5.4** — per-survey detection labels (`DetClass`) and the stopping criteria.
- **§6** — the photometric parameter vector as the code actually orders it
  (`0 u0, 1 tE, 2 fb0, 3 piE, 4 xi, 5 t0, 6 mbs0, 7 fb1, 8 mbs1` — Deviation 2), per-band
  fluxes, and the three-matrix-per-event structure.
- **Appendix** — the step-size convergence plot (`c3_step_sweep.png`).

And add what has happened since G3 was written:

- **Roman's real GBTDS schedule** (Deviation 16) — 10 seasons, `MISSION_START_DAY`, and the
  corrected `FoVRoman`. The old text describes none of this.
- **The lens population** (Deviation 22) — Kroupa (2001) plus remnants. Any number in the
  current draft predating this describes a MACHO population and must go.
- **The results as they stand** — F1 (+245 characterized events from combining), F2 (the
  gap-filling ordering, which is the *reverse* of what was predicted — Deviation 23.2), F2 in
  `piE` (Deviation 24), F3, F4 (Deviation 25). The gap-filling result is the headline; say
  plainly that it is a **yield** effect and not a precision effect, because that is what the
  data say and it is the more defensible claim.
- **Sampling** (Deviation 26) — the footprint stratification and the `w_area` weights, if the
  stratified run has happened by then.
- **Satellite parallax and astrometry** — H3 and H5, with H3's honest framing.
- **An open-items section** built from `OPEN_ITEMS.md`, including the `errlsstA` placeholder if
  H4 has not landed, and the fact that satellite parallax was absent from every result
  predating H1.

**Do the figures last**, and only from tables whose provenance stamp is correct — note that
every figure made from `test5.dat` carries a wrong `git_commit` in its footer (`PROGRESS.md`
§3), which must not be reproduced in a document going to collaborators.

---

## Ordering, and what depends on what

```
H1 (observer position)  ──┬─→ H2 (record it) ─→ H3 (experiment + figures) ─┐
                          └─→ (astrometric half of H1, needed by H5)       │
H4 (errRomanA)  ─────────────────────────────→ H5 (astrometric product) ───┼─→ H6 (whitepaper)
                                                                            │
E1a run (stratified sampling, already built) ──────────────────────────────┘
```

- **H1 before everything.** H2 and H3 are meaningless without it, and H5's astrometric parallax
  term depends on it.
- **H4 before H5.** An astrometric figure built on `errlsstA` is a figure about Rubin's error
  model wearing Roman's name.
- **H4's numbers are now sourced** (arXiv:1712.05420 and arXiv:2608.24998, see H4) and it no
  longer blocks on anything outside the repository.
- **H6 last**, and it wants the E1a production run to have happened, or it will quote
  sample-limited numbers to collaborators.
- **H1 before the E1a production run**, if only one run is affordable (§0.4).

## What is deliberately NOT in this phase

- **Roman's halo orbit around L2.** H1 uses the mean L2 offset. The halo orbit adds a periodic
  ~10⁵–10⁶ km wobble, i.e. a fractional correction to an effect that is already ~10⁻³. Record
  it in `OPEN_ITEMS.md` when H1 lands; revisit only if H3 finds the effect matters.
- **Step G2** (validation against Abrams et al. 2025). Unchanged, still wanted, still after
  this — and still requiring the reweighting caveat in the plan, or the comparison will look
  like a bug.
- **Step E1b** (`tE` stratification) and **Step E2** (stopping criteria). See `OPEN_ITEMS.md`.
