# ORIENTATION.md

Produced under `JOINT_FIT_REFACTOR_PLAN.md`, Step A1. Read-only pass over `Bulge.h`,
`Bulge_LSST.cpp`, `Lensing.cpp`, `helper.cpp` (all four files read in full — no other files).
Nothing in the source was changed to produce this document.

Line numbers refer to the file state as of this read (git commit `366bc5b`).

---

## 1. Struct inventory

For anything not self-evident I've traced every read/write site across the four files. Where I
could not find a use, or where the meaning is genuinely ambiguous, I've said so rather than guessed
— see the per-struct "unclear/dead" notes and the consolidated list in §6.

### `source` (`Bulge.h:213-258`) — one simulated background star

| Member | Scope | Meaning / units | Set in | Read in |
|---|---|---|---|---|
| `nums` | scalar | Index (0..`Num`-1) into the `Num=9500`-bin radial distance grid; `Ds = nums*step`. | `func_source` (`Lensing.cpp:139`) | `func_lens`, `optical_depth`, main's per-distance-bin histograms |
| `cl` | scalar | Besançon CMD "CL" column, copied verbatim from `CMD::cl_*`. Exact code meaning (evolutionary/luminosity class?) is not defined anywhere in these 4 files — only ever assigned, never read back or printed live (one `cout` reference is commented out, `Lensing.cpp:136`). **Unclear**, see §6. | `func_source` | nowhere active |
| `typ` | scalar | Besançon CMD "Type" column, same situation as `cl` — assigned, never read. **Unclear**, see §6. | `func_source` | nowhere active |
| `mass` | scalar | Source stellar mass [M☉], from CMD, for the "primary" star only (`k==1` in the blending loop). | `func_source` | not read elsewhere in these files |
| `logT` | scalar | log10(effective temperature) [log10 K], from CMD. | `func_source` | not read elsewhere |
| `age` | scalar | Stellar age [Gyr] (CHECK bounds in `read_cmd` are ≤10/13/14 Gyr per population, consistent with Gyr). | `func_source` | not read elsewhere |
| `ros` | scalar | Source radius in Einstein-radius units: `Rsun*xls/RE` (dimensionless finite-source parameter). | `func_lens` (`Lensing.cpp:287`) | **never read** — only appears inside a commented-out debug print (`Bulge_LSST.cpp:1313-1314`). Finite-source effects are computed but not used; the light curve (`lightcurve()`) uses the point-source magnification formula throughout. |
| `mus1, mus2` | scalar | Source proper-motion components [mas/day], sky-plane basis. | `vrel` (`Lensing.cpp:540-541`) | `lightcurve` (trajectory), `FisherM` (astrometric Fisher parameters `j=1,2`) |
| `mus` | scalar | `sqrt(mus1²+mus2²)` — source proper-motion magnitude [mas/day]. | `vrel` | normalizes `vsave` for `EventRecord`; `resu[14]` |
| `xv, yv, zv` | scalar | **Declared, never assigned or read anywhere in these 4 files.** Dead. | — | — |
| `Av` | scalar | Extinction in V-band [mag] at the source's distance along the current sightline, from `interpExtinctionAlongSightline`. | `func_source` | used locally in `func_source` to derive `Ai[i]` per filter |
| `SV_n1/2, LV_n1/2, VSun_n1/2` | scalar | Source/lens/Sun velocity components projected into the sky-plane rotated basis (see `deltal/deltas/deltao` below). [km/s] | `vrel` | `vrel` (combined into `Vt`, `mus1/2`, `mul1/2`) |
| `pos1b/2b`, `pos1c/2c` | scalar | Source sky-plane trajectory [mas]: `b` = straight motion + parallax (no lensing deflection), `c` = same + lensing deflection. | `lightcurve` | `FisherM` (astrometric derivatives), diagnostic file output |
| `def1a/2a`, `def1c/2c` | scalar | Lensing deflection [mas]: `a` = without-parallax model, `c` = with-parallax (the one actually added to `pos*c`). | `lightcurve` | diagnostic output, `def1p/def2p` (previous-step proper-motion-from-deflection calc, main loop) |
| `errM, errA` | scalar | Running/averaged photometric ("M") and astrometric ("A") error accumulators for the *LSST-only* detection pass: `errM` = mean fractional-flux photometric error, `errA` = mean astrometric error [mas]. | main loop (accum over LSST epochs, `Bulge_LSST.cpp:468-469`), then averaged `582-583` | `EventRecord` (`l->DeltaT/s->errA` — a detection-significance-like ratio, stored under the field name `DeltaT`, see §6) |
| `FWHM` | scalar | **Not** a filter FWHM — full-width-half-max *duration* of this event's light curve [days], from `func_lens` (`Lensing.cpp:301`). Same name as `Bulge.h`'s per-filter `FWHM[M]` constant array; different thing entirely. | `func_lens` | detection gate `s->FWHM < Tobs` (`Bulge_LSST.cpp:588`), `EventRecord` |
| `ut, ut0` | scalar | Instantaneous lens-source separation in Einstein-radius units: `ut` = full model (with parallax), `ut0` = without the `piE` term. Feeds magnification via `Astar = (ut²+2)/sqrt(ut²(ut²+4))`. | `lightcurve` | main loop magnification calc |
| `Astar` | scalar | Magnification at the current timestep for the full (with-parallax) model. | main loop (computed inline, not inside `lightcurve`) | `FisherM`, per-band magnitude calc |
| `xi` | scalar | Trajectory angle [radians] — direction of source's motion relative to the lens, computed once in `vrel` from the sign quadrant of `vls1,vls2`. One of the 5 photometric Fisher parameters. | `vrel` | `lightcurve`, `FisherM` |
| `ux, uy` | scalar | Source position relative to lens, in Einstein-radius units, including parallax. | `lightcurve` | feeds `ut`, `def1a/2a` |
| `Ds` | scalar | Source distance [kpc]. | `func_source` | pervasive |
| `TET, FI` | scalar | Internal Galactic-angle representation of the current sky position: `TET=(360-lon)/RAa`, `FI=lat/RAa` [radians]. | main loop, from `s->lon/lat` | `Disk_model`, `vrel`, `optical_depth` indirectly |
| `lat, lon` | scalar | Galactic latitude/longitude [deg] of the current field pointing — doubles as the outer-loop induction variables in `main()`. | main loop | pervasive |
| `vs` | scalar | Source total space velocity magnitude [km/s]. | `vrel` | `EventRecord` |
| `Nstart, Rostart` | scalar | Sightline-integrated total stellar number density [stars/deg²] and total mass density [M☉/deg²] out to `MaxD=12` kpc. | `Disk_model` | source-count normalization (`s->nstart` calc), `Gamma`/`Neven` |
| `Romaxs, Romins` | scalar | Running max/min of the per-shell mass-density-weighted quantity `Rostari[]`, used to bound the rejection-sampling envelope for drawing a source's distance. | `Disk_model` | `func_source` (rejection sampling `rho > s.Rostari[nums]`) |
| `nstart` | scalar | Population-corrected estimate of the number of *detectable* stars per deg² for this field (folds the per-distance-bin detection efficiency `nsdet/nssim` against the physical density profile). | main loop (`Bulge_LSST.cpp:722-726`) | `Neven` (expected yield) |
| `nstarti` | scalar | **Declared, never assigned or read.** Dead — likely a discarded intermediate of `nstart`'s calc (name collision with the `Nstari[]` vector below, see §6). | — | — |
| `od_thin/thick/bulge/halo, opt` | scalar | Per-component and total microlensing optical depth along the sightline (dimensionless, later ×1e6 for display). | `optical_depth` | `Gamma` (event-rate normalization) |
| `fb[2]` | array, **by telescope** | Blend fraction: `fb[0]`=LSST r, `fb[1]`=Roman F146. See §4. | `func_source` (`Lensing.cpp:169-170`) | `FisherM` (perturbed live during Fisher differencing!), `ErrorCal` |
| `mbs[2]` | array, **by telescope** | Baseline (unlensed) blended magnitude: `mbs[0]`=LSST r, `mbs[1]`=Roman F146. See §4. Note: `lens::mbs` is a *different*, unrelated array (see that struct's table) — same field name, two different meanings. | `func_source` | `FisherM` |
| `nssim[Num], nsdet[Num]` | vector, per distance-bin | Number of stars simulated / number passing the pre-selection "detectable" test, histogrammed by `nums`. Reset per sky position. | main loop | `s->nstart` calc |
| `rho_thin/thick/halo/bulge[Num]` | vector, per distance-bin | Mass density [M☉/pc³] of each Galactic component at that distance along the sightline. | `Disk_model` | `func_source`/`func_lens` component draws, `optical_depth` |
| `rho_stars[Num]` | vector | **Declared, allocated, never assigned or read.** Dead. | — | — |
| `Rostar0[Num]` | vector | Sum of the four component densities [M☉/pc³]. | `Disk_model` | `func_source`/`func_lens` rejection sampling |
| `Rostari[Num]` | vector | `Rostar0[i]` converted to [M☉/deg²] for this shell (×`x²·step`×solid-angle factor). | `Disk_model` | source-distance rejection sampling bound, `s->nstart` |
| `Nstari[Num]` | vector | Number-density analog of `Rostari[Num]` [stars/deg²] (mass→number via fixed component-average-mass divisors). | `Disk_model` | `s->Nstart` accumulation |
| `nsbl[M]` | vector, **by filter** | Number of stars blended into one seeing/PSF disk, per filter — computed from that filter's own `FWHM[i]` (`Bulge.h:99`), so Roman's F146 already gets a much smaller crowding radius (0.105″) than LSST's ugrizy (~0.6-1.2″). Relevant to Plan Step B3 — on this reading it looks like per-instrument PSF is *already* handled correctly, not a single shared radius, but that deserves the plan's own dedicated investigation rather than taking my word for it. | `func_source` | `func_source` (self, blending loop), diagnostics |
| `blend[M]` | vector, **by filter** | Fraction of aperture flux from the source star itself (vs. blended neighbours), per filter. | `func_source` | `s.fb[]` (telescope-indexed copy), main loop pre-selection/detection |
| `Fluxb[M]` | vector, **by filter** | Summed linear flux (Σ10^-0.4·Map over all blended stars) in the aperture, per filter. | `func_source` | `magb[]`, `blend[]` derivation |
| `magb[M]` | vector, **by filter** | Blended/composite baseline magnitude, per filter. | `func_source` | main loop peak-magnitude test, `s.mbs[]` |
| `Ai[M]` | vector, **by filter** | Extinction, per filter [mag]. | `func_source` | `Map[]` derivation |
| `Mab[M]` | vector, **by filter** | Absolute magnitude of the primary source star, per filter (from CMD). | `func_source` | `Map[]` derivation, consistency `CHECK` |
| `Map[M]` | vector, **by filter** | Apparent magnitude of the primary star (distance + extinction applied, *not* diluted by blending), per filter. | `func_source` | `EventRecord` (`s->Map[2]`, r-band only) |
| `struc` | scalar | Which Galactic component (`GalacticComponent` enum) the source belongs to. | `func_source` | `vrel` |

### `lens` (`Bulge.h:260-322`) — one simulated lens object for a given source

| Member | Scope | Meaning / units | Set in | Read in |
|---|---|---|---|---|
| `numl` | scalar | Distance-grid index for the lens (analogous to `source::nums`). | `func_lens` | `func_lens` (self, component draw) |
| `Ml, Dl` | scalar | Lens mass [M☉], lens distance [kpc]. Mass is drawn from one of 4 selectable power-law IMFs gated by the compile-time `IMnum` constant (`Bulge.h:53`, currently `2`). | `func_lens` | pervasive |
| `vl, Vt` | scalar | Lens space velocity magnitude [km/s]; `Vt` = lens-source *relative* transverse velocity [km/s], the quantity `tE` is actually derived from. | `vrel` | `EventRecord`, `tE` calc |
| `xls` | scalar | `Dl/Ds`, dimensionless lens-to-source distance ratio. | `func_lens` | `RE`/`tetE` calc, `vrel` |
| `u0` | scalar | Impact parameter [Einstein radii], one of the 5 photometric Fisher parameters. | `func_lens` (drawn `RandR(0.001,u0m)`) | `lightcurve`, `FisherM` |
| `A0` | scalar | Peak magnification (dimensionless), from `u0` via the standard point-lens formula. | `func_lens` | `mi1/mi2` calc |
| `mi1, mi2` | scalar | Apparent magnitude at `A0+1` / `A0-1` (rough peak-brightness bracket). Not used in the main detection/Fisher pipeline — only written to the `IMnum==1`-gated diagnostic file `BHLSSTMONTS.dat`. | `func_lens` | diagnostic output only |
| `rhomaxl` | scalar | Scratch: max of the rejection-sampling envelope for the lens-distance draw. | `func_lens` | `func_lens` (self) |
| `tE` | scalar | Einstein-radius crossing time [days] — `RE / (Vt in km/day)`. One of the 5 photometric Fisher parameters. | `func_lens` | pervasive; `FunctE` bins it |
| `RE` | scalar | Einstein radius [m internally, `/AU` at output/`EventRecord` time]. | `func_lens` | `tE`, `tetE` calc |
| `t0` | scalar | Time of closest approach / peak magnification [days from survey start], drawn uniformly over `[2, Tobs-2]` — i.e. the **full 10-year window including Roman's gaps**, confirming Plan Step D2's "investigate" item: `t0` sampling already does cover gaps as far as these 4 files show; the plan should still verify this against however `Tobs` interacts with the field's actual epoch coverage. **Not currently a free Fisher parameter** — this is exactly Plan Step C1's target. | `func_lens` | `lightcurve` |
| `murel` | scalar | Relative lens-source proper motion magnitude [mas/day]. | `vrel` | `tE`-adjacent diagnostics, `ErrorCal` |
| `DeltaT` | scalar | `sqrt(4+u0²)·tetE` [mas] — an astrometric-deflection amplitude proxy. | `func_lens` | `EventRecord` (stored as `DeltaT/s->errA`, a ratio — see §6 naming note) |
| `piE, pirel, tetE` | scalar | Microlensing parallax (dimensionless), relative parallax [mas], angular Einstein radius [mas]. `piE = pirel/tetE`. | `func_lens` | `lightcurve`, `FisherM` (both matrices) |
| `pos1, pos2` | scalar | Lens sky-plane trajectory [mas] (parallel to `source::pos1c/2c` but for the lens). Only used in the `IMnum==1` diagnostic dump. | `lightcurve` | diagnostic output only |
| `mul1, mul2, mul` | scalar | Lens proper-motion components/magnitude [mas/day]. | `vrel` | `lightcurve`, `ErrorCal` (`resu[11-13]`) |
| `betal, betas, deltal, deltas, deltao` | scalar | Spherical-geometry angles [radians] used inside `vrel`/`lightcurve` to rotate Galactocentric velocity components (R,T,Z) into the sky-plane basis, separately for the lens, source, and Sun. I'm confident of their *role* (basis-rotation angles feeding `SV_n1/2`, `LV_n1/2`, `VSun_n1/2`, and the annual-parallax basis vectors `Ve_n1/2`); I'm not confident I can state their exact geometric derivation (which triangle/convention) purely from the code — that likely needs the whitepaper's coordinate-system section, which is one of `JOINT_FIT_REFACTOR_PLAN.md`'s explicitly deferred items ("Tile geometry details"). | `vrel` | `vrel` (self), `lightcurve` |
| `Nhalo[2], Nself[2]` | array | **Never assigned anywhere in these 4 files.** Read (divided, printed) at `Bulge_LSST.cpp:707-708,718-719`. A commented-out legacy block (`Bulge_LSST.cpp:241-247`) shows these used to be populated by reading an external file (`EfLMC2B.dat`) that is no longer read. As `std::array<double,2>` members with no initializer in `lens`'s constructor init-list, they are **default-initialized to indeterminate values** — this is reading uninitialized memory, not just "always zero." Flagged prominently in §6. | *nowhere* | `main()` output (`fnEff`, `fnEffB` files) |
| `nstE[GG+1], ndtE[GG+1]` | vector, per tE-bin | Intended: simulated / detected-lensing event counts binned by true `tE`, reset per sky position. **Never incremented anywhere in these 4 files** — same "used to come from a file read that's now commented out" situation as `Nhalo/Nself`. See §6; this appears to make `EFF` (and hence `Gamma`, `Neven`) always evaluate to exactly 0, which would fail `CHECK(EFF > 0.0)` at `Bulge_LSST.cpp:872`. | *nowhere* | `NstE/NdtE` accumulation, `EFF` calc |
| `NstE[GG+1], NdtE[GG+1]` | vector | Field-cumulative counterparts of `nstE/ndtE` (never reset, so accumulate across sky positions) — inherits the same "always zero" issue. | main loop | `fnEff`/`fnEffB` output only |
| `NsMl/NdMl, Nspi/Ndpi, Nsu0/Ndu0, Nsmb/Ndmb, Nsfb/Ndfb, Nsmu/Ndmu` | vector, per bin | Same intended pattern as `nstE/ndtE` but for `Ml, piE, u0, mb, fb, murel` respectively. Same dead-histogram situation — none of these are incremented anywhere either. | *nowhere* | `fnEff`/`fnEffB` output only |
| `timn[coun], magn[coun], soux[coun], souy[coun], errm[coun], erra[coun]` | vector, per accepted epoch | The per-event light-curve buffer: time [days], magnitude [mag, native band of whichever telescope produced it — see §4], sky position x/y [mas], photometric/astrometric 1σ error. This is the buffer `Nx`/coun-sizing (Plan Step A2) concerns. | main loop (LSST and Roman branches) | `FisherM` |
| `tele[coun]` | vector, per accepted epoch | Telescope tag: `0`=LSST/Rubin, `1`=Roman/F146. The key architectural asset noted in the plan (§0.4). | main loop | `FisherM` (`tt = tele[i]` indexes `s.fb[tt]`/`s.mbs[tt]`) |
| `tEs, Mls, pis, u0s, mbs, fbs, mus` (all `GG+1`) | vector | The fixed bin-edge/grid-value arrays (via `make_grid`) for the histogram family above — constructed once from `Bulge.h`'s `*_min/*_max` constants. **Note:** `lens::mbs` (this grid array, size `GG+1`, magnitude-bin edges 15–26 mag) is a completely different variable from `source::mbs` (size-2, telescope-indexed baseline magnitude) — same short name, unrelated meaning, in two different structs. | constructed at `lens()` construction | `FuncMb` (bins `Map[2]` into this grid) |
| `struc` | scalar | Which Galactic component the lens belongs to. | `func_lens` | `vrel` |

### `astromet` (`Bulge.h:324-327`)

| Member | Meaning / units |
|---|---|
| `Ve_n1, Ve_n2` | Earth's orbital-motion basis vectors [dimensionless, radian-scale] at the current time, projected into the sky-plane basis — the geometric core of the *annual* parallax signal. Computed each `lightcurve()` call by differencing Earth's position at `timh` vs. `t=0`. |
| `ue_n1, ue_n2` | Same as `Ve_n1/2` but the *difference* (`int1`/`int2` in `lightcurve`) — this is what actually multiplies `piE` in `s.ux`/`s.uy`. |

### `covarian` (`Bulge.h:452-506`) — Fisher-matrix workspace, one instance reused across all events

| Member | Meaning / units |
|---|---|
| `sign` | **Declared, never used anywhere in these 4 files.** `invert_matrix` uses a local `int s` shadow instead. Dead member. |
| `flagi` | Intended as a per-event "was the Fisher matrix well-conditioned" flag. In practice: `FisherM` sets it via `co.flagi =+ 1;` (`Bulge_LSST.cpp:941`) — this parses as `co.flagi = (+1)`, i.e. an **unconditional reset to exactly 1**, not an increment. All the code that would later set it negative on a bad/singular inversion (`Bulge_LSST.cpp:1073, 1076, 1087, 1202`) is inside `/* */` block comments, i.e. **currently dead**. Net effect: `co->flagi` is always exactly `1` after any `FisherM` call, and every `if (co->flagi > 0)` gate downstream (`593, 781`) is currently a no-op — the conditioning check the plan's Step C4 wants doesn't currently do anything. Flagged in §6. |
| `deter` | Determinant of the Fisher matrix (from GSL LU decomposition), used only to detect an *exactly*-zero determinant and then nudge the diagonal by `1e-10` as a crude regularization (`invert_matrix`, `Bulge_LSST.cpp:1525-1564`). Not stored per-event; this is the ad hoc precursor to the proper condition-number diagnostic Plan Step C4 asks for. |
| `sigmul1, sigmul2` | Propagated 1σ uncertainty on the lens proper-motion components `mul1, mul2`, combining the astrometric Fisher uncertainty with a `piE`/`tE`-driven geometric term. Feeds `resu[11-13]`. |
| `f1, f2` | Scratch intermediates in the `sigmul1/2` calculation (combines `resu[5]` (tetE), `resu[1]` (tE), and a `tan(xi)`/`cot(xi)` cross-term with correlation `corr1`). |
| `magw` | Scratch: perturbed model magnitude at one finite-difference step, inside `FisherM`'s photometric derivative loop. Not read outside that loop. |
| `derm1f, derm2f, dera1f, derb1f, dera2f, derb2f` | Central-difference derivative estimates (averaged from the `[2]`-element `derm1/derm2/dera1/derb1/dera2/derb2` scratch arrays) — the actual entries that go into the Fisher matrices. |
| `diff` | Scratch: the signed finite-difference step (`±Delta1[j]·sign`) currently being applied to whichever parameter is being perturbed. |
| `resu[15]` (`nq=15`) | Output vector of normalized 1σ errors / derived quantities. Only a subset survive into `EventRecord`/output — see index table below. |
| `derm1, derm2, dera1, derb1, dera2, derb2, bb` (all size 2) | Two-sided (`h=0,1`) finite-difference scratch, one evaluation per side of the derivative. `bb` specifically holds the two `fb` step sizes (see below). |
| `Era[Nx], Erb[Ny]` | Absolute 1σ uncertainties on the photometric (`u0,tE,fb,piE,xi`) and astrometric (`tetE,mus1,mus2,piE`) Fisher parameter vectors — `sqrt(diag(inverse Fisher matrix))`. |
| `Delta1[Nx], Delta2[Ny]` | Finite-difference step sizes per parameter — legacy-tuned hardcoded constants (`Bulge_LSST.cpp:952-955, 1104-1107`). Explicit target of Plan Step C3. |
| `inputA/inverA` (`Nx×Nx`), `inputB/inverB` (`Ny×Ny`) | The photometric and astrometric Fisher matrices and their GSL-inverted covariance matrices. |
| `summA, summB` | `summB` is computed once, live, via `gsl_blas_dgemm` (`Bulge_LSST.cpp:1195`) as `inverB·inputB` (should be ≈identity if the inversion is good) — but the code that would actually check it against identity is commented out, so the result is computed and discarded. `summA`'s equivalent computation is entirely inside a block comment, so `summA` is allocated but never populated at all. Both are vestigial scaffolding for the conditioning check the plan's Step C4 wants built properly. |

`covarian::resu[]` index meaning (only indices actually forwarded to `EventRecord` are marked ✓):

| idx | Meaning | In `EventRecord`? |
|---|---|---|
| 0 | `σ(u0)/u0` | ✓ (`resu0`) |
| 1 | `σ(tE)/tE` | ✓ (`resu1`) |
| 2 | `σ(fb)/fb` (LSST r) | ✓ (`resu2`) |
| 3 | `σ(piE)/piE`, photometric — then overwritten by `MIN(resu[3], resu[8])` | ✓ (`resu3`) |
| 4 | `σ(xi)/xi` | ✗ — computed, dropped |
| 5 | `σ(tetE)/tetE` | ✓ (`resu5`) |
| 6 | `σ(mus1)/mus1` | ✗ — computed, dropped |
| 7 | `σ(mus2)/mus2` | ✗ — computed, dropped |
| 8 | `σ(piE)/piE`, astrometric branch (used only to `MIN()` into idx 3) | ✗ — computed, dropped (folded into 3) |
| 9 | Lens mass fractional error, `sqrt(resu3²+resu5²)` | ✓ (`resu9`) |
| 10 | Lens distance fractional error | ✓ (`resu10`) |
| 11 | `σ(mul1)/mul1` | ✗ — computed, dropped |
| 12 | `σ(mul2)/mul2` | ✗ — computed, dropped |
| 13 | Lens proper-motion fractional error | ✓ (`resu13`) |
| 14 | Source proper-motion fractional error | ✓ (`resu14`) |

Indices 4, 6, 7, 8, 11, 12 are computed by `ErrorCal` every time but never persisted past the local
`co->resu[]` array — whether that's intentional (only the combined/headline quantities matter) or an
oversight is worth confirming before Plan Step D1 extends `EventRecord`, since D1 explicitly wants
`σ(piE)`, `σ(tE)`, `σ(tetE)`, `σ(Ml)` per-survey, and the astrometric-branch `piE` (idx 8) and `xi`
(idx 4) are currently thrown away.

### Other structs (brief — not flagged as needing deep coverage by Step A1, but included for completeness)

- **`CMD`**: 4 parallel column-store tables (thin disk, bulge, thick disk, halo), each sized `N1..N4`, holding `logT, mass, Mab[M]` (absolute magnitude per filter), `typ, cl, age`. Populated once by `read_cmd()` from `CMD/components/*.dat`.
- **`extin` / `Sightline` / `ExtinctionProfile`**: `NFILES=2518` fixed sightlines, each with a `NROWS=3686`-row extinction-vs-distance profile. Populated once by `readBayestar()`.
- **`lsst`**: per-visit LSST/Rubin baseline (`Nl=7373` rows from `BulgeBaseline.dat`) plus a `Na=96`-row mag→astrometric-error lookup table (`sigmaA_LSST.txt`, despite the "photometric error" framing in some comments — confirmed by usage: `ls.mag/err` feeds `errlsstA`, the astrometric error function).
- **`roman`**: per-visit Roman baseline (`NlRoman=309084` rows from `RomanBaseline.dat`) plus a `NaRoman=123`-row mag→error lookup (`sigma_roman.txt`) used for photometric error only (`errRomanM`) — no astrometric error table exists yet for Roman (explicit TODO, `Bulge.h:439-441`).
- **`EventRecord`**: the flat one-row-per-detected-event output struct assembled in `main()`. **Naming trap**: three of its fields store a value that is *not* what the field name suggests — `vsave` actually stores `vsave/s->mus` (a normalized ratio, not the raw running-velocity sum), `DeltaT` actually stores `l->DeltaT/s->errA` (a detection-significance-like ratio, not the raw deflection amplitude), and `murel` actually stores `l->murel*year` (converted to mas/yr, not the raw mas/day value the `lens` struct member holds). See `Bulge_LSST.cpp:616-627` vs. the struct declaration `Bulge.h:526-537`. Flagged in §6.

---

## 2. Main loop structure (`main()`, `Bulge_LSST.cpp:107-931`)

1. **Setup** (`107-213`): allocate the 8 top-level state objects; read `BulgeBaseline.dat` (LSST, 7373 visits) into `ls`, `sigmaA_LSST.txt` (astrometric error table) into `ls`, `sigma_roman.txt` (photometric error table) into `ro`, `RomanBaseline.dat` (309,084 visits) into `ro`, run `readBayestar()` (2518 extinction sightlines) into `ex`, run `read_cmd()` (4 CMD population tables) into `cm`.
2. **Sky-position loop** (`287-929`, nested `for lon`/`for lat`, step `dd=0.02°`): **currently hardcoded to a 0.1°×0.1° smoke-test patch** (`lon∈[0.5,0.6]`, `lat∈[-1.0,-0.9]`) — the real production range (`l1-wid..l2+wid` × `b1-wid..b2+wid`, the full bulge footprint) is present but commented out (`286, 291`). A field is skipped entirely if `lon<lx and lat>bx` (`294-296`), which excludes one corner of the bounding box to approximate an L-shaped/cross-shaped survey footprint — I can see *that* it does this but not *why* that particular corner from these 4 files alone; likely tied to the plan's deferred "Tile geometry" item.
   - **Per-position epoch matching** (`301-307`): `matchVisibleEpochs()` called once for LSST (`FoV=1.75°` radius) and once for Roman (`FoVRoman=0.28°`), producing `ndd`/`nddR` visible-epoch counts and `minc`/`mincR` minimum cadences. See §5.
   - **Per-position setup** (`309-325`): reset per-position histograms, compute `s->TET/FI`, run `Disk_model()` (builds the radial density profile for this sightline), find the nearest of the 2518 extinction sightlines.
   - **Per-star detection loop** (`323-685`, `do { … } while (icon<20 or nlens<5 or nerr<1.0)` — the plan's Step E2 stopping-criteria floors, currently lowered from the production values `850/150/2` (visible commented-out at `684`) for fast local iteration):
     - Draw a source (`func_source`) and a lens for it (`func_lens`); compute optical depth (`optical_depth`).
     - **Pre-selection** (`354-365`): peak-magnification test in all `M=7` filters → `fdet`; accept if `fdet>1` (detectable in ≥2 filters) **and** a random draw is `≤ s->blend[2]` (LSST r-band blend fraction — this is the LSST-only pre-selection bias Plan Step B2 targets).
     - **Light-curve time loop** (`376-552`, adaptive `dt`): for each timestep, call `lightcurve()`, compute per-filter lensed/unlensed magnitudes, then two near-identical but separately-cadenced branches:
       - *LSST branch* (`402-475`): if `tim` lands on the next LSST visible epoch, test peak magnitude against that filter's `[satu,thre]`; if in range, compute photometric+astrometric error, simulate a noisy magnitude, accumulate `chi1/chi2/chi3` (lensed / no-parallax / baseline χ²) and `chi1a/chi2a/chi3a` (astrometric equivalents), push a datum into `l->timn/magn/errm/soux/souy/erra/tele[ndw]` tagged `tele=0`, update the 3-consecutive-3σ run test (`flag0/1/2 → flag_det`).
       - *Roman branch* (`486-521`): parallel structure, F146-only (`fiR=6` fixed), pushes into the same arrays tagged `tele=1`, but is **not yet wired into `chi1/chi2/chi3`/`flag_det`** (explicit TODO, `477-485`) — this is exactly the per-survey detection bookkeeping split Plan Step B1 asks for.
       - *Adaptive `dt`* (`523-550`): steps at the smaller of LSST's and Roman's local next-epoch cadence, so the loop doesn't step over Roman's dense high-cadence epochs.
     - **Detection decision** (`559-613`): `dchiL=|chi3-chi1|` (lensing effect), `dchiP=|chi2-chi1|` (parallax effect), `dchiA=|chi2a-chi1a|` (astrometric deflection effect — computed but not currently part of the detection gate). If `FWHM<Tobs and dchiL>2·ndw and flag_det>0 and ndw>10`: mark `FFG[0]=1` ("Lensing" detected), call `FisherM(ndw)`, and (since `co->flagi` is currently always `1`, see §1) always run `ErrorCal` and append to `LpLMC*.dat`.
     - Push one `EventRecord` per star (detected or not) into the in-memory `records` vector (`616-627`).
   - **Per-field aggregation** (`687-928`): fold the field's `nstE/ndtE` into cumulative `NstE/NdtE` (currently always-zero, §1/§6); replay `records` to accumulate summary sums (`741-787`); normalize into per-field means (`789-812`); write one summary row to `fnGam` (`MapLMC*.dat`); print diagnostics; run ~25 `CHECK()` assertions (`872-926`) that assume every summary statistic is strictly positive — these encode a load-bearing assumption that the field always produces enough detections, tied to the Step E2 stopping criteria.
3. Outer loop closes (`929`), `return 0` (`931`).

## 3. Filter convention

`M = 7` (`Bulge.h:82`). Confirmed index order, three independent ways:
- `Bulge.h:88-102`'s per-filter constant arrays (`gama, seeing, msky, Cm, Dci, km, sigma, thre, satu, FWHM, lambda_um`) are all comma-lists of 6 LSST values followed by 1 Roman value, with `lambda_um[6]=1.464 µm` matching Roman's F146 bandpass.
- `func_source` (`Lensing.cpp:169-172`): `s.fb[0]=s.blend[2]` (LSST r) `// r-LSST`, `s.fb[1]=s.blend[6]` (Roman F146) `// F146` — explicit comments confirming index 2 = r, index 6 = F146.
- `read_cmd()`'s column-order comment (`helper.cpp:222-223`): `Roman_F146 LSST_u LSST_g LSST_r LSST_i LSST_z LSST_y` mapped to `Mab[6],Mab[0..5]` respectively — same convention.

Cross-checked against the Python side: `CMD/BolometricCorrection.py`'s `FILTER_ORDER = ["LSST_u","LSST_g","LSST_r","LSST_i","LSST_z","LSST_y","Roman_F146"]` matches exactly (indices 0-5 LSST ugrizy, 6 Roman F146) — and that file auto-parses `thre[]`/`satu[]` out of `Bulge.h` specifically to prevent the drift the plan worries about, rather than hardcoding a second copy.

## 4. The two indexing systems

- **By filter** (`i` = `0..M-1`): `s.blend[i]`, `s.magb[i]`, `s.nsbl[i]`, `s.Ai[i]`, `s.Mab[i]`, `s.Map[i]`. Indices 0-5 = LSST `ugrizy`, index 6 = Roman F146.
- **By telescope** (`tt` = `0` or `1`): `s.fb[tt]`, `s.mbs[tt]`. `tt=0` = Rubin (fed from filter index 2, r-band), `tt=1` = Roman (fed from filter index 6, F146). Set once per star in `func_source` (`Lensing.cpp:169-172`), then **actively perturbed in place** during `FisherM`'s finite-difference loop (`s.fb[tt] += co.diff; …; s.fb[tt] -= co.diff;`) — `tt` there comes from `l.tele[i]`, the per-datum telescope tag, so the same `fb[]`/`mbs[]` pair is reused across both instruments' data within one Fisher-matrix accumulation.

Every datum in the per-event light-curve buffer (`l->timn/magn/errm/soux/souy/erra`) is tagged in the
parallel `l->tele[]` array: `0`=LSST/Rubin, `1`=Roman/F146. This is what lets `FisherM` pick the
right `s.fb[tt]`/`s.mbs[tt]` per datum, and is the hook any per-survey bookkeeping (Plan Phase B)
would key off.

## 5. Epoch-matching machinery

`matchVisibleEpochs()` (`Bulge_LSST.cpp:32-82`) is instrument-agnostic: given a sky position, a field-
of-view radius, and one instrument's full `l[]/b[]/tim[]` arrays, it walks the array once (assumes
pre-sorted-by-time, per-instrument, within-sky-cone) and fills `ct[]` with the indices (into that
instrument's own arrays) of every epoch within `fov` of the position, returning the count `ndd` and
the smallest observed gap `minCadence`. Tied timestamps (`cade<=0`, e.g. two overlapping adjacent
Roman fields sharing a schedule) are logged and skipped, keeping the first.

It's called twice per sky position (`303-304`): once as `("LSST", …, FoV=1.75°, ls->l/b/tim, Nl, ls->ct, minc)` → `ndd`, once as `("Roman", …, FoVRoman=0.28°, ro->l/b/tim, NlRoman, ro->ct, mincR)` → `nddR`.

Inside the time loop, each instrument gets its own cursor into its own `ct[]`: `gi`/`sq` for LSST
(`sq = ls->ct[gi]`, the current candidate epoch's row index in `ls`'s full arrays), `giR`/`sqR` for
Roman (same pattern against `ro`). A datum is only accepted for an instrument when `tim` has reached
that instrument's next `ct[]`-listed epoch (`tim >= ls->tim[sq]` / `tim >= ro->tim[sqR]`) *and* the
star's magnitude at that instant is within that instrument's `[satu,thre]` band; the cursor (`gi`/
`giR`) advances independently of whether the magnitude test passed, so a too-faint/too-bright epoch
is still "consumed" from the schedule.

`dt`, the adaptive light-curve timestep, is set every iteration to `min(cade, cadeR)` — `cade` = the
gap to LSST's *next* scheduled visible epoch if currently inside LSST's epoch window, else a flat 3-
day fallback; `cadeR` = the analogous quantity for Roman, or `cade` itself (i.e. "don't constrain")
if Roman isn't currently in-window. This is what lets the loop take Roman's ~12-minute steps during
its high-cadence seasons without silently skipping over them at LSST's coarser cadence, and without
permanently locking the whole light curve to Roman's cadence once Roman's season ends.

## 6. Variables/behaviors I could not confidently identify, or that look like bugs

Ranked roughly by how much they'd matter if wrong:

1. **`l->nstE`/`l->ndtE` (and the whole `Ns*/Nd*` histogram family) are never incremented anywhere
   in these 4 files.** They're zeroed per field (`314`), then at `690` `l->ndtE[i]` is *overwritten*
   (not incremented) with `ndtE[i]/(nstE[i]+eps)`, which is `0/eps = 0` since neither side was ever
   touched. `EFF` (`779`) accumulates `l->ndtE[gg]/(l.tE/year)` over every record — i.e. it
   accumulates zero every time — so `EFF` should evaluate to exactly `0.0` at the end of any field
   that reaches the summary stage, which would fail `CHECK(EFF > 0.0)` at `872`. A commented-out
   block at `241-247` shows these arrays used to be populated by reading an external file
   (`EfLMC2B.dat`) that is no longer read anywhere — my best guess is this is leftover plumbing from
   before that file-read was removed, not something that was ever re-wired to fill the histograms a
   different way. I have not run the binary to confirm this actually throws (no data files are
   present in this checkout — see `CLAUDE.md`'s note on `.gitignore`), so I'm flagging this as a
   high-confidence read of the code rather than a confirmed runtime observation. Worth confirming
   before doing anything else with the detection/yield numbers, since `Gamma` and `Neven` (expected
   event rate/yield) both derive from `EFF`.
2. **`l->Nhalo[2]`/`l->Nself[2]` are read (divided, printed) but never assigned anywhere**, and as
   `std::array<double,2>` members with no constructor initializer, they hold **indeterminate values**
   — this is undefined behavior (reading uninitialized memory), not merely "always zero." Same
   commented-out-file-read origin as item 1.
3. **`co.flagi =+ 1;` (`941`) parses as `= (+1)`, not `+= 1`.** Combined with all the code that would
   set it negative on a bad Fisher-matrix conditioning being commented out (`1073,1076,1087,1202`),
   `co->flagi` is currently always exactly `1` after any `FisherM` call — the `if (co->flagi > 0)`
   gates at `593`/`781` are currently unconditional. This directly bears on Plan Step C4
   ("replace crash-or-garbage behaviour on singular matrices with an explicit 'not characterizable'
   outcome") — right now nothing marks a singular/ill-conditioned event as such.
4. **`EventRecord` field-name/value mismatches**: `vsave` stores `vsave/s->mus`, `DeltaT` stores
   `l->DeltaT/s->errA`, `murel` stores `l->murel*year` — see the `EventRecord` row in §1. Anyone
   reading `r.vsave` expecting the raw quantity will get a normalized ratio instead.
5. **`s.cl`, `s.typ`** (Besançon CMD "CL"/"Type" columns): assigned from the CMD tables but never
   read back anywhere in these 4 files. I don't know what the numeric codes mean (evolutionary stage?
   luminosity class?) — `read_cmd`'s `CHECK(... <= 7)` / `CHECK(... <= 9.0)` bound them but don't
   document them. Not urgent (they're not used), but worth a one-line note if anyone later wants to
   use them.
6. **`s.ros`** (finite-source parameter): computed, never used — the light curve is point-source
   throughout. Confirms finite-source effects are out of scope for this simulation, but I can't tell
   from the code alone whether that's a deliberate scope decision or leftover from a fork that dropped
   finite-source modeling.
7. **Dead/unused members with no downstream effect**: `source::xv,yv,zv`, `source::rho_stars`,
   `source::nstarti`, `covarian::sign`. None of these affect behavior; listed for completeness in
   case any of them was meant to be wired up and simply wasn't.
8. **The `lon<lx and lat>bx` field-skip** (`294-296`): I can see it excludes one corner of the
   bounding box but not the geometric/survey-footprint reasoning behind exactly that corner. Likely
   answered by the whitepaper's Locations section (deferred in the plan).
9. **`betal/betas/deltal/deltas/deltao`** geometric derivation: I'm confident about their role
   (sky-plane basis-rotation angles for lens/source/Sun) but not confident I've correctly reconstructed
   the exact spherical-trig convention from the code alone — see §1.
