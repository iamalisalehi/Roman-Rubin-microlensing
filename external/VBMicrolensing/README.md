# VBMicrolensing (vendored)

[VBMicrolensing](https://github.com/valboz/VBMicrolensing) is a C++ microlensing library by Valerio
Bozza (University of Salerno), with Vito Saggese and Jiyuan Zhang, the successor of VBBinaryLensing.
It computes magnifications by single, binary and multiple lenses with finite, limb-darkened sources by
contour integration, image centroids, and light curves with parallax and orbital motion. The simulator
takes from it the finite-source (ESPL) magnification, the astrometric centroid and the annual parallax
from a tabulated Sun ephemeris.

## Version

| | |
|---|---|
| Upstream | https://github.com/valboz/VBMicrolensing |
| Version | v5.6.1 (tag `v5.6.1`) |
| Commit | `db4f95267d74ad5f4dc592ab1e6d46c05f6cc1b7`, 2026-10-08, "Merge pull request #76 from valboz/dev" |

## Files

Copied from that commit **unmodified except two local additions** (see
[Local modifications](#local-modifications)):

| Here | Upstream | What it is |
|---|---|---|
| `VBMicrolensingLibrary.cpp` | `VBMicrolensing/lib/VBMicrolensingLibrary.cpp` | the library |
| `VBMicrolensingLibrary.h` | `VBMicrolensing/lib/VBMicrolensingLibrary.h` | its header |
| `data/ESPL.tbl` | `VBMicrolensing/data/ESPL.tbl` | pre-computed uniform-disc point-lens magnification table (binary) |
| `data/SunEphemeris.txt` | `VBMicrolensing/data/SunEphemeris.txt` | JPL Horizons geocentric Sun ephemeris (DE441), daily, 1990-04-18 to 2050-07-12 |
| `LICENSE` | `LICENSE` | GNU LGPL v3 |

`COPYING` is not from upstream: it is the GNU GPL v3 text (https://www.gnu.org/licenses/gpl-3.0.txt),
which the LGPL incorporates and asks to accompany the library. Not taken: the Python bindings
(`lib/python_bindings.cpp`), the sample satellite and event-coordinate tables, the docs and examples.

To check the copies against upstream (no output = identical; both library files are expected to differ, see below):

```bash
git clone -q https://github.com/valboz/VBMicrolensing /tmp/vbm
git -C /tmp/vbm checkout -q db4f95267d74ad5f4dc592ab1e6d46c05f6cc1b7
for f in VBMicrolensingLibrary.cpp VBMicrolensingLibrary.h; do cmp /tmp/vbm/VBMicrolensing/lib/$f external/VBMicrolensing/$f; done
for f in ESPL.tbl SunEphemeris.txt; do cmp /tmp/vbm/VBMicrolensing/data/$f external/VBMicrolensing/data/$f; done
cmp /tmp/vbm/LICENSE external/VBMicrolensing/LICENSE
```

The `cmp` of the two library files reports a difference, by design. To see exactly what differs (`diff` prints the
lines with the files' CRLF line endings):

```bash
diff /tmp/vbm/VBMicrolensing/lib/VBMicrolensingLibrary.h external/VBMicrolensing/VBMicrolensingLibrary.h
diff /tmp/vbm/VBMicrolensing/lib/VBMicrolensingLibrary.cpp external/VBMicrolensing/VBMicrolensingLibrary.cpp
```

The `.cpp` diff is `>` lines only, every one inside a block opened by a comment `// LOCAL ADDITION 2 (2026-11-13)` and
closed by `// END LOCAL ADDITION 2`, or a bare brace or one-line call marked `// LOCAL ADDITION 2`. The header diff is
the two added blocks, the only changed upstream lines being one sentence of addition 1's own comment (it said
"the only change"). No upstream line is edited, re-indented or removed.

## Local modifications

Two additions, made for the Roman-Rubin simulator; the data files are as upstream. Addition 1 is in the header only,
addition 2 is in the header and, as added lines around upstream's, in the `.cpp`.

### Addition 1 (2026-10-11): `EarthVelocityOnSky`

| What | Where |
|---|---|
| `void EarthVelocityOnSky(double t0, double* v)`, an inline public method of `VBMicrolensing` | right after `AreCoordinatesSet()`, in a comment block marked "LOCAL ADDITION" |

It returns the Earth's velocity at `t0` projected on the sky at the target, along celestial South and West
in AU/day: the private member `vt0` that `ComputeParallax` builds for `t0_par = t0`, the geocentric frame
(Gould 2004) in which every `*Parallax` light curve is computed. The simulator draws events with
heliocentric proper motions and needs this velocity to convert them to that frame (mu_geo = mu_hel -
v_E,perp * pi_rel); without the method `vt0` is out of reach, and re-deriving it outside the library would
duplicate the ephemeris interpolation the light curve uses. The method invalidates the parallax caches
before and after and restores `t0_par_fixed` and `t0_par`, so it changes no result of any existing
function (`tests/vbm_test.cpp` passes unchanged). Only `src/events/vbm_model.cpp` calls it.

The LGPL v3 incorporates the GPL v3, whose section 5(a) asks a modified work to carry prominent notices
stating that it was modified and giving a relevant date: the comment block in the header and this section
are those notices. The method is under the library's licence, like the rest of the directory.

### Addition 2 (2026-11-13): the observers' positions once per sightline epoch

| What | Where |
|---|---|
| `void SetParallaxCache(bool on)` and `void SetObserverGrid(const double* t, int n)`, inline public methods; counters `pc_hits`, `pc_misses`; private storage and `PCClear`, `PCSame`, `PCFind` | header, right after addition 1, in a comment block marked "LOCAL ADDITION 2" |
| in `ComputeParallax`: `int pc_k = -1;`; an `if (...) {stored Ear} else {` around the Earth's interpolation at `t`; an `if (pc_k >= 0) {stored Spit} else {` around the satellite's search | `.cpp`, marked "LOCAL ADDITION 2" |
| `PCClear();` in `LoadSunTable` and in the two-argument `SetObjectCoordinates` | `.cpp`, one line each |

`ComputeParallax(t, t0)` works out at every call where the observers are at `t`: the Earth's ecliptic position, by
linear interpolation in the daily Sun table, and a satellite's, by a binary search of its table and interpolation. Neither
depends on the target or on the event, only on `t` and the tables, and a simulation evaluates the same epochs of a
sightline (about 10^5 of them) for every one of its events and again ~20 times per event for the Fisher matrices.
`SetObserverGrid` computes both once for the sorted epochs, by the very expressions `ComputeParallax` uses, and
stores them with the value those expressions would give `parallaxextrapolation`. When `t` equals a grid time bit
for bit (and the cache is on, `parallaxephemeris` is set, `t_in_HJD` is false) `ComputeParallax` writes the stored
doubles where it would have computed them and runs every other line, unchanged and in its order: the light-travel
term, `Ehel`, `Et`, the projection of the satellite on `rad` and `tang`, everything that involves `t0` or the
target. A time off the grid takes upstream's lines. The cache is off until `SetParallaxCache(true)`, so the library
without that call is upstream's. A table reload empties the grid; a change of target does not (nothing stored
depends on it). Nothing is approximated: same operands, same operations, same order, the same bits, and the cost
of a miss is a few comparisons. `tests/vbm_cache_test.cpp` (`make vbmcachetest`) checks it against an object that
never uses the cache, with `memcmp`, for both observers, both source models, with and without centroids, on
the real visit times, with off-grid times, a stale grid, a table reload and times outside the tables; the
adapter (`src/events/vbm_model.cpp`) turns the cache on and hands each sightline's grid in (`vbmSetObserverGrid`).
Bit identity assumes the compiler does not fuse the interpolation's multiply-add differently from `ComputeParallax`'s
own copy of it: not so for the flags in the Makefile (no `-march`, no `-ffast-math`), and the test is the check
whenever they change.

## Build and test

The `Makefile` compiles `VBMicrolensingLibrary.cpp` with its own flags (`-O2 -g -std=c++17 -w`: our
optimisation, never `-ffast-math`, and no warnings for code we do not edit) and links it into `./roman`
and the test programs. Our code includes the header as a system header (`-isystem`), so `-Wall -Wextra`
apply to our code only. The data paths are `PATH_VBM_ESPL_TABLE` and `PATH_VBM_SUN_TABLE` in
`config/parameters.h`; like every data path they are relative to the repository root.

`tests/vbm_test.cpp` pins the library against references that do not use it (`tests/vbm_reference.h`,
written by `tests/vbm_reference.py`): finite-source magnification against direct integrals over the
source disc, the point-source magnification against our own formula bit for bit, and the parallax
source track against astropy's ephemeris.

```bash
make vbmtest && ./vbmtest
```

## Updating

1. Copy the same five files from the newer upstream commit, then re-apply the local addition to the header
   (the block shown by the `diff` command above) and update its comment if upstream changed `ComputeParallax`.
2. Run `make vbmtest && ./vbmtest`, `make modeltest && ./modeltest` (which checks the Earth velocity against
   astropy) and the other tests. Do not regenerate `tests/vbm_reference.h` to
   make an update pass: the references are independent of the library, so a failure means its answers
   moved, which has to be understood first.
3. Record the new commit hash, date and version above.

## Licence

VBMicrolensing is distributed under the GNU Lesser General Public License v3 (`LICENSE`, which
incorporates the GNU GPL v3 in `COPYING`); that licence covers the files in this directory only. It is
compiled from source alongside the simulator, which is MIT-licensed, so anyone can modify these files and
rebuild with `make`, which is what the LGPL asks of a work that combines with the library.

## Citing

Upstream asks that scientific use cite the works relevant to the study. Those behind what the simulator
uses are:

- V. Bozza, MNRAS 408 (2010) 2188: the contour-integration algorithm and limb darkening;
- V. Bozza, E. Bachelet, F. Bartolic, T. Heintz, A. Hoag, M. Hundertmark, MNRAS 479 (2018) 5157: the
  extended-source point-lens methods (ESPLMag2) and the switch from point to extended source;
- V. Bozza, E. Khalouei, E. Bachelet, MNRAS 505 (2021) 126: astrometry and generalized limb darkening.

Upstream also lists V. Bozza, V. Saggese, G. Covone, P. Rota, J. Zhang, A&A 694 (2025) A219 (multiple
lenses) and J. Skowron, A. Gould, arXiv:1203.1034 (polynomial root finding), which apply to the
binary- and multiple-lens calculations.
