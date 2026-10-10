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

Copied from that commit **unmodified**:

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

To check that the files are unmodified (no output = identical):

```bash
git clone -q https://github.com/valboz/VBMicrolensing /tmp/vbm
git -C /tmp/vbm checkout -q db4f95267d74ad5f4dc592ab1e6d46c05f6cc1b7
for f in VBMicrolensingLibrary.cpp VBMicrolensingLibrary.h; do cmp /tmp/vbm/VBMicrolensing/lib/$f external/VBMicrolensing/$f; done
for f in ESPL.tbl SunEphemeris.txt; do cmp /tmp/vbm/VBMicrolensing/data/$f external/VBMicrolensing/data/$f; done
cmp /tmp/vbm/LICENSE external/VBMicrolensing/LICENSE
```

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

1. Copy the same five files, unmodified, from the newer upstream commit.
2. Run `make vbmtest && ./vbmtest` and the other tests. Do not regenerate `tests/vbm_reference.h` to
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
