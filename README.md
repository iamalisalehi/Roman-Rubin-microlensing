# Roman + Rubin microlensing forecasts for the Galactic bulge

A Monte Carlo simulator that forecasts what the Nancy Grace Roman Space Telescope's Galactic Bulge
Time Domain Survey (GBTDS) and the Vera C. Rubin Observatory's LSST will detect and measure in
microlensing toward the Galactic bulge, separately and in combination.

For each sightline it

1. draws source and lens stars from a Besançon Galaxy-model catalogue (or from Kroupa-IMF, remnant,
   black-hole or neutron-star lens populations), with 3D dust extinction at each star's distance;
2. builds per-band light curves on the two surveys' cadences: Rubin from an OpSim baseline database,
   Roman from the GBTDS season and field layout;
3. applies a detection test, then forecasts parameter precision with a Fisher-matrix analysis,
   photometric (u0, tE, piE, blending, ...) and astrometric (thetaE, source proper motion, piE),
   reported **separately for Roman alone, Rubin alone and the joint data set**.

## Requirements

- A C++17 compiler (`g++`), `make`, and the [GNU Scientific Library](https://www.gnu.org/software/gsl/)
  (`libgsl-dev` on Debian/Ubuntu, `gsl` on Homebrew, `module load gsl` on most clusters).
  `pipeline/pipeline.sh` can build GSL locally if none is installed (`GSL=build` in the config).
- Python ≥ 3.10 with the packages in [`requirements.txt`](requirements.txt). The pipeline's `setup`
  stage creates a virtual environment and installs them.
- `curl`, `tar`, `xz` for the downloads. Optionally Slurm for cluster runs.

## Quick start

```bash
git clone https://github.com/iamalisalehi/Roman-Rubin-microlensing.git
cd Roman-Rubin-microlensing
cp pipeline/config.example.sh my_run.sh     # set BESANCON_CATALOGUE, EXT_TABLES, POPULATIONS, ...
pipeline/pipeline.sh my_run.sh check        # lists what is missing; changes nothing
pipeline/pipeline.sh my_run.sh all          # setup, fetch, prep, sim, merge, analyze
```

The pipeline runs on a laptop or as a Slurm array. It downloads what it can (Rubin baseline, MIST
bolometric corrections, dust maps), builds the star lists, visit lists and extinction tables, splits
the sightline scan into chunks, merges the outputs and makes the figures. The Besançon catalogue has
to be requested from the [model's web form](https://model.obs-besancon.fr/); the parameters are in
[`CMD/Besancon/RERUN_bos10.md`](CMD/Besancon/RERUN_bos10.md). Full details:
[`pipeline/README.md`](pipeline/README.md).

To build and test only the C++ core:

```bash
make                                  # builds ./roman (run it from the repository root)
make fishertest && ./fishertest       # Fisher-matrix regression test, no data files needed
make extinctiontest && ./extinctiontest
make vbmtest && ./vbmtest             # VBMicrolensing against independent references
```

## Repository layout

| Path | Contents |
|---|---|
| `src/`, `include/` | C++ simulator, one module per directory: `sim/` (per-sightline and per-event Monte Carlo), `galaxy/` (density, kinematics, star lists, extinction), `events/` (sources, lenses, light curves), `surveys/` (footprints, schedules, noise models), `fisher/` (Fisher matrices), `run/` (I/O and command line) |
| `config/parameters.h` | Every model and survey choice, compile-time (edit, then `make`) |
| `config/data_products.h` | Generated from the data files by `tools/sync_data_products.py`; do not edit |
| `pipeline/` | End-to-end driver, example configuration, input validation, chunk merging |
| `CMD/` | Builds the star lists and lens mass-luminosity table from a Besançon catalogue |
| `Baseline/` | Builds the Rubin and Roman visit lists; vendored GBTDS field layout |
| `maps.py` | Builds the 3D extinction tables from the DECaPS and Marshall dust maps |
| `analysis/` | Yield, gap-filling and precision analyses and figures (all through `analysis/romanlib.py`) |
| `samples/` | Specifications of illustrative events whose full light curves `./roman` writes out |
| `tests/` | C++ regression and unit tests; data-quality checks for the generated inputs |
| `files/` | Small static inputs (noise curves); run outputs are written under `files/MONTLMC/files/` |
| `external/` | Third-party code compiled into the simulator, unmodified: [VBMicrolensing](external/VBMicrolensing/README.md) |
| `Whitepaper/` | LaTeX source of the write-up |

Data products are generated, not committed.

## Analysis conventions

The Monte Carlo does not draw events in proportion to the event rate, so pooled statistics must be
weighted. The analysis scripts require a sightline map (`--map`) or an explicit `--unweighted`, and
quote every weighted number with its Kish effective sample size. A Fisher entry of `-1.0` means
"not measured" and is never averaged.

## Acknowledgements

The simulator builds on earlier microlensing simulation code by Prof. Sedighe Sajadian, whom I thank
for her guidance and for the foundation it provided.

## License

MIT. See [`LICENSE`](LICENSE). The exception is `external/VBMicrolensing/`, which is VBMicrolensing
(Bozza et al.) under the GNU LGPL v3; see [its directory](external/VBMicrolensing/README.md).
