# Running the whole pipeline on another computer or a cluster

One driver, `pipeline/pipeline.sh`, takes you from a fresh `git clone` to the analysis plots. It is bash plus
three small Python helpers (no workflow framework). Every stage is idempotent: it records what it was run on
(input paths, sizes, modification times, options) in a small stamp file and skips itself when nothing changed.
`FORCE=1` redoes a stage.

```
pipeline/pipeline.sh CONFIG STAGE
```

`CONFIG` is a copy of `pipeline/config.example.sh` with your choices (every variable is commented there).
Stages: `check` `setup` `fetch` `prep` `sim` `merge` `analyze` `status`, and `all` = check, setup, fetch, prep,
sim, merge, analyze. `chunk POP K` is internal (what the schedulers call). The script can be started from any
directory.

## Your six steps, as commands

| # | You said | Command |
|---|---|---|
| 1 | Download a Besancon simulation | download it yourself (below); `BESANCON_CATALOGUE=<path>` in the config |
| 2 | Download a Rubin baseline | `fetch` downloads v5.1.0 if `OPSIM_DB` is missing; for another baseline put its `.db` at `OPSIM_DB` |
| 3 | Set parameters (populations, stride, ...) | edit the config: `POPULATIONS`, `STRIDE`, `STRIDE_ROMAN`, `EVENTS`, `LENSES`, ... |
| 4 | Run the preprocessing Python codes | `setup`, `fetch`, `prep` |
| 5 | Run the main C++ simulation | `sim` (then `merge`) |
| 6 | Run the analysis codes | `analyze` |

```bash
git clone <repo> && cd Roman
cp pipeline/config.example.sh my_run.sh            # edit: BESANCON_CATALOGUE, EXT_TABLES, POPULATIONS, ...
pipeline/pipeline.sh my_run.sh check               # what is missing? changes nothing
pipeline/pipeline.sh my_run.sh all                 # everything, on this machine
```

### What to download

1. **Besancon file** (not downloadable by script: it comes from the web form,
   <https://model.obs-besancon.fr/>, model 1612). The parameters that matter are in
   `CMD/Besancon/RERUN_bos10.md` (the full form block is there; do not change the rest of it). In short:
   - **ONE field**: (l, b) = (0.5, -1.4), distance 0 to 12 kpc, solid angle 0.05 deg^2 (about 10 M stars, 3.7 GB).
     The simulator takes one catalogue for every sightline; it draws star kinds from it, not positions.
   - **Extinction off** (Av = 0, 0 clouds): the simulator applies its own 3D dust at each star's distance.
   - **Kinematics off** (the simulator draws its own).
   - **No magnitude or colour limits** (every band -99 to 99; AbsMag [-10, 30[). A V limit silently removes the
     faint dwarfs and white dwarfs of the disc.
   - **Every error law = 0** (Teff, logg, [M/H], [a/Fe], age, mass, and the photometric errors). Besancon's default
     noise makes the catalogue mass a noisy label, and the lens-mass tables bin by it.
   - Age range [0, 15], mass [0, 90], spectral types O0.0 - D5.0.
   - The output must be the usual text file whose first line is a `#` header with the 38 column names (it must
     contain `Teff logg Pop Age Mass Mbol [M/H] [a/Fe] CL Typ`). `check` and `prep` verify this before spending
     minutes; a few damaged rows (wrong field count) are skipped by the builder, as in bos10.
2. **Rubin baseline**: an OpSim database (sqlite, table `observations`). `fetch` gets v5.1.0 (0.8 GB) from
   `OPSIM_URL` if the file is missing. For another baseline, download its `.db` and set `OPSIM_DB`; `check` verifies
   the table and the columns `Baseline/readbaselineBulge.py` selects, and prints the time span against the
   simulation clock (see "Another baseline" below).
3. **Extinction tables** (`EXT_TABLES`, required): path to a prebuilt `files/ext/ext_tables.dat` (47 MB; copy it
   with `ext_provenance.json` if that sits next to it), or `build`. The table depends only on the survey geometry
   (not on the Besancon file or the OpSim db), so copying one is the sensible default. `build` runs `maps.py`:
   `fetch` downloads the DECaPS mean map (7 GB) and Marshall's map into `dustmaps/`, then ~35 minutes of queries.
4. Fetched for you: MIST v2 bolometric-correction tables (`LSST.txz` 2.9 MB, `Roman.txz` 3.9 MB, from
   mist.science) into `CMD/Rubin`, `CMD/Roman`; the LSSTCam focal-plane map (`Baseline/lsstcam_fov/fov_map.txt`).
   Vendored in git: `files/sigmaA_LSST.txt`, `files/sigma_roman.txt`, `CMD/empirical/EEM_dwarf_UBVIJHK_colors_Teff.txt`.

## What each stage does

- **check**: tools (g++, make, GSL headers, awk, curl, tar, xz, python3, sbatch if `SCHEDULER=slurm`), your inputs, the
  venv and its modules, and every product as `current` / `missing` / `stale`. Exit status is non-zero if anything is
  missing or stale.
- **setup**: creates the venv named by `PYTHON` (default `.roman/bin/python`) and installs
  `requirements.txt` (numpy, scipy, pandas, astropy, matplotlib, dustmaps; minimum versions only). If the GSL
  headers are not found (`GSL=auto`, the default) or `GSL=build`, it also downloads GSL 2.8 from gnu.org and builds a
  static copy in `deps/gsl/` (~5 min, no root); the Makefile uses `deps/gsl/` whenever it exists.
- **fetch**: downloads the above. Never the Besancon file (it prints where the form parameters are).
- **prep**: star lists `CMD/components/{thin_disk,bulge,thick_disk,halo}.dat` and `lens_ml.dat` (from the Besancon
  file, `CATALOGUE_DWARFS`, `CATALOGUE_FILL`), `Baseline/BulgeBaseline.dat` (Rubin visits from the OpSim db),
  `Baseline/RomanBaseline.dat` (`ROMAN_MISSION_START`), the extinction table, `config/data_products.h` (row counts
  and mean masses, measured from the files), and finally `make roman`. The binary has the row counts compiled
  in, so a changed data file always means a rebuild; that is why prep does the steps in this order.
- **sim**: for each population, a dry run counts the sightlines N (2013 at the production grid), and the scan
  `[0, N)` is split into chunks of `CHUNK_SIZE` sightlines. Each chunk is an independent `./roman --start-index a
  --end-index b` in its own run directory `RUN_ROOT/<pop>/chunks/NNNN/` (shared inputs symlinked, its own outputs).
  The binary is copied once to `RUN_ROOT/bin/roman.<git describe>` so a later `make` cannot change a run in
  progress. A chunk is done when `./roman` exited 0: only then is a `DONE` file written. Re-running `sim` skips
  the done chunks and redoes the rest from scratch. `RUN_ROOT/pipeline_provenance.txt` records the config, the
  git commit and the md5 of every data product, at the first `sim`.
- **merge**: per population, refuses (listing the chunks) if any chunk lacks `DONE` or the ranges do not tile
  `[0, N)`; then concatenates every output file in chunk order into `RUN_ROOT/<pop>/` (dropping the repeated `#`
  header lines) and checks sizes and line counts. Because every sightline re-seeds from (seed, index), the
  result is byte for byte what one unchunked run writes for `test<tag>.dat`, `h3_pair.dat`, `MapLMC<tag>.dat`
  and `LpLMC<tag>.dat` (tested on the stub, both populations). **Two files are not exact: `EfLMC<tag>.dat` and
  `EfLMC<tag>B.dat`** hold a *running* detection efficiency (each block is accumulated over all sightlines so
  far, only percentages are written), which restarts at the start of every chunk; the merged file keeps its
  shape (one block per sightline) but a block is cumulative over its own chunk only, so its last block covers the
  last chunk, not the scan. Nothing in `analysis/` reads these files today; do not take an efficiency curve from a
  chunked run. `run.log` is the chunk logs in order (each chunk keeps its own `RUN TOTALS` block).
  `run_provenance.txt` is chunk 0's, with `start_index 0`, `end_index N` and `merged_chunks K` set and the
  end-of-run sightline counts and areas summed over the chunks (they equal the unchunked values).
  `KEEP_CHUNKS=0` deletes the chunks afterwards.
- **analyze** (`ANALYSES`): `yields` = y1 per population and y3 on their concatenation (y3's `eta` column needs a
  population named `bulge`; without one it is NaN, and a `besancon` run does not count as `bulge`); `figures` = per
  population f1 (`--fields-only`, on the full table: it also counts the non-detected draws), then, on a
  detection-only extract (`test<tag>_detJ.dat`: header kept, rows with `detJ == 1`, the column found from the
  header; the numbers equal those from the full table, only the printed N_eff counts detections), f2 (tE and
  piE), f3, f4 (footprint and `--scope all`), h5, and h3 (if `PAIR_SATELLITE=1`); then across populations p6
  (`--detections-only`) and p7. Each script's output goes to a `.log` next to its results in
  `RUN_ROOT/analysis/`. A failing script is reported (name and log) and the others still run; the stage then exits
  non-zero.
- **status**: chunks done / running / pending per population, CPU used, sizes.

## Laptop example

```bash
cp pipeline/config.example.sh my_run.sh
#   BESANCON_CATALOGUE=/data/besancon/my_field   EXT_TABLES=/data/ext/ext_tables.dat
#   POPULATIONS="bulge ns"   RUN_ROOT=runs/laptop_1   JOBS=3        (each ./roman uses ~0.7 GB of RAM)
pipeline/pipeline.sh my_run.sh all
pipeline/pipeline.sh my_run.sh status        # from another terminal
```

If it is interrupted (laptop closed, Ctrl-C, reboot), run the same command again: finished chunks are kept, the
chunk that was running is redone, and the earlier stages skip. A smoke test of the whole chain (about 10
minutes for two populations on two cores): a new `RUN_ROOT` and
`STUB=1 STRIDE=1 STRIDE_ROMAN=1 EVENTS=10 LENSES=2 MAXDRAWS=3000 CHUNK_SIZE=10` (the 0.1x0.1 deg patch is 36
sightlines; the stub fails with the production strides).

## Slurm example

The compute nodes must see the clone (a shared file system). `setup`, `fetch` and `prep` are serial and need
network, a few GB of RAM and ~30 minutes: run them on the login node or in an interactive allocation
(`srun --mem=6G --time=02:00:00 --pty bash`). Then `sim` submits the jobs.

```bash
#   SCHEDULER=slurm
#   SLURM_OPTS="--account=myproj --partition=normal --time=03:00:00 --mem=2G --cpus-per-task=1"
#   SLURM_OPTS_FINISH="--account=myproj --partition=normal --time=04:00:00 --mem=12G"
for s in check setup fetch prep; do pipeline/pipeline.sh my_run.sh $s || break; done
SUBMIT=0 pipeline/pipeline.sh my_run.sh sim   # only writes RUN_ROOT/slurm/<pop>.sbatch and finish.sbatch
pipeline/pipeline.sh my_run.sh sim            # submits: one array job per population over the pending chunks,
                                              # then a merge+analyze job with --dependency=afterok on the arrays
pipeline/pipeline.sh my_run.sh status
```

(`all` with `SCHEDULER=slurm` does the same and stops after submitting.) If a chunk fails the finish job never
starts; look at `RUN_ROOT/slurm/*.out` and `RUN_ROOT/<pop>/chunks/NNNN/run.log`, then run `sim` again (it resubmits
only the chunks that are not `DONE`) and submit `RUN_ROOT/slurm/finish.sbatch` yourself. Do not run `sim` while
its jobs are still queued (it refuses when it sees your `rm_*` jobs in `squeue`). Jobs read a snapshot of the
config (`RUN_ROOT/slurm/config.sh`), so editing your config afterwards does not change queued jobs.

## Time, memory, disk (production flags: `--events 300 --lenses 50 --stride 10 --stride-roman 5 --pair-satellite`)

- `./roman`: ~18-19 CPU-hours per population, one thread each, ~0.7 GB RAM each; 2013 sightlines, of which the
  ~160 on a Roman detector cost 4-6 CPU-min each and the others 0.1-0.3. Each chunk first reads the inputs (~25 s).
  Wall time = CPU time / `JOBS` (or the number of array tasks the cluster runs at once).
- `prep`: star lists ~8 min and 1.7 GB RAM; Rubin visit list ~3 min; the rest is seconds; the first `make` ~1-2 min.
  Downloads: ~7 MB MIST, 0.8 GB OpSim, 7 GB DECaPS only for `EXT_TABLES=build`.
- Disk: star lists 0.6 GB; per-event tables are large: about 6.0 GB (bulge), 4.1 GB (ns), 3.2 GB (bh) at
  production flags, and they exist twice (chunks and merged) unless `KEEP_CHUNKS=0`.
- `analyze`: y1 holds a population's draws in memory (bulge ~2.4 GB, ~33 min); p6/p7 hold two populations. Give
  `SLURM_OPTS_FINISH` ~12 GB.

## Where things go

```
RUN_ROOT/pipeline_provenance.txt        config, git commit, md5 of the inputs (first sim)
RUN_ROOT/bin/roman.<describe>           the binary the whole run uses
RUN_ROOT/<pop>/plan.txt                 N, chunk size, flags
RUN_ROOT/<pop>/chunks/NNNN/             one chunk: run.log, DONE, outputs
RUN_ROOT/<pop>/test<tag>.dat, h3_pair.dat, run.log, files/MONTLMC/files/*     the merged run (what analysis/ reads)
RUN_ROOT/analysis/{yields,figures}/     results, each script's .log next to its output
pipeline/logs/, pipeline/.stamps/       prep logs and stamps (ignored by git)
```

## Caveats

- **One clone holds one data configuration.** `prep` overwrites `CMD/components/`, `Baseline/*.dat`,
  `config/data_products.h` and the `./roman` binary, which has the row counts compiled in. To work with another
  Besancon file or OpSim baseline, use another clone (or finish and copy out what you need first). Changing only
  `POPULATIONS`, `STRIDE`, `EVENTS`, ... needs no new prep, only a new `RUN_ROOT`: `sim` refuses to reuse a
  `RUN_ROOT` whose plan used other flags.
- **Another OpSim baseline** should start on or before the simulation's day 0, MJD 61141.312 (2026-04-11, pinned in
  `Baseline/readbaselineBulge.py`; Roman's mission is placed on the same clock). `readbaselineBulge.py` keeps the
  visits within reach of a Roman field and silently drops those outside the clock `[0, Tobs]`
  (`Tobs` = 10 years = 3652.4 d in `config/parameters.h`); it prints how many (49 for v5.1.0, which starts 160 d before
  day 0 and ends 161 d before the clock does, so the last 161 days have no Rubin visits). A baseline that starts
  *after* day 0 is not rejected: Rubin simply has no visits at the start of the clock (an event peaking there is
  seen by Roman alone, or by nothing before Roman's mission starts at day 306). `check` prints what it finds.
- Chunking is exact for everything `./roman` appends to once per sightline or event (`test<tag>.dat`,
  `h3_pair.dat`, `MapLMC<tag>.dat`, `LpLMC<tag>.dat`); the `EfLMC` files are the exception (see `merge`). The
  `files/density/` scratch files and any `--dump-samples` directory are per-chunk and not merged.
- Queued Slurm tasks read `pipeline.sh` and the data products when they start: do not edit the script, rerun
  `prep` or change `CMD/`, `Baseline/`, `files/ext/` while jobs are queued or running.
- Tested on Linux (bash 5, GNU coreutils). The helper functions avoid bash-4-only syntax and try BSD `stat`/`md5`,
  but macOS is untested.
