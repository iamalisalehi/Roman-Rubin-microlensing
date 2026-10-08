# shellcheck shell=bash disable=SC2034
# config.example.sh -- settings for pipeline/pipeline.sh.   Copy it, edit the copy, then:
#
#     cp pipeline/config.example.sh my_run.sh
#     pipeline/pipeline.sh my_run.sh check
#
# This is a bash file that is `source`d from the repository root, so relative paths are relative to the
# repo root and a value may use $HOME or $(...). A deleted line falls back to the default shown.
# Do not put whitespace in any path.

# ---- step 1-2: the two downloads ------------------------------------------------------------------
# Besancon simulation (model 1612, ONE field, Av = 0, no magnitude limit: see pipeline/README.md).
# Not downloadable by script: it comes from the Besancon web form.
BESANCON_CATALOGUE=CMD/Besancon/bos10
# Rubin OpSim baseline database (sqlite). It must start on or before MJD 61141.312 (the simulation's day 0).
OPSIM_DB=Baseline/baseline_v5.1.0_10yrs.db
# `fetch` downloads OPSIM_DB from here if the file is missing (0.8 GB). Empty = never download.
OPSIM_URL=https://s3df.slac.stanford.edu/data/rubin/sim-data/sims_featureScheduler_runs5.1/baseline/baseline_v5.1.0_10yrs.db
# MIST v2 bolometric-correction tables (LSST.txz, Roman.txz) are fetched from this directory URL.
MIST_URL_BASE=https://mist.science/BC_tables/v2

# ---- extinction ------------------------------------------------------------------------------------
# MUST BE SET. Path of a prebuilt ext_tables.dat to copy into files/ext/ (with the ext_provenance.json
# next to it), or the word `build` to compute it with maps.py (downloads the 7 GB DECaPS map, ~35 min
# of queries). The table depends only on the survey geometry, not on the Besancon file or the OpSim db,
# so copying one is the normal choice.
EXT_TABLES=

# ---- step 3: what to simulate ----------------------------------------------------------------------
# Star-list builder (CMD/BolometricCorrection.py). Default: Besancon's own dwarfs, no synthetic dwarfs
# added. CATALOGUE_DWARFS: besancon | empirical.  CATALOGUE_FILL: none | kroupa | koshimoto.
CATALOGUE_DWARFS=besancon
CATALOGUE_FILL=none
# Day on the simulation clock (0 = 2026-04-11) at which Roman's season 0 starts.
ROMAN_MISSION_START=306

# Where this run's chunks, merged tables, logs and analysis go. One RUN_ROOT per set of options below.
RUN_ROOT=runs/pipeline_run1
# Lens populations, one ./roman run each: bulge bh ns besancon macho-* (`./roman --help` lists them).
POPULATIONS="bulge bh ns"
# Sightline grid: STRIDE x 0.02 deg outside Roman's footprint, STRIDE_ROMAN x 0.02 deg inside it.
# STRIDE_ROMAN must divide STRIDE. Larger = fewer sightlines = faster and coarser.
STRIDE=10
STRIDE_ROMAN=5
# Per-sightline targets: detected events (EVENTS) and Fisher-characterised events (LENSES).
EVENTS=300
LENSES=50
# Per-sightline cap on drawn stars. Empty = the binary's default (5e4).
MAXDRAWS=
# 1 = characterise every detection twice (Roman at L2 and at Earth) -> h3_pair.dat, the satellite-parallax figure.
PAIR_SATELLITE=1
SEED=42
# 1 = scan only the 0.1x0.1 deg stub patch (a smoke test; minutes instead of days). Needs STRIDE=1 and
# STRIDE_ROMAN=1 (36 sightlines); at the production strides the patch has no sightline in some Roman
# coverage class and ./roman refuses. Smoke-test values: EVENTS=10 LENSES=2 MAXDRAWS=3000 CHUNK_SIZE=10.
STUB=0
# Any further ./roman flags, passed to every population (e.g. "--dchi-det 300").
EXTRA_FLAGS=""

# ---- step 5: how the simulation is split and where it runs ----------------------------------------
# Sightlines per chunk. Chunks are independent runs whose concatenation equals the whole run exactly.
# 20 gives ~100 chunks per population at the default grid (2013 sightlines): a footprint sightline costs
# 4-6 CPU-min, so a chunk takes ~10 min on average and at most ~2 h, and the ~25 s of input reading per
# chunk is a few percent. Use a smaller value to spread work over more workers (cost is uneven over the
# index range, so have more chunks than workers).
CHUNK_SIZE=20
# local | slurm
SCHEDULER=local
# local: number of ./roman processes at once. Each needs ~0.7 GB of RAM (and one core).
JOBS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)
# slurm: words that go into #SBATCH lines of the chunk array jobs (one chunk per task, 1 core).
# ~2 GB and 2-3 h are comfortable at the default flags. No spaces inside an option.
SLURM_OPTS="--account=CHANGE_ME --partition=CHANGE_ME --time=03:00:00 --mem=2G --cpus-per-task=1"
# slurm: the same for the final merge + analysis job (holds a whole population's table in memory:
# the bulge y1 analysis peaks at ~2.5 GB; give it room and an hour or more). Default: SLURM_OPTS.
SLURM_OPTS_FINISH="--account=CHANGE_ME --partition=CHANGE_ME --time=04:00:00 --mem=12G --cpus-per-task=1"
# slurm: 1 = submit now, 0 = only write RUN_ROOT/slurm/*.sbatch and print how to submit them.
SUBMIT=1
# 1 = keep RUN_ROOT/<pop>/chunks after merging (the merged tables then exist twice on disk); 0 = delete them.
KEEP_CHUNKS=1

# ---- step 6: analysis ------------------------------------------------------------------------------
# yields  = absolute yields per population (y1) and yield vs abundance F (y3)
# figures = per population f1 f2 f3 f4 h5 (h3), then p6 and p7 across the populations
ANALYSES="yields figures"

# ---- Python ---------------------------------------------------------------------------------------
# The interpreter every script runs under. `setup` creates the venv if this path ends in /bin/python.
PYTHON=.roman/bin/python

# ---- C++ library ----------------------------------------------------------------------------------
# GSL (GNU Scientific Library), the one non-standard library ./roman links.
#   auto   = use the system's GSL; if its headers are missing, `setup` builds a static copy in deps/gsl
#   system = system GSL only (check fails if it is missing)
#   build  = always use the copy in deps/gsl (built once by `setup`, ~5 min, no root needed)
GSL=auto
