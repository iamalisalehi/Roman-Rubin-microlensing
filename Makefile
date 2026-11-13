# Compiler and flags
CXX = g++
# -O3 is byte-identical to -O2 (tests/identity_check.sh, 2026-11-13) and about 3% faster; without
# optimisation the Monte Carlo runs ~5x slower. GCC does not reorder floating-point arithmetic without
# -ffast-math, so the results are unchanged. NEVER add -ffast-math, -Ofast, -funsafe-math-optimizations,
# -fassociative-math or -freciprocal-math: they would let the compiler reassociate the chi-squared and
# Fisher sums and silently change the forecast. -flto was measured and gives nothing on top of -O3 (and
# slows the star-list parser), so it is not used. -g keeps CHECK() aborts backtraceable.
CXXFLAGS = -O3 -g -Wall -Wextra -std=c++17 -Iinclude -Iconfig
# CPU-specific flags: an opt-in mechanism, empty by default so the build runs on any x86-64. NOT recommended:
# ARCH_FLAGS="-mavx2 -ffp-contract=off" was measured head to head against plain `make pgo` (5 alternating pairs,
# sightline 795: bulge 0.991 with range 0.986-1.006, bh 0.991 with 0.970-1.003, both ranges crossing 1.0), so it
# gave no gain beyond noise. Check any flag set with tests/identity_check.sh. make does not track flags, so
# `make clean` before changing them.
ARCH_FLAGS ?=
CXXFLAGS += $(ARCH_FLAGS)
# Never let the compiler fuse a multiply and an add (no-op without FMA in the target, as in the default build):
# the VBM observer cache's two copies of the interpolation (vbm_model.o and the VBM object) and every sum then
# stay unfused whatever ARCH_FLAGS adds. Last, so that nothing in ARCH_FLAGS can switch it back on.
CXXFLAGS += -ffp-contract=off
# Profile-guided optimisation, driven by the `pgo` target below; PGO_MODE is internal and empty by default.
PGO_DIR ?= $(CURDIR)/build-pgo
override PGO_DIR := $(abspath $(PGO_DIR))
PGO_GEN_BIN = $(PGO_DIR)/roman-gen
PGO_GEN_FLAGS = -fprofile-generate=$(PGO_DIR)/profile -fprofile-update=single
PGO_USE_FLAGS = -fprofile-use=$(PGO_DIR)/profile -fprofile-partial-training -Wno-missing-profile
ifeq ($(PGO_MODE),gen)
PGO_FLAGS = $(PGO_GEN_FLAGS)
else ifeq ($(PGO_MODE),use)
PGO_FLAGS = $(PGO_USE_FLAGS)
endif
CXXFLAGS += $(PGO_FLAGS)
# The vendored VBMicrolensing header is included as a system header: our -Wall -Wextra apply to our
# code, not to third-party code we do not edit.
CXXFLAGS += -isystem external/VBMicrolensing
# Stamped into every run's provenance block so a result can be traced to its source.
# "unknown" outside a git checkout; a "-dirty" suffix marks edited sources.
GIT_COMMIT := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)$(shell git diff --quiet HEAD 2>/dev/null || echo -dirty)
CXXFLAGS += -DGIT_COMMIT='"$(GIT_COMMIT)"'
# A commit changes no source file, so make would keep the old stamp. .git_stamp is rewritten
# whenever the description differs, and build/run/outputs.o (which prints it) depends on it.
GIT_STAMP := .git_stamp
$(shell echo '$(GIT_COMMIT)' | cmp -s - $(GIT_STAMP) 2>/dev/null || echo '$(GIT_COMMIT)' > $(GIT_STAMP))

LDLIBS = -lgsl -lgslcblas -lm

# GSL from the system by default. GSL_PREFIX=<dir> points at another install; a static copy built
# by `pipeline/pipeline.sh CONFIG setup` in deps/gsl is picked up automatically.
GSL_PREFIX ?= $(wildcard $(CURDIR)/deps/gsl)
ifneq ($(GSL_PREFIX),)
CXXFLAGS += -I$(GSL_PREFIX)/include
LDFLAGS  += -L$(GSL_PREFIX)/lib
endif

# Target executable
TARGET = roman

# Sources: everything under src/, plus the vendored VBMicrolensing library. Objects mirror the
# tree under build/. LIB_OBJS is every object except main(); the test programs link it and supply
# their own main().
SRCS     = $(wildcard src/*.cpp src/*/*.cpp)
VBM_OBJ  = build/external/VBMicrolensing/VBMicrolensingLibrary.o
OBJS     = $(SRCS:src/%.cpp=build/%.o) $(VBM_OBJ)
LIB_OBJS = $(filter-out build/main.o,$(OBJS))

# Default target. The VBMicrolensing object is named first because it is the longest single compile:
# under make -j it then runs alongside the others instead of after them. Link order is unchanged.
all: $(VBM_OBJ) $(TARGET)

# Link step. Must be run from the repo root: every data path is a relative path.
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# -MMD -MP write build/<name>.d header dependencies (included below), so editing a header
# rebuilds exactly the objects that use it.
build/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

build/tests/%.o: tests/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# run/outputs.o prints GIT_COMMIT into run_provenance.txt; rebuild it when the stamp changes.
build/run/outputs.o: $(GIT_STAMP)

# VBMicrolensing (external/VBMicrolensing/, LGPL-3.0, copied from upstream with two local additions, see its README):
# finite-source magnification, astrometric centroids and ephemeris-based
# parallax. Its own flags: the same -O3 and -std=c++17 as ours and never -ffast-math (tests/vbm_test.cpp
# checks that its point-source magnification is bit-identical to ours), but -w: warnings in code we do
# not edit are noise.
VBM_CXXFLAGS = -O3 -g -std=c++17 -w $(ARCH_FLAGS) -ffp-contract=off $(PGO_FLAGS)

$(VBM_OBJ): external/VBMicrolensing/VBMicrolensingLibrary.cpp external/VBMicrolensing/VBMicrolensingLibrary.h
	@mkdir -p $(dir $@)
	$(CXX) $(VBM_CXXFLAGS) -c $< -o $@

# The header is included with -isystem, which -MMD leaves out of the .d files, and it holds inline code
# (LOCAL ADDITION 1 and 2): say here which objects depend on it.
build/events/vbm_model.o build/tests/vbm_test.o build/tests/vbm_cache_test.o build/tests/model_test.o: external/VBMicrolensing/VBMicrolensingLibrary.h

# ---------------------------------------------------------------------------
# Fisher-matrix regression fixture (tests/fisher_fixture.cpp)
#
# Hand-built synthetic events -> real FisherM/ErrorCal -> a diffable table of sigmas.
# Needs no data files:
#     make fishertest && ./fishertest > after.txt && diff before.txt after.txt
# Links every object except main.o, so it exercises the real FisherM, ErrorCal, evaluateModel
# and invert_matrix.
# ---------------------------------------------------------------------------
FIXTURE_TARGET = fishertest

$(FIXTURE_TARGET): build/tests/fisher_fixture.o $(LIB_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

fishertest-run: $(FIXTURE_TARGET)
	./$(FIXTURE_TARGET)

# ---------------------------------------------------------------------------
# CCM89 extinction-law unit test (tests/extinction_test.cpp)
#
# Pins A_V/A_V = 1 at V, a monotone fall from u to F146, and A_lambda/A_V at the seven
# survey bands. Needs no data files.
#     make extinctiontest && ./extinctiontest
# ---------------------------------------------------------------------------
EXT_TARGET = extinctiontest

# Links only galaxy/extinction.o, the one module it tests.
$(EXT_TARGET): build/tests/extinction_test.o build/galaxy/extinction.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# ---------------------------------------------------------------------------
# Astrometric error-model unit test (tests/noise_test.cpp)
#
# Pins errRomanA to McKinnon & van der Marel 2026's tabulated F146 curve and errlsstA to the
# kappa * FWHM / SNR form. Needs only the vendored files/sigma_roman.txt; run from the repo root:
#     make noisetest && ./noisetest
# ---------------------------------------------------------------------------
NOISE_TARGET = noisetest

$(NOISE_TARGET): build/tests/noise_test.o $(LIB_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# ---------------------------------------------------------------------------
# Input-reader unit test (tests/parse_test.cpp)
#
# Pins the buffered from_chars reader (src/util/parse.cpp) and read_cmd against std::istream: every
# token of the four CMD component files, the whole CMD struct against a verbatim copy of the istream
# read_cmd, and synthetic token forms. Reads the CMD star lists (~630 MB, about a minute of CPU; two
# CMD structs are read one after the other, ~700 MB each); run from the repo root:
#     make parsetest && ./parsetest
# ---------------------------------------------------------------------------
PARSE_TARGET = parsetest

$(PARSE_TARGET): build/tests/parse_test.o $(LIB_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# ---------------------------------------------------------------------------
# VBMicrolensing unit test (tests/vbm_test.cpp)
#
# Pins the vendored library against independent references: finite-source magnification against
# direct integrals over the source disc (uniform and limb-darkened), PSPLMag against our own
# point-source formula bit for bit, and the annual-parallax source track against astropy's
# ephemeris. The references are in tests/vbm_reference.h, written by tests/vbm_reference.py.
# Needs only the vendored files; run from the repo root:
#     make vbmtest && ./vbmtest
# ---------------------------------------------------------------------------
VBM_TARGET = vbmtest

# Links only the library it tests.
$(VBM_TARGET): build/tests/vbm_test.o $(VBM_OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# ---------------------------------------------------------------------------
# Test of the VBMicrolensing parallax cache (tests/vbm_cache_test.cpp)
#
# LOCAL ADDITION 2 of the vendored library: the observers' positions computed once per sightline epoch.
# Two library objects set up as the adapter does, one with the cache; every output array and the
# extrapolation flag of the light curves compared with memcmp, over three sightlines' epoch grids built
# from Baseline/BulgeBaseline.dat and RomanBaseline.dat (or a synthetic stand-in), with off-grid times, a stale grid, a table reload.
# Without those two git-ignored inputs it builds a synthetic grid from the committed tables (CI runs it that
# way); run from the repo root:
#     make vbmcachetest && ./vbmcachetest
# ---------------------------------------------------------------------------
VBMCACHE_TARGET = vbmcachetest

$(VBMCACHE_TARGET): build/tests/vbm_cache_test.o $(VBM_OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# ---------------------------------------------------------------------------
# Sky-geometry unit test of the light-curve model (tests/model_test.cpp)
#
# Pins the adapter to VBMicrolensing (src/events/vbm_model.cpp) and the Galactic axes of the kinematics
# against astropy: the ICRS direction of Galactic (l, b) and the rotation to (North, East); that the
# kinematics' projection axes are (e_l, e_b); the Earth's velocity at t0; Roman's L2 ephemeris as a
# satellite; and the light-curve model as the simulator calls it (evaluateModel, evaluateTrack): the side of
# the lens the source passes, the centroid, the heliocentric-to-geocentric conversion, the annual and
# satellite parallax. The references are in tests/model_reference.h, written by tests/model_reference.py.
# Links the real vrel and the adapter. Needs only committed files; run from the repo root:
#     make modeltest && ./modeltest
# ---------------------------------------------------------------------------
MODEL_TARGET = modeltest

$(MODEL_TARGET): build/tests/model_test.o $(LIB_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Verify config/data_products.h against the data files on disk. Not part of `all` (CI has no
# data files). To refresh the header, run `python3 tools/sync_data_products.py`.
check-data:
	python3 tools/sync_data_products.py --check

.PHONY: check-data

# Header dependencies written by -MMD
-include $(OBJS:.o=.d) build/tests/fisher_fixture.d build/tests/extinction_test.d build/tests/noise_test.d build/tests/vbm_test.d build/tests/vbm_cache_test.d build/tests/model_test.d build/tests/parse_test.d

# ---------------------------------------------------------------------------
# Profile-guided production build:  make clean && make -j2 pgo     (the production command)
#
# Always starts clean: builds an instrumented binary at $(PGO_GEN_BIN), never at ./roman, runs
# tools/pgo_train.sh with it (the stub plus bulge and bh sightline 795 at production settings, fixed seeds;
# profile in $(PGO_DIR)/profile, runs in $(PGO_DIR)/train), checks that a profile was written, then rebuilds
# roman and fishertest with it. ./roman does not exist until the very last step, so an interrupted or failed
# `make pgo` leaves none (the pipeline's `make -q roman` then fails, as it should). A trap and the failure
# handler also remove build/ and the instrumented binary; if even that is cut short (power loss), build/ may
# hold instrumented objects: a plain `make roman` then fails loudly at link time (undefined __gcov_* symbols).
# After an interrupted `make pgo`, run `make clean`.
# Output bytes do not change (checked with tests/identity_check.sh). Against the -O3 default the build is
# about x0.955 in CPU time (sightline 795, bulge; 0.923 against 0.966 relative to the -O2 base); outside the
# training sightlines x0.959 (bulge 796) and x0.964 (ns 795).
# ---------------------------------------------------------------------------
pgo:
	rm -rf build $(TARGET) $(PGO_DIR)/profile $(PGO_DIR)/train $(PGO_GEN_BIN)
	mkdir -p $(PGO_DIR)
	trap 'rm -rf build $(PGO_GEN_BIN); exit 130' INT TERM HUP; \
	fail() { rm -rf build $(TARGET) $(PGO_GEN_BIN); exit 1; }; \
	$(MAKE) PGO_MODE=gen TARGET=$(PGO_GEN_BIN) $(PGO_GEN_BIN) || fail; \
	tools/pgo_train.sh $(PGO_GEN_BIN) $(PGO_DIR)/train || fail; \
	[ -n "$$(find $(PGO_DIR)/profile -name '*.gcda' 2>/dev/null | head -n 1)" ] \
	    || { echo "make pgo: training wrote no profile (*.gcda) in $(PGO_DIR)/profile" >&2; fail; }; \
	rm -rf build $(PGO_GEN_BIN); \
	$(MAKE) PGO_MODE=use $(TARGET) $(FIXTURE_TARGET) || fail

.PHONY: pgo

# Clean: only what this Makefile made (PGO_DIR may be set by the user, so the directory itself goes only if empty).
clean:
	rm -rf build $(PGO_DIR)/profile $(PGO_DIR)/train $(PGO_GEN_BIN) $(TARGET) $(FIXTURE_TARGET) $(EXT_TARGET) $(NOISE_TARGET) $(VBM_TARGET) $(VBMCACHE_TARGET) $(MODEL_TARGET) $(PARSE_TARGET)
	-rmdir $(PGO_DIR) 2>/dev/null
