# Compiler and flags
CXX = g++
# -O2 is required: without it the Monte Carlo runs ~5x slower. Results are unchanged
# because GCC does not reorder floating-point arithmetic without -ffast-math; -O3 gives no
# further gain. Do NOT add -ffast-math: it would let the compiler reassociate the chi-squared
# and Fisher sums and silently change the forecast. -g keeps CHECK() aborts backtraceable.
CXXFLAGS = -O2 -g -Wall -Wextra -std=c++17 -Iinclude -Iconfig
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

# VBMicrolensing (external/VBMicrolensing/, LGPL-3.0, copied from upstream with one local addition to the
# header, see its README): finite-source magnification, astrometric centroids and ephemeris-based
# parallax. Its own flags: the same -O2 and -std=c++17 as ours and never -ffast-math (tests/vbm_test.cpp
# checks that its point-source magnification is bit-identical to ours), but -w: warnings in code we do
# not edit are noise.
VBM_CXXFLAGS = -O2 -g -std=c++17 -w

$(VBM_OBJ): external/VBMicrolensing/VBMicrolensingLibrary.cpp external/VBMicrolensing/VBMicrolensingLibrary.h
	@mkdir -p $(dir $@)
	$(CXX) $(VBM_CXXFLAGS) -c $< -o $@

# ---------------------------------------------------------------------------
# Fisher-matrix regression fixture (tests/fisher_fixture.cpp)
#
# Hand-built synthetic events -> real FisherM/ErrorCal -> a diffable table of sigmas.
# Needs no data files:
#     make fishertest && ./fishertest > after.txt && diff before.txt after.txt
# Links every object except main.o, so it exercises the real FisherM, ErrorCal, lightcurve
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
# Sky-geometry unit test of the light-curve model (tests/model_test.cpp)
#
# Pins the adapter to VBMicrolensing (src/events/vbm_model.cpp) and the Galactic axes of the kinematics
# against astropy: the ICRS direction of Galactic (l, b) and the rotation to (North, East); that the
# kinematics' projection axes are (e_l, e_b); the Earth's velocity at t0; and Roman's L2 ephemeris as a
# satellite. The references are in tests/model_reference.h, written by tests/model_reference.py. Links the
# real vrel and the adapter. Needs only committed files; run from the repo root:
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
-include $(OBJS:.o=.d) build/tests/fisher_fixture.d build/tests/extinction_test.d build/tests/noise_test.d build/tests/vbm_test.d build/tests/model_test.d

# Clean
clean:
	rm -rf build $(TARGET) $(FIXTURE_TARGET) $(EXT_TARGET) $(NOISE_TARGET) $(VBM_TARGET) $(MODEL_TARGET)
