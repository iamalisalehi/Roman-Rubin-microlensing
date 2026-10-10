# Compiler and flags
CXX = g++
# -O2 is required: without it the Monte Carlo runs ~5x slower. Results are unchanged
# because GCC does not reorder floating-point arithmetic without -ffast-math; -O3 gives no
# further gain. Do NOT add -ffast-math: it would let the compiler reassociate the chi-squared
# and Fisher sums and silently change the forecast. -g keeps CHECK() aborts backtraceable.
CXXFLAGS = -O2 -g -Wall -Wextra -std=c++17 -Iinclude -Iconfig
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

# Sources: everything under src/. Objects mirror that tree under build/. LIB_OBJS is every
# object except main(); the test programs link it and supply their own main().
SRCS     = $(wildcard src/*.cpp src/*/*.cpp)
OBJS     = $(SRCS:src/%.cpp=build/%.o)
LIB_OBJS = $(filter-out build/main.o,$(OBJS))

# Default target
all: $(TARGET)

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

# Verify config/data_products.h against the data files on disk. Not part of `all` (CI has no
# data files). To refresh the header, run `python3 tools/sync_data_products.py`.
check-data:
	python3 tools/sync_data_products.py --check

.PHONY: check-data

# Header dependencies written by -MMD
-include $(OBJS:.o=.d) build/tests/fisher_fixture.d build/tests/extinction_test.d build/tests/noise_test.d

# Clean
clean:
	rm -rf build $(TARGET) $(FIXTURE_TARGET) $(EXT_TARGET) $(NOISE_TARGET)
