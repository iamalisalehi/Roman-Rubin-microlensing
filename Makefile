# Compiler and flags
CXX = g++
# -O2 is not optional: without it the Monte Carlo runs ~5x slower (measured: 257 vs
# 1280 events per 300 s on the same seeded run). It is safe to turn on because GCC
# does NOT reorder floating-point arithmetic without -ffast-math, so results are
# unchanged -- verified bit-identical against an -O0 build across EfLMC2.dat,
# EfLMC2B.dat, LpLMC2.dat, test2.dat, stdout and the full fishertest table.
# -O3 was measured too and gives no further gain. Do NOT add -ffast-math: it would
# let the compiler reassociate the chi-squared and Fisher sums and silently change
# the forecast. -g is kept so the CHECK() aborts still produce a usable backtrace.
CXXFLAGS = -O2 -g -Wall -Wextra -std=c++17
LDLIBS = -lgsl -lgslcblas -lm

# Target executable
TARGET = roman

# Source files
SRCS = Bulge_LSST.cpp Lensing.cpp helper.cpp
OBJS = $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Link step
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

# Compile each .cpp into .o
%.o: %.cpp Bulge.h
	$(CXX) $(CXXFLAGS) -c $<

# ---------------------------------------------------------------------------
# Fisher-matrix regression fixture (tests/fisher_fixture.cpp)
#
# Hand-built synthetic events -> real FisherM/ErrorCal -> a diffable table of sigmas.
# Runs in seconds and needs NO data files, unlike ./roman. Use it to check that a change
# to FisherM moves sigma the way you expect:
#     make fishertest && ./fishertest > after.txt && diff before.txt after.txt
#
# Bulge_LSST.cpp is recompiled here with -DFISHER_FIXTURE_BUILD, which drops its main()
# so the fixture can supply its own while still linking against FisherM, ErrorCal,
# lightcurve and invert_matrix.
# ---------------------------------------------------------------------------
FIXTURE_TARGET = fishertest

$(FIXTURE_TARGET): tests/fisher_fixture.cpp Bulge_LSST.cpp Lensing.cpp helper.cpp Bulge.h
	$(CXX) $(CXXFLAGS) -DFISHER_FIXTURE_BUILD -I. -o $@ \
	    tests/fisher_fixture.cpp Bulge_LSST.cpp Lensing.cpp helper.cpp $(LDLIBS)

fishertest-run: $(FIXTURE_TARGET)
	./$(FIXTURE_TARGET)

# Clean
clean:
	rm -f $(OBJS) $(TARGET) $(FIXTURE_TARGET)
