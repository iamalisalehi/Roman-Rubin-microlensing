# Compiler and flags
CXX = g++
CXXFLAGS = -g -Wall -Wextra -std=c++17
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
