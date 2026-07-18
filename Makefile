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

# Clean
clean:
	rm -f $(OBJS) $(TARGET)
