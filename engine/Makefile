# Portable Makefile for the alloy chemical-ordering Monte Carlo engine.
#
#   make            -> build ./alloy_mc
#   make CXX=mpicxx -> override the MPI C++ compiler wrapper
#   make debug      -> build with sanitizers and assertions
#   make clean
#
# On an HPC system, `module load` your MPI stack first so that mpic++/mpicxx
# resolves to the system MPI (see references/hpc.md).

# Default to the MPI wrapper, but honor an explicit env/command-line CXX.
# (`?=` alone is defeated by make's built-in CXX, so test its origin.)
ifeq ($(origin CXX),default)
CXX := mpic++
endif
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra
INCLUDES  = -Iinclude
LDFLAGS  ?=

SRC  = src/config.cpp src/lattice.cpp src/alloy_model.cpp \
       src/parallel_tempering.cpp src/main.cpp
OBJ  = $(SRC:.cpp=.o)
BIN  = alloy_mc

# Serial test binary (no MPI): links everything except main/parallel_tempering.
TEST_BIN = test_consistency
TEST_OBJ = src/config.o src/lattice.o src/alloy_model.o

.PHONY: all debug test clean
all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) -o $@

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

debug: CXXFLAGS = -O0 -g -std=c++17 -Wall -Wextra -fsanitize=address,undefined
debug: clean $(BIN)

# Build the serial test with a plain C++ compiler (no MPI needed).
$(TEST_BIN): tests/test_consistency.cpp $(TEST_OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDES) tests/test_consistency.cpp $(TEST_OBJ) -o $@

test: $(TEST_BIN)
	@echo "Run from a directory containing coupling.input and a control file."

clean:
	rm -f $(OBJ) $(BIN) $(TEST_BIN)
