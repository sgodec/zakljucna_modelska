# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -flto -Wall -Wextra -fopenmp 
LDFLAGS = -fopenmp -llemon -lboost_math_c99

# Directories
LEMON_INCLUDE = /usr/local/include
LEMON_LIB = /usr/local/lib
INCLUDE_DIR = include
SRC_DIR = src

SOURCES = $(SRC_DIR)/main.cpp $(SRC_DIR)/solve_gridgraph.cpp 
OBJECTS = $(SOURCES:.cpp=.o)
EXEC = solve_gridgraph

all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ -L$(LEMON_LIB) $(LDFLAGS)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(LEMON_INCLUDE) -I$(INCLUDE_DIR) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXEC)

$(SRC_DIR)/main.o: $(INCLUDE_DIR)/solve_gridgraph.h 
$(SRC_DIR)/solve_gridgraph.o: $(INCLUDE_DIR)/solve_gridgraph.h

.PHONY: all clean
