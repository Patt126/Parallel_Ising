# Use mpicxx as the compiler
CXX = mpicxx
CXXFLAGS = -std=c++17 -fopenmp
TARGET = REMprogram

SRCS = main.cpp ../Lattice/SquareLattice.cpp Rem.cpp
OBJS = $(SRCS:.cpp=.o)

# Main target
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Rule to compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(TARGET) $(OBJS)

# Add a target for running the program with mpirun
run: $(TARGET)
	mpirun ./$(TARGET)

