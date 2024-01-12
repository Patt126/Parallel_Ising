CXX = g++
CXXFLAGS = -std=c++17 -fopenmp
SRC_DIR = .
LATTICE_DIR = ../../Lattice
SRC = $(SRC_DIR)/main.cpp $(SRC_DIR)/SlidingWindow.cpp $(LATTICE_DIR)/SquareLattice.cpp
TARGET = simulation

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)

