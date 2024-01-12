CXX = clang++
CXXFLAGS = -std=c++17 
SRC = main.cpp SquareLattice.cpp AutostopMontecarlo.cpp
TARGET = Autoprogram

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(TARGET)

