CXX = mpic++
CXXFLAGS = -std=c++17 -fopenmp
MPIFLAGS = -lmpi
SRCDIR = ./Lattice  ./Metropolis/DomainDecomposition ./Metropolis
INCLUDEDIRS = -I./Lattice  -I./Metropolis/DomainDecomposition   -I./Metropolis

SOURCES = main.cpp \
          $(foreach dir, $(SRCDIR), $(wildcard $(dir)/*.cpp))

OBJECTS = $(SOURCES:.cpp=.o)

TARGET = simulation

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(MPIFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDEDIRS) -c $< -o $@

.PHONY: clean

clean:
	rm -f $(OBJECTS) $(TARGET)
