CXX = g++
CXXFLAGS = -std=c++17 -fopenmp
SRCDIR = ./Lattice ./Metropolis/SerialMetropolis ./Metropolis/DomainDecomposition ./Metropolis/SlidingWindow ./Metropolis
INCLUDEDIRS = -I./Lattice -I./Metropolis/SerialMetropolis -I./Metropolis/DomainDecomposition -I./Metropolis/SlidingWindow  -I./Metropolis

SOURCES = main.cpp \
          $(foreach dir, $(SRCDIR), $(wildcard $(dir)/*.cpp))

OBJECTS = $(SOURCES:.cpp=.o)

TARGET = simulation

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDEDIRS) -c $< -o $@

.PHONY: clean

clean:
	rm -f $(OBJECTS) $(TARGET)

