CXX = clang++
CXXFLAGS = -std=c++17 -fopenmp
SRCDIR = ./Lattice ./Metropolis/AutostopSerial ./Metropolis/AutostopDomainDecomposition 
INCLUDEDIRS = -I./Lattice -I./Metropolis/AutostopSerial -I./Metropolis/AutostopDomainDecomposition 

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