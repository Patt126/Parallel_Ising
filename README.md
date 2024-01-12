# Benchmarking Strategies for Parallel 2D Ising Model Simulation     
This project implements the Ising model simulation using C++ and offers parallel computation using  OpenMP (for multicore CPUs),. It models the magnetic dipole moments of atomic "spins" on a lattice, simulating phase transitions in ferromagnetic materials.

## Features

- 2D Ising Model simulation with periodic boundary conditions.
- Serial computation with C++ with  Metropolis-Hastings alghoritm.
-  parallel computation utilizing OpenMP with different parallelization techniques.
- Configurable parameters such as lattice size, temperature range, and number of iterations.
- Analysis of physical properties like total energy, magnetization, correlation lenght and specific heat. Following their evolution with temperature around the critical one.
- Performance evaluation for benchmarking with fixed number of iteration.
- Convergence criterion implementation: exact solution comparison and fluctuation analisys.
- Visualization

# Parallel Algorithm

## Domain Decomposition
The code implements a parallel Metropolis algorithm in which lattice is divided into blocks,. Unique to this approach is the use of atomic updates for boundary sites, ensuring thread safety and consistency during parallelized simulations.

## Sliding Window
Lattice is divided into blocks, and it is implemented a Sliding Window Metropolis algorithm. In this approach, spins are updated only if they are not at block blundary. A periodic translation of the lattice matrix is performed, ensuring proper visiting and updates of boundary spins.

## Checkboard
The lattice is divided into two alternating sublattices resembling a checkerboard, and the simulation progresses by updating spins on each sublattice in alternating steps. 

## Prerequisites

- C++ Compiler with C++17 or later support.
- OpenMP for multicore CPU parallelization.
- Gnuplot-iostream Library for results visualisation. You can find the Gnuplot-iostream library [Here](http://stahlke.org/dan/gnuplot-iostream/). Follow the installation instructions provided on the website.



## Building and Running the Simulation
Just move to the folder corresponding to a simulation method (autoconvergence or there are several inside metropolis).
```bash
cd path/to/folder
make
./simulation
```
Simulation start, at the end file relative to performance and results are inside the same directory.
```bash
make clean 
```
To remove the executable.



## Visualize Results

```bash
cd Visualisation
g++ -std=c++17 -o visualisation visualisation.cpp
./visualisation
```
After entering the file path, the program will display an interactive menu:

-Option 1: Plot Energy
-Option 2: Plot Magnetization over Temperature
-Option 3: Plot Specific Heat
-Option 4: Evaluate Critical Temperature
-Option 0: Exit
Choose an option by entering the corresponding number. The program will execute the selected operation.





