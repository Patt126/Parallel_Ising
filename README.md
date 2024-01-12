## Benchmarking Strategies for Parallel 2D Ising Model Simulation     
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

## Parallel alghoritm
- DOMAIN DECOMPOSITION: the code implements a parallel Metropolis algorithm in which lattice is divided into blocks,. Unique to this approach is the use of atomic updates for boundary sites, ensuring thread safety and consistency during parallelized simulations.
- SLIDING WINDOW: lattice is divided into blocks, and it is implemented a Sliding Window Metropolis algorithm. In this approach, spins are updated only if they are not at block blundary. A periodic translation of the lattice matrix is performed, ensuring proper visiting and updates of boundary spins.

## Prerequisites

- C++ Compiler with C++17 or later support.
- OpenMP for multicore CPU parallelization.


## Building the Simulation



## Running the Simulation




