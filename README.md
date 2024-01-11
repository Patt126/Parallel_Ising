## Benchmarking Strategies for Parallel 2D Ising Model Simulation     
This project implements the Ising model simulation using C++ and offers parallel computation using MPI (for distributed systems), OpenMP (for multicore CPUs), and CUDA (for NVIDIA GPUs). It models the magnetic dipole moments of atomic "spins" on a lattice, simulating phase transitions in ferromagnetic materials.

## Features

- 2D Ising Model simulation with periodic boundary conditions.
- Serial computation with C++ with three differend algorithms, Metropolis-Hastings, Wolff and Swendesn-Wang
- (Hybrid) parallel computation utilizing MPI, OpenMP, and CUDA.
- Configurable parameters such as lattice size, temperature range, and number of iterations.
- Analysis of physical properties like total energy and magnetization.

## Prerequisites

- C++ Compiler with C++17 or later support.
- OpenMP for multicore CPU parallelization.
- MPI implementation (e.g., MPICH or OpenMPI) for distributed system parallelization.
- CUDA Toolkit and NVIDIA GPU for GPU-accelerated execution. The Code has been developed and tested with CUDA compute capability higher than 6.1 and CUDA version higher than 11.

## Building the Simulation

**For CPU-based Parallelization with OpenMP:**
```bash
g++ -std=c++17 -fopenmp ising_simulation.cpp -o ising_simulation
```

**For GPU-accelerated Parallelization with CUDA:**
```bash
nvcc -std=c++17 ising_simulation.cu -o ising_simulation
```

**For Distributed System Parallelization with MPI:**
```bash
mpicxx -std=c++17 -fopenmp ising_simulation_mpi.cpp -o ising_simulation_mpi
```
or
```bash
mpic++ -std=c++17 -fopenmp ising_simulation_mpi.cpp -o ising_simulation_mpi
```

## Running the Simulation

**For OpenMP and CUDA Versions:**
```bash
./ising_simulation
```

**For the MPI Version:**
```bash
mpirun -np <number_of_processes> ./ising_simulation_mpi
```
Replace `<number_of_processes>` with the desired number of MPI processes.

## Notes on MPI Version

- The MPI version of the simulation is designed for distributed systems and can be run on multiple nodes or processors.
- Ensure MPI is correctly installed and configured on your system or computing cluster.
- The MPI version can be combined with OpenMP to leverage multi-threading on each node in addition to distributed computing. 

## Notes on CUDA Folder

The CUDA folder contains three different versions. The fastest among them is named `ising_simulation`. There's another implementation that attempts to utilize shared memory, but it's not the quickest.

The output of the program presents the results of the simulation, including the system's energy and magnetization at each temperature.

## Hands On

In the Hands On section, the same Metropolis algorithm used for the 2D Ising Simulation is employed to write parallel code. This code calculates Monte Carlo Integrals on rectangular and spherical domains in any number of dimensions.

