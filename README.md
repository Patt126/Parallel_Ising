# Benchmarking Strategies for Parallel 2D Ising Model Simulation     
This project implements the Ising model simulation using C++ and offers parallel computation using  OpenMP (for multicore CPUs),. It models the magnetic dipole moments of atomic "spins" on a lattice, simulating phase transitions in ferromagnetic materials.

## Features

- 2D Ising Model simulation with periodic boundary conditions.
- Serial computation with C++ with  Metropolis-Hastings alghoritm.
-  parallel computation utilizing OpenMP with different parallelization techniques.
- Configurable parameters such as lattice size, temperature range, and number of iterations.
- Automated simulation specifying minimum and maximum lattice size. Allows to get performance and Physical results, stored in dedicated folders 
- Analysis of physical properties like total energy, magnetization, correlation lenght and specific heat. Following their evolution with temperature around the critical one.
- Performance evaluation for benchmarking with fixed number of iteration.
- Convergence criterion implementation: exact solution comparison and fluctuation analisys.
- Visualization

## Update

- In the latest simulation, conducted on a high-performance computer equipped with a substantial number of cores, we were able to explore larger lattices. However, we observed that the approach of pre-storing and manipulating a fixed sequence of random numbers proved to be resource-intensive in terms of memory usage. Consequently, I reverted to a simpler strategy where a random number is generated on-the-fly when needed. While this incurs a higher computational cost, it significantly improves memory performance. In the code, I have ve preserved commented or unused sections from the previous implementation as a form of documentation to trace the development process.
- Autostop code has still some bug about storing results

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
The code is provided with a main program to execute a BenchMark of three of the four methods.
To Build and run in the main directory:  
```bash
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

1. Plot Energy
2. Plot Magnetization over Temperature
3. Plot Specific Heat
4. Evaluate Critical Temperature
5. Exit
   
Choose an option by entering the corresponding number. The program will execute the selected operation.


## Note

- The codes inside Metropolis directory are meant for benchmarking, for that reason they consider, for a given lattice widht $L$, a fixed number of iteration obtained as $L^{4.25}$. This value is obtained from several trial with autostop algorithm, and provide an upper bound to the tipical number of convergence steps. For that reason, this implementation are much slower to autostopping one but are usefull to perform a deterministic simulation in terms of steps, in order to evaluate speed up. 
- In parallel implementation, in order to optimize cuncurrent execution, the lattice is divided in a number of blocks corresponding to the number of thread. With the exclusion of checkboard in which ,considering the alternation of the two sublattice, the number of blocks is twice the number of threads. Be carefull when choosing the number of thread, remembering that it must be a perfect square in order to cover with square blocks a square lattice.
**N.B.** for that reason I warmly suggest to work with a even power of two in term of lattice width and Thread Number.

