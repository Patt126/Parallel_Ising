# Parallel 2D Ising Model Simulation and Benchmarking

This project implements several parallel strategies for simulating the 2D Ising model using **C++**, integrating **MPI**, **OpenMP**, and **CUDA** backends. It models the magnetic dipole moments of atomic spins on a lattice and reproduces the phase transition of ferromagnetic materials. The repository unifies the work carried out **together** by me and colleague Luigi Pagani who explored complementary parallelization approaches.

---

## Features

* 2D Ising Model with **periodic boundary conditions**.
* Multiple simulation algorithms:

  * **Metropolis–Hastings** (serial, OpenMP,MPI, CUDA).
  * **Wolff** and **Swendsen–Wang** (cluster updates, OpenMP).
  * **Replica Exchange Method (REM)** combining **MPI + OpenMP**.
* Parallelization schemes:

  * **Domain decomposition**.
  * **Checkerboard** (red–black) updates.
  * **Sliding window** parallel updates.
* Benchmarking framework to measure scaling and convergence.
* Configurable parameters: lattice size, temperature range, iteration count.
* Visualization utilities (gnuplot-iostream).
* Physical observables: total energy, magnetization, specific heat, correlation length.

---

## Algorithms Overview

### 1. Metropolis–Hastings

Local spin-flip dynamics with Boltzmann acceptance. Implemented in serial, OpenMP, and CUDA variants. Parallelization via checkerboard and domain decomposition.

### 2. Cluster Algorithms

* **Wolff Algorithm**: builds and flips single clusters.
* **Swendsen–Wang Algorithm**: identifies all clusters per step using parallel labeling.
  Both reduce critical slowing down and are OpenMP-parallelized.

### 3. Replica Exchange (MPI + OpenMP Hybrid)

Each MPI process simulates a replica at a distinct temperature; periodic replica exchanges improve sampling. Each replica internally uses OpenMP-parallelized Swendsen–Wang dynamics.

---

## Prerequisites

* **C++17** or later compiler
* **OpenMP** (for CPU parallelization)
* **MPI** (e.g., MPICH, OpenMPI)
* **CUDA 12.1+** (for GPU version)
* **gnuplot-iostream** library for visualization

---

## Build and Run

### Serial / OpenMP / MPI

```bash
make
./simulation
make clean
```

### MPI version

```bash
cd MPI
make
mpirun -np <num_processes> ./simulation
```

Example:
`L_TOTAL=256`, `num_processes=4`, `NUMTHREAD=16`.

### CUDA version

```bash
cd Test/MetropolisCuda
nvcc -std=c++17 MetropolisCuda.cu -o MCprogram
./MCprogram
```

### Cluster Algorithms

```bash
cd src/Wolff
make && ./Wprogram && make clean
cd src/SwendsenWangParallel
make && ./SWprogram && make clean
```

---

## Visualization

```bash
cd Visualisation
g++ -std=c++17 -o visualisation visualisation.cpp
./visualisation
```

Menu options allow plotting energy, magnetization, specific heat, and critical temperature.

---

## Benchmarking

Fixed-iteration benchmarks based on ( L^{4.25} ) scaling to evaluate speed-up across methods and thread counts.
Deterministic step counts ensure consistent performance comparisons.

---

## Results

### Physical Validation

All implementations reproduce the **2D Ising model phase transition** at the expected critical temperature ( T_c \approx 2.27,K ). Magnetization follows Onsager’s analytical curve:

* Ordered ferromagnetic phase for ( T < T_c ).
* Disordered paramagnetic phase for ( T > T_c ).
* Sharp magnetization drop near ( T_c ).

### Performance and Scaling

#### Domain Decomposition (OpenMP)

* Most balanced and stable parallel scheme.
* Achieved **speed-up up to 2.8×** on 64 threads.
* Scales consistently with thread count, limited only by atomic operations on boundaries.
* Accuracy comparable to the serial Metropolis implementation.

#### Sliding Window and Checkerboard (CPU)

* Less efficient due to synchronization and boundary update overhead.
* Suitable mainly for conceptual benchmarking.

#### Checkerboard (CUDA)

* Highly efficient on GPU; fine-grained parallelism eliminates CPU synchronization bottlenecks.
* Demonstrates major throughput gains per lattice site, constrained only by GPU memory.

#### Cluster Algorithms (Wolff, Swendsen–Wang)

* Reduce critical slowing down near ( T_c ).
* Swendsen–Wang parallelized version achieved **speed-up ≈1.5×** with OpenMP before reaching synchronization limits.
* Wolff algorithm remains more accurate but sequential due to cluster dependency.

#### Replica Exchange (MPI + OpenMP)

* Allows independent simulation replicas at different temperatures with periodic exchange.
* Strong scaling across nodes for large systems; efficiency drops for small lattices due to communication latency.

#### GPU vs CPU Comparison (N = 256×256)

| Algorithm                     | Speed-Up / Efficiency | Observation                            |
| ----------------------------- | --------------------- | -------------------------------------- |
| Domain Decomposition          | 2.8×                  | Optimal CPU scaling.                   |
| Sliding Window                | <1×                   | Synchronization overhead.              |
| Checkerboard (CPU)            | <1×                   | Thread contention.                     |
| Checkerboard (CUDA)           | >3×                   | Highest throughput on GPU.             |
| Wolff                         | Sequential            | Best physical convergence.             |
| Swendsen–Wang (OpenMP)        | 1.5×                  | Moderate improvement.                  |
| Replica Exchange (MPI+OpenMP) | Scalable              | Communication-bound on small lattices. |

### Summary

* **Most efficient CPU method:** Domain Decomposition (OpenMP).
* **Most physically accurate:** Cluster algorithms (Wolff, Swendsen–Wang).
* **Best large-scale scaling:** Replica Exchange (MPI+OpenMP).
* **Highest raw performance:** Checkerboard CUDA Metropolis.

All methods collectively confirm the same physical results, validating correctness and scalability across architectures.

---

## Notes

* Use even powers of two for lattice size and thread count to maintain square domain coverage.
* For OpenMP domain decomposition, ensure `sqrt(NUMTHREAD)` divides lattice width.
* For MPI + OpenMP hybrid runs, adapt `--map-by` options according to hardware topology.
