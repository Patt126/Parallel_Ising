/**
 * @file Rem.cpp
 * Implementation of the Replica Exchange Monte Carlo (REM) algorithm for simulating phase transitions.
 */
#include "Rem.h"
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <random>
#include <omp.h> 
#include <iomanip>  

std::vector<std::mt19937> rngs;

/**
 * @brief Constructor for the Rem class.
 * Initializes the lattice model and simulation parameters.
 * @param interactionStrength Strength of the interaction in the lattice model.
 * @param latticeSize Size of the square lattice.
 * @param T_MIN Minimum temperature for the simulation.
 * @param T_MAX Maximum temperature for the simulation.
 * @param T_STEP Temperature step for the simulation.
 * @param IT Number of iterations for each temperature.
 */

Rem::Rem(float interactionStrength, int latticeSize,  float T_MIN, float T_MAX, float T_STEP, long int IT):
     lattice(interactionStrength, latticeSize),  
        MagnetizationResults(),
        Temperatures(),
        T_STEP(T_STEP) ,
        T_MIN(T_MIN),
        T_MAX(T_MAX),
        L(latticeSize),
        N(latticeSize * latticeSize),
        IT(IT)
{      

    // Initialize MagnetizationResults with std::make_unique
    MagnetizationResults = std::make_unique<std::vector<float> >();
    // Initialize temperatures with std::make_unique
    Temperatures = std::make_unique<std::vector<float> >();                     

}

/**
 * @brief Simulates the phase transition process across a range of temperatures using the REM method.
 * This function divides the temperature range among MPI processes and performs simulations in parallel.
 * @param argc Argument count from the command line.
 * @param argv Argument vector from the command line.
 */

void Rem::simulate_phase_transition(int argc, char* argv[]) {
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Seed RNG uniquely for each process
    srand(static_cast<unsigned>(time(nullptr)) + world_rank);

    // Calculate the start temperature for each process
    float T_start = T_MIN + world_rank * 0.4f;
    float T_end = T_MAX + world_rank*0.4f;

    omp_set_num_threads(omp_get_max_threads());
    int num_threads = omp_get_max_threads();
    rngs.resize(num_threads);

    // Unique seed for each thread
    std::random_device rd;
    for (int i = 0; i < num_threads; ++i) {
        rngs[i].seed(rd() + i);
    }

    for (float T = T_start; T <= T_end; T += T_STEP) {
        float P = 1 - exp(-2 * lattice.get_interaction_energy() / T);
        simulate(P, lattice.get_lattice(), world_size, T, lattice.J, world_rank);       
        Temperatures->emplace_back(T);

        // Calculate and store magnetization per site
        float m = calculate_magnetization_per_site(lattice.get_lattice());
        MagnetizationResults->emplace_back(fabs(m)); 

        // Restore lattice to a random state for next iteration
        lattice.restore_random_lattice();
        
        MPI_Barrier(MPI_COMM_WORLD);
/*
        if (world_rank == 0) {
            std::cout << "Magnetization at T=" << T << ": " << m << std::endl;

        } */
    } 
}


/**
 * @brief Evaluates the energy of the given lattice configuration.
 * @param lattice A reference to the lattice vector.
 * @param J Interaction strength parameter.
 * @return The calculated energy of the lattice.
 */


//To correctly evaluate energy
float Rem::evaluate(std::vector<int>& lattice, float J) {
    int sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < N; i++) {
        //energy contibute
        if (i >= L) { //NO FIRST ROW
            sum += lattice[i - L] * lattice[i] * 2;
        }
        if (i % L != 0) { //NO FIRST COLUMN
            sum += lattice[i - 1] * lattice[i] * 2;
        }

        if (i >= L * (L - 1)) { //LAST ROW
            sum += lattice[i - L * (L - 1)] * lattice[i] * 2; //times 2 two count also contribute where i = 0
        }
        if ((i + 1) % L == 0) { //LAST COLUMN
            sum += lattice[i - (L - 1)] * lattice[i] * 2;
        }
    }
    return -J * sum;
}

/**
 * @brief Calculates the magnetization per site of the lattice.
 * @param lattice A reference to the lattice vector.
 * @return The average magnetization per site.
 */


float Rem::calculate_magnetization_per_site(std::vector<int>& lattice) {
    int total_spin = std::accumulate(lattice.begin(), lattice.end(), 0);
    return static_cast<float>(total_spin) / N;
}

/**
 * @brief Finds the root of the set that element x belongs to.
 * Uses path compression for efficiency.
 * @param x Element to find the set for.
 * @param parent The parent vector representing the disjoint set forest.
 * @return The root of the set.
 */

int Rem::find_set(int x, std::vector<int>& parent) {
    while (x != parent[x]) {
        parent[x] = parent[parent[x]];  // Path compression
        x = parent[x];
    }
    return x;
}

/**
 * @brief Unions the sets containing elements x and y.
 * Uses rank to keep the tree flat.
 * @param x First element.
 * @param y Second element.
 * @param parent The parent vector representing the disjoint set forest.
 * @param rank The rank vector used for union by rank.
 */


void Rem::union_sets(int x, int y, std::vector<int>& parent, std::vector<int>& rank) {
    x = find_set(x, parent);
    y = find_set(y, parent);
    if (x != y) {
        if (rank[x] < rank[y]) {
            std::swap(x, y);
        }
        parent[y] = x;
        if (rank[x] == rank[y]) {
            rank[x]++;
        }
    }
}

/**
 * @brief Performs a single simulation step of the Monte Carlo algorithm.
 * @param lattice The lattice to simulate on.
 * @param P Probability threshold for forming bonds.
 * @param parent The parent vector representing the disjoint set forest.
 * @param rank The rank vector used for union by rank.
 */

void Rem::simulate_step(std::vector<int>& lattice, float P, std::vector<int>& parent, std::vector<int>& rank) {
    
    // Reset parent and rank for each site
    std::fill(parent.begin(), parent.end(), 0);
    std::fill(rank.begin(), rank.end(), 0);
    std::iota(parent.begin(), parent.end(), 0); // Set each site as its own parent initially

    // Form clusters (Parallelized)
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        int thread_id = omp_get_thread_num();
        std::mt19937& local_rng = rngs[thread_id];
        std::uniform_real_distribution<float> dist(0.0, 1.0);

        int right = (i % L == L - 1) ? i - (L - 1) : i + 1;
        int down = (i + L) % N;

        if (lattice[i] == lattice[right] && dist(local_rng) < P) {
            #pragma omp critical
            union_sets(i, right, parent, rank);
        }
        if (lattice[i] == lattice[down] && dist(local_rng) < P) {
            #pragma omp critical
            union_sets(i, down, parent, rank);
        }
    }

    // Determine flip decisions (Single-threaded)
    std::vector<int> flip_decision(N);
    std::mt19937 rng; // Non-thread local RNG for single-threaded operation
    std::uniform_int_distribution<int> dist(0, 1);
    for (int i = 0; i < N; ++i) {
        flip_decision[i] = dist(rng);
    }

    // Flip clusters (Single-threaded)
    for (int i = 0; i < N; ++i) {
        if (flip_decision[find_set(i, parent)]) {
            lattice[i] *= -1;
        }
    }
}


/**
 * @brief Runs the full Monte Carlo simulation for a given temperature.
 * @param P Probability threshold for forming bonds.
 * @param lattice The lattice to simulate on.
 * @param world_size The size of the MPI world.
 * @param T Current temperature.
 * @param J Interaction strength parameter.
 * @param world_rank The rank of the MPI process.
 */


void Rem::simulate(float P, std::vector<int>& lattice,int world_size, float& T, float J,int world_rank) {
    using namespace std;
    vector<int> parent(N), rank(N, 0);
    for (int i = 1; i < IT; ++i) {
        simulate_step(lattice, P, parent, rank);
                 // Attempt to exchange configuration with neighboring replica every 5 iterations
        if (i % 5 == 0 && world_size > 1) {
            attempt_replica_exchange(lattice, T, J,world_rank, world_size, MPI_COMM_WORLD);
        }
    }

}

/**
 * @brief Attempts to exchange lattice configurations with a neighboring replica.
 * @param lattice The lattice to attempt exchange with.
 * @param T Current temperature.
 * @param J Interaction strength parameter.
 * @param world_rank The rank of the MPI process.
 * @param world_size The size of the MPI world.
 * @param comm The MPI communicator.
 */

void Rem::attempt_replica_exchange(std::vector<int>& lattice, float& T, float J, int world_rank, int world_size, MPI_Comm comm) {
    int partner;
    if (world_rank % 2 == 0) {
        partner = world_rank + 1;
    } else {
        partner = world_rank - 1;
    }

    // Boundary condition for the replicas
    if (partner < 0 || partner >= world_size) return;

    float my_energy = evaluate(lattice,J);
    float partner_energy;
    float my_T = T, partner_T;

    // Exchange temperatures with the partner
    MPI_Sendrecv(&my_T, 1, MPI_FLOAT, partner, 0, &partner_T, 1, MPI_FLOAT, partner, 0, comm, MPI_STATUS_IGNORE);

    // Only exchange energies and attempt swapping if temperatures are different
    if (my_T != partner_T) {
        MPI_Sendrecv(&my_energy, 1, MPI_FLOAT, partner, 1, &partner_energy, 1, MPI_FLOAT, partner, 1, comm, MPI_STATUS_IGNORE);

        // Metropolis criterion
        float delta = (1.0/my_T - 1.0/partner_T) * (partner_energy - my_energy);
        std::mt19937& rng = rngs[0];  // Use the first RNG in the vector
        std::uniform_real_distribution<float> dist(0.0, 1.0);
        if (delta < 0 || exp(-delta) > dist(rng)) {
            // Accept the swap
            //std::cout<<"Exchange Successfull" <<std::endl; Uncomment to see the method in action, more likely to happen in simulations at low temperatire. 
            //It always print it twice for both the process exchanging
            std::vector<int> partner_lattice(N);  // Assuming N is the number of lattice sites
            MPI_Sendrecv(lattice.data(), N, MPI_INT, partner, 2,
                         partner_lattice.data(), N, MPI_INT, partner, 2, comm, MPI_STATUS_IGNORE);
            lattice.swap(partner_lattice);
        }
    }
}

/**
 * @brief Stores the results of the simulation to a file.
 */

void Rem::store_results_to_file() const{
    std::ofstream outFile("result_" + std::to_string(N) + ".txt");

    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
        return;
    }

    outFile << " T  M" << std::endl;

    std::size_t numResults = MagnetizationResults->size();
    for (std::size_t i = 0; i < numResults; ++i) {
        // Set fixed point notation and precision to 2 for temperature values
        outFile << std::fixed << std::setprecision(2) << (*Temperatures)[i] << "   ";
        // Set precision to 4 for magnetization values
        outFile << std::setprecision(4) << (*MagnetizationResults)[i] << std::endl;
    }

    outFile.close();
}


