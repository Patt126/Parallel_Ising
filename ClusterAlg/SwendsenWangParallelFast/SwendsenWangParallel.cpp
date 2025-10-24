/**
 * @file SwendsenWangParallel.cpp
 * Implementation of the Swendsen-Wang algorithm in parallel for simulating phase transitions.
 */

#include "SwendsenWangParallel.h"
//#include <cstdlib>
#include <fstream>
#include <cmath>
#include <random>
#include <omp.h> 
#include <vector>
#include <iomanip>
#include <mutex>
std::vector<std::mt19937> rngs;

/**
 * @brief Constructor for SwendsenWangParallel.
 * Initializes the lattice model and simulation parameters.
 * @param interactionStrength Strength of the interaction in the lattice model.
 * @param latticeSize Size of the square lattice.
 * @param T_MIN Minimum temperature for the simulation.
 * @param T_MAX Maximum temperature for the simulation.
 * @param T_STEP Temperature step for the simulation.
 * @param IT Number of iterations for each temperature.
 */


SwendsenWangParallel::SwendsenWangParallel(float interactionStrength, int latticeSize,  float T_MIN, float T_MAX, float T_STEP, long int IT):
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
 * @brief Simulates the phase transition across a range of temperatures.
 * This function iterates over a range of temperatures, performing the Swendsen-Wang algorithm at each step and recording the magnetization.
 */

void SwendsenWangParallel::simulate_phase_transition() {
 
    float m = 0;
    float T = T_MIN;
    // Set up the number of threads for OpenMP
    omp_set_num_threads(omp_get_max_threads());
    int num_threads = omp_get_max_threads();
    rngs.resize(num_threads); // Seed each RNG with a unique seed

    unsigned int seed = static_cast<unsigned int>(time(nullptr));
    for (int i = 0; i < num_threads; ++i) {
        rngs[i].seed(seed + i);
    }
    
    while (T < T_MAX) {
        float P = 1 - exp(-2 * lattice.get_interaction_energy() / T);
        simulate(P, lattice.get_lattice());       
        Temperatures->emplace_back(T);
        T += T_STEP;
        m =  calculate_magnetization_per_site(lattice.get_lattice());
        MagnetizationResults->emplace_back(fabs(m)); 
        lattice.restore_random_lattice();
        //std::cout<<N<<std::endl;
    }
}

/**
 * @brief Calculates the magnetization per site of the lattice.
 * @param lattice A reference to the lattice vector.
 * @return The average magnetization per site.
 */

float SwendsenWangParallel::calculate_magnetization_per_site(std::vector<int>& lattice) {
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

int SwendsenWangParallel::find_set(int x, std::vector<int>& parent) {
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



void SwendsenWangParallel::union_sets(int x, int y, std::vector<int>& parent, std::vector<int>& rank) {
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
 * @brief Performs a single simulation step of the Swendsen-Wang algorithm.
 * @param lattice The lattice to simulate on.
 * @param P Probability threshold for forming bonds.
 * @param parent The parent vector representing the disjoint set forest.
 * @param rank The rank vector used for union by rank.
 */

void SwendsenWangParallel::simulate_step(std::vector<int>& lattice, float P, std::vector<int>& parent, std::vector<int>& rank) {
    

    // Reset parent and rank for each site
    std::fill(parent.begin(), parent.end(), 0);
    std::fill(rank.begin(), rank.end(), 0);
    std::iota(parent.begin(), parent.end(), 0); // Set each site as its own parent initially

    // Form clusters
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        int thread_id = omp_get_thread_num();
        std::mt19937& local_rng = rngs[thread_id];
        std::uniform_real_distribution<float> dist(0.0, 1.0);

        int right = (i % L == L - 1) ? i - (L - 1) : i + 1;
        int down = (i + L) % N;

        if (lattice[i] == lattice[right] && dist(local_rng) < P) {
            union_sets(i, right, parent, rank);
        }
        if (lattice[i] == lattice[down] && dist(local_rng) < P) {
            union_sets(i, down, parent, rank);
        }
    }

    std::vector<int> flip_decision(N);
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        int thread_id = omp_get_thread_num();
        std::mt19937& local_rng = rngs[thread_id];
        std::uniform_int_distribution<int> dist(0, 1);

        flip_decision[i] = dist(local_rng);
    }

    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        if (flip_decision[find_set(i, parent)]) {
            lattice[i] *= -1;
        }
    }
}


/**
 * @brief Runs the full Swendsen-Wang simulation for a given temperature.
 * @param P Probability threshold for forming bonds.
 * @param lattice The lattice to simulate on.
 */

void SwendsenWangParallel::simulate(float P, std::vector<int>& lattice) {
    using namespace std;
    vector<int> parent(N), rank(N, 0);
    for (int i = 1; i < IT; ++i) {
        simulate_step(lattice, P, parent, rank);
        }
}

/**
 * @brief Stores the results of the simulation to a file.
 * Results include the magnetization and temperature for each step of the simulation.
 */

void SwendsenWangParallel::store_results_to_file() const {
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
