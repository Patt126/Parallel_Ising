#include "Wolff.h"
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <random>
#include <unordered_set>
#include <queue>
#include <vector>
#include <iomanip>



/**
 * @file Wolff.cpp
 * @brief Implementation of the Wolff algorithm for simulating phase transitions in a spin lattice.
 */

/**
 * @brief Constructor for the Wolff simulation.
 * @param interactionStrength The interaction strength J in the model.
 * @param latticeSize The size of one side of the square lattice.
 * @param T_MIN The minimum temperature to simulate.
 * @param T_MAX The maximum temperature to simulate.
 * @param T_STEP The step size in temperature for the simulation.
 * @param IT The number of iterations for the simulation.
 */
Wolff::Wolff(float interactionStrength, int latticeSize,  float T_MIN, float T_MAX, float T_STEP, long int IT):
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
 * @brief Runs the simulation to study the phase transition over a range of temperatures.
 */

void Wolff::simulate_phase_transition() {
    using namespace std;
    unsigned seed = time(0);
    srand(seed);
    float M = 0;
    float  energy =0;
    float T = T_MIN;
    using namespace std;
    std::vector<int> rand_vect;
    create_rand_vect(rand_vect);
    while (T < T_MAX) {
        simulate(T, lattice.get_lattice(),energy,  M,rand_vect,lattice.J);       
        Temperatures->emplace_back(T);
        T += T_STEP;
        M =  calculate_magnetization_per_site(lattice.get_lattice());
        MagnetizationResults->emplace_back(fabs(M)); 
        lattice.restore_random_lattice();
        //std::cout<<M<<std::endl;
    }

}

/**
 * @brief Calculates the magnetization per site of the lattice.
 * @param lattice A reference to the vector representing the spin lattice.
 * @return The average magnetization per site.
 */

float Wolff::calculate_magnetization_per_site(std::vector<int>& lattice) {
    int total_spin = std::accumulate(lattice.begin(), lattice.end(), 0);
    return static_cast<float>(total_spin) / N;
}

/**
 * @brief Adds a spin to the cluster if it meets the criteria.
 * @param lattice The spin lattice.
 * @param cluster A set containing indexes of spins already in the cluster.
 * @param spin_queue A queue of spins to be processed.
 * @param P The probability threshold for adding a neighbor spin to the cluster.
 * @param i The index of the current spin.
 */




void Wolff::add_to_cluster(std::vector<int>& lattice, std::unordered_set<int>& cluster, std::queue<int>& spin_queue, float P, int i) {
    if (cluster.find(i) == cluster.end()) {
        cluster.insert(i);
        std::vector<int> neighbors = {
            (i + L) % N,
            (i - L + N) % N,
            (i + 1) % L + (i / L) * L,
            (i - 1 + L) % L + (i / L) * L
        };
        for (int neighbor : neighbors) {
            if (lattice[neighbor] == lattice[i] && (rand() / (float)RAND_MAX) < P) {
                if (cluster.find(neighbor) == cluster.end()) {
                    spin_queue.push(neighbor);
                }
            }
        }
    }
}

/**
 * @brief Updates the lattice for a given temperature T using the Wolff algorithm.
 * @param lattice The spin lattice.
 * @param T The current temperature.
 * @param J The interaction strength.
 */

void Wolff::update(std::vector<int>& lattice, float T,  float J) {
    float P = 1 - exp(-2 * J / T);
    int seed = rand() % N;
    std::unordered_set<int> cluster;
    std::queue<int> spin_queue;
    spin_queue.push(seed);
    while (!spin_queue.empty()) {
        int current = spin_queue.front();
        spin_queue.pop();
        add_to_cluster(lattice, cluster, spin_queue, P, current);
    }
    for (int i : cluster) {
        lattice[i] *= -1;
    }
}
/**
 * @brief Simulates the lattice for a single temperature.
 * @param T The temperature at which to simulate.
 * @param lattice The spin lattice.
 * @param energy A reference to the energy variable.
 * @param M A reference to the magnetization variable.
 * @param rand_vect A vector of random numbers.
 * @param J The interaction strength.
 * @return The final magnetization per site after simulation.
 */

float Wolff::simulate(float T,std::vector <int> & lattice, float& energy, float& M,std::vector<int>& rand_vect,  float J) {
    //using namespace std;
    //vector<float> prob(2);
    //prob[0] = exp(-4 * J/T);
    //prob[1] = exp(-8 * J/T);
    std::vector<float> energy_vec(1);
    energy_vec[0] =  energy;
    std::vector<float> m(1);
    m[0] = (float)M/N;
    std::vector<int> t_axis(1);
    t_axis[0] = 0;
    int n;
    for (int i = 1; i < IT; i++) {
        update(lattice, T, J);
        if (i % N == 0) {
            m.push_back((float) M / N);
            energy_vec.push_back(energy);
            t_axis.push_back(i / N);
        }
    }
    return m.back();
}
/**
 * @brief Creates a vector of random integers.
 * @param rand_vect_0 A reference to the vector where random integers will be stored.
 */

void Wolff::create_rand_vect(std::vector<int>& rand_vect_0) {
    int i;
    for (i = 0; i < IT; i++) {
        rand_vect_0.push_back(rand() % N );
    }
}



/**
 * @brief Stores the results of the simulation to a file.
 */

void Wolff::store_results_to_file() const {
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

