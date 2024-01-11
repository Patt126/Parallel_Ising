/**
 * @file main.cpp
 * @brief Benchmarking different implementations of Ising Monte Carlo.
 */

#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>
#include "Metropolis/DomainDecomposition/DomainDecomposition.h"

/**
 * @brief Stores the performance results to a file.
 * @param time_results The vector of time results.
 * @param size The vector of sizes.
 * @param filename The name of the file to store the results.
 */
void store_performance_to_file(std::vector<float> time_results, std::vector<int> size, std::string filename) {
    std::ofstream file;
    file.open("Performance_" + filename + ".txt");

    // Write time results and sizes to the file
    for (int i = 0; i < time_results.size(); i++) {
        file << size[i] << " " << time_results[i] << std::endl;
    }

    file.close(); // Close the file
}

/**
 * @brief Main function for benchmarking a single Ising Monte Carlo implementations.
 * @return 0 on successful execution.
 */
int main() {
    // User input variables
    int L_MIN, L_MAX, NUMTHREAD;
    float T_MIN, T_MAX, T_STEP, interactionStrength;
    long int IT;

    // Taking user input for variables
    std::cout << "Enter minimum L: ";
    std::cin >> L_MIN;

    std::cout << "Enter maximum L: "; 
    std::cin >> L_MAX;

    std::cout << "Enter minimum temperature (T_MIN): ";
    std::cin >> T_MIN;

    std::cout << "Enter maximum temperature (T_MAX): ";
    std::cin >> T_MAX;

    std::cout << "Enter temperature step (T_STEP): ";
    std::cin >> T_STEP;

    std::cout << "Enter interaction strength: ";
    std::cin >> interactionStrength;

    std::cout << "Enter number of threads \n"
              << "(please for organization reason enter a perfect square and a number that allows covering the whole lattice in blocks): ";
    std::cin >> NUMTHREAD;

    // Vector to store time results and sizes
    std::vector<float> time_results;
    std::vector<int> size;

    // Serial Metropolis
    for (int L = L_MIN; L <= L_MAX; L *= 2) {
        IT = pow(L, 4.4);
        std::cout << "Simulation start for L = " << L << std::endl;
        DomainDecomposition simulation(interactionStrength, L, NUMTHREAD, T_MIN, T_MAX, T_STEP, IT);

        auto start = std::chrono::high_resolution_clock::now();

        simulation.simulate_phase_transition();

        auto stop = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

        time_results.emplace_back(duration.count());
        size.emplace_back(L * L);
        simulation.store_results_to_file();
    }

    // Store performance results for Serial Metropolis
    store_performance_to_file(time_results, size, "DomainDecomposition");
    time_results.clear();
    size.clear();
}