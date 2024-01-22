#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>
#include <sstream>
#include "DomainDecomposition.h"
#include "AutostopMontecarlo.h"

void store_performance_to_file(std::vector<float> time_results,int L , std::string filename) {
    std::ofstream file;
    file.open("./Performance/" + filename + "_L_" + std::to_string(L) + ".txt", std::ios::trunc);

    // Write time results and sizes to the file for the current method
    file <<"Serial: " <<time_results[0] << " s"<< std::endl;
    file <<"Domain Decompositon: " <<time_results[1] << " s"<< std::endl;
    file.close();
}

/**
 * @brief Main function for Benchmarking Ising Monte Carlo implementations.
 * @return 0 on successful execution.
 */
int main() {
    // User input variables
    int L_MIN =64;
    int L_MAX = 128;
    int NUMTHREAD = 16;
    float T_MIN =0.1;
    float T_MAX = 0.6;
    float T_STEP = 0.1;
    float interactionStrength = 1.0;
    float tolerance = 0.05;
    std::string filename;

   /* // Taking user input for variables
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

    std::cout << "Enter % of not corrected aligned spin tollerated: ";
    std::cin >> tolerance;

     std::cout << "Enter filename:";
    std::cin >> filename;
    */
    std::vector<float> time_results;

    // Loop over lattice sizes
    for (int L = L_MIN; L <= L_MAX; L *= 2) {
        std::cout << "Simulation start for L = " << L << std::endl;

        {
            // Serial Metropolis simulation
            AutostopMontecarlo serialSimulation(interactionStrength, L, T_MIN, T_MAX, T_STEP, tolerance);

            auto startSerial = std::chrono::high_resolution_clock::now();
            serialSimulation.simulate_phase_transition();
            auto stopSerial = std::chrono::high_resolution_clock::now();
            serialSimulation.store_results_to_file();
            auto durationSerial = std::chrono::duration_cast<std::chrono::seconds>(stopSerial - startSerial);

            time_results.push_back(durationSerial.count());
        }

        {
            // Domain Decomposition simulation
            DomainDecomposition domainSimulation(interactionStrength, L, NUMTHREAD, T_MIN, T_MAX, T_STEP, tolerance);

            auto startDomain = std::chrono::high_resolution_clock::now();
            domainSimulation.simulate_phase_transition();
            auto stopDomain = std::chrono::high_resolution_clock::now();
            domainSimulation.store_results_to_file();
            auto durationDomain = std::chrono::duration_cast<std::chrono::seconds>(stopDomain - startDomain);

            time_results.push_back(durationDomain.count());

        }

    
        // Store performance results for all methods after completing all iterations
    store_performance_to_file(time_results, L, filename);
    time_results.clear();
    }

    

    return 0;
}
