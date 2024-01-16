// AutostopMontecarlo.tpp

#include "AutostopMontecarlo.h"
#include <cstdlib>
#include <fstream>


AutostopMontecarlo::AutostopMontecarlo(float interactionStrength, int latticeSize,  float Tolerance, float T_MIN, float T_MAX, float T_STEP):
     lattice(interactionStrength, latticeSize),  
        randVect(),
        energyResults(),
        magnetizationResults(),
        monteCarloStepsResults(),
        temperatures(),
        tolerance(Tolerance),  // Initialize tolerance with the specified value
        T_STEP(T_STEP) ,
        T_MIN(T_MIN),
        T_MAX(T_MAX),
        L(latticeSize),
        N(latticeSize * latticeSize)
{   
    randVect = std::make_unique<std::vector<int> >(N);
    // Initialize randVect with createRandVect function
    create_rand_vector();
    // Initialize energyResults with std::make_unique       
    energyResults = std::make_unique<std::vector<float> >();
    // Initialize magnetizationResults with std::make_unique
    magnetizationResults = std::make_unique<std::vector<float> >();
    // Initialize monteCarloStepsResults with std::make_unique
    monteCarloStepsResults = std::make_unique<std::vector<int> >();
    // Initialize temperatures with std::make_unique
    temperatures = std::make_unique<std::vector<float> >();                     

}


void AutostopMontecarlo::simulate_phase_transition() {
    int deltaE ;
    int deltaM ;
    int step;
    float m = 0;
    float mexact = 1;
    float T = T_MIN;

    float error = 0;
    std::array<float, 2> prob;
    
    while (T < T_MAX) {
        prob[0] = std::exp(-4 * lattice.get_interaction_energy() / T);
        prob[1] = std::exp(-8 * lattice.get_interaction_energy() / T);
        step = 0;
        mexact = m_exact(T); //evaluate exact magnetization
        m = static_cast<float>(lattice.get_magnetization()) / N;
        error = std::abs(std::abs(m) - mexact); //evaluate error
        while (error > tolerance ) { //continue simulation until convergence criterion is met 
            deltaE = 0;
            deltaM = 0; 
            create_rand_vector();
            simulate_step(prob, lattice.get_lattice(), deltaM, deltaE);
            lattice.increment_magnetization(deltaM);
            lattice.increment_energy(deltaE);
            m = static_cast<float>(lattice.get_magnetization()) / N;
            error = std::abs(std::abs(m) - mexact);
            step++;
        
        }
        //record simulation results for the current temperature
        temperatures->emplace_back(T);
        T += T_STEP;
        monteCarloStepsResults->emplace_back(step);
        energyResults->emplace_back(1);
        magnetizationResults->emplace_back(abs(m));
        step = 0;
        lattice.restore_random_lattice(); //restore the lattice to its initial state for the next temperature
    }
}





void AutostopMontecarlo::flip(std::vector<int>& lattice, std::array<float, 2>& prob, int site, int& M, int& E) {
    int sum = 0;

    if (site < L) {
        sum += lattice[site + L * (L - 1)];
    } else {
        sum += lattice[site - L];
    }
    if (site % L == 0) {
        sum += lattice[site + (L - 1)];
    } else {
        sum += lattice[site - 1];
    }

    if (site >= L * (L - 1)) {
        sum += lattice[site - L * (L - 1)];
    } else {
        sum += lattice[site + L];
    }
    if ((site + 1) % L == 0) {
        sum += lattice[site - (L - 1)];
    } else {
        sum += lattice[site + 1];
    }
    int delta = 2 * sum * lattice[site];
    if (delta <= 0) {
        lattice[site] = -lattice[site];
    } else if (delta == 4) {
        float rnd = (rand() % 10000) / 1e4;
        if (rnd < prob[0]) {
            lattice[site] = -lattice[site];
        } else {
            return;
        }
    } else if (delta == 8) {
        float rnd = (rand() % 10000) / 1e4;
        if (rnd < prob[1]) {
            lattice[site] = -lattice[site];
        } else {
            return;
        }
    }
    M += 2 * lattice[site];
}


int AutostopMontecarlo::simulate_step(std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E,int offset ) {
    for (unsigned long int i = 0; i < (N); i++) {
        int n = (*randVect)[i];
        if (n != -1) {
            flip(lattice, prob, n, M, E);
        }
    }
    return 0;
}


void AutostopMontecarlo::store_results_to_file() const {

    std::string filePath = "Results/Serial/Result_" + std::to_string(L) + ".txt";
    // Open the file for writing
    std::ofstream outFile(filePath);

    // Check if the file is open
    if (!outFile.is_open()) {
        // Handle error: unable to open the file
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
        return;
    }

    // Write column headers
    outFile << "E  M  T " << std::endl;

    // Determine the number of results to write
    std::size_t numResults = energyResults->size(); //they all have same lenght
    // Write results to the file
    for (std::size_t i = 0; i < numResults; ++i) {
        // Write data for each row
        outFile << (*energyResults)[i] << " " << (*magnetizationResults)[i] << " " << (*temperatures)[i] <<(*monteCarloStepsResults)[i]  << std::endl;

    }

    // Close the file
    outFile.close();
}


