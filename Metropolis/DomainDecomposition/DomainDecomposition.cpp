
#include "DomainDecomposition.h"
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include <random>


DomainDecomposition::DomainDecomposition(float interactionStrength, int latticeSize ,int NUMTHREAD , float T_MIN, float T_MAX, float T_STEP, long int IT)
: lattice(interactionStrength, latticeSize),  
        EnergyResults(),
        MagnetizationResults(),
        Temperatures(),
        T_STEP(T_STEP) ,
        T_MIN(T_MIN),
        T_MAX(T_MAX),
        L(latticeSize),
        N(latticeSize * latticeSize),
        IT(IT),
        NUMTHREAD(NUMTHREAD),
        A(),
        ThreadStart(),
        dist(0.0, 1.0)
{   
    ThreadStart.resize(NUMTHREAD);
    // Initialize A with set_block_size function
    A = set_block_size();
    if (A == -1) {
    std::cerr << "Error: set_block_size failed\n";
    exit(EXIT_FAILURE);
    }
}


void DomainDecomposition::simulate_phase_transition() {
    int E_loc ;
    int M_loc ;
    int deltaE ;
    int deltaM ;
    float m = 0;
    float T = T_MIN;
    int M = lattice.get_magnetization();
    std::array<float, 2> prob;
    

#pragma omp parallel 
    {
        #pragma omp single nowait
        {
             while (T < T_MAX) {
                // Initialize prob for metropolis algorithm at temperature T
                prob[0] = std::exp(-4 * lattice.get_interaction_energy() / T);
                prob[1] = std::exp(-8 * lattice.get_interaction_energy() / T);
                E_loc = 0; //local variables for each thread
                M_loc = 0; 
                deltaE = 0; //shared variables for update
                deltaM = 0;
                //create a random vector for each temperature
                    for(int taskNum = 0; taskNum < NUMTHREAD; taskNum++ ){
                        // create and assign a task to the thread
                        #pragma omp task  firstprivate(M_loc, E_loc) shared(deltaM,deltaE)
                        {

                            //perform the simulation for a block
                            simulate_step(prob,lattice.get_lattice(), M_loc, E_loc, ThreadStart[taskNum]);
                            #pragma omp atomic update //update shared variables
                            deltaM += M_loc;
                            #pragma omp atomic update
                            deltaE += E_loc;
                        }
                    }
                    #pragma omp taskwait
                  
                  
            //update lattice variables
            lattice.increment_magnetization(deltaM);
            lattice.increment_energy(deltaE*lattice.get_interaction_energy());
            m = static_cast<float>(lattice.get_magnetization()) / N;
            
            //store resulabsts
            Temperatures.emplace_back(T);
            T += T_STEP;
            EnergyResults.emplace_back(lattice.get_energy());
            MagnetizationResults.emplace_back(fabs(m));
            lattice.restore_random_lattice();
             }
        }
    }
}






// organize division of the lattice in blocks
int DomainDecomposition::set_block_size() {
    int THREADPERSIDE = sqrt(NUMTHREAD);
    if(THREADPERSIDE*THREADPERSIDE != NUMTHREAD){
        std::cerr << "Error: Number of threads must be a perfect square." << std::endl;
        return -1;
    }
    else{
        if (L % THREADPERSIDE == 0) {
            int A_loc = L / THREADPERSIDE; // = larghezze di un blocco
        for(int i=0;i<  NUMTHREAD;i++){
                ThreadStart[i] = (floor(i/THREADPERSIDE)*A_loc*L + i%THREADPERSIDE*A_loc); //thread at which each block start
            }
            return A_loc;
        }
        else {
            std::cerr << "Error: Unable to fill the line with the given number of blocks." << std::endl;
            return -1;
        }
    }
}
//flip the spin of the site and update the energy and magnetization 
//atomic version of the flip function, used to flip boundary sites.
void DomainDecomposition::atomic_flip(std::vector<int>& lattice, std::array<float, 2>& prob, const int& site, int& M, int& E,std::mt19937& rng_private) {
    int sum = 0;

    if (site < L) {
        sum += lattice[site+L*(L-1)];
    }
    else {
        sum += lattice[site-L];
    }
    if (site % L == 0) {
        sum += lattice[site + (L - 1)];
    }
    else {
        sum += lattice[site - 1];
    }

    if (site >= L*(L - 1)) {
        sum += lattice[site - L*(L-1)];
    }
    else {
        sum += lattice[site + L];
    }
    if ((site+1) % L == 0) {
        sum += lattice[site - (L-1)];
    }
    else {
        sum += lattice[site + 1];
    }
    int delta = 2*sum*lattice[site];
    if (delta <= 0) {
#pragma omp atomic write
        lattice[site] = -lattice[site];
    }

    else if (delta == 4) {
        float rnd = dist(rng_private);
        if (rnd < prob[0] ){
#pragma omp atomic write
            lattice[site] = -lattice[site];
        }
        else{
            return ;
        }
    }
    else if (delta==8){
        float rnd = dist(rng_private);
        if (rnd < prob[1]) {
#pragma omp atomic write
            lattice[site] = -lattice[site];
        }
        else{
            return ;
        }
    }
    M += 2 * lattice[site];
    E += delta;

}
//flip the spin of the site and update the energy and magnetization 
//implement the metropolis algorithm described in the report
void DomainDecomposition::flip(std::vector<int>& lattice, std::array<float, 2>& prob, const int& site, int& M, int& E, std::mt19937& rng_private) {
    int sum = 0; //we arleady know that we are not out of border
    sum += lattice[site - L];
    sum += lattice[site - 1];
    sum += lattice[site + L];
    sum += lattice[site + 1];

    int delta = 2 * sum * lattice[site];
    if (delta <= 0) {
        lattice[site] = -lattice[site];
    } else if (delta == 4) {
        float rnd = dist(rng_private);
        if (rnd < prob[0]) {
            lattice[site] = -lattice[site];
        } else {
            return;
        }
    } else if (delta == 8) {
        float rnd = dist(rng_private);
        if (rnd < prob[1]) {
            lattice[site] = -lattice[site];
        } else {
            return;
        }
    }
    M += 2 * lattice[site];
    E += delta;
}

// simulate the operation for one block of the lattice each Temperature
void DomainDecomposition::simulate_step ( std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E, const int& offset) {
    int n;
    //create a private random vector and rng engine for each thread
    std::random_device rd;
    std::mt19937 rng_private(rd());
    std::uniform_int_distribution<int> dist_private(0, A - 1);



    int r ,c ;
    for (unsigned long int i = 0; i < ceil((IT/NUMTHREAD));i++) {
        int r = dist_private(rng_private);
        int c = dist_private(rng_private);
        n = r * L + c;
        //if boundary flip atomically
        if (r == 0 || c == 0 || r == A - 1 || c == A - 1) {
            
            atomic_flip(lattice, prob, n + offset, M, E, rng_private);
        }
        else{
            flip(lattice, prob, n + offset, M, E,  rng_private);
        }
    }
}

//store results in a file
void DomainDecomposition::store_results_to_file() const {

 std::string filePath = "./Results/Domain/Result_" + std::to_string(L) + ".txt";
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
    std::size_t numResults = EnergyResults.size(); //they all have same lenght
    // Write results to the file
    for (std::size_t i = 0; i < numResults; ++i) {
        // Write data for each row
        outFile << EnergyResults[i] << " " << MagnetizationResults[i] << " " << Temperatures[i]  << std::endl;

    }

    // Close the file
    outFile.close();
}


