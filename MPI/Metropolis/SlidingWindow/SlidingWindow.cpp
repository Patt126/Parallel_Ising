
#include "SlidingWindow.h"
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include <random>


SlidingWindow::SlidingWindow(float interactionStrength, int latticeSize ,int NUMTHREAD , float T_MIN, float T_MAX, float T_STEP, long int IT)
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
        NumSlide(),
        dist(0.0, 1.0)
{   
    
    ThreadStart.resize(NUMTHREAD);
    // Initialize A with set_block_size function
    A = set_block_size();
    NumSlide = ceil(4*A);
    if (A == -1) {
    std::cerr << "Error: set_block_size failed\n";
    exit(EXIT_FAILURE);
    }
    omp_set_num_threads(NUMTHREAD);

}


void SlidingWindow::simulate_phase_transition() {
    int E_loc ;
    int M_loc ;
    int deltaE ;
    int deltaM ;
    int step;
    float m = 0;
    float T = T_MIN;
    int M = lattice.get_magnetization();
    float error = 0;
    std::array<float, 2> prob;
    

#pragma omp parallel 
    {
        
        /*std::vector<int> randVector_private;

        randVector_private.reserve(ceil(IT/(NUMTHREAD*NumSlide)));
        #pragma omp barrier*/
        
    
        #pragma omp single nowait
        {
            std::cout<<"sliding window"<<std::endl;
             while (T < T_MAX) {

                // Initialize prob for the current temperature
                prob[0] = std::exp(-4 * lattice.get_interaction_energy() / T);
                prob[1] = std::exp(-8 * lattice.get_interaction_energy() / T);
                
                deltaE = 0; //shared variable
                deltaM = 0;
                for(int slide = 0; slide < NumSlide; slide++ ){
                    E_loc = 0; //private variable
                    M_loc = 0; 
                    deltaM = 0;
                    for(int taskNum = 0; taskNum < NUMTHREAD; taskNum++ ){
                        #pragma omp task  firstprivate(M_loc, E_loc) shared(deltaM,deltaE)
                        {
                            //create_rand_vector(randVector_private, rng_private);
                            simulate_step(prob,lattice.get_lattice(), M_loc, E_loc, ThreadStart[taskNum]);
                            #pragma omp atomic update
                            deltaM += M_loc;
                            #pragma omp atomic update
                            deltaE += E_loc;
                        }
                    }
                    #pragma omp taskwait

                    lattice.increment_magnetization(deltaM);
                    m = static_cast<float>(lattice.get_magnetization()) / N;

                    //periodically translate the matrix to update freezed boundaries
                    translate_matrix(lattice.get_lattice());  
                }      
            // Update the magnetization and energy of lattice
            lattice.increment_magnetization(deltaM);
            lattice.increment_energy(deltaE*lattice.get_interaction_energy());
            m = static_cast<float>(lattice.get_magnetization()) / N;
            
            //store the results
            Temperatures.emplace_back(T);
            T += T_STEP;
            EnergyResults.emplace_back(lattice.get_energy());
            MagnetizationResults.emplace_back(abs(m));
            lattice.restore_random_lattice();
        }
        }
    }
}



//organice division of the lattice in blocks
int SlidingWindow::set_block_size() {
    int THREADPERSIDE = sqrt(NUMTHREAD);
    if(THREADPERSIDE*THREADPERSIDE != NUMTHREAD){
        std::cerr << "Error: Number of threads must be a perfect square." << std::endl;
        return -1;
    }
    else{
        if (L % THREADPERSIDE == 0) {
            int A_loc = L / THREADPERSIDE; // = single block width
        for(int i=0;i<  NUMTHREAD;i++){
                ThreadStart[i] = (floor(i/THREADPERSIDE)*A_loc*L + i%THREADPERSIDE*A_loc); //index at which each block start
            }
            return A_loc;
        }
        else {
            std::cerr << "Error: Unable to fill the line with the given number of blocks." << std::endl;
            return -1;
        }
    }
}

//translate the matrix by a vector (2,2) with periodic boundary conditions
void SlidingWindow::translate_matrix(std::vector<int>& inputMatrix) {
    std::vector<int> localCopy(N);
    std::copy(inputMatrix.begin(), inputMatrix.end(), localCopy.begin());

    for (int i = 0; i < N; ++i) {
        // Calculate the new indices with periodic boundary conditions
        int newRow = (i / L + 2) % L; //consider that being i and L integers, (i/L) is the floor of the division
        int newCol = (i % L + 2) % L;

        // Translate the matrix element to the new index
        int newIndex = newRow * L + newCol;
        inputMatrix[newIndex] = localCopy[i];
    }
}

//flip the spin of the site and update the energy and magnetization 
//implement the metropolis algorithm described in the report
void SlidingWindow::flip(std::vector<int>& lattice, std::array<float, 2>& prob, const int& site, int& M, int& E, std::mt19937& rng_private) {
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

// simulate the operation for one block freezing boundary flips
void SlidingWindow::simulate_step (std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E, const int& offset) {
    int n;
    std::random_device rd;
    std::mt19937 rng_private(rd());
    std::uniform_int_distribution<int> dist_private(0, A - 1);
    int r ,c ;
    for (unsigned long int i = 0; i < ceil((IT/(NUMTHREAD*NumSlide)));i++) {
        int r = dist_private(rng_private);
        int c = dist_private(rng_private);
        
        //if not boundary
        if (r  != 0 && c != 0 && r != A-1 && c != A-1) {
            n = r * L + c;
            flip(lattice, prob, n + offset, M, E, rng_private);
            }
    }
}



//store the results in a file
void SlidingWindow::store_results_to_file() const {

   std::string filePath = "./Results/Sliding/Result_" + std::to_string(L) + ".txt";
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



