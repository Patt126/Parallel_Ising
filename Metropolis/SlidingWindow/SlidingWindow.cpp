
#include "SlidingWindow.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include <random>


SlidingWindow::SlidingWindow(float interactionStrength, int latticeSize ,int NUMTHREAD , float T_MIN, float T_MAX, float T_STEP, long int IT)
: lattice(interactionStrength, latticeSize),  
        RandVect(),
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
        rng(std::random_device{}()),
        dist(0.0, 1.0)
{   
    // Initialize tStart with std::make_unique
    ThreadStart = std::make_unique<std::vector<int> >();
    ThreadStart->resize(NUMTHREAD);
    // Initialize A with set_block_size function
    A = set_block_size();
    NumSlide = ceil(A/2);
    if (A == -1) {
    std::cerr << "Error: set_block_size failed\n";
    exit(EXIT_FAILURE);
    }
    omp_set_num_threads(NUMTHREAD);
    // Initialize RandVect with std::make_unique
    RandVect = std::make_unique<std::vector<int> >();
    RandVect->resize(0);
    // Initialize EnergyResults with std::make_unique       
    EnergyResults = std::make_unique<std::vector<float> >();
    // Initialize MagnetizationResults with std::make_unique
    MagnetizationResults = std::make_unique<std::vector<float> >();
    // Initialize Temperatures with std::make_unique
    Temperatures = std::make_unique<std::vector<float> >();


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
             while (T < T_MAX) {
                // Initialize prob for the current temperature
                prob[0] = std::exp(-4 * lattice.get_interaction_energy() / T);
                prob[1] = std::exp(-8 * lattice.get_interaction_energy() / T);
                
                deltaE = 0; //shared variable
                deltaM = 0;
                for(int slide = 0; slide < NumSlide; slide++ ){
                    E_loc = 0; //private variable
                    M_loc = 0; 
                    for(int taskNum = 0; taskNum < NUMTHREAD; taskNum++ ){
                        #pragma omp task  firstprivate(M_loc, E_loc) shared(deltaM,deltaE)
                        {
                            //create_rand_vector(randVector_private, rng_private);
                            simulate_step(prob,lattice.get_lattice(), M_loc, E_loc, (*ThreadStart)[taskNum]);
                            #pragma omp atomic update
                            deltaM += M_loc;
                            #pragma omp atomic update
                            deltaE += E_loc;
                        }
                    }
                    #pragma omp taskwait
                    //periodically translate the matrix to update freezed boundaries
                    translate_matrix(lattice.get_lattice());  
                }      
            // Update the magnetization and energy of lattice
            lattice.increment_magnetization(deltaM);
            lattice.increment_energy(deltaE*lattice.get_interaction_energy());
            m = static_cast<float>(lattice.get_magnetization()) / N;
            
            //store the results
            Temperatures->emplace_back(T);
            T += T_STEP;
            EnergyResults->emplace_back(lattice.get_energy());
            MagnetizationResults->emplace_back(abs(m));
            lattice.restore_random_lattice();
        }
        }
    }
}


void SlidingWindow::create_rand_vector() {
    RandVect->clear();
    for (int j = 0;j<(ceil(IT/(NUMTHREAD*NumSlide)));j++) {
        int r = static_cast<int>(dist(rng) * A);
        int c = static_cast<int>(dist(rng) * A);
        if(r%A != 0 && c%A!=0 ){ //not last or first row or column
            RandVect->emplace_back(r * L + c);
        }
        else{
            RandVect->emplace_back(-1); //will mean do nothing
        }
            
 }
}

//Use for create in parallel a random vector of size IT/(NUMTHREAD*NumSlide) with random sites to flip
void SlidingWindow::create_rand_vector(std::vector<int>& randVector, std::mt19937& rng_private) {
    for (int i = 0; i < ceil(IT/(NUMTHREAD*NumSlide)); i++) {
        int r = static_cast<int>(dist(rng_private) * A);
        int c = static_cast<int>(dist(rng_private) * A);
        if(r%A != 0 && c%A!=0 ){ //not last or first row or column
            randVector.emplace_back(r * L + c);
        }
        else{
            randVector.emplace_back(-1); //will mean do nothing
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
                (*ThreadStart)[i] = (floor(i/THREADPERSIDE)*A_loc*L + i%THREADPERSIDE*A_loc); //index at which each block start
            }
            return A_loc;
        }
        else {
            std::cerr << "Error: Unable to fill the line with the given number of blocks." << std::endl;
            return -1;
        }
    }
}

//translate the matrix by a vector (1,1) with periodic boundary conditions
void SlidingWindow::translate_matrix(std::vector<int>& inputMatrix) {
    std::vector<int> localCopy(N);
    std::copy(inputMatrix.begin(), inputMatrix.end(), localCopy.begin());
//#pragma omp taskloop grainsize(N/NUMTHREAD)
    for (int i = 0; i < N; ++i) {
        // Check if the new indices are within bounds
        if ((i + L) < N) { // Check if the new index is not on the bottom boundary
            if ((i + 1) % L != 0) // Check if the new index is not on the right boundary
                inputMatrix[i + L + 1] = localCopy[i];
            else // If it is, then the new index is the left boundary
                inputMatrix[i + 1] = localCopy[i];
        } else if ((i + 1) % L != 0) { // Check if the new index is not on the right boundary but is on the bottom boundary
            inputMatrix[i + L + 1 - N] = localCopy[i];
        } else { //if it is the last corner move to beginning corner
            inputMatrix[0] = localCopy[i];
        }
    }
    
}
//flip the spin of the site and update the energy and magnetization
void SlidingWindow::flip(std::vector<int>& lattice, std::array<float, 2>& prob, int site, int& M, int& E) {
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
        float rnd = dist(rng);
        if (rnd < prob[0]) {
            lattice[site] = -lattice[site];
        } else {
            return;
        }
    } else if (delta == 8) {
        float rnd = dist(rng);
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
void SlidingWindow::simulate_step (std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E, int offset) {
    int n;
    std::random_device rd;
    std::mt19937 rng_private(rd());
    std::uniform_int_distribution<int> dist_private(0, A - 1);
    /*auto randVector_private = std::make_unique<std::vector<int>>();
    randVector_private->reserve(ceil(IT/(NUMTHREAD*NumSlide)));
    create_rand_vector(*randVector_private, *rng_private);*/
    int r ,c ;
    for (unsigned long int i = 0; i < ceil((IT/NUMTHREAD));i++) {
        int r = dist_private(rng_private);
        int c = dist_private(rng_private);
        
        //if not boundary
        if (r  != 0 && c != 0 && r != A-1 && c != A-1) {
            n = r * L + c;
            flip(lattice, prob, n + offset, M, E);
            }
    }
}

void SlidingWindow::simulate_step (std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E, int offset,std::vector<int>& randVector_private) {
    int n;
    for (unsigned long int i = 0; i < ceil(IT/(NUMTHREAD*NumSlide));i++) {
        n = randVector_private[i];
        //if (not boundary) 
        if (n != -1){
            flip(lattice, prob, n + offset, M, E);
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
    std::size_t numResults = EnergyResults->size(); //they all have same lenght
    // Write results to the file
    for (std::size_t i = 0; i < numResults; ++i) {
        // Write data for each row
        outFile << (*EnergyResults)[i] << " " << (*MagnetizationResults)[i] << " " << (*Temperatures)[i]  << std::endl;

    }

    // Close the file
    outFile.close();
}



