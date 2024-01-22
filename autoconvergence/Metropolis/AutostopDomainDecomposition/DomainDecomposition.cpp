
#include "DomainDecomposition.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include <random>
#include <numeric>


DomainDecomposition::DomainDecomposition(float interactionStrength, int latticeSize ,int NUMTHREAD , float T_MIN, float T_MAX, float T_STEP, float tolerance)
: lattice(interactionStrength, latticeSize),  
        EnergyResults(),
        MagnetizationResults(),
        Temperatures(),
        MontecarloStepResults(),
        T_STEP(T_STEP) ,
        T_MIN(T_MIN),
        T_MAX(T_MAX),
        L(latticeSize),
        N(latticeSize * latticeSize),
        tolerance(tolerance),
        NUMTHREAD(NUMTHREAD),
        mExact(),
        A(),
        continue_flag(true),
        ThreadStart(),
        rngPrivate(NUMTHREAD),
        dist(0.0, 1.0)
        
{   
    // Initialize tStart with std::make_unique
    ThreadStart.resize(NUMTHREAD);
     // Initialize A with set_block_size function
    A = set_block_size();
    if (A == -1) {
    std::cerr << "Error: set_block_size failed\n";
    exit(EXIT_FAILURE);
    }
    omp_set_num_threads(NUMTHREAD);

     for (int i = 0; i < NUMTHREAD; ++i) {
            // Choose a seed value, for example, based on the current time
            std::random_device rd;
            std::seed_seq seedSequence{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
            rngPrivate[i].seed(seedSequence);
        }

}


void DomainDecomposition::simulate_phase_transition() {
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
    std::vector<int> iterations;
    iterations.resize(NUMTHREAD);

#pragma omp parallel 
    {

        #pragma omp single nowait
        {
             while (T < T_MAX) {
            
                // Initialize prob for metropolis algorithm at temperature T
                prob[0] = std::exp(-4 * lattice.get_interaction_energy() / T);
                prob[1] = std::exp(-8 * lattice.get_interaction_energy() / T);
                mExact = m_exact(T);
                deltaE = 0; //shared variables for update
                step = 0;
                M = lattice.get_magnetization();
                continue_flag = true;
                std::fill(iterations.begin(), iterations.end(), 0);
                m = static_cast<float>(lattice.get_magnetization()) / N;
                error = std::abs(std::abs(m) - mExact); //evaluate error
                while (error > tolerance ) { //continue simulation until convergence criterion is met 
                //create a random vector for each temperature
                    E_loc = 0; 
                    M_loc = 0; 
                    for(int taskNum = 0; taskNum < NUMTHREAD; taskNum++ ){
                        // create and assign a task to the thread
                        #pragma omp task  firstprivate(E_loc,M_loc) shared(M,iterations,rngPrivate)
                        {
                            //perform the simulation for a block
                            iterations[taskNum] += simulate_step(prob,lattice.get_lattice(), M_loc, E_loc, ThreadStart[taskNum], taskNum);
                        #pragma omp atomic update
                            M += M_loc;    
                        }
                    }
                    #pragma omp taskwait
                    m = static_cast<float>(M) / N;
                    error = std::abs(std::abs(m) - mExact);
                    continue_flag = true;
                }      
            step  = std::accumulate(iterations.begin(), iterations.end(), 0);  
            //update lattice variables
            lattice.increment_magnetization(M-lattice.get_magnetization());
            lattice.increment_energy(deltaE*lattice.get_interaction_energy());
            m = static_cast<float>(lattice.get_magnetization()) / N;
            //store results
            Temperatures.emplace_back(T);
            std::cout << "T = " << T << " m = " << m << " error = " << error << " step = " << step/N << std::endl;
            T += T_STEP;
            EnergyResults.emplace_back(lattice.get_energy());
            MagnetizationResults.emplace_back(abs(m));
            lattice.restore_random_lattice();
            MontecarloStepResults.emplace_back(step/N);
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
void DomainDecomposition::atomic_flip(std::vector<int>& lattice, std::array<float, 2>& prob, const int& site, int& M_loc, int& E, std::mt19937& rngPrivate) {
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
    if (delta <= 0 ) {
#pragma omp atomic write
        lattice[site] = -lattice[site];
    }

    else if (delta == 4 ) {
        float rnd = dist(rngPrivate);
        if (rnd < prob[0]  ){
#pragma omp atomic write
            lattice[site] = -lattice[site];
        }
        else{
            return ;
        }
    }
    else if (delta==8 ){
        float rnd = dist(rngPrivate);
        if (rnd < prob[1]) {
#pragma omp atomic write
            lattice[site] = -lattice[site];
        }
        else{
            return ;
        }
    }
    M_loc += 2 * lattice[site];
    E += delta;
}
//flip the spin of the site and update the energy and magnetization 
//implement the metropolis algorithm described in the report
void DomainDecomposition::flip(std::vector<int>& lattice, std::array<float, 2>& prob, const int& site, int& M_loc, int& E, std::mt19937& rngPrivate) {
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
        lattice[site] = -lattice[site];
    }

    else if (delta == 4 ) {
        float rnd = dist(rngPrivate);
        if (rnd < prob[0]  ){
            lattice[site] = -lattice[site];
        }
        else{
            return ;
        }
    }
    else if (delta==8 ){
        float rnd = dist(rngPrivate);
        if (rnd < prob[1] ) {
            lattice[site] = -lattice[site];

        }
        else{
            return ;
        }
    }
    M_loc += 2 * lattice[site];
    E += delta;
    
}

// simulate the operation for one block of the lattice each Temperature
int DomainDecomposition::simulate_step (std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E, const int& offset) {
    int n;
    //create a private random vector and rng engine for each thread
    std::random_device rd;
    std::mt19937 rngPrivate(rd());
    std::uniform_int_distribution<int> dist_private(0, A - 1);

    int local_M = 0;
    float error = 0;
    float m ;
    int i = 0;
    int r ,c ;
    while (continue_flag){
        int r = dist_private(rngPrivate);
        int c = dist_private(rngPrivate);
        n = r * L + c;
        //if boundary flip atomically
        if ((r == 0 || c == 0 || r == A - 1 || c == A - 1)) {
            
            atomic_flip(lattice, prob, n + offset, local_M, E, rngPrivate);
        }
        else {
            flip(lattice, prob, n + offset, local_M, E, rngPrivate);
        }
        #pragma omp atomic update
            M += local_M;
        if (i % (N/NUMTHREAD) == 0 && continue_flag){
            #pragma omp barrier
            #pragma omp critical
            {

            m = static_cast<float>(M) / N;
            error = std::abs(std::abs(m) - mExact); 
            if (error < tolerance ){
            #pragma omp atomic write
                continue_flag = false;
            }
            }
            
        }
        i++;
    }
    return i;
}

// simulate the operation for one block of the lattice each Temperature



   



int DomainDecomposition::simulate_step (std::array<float, 2> prob, std::vector<int>& lattice, int& M_loc, int& E_loc,const int& offset,const int& numBlock) {
        int n;
        std::uniform_int_distribution<int>dist_A(0, A - 1);
        int i = 0;
        int r ,c ;
        while (continue_flag){
                int r = dist_A(rngPrivate[numBlock]);
                int c = dist_A(rngPrivate[numBlock]);
                n = r * L + c;
                //if boundary flip atomically
                if ((r == 0 || c == 0 || r == A - 1 || c == A - 1) ) {
                    atomic_flip(lattice, prob, n + offset, M_loc ,E_loc, rngPrivate[numBlock]);
                }
                else {
                    flip(lattice, prob, n + offset, M_loc, E_loc, rngPrivate[numBlock]);
                }
            if ( (numBlock == 0) && (i % (2*N) == 0)  && (i!=0)){
                continue_flag = false;
            }
            i++;
        }
        
        return i;
    }



//store results in a file
void DomainDecomposition::store_results_to_file() const {

    std::string filePath = "Results/Domain/Result_" + std::to_string(L) + ".txt";
    // Open the file for writing
    std::ofstream outFile(filePath);

    // Check if the file is open
    if (!outFile.is_open()) {
        // Handle error: unable to open the file
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
        return;
    }

    // Write column headers
    outFile << "E  M  T  N " << std::endl;

    // Determine the number of results to write
    std::size_t numResults = EnergyResults.size(); //they all have same lenght
    // Write results to the file
    for (std::size_t i = 0; i < numResults; ++i) {
        // Write data for each row
        outFile << EnergyResults[i] << " " << MagnetizationResults[i] << " " << Temperatures[i] << " " << MontecarloStepResults[i] << std::endl;
    }

    // Close the file
    outFile.close();
}