#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>
#include <sstream>
#include <mpi.h>
#include "Metropolis/MPIDomainDecomposition/MPIDomainDecomposition.h"


void store_performance_to_file(std::vector<float> time_results,int L , std::string filename) {
    std::ofstream file;
    file.open("./Performance/" + filename + "_L_" + std::to_string(L) + ".txt", std::ios::trunc);

    file <<"Domain Decompositon: " <<time_results[0] << " s"<< std::endl;

    file.close();
}

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int L_TOTAL = 64; //total side of the complessive square lattice, must be divide by MPISIZE
    int NUMTHREAD = 4; // Number of threads per MPI process
    float T_MIN = 0.1; 
    float T_MAX = 0.2;
    float T_STEP = 0.1; 
    float interactionStrength = 1.0;
    long int IT;
    std::string filename = "_test_1";
    std::vector<float> time_results;

    //check the possibility of cover the lattice with the thread and process indicated
    if (world_rank == 0) {
    int GLOBAL_NUMTHREAD = world_size*NUMTHREAD;
    int THREADPERSIDE = sqrt(GLOBAL_NUMTHREAD);
   
    if(THREADPERSIDE*THREADPERSIDE != GLOBAL_NUMTHREAD){
        std::cerr << "Error: Global number of threads must be a perfect square." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    else if (L_TOTAL % THREADPERSIDE != 0) {
        std::cerr << "Error: Unable to fill a row with the given number of blocks." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    else if(THREADPERSIDE % world_size != 0) {
        std::cerr << "Error: Number of block in a side must be a multiple of the number of processes." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }



    }

    // Determine the width and heigth of the rectangular lattice for each MPI process
    int width = L_TOTAL;
    int heigth = L_TOTAL/world_size;
    
    // Loop over lattice sizes
        IT = ceil(pow(L_TOTAL, 4.7));

        // instantiate for each process an object
        MPIDomainDecomposition domainSimulation(interactionStrength, width,heigth, NUMTHREAD,world_size, T_MIN, T_MAX, T_STEP, IT);
        auto startDomain = std::chrono::high_resolution_clock::now();
        domainSimulation.simulate_phase_transition();
        auto stopDomain = std::chrono::high_resolution_clock::now();
        
        auto durationDomain = std::chrono::duration_cast<std::chrono::seconds>(stopDomain - startDomain);
        time_results.push_back(durationDomain.count());

        //store results
        domainSimulation.store_results_to_file();
        
        if(world_rank == 0)
        {
            store_performance_to_file(time_results, L_TOTAL, filename);
            time_results.clear();
        }



    
    MPI_Finalize();    

    return 0;
}