
#include "DomainDecomposition.h"
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include <random>
#include<mpi.h>


DomainDecomposition::DomainDecomposition(float interactionStrength, int widht,int heigth ,int NUMTHREAD,int MPISIZE , float T_MIN, float T_MAX, float T_STEP, long int IT)
: lattice(interactionStrength, widht, heigth),  
        EnergyResults(),
        MagnetizationResults(),
        Temperatures(),
        T_STEP(T_STEP) ,
        T_MIN(T_MIN),
        T_MAX(T_MAX),
        W(widht),
        H(heigth),
        N(MPISIZE * widht * heigth),
        localN(widht*heigth),
        IT(IT),
        NUMTHREAD(NUMTHREAD),
        MPISIZE(MPISIZE),
        numSlide(heigth/2),
        A(),
        ThreadStart(),
        dist(0.0, 1.0)
{   
    ThreadStart.resize(NUMTHREAD);
    // Initialize A with set_block_size function
    A = set_block_size();    
    localN = N/MPISIZE;

}


void DomainDecomposition::simulate_phase_transition() {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(static_cast<unsigned>(time(nullptr)) + rank);

    // Scatter the relevant part of the linear vector among MPI processes
   
    int E_loc ;
    int M_loc ;
    int deltaE ;
    int deltaM ;
    int totalDeltaM = 0;
    int totalDeltaE = 0;
    float m = 0;
    float T = T_MIN;
    int M = lattice.get_magnetization();
    std::array<float, 2> prob;
    

#pragma omp parallel num_threads(NUMTHREAD/MPISIZE)
    {
        #pragma omp single nowait
        {
             int world_size, world_rank;
            MPI_Comm_size(MPI_COMM_WORLD, &world_size);
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
             while (T < T_MAX) {
                // Initialize prob for metropolis algorithm at temperature T
                prob[0] = std::exp(-4 * lattice.get_interaction_energy() / T);
                prob[1] = std::exp(-8 * lattice.get_interaction_energy() / T);
                E_loc = 0; //local variables for each thread
                M_loc = 0; 
                deltaE = lattice.get_energy(); //shared variables for update
                deltaM = lattice.get_magnetization();
                for(int slide = 0; slide<numSlide; slide++){
                    for(int taskNum = 0; taskNum < NUMTHREAD; taskNum++ ){
                        // create and assign a task to the thread
                        #pragma omp task  firstprivate(M_loc, E_loc) shared(deltaM,deltaE)
                        {
                           
                            simulate_step(prob,lattice.get_lattice(), M_loc, E_loc, ThreadStart[taskNum]);
                            
                            #pragma omp atomic update //update shared variables
                            deltaM += M_loc;
                            #pragma omp atomic update
                            deltaE += E_loc;
                        }
                    }
                    #pragma omp taskwait
                    exchange_rows();
                 } 
                  
            totalDeltaE = 0;
            totalDeltaM = 0;
            MPI_Allreduce(&deltaM, &totalDeltaM, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&deltaE, &totalDeltaE, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            //store results
            Temperatures.emplace_back(T);
            EnergyResults.emplace_back(lattice.get_interaction_energy()*totalDeltaE);
            m = static_cast<float>(totalDeltaM) / N;
            MagnetizationResults.emplace_back(fabs(m));

             
             T += T_STEP;
             lattice.restore_random_lattice();
             

        }
    }
}
}

void DomainDecomposition::print_full_lattice() {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> localLattice = lattice.get_lattice();


    std::vector<MPI_Request> recvRequests(size - 1);

    if (rank != 0) {
        // Processes except the first send their local lattice to the first process
        MPI_Request sendRequest;
        MPI_Isend(localLattice.data(), localN, MPI_INT, 0, 0, MPI_COMM_WORLD, &sendRequest);
    }
    

    // First process collects data from all processes
    else {
        //  total lattice
        std::vector<int> totalLattice(N, 0);

        // Copy the local lattice of the first process i
        std::copy(localLattice.begin(), localLattice.end(), totalLattice.begin());

        // Receive and store the lattices from other processes
        for (int source = 1; source < size; ++source) {
            //  non-blocking receive
            MPI_Irecv(totalLattice.data() + source * localN, localN, MPI_INT, source, 0, MPI_COMM_WORLD, &recvRequests[source - 1]);
        }

        // Wait for all non-blocking receives to complete and get status
        //MPI_Waitall(size - 1, recvRequests.data(), MPI_STATUS_IGNORE);

        // Print the full lattice (only the first process should print)
        for (int i = 0; i < W; ++i) {
            for (int j = 0; j < H*MPISIZE; ++j) { //often W=H*MPISIZE
                int index = i * W + j;
                if(totalLattice[index] == -1){
                    std::cout<<"o ";
                }
                else{
                    std::cout<<"x ";
                }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


}


// organize division of the lattice in blocks
int DomainDecomposition::set_block_size() { 
    int A_loc = sqrt(localN / NUMTHREAD);
    int THREADPERSIDE = W/A_loc;
    for(int i=0;i<  NUMTHREAD;i++){
                ThreadStart[i] = (floor(i/THREADPERSIDE)*A_loc*W + i%THREADPERSIDE*A_loc); //thread at which each block start
            }
    return A_loc;
  
}
//flip the spin of the site and update the energy and magnetization 
//atomic version of the flip function, used to flip boundary sites.
void DomainDecomposition::atomic_flip(std::vector<int>& lattice, std::array<float, 2>& prob, const int& site, int& M, int& E,std::mt19937& rng_private) {
    int sum = 0;

    if (site > W && site < W*(H-1)) { //if not boundary of the whole rectangle
        sum += lattice[site-W];
        sum += lattice[site + W];
        if (site % W == 0) {
        sum += lattice[site + (W - 1)];
    }
        else {
            sum += lattice[site - 1];
        }

        if ((site+1) % W == 0) {
            sum += lattice[site - (W-1)];
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
    
   
    

}
//flip the spin of the site and update the energy and magnetization 
//implement the metropolis algorithm described in the report
void DomainDecomposition::flip(std::vector<int>& lattice, std::array<float, 2>& prob, const int& site, int& M, int& E, std::mt19937& rng_private) {
    int sum = 0; //note that here when we call flip and not atomicflip we arleady know that is not border

    sum += lattice[site - W];
    sum += lattice[site - 1];
    sum += lattice[site + W];
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
    for (unsigned long int i = 0; i < ceil((IT/(NUMTHREAD*numSlide)));i++) {
        int r = dist_private(rng_private);
        int c = dist_private(rng_private);
        n = r * W + c;
        //if boundary flip atomically
        if (r == 0 || c == 0 || r == A - 1 || c == A - 1) {
            
            atomic_flip(lattice, prob, n + offset, M, E, rng_private);
        }
        else{
            flip(lattice, prob, n + offset, M, E,  rng_private);
        }
    }
}


void DomainDecomposition::exchange_rows() {
    // Determine the ranks of the adjacent processes
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int prev_rank = (world_rank - 1 + world_size) % world_size;
    int next_rank = (world_rank + 1) % world_size;

    // Allocate memory for the send and receive buffers
    std::vector<int> send_buffer(2 * W);
    std::vector<int> recv_buffer(2 * W);

    // Copy the last two rows to the send buffer
    std::copy(lattice.get_lattice().end() - 2 * W, lattice.get_lattice().end(), send_buffer.begin());

    // Use MPI_Sendrecv to simultaneously send the last two rows to the next process
    // and receive two rows from the previous process
    MPI_Sendrecv(send_buffer.data(), 2 * W, MPI_INT, next_rank, 0,
                 recv_buffer.data(), 2 * W, MPI_INT, prev_rank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Move all rows in the lattice down by two rows
    std::copy(lattice.get_lattice().begin(), lattice.get_lattice().end() - 2 * W, lattice.get_lattice().begin() + 2 * W);

    // Copy the received rows to the beginning of the lattice
    std::copy(recv_buffer.begin(), recv_buffer.end(), lattice.get_lattice().begin());
}


//store results in a file
void DomainDecomposition::store_results_to_file() const {
    std::string filePath = "./Results/Domain/Result_" + std::to_string(W) + ".txt";
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


