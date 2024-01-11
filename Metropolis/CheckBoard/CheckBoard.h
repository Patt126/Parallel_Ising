#ifndef MY_PROJECT_MONTE_CARLO_SIMULATION_H
#define MY_PROJECT_MONTE_CARLO_SIMULATION_H

#include <array>
#include <vector>
#include <memory>
#include <random>
#include "../AbstractMonteCarloSimulation.h"
#include "../../Lattice/SquareLattice.h"

/**
 * @brief Class representing a Monte Carlo simulation using a checkerboard pattern.
 *
 * This class inherits from AbstractMonteCarloSimulation and provides specific
 * implementations for the simulation steps using a checkerboard pattern, consisting
 * of two sublattices, like black and white in a chessboard. The simulation is performed alternating
 * between the two sublattices. The simulation is performed concurrently using OpenMP tasks.
 * This help when introducing convergence criteria and to menage alternation between sublattices 
 * and synchronization.
 */
class CheckBoard : public AbstractMonteCarloSimulation {
public:
    /**
     * @brief Constructor for the CheckBoard class.
     *
     * @param interactionStrength The interaction strength parameter.
     * @param latticeSize The size of the square lattice.
     * @param NUMTHREAD The number of threads to use in the simulation.
     * @param T_MIN The minimum temperature.
     * @param T_MAX The maximum temperature.
     * @param T_STEP The temperature step size.
     * @param IT The number of Monte Carlo steps.
     */
    CheckBoard(float interactionStrength, int latticeSize, int NUMTHREAD, float T_MIN, float T_MAX, float T_STEP, long int IT);

    /**
     * @brief Perform the simulation of a phase transition using a checkerboard pattern.
     */
    void simulate_phase_transition() override;

    /**
     * @brief Store the results of the simulation to a file.
     */
    void store_results_to_file() const override;

protected:
    /**
     * @brief Create a random vector for the simulation.
     */
    void create_rand_vector() override;

    /**
     * @brief Flip a spin at a given lattice site during the simulation.
     *
     * @param lattice Reference to the lattice vector.
     * @param prob Probability array for the simulation.
     * @param site Index of the lattice site to flip.
     * @param M Reference to the magnetization variable.
     * @param E Reference to the energy variable.
     */
    void flip(std::vector<int>& lattice, std::array<float, 2>& prob, int site, int& M, int& E);

    /**
     * @brief Simulate a step for a single block of the checkerboard,
     * number of iterations is given by NumFlipPerBlock. It is a parameter
     * to be tuned.
     *
     * @param prob Probability array for the simulation.
     * @param lattice Reference to the lattice vector.
     * @param M Reference to the magnetization variable.
     * @param E Reference to the energy variable.
     * @param offset Offset for the starting point of each block.
     */
    void simulate_step(std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E, int offset) override;

    /**
     * @brief Set the block widht for the checkerboard pattern.
     *
     * @return The block size.
     */
    int set_block_size();

private:
    SquareLattice lattice;  // Use SquareLattice as a private member
    std::unique_ptr<std::vector<int>> RandVect;  // Random vector of size N.
    std::unique_ptr<std::vector<float>> EnergyResults;  // Vector to store energy results.
    std::unique_ptr<std::vector<float>> MagnetizationResults;  // Vector to store magnetization results.
    std::unique_ptr<std::vector<float>> Temperatures;  // Vector to store temperature visited
    std::unique_ptr<std::vector<int>> ThreadStart;  // Vector to store Monte Carlo steps
    float T_MIN;
    float T_MAX;
    float T_STEP;
    const int L;
    const int N;
    const long int IT;
    const int NUMTHREAD;
    int A;
    int NumFlipPerBlock;
    std::mt19937 rng;  // Mersenne Twister 19937 generator
    std::uniform_real_distribution<double> dist;  // Uniform distribution in [0, 1)
};

#endif // MY_PROJECT_MONTE_CARLO_SIMULATION_H
