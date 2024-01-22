/**
 * @file SlidingWindow.h
 * @brief Header file for the SlidingWindow class.
 */

#ifndef MY_PROJECT_SLIDING_WINDOW_H
#define MY_PROJECT_SLIDING_WINDOW_H

#include <array>
#include <vector>
#include <memory>
#include <random>
#include "AbstractMonteCarloSimulation.h"
#include "SquareLattice.h"

/**
 * @brief Class representing a Monte Carlo simulation using a sliding window approach.
 * Lattice is divided in blocks of widht A, and each block is simulated concurrently, BUT, if
 * the flipping site is on the boundary of the block, the flip is not performed.
 * Periodocally the lattice is translated of A sites, so that each site is visited..
 *
 * This class inherits from AbstractMonteCarloSimulation and provides specific
 * implementations for the simulation steps using a sliding window approach.
 */
class SlidingWindow : public AbstractMonteCarloSimulation {
public:
    /**
     * @brief Constructor for the SlidingWindow class.
     *
     * @param interactionStrength The spins interaction strength parameter.
     * @param latticeSize The size of the square lattice.
     * @param NUMTHREAD The number of threads to use in the simulation.
     * @param T_MIN The minimum temperature.
     * @param T_MAX The maximum temperature.
     * @param T_STEP The temperature step size.
     * @param IT The number of Monte Carlo steps.
     */
    SlidingWindow(float interactionStrength, int latticeSize, int NUMTHREAD, float T_MIN, float T_MAX, float T_STEP, long int IT);

    /**
     * @brief Perform the simulation of a phase transition using a sliding window approach.
     */
    void simulate_phase_transition() override;

    /**
     * @brief Store the results of the simulation to a file.
     */
    void store_results_to_file() const override;

protected:
 

    /**
     * @brief Flip a spin at a given lattice site during the simulation.
     *
     * @param lattice Reference to the lattice vector.
     * @param prob Probability array for the simulation.
     * @param site Index of the lattice site to flip.
     * @param magnetization Reference to the magnetization variable.
     * @param energy Reference to the energy variable.
     * @param rng_private Random number generator.
     */
    void flip(std::vector<int>& lattice, std::array<float, 2>& prob, const int& site, int& magnetization, int& energy, std::mt19937& rng_private);


    /**
     * @brief Simulate a step in the Monte Carlo simulation using a sliding window approach. The
     * number of iterations is given by NumFlipPerBlock. It is a parameter to be tuned, depending on how
     * often the lattice is translated. This affect results and performance.
     *
     * @param prob Probability array for the simulation.
     * @param lattice Reference to the lattice vector.
     * @param M Reference to the magnetization variable.
     * @param E Reference to the energy variable.
     * @param offset Offset for the starting point of each block.
     */
    void simulate_step(std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E, const int& offset) override;
    /**
     * @brief Translate the matrix in a sliding window manner.
     *
     * @param lattice Reference to the lattice vector.
     */
    void translate_matrix(std::vector<int>& lattice);

    /**
     * @brief Set the block width for the sliding window approach.
     *
     * @return The block size.
     */
    int set_block_size();

private:
    SquareLattice lattice;  /**< Square lattice used in the simulation. */
    std::vector<float> EnergyResults;  /**< Vector to store energy results. */
    std::vector<float> MagnetizationResults;  /**< Vector to store magnetization results. */
    std::vector<float>Temperatures;  /**< Vector to store temperature visited. */
    std::vector<int> ThreadStart;  /**< Vector to store Monte Carlo steps. */
    const float T_MIN;
    const float T_MAX;
    const float T_STEP;
    const int L;
    const int N;
    const long int IT;
    const int NUMTHREAD;
    int A;
    int NumSlide;
    std::uniform_real_distribution<double> dist;  /**< Uniform distribution in [0, 1). */
};

#endif // MY_PROJECT_SLIDING_WINDOW_H
