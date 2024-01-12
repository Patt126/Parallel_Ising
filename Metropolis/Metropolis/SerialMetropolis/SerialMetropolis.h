/**
 * @file SerialMetropolis.h
 * @brief Header file for the SerialMetropolis class.
 */

#ifndef MY_PROJECT_SERIAL_METROPOLIS_H
#define MY_PROJECT_SERIAL_METROPOLIS_H

#include <array>
#include <vector>
#include <memory>
#include <random>
#include "../AbstractMonteCarloSimulation.h"
#include "../../Lattice/SquareLattice.h"

/**
 * @brief Class representing a Monte Carlo simulation using the serial Metropolis algorithm.
 *
 * This class inherits from AbstractMonteCarloSimulation and provides specific
 * implementations for the simulation steps using the serial Metropolis algorithm.
 */
class SerialMetropolis : public AbstractMonteCarloSimulation {
public:
    /**
     * @brief Constructor for the SerialMetropolis class.
     *
     * @param interactionStrength The interaction strength parameter.
     * @param latticeSize The size of the square lattice.
     * @param T_MIN The minimum temperature.
     * @param T_MAX The maximum temperature.
     * @param T_STEP The temperature step size.
     * @param IT The number of Monte Carlo steps.
     */
    SerialMetropolis(float interactionStrength, int latticeSize, float T_MIN, float T_MAX, float T_STEP, long int IT);

    /**
     * @brief Perform the simulation of a phase transition using the serial Metropolis algorithm.
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
     * @brief Simulate a step in the Monte Carlo simulation using the serial Metropolis algorithm.
     *
     * @param prob Probability array for the simulation.
     * @param lattice Reference to the lattice vector.
     * @param M Reference to the magnetization variable.
     * @param E Reference to the energy variable.
     * @param offset Offset for the starting point of each block.
     */
    void simulate_step(std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E, int offset = 0) override;

private:
    SquareLattice lattice;  /**< Square lattice used in the simulation. */
    std::unique_ptr<std::vector<int>> RandVect;  /**< Random vector of size N. */
    std::unique_ptr<std::vector<float>> EnergyResults;  /**< Vector to store energy results. */
    std::unique_ptr<std::vector<float>> MagnetizationResults;  /**< Vector to store magnetization results. */
    std::unique_ptr<std::vector<float>> Temperatures;  /**< Vector to store temperature visited. */
    std::unique_ptr<std::vector<int>> ThreadStart;  /**< Vector to store Monte Carlo steps. */
    float T_MIN;
    float T_MAX;
    float T_STEP;
    const int L;
    const int N;
    const long int IT;
    std::mt19937 rng;  /**< Mersenne Twister 19937 generator. */
    std::uniform_real_distribution<double> dist;  /**< Uniform distribution in [0, 1). */
};

#endif // MY_PROJECT_SERIAL_METROPOLIS_H
