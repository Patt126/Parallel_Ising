#ifndef MY_PROJECT_MONTE_CARLO_SIMULATION_H
#define MY_PROJECT_MONTE_CARLO_SIMULATION_H

#include <array>
#include <vector>
#include <cmath>
#include <memory>
#include <random>
#include "AbstractAUTOMonteCarloSimulation.h"
#include "SquareLattice.h"

/**
 * @brief Class representing a Monte Carlo simulation with an automatic stopping criterion.
 *
 * This class inherits from AbstractAUTOMonteCarloSimulation and provides specific
 * implementations for the simulation steps with an automatic stopping criterion
 * based on exact solution or alternatively on the fluctuation of magnetization. The simulation results are stored for analysis.
 */
class AutostopMontecarlo : public AbstractAUTOMonteCarloSimulation {
public:
    /**
     * @brief Constructor for the AutostopMontecarlo class.
     *
     * @param interactionStrength The interaction strength parameter.
     * @param latticeSize The size of the square lattice.
     * @param Tolerance Tolerance value for the automatic stopping criterion.
     * @param T_MIN The minimum temperature.
     * @param T_MAX The maximum temperature.
     * @param T_STEP The temperature step size.
     */
    AutostopMontecarlo(float interactionStrength, int latticeSize, float T_MIN, float T_MAX, float T_STEP,float tolerance);

    /**
     * @brief Perform the simulation of a phase transition with an automatic stopping criterion.
     */
    void simulate_phase_transition() override;

    /**
     * @brief Store the results of the simulation to a file.
     */
    void store_results_to_file() const override;

protected:
  
    /**
     * @brief Calculate the exact magnetization at a given temperature from Onsager formula.
     *
     * @param T The temperature.
     * @return The exact magnetization at the given temperature.
     */
    float m_exact(float T) const {
        return std::pow((1.0 - std::pow(std::sinh(2 * lattice.get_interaction_energy() / T), -4)), 1.0 / 8.0);
    }

    /**
     * @brief Flip a spin at a given lattice site during the simulation.
     *
     * @param lattice Reference to the lattice vector.
     * @param prob Probability array for the simulation.
     * @param site Index of the lattice site to flip.
     * @param M Reference to the magnetization variable.
     * @param E Reference to the energy variable.
     */
    void flip(std::vector<int>& lattice, std::array<float, 2>& prob, int site, int& M, int& E) ;

    /**
     * @brief Simulate a step for the entire lattice in the Monte Carlo simulation.
     *
     * @param prob Probability array for the simulation.
     * @param lattice Reference to the lattice vector.
     * @param M Reference to the magnetization variable.
     * @param E Reference to the energy variable.
     */
    int simulate_step(std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E,const int& offset = 0 ) override;

private:

    SquareLattice lattice;  // Use SquareLattice as a private member
    std::vector<int> randVect;               // Random vector of size N.
    std::vector<float> energyResults;        // Vector to store energy results.
    std::vector<float> magnetizationResults; // Vector to store magnetization results.
    std::vector<int> monteCarloStepsResults; // Vector to store number of Monte Carlo steps performed.
    std::vector<float> temperatures;         // Vector to store temperatures visited.
    float T_MIN;
    float T_MAX;
    float tolerance; // Tolerance value.
    float T_STEP;
    const int L;
    const int N;
    std::mt19937 rng; // Random number generator.
    std::uniform_real_distribution<double> dist; // Uniform distribution.
};

#endif // MY_PROJECT_MONTE_CARLO_SIMULATION_H
