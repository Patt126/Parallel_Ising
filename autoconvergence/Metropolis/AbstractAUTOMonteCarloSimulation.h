#ifndef ABSTRACT_MONTE_CARLO_SIMULATION_H
#define ABSTRACT_MONTE_CARLO_SIMULATION_H

#include <array>
#include <vector>

/**
 * @brief Abstract base class for Monte Carlo simulations.
 *
 * This class defines the interface for Monte Carlo simulations and
 * declares pure virtual functions that must be implemented by derived classes.
 */
class AbstractAUTOMonteCarloSimulation {
public:
    /**
     * @brief Perform the simulation of a phase transition. performimg 
     * the simulation for a range of temperatures.
     * 
     * Pure virtual function to be implemented by derived classes.
     */
    virtual void simulate_phase_transition() = 0;

protected:

    /**
     * @brief Simulate a step in the Monte Carlo simulation. 
     * A step correspong in general to multiple Monte Carlo steps.
     *
     * @param prob Probability array governed by Boltzman factor for
     * detailed balance of flipping probability.
     * @param lattice Reference to the lattice of spin for the simulation.
     * @param M Reference to the magnetization variable.
     * @param E Reference to the energy variable.
     * @param offset When working in parallel the offset indicate the starting
     * point in which each block starts.
     *
     * Pure virtual function to be implemented by derived classes.
     */
    virtual int simulate_step(std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E, const int& offset ) = 0;

    /**
     * @brief Store the results of Energy and Magnetization of the simulation to a file.
     *
     * Pure virtual function to be implemented by derived classes.
     */
    virtual void store_results_to_file() const = 0;

};

#endif // ABSTRACT_MONTE_CARLO_SIMULATION_H

