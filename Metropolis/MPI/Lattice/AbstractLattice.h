// AbstractLattice.h

#ifndef ABSTRACTLATTICE_H
#define ABSTRACTLATTICE_H

/**
 * @file AbstractLattice.h
 * @brief Declaration of the AbstractLattice class.
 */

/**
 * @class AbstractLattice
 * @brief Abstract base class for representing a lattice in the Ising Monte Carlo simulation.
 * The lattice is represented as a vector of integers, where each integer represents a spin.
 */
class AbstractLattice {
public:
    /**
     * @brief Pure virtual function to initialize a random lattice
     * corrisponding to a situation above critical temperature.
     */
    virtual void initialize() = 0;

    /**
     * @brief Pure virtual function to evaluate the energy of the lattice.
     * @return The energy of the lattice as a floating-point value.
     */
    virtual float evaluate_energy() const = 0;

    /**
     * @brief Pure virtual function to print the lattice.
     * @details This function is intended for debugging or visualization purposes.
     */
    virtual void print_lattice() const = 0;

    /**
     * @brief Pure virtual function to get the interaction energy of the lattice's spins.
     * @return The interaction energy of the lattice's spins as a floating-point value.
     */
    virtual float get_interaction_energy() const = 0;

    /**
     * @brief Virtual destructor to ensure proper cleanup in derived classes.
     */
    virtual ~AbstractLattice() = default;
};

#endif // ABSTRACTLATTICE_H
