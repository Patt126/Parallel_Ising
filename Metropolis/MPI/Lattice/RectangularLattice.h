#ifndef RectangularLATTICE_H
#define RectangularLATTICE_H

#include "AbstractLattice.h"
#include <iostream>
#include <memory>
#include <vector>

/**
 * @file RectangularLattice.h
 * @brief Declaration of the RectangularLattice class.
 */

/**
 * @class RectangularLattice
 * @brief Represents a 2D rectangular lattice in the Ising Monte Carlo simulation. The lattice is organized
 * in a linear fashion, with each element of the vector representing a spin. Boundary conditions are
 * menaged consider Born-Von Karman boundary conditions, namely periodic boundary conditions.
 * @details Inherited from AbstractLattice, provides concrete implementations of its functions.
 */
class RectangularLattice : public AbstractLattice {
public:
    /**
     * @brief Constructor for RectangularLattice.
     * @param interactionStrength The interaction strength for the spins interaction.
     * @param widht The widht of the rectangular lattice.
     * @param heigth the heigth of the rectangular lattice
     */
    RectangularLattice(float interactionStrength, int widht, int heigth);

    /**
     * @brief Overrides the AbstractLattice function to print the lattice.
     */
    void print_lattice() const override;

    /**
     * @brief Graph the lattice using Gnuplot.
     */
    //void graph_lattice() const;

    /**
     * @brief Overrides the AbstractLattice function to evaluate the energy of the lattice.
     * @return The energy of the lattice as a floating-point value.
     */
    float evaluate_energy() const override;

    /**
     * @brief Overrides the AbstractLattice function to initialize the lattice.
     */
    void initialize() override;

    /**
     * @brief Gets a reference to the lattice.
     * @return A reference to the lattice vector.
     */
    std::vector<int>& get_lattice();

    /**
     * @brief Overrides the AbstractLattice function to get the interaction energy.
     * @return The interaction energy of the lattice as a floating-point value.
     */
    float get_interaction_energy() const override {
        return J;
    }

    /**
     * @brief Gets the magnetization of the lattice.
     * @return The magnetization of the lattice as a floating-point value.
     */
    float get_magnetization() const {
        return M;
    }

    /**
     * @brief Increments the magnetization of the lattice.
     * @param Madd The amount to increment the magnetization by.
     */
    void increment_magnetization(int Madd) {
        M += Madd;
    }

    /**
     * @brief Increments the energy of the lattice.
     * @param Eadd The amount to increment the energy by.
     */
    void increment_energy(int Eadd) {
        E += Eadd;
    }
    float get_energy() const {
        return E;
    }
    /**
     * @brief Restores the lattice to a randomly initialized state.
     * usefull after a simulation for one temperature to start a new one.
     */
    void restore_random_lattice() {
        lattice.assign(randomLattice.begin(), randomLattice.end());
        M = M_rand;
        E = E_rand;
    }

private:
    const int H; ///< heigth of the Rectangular lattice.
    const int W; ///<width of the Rectangular lattice
    const float J; ///< Interaction strength for the lattice.
    const int N; ///< Total number of lattice sites.
    int M; ///< Magnetization of the lattice.
    float E; ///< Energy of the lattice.
    int M_rand; ///< Magnetization of the randomly initialized lattice.
    float E_rand; ///< Energy of the randomly initialized lattice.
    std::vector<int> lattice; ///< Pointer to the lattice vector.
    std::vector<int> randomLattice; ///< Pointer to the randomly initialized lattice vector.
};

#endif // RectangularLATTICE_H
