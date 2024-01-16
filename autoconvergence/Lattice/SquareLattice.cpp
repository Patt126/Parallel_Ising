#include "SquareLattice.h"
//#include <gnuplot-iostream.h>


/**
 * @file SquareLattice.cpp
 * @brief Implementation of the SquareLattice class.
 */

/**
 * @brief Constructor for the SquareLattice class.
 * @param interactionStrength The interaction strength for the lattice.
 * @param latticeSize The size of the square lattice.
 */
SquareLattice::SquareLattice(float interactionStrength, int latticeSize)
    : AbstractLattice(),
      L(latticeSize),
      N(latticeSize * latticeSize),
      J(interactionStrength),
      lattice(),
      randomLattice(),
      M_rand(),
      E_rand(),
      M(),
      E()
{
    lattice = std::make_unique<std::vector<int>>(N); // Automatically reserved and initialized to 0
    randomLattice = std::make_unique<std::vector<int>>(N); // Automatically reserved and initialized to 0
    initialize();
    restore_random_lattice();
}

/**
 * @brief Overrides the AbstractLattice function to print the lattice.
 */
void SquareLattice::print_lattice() const {
    for (int i = 0; i < N; i++) {
        if (i % L == 0)
            std::cout << std::endl;
        if ((*lattice)[i] == -1)
            std::cout << "o ";
        else
            std::cout << "x ";
    }
    std::cout << std::endl;
}

/*void graph_lattice() const {
        // Open a pipe to Gnuplot
        FILE* gp = popen("gnuplot -persist", "w");

        // Check if the pipe is open
        if (!gp) {
            std::cerr << "Error: Unable to open a pipe to Gnuplot." << std::endl;
            return;
        }

        // Write Gnuplot commands to the pipe
        fprintf(gp, "set title 'Square Lattice'\n");
        fprintf(gp, "set xrange [0:%d]\n", L);
        fprintf(gp, "set yrange [0:%d]\n", L);
        fprintf(gp, "unset xtics\n");
        fprintf(gp, "unset ytics\n");
        fprintf(gp, "set size square\n");
        fprintf(gp, "plot '-' matrix with image notitle\n");

        // Write lattice data to the pipe
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                fprintf(gp, "%d ", (*lattice)[i * L + j]);
            }
            fprintf(gp, "\n");
        }

        // End the Gnuplot plot command
        fprintf(gp, "e\n");

        // Close the pipe
        fclose(gp);
    }
*/

/**
 * @brief Overrides the AbstractLattice function to evaluate the energy of the lattice.
 * @return The energy of the lattice as a floating-point value.
 */
float SquareLattice::evaluate_energy() const {
    int sum = 0;
    for (int i = 0; i < N; i++) {
        if (i >= L) // NO FIRST ROW
            sum += (*lattice)[i - L] * (*lattice)[i] * 2;
        if (i % L != 0) // NO FIRST COLUMN
            sum += (*lattice)[i - 1] * (*lattice)[i] * 2;
        if (i >= L * (L - 1)) // LAST ROW
            sum += (*lattice)[i - L * (L - 1)] * (*lattice)[i] * 2;
        if ((i + 1) % L == 0) // LAST COLUMN
            sum += (*lattice)[i - (L - 1)] * (*lattice)[i] * 2;
    }
    return -J * sum;
}

/**
 * @brief Overrides the AbstractLattice function to initialize the lattice.
 */
void SquareLattice::initialize() {
    int sum = 0;
    M_rand = 0;

    for (int i = 0; i < N; i++) {
        int k = rand() % 2;
        if (k == 0) {
            (*randomLattice)[i] = -1;
            M_rand -= 1;
        } else {
            (*randomLattice)[i] = 1;
            M_rand += 1;
        }

        if (i >= L) // NO FIRST ROW
            sum += (*randomLattice)[i - L] * (*randomLattice)[i] * 2;
        if (i % L != 0) // NO FIRST COLUMN
            sum += (*randomLattice)[i - 1] * (*randomLattice)[i] * 2;
        if (i >= L * (L - 1)) // LAST ROW
            sum += (*randomLattice)[i - L * (L - 1)] * (*randomLattice)[i] * 2;
        if ((i + 1) % L == 0) // LAST COLUMN
            sum += (*randomLattice)[i - (L - 1)] * (*randomLattice)[i] * 2;
    }

    E_rand = -J * sum;
}

/**
 * @brief Gets a reference to the lattice.
 * @return A reference to the lattice vector.
 */
std::vector<int>& SquareLattice::get_lattice() {
    return *lattice;
}
