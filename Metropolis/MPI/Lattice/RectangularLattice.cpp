#include "RectangularLattice.h"

RectangularLattice::RectangularLattice(float interactionStrength, int width, int heigth)
    : AbstractLattice(),
      W(width),
      H(heigth),
      N(width * heigth),
      J(interactionStrength),
      lattice(),
      randomLattice(),
      M_rand(),
      E_rand(),
      M(),
      E()
{
    randomLattice.resize(N);
    lattice.resize(N);
    initialize();
    restore_random_lattice();
}

void RectangularLattice::print_lattice() const {
    for (int i = 0; i < N; i++) {
        if (i % W == 0)
            std::cout << std::endl;
        if (lattice[i] == -1)
            std::cout << "o ";
        else
            std::cout << "x ";
    }
    std::cout << std::endl;
}

float RectangularLattice::evaluate_energy() const {
    int sum = 0;
    for (int i = 0; i < N; i++) {
        if (i >= W) // NO FIRST ROW
            sum += lattice[i - W] * lattice[i] * 2;
        if (i % W != 0) // NO FIRST COLUMN
            sum += lattice[i - 1] * lattice[i] * 2;
        if (i >= W * (H - 1)) // LAST ROW
            sum += lattice[i - W * (H - 1)] * lattice[i] * 2;
        if ((i + 1) % W == 0) // LAST COLUMN
            sum += lattice[i - (W - 1)] * lattice[i] * 2;
    }
    return -J * sum;
}

void RectangularLattice::initialize() {
    int sum = 0;
    M_rand = 0;

    for (int i = 0; i < N; i++) {
        int k = rand() % 2;
        if (k == 0) {
            randomLattice[i] = -1;
            M_rand -= 1;
        } else {
            randomLattice[i] = 1;
            M_rand += 1;
        }

        if (i >= W) // NO FIRST ROW
            sum += randomLattice[i - W] * randomLattice[i] * 2;
        if (i % W != 0) // NO FIRST COLUMN
            sum += randomLattice[i - 1] * randomLattice[i] * 2;
        if (i >= W * (H - 1)) // LAST ROW
            sum += randomLattice[i - W * (H - 1)] * randomLattice[i] * 2;
        if ((i + 1) % W == 0) // LAST COLUMN
            sum += randomLattice[i - (W - 1)] * randomLattice[i] * 2;
    }

    E_rand = -J * sum;
}

std::vector<int>& RectangularLattice::get_lattice() {
    return lattice;
}
