/**
 * @brief Simulate the phase transition using a convergence criterion based on the fluctuation of magnetization.
 *
 * This method performs a Monte Carlo simulation over a range of temperatures, utilizing a convergence criterion
 * based on the standard deviation of the last few magnetization values. Simulation results are recorded for analysis.
 */
void AutostopMontecarlo::simulate_phase_transition() {
    int deltaE;
    int deltaM;
    int step;
    float m = 0;
    float T = T_MIN;

    float error = 0;
    std::array<float, 2> prob;

    // Parameters for convergence check
    int windowSize = 10;  // Number of last magnetization values to consider
    float convergenceThreshold = 1e-5;

    // Loop over temperature range
    while (T < T_MAX) {
        // Calculate transition probabilities
        prob[0] = std::exp(-4 * lattice.getInteractionEnergy() / T);
        prob[1] = std::exp(-8 * lattice.getInteractionEnergy() / T);

        step = 0;
        m = static_cast<float>(lattice.getMagnetization()) / N;

        // Continue simulation until convergence criterion is met or a maximum number of steps is reached
        bool simulationContinue = true;
        while (step < maxSimulationSteps && simulationContinue) {
            // Perform a Monte Carlo step
            deltaE = 0;
            deltaM = 0;
            create_rand_vector();
            simulate_step(prob, lattice.getLattice(), deltaM, deltaE);
            lattice.incrementMagnetization(deltaM);
            lattice.incrementEnergy(deltaE);
            m = static_cast<float>(lattice.getMagnetization()) / N;

            // Update error and step count
            error = std::abs(m - mexact);
            step++;

            // Check for convergence based on fluctuation of last magnetization values
            if (step >= windowSize) {
                float lastMagnetizationsSD = calculateLastMagnetizationsSD(windowSize);
                if (lastMagnetizationsSD < convergenceThreshold) {
                    simulationContinue = false;  // Set flag to end simulation
                }
            }
        }

        // Record simulation results for the current temperature
        temperatures->emplace_back(T);
        T += T_STEP;
        monteCarloStepsResults->emplace_back(step);
        energyResults->emplace_back(1);  // Placeholder value, update as needed
        magnetizationResults->emplace_back(abs(m));
        step = 0;

        // Restore the lattice to its initial state for the next temperature
        lattice.restoreRandomLattice();
    }
}
