#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot-iostream.h>

void readResultsFromFile(const std::string& filePath, std::vector<float>& M_result, std::vector<float>& E_result, std::vector<float>& T_value) {
    // Open the file for reading
    std::ifstream inFile(filePath);

    // Check if the file is open
    if (!inFile.is_open()) {
        // Handle error: unable to open the file
        std::cerr << "Error: Unable to open the file for reading." << std::endl;
        return;
    }

    // Variables to store data from each line
    float energy, magnetization, temperature;

    // Skip the header line
    std::string headerLine;
    std::getline(inFile, headerLine);

    // Read data from each line and store in vectors
    while (inFile >> energy >> magnetization >> temperature) {
        E_result.emplace_back(energy);
        M_result.emplace_back(magnetization);
        T_value.emplace_back(temperature);
    }

    // Close the file
    inFile.close();
}

// Function to calculate the specific heat (C_V) for the entire range of temperatures
void calculateSpecificHeat(const std::vector<float>& E_result, const std::vector<float>& T_value, std::vector<float>& specificHeatValues) {
    // Clear the output vector
    specificHeatValues.clear();

    // Check if the input vectors have the same size
    if (E_result.size() != T_value.size()) {
        std::cerr << "Error: Input vectors have different sizes." << std::endl;
        return; // Return without modifying the output vector
    }

    // Calculate specific heat for each temperature
    for (size_t i = 0; i < T_value.size(); ++i) {
        // Calculate the mean of the energy
        float meanEnergy = 0.0;
        for (float energy : E_result) {
            meanEnergy += energy;
        }
        meanEnergy /= static_cast<float>(E_result.size());

        // Calculate the variance of the energy
        float varianceEnergy = 0.0;
        for (float energy : E_result) {
            varianceEnergy += std::pow(energy - meanEnergy, 2);
        }
        varianceEnergy /= static_cast<float>(E_result.size());

        // Calculate the specific heat (C_V) for the current temperature
        float specificHeat = varianceEnergy / (std::pow(meanEnergy, 2) * std::pow(T_value[i], 2));

        // Store the specific heat value in the output vector
        specificHeatValues.emplace_back(specificHeat);
    }
}

// Function to estimate critical temperature from magnetization trend, where derivative changes sign
float estimateCriticalTemperature(const std::vector<float>& magnetization, const std::vector<float>& temperatures) {
    // Check if input vectors have the same size
    if (magnetization.size() != temperatures.size() || magnetization.size() < 3) {
        std::cerr << "Error: Invalid input vectors." << std::endl;
        return -1.0f; // Return an invalid value
    }

    // Calculate the derivative of magnetization with respect to temperature
    std::vector<float> magnetizationDerivative(magnetization.size() - 1);
    for (size_t i = 0; i < magnetization.size() - 1; ++i) {
        magnetizationDerivative[i] = (magnetization[i + 1] - magnetization[i]) / (temperatures[i + 1] - temperatures[i]);
    }

    // Find the index where the derivative changes sign
    size_t inflectionIndex = 0;
    while (inflectionIndex < magnetizationDerivative.size() && magnetizationDerivative[inflectionIndex] >= 0) {
        ++inflectionIndex;
    }

    // Check if the inflection point is at the boundaries
    if (inflectionIndex == 0 || inflectionIndex == magnetizationDerivative.size()) {
        std::cerr << "Warning: Inflection point found at boundary. Critical temperature estimation may be less accurate." << std::endl;
    }

    // Estimate the critical temperature corresponding to the inflection point
    float criticalTemperature = temperatures[inflectionIndex];

    return criticalTemperature;
}


void plotEnergy(const std::vector<float>& E_result, const std::vector<float>& T_value) {
    Gnuplot gp;

    // Set the title of the plot
    gp << "set title 'Energy vs Temperature'\n";

    // Set the labels for the x and y axes
    gp << "set xlabel 'Temperature'\n";
    gp << "set ylabel 'Energy'\n";

    // Plot the energy values
    gp << "plot '-' with lines title 'Energy'\n";
    for (size_t i = 0; i < E_result.size(); ++i) {
        gp << T_value[i] << " " << E_result[i] << "\n";
    }
    gp << "e\n";

    // Pause to keep the plot window open
    std::cout << "Press enter to exit...\n";
    std::cin.get();
}


void plotMagnetizationOverTemperature(const std::vector<float>& M_result, const std::vector<float>& T_value) {
    // Check if input vectors have the same size
    if (M_result.size() != T_value.size()) {
        std::cerr << "Error: Input vectors must have the same size." << std::endl;
        return;
    }
    Gnuplot gp;
    gp << "plot '-' with lines title 'Magnetization'\n";
    gp << "set xlabel 'Temperature'\n";
    gp << "set ylabel 'Magnetization'\n";

    for (size_t i = 0; i < M_result.size(); ++i) {
        gp << T_value[i] << " " << M_result[i] << "\n";
    }
    gp << "e\n";

    // Wait for the user to close the plot
    std::cout << "Press enter to exit." << std::endl;
    std::cin.get();
}

void plotSpecificHeat(const std::vector<float>& specificHeat, const std::vector<float>& T_value){
    Gnuplot gp;
    gp << "set title 'Specific Heat vs Temperature'\n";
    gp << "set xlabel 'Temperature'\n";
    gp << "set ylabel 'Specific Heat'\n";

 
    // Plot the specific heat values
    gp << "plot '-' with lines title 'Specific Heat'\n";
    for (size_t i = 0; i < specificHeat.size(); ++i) {
        gp << T_value[i] << " " << specificHeat[i] << "\n";
    }
    gp << "e\n";

    // Pause to keep the plot window open
    std::cout << "Press enter to exit...\n";
    std::cin.get();
}


int main() {
    // Variables to store the file path and the result vectors
    std::string filePath;
    std::vector<float> M_result;
    std::vector<float> E_result ;
    std::vector<float> T_value ;
    std::vector<float> specificHeatValues;

    // Prompt the user to enter the file path
    std::cout << "Enter the path of the result file: ";
    std::cin >> filePath;

    // Call the function to read results from the file
    readResultsFromFile(filePath, M_result, E_result, T_value);

    // Choice menu
    int choice;
    do {
        std::cout << "Choose an option:" << std::endl;
        std::cout << "1. Plot energy" << std::endl;
        std::cout << "2. Plot magnetization over temperature" << std::endl;
        std::cout << "3. Plot specific heat" << std::endl;
        std::cout << "4. Evaluate critical temperature" << std::endl;
        std::cout << "0. Exit" << std::endl;
        std::cout << "Enter your choice: ";
        std::cin >> choice;

        switch (choice) {
            case 1:
                plotEnergy(E_result, T_value);
                break;
            case 2:
                plotMagnetizationOverTemperature( M_result,  T_value);
                break;
            case 3:
                calculateSpecificHeat( E_result,  T_value,  specificHeatValues);
                plotSpecificHeat( specificHeatValues,  T_value);
                break;
            case 4:
                float criticalTemperature = estimateCriticalTemperature( M_result,  T_value);
                std::cout << "Critical temperature: " << criticalTemperature << std::endl;
                break;
            case 0:
                std::cout << "Exiting..." << std::endl;
                break;
            default:
                std::cout << "Invalid choice. Please try again." << std::endl;
                break;
        }
    } while (choice != 0);

    return 0;
}
