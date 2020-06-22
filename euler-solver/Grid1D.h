#pragma once

#include <vector>
#include <string>
#include <fstream>

struct Grid1D
{
private:
    // Output files
    std::ofstream VelocitySaveFile;
    std::ofstream DensitySaveFile;
    std::ofstream PressureSaveFile;
    std::ofstream inEnergySaveFile;

public:
    // Properties of the grid
    size_t numGhostCells;
    size_t numRealCells;
    size_t numTotCells;

    // Grid values
    std::vector<double> velocity;
    std::vector<double> density;
    std::vector<double> pressure;
    std::vector<double> inEnergy; //Specific Internal Energy

    // Get the conserved variables
    double ComputeMomentumElement(size_t const &i);
    double ComputeTotalEnergyElement(size_t const &i);
    std::vector<double> ComputeMomentumVec();
    std::vector<double> ComputeTotalEnergyVec();

    // Save the array
    void SaveState();

    // Constructors and Destructor =============================================
    Grid1D(size_t const &reals,
           size_t const &ghosts,
           std::string const &saveloc);
    Grid1D(size_t const &reals,  // This constructor is for instances that will
           size_t const &ghosts);  // never be saved
    ~Grid1D();
};


