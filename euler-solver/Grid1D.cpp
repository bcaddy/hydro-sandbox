#include <cmath>
#include <iostream>
using std::cout;
using std::endl;

#include "Grid1D.h"

// =============================================================================
double Grid1D::ComputeMomentumElement(size_t const &i)
{
    return velocity[i] * density[i];
}
// =============================================================================

// =============================================================================
double Grid1D::ComputeTotalEnergyElement(size_t const &i)
{
    return inEnergy[i] + 0.5 * density[i] * std::pow(velocity[i], 2);
}
// =============================================================================

// =============================================================================
std::vector<double> Grid1D::ComputeMomentumVec()
{
    std::vector<double> Momentum(numTotCells);

    for (size_t i = 0; i < numTotCells; i++)
    {
        Momentum[i] = ComputeMomentumElement(i);
    }

    return Momentum;
};
// =============================================================================

// =============================================================================
std::vector<double> Grid1D::ComputeTotalEnergyVec()
{
    std::vector<double> TotEnergy(numTotCells);

    for (size_t i = 0; i < numTotCells; i++)
    {
        TotEnergy[i] = ComputeTotalEnergyElement(i);
    }

    return TotEnergy;
};
// =============================================================================

// =============================================================================
void Grid1D::UpdateBoundaries()
    // Set boundary conditions (periodic)
{
    for (size_t j = 0; j < numGhostCells; j++)
    {
        // Compute indices
        size_t const rightReal  = -(2 * numGhostCells - j);
        size_t const rightGhost = -(numGhostCells - j);
        size_t const leftReal   = j + numGhostCells;
        size_t const leftGhost  = j;

        // Update Velocity BC's
        velocity[leftGhost] = velocity.end()[rightReal];
        velocity.end()[rightGhost] = velocity[leftReal];

        // Update Density BC's
        density[leftGhost] = density.end()[rightReal];
        density.end()[rightGhost] = density[leftReal];

        // Update Pressure BC's
        pressure[leftGhost] = pressure.end()[rightReal];
        pressure.end()[rightGhost] = pressure[leftReal];

        // Update Internal Energy BC's
        inEnergy[leftGhost] = inEnergy.end()[rightReal];
        inEnergy.end()[rightGhost] = inEnergy[leftReal];
    }
}
// =============================================================================

// =============================================================================
void Grid1D::SaveState()
{
    if (VelocitySaveFile.is_open() &&
        DensitySaveFile.is_open() &&
        PressureSaveFile.is_open() &&
        inEnergySaveFile.is_open())
    {
        VelocitySaveFile << velocity[numGhostCells];
        DensitySaveFile  << density[numGhostCells];
        PressureSaveFile << pressure[numGhostCells];
        inEnergySaveFile << inEnergy[numGhostCells];

        for (size_t i = numGhostCells + 1; i < (numTotCells - numGhostCells); i++)
        {
            VelocitySaveFile << "," << velocity[i];
            DensitySaveFile  << "," << density[i];
            PressureSaveFile << "," << pressure[i];
            inEnergySaveFile << "," << inEnergy[i];
        }
        VelocitySaveFile << endl;
        DensitySaveFile  << endl;
        PressureSaveFile << endl;
        inEnergySaveFile << endl;
    }
    else
    {
        cout << "File not open";
    }
}
// =============================================================================

// =============================================================================
Grid1D::Grid1D(size_t const &reals,
               size_t const &ghosts,
               std::string const &saveDir)
{
    numGhostCells = ghosts;
    numRealCells = reals;
    numTotCells = 2*numGhostCells + numRealCells;

    velocity.reserve(numTotCells);
    density.reserve(numTotCells);
    pressure.reserve(numTotCells);
    inEnergy.reserve(numTotCells);

    VelocitySaveFile.open(saveDir + "Velocity.csv");
    DensitySaveFile.open(saveDir + "Density.csv");
    PressureSaveFile.open(saveDir + "Pressure.csv");
    inEnergySaveFile.open(saveDir + "Internal-Energy.csv");
}

Grid1D::Grid1D(size_t const &reals,
               size_t const &ghosts)
{
    numGhostCells = ghosts;
    numRealCells = reals;
    numTotCells = 2*numGhostCells + numRealCells;

    velocity.reserve(numTotCells);
    density.reserve(numTotCells);
    pressure.reserve(numTotCells);
    inEnergy.reserve(numTotCells);
}
Grid1D::~Grid1D()
{
    VelocitySaveFile.close();
    DensitySaveFile.close();
    PressureSaveFile.close();
    inEnergySaveFile.close();
}
// =============================================================================