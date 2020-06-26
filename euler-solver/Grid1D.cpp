/*******************************************************************************
 * \file Grid1D.cpp
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains the implementation of the Grid1D struct.
 *
 * \date 2020-06-24
 * 
 * \copyright Copyright (c) 2020
 * 
 * \details
 ******************************************************************************/

#include <cmath>
#include <iostream>
using std::cout;
using std::endl;

#include "Grid1D.h"

// =============================================================================
double Grid1D::computeMomentumElement(size_t const &i)
{
    return velocity[i] * density[i];
}
// =============================================================================

// =============================================================================
double Grid1D::computeTotalSpecEnergyElement(size_t const &i)
{
    return siEnergy[i] + 0.5 * std::pow(velocity[i], 2);
}
// =============================================================================

// =============================================================================
std::vector<double> Grid1D::computeMomentumVec()
{
    std::vector<double> Momentum(numTotCells);

    for (size_t i = 0; i < numTotCells; i++)
    {
        Momentum[i] = computeMomentumElement(i);
    }

    return Momentum;
};
// =============================================================================

// =============================================================================
std::vector<double> Grid1D::computeTotalSpecEnergyVec()
{
    std::vector<double> TotEnergy(numTotCells);

    for (size_t i = 0; i < numTotCells; i++)
    {
        TotEnergy[i] = computeTotalSpecEnergyElement(i);
    }

    return TotEnergy;
};
// =============================================================================

// =============================================================================
void Grid1D::updateBoundaries()
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
        siEnergy[leftGhost] = siEnergy.end()[rightReal];
        siEnergy.end()[rightGhost] = siEnergy[leftReal];
    }
}
// =============================================================================

// =============================================================================
void Grid1D::saveState()
{
    if (_velocitySaveFile.is_open() &&
        _densitySaveFile.is_open() &&
        _pressureSaveFile.is_open() &&
        _siEnergySaveFile.is_open())
    {
        _velocitySaveFile << velocity[numGhostCells];
        _densitySaveFile  << density[numGhostCells];
        _pressureSaveFile << pressure[numGhostCells];
        _siEnergySaveFile << siEnergy[numGhostCells];

        for (size_t i = numGhostCells; i < (numTotCells - numGhostCells); i++)
        {
            _velocitySaveFile << "," << velocity[i];
            _densitySaveFile  << "," << density[i];
            _pressureSaveFile << "," << pressure[i];
            _siEnergySaveFile << "," << siEnergy[i];
        }
        _velocitySaveFile << endl;
        _densitySaveFile  << endl;
        _pressureSaveFile << endl;
        _siEnergySaveFile << endl;
    }
    else
    {
        cout << "File not open";
    }
}
// =============================================================================

// =============================================================================
// Constructors and destructors along with Init function
// =============================================================================
void Grid1D::init(size_t const &reals,
                  size_t const &ghosts,
                  std::string const &saveDir)
{
    numGhostCells = ghosts;
    numRealCells = reals;
    numTotCells = 2 * numGhostCells + numRealCells;

    velocity.reserve(numTotCells);
    density.reserve(numTotCells);
    pressure.reserve(numTotCells);
    siEnergy.reserve(numTotCells);

    if (saveDir != "no saving")
    {
        _velocitySaveFile.open(saveDir + "Velocity.csv");
        _densitySaveFile.open(saveDir + "Density.csv");
        _pressureSaveFile.open(saveDir + "Pressure.csv");
        _siEnergySaveFile.open(saveDir + "Internal-Energy.csv");
    }
    
}

Grid1D::Grid1D(size_t const &reals,
               size_t const &ghosts,
               std::string const &saveDir)
{
    init(reals, ghosts, saveDir);
}

Grid1D::Grid1D(size_t const &reals,
               size_t const &ghosts)
{
    init(reals, ghosts, "no saving");
}

Grid1D::~Grid1D()
{
    _velocitySaveFile.close();
    _densitySaveFile.close();
    _pressureSaveFile.close();
    _siEnergySaveFile.close();
}
// =============================================================================