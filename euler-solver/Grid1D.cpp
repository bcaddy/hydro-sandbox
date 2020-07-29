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


#include <iostream>
#include <stdexcept>

#include "Grid1D.h"

// =============================================================================
void Grid1D::updateBoundaries(double const &gamma)
{
    if (boundaryConditionKind == "sod")
    {
        for (size_t j = 0; j < numGhostCells; j++)
        {
            // Compute indices
            int rightGhost = numTotCells - numGhostCells + j;
            int leftGhost  = j;

            // Update Density BC's
            density[leftGhost]  = 1.0;
            density[rightGhost] = 0.125;

            // Update Momentum BC's
            momentum[leftGhost]  = 0.0;
            momentum[rightGhost] = 0.0;

            // Update Energy BC's
            energy[leftGhost]  = 1. / (gamma - 1.);
            energy[rightGhost] = 0.1 / (gamma - 1);
        }

    }
    else if (boundaryConditionKind == "periodic")
    {
        // Set boundary conditions (periodic)
        for (size_t j = 0; j < numGhostCells; j++)
        {
            // Compute indices
            size_t rightReal  = numTotCells - (2 * numGhostCells) + j;
            size_t rightGhost = numTotCells - numGhostCells + j;
            size_t leftReal   = j + numGhostCells;
            size_t leftGhost  = j;

            // Update Density BC's
            density[leftGhost]  = density[rightReal];
            density[rightGhost] = density[leftReal];

            // Update Momentum BC's
            momentum[leftGhost]  = momentum[rightReal];
            momentum[rightGhost] = momentum[leftReal];

            // Update Energy BC's
            energy[leftGhost]  = energy[rightReal];
            energy[rightGhost] = energy[leftReal];
        }
    }
    else if (boundaryConditionKind == "pass")
    {
        ; // Pass
    }
    else
    {
        throw std::invalid_argument("Invalid kind of boundary conditions");
    }

}
// =============================================================================

// =============================================================================
void Grid1D::saveState()
{
    if (_densitySaveFile.is_open() and
        _momentumSaveFile.is_open() and
        _energySaveFile.is_open())
    {
        _densitySaveFile  << density[numGhostCells];
        _momentumSaveFile << momentum[numGhostCells];
        _energySaveFile << energy[numGhostCells];

        for (size_t i = numGhostCells + 1;
             i < (numTotCells - numGhostCells);
             i++)
        {
            _densitySaveFile  << "," << density[i];
            _momentumSaveFile << "," << momentum[i];
            _energySaveFile << "," << energy[i];
        }
        _densitySaveFile  << std::endl;
        _momentumSaveFile << std::endl;
        _energySaveFile << std::endl;
    }
    else
    {
        std::cout << "File not open";
    }
}
// =============================================================================

// =============================================================================
// Constructors and destructors along with Init function
// =============================================================================
void Grid1D::init(size_t const &reals,
                  size_t const &ghosts,
                  std::string const &saveDir,
                  std::string const &boundaryConditions)
{
    numGhostCells = ghosts;
    numRealCells = reals;
    numTotCells = 2 * numGhostCells + numRealCells;
    boundaryConditionKind = boundaryConditions;

    density.resize(numTotCells);
    momentum.resize(numTotCells);
    energy.resize(numTotCells);

    if (saveDir != "no saving")
    {
        _densitySaveFile.open(saveDir + "Density.csv");
        _momentumSaveFile.open(saveDir + "Momentum.csv");
        _energySaveFile.open(saveDir + "Energy.csv");
    }

}

Grid1D::Grid1D(size_t const &reals,
               size_t const &ghosts,
               std::string const &saveDir,
               std::string const &boundaryConditions)
{
    init(reals, ghosts, saveDir, boundaryConditions);
}

Grid1D::Grid1D(size_t const &reals,
               size_t const &ghosts)
{
    init(reals, ghosts, "no bc", "no saving");
}

Grid1D::~Grid1D()
{
    // Close the file if it's open
    if (_densitySaveFile.is_open() and
        _momentumSaveFile.is_open() and
        _energySaveFile.is_open())
    {
        _densitySaveFile.close();
        _momentumSaveFile.close();
        _energySaveFile.close();
    }
}
// =============================================================================