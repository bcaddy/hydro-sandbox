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
void Grid1D::updateBoundaries()
{
    if (boundaryConditionKind == "sod")
    {
        for (size_t j = 0; j < numGhostCells; j++)
        {
            // Compute indices
            int rightGhost = numTotCells - numGhostCells + j;
            int leftGhost  = j;

            // Update Velocity BC's
            velocity[leftGhost]  = 0.0;
            velocity[rightGhost] = 0.0;

            // Update Density BC's
            density[leftGhost]  = 1.0;
            density[rightGhost] = 0.125;

            // Update Pressure BC's
            pressure[leftGhost]  = 1.0;
            pressure[rightGhost] = 0.1;
        }

    }
    else if (boundaryConditionKind == "periodic")
    {
        // Set boundary conditions (periodic)
        for (size_t j = 0; j < numGhostCells; j++)
        {
            // Compute indices
            int rightReal  = numTotCells - (2 * numGhostCells) + j;
            int rightGhost = numTotCells - numGhostCells + j;
            int leftReal   = j + numGhostCells;
            int leftGhost  = j;

            std::cout
            << rightReal << " "
            << leftGhost << " "
            << leftReal << " "
            << rightGhost<<std::endl;


            // Update Velocity BC's
            velocity[leftGhost] = velocity[rightReal];
            velocity[rightGhost] = velocity[leftReal];

            // Update Density BC's
            density[leftGhost] = density[rightReal];
            density[rightGhost] = density[leftReal];

            // Update Pressure BC's
            pressure[leftGhost] = pressure[rightReal];
            pressure[rightGhost] = pressure[leftReal];
        }
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
    if (_velocitySaveFile.is_open() &&
        _densitySaveFile.is_open() &&
        _pressureSaveFile.is_open())
    {
        _velocitySaveFile << velocity[numGhostCells];
        _densitySaveFile  << density[numGhostCells];
        _pressureSaveFile << pressure[numGhostCells];

        for (size_t i = numGhostCells; i < (numTotCells - numGhostCells); i++)
        {
            _velocitySaveFile << "," << velocity[i];
            _densitySaveFile  << "," << density[i];
            _pressureSaveFile << "," << pressure[i];
        }
        _velocitySaveFile << std::endl;
        _densitySaveFile  << std::endl;
        _pressureSaveFile << std::endl;
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

    velocity.reserve(numTotCells);
    density.reserve(numTotCells);
    pressure.reserve(numTotCells);

    if (saveDir != "no saving")
    {
        _velocitySaveFile.open(saveDir + "Velocity.csv");
        _densitySaveFile.open(saveDir + "Density.csv");
        _pressureSaveFile.open(saveDir + "Pressure.csv");
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
    _velocitySaveFile.close();
    _densitySaveFile.close();
    _pressureSaveFile.close();
}
// =============================================================================