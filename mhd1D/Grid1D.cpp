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
    if (boundaryConditionKind == "periodic")
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

            // Update Energy BC's
            energy[leftGhost]  = energy[rightReal];
            energy[rightGhost] = energy[leftReal];

            for (size_t i = 0; i < 3; i++)
            {
                // Update Momentum BC's
                momentum[leftGhost][i]  = momentum[rightReal][i];
                momentum[rightGhost][i] = momentum[leftReal][i];

                // Update Magnetic field BC's
                magnetic[leftGhost][i]  = magnetic[rightReal][i];
                magnetic[rightGhost][i] = magnetic[leftReal][i];
            }
        }
    }
    else if (boundaryConditionKind == "pass")
    {
        ; // Pass
    }
    else
    {
        std::cout << std::endl <<  "Given boundary Condition type: " <<
            boundaryConditionKind << std::endl;
        throw std::invalid_argument("Invalid kind of boundary conditions");
    }

}
// =============================================================================

// =============================================================================
void Grid1D::saveState()
{
    if (_densitySaveFile.is_open() and
        _momentumXSaveFile.is_open() and
        _momentumYSaveFile.is_open() and
        _momentumZSaveFile.is_open() and
        _magneticXSaveFile.is_open() and
        _magneticYSaveFile.is_open() and
        _magneticZSaveFile.is_open() and
        _energySaveFile.is_open())
    {

        _densitySaveFile   << density[numGhostCells];
        _momentumXSaveFile << momentum[numGhostCells][0];
        _momentumYSaveFile << momentum[numGhostCells][1];
        _momentumZSaveFile << momentum[numGhostCells][2];
        _magneticXSaveFile << magnetic[numGhostCells][0];
        _magneticYSaveFile << magnetic[numGhostCells][1];
        _magneticZSaveFile << magnetic[numGhostCells][2];
        _energySaveFile    << energy[numGhostCells];

        for (size_t i = numGhostCells + 1;
             i < (numTotCells - numGhostCells);
             i++)
        {

            _densitySaveFile   << "," << density[i];
            _momentumXSaveFile << "," << momentum[i][0];
            _momentumYSaveFile << "," << momentum[i][1];
            _momentumZSaveFile << "," << momentum[i][2];
            _magneticXSaveFile << "," << magnetic[i][0];
            _magneticYSaveFile << "," << magnetic[i][1];
            _magneticZSaveFile << "," << magnetic[i][2];
            _energySaveFile    << "," << energy[i];
        }

        _densitySaveFile   << std::endl;
        _momentumXSaveFile << std::endl;
        _momentumYSaveFile << std::endl;
        _momentumZSaveFile << std::endl;
        _magneticXSaveFile << std::endl;
        _magneticYSaveFile << std::endl;
        _magneticZSaveFile << std::endl;
        _energySaveFile    << std::endl;
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
Grid1D::Grid1D(size_t const &reals,
               size_t const &ghosts,
               std::string const &boundaryConditions,
               std::string const &saveDir)
    :
    numGhostCells(ghosts),
    numRealCells(reals),
    numTotCells(2*numGhostCells + numRealCells),
    boundaryConditionKind(boundaryConditions),
    density(numTotCells),
    momentum(numTotCells, std::vector<double> (3, 0)),
    magnetic(numTotCells, std::vector<double> (3, 0)),
    energy(numTotCells)
{
    // Open files if saving is enabled
    if (saveDir != "0")
    {
        _densitySaveFile.open(saveDir + "Density.csv");

        _momentumXSaveFile.open(saveDir + "MomentumX.csv");
        _momentumYSaveFile.open(saveDir + "MomentumY.csv");
        _momentumZSaveFile.open(saveDir + "MomentumZ.csv");

        _magneticXSaveFile.open(saveDir + "MagneticX.csv");
        _magneticYSaveFile.open(saveDir + "MagneticY.csv");
        _magneticZSaveFile.open(saveDir + "MagneticZ.csv");

        _energySaveFile.open(saveDir + "Energy.csv");
    }
}

Grid1D::~Grid1D()
{
    // Close the file if it's open
    if (_densitySaveFile.is_open())   _densitySaveFile.close();
    if (_momentumXSaveFile.is_open()) _momentumXSaveFile.close();
    if (_momentumYSaveFile.is_open()) _momentumYSaveFile.close();
    if (_momentumZSaveFile.is_open()) _momentumZSaveFile.close();
    if (_magneticXSaveFile.is_open()) _magneticXSaveFile.close();
    if (_magneticYSaveFile.is_open()) _magneticYSaveFile.close();
    if (_magneticZSaveFile.is_open()) _magneticZSaveFile.close();
    if (_energySaveFile.is_open())    _energySaveFile.close();

}
// =============================================================================