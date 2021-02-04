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
#include "mhdUtilities.h"

using namespace mhdUtilities;

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
    else if (boundaryConditionKind == "bwShockTube")
    {
        // Choose Values
        double denL  = 1.0,
        presL = 1.0,
        denR  = 0.125,
        presR = 0.1;

        std::vector<double> velL{0.0, 0.0, 0.0},
                            bL{0.75, 1.0, 0.0},
                            velR{0.0, 0.0, 0.0},
                            bR{0.75, -1.0, 0.0};

        for (size_t iL = 0; iL < numGhostCells; iL++)
        {
            // Compute right most index
            size_t iR = numTotCells - 1 - iL;

            // Set the left ghost cells
            density[iL]     = denL;
            momentum[iL][0] = computeMomentum(velL[0], denL);
            momentum[iL][1] = computeMomentum(velL[1], denL);
            momentum[iL][2] = computeMomentum(velL[2], denL);
            magnetic[iL][0] = bL[0];
            magnetic[iL][1] = bL[1];
            magnetic[iL][2] = bL[2];
            energy[iL]      = computeEnergy(presL, denL, velL, bL, gamma);

            density[iR]     = denR;
            momentum[iR][0] = computeMomentum(velL[0], denR);
            momentum[iR][1] = computeMomentum(velL[1], denR);
            momentum[iR][2] = computeMomentum(velL[2], denR);
            magnetic[iR][0] = bR[0];
            magnetic[iR][1] = bR[1];
            magnetic[iR][2] = bR[2];
            energy[iR]      = computeEnergy(presR, denR, velR, bR, gamma);
        }
    }
    else if (boundaryConditionKind == "dwShockTube")
    {
        // Choose Values
        double coef = 1. / std::sqrt(4. * M_PI);
        double denL  = 1.08,
               presL = 0.95,
               denR  = 1.,
               presR = 1.;

        std::vector<double> velL{1.2, 0.01, 0.5},
                            bL{4. * coef, 3.6 * coef, 2.0 * coef},
                            velR{0.0, 0.0, 0.0},
                            bR{4.0 * coef, 4.0 * coef, 2.0 * coef};

        for (size_t iL = 0; iL < numGhostCells; iL++)
        {
            // Compute right most index
            size_t iR = numTotCells - 1 - iL;

            // Set the left ghost cells
            density[iL]     = denL;
            momentum[iL][0] = computeMomentum(velL[0], denL);
            momentum[iL][1] = computeMomentum(velL[1], denL);
            momentum[iL][2] = computeMomentum(velL[2], denL);
            magnetic[iL][0] = bL[0];
            magnetic[iL][1] = bL[1];
            magnetic[iL][2] = bL[2];
            energy[iL]      = computeEnergy(presL, denL, velL, bL, gamma);

            density[iR]     = denR;
            momentum[iR][0] = computeMomentum(velL[0], denR);
            momentum[iR][1] = computeMomentum(velL[1], denR);
            momentum[iR][2] = computeMomentum(velL[2], denR);
            magnetic[iR][0] = bR[0];
            magnetic[iR][1] = bR[1];
            magnetic[iR][2] = bR[2];
            energy[iR]      = computeEnergy(presR, denR, velR, bR, gamma);
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
    density(numTotCells),
    momentum(numTotCells, std::vector<double> (3, 0)),
    magnetic(numTotCells, std::vector<double> (3, 0)),
    energy(numTotCells),
    boundaryConditionKind(boundaryConditions)
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