/*!
 * \file Simulation1D.cpp
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Implementation of Simulation1D Class
 * \version 0.1
 * \date 2020-07-14
 *
 * \copyright Copyright (c) 2020
 *
 */

#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "Simulation1D.h"

// =============================================================================
void Simulation1D::_setInitialConditions(std::string const &initialConditionsKind)
{
    // Set up Sod initial conditions
    // TODO Add other kinds of initial conditions

    if (initialConditionsKind == "sod")
    {
        size_t const half  = grid.numTotCells / 2;

        // Iterate over just the real cells on the left side
        for (size_t i = grid.numGhostCells;
            i < half;
            i++)
        {
            grid.density[i]  = 1.;
            grid.momentum[i] = 0.;
            grid.energy[i] = 1. / (_gamma - 1.) ;
        }

        // Iterate over the real cells on the right side
        for (size_t i = half;
            i < (grid.numTotCells - grid.numGhostCells);
            i++)
        {
            grid.density[i]  = 0.125;  // 1/8th
            grid.momentum[i] = 0.;
            grid.energy[i] = 0.1 / (_gamma - 1);
        }
    }
    else
    {
        throw std::invalid_argument("Invalid kind of initial conditions");
    }
}
// =============================================================================

// =============================================================================
double Simulation1D::_slope(std::array<double, _arraySize> const &primitive,
                            size_t const &idx)
{
    // MC Limiter

    // Declare variables
    double outValue;
    double leftDerive, rightDerive, leftRightProduct;

    // Compute the derivatives
    leftDerive =  (primitive[idx] - primitive[idx-1]) / _deltaX;
    rightDerive = (primitive[idx + 1] - primitive[idx]) / _deltaX;
    leftRightProduct = leftDerive * rightDerive;

    // Choose what value to output
    if ((std::abs(leftDerive) < std::abs(rightDerive)) && (leftRightProduct > 0.))
    {
        outValue = leftDerive;
    }
    else if ((std::abs(leftDerive) > std::abs(rightDerive)) && (leftRightProduct > 0.))
    {
        outValue = rightDerive;
    }
    else
    {
        outValue = 0.;
    }

    return outValue;
}
// =============================================================================

// =============================================================================
void Simulation1D::_computeEigens(size_t const &idx,
                                  std::vector<double> &eigVal,
                                  std::vector<std::vector<double> > &rEigVec,
                                  std::vector<std::vector<double> > &lEigVec)
{
    // Compute the eigenvalues and vectors

    // first we have to find the speed of sound c
    double c = std::sqrt(_gamma * _pressure[idx] / _density[idx]);

    // Compute a couple of common terms
    double const cOverDensity = c/grid.density[idx];
    double const cSquared = c*c;

    // Eigenvalues are lambda^-, lambda^0, lambda^+ in that order
    eigVal[0] = _velocity[idx] - c;
    eigVal[1] = _velocity[idx];
    eigVal[2] = _velocity[idx] + c;

    // Right Eigenvectors
    rEigVec[0][0] = 1.;            rEigVec[0][1] = 1.; rEigVec[0][2] = 1.;
    rEigVec[1][0] = -cOverDensity; rEigVec[1][1] = 0.; rEigVec[1][2] = cOverDensity;
    rEigVec[2][0] = cSquared;      rEigVec[2][1] = 0.; rEigVec[2][2] = cSquared;

    // Left Eigenvectors
    lEigVec[0][0] = 0.; lEigVec[0][1] = -0.5/cOverDensity; lEigVec[0][2] =  0.5/cSquared;
    lEigVec[1][0] = 1.; lEigVec[1][1] =  0.;               lEigVec[1][2] = -1./cSquared;
    lEigVec[2][0] = 0.; lEigVec[2][1] =  0.5/cOverDensity; lEigVec[2][2] =  0.5/cSquared;
}
// =============================================================================

// =============================================================================
void Simulation1D::setPrimitives(std::string const &operation)
{
    if (operation == "reset")
    {
        // Reset current index
        _currentIndex = grid.numGhostCells;

        // Set all the values to the first few cells
        for (size_t i = 0; i < _arraySize; i++)
        {
            _density[i]  = grid.density[i];
            _velocity[i] = grid.momentum[i] / grid.density[i];
            _pressure[i] = (_gamma - 1) * (grid.energy[i]
                           - 0.5 * std::pow(grid.momentum[i], 2));
        }
    }
    else if (operation == "update")
    {
        // Shift all but the last value one to the left
        for (size_t i = 0; i < (_arraySize - 1); i++)
        {
            _density[i]  = _density[i+1];
            _velocity[i] = _velocity[i+1];
            _pressure[i] = _pressure[i+1];
        }

        // Set the final array elements
        _currentIndex++;
        _density[_arraySize]  = grid.density[_currentIndex + 2];
        _velocity[_arraySize] = grid.momentum[_currentIndex + 2] / grid.density[_currentIndex + 2];
        _pressure[_arraySize] = (_gamma - 1) * (grid.energy[_currentIndex + 2]
                       - 0.5 * std::pow(grid.momentum[_currentIndex + 2], 2));
    }
    else
    {
        throw std::invalid_argument("Invalid value for how to update primitive"
                                    "variable arrays in Simulation1D::setPrimitives");
    }
}
// =============================================================================

// =============================================================================
void Simulation1D::computeTimeStep()
{
    // Find the maximum speed in the simulation
    double vMaxTemp, vMax = 0.;

    for (size_t i = grid.numGhostCells + 1;
         i < (grid.numTotCells - grid.numGhostCells);
         i++)
    {
        double velocity = grid.momentum[i]/grid.density[i];
        double pressure = (_gamma - 1) * (grid.energy[i]
                          - 0.5 * std::pow(grid.momentum[i], 2));

        // Compute the maximum wave speed
        vMaxTemp = std::abs(velocity) +
                   std::sqrt(_gamma * pressure / grid.density[i]);
        if (vMax < vMaxTemp)
        {
            vMax = vMaxTemp;
        }
    }

    _timeStep = _cflNum * _deltaX / vMax;
}
// =============================================================================

// =============================================================================
// This function should work on the ith element and only one primitive at a time
void Simulation1D::interfaceStates(std::string const &side,
                                   std::vector<double> &leftSideOfInterface,
                                   std::vector<double> &rightSideOfInterface)
{
    // If we're computing the right (i+1/2) state then we set the index to i but
    // if we're computing the left (i-1/2) state then set the index to i-1
    size_t idx;
    if (side == "right")
    {
        idx = _center;
    }
    else if (side == "left")
    {
        idx = _center - 1;
    }
    else
    {
        throw std::invalid_argument("Invalid value for which interface to compute.");
    }

    // Some common terms that I don't want to compute multiple times
    double const dtOverDx = _timeStep / _deltaX;

    // Declare eigenvalues and eigenvectors vectors
    std::vector<double> eigVal(3);
    std::vector<std::vector<double>> rEigVec(3, std::vector<double>(3));
    std::vector<std::vector<double>> lEigVec(3, std::vector<double>(3));

    // Resize output vectors to make sure they're the right size
    leftSideOfInterface.resize(3);
    rightSideOfInterface.resize(3);

    // ===== Compute the left side of the interface ============================
    // Compute eigenvalues and eigenvectors
    _computeEigens(idx, eigVal, rEigVec, lEigVec);

    // Compute the slopes. The order is density, velocity, pressure
    std::vector<double> slopes({_slope(_density, idx),
                                _slope(_velocity, idx),
                                _slope(_pressure, idx)});

    // Compute lEigVec^nu dot slope
    std::vector<double> lEigVecDotSlope(3,0);
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            lEigVecDotSlope[i] += lEigVec[j][i] * slopes[j];
        }
    }

    // Compute the reference state
    double coef = 0.5 * (1 - dtOverDx * std::max(0., eigVal[2]));
    std::vector<double> refState({grid.density[idx] + coef*slopes[0],
                                  grid.density[idx] + coef*slopes[0],
                                  grid.density[idx] + coef*slopes[0]});

    // To find the left side of the interface state we first compute the sum in
    // the interface equation
    std::vector<double> sum(3, 0);
    for (size_t nu = 0; nu < 3; nu++) // loop over the elements in the sum
    {
        if (eigVal[nu] >= 0.)
        {
            for (size_t j = 0; j < 3; j++) // loop over primitives
            {
                sum[j] += (std::max(eigVal[2], 0.) - eigVal[nu])
                          * lEigVecDotSlope[nu]
                          * rEigVec[j][nu];
            }
        }
    }

    // Now we compute the left side of the interface
    for (size_t i = 0; i < 3; i++)
    {
        leftSideOfInterface[i] = refState[i] + 0.5 * dtOverDx * sum[i];
    }
    // ===== End computing the left side of the interface ======================

    // ===== Compute the right side of the interface ===========================

    // First increase the index to move to working on the right side of the interface
    idx++;

    // Compute eigenvalues and eigenvectors
    _computeEigens(idx, eigVal, rEigVec, lEigVec);

    // Compute the slopes. The order is density, velocity, pressure
    slopes.assign({_slope(_density, idx),
                   _slope(_velocity, idx),
                   _slope(_pressure, idx)});

    // Compute lEigVec^nu dot slope
    lEigVecDotSlope.assign(3, 0.); // zero out the vector
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            lEigVecDotSlope[i] += lEigVec[j][i] * slopes[j];
        }
    }

    // Compute the reference state
    coef = 0.5 * (1 + dtOverDx * std::min(0., eigVal[0]));
    refState.assign({grid.density[idx] - coef * slopes[0],
                     grid.density[idx] - coef * slopes[0],
                     grid.density[idx] - coef * slopes[0]});

    // To find the right side of the interface state we first compute the sum in
    // the interface equation
    sum.assign(3, 0);
    for (size_t nu = 0; nu < 3; nu++) // loop over the elements in the sum
    {
        if (eigVal[nu] <= 0.)
        {
            for (size_t j = 0; j < 3; j++) // loop over primitives
            {
                sum[j] += (std::min(eigVal[0], 0.) - eigVal[nu]) * lEigVecDotSlope[nu] * rEigVec[j][nu];
            }
        }
    }

    // Now we compute the left side of the interface
    for (size_t i = 0; i < 3; i++)
    {
        rightSideOfInterface[i] = refState[i] + 0.5 * dtOverDx * sum[i];
    }
    // ===== End computing the right side of the interface =====================
}
// =============================================================================

// =============================================================================
void Simulation1D::solveRiemann(double const &densityR,
                                double const &velocityR,
                                double const &pressureR,
                                double const &densityL,
                                double const &velocityL,
                                double const &pressureL,
                                double const &energy,
                                double const &posOverT,
                                double &energyFlux,
                                double &momentumFlux,
                                double &densityFlux)
{
    _riemannSolver.riemannMain(densityR,
                               velocityR,
                               pressureR,
                               densityL,
                               velocityL,
                               pressureL,
                               energy,
                               posOverT,
                               energyFlux,
                               momentumFlux,
                               densityFlux);
}
// =============================================================================

// =============================================================================
// Performe the conservative update
void Simulation1D::conservativeUpdate(size_t const &idxInput,
                                      double const &densityFluxL,
                                      double const &momentumFluxL,
                                      double const &energyFluxL,
                                      double const &densityFluxR,
                                      double const &momentumFluxR,
                                      double const &energyFluxR)
{
    _tempGrid.density[idxInput]  = grid.density[idxInput]
                                   + (_timeStep / _deltaX)
                                   * (densityFluxL - densityFluxR);

    _tempGrid.momentum[idxInput] = grid.momentum[idxInput]
                                   + (_timeStep / _deltaX)
                                   * (momentumFluxL - momentumFluxR);

    _tempGrid.energy[idxInput]   = grid.energy[idxInput]
                                   + (_timeStep / _deltaX)
                                   * (energyFluxL - energyFluxR);
}
// =============================================================================

// =============================================================================
void Simulation1D::updateGrid()
{
    // Copy every real element in Simulation1D::_tempGrid to Simulation1D::grid
    for (size_t i = grid.numGhostCells;
             i < (grid.numTotCells-grid.numGhostCells);
             i++)
    {
        grid.density[i]  = _tempGrid.density[i];
        grid.momentum[i] = _tempGrid.momentum[i];
        grid.energy[i]   = _tempGrid.energy[i];
    }
}
// =============================================================================

// =============================================================================
// Constructor
Simulation1D::Simulation1D(double const &physicalLength,
                           double const &gamma,
                           double const &CFL,
                           size_t const &reals,
                           size_t const &ghosts,
                           std::string const &initialConditionsKind,
                           std::string const &boundaryConditions,
                           std::string const &saveDir)

    // Start by initializing all the const member variables
    : _physLen(physicalLength),
      _cflNum(CFL),
      _deltaX(_physLen / static_cast<double>(reals)),
      _gamma(gamma),
      _riemannSolver(gamma),
      currentTime(0.0)
{
    // Now we initialize the grids.
    grid.init(reals, ghosts, saveDir, boundaryConditions);
    _tempGrid.init(reals, ghosts, "no saving", "no bc");

    // And set the initial conditions
    _setInitialConditions(initialConditionsKind);
}
// =============================================================================