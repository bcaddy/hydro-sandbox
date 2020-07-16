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

    size_t const half  = grid.numTotCells / 2;

    // Iterate over just the real cells on the left side
    for (size_t i = grid.numGhostCells;
         i < half;
         i++)
    {
        grid.velocity[i] = 0.;
        grid.density[i]  = 1.;
        grid.pressure[i] = 1.;
    }

    // Iterate over the real cells on the right side
    for (size_t i = half;
         i < (grid.numTotCells - grid.numGhostCells);
         i++)
    {
        grid.velocity[i] = 0;
        grid.density[i]  = 0.125;  // 1/8th
        grid.pressure[i] = 0.1;
    }
}
// =============================================================================

// =============================================================================
double Simulation1D::_slope(std::vector<double> const &primitive,
                            size_t const &i)
{
    // MC Limiter

    // Declare variables
    double outValue;
    double leftDerive, rightDerive, leftRightProduct;

    // Compute the derivatives
    leftDerive =  (primitive[i] - primitive[i-1]) / _deltaX;
    rightDerive = (primitive[i + 1] - primitive[i]) / _deltaX;
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
    double c = std::sqrt(_gamma * grid.pressure[idx] / grid.density[idx]);

    // Compute a couple of common terms
    double const cOverDensity = c/grid.density[idx];
    double const cSquared = c*c;

    // Eigenvalues are lambda^-, lambda^0, lambda^+ in that order
    eigVal[0] = grid.velocity[idx] - c;
    eigVal[1] = grid.velocity[idx];
    eigVal[2] = grid.velocity[idx] + c;

    // Right Eigenvector
    rEigVec[0][0] = 1.;            rEigVec[0][1] = 1.; rEigVec[0][2] = 1.;
    rEigVec[1][0] = -cOverDensity; rEigVec[1][1] = 0.; rEigVec[1][2] = cOverDensity;
    rEigVec[2][0] = cSquared;      rEigVec[2][1] = 0.; rEigVec[2][2] = cSquared;

    // Left Eigenvector
    lEigVec[0][0] = 0.; lEigVec[0][1] = -0.5/cOverDensity; lEigVec[0][2] =  0.5/cSquared;
    lEigVec[1][0] = 1.; lEigVec[1][1] =  0.;               lEigVec[1][2] = -1./cSquared;
    lEigVec[2][0] = 0.; lEigVec[2][1] =  0.5/cOverDensity; lEigVec[2][2] =  0.5/cSquared;
}
// =============================================================================

// =============================================================================
void Simulation1D::computeTimeStep()
{
    // I don't want to compute this in every iteration
    double const cflTimesDeltaX = _cflNum * _deltaX;

    // Set the timestep to the value determined by the first real cell
    _timeStep = cflTimesDeltaX / std::abs(grid.velocity[grid.numGhostCells]);

    // Go through the entire grid, compute the time step for each cell, and
    // choose the smallest one by setting _timeStep equal to it.
    for (size_t i = grid.numGhostCells + 1;
         i < (grid.numTotCells - grid.numGhostCells);
         i++)
    {
        double deltatTemp = cflTimesDeltaX / std::abs(grid.velocity[i]);
        if (_timeStep > deltatTemp)
        {
            _timeStep = deltatTemp;
        }
    }
    ////////////////////////////////////////////////////////////////////////////
    // Find the maximum speed in the simulation
    double vMaxTemp, vMax = 0.;

    for (size_t i = grid.numGhostCells + 1;
         i < (grid.numTotCells - grid.numGhostCells);
         i++)
    {
        // TODO change this to call Riemann sound speed function
        vMaxTemp = std::abs(grid.velocity[i]) +
                   std::sqrt(_gamma*grid.pressure[i] / grid.density[i]);
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
void Simulation1D::interfaceStates(size_t const &idxInput,
                                   std::string const &side,
                                   std::vector<double> &leftSideOfInterface,
                                   std::vector<double> &rightSideOfInterface)
{
    // If we're computing the right (i+1/2) state then we set the index to i but
    // if we're computing the left (i-1/2) state then set the index to i-1
    size_t idx;
    if (side == "right")
    {
        idx = idxInput;
    }
    else if (side == "left")
    {
        idx = idxInput - 1;
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
    std::vector<double> slopes({_slope(grid.density, idx),
                                _slope(grid.velocity, idx),
                                _slope(grid.pressure, idx)});

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
    slopes.assign({_slope(grid.density, idx),
                   _slope(grid.velocity, idx),
                   _slope(grid.pressure, idx)});

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
                                double const &posOverT,
                                double &energyFlux,
                                double &momentumFlux,
                                double &massFlux)
{
    _riemannSolver.riemannMain(densityR,
                               velocityR,
                               pressureR,
                               densityL,
                               velocityL,
                               pressureL,
                               posOverT,
                               energyFlux,
                               momentumFlux,
                               massFlux);
}
// =============================================================================

// =============================================================================
// Performe the conservative update
void Simulation1D::conservativeUpdate(size_t const &idxInput,
                                      double const &massFluxL,
                                      double const &momentumFluxL,
                                      double const &energyFluxL,
                                      double const &massFluxR,
                                      double const &momentumFluxR,
                                      double const &energyFluxR)
{
    // Density update
    _tempGrid.density[idxInput] = grid.density[idxInput]
                           + (_timeStep / _deltaX) * ( massFluxL - massFluxR );

    // Velocity update
    _tempGrid.velocity[idxInput] = (1 / _tempGrid.density[idxInput]) * (
                     (grid.density[idxInput] * grid.velocity[idxInput] )
                     + (_timeStep / _deltaX) * ( momentumFluxL - momentumFluxR )
                     );

    // Pressure update. This one is more complicated since I have to project into
    // energy space, compute the update, and project back.
    double currentEnergy = (grid.pressure[idxInput] / ((_gamma - 1)*grid.density[idxInput]))
                           + 0.5 * std::pow(grid.velocity[idxInput], 2);
    double nextEnergy = (1 / _tempGrid.density[idxInput]) * (
                        (grid.density[idxInput] * currentEnergy )
                        + (_timeStep / _deltaX) * ( energyFluxL - energyFluxR )
                        );
    _tempGrid.pressure[idxInput] = (_gamma - 1)
                                   * (_tempGrid.density[idxInput] * nextEnergy
                                   + 0.5 * _tempGrid.density[idxInput]
                                   * std::pow( _tempGrid.velocity[idxInput], 2)
                                   );
}
// =============================================================================

// =============================================================================
void Simulation1D::updateGrid()
{
    // Copy every real element in Simulation1D::_tempGrid to Simulation1D::grid
    grid.velocity = _tempGrid.velocity;
    grid.density  = _tempGrid.density;
    grid.pressure = _tempGrid.pressure;
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
    grid.init(reals, ghosts, saveDir);
    _tempGrid.init(reals, ghosts, "no saving");

    // And set the initial conditions
    _setInitialConditions(initialConditionsKind);
}
// =============================================================================