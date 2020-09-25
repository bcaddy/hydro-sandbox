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
#include <numeric>

#include "Simulation1D.h"
#include "HydroUtilities.h"

using namespace HydroUtilities;

// =============================================================================
void Simulation1D::_setInitialConditions(std::string const &initialConditionsKind)
{
    if (initialConditionsKind == "sod")
    {
        size_t const half  = grid.numTotCells / 2;

        // Iterate over just the real cells on the left side
        for (size_t i = grid.numGhostCells;
            i < half;
            i++)
        {
            grid.density[i]  = 1.;
            grid.momentum[i] = computeMomentum(0.0, 1.0);
            grid.energy[i]   = computeEnergy(1.0, 1.0, 0.0, _gamma);
        }

        // Iterate over the real cells on the right side
        for (size_t i = half;
            i < (grid.numTotCells - grid.numGhostCells);
            i++)
        {
            grid.density[i]  = 0.1;
            grid.momentum[i] = computeMomentum(0.0, 0.1);;
            grid.energy[i]   = computeEnergy(0.1, 0.1, 0.0, _gamma);
        }
    }
    else if (initialConditionsKind == "indexCheck")
    {
        for (size_t i = grid.numGhostCells;
             i < (grid.numTotCells-grid.numGhostCells);
             i++)
        {
            double dIdx = static_cast<double>(i) + 0.001;
            grid.density[i]  = dIdx;
            grid.momentum[i] = computeMomentum(dIdx, dIdx);
            grid.energy[i]   = computeEnergy(dIdx, dIdx, dIdx, _gamma);
        }

    }
    else if (initialConditionsKind == "advectionStep")
    {
        double pressure = 1.0, velocity = 1.0;
        double denLow = 1.0, denHigh = 2.0;

        for (size_t i = grid.numGhostCells;
             i < (grid.numTotCells-grid.numGhostCells);
             i++)
        {
            if ( (grid.numTotCells/4 < i) and (i < grid.numTotCells/2) )
            // if ((i > 3*grid.numRealCells/4) )
            {
                // We're in the high pressure region
                grid.density[i] = denHigh;
            }
            else
            {
                // We're in the low pressure region
                grid.density[i] = denLow;
            }

            // Set the other conserved variables
            grid.momentum[i] = computeMomentum(velocity, grid.density[i]);
            grid.energy[i]   = computeEnergy(pressure, grid.density[i], velocity, _gamma);
        }
    }
    else if (initialConditionsKind == "advectionGauss")
    {
        double const pressure = 1.0, velocity = 1.0;
        double const denLow = 10.0, amplitude = 1.0;

        for (size_t i = grid.numGhostCells;
             i < (grid.numTotCells-grid.numGhostCells);
             i++)
        {

            double x = 4. * (static_cast<double>(i) / static_cast<double>(grid.numTotCells)) - 2.;

            grid.density[i] = denLow + amplitude * std::exp( -std::pow(x,2));

            // Set the other conserved variables
            grid.momentum[i] = computeMomentum(velocity, grid.density[i]);
            grid.energy[i]   = computeEnergy(pressure, grid.density[i], velocity, _gamma);
        }
    }
    else
    {
        throw std::invalid_argument("Invalid kind of initial conditions");
    }
}
// =============================================================================

// =============================================================================
double Simulation1D::_slope(std::vector<double> const &primitive,
                            size_t const &idx)
{
    if (limiterKind == "zeroSlope")  // Always return zero slope
    {
        return 0.0;
    }
    else if (limiterKind == "centerDiff")  // Use centered diff
    {
        return 0.5 * (primitive[idx+1] - primitive[idx-1]);
    }
    else if (limiterKind == "minMod") // Use MinMod limiter
    {
        // Declare variables
        double outValue;
        double leftDerive, rightDerive, leftRight;

        // Compute the derivatives
        leftDerive = (primitive[idx] - primitive[idx-1]);
        rightDerive = (primitive[idx+1] - primitive[idx]);
        leftRight = leftDerive * rightDerive;

        // Choose what value to output
        if ((std::abs(leftDerive) < std::abs(rightDerive)) && (leftRight > 0.))
        {
            outValue = leftDerive;
        }
        else if ((std::abs(leftDerive) > std::abs(rightDerive)) && (leftRight > 0.))
        {
            outValue = rightDerive;
        }
        else
        {
            outValue = 0.0;
        }

        return outValue;
    }
    else if (limiterKind == "MC")
    {
        double xi = (primitive[idx+1] - primitive[idx]) * (primitive[idx] - primitive[idx-1]);

        if (xi > 0)
        {
            double centerDif, forwardDif, backwardDif, sign;
            centerDif   = 0.5 * std::abs(primitive[idx+1] - primitive[idx-1]);
            forwardDif  = 2.0 * std::abs(primitive[idx+1] - primitive[idx]);
            backwardDif = 2.0 * std::abs(primitive[idx]   - primitive[idx-1]);

            sign = std::copysign(1.0, primitive[idx+1]-primitive[idx-1]);// equivalent to sign(a-b)

            return sign * std::min({centerDif, forwardDif, backwardDif});
        }
        else
        {
            return 0.;
        }
    }
    else
    {
        throw std::invalid_argument("Invalid kind of slope limiter chosen.");
    }

}
// =============================================================================

// =============================================================================
void Simulation1D::_piecewiseLinearReconstruction()
{
    // Some common terms that I don't want to compute multiple times
    double const dtOverDx = _timeStep / _deltaX;

    // Compute all the primitive values
    std::vector<double> velocity(grid.numTotCells), pressure(grid.numTotCells);
    for (size_t i = 0;
         i < grid.numTotCells;
         i++)
    {
        velocity[i] = computeVelocity(grid.momentum[i],
                                      grid.density[i]);
        pressure[i] = computePressure(grid.energy[i],
                                      grid.density[i],
                                      velocity[i],
                                      _gamma);
    }

    // Loop over the entire grid and reconstruct the state within each cell then
    // compute the i-1/2,R and i+1/2,L states
    for (size_t i = (grid.numGhostCells - 1);
         i < (grid.numTotCells - (grid.numGhostCells - 1) );
         i++)
    {
        // ===== Compute the eigenvalues and eigenvectors ======================
        // Compute sound speed
        double c = soundSpeed(pressure[i], grid.density[i], _gamma);

        // Compute eigenvalues
        std::vector<double> eigVal{velocity[i] - c,
                                   velocity[i],
                                   velocity[i] + c};

        // Compute Left eigenvectors
        double const cSquared = c*c;
        std::vector<std::vector<double>> eigVecL
                        {{0.0,  -0.5 * grid.density[i] / c,   0.5 / cSquared },
                         {1.0,                         0.0,  -1.0 / cSquared },
                         {0.0,   0.5 * grid.density[i] / c,   0.5 / cSquared }};

        // Compute Right eigenvectors
        std::vector<std::vector<double>> eigVecR
                        {{                 1.0,  1.0,                  1.0},
                         {-c / grid.density[i],  0.0,  c / grid.density[i]},
                         {            cSquared,  0.0,             cSquared}};
        // ===== Finished computing the eigenvalues and eigenvectors ===========

        // ===== Compute the reference states ==================================
        // Slopes
        std::vector<double> slopes {_slope(grid.density, i),
                                    _slope(velocity,     i),
                                    _slope(pressure,     i)};

        // Coefficients
        double coefL = 0.5 * (1.0 - dtOverDx * std::max(eigVal[2], 0.0));
        double coefR = 0.5 * (1.0 + dtOverDx * std::min(eigVal[0], 0.0));

        // Compute the reference states
        std::vector<double> refStateL(3), refStateR(3);

        // Left
        refStateL[0] = grid.density[i] + coefL * slopes[0];
        refStateL[1] = velocity[i]     + coefL * slopes[1];
        refStateL[2] = pressure[i]     + coefL * slopes[2];

        // Right
        refStateR[0] = grid.density[i] - coefR * slopes[0];
        refStateR[1] = velocity[i]     - coefR * slopes[1];
        refStateR[2] = pressure[i]     - coefR * slopes[2];
        // ===== Finished computing the reference states =======================

        // ===== Compute the Vhat functions ====================================
        std::vector<double> betaL(3), betaR(3);
        for (size_t j = 0; j < 3; j++)
        {
            double sum = std::inner_product(eigVecL[j].begin(),
                                            eigVecL[j].end(),
                                            slopes.begin(),
                                            0);

            betaL[j] = (dtOverDx/4.) * (eigVal[2] - eigVal[j])
                       * (std::copysign(1.0, eigVal[j]) + 1.0) * sum;
            betaR[j] = (dtOverDx/4.) * (eigVal[0] - eigVal[j])
                       * (1.0 - std::copysign(1.0, eigVal[j])) * sum;
        }
        // ===== Finished computing the Vhat functions =========================

        // ===== Compute the i+1/2,L and i-1/2,R interface states ==============
        // i+1/2,L state
        _densityInterfaceL[i+1]  = refStateL[0] + std::inner_product(betaL.begin(), betaL.end(), eigVecR[0].begin(), 0);
        _velocityInterfaceL[i+1] = refStateL[1] + std::inner_product(betaL.begin(), betaL.end(), eigVecR[1].begin(), 0);
        _pressureInterfaceL[i+1] = refStateL[2] + std::inner_product(betaL.begin(), betaL.end(), eigVecR[2].begin(), 0);

        // i-1/2,R state
        _densityInterfaceR[i]  = refStateR[0] + std::inner_product(betaR.begin(), betaR.end(), eigVecR[0].begin(), 0);
        _velocityInterfaceR[i] = refStateR[1] + std::inner_product(betaR.begin(), betaR.end(), eigVecR[1].begin(), 0);
        _pressureInterfaceR[i] = refStateR[2] + std::inner_product(betaR.begin(), betaR.end(), eigVecR[2].begin(), 0);
        // ===== Finished computing the i+1/2,L and i-1/2,R interface states ===
    }
}
// =============================================================================

// =============================================================================
// Implement the piecewise constant reconstruction of the interface states
void Simulation1D::_piecewiseConstantReconstruction()
{
    // Loop through every element of the grid and compute the interface states.
    // Note that we have to go 1 element farther than usual since we need the
    // interface states on both sides of the real part of the grid
    for (size_t i = grid.numGhostCells;
         i < (grid.numTotCells - grid.numGhostCells + 1);
         i++)
    {
        // Compute the left interfaces
        _densityInterfaceL[i] = grid.density[i-1];
        _velocityInterfaceL[i] = computeVelocity(grid.momentum[i-1],
                                                 grid.density[i-1]);
        _pressureInterfaceL[i] = computePressure(grid.energy[i-1],
                                                 grid.density[i-1],
                                                 _velocityInterfaceL[i],
                                                 _gamma);

        // Compute the right interfaces
        _densityInterfaceR[i] = grid.density[i];
        _velocityInterfaceR[i] = computeVelocity(grid.momentum[i],
                                                 grid.density[i]);
        _pressureInterfaceR[i] = computePressure(grid.energy[i],
                                                 grid.density[i],
                                                 _velocityInterfaceR[i],
                                                 _gamma);
    }
}
// =============================================================================

// =============================================================================
void Simulation1D::computeTimeStep()
{
    // Find the maximum speed in the simulation
    double vMaxTemp, vMax = 0.;

    for (size_t i = grid.numGhostCells;
         i < (grid.numTotCells - grid.numGhostCells);
         i++)
    {
        double velocity = computeVelocity(grid.momentum[i], grid.density[i]);
        double pressure = computePressure(grid.energy[i], grid.density[i], velocity, _gamma);

        // Compute the maximum wave speed
        vMaxTemp = std::abs(velocity) +
                   soundSpeed(pressure, grid.density[i], _gamma);
        if (vMax < vMaxTemp)
        {
            vMax = vMaxTemp;
        }
    }

    _timeStep = _cflNum * _deltaX / vMax;
}
// =============================================================================

// =============================================================================
void Simulation1D::interfaceStates()
{
    // Choose which interface reconstruction technique to use
    if (reconstructionKind == "PCM")
    {
        _piecewiseConstantReconstruction();
    }
    else if (reconstructionKind == "PLM")
    {
        _piecewiseLinearReconstruction();
    }
    else
    {
        throw std::invalid_argument("Invalid kind of initial conditions");
    }

}
// =============================================================================

// =============================================================================
void Simulation1D::solveRiemann()
{
    // Loop through the grid and solve the Riemann problem at each interface
    for (size_t i = grid.numGhostCells;
    i < (grid.numTotCells - grid.numGhostCells + 1);
    i++)
    {
        _riemannSolver.riemannMain(_densityInterfaceL[i],
                                   _velocityInterfaceL[i],
                                   _pressureInterfaceL[i],
                                   _densityInterfaceR[i],
                                   _velocityInterfaceR[i],
                                   _pressureInterfaceR[i],
                                   0.0,
                                   _energyFlux[i],
                                   _momentumFlux[i],
                                   _densityFlux[i]);
    }
}
// =============================================================================

// =============================================================================
// Performe the conservative update
void Simulation1D::conservativeUpdate()
{
    for (size_t i = grid.numGhostCells;
         i < (grid.numTotCells - grid.numGhostCells);
         i++)
    {
        grid.density[i]  = grid.density[i]
                            + (_timeStep / _deltaX)
                            * (_densityFlux[i] - _densityFlux[i+1]);

        grid.momentum[i] = grid.momentum[i]
                            + (_timeStep / _deltaX)
                            * (_momentumFlux[i] - _momentumFlux[i+1]);

        grid.energy[i]   = grid.energy[i]
                            + (_timeStep / _deltaX)
                            * (_energyFlux[i] - _energyFlux[i+1]);
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
                           std::string const &reconstructionKind,
                           std::string const &limiterKind,
                           std::string const &boundaryConditions,
                           std::string const &saveDir)

    // Start by initializing all the const member variables
    : _physLen(physicalLength),
      _cflNum(CFL),
      _deltaX(_physLen / static_cast<double>(reals)),
      _gamma(gamma),
      _riemannSolver(gamma),
      grid(reals, ghosts, saveDir, boundaryConditions),
      currentTime(0.0),
      reconstructionKind(reconstructionKind),
      limiterKind(limiterKind)
{
    // Resize all the arrays
    _densityInterfaceL.resize(grid.numTotCells);
    _velocityInterfaceL.resize(grid.numTotCells);
    _pressureInterfaceL.resize(grid.numTotCells);
    _densityInterfaceR.resize(grid.numTotCells);
    _velocityInterfaceR.resize(grid.numTotCells);
    _pressureInterfaceR.resize(grid.numTotCells);
    _densityFlux.resize(grid.numTotCells);
    _momentumFlux.resize(grid.numTotCells);
    _energyFlux.resize(grid.numTotCells);

    // Set the initial conditions
    _setInitialConditions(initialConditionsKind);
}
// =============================================================================