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
#include <memory>

#include "Simulation1D.h"
#include "HydroUtilities.h"
#include "HllcRiemannSolver.h"

using namespace HydroUtilities;

// =============================================================================
void Simulation1D::_setInitialConditions(std::string const &initialConditionsKind)
{
    if (initialConditionsKind == "dwShockTube")
    {
        size_t const half  = grid.numTotCells / 2;
        double coef = 1. / std::sqrt(4. * M_PI);
        double denL  = 1.08,
               presL = 0.95,
               denR  = 1.,
               presR = 1.;

        std::vector<double> velL{1.2, 0.01, 0.5},
                            bL{4. * coef, 3.6 * coef, 2.0 * coef},
                            velR{0.0, 0.0, 0.0},
                            bR{4.0 * coef, 4.0 * coef, 2.0 * coef};

        // Iterate over just the real cells on the left side
        for (size_t i = grid.numGhostCells;
            i < half;
            i++)
        {
            grid.density[i]     = denL;
            grid.momentum[i][0] = computeMomentum(velL[0], denL);
            grid.momentum[i][1] = computeMomentum(velL[1], denL);
            grid.momentum[i][2] = computeMomentum(velL[2], denL);
            grid.magnetic[i][0] = bL[0];
            grid.magnetic[i][1] = bL[1];
            grid.magnetic[i][2] = bL[2];
            grid.energy[i]      = computeEnergy(presL, denL, velL, bL, _gamma);
        }

        // Iterate over the real cells on the right side
        for (size_t i = half;
            i < (grid.numTotCells - grid.numGhostCells);
            i++)
        {
            grid.density[i]     = denR;
            grid.momentum[i][0] = computeMomentum(velL[0], denR);
            grid.momentum[i][1] = computeMomentum(velL[1], denR);
            grid.momentum[i][2] = computeMomentum(velL[2], denR);
            grid.magnetic[i][0] = bR[0];
            grid.magnetic[i][1] = bR[1];
            grid.magnetic[i][2] = bR[2];
            grid.energy[i]      = computeEnergy(presR, denR, velR, bR, _gamma);
        }
    }
    else if (initialConditionsKind == "bwShockTube")
    {
        size_t const half  = grid.numTotCells / 2;
        double denL  = 1.0,
               presL = 1.0,
               denR  = 0.125,
               presR = 0.1;

        std::vector<double> velL{0.0, 0.0, 0.0},
                            bL{0.75, 1.0, 0.0},
                            velR{0.0, 0.0, 0.0},
                            bR{0.75, -1.0, 0.0};

        // Iterate over just the real cells on the left side
        for (size_t i = grid.numGhostCells;
            i < half;
            i++)
        {
            grid.density[i]     = denL;
            grid.momentum[i][0] = computeMomentum(velL[0], denL);
            grid.momentum[i][1] = computeMomentum(velL[1], denL);
            grid.momentum[i][2] = computeMomentum(velL[2], denL);
            grid.magnetic[i][0] = bL[0];
            grid.magnetic[i][1] = bL[1];
            grid.magnetic[i][2] = bL[2];
            grid.energy[i]      = computeEnergy(presL, denL, velL, bL, _gamma);
        }

        // Iterate over the real cells on the right side
        for (size_t i = half;
            i < (grid.numTotCells - grid.numGhostCells);
            i++)
        {
            grid.density[i]     = denR;
            grid.momentum[i][0] = computeMomentum(velL[0], denR);
            grid.momentum[i][1] = computeMomentum(velL[1], denR);
            grid.momentum[i][2] = computeMomentum(velL[2], denR);
            grid.magnetic[i][0] = bR[0];
            grid.magnetic[i][1] = bR[1];
            grid.magnetic[i][2] = bR[2];
            grid.energy[i]      = computeEnergy(presR, denR, velR, bR, _gamma);
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
void Simulation1D::_piecewiseLinearReconstruction(Grid1D const &workingGrid)
{
    // Compute all the primitive values
    std::vector<double> velocity(grid.numTotCells), pressure(grid.numTotCells);
    for (size_t i = 1;
         i < grid.numTotCells;
         i++)
    {
        velocity[i] = computeVelocity(workingGrid.momentum[i],
                                      workingGrid.density[i]);
        pressure[i] = computePressure(workingGrid.energy[i],
                                      workingGrid.density[i],
                                      velocity[i],
                                      _gamma);
    }

    // Loop over the entire grid and reconstruct the state within each cell then
    // compute the i-1/2,R and i+1/2,L states
    for (size_t i = 2;
         i < grid.numTotCells - 2 ;
         i++)
    {
        // Slopes
        std::vector<double> slopes {_slope(workingGrid.density, i),
                                    _slope(velocity,     i),
                                    _slope(pressure,     i)};

        // ===== Compute the i+1/2,L and i-1/2,R interface states ==============
        // i+1/2,L state
        _interfaceL.density[i+1]  = workingGrid.density[i] +  0.5 * slopes[0];
        _interfaceL.velocity[i+1] = velocity[i]            +  0.5 * slopes[1];
        _interfaceL.pressure[i+1] = pressure[i]            +  0.5 * slopes[2];

        // i-1/2,R state
        _interfaceR.density[i]  = workingGrid.density[i] - 0.5 * slopes[0];
        _interfaceR.velocity[i] = velocity[i]            - 0.5 * slopes[1];
        _interfaceR.pressure[i] = pressure[i]            - 0.5 * slopes[2];
        // ===== Finished computing the i+1/2,L and i-1/2,R interface states ===
    }
}
// =============================================================================

// =============================================================================
// Implement the piecewise constant reconstruction of the interface states
void Simulation1D::_piecewiseConstantReconstruction(Grid1D const &workingGrid)
{
    // Loop through every element of the grid and compute the interface states.
    // Note that we have to go 1 element farther than usual since we need the
    // interface states on both sides of the real part of the grid
    for (size_t i = 1;
         i < grid.numTotCells;
         i++)
    {
        // Compute the left interfaces
        _interfaceL.density[i] = workingGrid.density[i-1];
        _interfaceL.velocity[i] = computeVelocity(workingGrid.momentum[i-1],
                                                 workingGrid.density[i-1]);
        _interfaceL.pressure[i] = computePressure(workingGrid.energy[i-1],
                                                 workingGrid.density[i-1],
                                                 _interfaceL.velocity[i],
                                                 _gamma);

        // Compute the right interfaces
        _interfaceR.density[i] = workingGrid.density[i];
        _interfaceR.velocity[i] = computeVelocity(workingGrid.momentum[i],
                                                 workingGrid.density[i]);
        _interfaceR.pressure[i] = computePressure(workingGrid.energy[i],
                                                 workingGrid.density[i],
                                                 _interfaceR.velocity[i],
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
void Simulation1D::interfaceStates(std::string const &algoStep)
{
    if (algoStep == "first reconstruction")
    {
        _piecewiseConstantReconstruction(grid);
    }
    else if (algoStep == "second reconstruction")
    {
            // Choose which interface reconstruction technique to use
        if (reconstructionKind == "PCM")
        {
            _piecewiseConstantReconstruction(_gridHalfTime);
        }
        else if (reconstructionKind == "PLM")
        {
            _piecewiseLinearReconstruction(_gridHalfTime);
        }
        else
        {
            throw std::invalid_argument("Invalid kind of initial conditions");
        }
    }
    else
    {
        throw std::invalid_argument("Invalid option for Simulation1D::interfaceStates");
    }
}
// =============================================================================

// =============================================================================
void Simulation1D::solveRiemann()
{
    // Loop through the grid and solve the Riemann problem at each interface
    for (size_t i = 1;
    i < grid.numTotCells;
    i++)
    {
        _riemannSolver->riemannMain(_interfaceL.density[i],
                                    _interfaceL.velocity[i],
                                    _interfaceL.pressure[i],
                                    _interfaceR.density[i],
                                    _interfaceR.velocity[i],
                                    _interfaceR.pressure[i],
                                    _flux.density[i],
                                    _flux.momentum[i],
                                    _flux.energy[i]);
    }
}
// =============================================================================

// =============================================================================
// Performe the conservative update
void Simulation1D::conservativeUpdate(std::string const &timeChoice)
{
    // Choose which grid to update
    std::unique_ptr<Grid1D> sourceGrid;
    std::unique_ptr<Grid1D> destinationGrid;
    double localTimeStep;
    if (timeChoice == "half time update")
    {
        sourceGrid      = std::unique_ptr<Grid1D>(&grid);
        destinationGrid = std::unique_ptr<Grid1D>(&_gridHalfTime);
        localTimeStep = 0.5 * _timeStep;
    }
    else if (timeChoice == "full time update")
    {
        sourceGrid      = std::unique_ptr<Grid1D>(&grid);
        destinationGrid = std::unique_ptr<Grid1D>(&grid);
        localTimeStep = _timeStep;
    }
    else
    {
        throw std::invalid_argument("Invalid option for Simulation1D::conservativeUpdate");
    }

    for (size_t i = 1;
         i < grid.numTotCells - 1;
         i++)
    {
        destinationGrid->density[i]  = sourceGrid->density[i]
                                      + (localTimeStep / _deltaX)
                                      * (_flux.density[i] - _flux.density[i+1]);

        destinationGrid->momentum[i] = sourceGrid->momentum[i]
                                      + (localTimeStep / _deltaX)
                                      * (_flux.momentum[i] - _flux.momentum[i+1]);

        destinationGrid->energy[i]   = sourceGrid->energy[i]
                                      + (localTimeStep / _deltaX)
                                      * (_flux.energy[i] - _flux.energy[i+1]);
    }

    // Release pointers
    sourceGrid.release();
    destinationGrid.release();
}
// =============================================================================

// =============================================================================
// Constructor
Simulation1D::Simulation1D(double const &physicalLength,
                           double const &gamma,
                           double const &CFL,
                           size_t const &reals,
                           std::string const &initialConditionsKind,
                           std::string const &reconstructionKind,
                           std::string const &limiterKind,
                           std::string const &riemannSolverKind,
                           std::string const &boundaryConditions,
                           std::string const &saveDir)

    // Start by initializing all the const member variables
    : _physLen(physicalLength),
      _cflNum(CFL),
      _deltaX(_physLen / static_cast<double>(reals)),
      _gamma(gamma),
      _interfaceL(reals, _numGhosts),
      _interfaceR(reals, _numGhosts),
      _flux(reals, _numGhosts),
      _gridHalfTime(reals, _numGhosts, boundaryConditions),
      grid(reals, _numGhosts, boundaryConditions, saveDir),
      currentTime(0.0),
      reconstructionKind(reconstructionKind),
      limiterKind(limiterKind),
      riemannSolverKind(riemannSolverKind)
{
    // Set the initial conditions
    _setInitialConditions(initialConditionsKind);

    // Choose the Riemann Solver
    if (riemannSolverKind == "HLLC")
    {
        _riemannSolver = std::unique_ptr<RiemannSolver>(new HllcRiemannSolver(_gamma));
    }
    else
    {
        throw std::invalid_argument("Invalid kind of Riemann Solver");
    }
}
// =============================================================================