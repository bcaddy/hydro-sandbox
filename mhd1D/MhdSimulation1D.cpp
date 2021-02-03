/*!
 * \file MhdSimulation1D.cpp
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Implementation of MhdSimulation1D Class
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

#include "MhdSimulation1D.h"
#include "mhdUtilities.h"
#include "HlldRiemannSolver.h"

using namespace mhdUtilities;

// =============================================================================
void MhdSimulation1D::_setInitialConditions(std::string const &initialConditionsKind)
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

        // Iterate over all cells on the left side
        for (size_t i = 0;
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
            i < grid.numTotCells;
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
        for (size_t i = 0;
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
            i < grid.numTotCells;
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
    else
    {
        throw std::invalid_argument("Invalid kind of initial conditions");
    }
}
// =============================================================================

// =============================================================================
void MhdSimulation1D::_centeredMagneticField(Grid1D const &activeGrid)
{
    // Loop through every element of the grid and compute the cell centered
    // electric fields. Note that we have to go 1 element farther than usual
    // since we need the centered state on both sides of the real part of the
    // grid
    for (size_t i = 1;
         i < grid.numTotCells;
         i++)
    {
        for (int j = 0; j < 3; j++)
        {
            _magCentered[i][j] = 0.5 * (activeGrid.magnetic[i][j]
                                      + activeGrid.magnetic[i+1][j]);
        }

    }
}
// =============================================================================

// =============================================================================
double MhdSimulation1D::_slope(std::vector<double> const &primitive,
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
double MhdSimulation1D::_ctSlope(double const &centerL,
                                 double const &faceL,
                                 double const &centerR,
                                 double const &faceR,
                                 double const &velocity)
{
    // Upwinding
    if (velocity > 0.0)
    {
        // Return the slope on the left side
        return (centerL - faceL);
    }
    else if (velocity < 0.0)
    {
        // Return the slope on the right side
        return (centerR - faceR);
    }
    else
    {
        // Return the average of the left and right side slopes
        return 0.5 * ((centerL - faceL) + (centerR - faceR));
    }
}
// =============================================================================

// =============================================================================
void MhdSimulation1D::_piecewiseLinearReconstruction(Grid1D const &workingGrid)
{
    // Compute all the primitive values
    std::vector<double> velocityX(grid.numTotCells),
                        velocityY(grid.numTotCells),
                        velocityZ(grid.numTotCells),
                        magneticX(grid.numTotCells),
                        magneticY(grid.numTotCells),
                        magneticZ(grid.numTotCells),
                        pressure(grid.numTotCells);
    for (size_t i = 1;
         i < grid.numTotCells;
         i++)
    {
        std::vector<double> tempVec = computeVelocity(workingGrid.momentum[i],
                                                     workingGrid.density[i]);
        pressure[i] = computePressure(workingGrid.energy[i],
                                      workingGrid.density[i],
                                      tempVec,
                                      workingGrid.magnetic[i],
                                      _gamma);

        // Assign velocity and magnetic field vectors. Yes I know this is hacky
        // but I don't want to have to change the slope function
        velocityX[i] = tempVec[0];
        velocityY[i] = tempVec[1];
        velocityZ[i] = tempVec[2];
        magneticX[i] = _magCentered[i][0];
        magneticY[i] = _magCentered[i][1];
        magneticZ[i] = _magCentered[i][2];
    }

    // Loop over the entire grid and reconstruct the state within each cell then
    // compute the i-1/2,R and i+1/2,L states
    for (size_t i = 2;
         i < grid.numTotCells - 2 ;
         i++)
    {
        // Slopes
        std::vector<double> slopes {_slope(workingGrid.density, i),
                                    _slope(velocityX,     i),
                                    _slope(velocityY,     i),
                                    _slope(velocityZ,     i),
                                    _slope(magneticX,     i),
                                    _slope(magneticY,     i),
                                    _slope(magneticZ,     i),
                                    _slope(pressure,      i)};

        // ===== Compute the i+1/2,L and i-1/2,R interface states ==============
        // i+1/2,L state
        _interfaceL.density[i+1]     = workingGrid.density[i] +  0.5 * slopes[0];
        _interfaceL.velocity[i+1][0] = velocityX[i]           +  0.5 * slopes[1];
        _interfaceL.velocity[i+1][1] = velocityY[i]           +  0.5 * slopes[2];
        _interfaceL.velocity[i+1][2] = velocityZ[i]           +  0.5 * slopes[3];
        _interfaceL.magnetic[i+1][0] = magneticX[i]           +  0.5 * slopes[4];
        _interfaceL.magnetic[i+1][1] = magneticY[i]           +  0.5 * slopes[5];
        _interfaceL.magnetic[i+1][2] = magneticZ[i]           +  0.5 * slopes[6];
        _interfaceL.pressure[i+1]    = pressure[i]            +  0.5 * slopes[7];

        // i-1/2,R state
        _interfaceR.density[i]     = workingGrid.density[i] -  0.5 * slopes[0];
        _interfaceR.velocity[i][0] = velocityX[i]           -  0.5 * slopes[1];
        _interfaceR.velocity[i][1] = velocityY[i]           -  0.5 * slopes[2];
        _interfaceR.velocity[i][2] = velocityZ[i]           -  0.5 * slopes[3];
        _interfaceR.magnetic[i][0] = magneticX[i]           -  0.5 * slopes[4];
        _interfaceR.magnetic[i][1] = magneticY[i]           -  0.5 * slopes[5];
        _interfaceR.magnetic[i][2] = magneticZ[i]           -  0.5 * slopes[6];
        _interfaceR.pressure[i]    = pressure[i]            -  0.5 * slopes[7];
        // ===== Finished computing the i+1/2,L and i-1/2,R interface states ===
    }
}
// =============================================================================

// =============================================================================
// Implement the piecewise constant reconstruction of the interface states
void MhdSimulation1D::_piecewiseConstantReconstruction(Grid1D const &workingGrid)
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
                                                 workingGrid.magnetic[i-1],
                                                 _gamma);
        _interfaceL.magnetic[i][0] = workingGrid.magnetic[i][0];
        _interfaceL.magnetic[i][1] = _magCentered[i-1][1];
        _interfaceL.magnetic[i][2] = _magCentered[i-1][2];

        // Compute the right interfaces
        _interfaceR.density[i] = workingGrid.density[i];
        _interfaceR.velocity[i] = computeVelocity(workingGrid.momentum[i],
                                                 workingGrid.density[i]);
        _interfaceR.pressure[i] = computePressure(workingGrid.energy[i],
                                                 workingGrid.density[i],
                                                 _interfaceR.velocity[i],
                                                 workingGrid.magnetic[i],
                                                 _gamma);
        _interfaceR.magnetic[i][0] = workingGrid.magnetic[i][0];
        _interfaceR.magnetic[i][1] = _magCentered[i][1];
        _interfaceR.magnetic[i][2] = _magCentered[i][2];
    }
}
// =============================================================================

// =============================================================================
void MhdSimulation1D::computeTimeStep()
{
    // Compute all the cell centered magnetic fields
    _centeredMagneticField(grid);

    // Find the minimum crossing time in the simulation
    double tMin = __DBL_MAX__;

    for (size_t i = grid.numGhostCells;
         i < (grid.numTotCells - grid.numGhostCells);
         i++)
    {
        std::vector<double> velocity = computeVelocity(grid.momentum[i], grid.density[i]);

        double pressure = computePressure(grid.energy[i],
                                          grid.density[i],
                                          velocity,
                                          _magCentered[i],
                                          _gamma);

        double magSonicSpeed = magnetosonicSpeed(pressure,
                                                 grid.density[i],
                                                 _magCentered[i],
                                                 _gamma);

        // Compute the minimum time to cross a cell (ie the maximum wave speed)
        std::vector<double> crossTimes(3);
        crossTimes[0] = _deltaX / (std::abs(velocity[0]) + magSonicSpeed);
        crossTimes[1] = _deltaY / (std::abs(velocity[1]) + magSonicSpeed);
        crossTimes[2] = _deltaZ / (std::abs(velocity[2]) + magSonicSpeed);

        // Check if the value is smaller than the stored value or not
        double localMin = *std::min_element(crossTimes.begin(), crossTimes.end());
        if (localMin < tMin)
        {
            tMin = localMin;
        }
    }
    // compute the time step
    _timeStep = _cflNum * tMin;
}
// =============================================================================

// =============================================================================
void MhdSimulation1D::interfaceStates(std::string const &algoStep)
{
    if (algoStep == "first reconstruction")
    {
        _piecewiseConstantReconstruction(grid);
    }
    else if (algoStep == "second reconstruction")
    {
        // Compute the centered magnetic field
        _centeredMagneticField(_gridHalfTime);

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
            throw std::invalid_argument("Invalid order of reconstruction");
        }
    }
    else
    {
        throw std::invalid_argument("Invalid option for MhdSimulation1D::interfaceStates");
    }
}
// =============================================================================

// =============================================================================
void MhdSimulation1D::solveRiemann()
{
    // Loop through the grid and solve the Riemann problem at each interface
    for (size_t i = 1;
    i < grid.numTotCells;
    i++)
    {
        _riemannSolver->riemannMain(_interfaceL.density[i],
                                    _interfaceL.velocity[i],
                                    _interfaceL.pressure[i],
                                    _interfaceL.magnetic[i],
                                    _interfaceR.density[i],
                                    _interfaceR.velocity[i],
                                    _interfaceR.pressure[i],
                                    _interfaceR.magnetic[i],
                                    _flux.density[i],
                                    _flux.momentum[i],
                                    _flux.magnetic[i],
                                    _flux.energy[i]);

        // Capture the value of the velocities
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                _ctVelocities[i][j][k][0] = _riemannSolver->getVelocityState();
            }
        }
    }

}
// =============================================================================

// =============================================================================
void MhdSimulation1D::ctElectricFields(Grid1D const &activeGrid)
{
    // First we need to compute the cell centered electric fields
    // =========================================================================
    // Declare a 4D vector of shape Nx3x3x3 to store the 3D centered field in each cell
    std::vector<std::vector<std::vector<std::vector<double>>>>
    electricCentered(grid.numTotCells, std::vector<std::vector<std::vector<double>>>(
                      3, std::vector<std::vector<double>>(
                      3, std::vector<double>(
                      3, 0.0))));

    for (size_t i = 0; i < grid.numTotCells; i++)
    {
        // Compute the electric field using a cross product
        std::vector<double> eRef(3, 0.0), velocity(3, 0.0);
        velocity = computeVelocity(activeGrid.momentum[i], activeGrid.density[i]);
        eRef[0] = velocity[2] * activeGrid.magnetic[i][1] - velocity[1] * activeGrid.magnetic[i][2];
        eRef[1] = velocity[0] * activeGrid.magnetic[i][2] - velocity[2] * activeGrid.magnetic[i][0];
        eRef[2] = velocity[1] * activeGrid.magnetic[i][0] - velocity[0] * activeGrid.magnetic[i][1];

        // Now assign the values to the correct parts of the 4D vector
        for (size_t j = 0; j < 3; j++)
        {
            for (size_t k = 0; k < 3; k++)
            {
                for (size_t m = 0; m < 3; m++)
                {
                    electricCentered[i][j][k][m] = eRef[m];
                }
            }
        }
    }
    // Finished computing the centered electric fields
    // =========================================================================

    // Create a virtual grid for the face averaged magnetic fluxes
    // =========================================================================
    // declare the 5d vector used to hold the magnetic fluxes. Indices are as follows
    // [x position]
    // [y position]
    // [z position]
    // [which face, i-1/2, j-1/2, or z-1/2 in that order]
    // [the flux in each direction, x, y, and z respectively]
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>
    magFlux(grid.numTotCells, std::vector<std::vector<std::vector<std::vector<double>>>>(
            3, std::vector<std::vector<std::vector<double>>>(
            3, std::vector<std::vector<double>>(
            3, std::vector<double>(
            3, 0.0)))));

    for (size_t i = 0; i < grid.numTotCells; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            for (size_t k = 0; k < 3; k++)
            {
                for (size_t m = 0; m < 3; m++)
                {
                    magFlux[i][j][k][0][m] = _flux.magnetic[i][m];
                }
            }
        }
    }
    // =========================================================================


    // Then iterate over that grid as to compute all the CT electrid fields as
    // expected in a full 3D simulation
    // =========================================================================
    for (size_t i = 1; i < grid.numTotCells; i++)   // Loop in x-direction
    {
        for (int j = 1; j < 3; j++)  // Loop in y-direction
        {
            for (int k = 1; k < 3; k++)  // Loop in z-direction
            {
                for (int m = 0; m < 3; m++)  // Loop over vector elements
                {
                    // Compute the other two indices that will be needed
                    int m1 = _mod3(m+1), m2 = _mod3(m+2);

                    // All our offset arrays
                    std::vector<int> m2Offset(3, 0), m1Offset(3, 0);

                    /// Compute the offsets
                    /// \todo: consider replacing with vectors that just are the
                    /// the indices. That way I don't have to compute it over
                    /// and over again in indices. This level of optimization is
                    /// worth considering but is probably at the level of
                    /// considering cache per ALU
                    m1Offset[m1] = -1;
                    m2Offset[m2] = -1;

                    // Compute the first term, the sum of surrounding faces
                    double firstTerm = 0.25 * ( magFlux[i][j][k][m1][m]
                                              + magFlux[i][j][k][m2][m]
                                              + magFlux[i + m2Offset[0]][j + m2Offset[0]][k + m2Offset[0]][m1][m]
                                              + magFlux[i + m1Offset[1]][j + m1Offset[1]][k + m1Offset[1]][m2][m]
                                              );
                    // The slopes in the m1 direction
                    double secondTerm = 0.25 * (// The -1/4 slopes
                                                _ctSlope(electricCentered[i][k][k][m],
                                                         magFlux[i][j][k][m1][m],
                                                         electricCentered[i + m2Offset[0]][j + m2Offset[0]][k + m2Offset[0]][m],
                                                         magFlux[i + m2Offset[0]][j + m2Offset[0]][k + m2Offset[0]][m1][m],
                                                         _ctVelocities[i][j][k][m2][m2])
                                                // The -3/4 slopes
                                              - _ctSlope(electricCentered[i + m1Offset[1]][j + m1Offset[1]][k + m1Offset[1]][m],
                                                         magFlux[i][j][k][m1][m],
                                                         electricCentered[i + m1Offset[1] + m2Offset[0]][j + m1Offset[1] + m2Offset[1]][k + m1Offset[1] + m2Offset[2]][m],
                                                         magFlux[i + m2Offset[0]][j + m2Offset[0]][k + m2Offset[0]][m1][m],
                                                         _ctVelocities[i + m2Offset[0]][j + m2Offset[0]][k + m2Offset[0]][m2][m2]));

                    // The slopes in the m2 directions
                    double thirdTerm = 0.25 * (// The -1/4 slopes
                                                _ctSlope(electricCentered[i + m1Offset[1]][j + m1Offset[1]][k + m1Offset[1]][m],
                                                         magFlux[i + m1Offset[1]][j + m1Offset[1]][k + m1Offset[1]][m2][m],
                                                         electricCentered[i][j][k][m],
                                                         magFlux[i][j][k][m2][m],
                                                         _ctVelocities[i][j][k][m1][m1])
                                                // The -3/4 slopes
                                              - _ctSlope(electricCentered[i + m1Offset[1] + m2Offset[0]][j + m1Offset[1] + m2Offset[1]][k + m1Offset[1] + m2Offset[2]][m],
                                                         magFlux[i + m1Offset[1]][j + m1Offset[1]][k + m1Offset[1]][m2][m],
                                                         electricCentered[i + m2Offset[0]][j + m2Offset[0]][k + m2Offset[0]][m],
                                                         magFlux[i][j][k][m2][m],
                                                         _ctVelocities[i + m2Offset[0]][j + m2Offset[0]][k + m2Offset[0]][m1][m1]));

                    // Now we fill in the array of edge values
                    _edgeFields[i][j][k][m] = firstTerm + secondTerm + thirdTerm;
                }
            }
        }
    }
    // Done computing the CT electric fields
    // =========================================================================
}
// =============================================================================

// =============================================================================
// Performe the conservative update
void MhdSimulation1D::conservativeUpdate(std::string const &timeChoice)
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
        throw std::invalid_argument("Invalid option for MhdSimulation1D::conservativeUpdate");
    }

    for (size_t i = 1;
         i < grid.numTotCells - 1;
         i++)
    {
        destinationGrid->density[i]  = sourceGrid->density[i]
                                      + (localTimeStep / _deltaX)
                                      * (_flux.density[i] - _flux.density[i+1]);

        destinationGrid->energy[i]   = sourceGrid->energy[i]
                                      + (localTimeStep / _deltaX)
                                      * (_flux.energy[i] - _flux.energy[i+1]);

        // Update momentum
        // todo this might be better in a loop but idgaf right now
        destinationGrid->momentum[i][0] = sourceGrid->momentum[i][0]
                                            + (localTimeStep / _deltaX)
                                            * (_flux.momentum[i][0] - _flux.momentum[i+1][0]);
        destinationGrid->momentum[i][1] = sourceGrid->momentum[i][1]
                                            + (localTimeStep / _deltaY)
                                            * (_flux.momentum[i][1] - _flux.momentum[i+1][1]);
        destinationGrid->momentum[i][2] = sourceGrid->momentum[i][2]
                                            + (localTimeStep / _deltaZ)
                                            * (_flux.momentum[i][2] - _flux.momentum[i+1][2]);

        // Update magnetic field using CT fields
        // todo this might be better in a loop but idgaf right now
        destinationGrid->magnetic[i][0] = sourceGrid->magnetic[i][0]
                                            - (localTimeStep / _deltaY)
                                            * (_edgeFields[i][2][1][2] - _edgeFields[i][1][1][2])
                                            + (localTimeStep / _deltaZ)
                                            * (_edgeFields[i][1][2][1] - _edgeFields[i][1][1][1]);

        destinationGrid->magnetic[i][2] = sourceGrid->magnetic[i][1]
                                            - (localTimeStep / _deltaX)
                                            * (_edgeFields[i+1][1][1][2] - _edgeFields[i][1][1][2])
                                            + (localTimeStep / _deltaZ)
                                            * (_edgeFields[i][1][2][0] - _edgeFields[i][1][1][0]);

        destinationGrid->magnetic[i][2] = sourceGrid->magnetic[i][2]
                                            - (localTimeStep / _deltaX)
                                            * (_edgeFields[i+1][1][1][1] - _edgeFields[i][1][1][1])
                                            + (localTimeStep / _deltaY)
                                            * (_edgeFields[i][2][1][0] - _edgeFields[i][1][1][0]);
    }

    // Release pointers
    sourceGrid.release();
    destinationGrid.release();
}
// =============================================================================

// =============================================================================
// Constructor
MhdSimulation1D::MhdSimulation1D(double const &physicalLength,
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
      _deltaY(_deltaX),
      _deltaZ(_deltaX),
      _gamma(gamma),
      _magCentered(2 * _numGhosts + reals, std::vector<double> (3, 0)),
      _interfaceL(reals, _numGhosts),
      _interfaceR(reals, _numGhosts),
      _flux(reals, _numGhosts),
      _gridHalfTime(reals, _numGhosts, boundaryConditions),
      _edgeFields(2 * _numGhosts + reals, std::vector<std::vector<std::vector<double>>>(
                  3, std::vector<std::vector<double>>(
                  3, std::vector<double>(
                  3, 0.0)))),
      _ctVelocities(2 * _numGhosts + reals, std::vector<std::vector<std::vector<std::vector<double>>>>(
                    3, std::vector<std::vector<std::vector<double>>>(
                    3, std::vector<std::vector<double>>(
                    3, std::vector<double>(
                    3, 0.0))))),
      grid(reals, _numGhosts, boundaryConditions, saveDir),
      currentTime(0.0),
      reconstructionKind(reconstructionKind),
      limiterKind(limiterKind),
      riemannSolverKind(riemannSolverKind)
{
    // Set the initial conditions
    _setInitialConditions(initialConditionsKind);

    // Choose the Riemann Solver
    if (riemannSolverKind == "HLLD")
    {
        _riemannSolver = std::unique_ptr<MhdRiemannSolver>(new HlldRiemannSolver(_gamma));
    }
    else
    {
        throw std::invalid_argument("Invalid kind of Riemann Solver");
    }
}
// =============================================================================