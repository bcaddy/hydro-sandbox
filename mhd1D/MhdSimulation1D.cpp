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
#include <iostream>

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

        std::vector<double> velL={1.2, 0.01, 0.5},
                            bL={4. * coef, 3.6 * coef, 2.0 * coef},
                            velR={0.0, 0.0, 0.0},
                            bR={4.0 * coef, 4.0 * coef, 2.0 * coef};

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
            grid.momentum[i][0] = computeMomentum(velR[0], denR);
            grid.momentum[i][1] = computeMomentum(velR[1], denR);
            grid.momentum[i][2] = computeMomentum(velR[2], denR);
            grid.magnetic[i][0] = bR[0];
            grid.magnetic[i][1] = bR[1];
            grid.magnetic[i][2] = bR[2];
            grid.energy[i]      = computeEnergy(presR, denR, velR, bR, _gamma);
        }
        // Set the last face in the magnetic field
        grid.magnetic[grid.numTotCells][0] = bR[0];
        grid.magnetic[grid.numTotCells][1] = bR[1];
        grid.magnetic[grid.numTotCells][2] = bR[2];
    }
    else if (initialConditionsKind == "bwShockTube")
    {
        size_t const half  = grid.numTotCells / 2;
        double denL  = 1.0,
               presL = 1.0,
               denR  = 0.125,
               presR = 0.1;

        std::vector<double> velL={0.0, 0.0, 0.0},
                            bL={0.75, 1.0, 0.0},
                            velR={0.0, 0.0, 0.0},
                            bR={0.75, -1.0, 0.0};

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
            grid.momentum[i][0] = computeMomentum(velR[0], denR);
            grid.momentum[i][1] = computeMomentum(velR[1], denR);
            grid.momentum[i][2] = computeMomentum(velR[2], denR);
            grid.magnetic[i][0] = bR[0];
            grid.magnetic[i][1] = bR[1];
            grid.magnetic[i][2] = bR[2];
            grid.energy[i]      = computeEnergy(presR, denR, velR, bR, _gamma);
        }
        // Set the last face in the magnetic field
        grid.magnetic[grid.numTotCells][0] = bR[0];
        grid.magnetic[grid.numTotCells][1] = bR[1];
        grid.magnetic[grid.numTotCells][2] = bR[2];
    }
    else if (initialConditionsKind == "chollaSodShockTube")
    {
        size_t const half  = grid.numTotCells / 2;
        double denL  = 1.0,
               presL = 1.0,
               denR  = 0.1,
               presR = 0.1;

        std::vector<double> velL={0.0, 0.0, 0.0},
                            bL={0.0, 0.0, 0.0},
                            velR={0.0, 0.0, 0.0},
                            bR={0.0, 0.0, 0.0};

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
            grid.momentum[i][0] = computeMomentum(velR[0], denR);
            grid.momentum[i][1] = computeMomentum(velR[1], denR);
            grid.momentum[i][2] = computeMomentum(velR[2], denR);
            grid.magnetic[i][0] = bR[0];
            grid.magnetic[i][1] = bR[1];
            grid.magnetic[i][2] = bR[2];
            grid.energy[i]      = computeEnergy(presR, denR, velR, bR, _gamma);
        }
        // Set the last face in the magnetic field
        grid.magnetic[grid.numTotCells][0] = bR[0];
        grid.magnetic[grid.numTotCells][1] = bR[1];
        grid.magnetic[grid.numTotCells][2] = bR[2];
    }
    else if (initialConditionsKind.substr(0,10) == "singleWave")
    {
        // Setup the background state
        double backgroundPres = 3./5., backgroundDen = 1.;
        std::vector<double> backgroundVel = {0.,0.,0.},
                            backgroundMag = {1., std::sqrt(2.), 0.5};

        double backgroundEnergy = computeEnergy(backgroundPres, backgroundDen, backgroundVel, backgroundMag, _gamma);
        std::vector<double> backgroundMom = {computeMomentum(backgroundVel[0], backgroundDen),
                                             computeMomentum(backgroundVel[1], backgroundDen),
                                             computeMomentum(backgroundVel[2], backgroundDen)};

        // Set the wave amplitude
        double amp = 0.1;//1.E-6;

        // Choose left or right moving wave
        double lrSign = (initialConditionsKind.substr(11,1) == "R")? 1.: -1.;

        // Choose the correct right eigenvector
        double rightVecEnergy, rightVecDen;
        std::vector<double> rightVecMom(3), rightVecMag(3);
        if (initialConditionsKind.substr(10,1) == "C")
        {   // The contact/entropy wave
            rightVecDen    = 1.;
            rightVecMom    = {lrSign, 0., 0.};
            rightVecMag    = {0., 0., 0.};
            rightVecEnergy = 0.5;

            backgroundVel[0] = lrSign;
            backgroundMom[0] = computeMomentum(backgroundVel[0], backgroundDen);
            backgroundEnergy = computeEnergy(backgroundPres, backgroundDen, backgroundVel, backgroundMag, _gamma);
        }
        else if (initialConditionsKind.substr(10,1) == "F")
        {   // The fast magnetosonic wave
            double coef = 1. / (6. * std::sqrt(5.));

            rightVecDen    = coef * 6.;
            rightVecMom    = {coef * lrSign * 12.,
                              coef * (-lrSign) * 4. * std::sqrt(2.),
                              coef * (-lrSign) * 2.};
            rightVecMag    = {0.,
                              coef * 8. * std::sqrt(2.),
                              coef * 4.};
            rightVecEnergy = coef * 27.;
        }
        else if (initialConditionsKind.substr(10,1) == "S")
        {   // The slow magnetosonic wave
            double coef = 1. / (6. * std::sqrt(5.));

            rightVecDen    = coef * 12.;
            rightVecMom    = {coef * lrSign * 6.,
                              coef * lrSign * 8. * std::sqrt(2.),
                              coef * lrSign * 4.};
            rightVecMag    = {0.,
                              -coef * 4. * std::sqrt(2.),
                              -coef * 2.};
            rightVecEnergy = coef * 9.;
        }
        else if (initialConditionsKind.substr(10,1) == "A")
        {  // The Alfven wave
            double coef = 1. / 3.;

            rightVecDen    = 0.;
            rightVecMom    = {0.,
                              coef * lrSign,
                              coef * (-lrSign) * 2. * std::sqrt(2.)};
            rightVecMag    = {0.,
                              coef * (-1.),
                              coef * 2. * std::sqrt(2.)};
            rightVecEnergy = 0.;
        }
        else
        {
            throw std::invalid_argument("Invalid kind of initial conditions in singleWave");
        }

        // Choose hydro only wave
        if (initialConditionsKind.substr(12,1) == "H")
        {
            backgroundMag = {0., 0., 0.};
            rightVecMag = {0., 0., 0.};
        }

        // Now we actually compute the initial conditions
        double const twoPi = 2.* M_PI; // just to save some compute time
        double lFacePosition = 0.;
        double offset = 0.25;
        for (size_t i = grid.numGhostCells;
            i < (grid.numTotCells - grid.numGhostCells); i++)
        {
            // Compute the positions and the sines required
            lFacePosition  = (i - grid.numGhostCells) * _deltaX + offset;
            double centerPosition = lFacePosition + _deltaX/2.;

            double faceSine   = std::sin(twoPi * lFacePosition);
            double centerSine = std::sin(twoPi * centerPosition);

            // Set the state at this grid point
            grid.density[i]     = backgroundDen    + amp * rightVecDen    * centerSine;
            grid.momentum[i][0] = backgroundMom[0] + amp * rightVecMom[0] * centerSine;
            grid.momentum[i][1] = backgroundMom[1] + amp * rightVecMom[1] * centerSine;
            grid.momentum[i][2] = backgroundMom[2] + amp * rightVecMom[2] * centerSine;
            grid.magnetic[i][0] = backgroundMag[0] + amp * rightVecMag[0] * faceSine;
            grid.magnetic[i][1] = backgroundMag[1] + amp * rightVecMag[1] * faceSine;
            grid.magnetic[i][2] = backgroundMag[2] + amp * rightVecMag[2] * faceSine;
            grid.energy[i]      = backgroundEnergy + amp * rightVecEnergy * centerSine;
        }
        // lastly compute the final magnetic field face
        double faceSine = std::sin(twoPi * (lFacePosition + _deltaX));
        grid.magnetic[grid.numTotCells - grid.numGhostCells][0] = backgroundMag[0] + amp * rightVecMag[0] * faceSine;
        grid.magnetic[grid.numTotCells - grid.numGhostCells][1] = backgroundMag[1] + amp * rightVecMag[1] * faceSine;
        grid.magnetic[grid.numTotCells - grid.numGhostCells][2] = backgroundMag[2] + amp * rightVecMag[2] * faceSine;
    }
    else if (initialConditionsKind == "indexCheck")
    {
        for (size_t i = grid.numGhostCells;
             i < (grid.numTotCells - grid.numGhostCells); i++)
        {
            // Set the state at this grid point
            grid.density[i]     = double(i);
            grid.momentum[i][0] = double(i);
            grid.momentum[i][1] = double(i);
            grid.momentum[i][2] = double(i);
            grid.magnetic[i][0] = double(i);
            grid.magnetic[i][1] = double(i);
            grid.magnetic[i][2] = double(i);
            grid.energy[i]      = computeEnergy(double(i), double(i), grid.momentum[i], grid.magnetic[i], _gamma);
        }
    }
    else
    {
        throw std::invalid_argument("Invalid kind of initial conditions");
    }
}
// =============================================================================

// =============================================================================
void MhdSimulation1D::_centeredMagneticField(Grid1D const &workingGrid)
{
    // Loop through every element of the grid and compute the cell centered
    // electric fields. Note that we have to go 1 element farther than usual
    // since we need the centered state on both sides of the real part of the
    // grid
    for (size_t i = 0;
         i < grid.numTotCells;
         i++)
    {
        for (int j = 0; j < 3; j++)
        {
            _magCentered[i][j] = 0.5 * (workingGrid.magnetic[i][j]
                                      + workingGrid.magnetic[i+1][j]);
        }

    }

    /// \note The issue is tha we're missing the right most face-centered
    /// magnetic field. In general to set that we need the boundary conditions
    /// or we just set the average in the cell to the value at the face. This
    /// last option will effectively break periodic boundary conditions.
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
double MhdSimulation1D::_ctSlope(double const &centerNeg,
                                 double const &faceNeg,
                                 double const &centerPos,
                                 double const &facePos,
                                 double const &velocity)
{
    // Upwinding
    if (velocity > 0.0)
    {
        // Return the slope on the left side
        return (centerNeg - faceNeg);
    }
    else if (velocity < 0.0)
    {
        // Return the slope on the right side
        return (centerPos - facePos);
    }
    else
    {
        // Return the average of the left and right side slopes
        return 0.5 * ((centerNeg - faceNeg) + (centerPos - facePos));
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
                                    _slope(magneticY,     i),
                                    _slope(magneticZ,     i),
                                    _slope(pressure,      i)};

        // ===== Compute the i+1/2,L and i-1/2,R interface states ==============
        // i+1/2,L state
        _interfaceL.density[i+1]     = workingGrid.density[i] +  0.5 * slopes[0];
        _interfaceL.velocity[i+1][0] = velocityX[i]           +  0.5 * slopes[1];
        _interfaceL.velocity[i+1][1] = velocityY[i]           +  0.5 * slopes[2];
        _interfaceL.velocity[i+1][2] = velocityZ[i]           +  0.5 * slopes[3];
        _interfaceL.magnetic[i+1][0] = workingGrid.magnetic[i+1][0];
        _interfaceL.magnetic[i+1][1] = magneticY[i]           +  0.5 * slopes[4];
        _interfaceL.magnetic[i+1][2] = magneticZ[i]           +  0.5 * slopes[5];
        _interfaceL.pressure[i+1]    = pressure[i]            +  0.5 * slopes[6];

        // i-1/2,R state
        _interfaceR.density[i]     = workingGrid.density[i] -  0.5 * slopes[0];
        _interfaceR.velocity[i][0] = velocityX[i]           -  0.5 * slopes[1];
        _interfaceR.velocity[i][1] = velocityY[i]           -  0.5 * slopes[2];
        _interfaceR.velocity[i][2] = velocityZ[i]           -  0.5 * slopes[3];
        _interfaceR.magnetic[i][0] = workingGrid.magnetic[i][0];
        _interfaceR.magnetic[i][1] = magneticY[i]           -  0.5 * slopes[4];
        _interfaceR.magnetic[i][2] = magneticZ[i]           -  0.5 * slopes[5];
        _interfaceR.pressure[i]    = pressure[i]            -  0.5 * slopes[6];
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
            for (int k = 0; k < 3; k++)
            {
                _ctVelocities[i][j][k][0] = _riemannSolver->getVelocityState();
            }
        }
    }

}
// =============================================================================

// =============================================================================
void MhdSimulation1D::ctElectricFields(std::string const &timeChoice)
{
    // =========================================================================
    // We need to compute the cell centered electric fields
    // =========================================================================
    // First we choose the working grid
    std::unique_ptr<Grid1D> workingGrid;
    if (timeChoice == "Half time")
    {
        workingGrid = std::unique_ptr<Grid1D>(&grid);
    }
    else if (timeChoice == "Full time")
    {
        workingGrid = std::unique_ptr<Grid1D>(&_gridHalfTime);
    }
    else
    {
        throw std::invalid_argument("Invalid option for MhdSimulation1D::ctElectricFields");
    }
    // =========================================================================
    // Finish choosing grid
    // =========================================================================

    // =========================================================================
    // Declare a 4D vector of shape Nx3x3x3 to store the 3D centered field in
    // each cell
    // =========================================================================
    // declare the 4D vector used to hold the centered electric field. Indices
    // are as follows
    // [x position]
    // [y position]
    // [z position]
    // [the magnetic field in each direction, x, y, and z respectively]
    stdVector4D electricCentered(grid.numTotCells, stdVector3D(
                                                3, stdVector2D(
                                                3, stdVector1D(
                                                3, 0.0))));

    // Loop through the cells to compute all the electric fields then assign
    // them to their correct spots
    for (size_t i = 0; i < grid.numTotCells; i++)
    {
        // Compute the electric field using a cross product
        stdVector1D eRef(3, 0.0), velocity(3, 0.0);

        velocity = computeVelocity(workingGrid->momentum[i], workingGrid->density[i]);

        eRef[0]  = velocity[2] * workingGrid->magnetic[i][1] - velocity[1] * workingGrid->magnetic[i][2];
        eRef[1]  = velocity[0] * workingGrid->magnetic[i][2] - velocity[2] * workingGrid->magnetic[i][0];
        eRef[2]  = velocity[1] * workingGrid->magnetic[i][0] - velocity[0] * workingGrid->magnetic[i][1];


        // Now assign the values to the correct parts of the 4D vector
        for (size_t j = 0; j < 3; j++)
        {
            for (size_t k = 0; k < 3; k++)
            {
                electricCentered[i][j][k] = eRef;
            }
        }
    }
    // We can release the workingGrid pointer now
    workingGrid.release();
    // =========================================================================
    // Finished computing the centered electric fields
    // =========================================================================

    // =========================================================================
    // Create a virtual grid for the face averaged magnetic fluxes
    // =========================================================================
    // declare the 5D vector used to hold the magnetic fluxes. Indices are as follows
    // [x position]
    // [y position]
    // [z position]
    // [which face, i-1/2, j-1/2, or z-1/2 in that order]
    // [the flux in each direction, x, y, and z respectively]
    stdVector5D magFlux(grid.numTotCells, stdVector4D(
                                       3, stdVector3D(
                                       3, stdVector2D(
                                       3, stdVector1D(
                                       3, 0.0)))));

    for (size_t i = 0; i < grid.numTotCells; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            for (size_t k = 0; k < 3; k++)
            {
                magFlux[i][j][k][0] = _flux.magnetic[i];
            }
        }
    }
    // =========================================================================
    // Finished virtual magnetic flux grid
    // =========================================================================

    // =========================================================================
    // Then iterate over that grid as to compute all the CT electrid fields as
    // expected in a full 3D simulation
    // =========================================================================
    for (size_t i = 1; i < grid.numTotCells; i++)   // Loop in x-direction
    {
        for (size_t j = 1; j < 3; j++)  // Loop in y-direction
        {
            for (size_t k = 1; k < 3; k++)  // Loop in z-direction
            {
                for (size_t m = 0; m < 3; m++)  // Loop over vector elements
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
                    double firstTerm = ( magFlux[i][j][k][m1][m]
                                       + magFlux[i][j][k][m2][m]
                                       + magFlux[i + m2Offset[0]][j + m2Offset[1]][k + m2Offset[2]][m1][m]
                                       + magFlux[i + m1Offset[0]][j + m1Offset[1]][k + m1Offset[2]][m2][m]);

                    // The slopes in the m1 direction
                    double secondTerm = // The -1/4 slopes
                                        _ctSlope(electricCentered[i + m2Offset[0]][j + m2Offset[1]][k + m2Offset[2]][m],
                                                 magFlux[i + m2Offset[0]][j + m2Offset[1]][k + m2Offset[2]][m1][m],
                                                 electricCentered[i][j][k][m],
                                                 magFlux[i][j][k][m1][m],
                                                 _ctVelocities[i][j][k][m2][m2])
                                        // The -3/4 slopes
                                        + _ctSlope(electricCentered[i + m1Offset[0] + m2Offset[0]][j + m1Offset[1] + m2Offset[1]][k + m1Offset[2] + m2Offset[2]][m],
                                                   magFlux[i + m2Offset[0]][j + m2Offset[1]][k + m2Offset[2]][m1][m],
                                                   electricCentered[i + m1Offset[0]][j + m1Offset[1]][k + m1Offset[2]][m],
                                                   magFlux[i][j][k][m1][m],
                                                   _ctVelocities[i + m1Offset[0]][j + m1Offset[1]][k + m1Offset[2]][m2][m2]);

                    // The slopes in the m2 directions
                    double thirdTerm =  // The -1/4 slopes
                                        _ctSlope(electricCentered[i + m1Offset[0]][j + m1Offset[1]][k + m1Offset[2]][m],
                                                 magFlux[i + m1Offset[0]][j + m1Offset[1]][k + m1Offset[2]][m2][m],
                                                 electricCentered[i][j][k][m],
                                                 magFlux[i][j][k][m2][m],
                                                 _ctVelocities[i][j][k][m1][m1])
                                        // The -3/4 slopes
                                        + _ctSlope(electricCentered[i + m1Offset[0] + m2Offset[0]][j + m1Offset[1] + m2Offset[1]][k + m1Offset[2] + m2Offset[2]][m],
                                                   magFlux[i + m1Offset[0]][j + m1Offset[1]][k + m1Offset[2]][m2][m],
                                                   electricCentered[i + m2Offset[0]][j + m2Offset[1]][k + m2Offset[2]][m],
                                                   magFlux[i][j][k][m2][m],
                                                   _ctVelocities[i + m2Offset[0]][j + m2Offset[1]][k + m2Offset[2]][m1][m1]);

                    // Now we fill in the array of edge values
                    _edgeFields[i][j][k][m] = 0.25 * (firstTerm + secondTerm + thirdTerm);
                }
            }
        }
    }
    // =========================================================================
    // Done computing the CT electric fields
    // =========================================================================
}
// =============================================================================

// =============================================================================
// Performe the conservative update
void MhdSimulation1D::conservativeUpdate(std::string const &timeChoice)
{
    // Choose which grid to update
    std::unique_ptr<Grid1D> destinationGrid;
    double localTimeStep;
    if (timeChoice == "half time update")
    {
        destinationGrid = std::unique_ptr<Grid1D>(&_gridHalfTime);
        localTimeStep = 0.5 * _timeStep;
    }
    else if (timeChoice == "full time update")
    {
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
        destinationGrid->density[i]  = grid.density[i]
                                      + (localTimeStep / _deltaX)
                                      * (_flux.density[i] - _flux.density[i+1]);

        destinationGrid->energy[i]   = grid.energy[i]
                                      + (localTimeStep / _deltaX)
                                      * (_flux.energy[i] - _flux.energy[i+1]);

        // Update momentum and magnetic field
        std::vector<double> deltas = {_deltaX, _deltaY, _deltaZ};
        for (size_t j = 0; j < 3; j++)
        {
            // Update momentum
            destinationGrid->momentum[i][j] = grid.momentum[i][j]
                                                + (localTimeStep / _deltaX)
                                                * (_flux.momentum[i][j] - _flux.momentum[i+1][j]);

            // Update magnetic field
            // Compute the modulos
            int j1 = _mod3(j+1), j2 = _mod3(j+2);
            std::vector<size_t> j2Offset(3, 0), j1Offset(3, 0);
            j1Offset[j1] = 1;
            j2Offset[j2] = 1;

            destinationGrid->magnetic[i][j] = grid.magnetic[i][j]
                                                + (localTimeStep / deltas[j2])
                                                * (_edgeFields[i+j2Offset[0]][1+j2Offset[1]][1+j2Offset[2]][j1] - _edgeFields[i][1][1][j1])
                                                - (localTimeStep / deltas[j1])
                                                * (_edgeFields[i+j1Offset[0]][1+j1Offset[1]][1+j1Offset[2]][j2] - _edgeFields[i][1][1][j2]);
        }
    }
    // Release pointer
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
      _magCentered(2 * _numGhosts + reals, stdVector1D (3, 0)),
      _interfaceL(reals, _numGhosts),
      _interfaceR(reals, _numGhosts),
      _flux(reals, _numGhosts),
      _gridHalfTime(reals, _numGhosts, boundaryConditions),
      _edgeFields(2 * _numGhosts + reals, stdVector3D(
                                       3, stdVector2D(
                                       3, stdVector1D(
                                       3, 0.0)))),
      _ctVelocities(2 * _numGhosts + reals, stdVector4D(
                                         3, stdVector3D(
                                         3, stdVector2D(
                                         3, stdVector1D(
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