/*!
 * \file HlldRiemannSolver.cpp
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Implementation of the HlldRiemannSolver class
 * \version 0.1
 * \date 2020-10-09
 *
 * \copyright Copyright (c) 2020
 *
 */

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "HlldRiemannSolver.h"
#include "mhdUtilities.h"

using namespace mhdUtilities;

// =============================================================================
void HlldRiemannSolver::riemannMain(double const &densityL,
                                    std::vector<double> const &velocityL,
                                    double const &pressureL,
                                    std::vector<double> const &magneticL,
                                    double const &densityR,
                                    std::vector<double> const &velocityR,
                                    double const &pressureR,
                                    std::vector<double> const &magneticR,
                                    double &densityFlux,
                                    std::vector<double> &momentumFlux,
                                    std::vector<double> &magneticFlux,
                                    double &energyFlux,
                                    double const &posOverT)
{
    // Copy arguments to member variables
    _densityL  = densityL;
    _velocityL = velocityL;
    _pressureL = std::max(pressureL, 1.0E-20);
    _pressureTotL = computeTotalPressure(_pressureL, magneticL);
    _magneticL = magneticL;
    _energyL   = computeEnergy(pressureL, _densityL, _velocityL, _magneticL, _gamma);
    _densityR  = densityR;
    _velocityR = velocityR;
    _pressureR = std::max(pressureR, 1.0E-20);
    _pressureTotR = computeTotalPressure(_pressureR, magneticR);
    _magneticR = magneticR;
    _energyR   = computeEnergy(pressureR, _densityR, _velocityR, _magneticR, _gamma);

    // Compute the wave speeds. This compute the _sL, _sM, and _sR approximate
    // wave speeds
    _computeWaveSpeeds();

    // Now we need to figure out which state we're in
    // According to athena this can be done as
    //   if (_sL >= 0.0) {
    //     return Fl if flow is supersonic
    //   } else if (_sR <= 0.0) {
    //     return Fr if flow is supersonic
    //   } else if (_sStarL >= 0.0) {
    //     return Fl*
    //   } else if (_sM >= 0.0) {
    //     return Fl**
    //   } else if (_sStarR > 0.0) {
    //     return Fr**
    //   } else {
    //     return Fr*
    //   }
    if (posOverT < _sL)
    {
        _computeStandardFluxes(_densityL, _velocityL, _pressureTotL, _magneticL, _energyL,
                               densityFlux, momentumFlux, magneticFlux, energyFlux);
    }
    else if ( (_sL <= posOverT) and (posOverT < _sStarL) )
    {
        _computeStarFluxes(_densityL, _velocityL, _pressureL, _pressureTotL, _magneticL, _energyL, _sL, _densityStarL,
                           densityFlux, momentumFlux, magneticFlux, energyFlux);
    }
    else if ( (_sStarL <= posOverT) and (posOverT < _sM) )
    {
        _computeDblStarFluxes(_magneticL, _sStarL, _densityStarL, -1.0,
                               densityFlux, momentumFlux, magneticFlux, energyFlux);
    }
    else if ( (_sM <= posOverT) and (posOverT < _sStarR) )
    {
        _computeDblStarFluxes(_magneticR, _sStarR, _densityStarR, 1.0,
                               densityFlux, momentumFlux, magneticFlux, energyFlux);
    }
    else if ( (_sStarR <= posOverT) and (posOverT <= _sR) )
    {
        _computeStarFluxes(_densityR, _velocityR, _pressureR, _pressureTotR, _magneticR, _energyR, _sR, _densityStarR,
                           densityFlux, momentumFlux, magneticFlux, energyFlux);
    }
    else if (_sR < posOverT)
    {
        _computeStandardFluxes(_densityR, _velocityR, _pressureTotR, _magneticR, _energyR,
                               densityFlux, momentumFlux, magneticFlux, energyFlux);
    }
    else
    {
        throw std::runtime_error("Error in HLLD Riemann Solver: No state chosen");
    }
}
// =============================================================================


// =============================================================================
void HlldRiemannSolver::_computeStandardFluxes(double const &density,
                                               std::vector<double> const &velocity,
                                               double const &pressureTot,
                                               std::vector<double> const &magnetic,
                                               double const &energy,
                                               double &densityFlux,
                                               std::vector<double> &momentumFlux,
                                               std::vector<double> &magneticFlux,
                                               double &energyFlux)
{
    // Compute the regular fluxes for the left or right states
    // Note that a lot of these terms are momentum and so I shouldn't recompute
    // them
    densityFlux     = density  * velocity[0];

    momentumFlux[0] = density  * velocity[0] * velocity[0] + pressureTot - magnetic[0] * magnetic[0];
    momentumFlux[1] = density  * velocity[0] * velocity[1] - magnetic[0] * magnetic[1];
    momentumFlux[2] = density  * velocity[0] * velocity[2] - magnetic[0] * magnetic[2];

    magneticFlux[0] = 0.;
    magneticFlux[1] = magnetic[1] * velocity[0] - magnetic[0] * velocity[1];
    magneticFlux[2] = magnetic[2] * velocity[0] - magnetic[0] * velocity[2];

    // Group transverse terms for FP associative symmetry
    energyFlux = velocity[0] * (energy + pressureTot) - magnetic[0]
                 * (velocity[0] * magnetic[0] + ((velocity[1] * magnetic[1])  + (velocity[2] * magnetic[2])));

    // Set member variables to the current state for retrieval if needed
    _densityState     = density;
    _velocityState    = velocity;
    _magneticState    = magnetic;
    _pressureState    = pressureTot - 0.5 * (magnetic[0] * magnetic[0] + ((magnetic[1] * magnetic[1])  + (magnetic[2] * magnetic[2])));
    _pressureTotState = pressureTot;
    _energyState      = energy;
}
// =============================================================================


// =============================================================================
void HlldRiemannSolver::_computeVandBStar(double const &density,
                                          std::vector<double> const &velocity,
                                          double const &pressure,
                                          std::vector<double> const &magnetic,
                                          double const &sSide,
                                          std::vector<double> &velocityStar,
                                          std::vector<double> &magneticStar)
{
    // Set known quantities
    velocityStar[0] = _sM;
    magneticStar[0] = magnetic[0];

    // Compute the magnetosonic speed
    double magSonicSide = magnetosonicSpeed(pressure, density, magnetic, _gamma);

    // Check for divide-by-zero errors
    _divZero = (_sM == velocity[0])
                and ((sSide == velocity[0] + magSonicSide) or (sSide == velocity[0] - magSonicSide))
                and (magnetic[0] * magnetic[0] >= _gamma * pressure)
                and (magnetic[1] == 0.)
                and (magnetic[2] == 0.);
    if ( _divZero)
    {
        // This is the divide-by-zero case
        velocityStar = velocity;
        magneticStar[1] = 0.;
        magneticStar[2] = 0.;
        std::cout << "Degenerate case" << std::endl;
    }
    else
    {
        // This is the standard case where we have no divide-by-zero errors

        // Compute the denominator that is used for both velocity and magnetic
        // field computations
        double const denom = density * (sSide - velocity[0]) * (sSide - _sM) - (magnetic[0] * magnetic[0]);

        // Compute the velocity and magnetic field in the star state
        double coef = magnetic[0]  * (_sM - velocity[0]) / denom;
        velocityStar[1] = velocity[1] - magnetic[1] * coef;
        velocityStar[2] = velocity[2] - magnetic[2] * coef;

        double tmpPower = (sSide - velocity[0]);
        tmpPower = tmpPower * tmpPower;
        coef = (density * tmpPower - (magnetic[0] * magnetic[0])) / denom;
        magneticStar[1] = magnetic[1] * coef;
        magneticStar[2] = magnetic[2] * coef;
    }
}
// =============================================================================


// =============================================================================
void HlldRiemannSolver::_computeStarFluxes(double const &density,
                                           std::vector<double> const &velocity,
                                           double const &pressure,
                                           double const &pressureTot,
                                           std::vector<double> const &magnetic,
                                           double const &energy,
                                           double const &sSide,
                                           double const &densityStarSide,
                                           double &densityFlux,
                                           std::vector<double> &momentumFlux,
                                           std::vector<double> &magneticFlux,
                                           double &energyFlux)
{
    // Compute the star state (between the fast magnetosonic wave and the Alfven wave)

    // We already know some of the star state values
    //  1. density has already been found and is in densityStarSide
    //  2. v_x^* is S_M

    // Declare all the star state variables
    std::vector<double> velocityStar(3), magneticStar(3);
    double densityStar, pressureTotStar, energyStar;

    // Compute the velocity and magnetic field in the star region
    _computeVandBStar(density, velocity, pressure, magnetic, sSide,
                      velocityStar,magneticStar);

    // Check for divide-by-zero errors, the _divZero member variables was
    // already compute in _computeVandBStar
    if (_divZero)
    {
        // This is the divide-by-zero case
        densityStar = density;
        pressureTotStar = pressureTot;
    }
    else
    {
        // This is the standard case where we have no divide-by-zero errors
        densityStar = densityStarSide;  // Density in the star state
        // pressureTotStar =  // Total pressure in the star state
        //     ( _densityR * _pressureTotL * (_sR - _velocityR[0])
        //     - _densityL * _pressureTotR * (_sL - _velocityL[0])
        //     + _densityL * _densityR
        //         * (_sR - _velocityR[0])
        //         * (_sL - _velocityL[0])
        //         * (_velocityR[0] - _velocityL[0]))
        //     /
        //     (_densityR * (_sR - _velocityR[0]) - _densityL * (_sL - _velocityL[0]));

        // Athena version of the equation. They claim that eq 41 in MK2005 becomes eq 23
        double pressureTotStarL = _pressureTotL + _densityL * (_sL - _velocityL[0]) * (_sM - _velocityL[0]);
        double pressureTotStarR = _pressureTotR + _densityR * (_sR - _velocityR[0]) * (_sM - _velocityR[0]);
        pressureTotStar  = 0.5 * (pressureTotStarL + pressureTotStarR);
    }

    // Now we'll compute the energy in the state
    energyStar =
        ( energy * (sSide - velocity[0])
        - pressureTot * velocity[0]
        + pressureTotStar * _sM
        + magnetic[0] * ((velocity[0] * magnetic[0] + ((velocity[1] * magnetic[1])  + (velocity[2] * magnetic[2])))
                      -  (velocityStar[0] * magneticStar[0] + ((velocityStar[1] * magneticStar[1])  + (velocityStar[2] * magneticStar[2])))))
        / (sSide - _sM);

    // Compute the standard flux
    double densityFluxSide, energyFluxSide;
    std::vector<double> momentumFluxSide(3), magneticFluxSide(3);
    _computeStandardFluxes(density, velocity, pressureTot, magnetic, energy,
                           densityFluxSide, momentumFluxSide, magneticFluxSide, energyFluxSide);

    // Now we actually compute the HLLD Star State Flux
    for (size_t i = 0; i < 3; i++)
    {
        momentumFlux[i] = momentumFluxSide[i] + sSide * (computeMomentum(velocityStar[i], densityStar) - computeMomentum(velocity[i], density));
        magneticFlux[i] = magneticFluxSide[i] + sSide * (magneticStar[i] - magnetic[i]);
    }
    densityFlux = densityFluxSide + sSide * (densityStar - density);
    energyFlux  = energyFluxSide  + sSide * (energyStar  - energy);
    magneticFlux[0] = 0.0;  // The flux in the x-direction is always zero

    // Set member variables to the current state for retrieval if needed
    _densityState     = densityStar;
    _velocityState    = velocityStar;
    _magneticState    = magneticStar;
    _pressureState    = pressureTotStar - 0.5 * (magnetic[0] * magnetic[0] + ((magnetic[1] * magnetic[1])  + (magnetic[2] * magnetic[2])));
    _pressureTotState = pressureTotStar;
    _energyState      = energyStar;
}
// =============================================================================


// =============================================================================
void HlldRiemannSolver::_computeDblStarFluxes(std::vector<double> const &magnetic,
                                              double const &sStarSide,
                                              double const &densityStarSide,
                                              double const &sideSign,
                                              double &densityFlux,
                                              std::vector<double> &momentumFlux,
                                              std::vector<double> &magneticFlux,
                                              double &energyFlux)
{
    // Compute the double star state, the region between the Alfven wave and the
    // contact discontinuity

    // We already know some of the star state values
    //  1. densityDblStarSide = densityStarSide
    //  2. v_x^** = v_x^* = S_M
    //  3. pressureTotDblStar = pressureTotStar
    //  4. B_x^** = B_x
    // All that is left is v_y, v_z, B_y, B_z, and energy

    // Lets declare all our variables to store the double star state
    double densityDblStar, pressureTotDblStar, energyDblStar;
    std::vector<double> velocityDblStar(3), magneticDblStar(3);

    // We will also need to have all the values for both star states computed to
    // compute the double star state. Once the _computeStarFluxes function is
    // finished then the star state will be stored in the _"primitive"Star
    // variables and we can use them directly.
    double densityFluxStarSide, energyFluxStarSide;
    std::vector<double> velocityStarL(3), magneticStarL(3),
                        velocityStarR(3), magneticStarR(3),
                        momentumFluxStarSide(3), magneticFluxStarSide(3);
    // todo:    instead of this complicated bs try calling multiple functions,
    // todo:    one for each component in the IF statement where I choose which
    // todo:    state I am in
    if (sideSign == 1.0)
    {
        // We're on the right side

        // Compute the off (left) side velocity and magnetic field
        _computeVandBStar(_densityL, _velocityL, _pressureL, _magneticL, _sL,
                          velocityStarL, magneticStarL);

        // Compute the on (right) side state and flux
        _computeStarFluxes(_densityR, _velocityR, _pressureR, _pressureTotR, _magneticR, _energyR, _sR, _densityStarR,
                           densityFluxStarSide, momentumFluxStarSide, magneticFluxStarSide, energyFluxStarSide);
        velocityStarR = _velocityState;
        magneticStarR = _magneticState;
    }
    else if (sideSign == -1.0)
    {
        // We're on the left side
        // Compute the off (right) side velocity and magnetic field
        _computeVandBStar(_densityR, _velocityR, _pressureR, _magneticR, _sR,
                          velocityStarR, magneticStarR);

        // Compute the on (left) side state and flux
        _computeStarFluxes(_densityL, _velocityL, _pressureL, _pressureTotL, _magneticL, _energyL, _sL, _densityStarL,
                           densityFluxStarSide, momentumFluxStarSide, magneticFluxStarSide, energyFluxStarSide);
        velocityStarL = _velocityState;
        magneticStarL = _magneticState;
    }
    else
    {
        throw std::runtime_error("Error in HLLD Riemann Solver: Incorrect value for sideSign in HlldRiemannSolver::_computeDblStarFluxes");
    }

    // Now that we have all the prerequisites we can start computing the double
    // star state starting with the known quantities
    densityDblStar     = densityStarSide;
    velocityDblStar[0] = _sM;
    pressureTotDblStar = _pressureTotState;
    magneticDblStar[0] = magnetic[0];

    // There are several terms that are common to many double star state
    // computations so we'll assign variables to them.
    double const sqrtDenStarL  = std::sqrt(_densityStarL);
    double const sqrtDenStarR  = std::sqrt(_densityStarR);
    double const signMagneticX = (magnetic[0] >= 0.)? 1.0: -1.0;
    double const denom         = sqrtDenStarL + sqrtDenStarR;

    // Compute the double star state velocity and magnetic field
    for (size_t i = 1; i < 3; i++)
    {
        velocityDblStar[i] =
            (velocityStarL[i] * sqrtDenStarL
            + velocityStarR[i] * sqrtDenStarR
            + signMagneticX * (magneticStarR[i] - magneticStarL[i]))
            / denom;

        magneticDblStar[i] =
            (magneticStarR[i] * sqrtDenStarL
            + magneticStarL[i] * sqrtDenStarR
            + signMagneticX * sqrtDenStarL * sqrtDenStarR * (velocityStarR[i] - velocityStarL[i]))
            / denom;
    }

    // Compute the double star state energy
    double sqrtDenSide  = (sideSign > 0.)? sqrtDenStarR: sqrtDenStarL;
    std::vector<double> *velocityStarSide = (sideSign > 0.)? &velocityStarR: &velocityStarL;
    std::vector<double> *magneticStar = (sideSign > 0.)? &magneticStarR: &magneticStarL;

    energyDblStar = _energyState + sideSign * sqrtDenSide * signMagneticX
        * (((*velocityStarSide)[0] * (*magneticStar)[0] + (((*velocityStarSide)[1] * (*magneticStar)[1])  + ((*velocityStarSide)[2] * (*magneticStar)[2])))
        - (velocityDblStar[0] * magneticDblStar[0] + ((velocityDblStar[1] * magneticDblStar[1])  + (velocityDblStar[2] * magneticDblStar[2]))));

    // Compute the double star state HLLD Fluxes
    densityFlux = densityFluxStarSide + sStarSide * (densityDblStar - densityStarSide);
    energyFlux = energyFluxStarSide + sStarSide * (energyDblStar - _energyState);
    for (size_t i = 0; i < 3; i++)
    {
        momentumFlux[i] = momentumFluxStarSide[i] + sStarSide * (computeMomentum(velocityDblStar[i], densityDblStar) - computeMomentum(_velocityState[i], _densityState));
        magneticFlux[i] = magneticFluxStarSide[i] + sStarSide * (magneticDblStar[i] - _magneticState[i]);
    }
    magneticFlux[0] = 0.0;

    // Set member variables to the current state for retrieval if needed
    _densityState     = densityDblStar;
    _velocityState    = velocityDblStar;
    _magneticState    = magneticDblStar;
    _pressureState    = pressureTotDblStar - 0.5 * (magnetic[0] * magnetic[0] + ((magnetic[1] * magnetic[1])  + (magnetic[2] * magnetic[2])));
    _pressureTotState = pressureTotDblStar;
    _energyState      = energyDblStar;
}
// =============================================================================


// =============================================================================
void HlldRiemannSolver::_computeWaveSpeeds()
{
    // Compute the fast magnetosonic wave speeds
    double magSonicL, magSonicR;  // The speeds of the left and right fast magnetosonic waves
    magSonicL = magnetosonicSpeed(_pressureL, _densityL, _magneticL, _gamma);
    magSonicR = magnetosonicSpeed(_pressureR, _densityR, _magneticR, _gamma);

    // Compute the S_L and S_R wave speeds.
    _sL = std::min(_velocityL[0], _velocityR[0]) - std::max(magSonicL, magSonicR);
    _sR = std::max(_velocityL[0], _velocityR[0]) + std::max(magSonicL, magSonicR);

    // Athena's version
    // _sL = std::min(_velocityL[0] - magSonicL, _velocityR[0] - magSonicR);
    // _sR = std::max(_velocityL[0] + magSonicL, _velocityR[0] + magSonicR);


    // Compute the S_M wave speed, note that the density*velocity terms are just
    // momentum, might want to pass that in. Also the velocity differences are
    // used twice in _sM and once in _densityStarK where K=L,R
    _sM = // Numerator
          ( _densityR * _velocityR[0] * (_sR - _velocityR[0])
          - _densityL * _velocityL[0] * (_sL - _velocityL[0])
          + (_pressureTotL - _pressureTotR))
          /
          // Denominator
          ( _densityR * (_sR - _velocityR[0])
          - _densityL * (_sL - _velocityL[0]));

    // Compute the densities in the star state
    _densityStarL = _densityL * (_sL - _velocityL[0]) / (_sL - _sM);
    _densityStarR = _densityR * (_sR - _velocityR[0]) / (_sR - _sM);

    // Compute the S_L^* and S_R^* wave speeds
    _sStarL = _sM - std::abs(_magneticL[0]) / std::sqrt(_densityStarL);
    _sStarR = _sM + std::abs(_magneticR[0]) / std::sqrt(_densityStarR);
}
// =============================================================================