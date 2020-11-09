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

#include "HlldRiemannSolver.h"
#include "mhdUtilities.h"

using namespace mhdUtilities;

// =============================================================================
void HlldRiemannSolver::riemannMain(double const &densityL,
                                    std::vector<double> const &velocityL,
                                    double const &pressureL,
                                    std::vector<double> magneticL,
                                    double const &densityR,
                                    std::vector<double> const &velocityR,
                                    double const &pressureR,
                                    std::vector<double> magneticR,
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
    if (posOverT < _sL)
    {
        _computeStandardFluxes(_densityL, _velocityL, _pressureTotL, _magneticL, _energyL,
                               densityFlux, momentumFlux, magneticFlux, energyFlux);
    }
    else if ( (_sL <= posOverT) and (posOverT <= _sStarL) )
    {
        _computeStarFluxes(_densityL, _velocityL, _pressureL, _pressureTotL, _magneticL, _energyL, _sL, _densityStarL,
                           densityFlux, momentumFlux, magneticFlux, energyFlux);
    }
    else if ( (_sStarL <= posOverT) and (posOverT <= _sM) )
    {
        _computeStarStarFluxes(_densityL, _velocityL, _pressureTotL, _magneticL, _energyL, _sStarL,
                               densityFlux, momentumFlux, magneticFlux, energyFlux);
    }
    else if ( (_sM <= posOverT) and (posOverT <= _sStarR) )
    {
        _computeStarStarFluxes(_densityR, _velocityR, _pressureTotR, _magneticR, _energyR, _sStarR,
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
    densityFlux     = density  * velocity[0];

    momentumFlux[0] = density  * velocity[0] * velocity[0] + pressureTot - magnetic[0] * magnetic[0];
    momentumFlux[1] = density  * velocity[0] * velocity[1] - magnetic[0] * magnetic[1];
    momentumFlux[3] = density  * velocity[0] * velocity[3] - magnetic[0] * magnetic[3];

    magneticFlux[0] = 0.;
    magneticFlux[1] = magnetic[1] * velocity[0] - magnetic[0] * velocity[1];
    magneticFlux[2] = magnetic[2] * velocity[0] - magnetic[0] * velocity[2];

    energyFlux = velocity[0] * (energy + pressureTot) - magnetic[0]
                 * (std::inner_product(velocity.begin(), velocity.end(), magnetic.begin(), 0));

    // Set member variables to the current state for retrieval if needed
    _densityState  = density;
    _velocityState = velocity;
    _magneticState = magnetic;
    _pressureState = pressureTot - 0.5 * std::inner_product(magnetic.begin(), magnetic.end(), magnetic.begin(), 0);
    _pressureTotState = pressureTot;
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
    // First we must compute the standard flux
    double densityFluxSide, energyFluxSide;
    std::vector<double> momentumFluxSide(3), magneticFluxSide(3);
    _computeStandardFluxes(density, velocity, pressureTot, magnetic, energy,
                           densityFluxSide, momentumFluxSide, magneticFluxSide, energyFluxSide);

    // At this point the various State member variables are the state in the L
    // or R state so we can use them to compute the L* or R* state

    // Compute the star state (between the fast magnetosonic wave and the Alfven wave)
        // density has already been found and is in densityStarSide
        // v_x^* is S_M
    // Compute the velocity and magnetic field
    std::vector<double> velocityStar(3), magneticStar(3);
    double densityStar, pressureTotStar, energyStar;
    velocityStar[0] = _sM;
    magneticStar[0] = magnetic[0];

    double magSonicSide = magnetosonicSpeed(pressure, density, magnetic, _gamma);

    // Check for divide by zero errors
    if ( (_sM == velocity[0])
         and ((sSide == velocity[0] + magSonicSide) or (sSide == velocity[0] - magSonicSide))
         and (magnetic[0] * magnetic[0] >= _gamma * pressure)
         and (magnetic[1] == 0.)
         and (magnetic[2] == 0.))
    {
        // This is the divide by zero case
        densityStar = density;
        velocityStar = velocity;
        pressureTotStar = pressureTot;
        magneticStar[1] = 0.;
        magneticStar[2] = 0.;
    }
    else
    {
        densityStar = densityStarSide;
        pressureTotStar =
            ( _densityR * _pressureTotL * (_sR - _velocityR[0])
            - _densityL * _pressureTotR * (_sL - _velocityL[0])
            + _densityL * _densityR
                * (_sR - _velocityR[0])
                * (_sL - _velocityL[0])
                * (_velocityR[0] - _velocityL[0]))
            /
            (_densityR * (_sR - _velocityR[0]) - _densityL * (_sL - _velocityL[0]));

        double denom = density * (sSide - velocity[0]) * (sSide - _sM) - magnetic[0] * magnetic[0];

        velocityStar[1] = velocity[1] - magnetic[0] * magnetic[1] * (_sM - velocity[0]) / denom;
        velocityStar[2] = velocity[2] - magnetic[0] * magnetic[2] * (_sM - velocity[0]) / denom;
        magneticStar[1] = magnetic[1] * (density * std::pow(sSide - velocity[0], 2) - magnetic[0] * magnetic[0]) / denom;
        magneticStar[2] = magnetic[2] * (density * std::pow(sSide - velocity[0], 2) - magnetic[0] * magnetic[0]) / denom;
    }
    // Now we'll work on computing the pressure and energy in the state
    energyStar =
        ( energy * (sSide - velocity[0])
        - pressureTot * velocity[0]
        + pressureTotStar * _sM
        + magnetic[0] * (std::inner_product(velocity.begin(), velocity.end(), magnetic.begin(), 0)
                    - std::inner_product(velocityStar.begin(), velocityStar.end(), magneticStar.begin(), 0)))
        / (sSide - _sM);

    // Now we actually compute the HLLD Flux. Note that the _xxxxState variables
    // are all from the standard flux function
    densityFlux = densityFluxSide + sSide * (densityStar - density);
    energyFlux  = energyFluxSide  + sSide * (energyStar  - energy);
    for (size_t i = 0; i < 3; i++)
    {
        momentumFlux[i] = momentumFluxSide[i] + sSide * (velocityStar[i] - densityStar);
        magneticFlux[i] = magneticFluxSide[i] + sSide * (magneticStar[i] - magnetic[i]);
    }
    magneticFlux[0] = 0.0;

    // Set member variables to the current state for retrieval if needed
    _densityState  = densityStar;
    _velocityState = velocityStar;
    _magneticState = magneticStar;
    _pressureState = pressureTotStar - 0.5 * std::inner_product(magneticStar.begin(), magneticStar.end(), magneticStar.begin(), 0);
    _pressureTotState = pressureTotStar;
}
// =============================================================================


// =============================================================================
void HlldRiemannSolver::_computeWaveSpeeds()
{
    // Compute the fast magnetosonic wave speeds
    double magSonicL, magSonicR;  // The speeds of the left and right fast magnetosonic waves
    magSonicL = magnetosonicSpeed(_pressureL, _densityL, _magneticL, _gamma);
    magSonicR = magnetosonicSpeed(_pressureR, _densityR, _magneticR, _gamma);

    // Compute the S_L and S_R wave speeds
    _sL = std::min(_velocityL[0], _velocityR[0]) - std::max(magSonicL, magSonicR);
    _sL = std::max(_velocityL[0], _velocityR[0]) + std::max(magSonicL, magSonicR);

    // Compute the S_M wave speed
    _sM = // Numerator
          ( _densityR * _velocityR[0] * (_sR - _velocityR[0])
          - _densityL * _velocityL[0] * (_sL - _velocityL[0])
          + _pressureTotL - _pressureTotR)
          /
          // Denominator
          ( _densityR * (_sR - _velocityR[0])
          - _densityL * (_sL - _velocityL[0]));

    // Compute the densities in the star state
    _densityStarL = _densityL * (_sL - _velocityL[0]) / (_sL - _sM);
    _densityStarL = _densityR * (_sR - _velocityR[0]) / (_sR - _sM);

    // Compute the S_L^* and S_R^* wave speeds
    _sStarL = _sM - std::abs(_magneticL[0]) / std::sqrt(_densityStarL);
    _sStarL = _sM + std::abs(_magneticR[0]) / std::sqrt(_densityStarR);
}
// =============================================================================