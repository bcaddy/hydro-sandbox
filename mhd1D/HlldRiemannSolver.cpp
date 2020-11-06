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
                                    double &momentumFlux,
                                    double &energyFlux,
                                    double const &posOverT)
{
    // Copy arguments to member variables
    _densityL  = densityL;
    _velocityL = velocityL;
    _pressureL = std::max(pressureL, 1.0E-20);
    _magneticL = magneticL;
    _energyL   = computeEnergy(_pressureL, _densityL, _velocityL, _magneticL, _gamma);
    _densityR  = densityR;
    _velocityR = velocityR;
    _pressureR = std::max(pressureR, 1.0E-20);
    _magneticR = magneticR;
    _energyR   = computeEnergy(_pressureR, _densityR, _velocityR, _magneticR, _gamma);

    // Compute the wave speeds. This compute the _sL, _sM, and _sR approximate
    // wave speeds
    _computeWaveSpeeds();

    // Now we need to figure out which state we're in
    if (posOverT < _sL)
    {
        _computeStandardFluxes(_densityL, _velocityL, _pressureL, _energyL,
                               densityFlux, momentumFlux, energyFlux);
    }
    else if ( (_sL <= posOverT) and (posOverT <= _sM) )
    {
        _computeStarFluxes(_densityL, _velocityL, _pressureL, _energyL, _sL,
                               densityFlux, momentumFlux, energyFlux);
    }
    else if ( (_sM <= posOverT) and (posOverT <= _sR) )
    {
        _computeStarFluxes(_densityR, _velocityR, _pressureR, _energyR, _sR,
                               densityFlux, momentumFlux, energyFlux);
    }
    else if (_sR < posOverT)
    {
        _computeStandardFluxes(_densityR, _velocityR, _pressureR, _energyR,
                               densityFlux, momentumFlux, energyFlux);
    }
    else
    {
        throw std::runtime_error("Error in HLLD Riemann Solver: No state chosen");
    }
}
// =============================================================================


// =============================================================================
void HlldRiemannSolver::_computeStandardFluxes(double const &density,
                                               double const &velocity,
                                               double const &pressure,
                                               double const &energy,
                                               double &densityFlux,
                                               double &momentumFlux,
                                               double &energyFlux)
{
    // Compute the regular fluxes for the left or right states
    densityFlux  = density  * velocity;
    momentumFlux = density  * std::pow(velocity, 2) + pressure;
    energyFlux   = velocity * (energy + pressure);

    // Set member variables to the current state for retrieval if needed
    _densityState  = density;
    _velocityState = velocity;
    _pressureState = pressure;
}
// =============================================================================


// =============================================================================
void HlldRiemannSolver::_computeStarFluxes(double const &density,
                                           double const &velocity,
                                           double const &pressure,
                                           double const &energy,
                                           double const &sSide,
                                           double &densityFlux,
                                           double &momentumFlux,
                                           double &energyFlux)
{
    // Compute the state in the star region
    double densityStar  = density * (sSide - velocity) / (sSide - _sM);
    double pressureStar = _pressureL
                          + _densityL * (_velocityL - _sL) * (_velocityL - _sM);
    double momentumStar = (
                            density * velocity * (sSide - velocity)
                            + (pressureStar - pressure)
                          )
                          /
                          (sSide - _sM);
    double energyStar   = (
                            energy * (sSide - velocity)
                            - pressure * velocity
                            + pressureStar * _sM
                          )
                          /
                          (sSide - _sM);

    // Compute the standard flux
    double sideDensityFlux, sideMomentumFlux, sideEnergyFlux;
    _computeStandardFluxes(density, velocity, pressure, energy,
                           sideDensityFlux, sideMomentumFlux, sideEnergyFlux);

    // Compute and return the fluxes
    densityFlux  = sideDensityFlux  + sSide * (densityStar  - density);
    momentumFlux = sideMomentumFlux + sSide * (momentumStar - computeMomentum(velocity, density));
    energyFlux   = sideEnergyFlux   + sSide * (energyStar   - energy);

    // Set member variables to the current state for retrieval if needed
    _densityState  = densityStar;
    _velocityState = momentumStar/densityStar;
    _pressureState = pressureStar;
}
// =============================================================================


// =============================================================================
void HlldRiemannSolver::_computeWaveSpeeds()
{
    // Compute the fast magnetosonic wave speeds
    double magSonicL, magSonicR;  // The speeds of the left and right fast magnetosonic waves
    double magSonicL = magnetosonicSpeed(_pressureL, _densityL, _magneticL, _gamma);
    double magSonicR = magnetosonicSpeed(_pressureR, _densityR, _magneticR, _gamma);

    // Compute the S_L and S_R wave speeds
    _sL = std::min(_velocityL[0], _velocityR[0]) - std::max(magSonicL, magSonicR);
    _sL = std::max(_velocityL[0], _velocityR[0]) + std::max(magSonicL, magSonicR);

    // Compute the S_M wave speed
    _sM = // Numerator
          ( _densityR * _velocityR[0] * (_sR - _velocityR[0])
          - _densityL * _velocityL[0] * (_sL - _velocityL[0])
          + _pressureL - _pressureR)
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