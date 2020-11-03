/*!
 * \file HllcRiemannSolver.cpp
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Implementation of the HllcRiemannSolver class
 * \version 0.1
 * \date 2020-10-09
 *
 * \copyright Copyright (c) 2020
 *
 */

#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "HllcRiemannSolver.h"
#include "mhdUtilities.h"

using namespace mhdUtilities;

// =============================================================================
void HllcRiemannSolver::riemannMain(double const &densityL,
                                    double const &velocityL,
                                    double const &pressureL,
                                    double const &densityR,
                                    double const &velocityR,
                                    double const &pressureR,
                                    double &densityFlux,
                                    double &momentumFlux,
                                    double &energyFlux,
                                    double const &posOverT)
{
    // Copy arguments to member variables
    _densityL  = densityL;
    _velocityL = velocityL;
    _pressureL = std::max(pressureL, 1.0E-20);
    _energyL   = computeEnergy(_pressureL, _densityL, _velocityL, _gamma);
    _densityR  = densityR;
    _velocityR = velocityR;
    _pressureR = std::max(pressureR, 1.0E-20);
    _energyR   = computeEnergy(_pressureR, _densityR, _velocityR, _gamma);

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
        throw std::runtime_error("Error in HLLC Riemann Solver: No state chosen");
    }
}
// =============================================================================


// =============================================================================
void HllcRiemannSolver::_computeStandardFluxes(double const &density,
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
void HllcRiemannSolver::_computeStarFluxes(double const &density,
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
void HllcRiemannSolver::_computeWaveSpeeds()
{
    // Square root of density ratios
    double denRatio = std::sqrt(_densityR/_densityL);

    // Compute the enthalpies
    double hL = (_energyL + _pressureL) / _densityL;
    double hR = (_energyR + _pressureR) / _densityR;

    double hTilde = (hL + hR * denRatio) / (1 + denRatio);

    // Compute the tilde velocity and sound speed
    double velTilde = (_velocityL + _velocityR * denRatio) / (1 + denRatio);
    double cTilde = std::sqrt( (_gamma - 1)
                               * (hTilde - 0.5 * std::pow(velTilde, 2)) );

    // Compute the S_L and S_R wave speeds
    _sL = std::min( _velocityL - soundSpeed(_pressureL, _densityL, _gamma),
                    velTilde - cTilde);
    _sR = std::max( _velocityR + soundSpeed(_pressureR, _densityR, _gamma),
                    velTilde + cTilde);

    // Compute the S_M wave speed
    _sM = // Numerator
          ( _densityR * _velocityR * (_sR - _velocityR)
          - _densityL * _velocityL * (_sL - _velocityL)
          + _pressureL - _pressureR)
          /
          // Denominator
          ( _densityR * (_sR - _velocityR)
          - _densityL * (_sL - _velocityL));
}
// =============================================================================