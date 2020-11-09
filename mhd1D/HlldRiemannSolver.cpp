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
        _computeStarFluxes(_densityL, _velocityL, _pressureTotL, _magneticL, _energyL, _sL,
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
        _computeStarFluxes(_densityR, _velocityR, _pressureTotR, _magneticR, _energyR, _sR,
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
// =============================================================================


// =============================================================================
void HlldRiemannSolver::_computeStarFluxes(double const &density,
                                           std::vector<double> const &velocity,
                                           double const &pressureTot,
                                           std::vector<double> const &magnetic,
                                           double const &energy,
                                           double const &sSide,
                                           double &densityFlux,
                                           std::vector<double> &momentumFlux,
                                           std::vector<double> &magneticFlux,
                                           double &energyFlux)
{
    // First we must compute the standard flux
    double sideDensityFlux, sideEnergyFlux;
    std::vector<double> sideMomentumFlux(3), sideMagneticFlux(3);
    _computeStandardFluxes(density, velocity, pressureTot, magnetic, energy,
                           sideDensityFlux, sideMomentumFlux, sideMagneticFlux, sideEnergyFlux);

    // At this point the various State member variables are the state in the L
    // or R state so we can use them to compute the L* or R* state

    // Compute the star state (between the fast magnetosonic wave and the Alfven wave)
    



    // Compute and return the fluxes
    /// \todo update this
    // densityFlux  = ;//??
    // momentumFlux = ;//??
    // energyFlux   = ;//??

    // Set member variables to the current state for retrieval if needed
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