/*!
 * \file RiemannSolver.cpp
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Implementation of the RiemannSolver class
 * \version 0.1
 * \date 2020-07-14
 * 
 * \copyright Copyright (c) 2020
 * 
 */

#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "RiemannSolver.h"

// =============================================================================
void RiemannSolver::riemannMain(double const &densityR,
                                double const &velocityR,
                                double const &pressureR,
                                double const &densityL,
                                double const &velocityL,
                                double const &pressureL,
                                double const &position,
                                double const &time,
                                double const &energyFlux,
                                double const &momentumFlux,
                                double const &massFlux)
{
    // Copy arguments to member variables
    _densityR     = densityR;
    _velocityR    = velocityR;
    _pressureR    = pressureR;
    _densityL     = densityL;
    _velocityL    = velocityL;
    _pressureL    = pressureL;
    _position     = position;
    _time         = time;
    _energyFlux   = energyFlux;
    _momentumFlux = momentumFlux;
    _massFlux     = massFlux;


    // Compute the sound speeds
    _cR = std::sqrt(_gamma * pressureR / densityR);
    _cL = std::sqrt(_gamma * pressureL / densityL);

    // Check for a vacuum
    if ((2 / (_gamma - 1)) * (_cL + _cR) >= (_velocityR - _velocityL))
    {
        throw std::runtime_error("Vacuum Detected. Exiting.");
    }

    // Figure out the pressure in the star region. This requires a good initial
    // guess then refining that guess with Newton-Raphson iterations
    _pressureStar = _computePressureStar();
}
// =============================================================================

// =============================================================================
double RiemannSolver::_computePressureStar()
{
    ;
}
// =============================================================================

// =============================================================================
double RiemannSolver::_guessPressureStar()
{
    // Compute min and maximum pressures of the two sides of the interface
    double pMin = std::min(_pressureL, _pressureR);
    double pMax = std::max(_pressureL, _pressureR);

    // First compute the primitive variable approximation
    double pPrim = 0.5   * (_pressureL + _pressureR) -
                   0.125 * (_velocityR - _velocityL) 
                         * (_densityL + _densityR) 
                         * (_cL + _cR);
    // Make sure it's not negative
    pPrim = std::max(_tol, pPrim);

    // Check to see if we should use the primitive variable approximation or not
    // TODO Stopped at cell with equation 9
    if ( ((pMax/pMin) <= 2.0) && (pMin < ) && () )
    {
        /* code */
    }
    
}
// =============================================================================

// =============================================================================
// Constructor
RiemannSolver::RiemannSolver(double const &gamma)

    // Start by initializing all the const member variables
    : _gamma(gamma)
{
}
// =============================================================================