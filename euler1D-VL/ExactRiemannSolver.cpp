/*!
 * \file ExactRiemannSolver.cpp
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Implementation of the ExactRiemannSolver class
 * \version 0.1
 * \date 2020-07-14
 *
 * \copyright Copyright (c) 2020
 *
 */

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "HydroUtilities.h"
#include "ExactRiemannSolver.h"

using namespace HydroUtilities;

// =============================================================================
void ExactRiemannSolver::riemannMain(double const &densityL,
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
    _densityL     = densityL;
    _velocityL    = velocityL;
    _pressureL    = std::max(pressureL, 1.0E-20);
    _densityR     = densityR;
    _velocityR    = velocityR;
    _pressureR    = std::max(pressureR, 1.0E-20);
    _posOverT     = posOverT;

    // Compute the sound speeds
    _cL = soundSpeed(pressureL, densityL, _gamma);
    _cR = soundSpeed(pressureR, densityR, _gamma);

    // Check for a vacuum
    if ((2 / (_gamma - 1)) * (_cL + _cR) <= (_velocityR - _velocityL))
    {
        throw std::runtime_error("Vacuum Detected. Exiting.");
    }

    // Figure out the pressure in the star region. This requires a good initial
    // guess then refining that guess with Newton-Raphson iterations
    _pressureStar = _computePressureStar();

    // Compute the velocity in the star region
    {
        double fL = 0, fR = 0, df = 0;
        _pressureFunctions(_pressureStar, _pressureL, _densityL, _cL, fL, df);
        _pressureFunctions(_pressureStar, _pressureR, _densityR, _cR, fR, df);
        _velocityStar = 0.5 * (_velocityL + _velocityR + fR - fL);
    }

    // Determine if we're in the L/L_* state or the R/R_* state
    if (_velocityStar >= _posOverT)
    {
        // We're in the L or L_* state
        if (_pressureStar > _pressureL)
        {
            // The Left non-linear wave is a shockwave
            double shockSpeed = _shockSpeed("left");

            // Decide between L and L_* state
            if (shockSpeed >= _posOverT)
            {
                // We're in the L state
                _pressureState = _pressureL;
                _velocityState = _velocityL;
                _densityState  = _densityL;
            }
            else
            {
                // We're in the L_* state
                _pressureState = _pressureStar;
                _velocityState = _velocityStar;
                _densityState  = _densityShock("left");
            }

        }
        else
        {
            // The Left non-linear wave is a rarefaction
            double rareSpeedHead, rareSpeedTail, cRare;
            cRare = _cL * std::pow(_pressureStar/_pressureL , (_gamma - 1)/(2 * _gamma));
            rareSpeedTail = _velocityStar - cRare;
            rareSpeedHead = _velocityL - _cL;

            if (rareSpeedHead >= _posOverT)
            {
                // We're in the L state
                _pressureState = _pressureL;
                _velocityState = _velocityL;
                _densityState  = _densityL;
            }
            else if (rareSpeedTail < _posOverT)
            {
                // We're in the L_* state
                _pressureState = _pressureStar;
                _velocityState = _velocityStar;
                _densityState  = _densityRare("left");
            }
            else
            {
                // We're somewhere in the fan itself
                _velocityState = (2 / (_gamma + 1))
                                 * (
                                 _cL
                                 + _velocityL * ((_gamma - 1) / 2)
                                 + _posOverT);

                double coef = ( (2 / (_gamma + 1))
                              * (
                              _cL
                              + ((_gamma - 1) / 2)
                              * ( _velocityL - _posOverT) ) )
                              /
                              _cL;
                _pressureState = _pressureL * std::pow(coef, 2 * _gamma / (_gamma - 1));
                _densityState  = _densityL  * std::pow(coef, 2 / (_gamma - 1));
            }
        }

    }
    else
    {
        // We're in the R or R_* state
        if (_pressureStar > _pressureR)
        {
            // The Right non-linear wave is a shockwave
            double shockSpeed = _shockSpeed("right");

            // Decide between R and R_* state
            if (shockSpeed <= _posOverT)
            {
                // We're in the R state
                _pressureState = _pressureR;
                _velocityState = _velocityR;
                _densityState  = _densityR;
            }
            else
            {
                // We're in the R_* state
                _pressureState = _pressureStar;
                _velocityState = _velocityStar;
                _densityState  = _densityShock("right");
            }

        }
        else
        {
            // The Right non-linear wave is a rarefaction
            double rareSpeedHead, rareSpeedTail, cRare;
            cRare = _cR * std::pow(_pressureStar/_pressureR , (_gamma - 1)/(2 * _gamma));
            rareSpeedTail = _velocityStar + cRare;
            rareSpeedHead = _velocityR + _cR;

            if (rareSpeedHead <= _posOverT)
            {
                // We're in the R state
                _pressureState = _pressureR;
                _velocityState = _velocityR;
                _densityState  = _densityR;
            }
            else if (rareSpeedTail >= _posOverT)
            {
                // We're in the R_* state
                _pressureState = _pressureStar;
                _velocityState = _velocityStar;
                _densityState  = _densityRare("right");
            }
            else
            {
                // We're somewhere in the fan itself
                _velocityState = (2 / (_gamma + 1))
                                 * (
                                 -_cR
                                 + ((_gamma - 1) / 2) * _velocityR
                                 + _posOverT);

                double coef = ((2.0 / (_gamma + 1.0))
                              * (
                              _cR
                              - ((_gamma - 1.0) / 2.0)
                              * (_velocityR - _posOverT) ))
                              /
                              _cR;
                _pressureState = _pressureR * std::pow(coef, 2.0 * _gamma / (_gamma - 1.0));
                _densityState  = _densityR  * std::pow(coef, 2.0 / (_gamma - 1.0));;
            }
        }
    }

    // Compute and return the fluxes
    double energyState = computeEnergy(_pressureState, _densityState, _velocityState, _gamma);

    densityFlux = _densityState * _velocityState;
    momentumFlux = _densityState * std::pow(_velocityState, 2) + _pressureState;
    energyFlux = _velocityState * (energyState + _pressureState);
}
// =============================================================================

// =============================================================================
double ExactRiemannSolver::_computePressureStar()
{
    // Guess a value for the pressure in the star region
    double pStar = _guessPressureStar();

    // Perform the Newton-Raphson iterations
    size_t i = 0, maxIters = 20;
    double pTemp, fL = 0, fR = 0, dfL = 0, dfR = 0;

    while (true)
    {
        // Compute pressure functions
        _pressureFunctions(pStar, _pressureL, _densityL, _cL, fL, dfL);
        _pressureFunctions(pStar, _pressureR, _densityR, _cR, fR, dfR);

        // Compute new value of pStar
        pTemp = pStar - ( (fR + fL + _velocityR - _velocityL) / (dfL + dfR) );

        // Check for positivity and for convergence
        if (pTemp < 0.0)
        {
            pTemp = _tol;
        }
        else if ( 2.0 * (std::abs(pTemp - pStar) / (pTemp + pStar)) <= _tol)
        {
            // Change is below tolerance so we're done
            return pTemp;
        }

        // Check for iteration limit
        i++; // Increment counter
        if (i >= maxIters)
        {
            std::cout << "Max iterations reached in Newton-Raphson iterations"
                      << std::endl;
            return pTemp;
        }

        // Copy pTemp to pStar
        pStar = pTemp;
    }



}
// =============================================================================

// =============================================================================
double ExactRiemannSolver::_guessPressureStar()
{
    // Compute min and maximum pressures of the two sides of the interface
    double pMin = std::min(_pressureL, _pressureR);
    double pMax = std::max(_pressureL, _pressureR);

    // First compute the primitive variable approximation
    double pPrim = 0.5   * (_pressureL + _pressureR)
                   +
                   0.125 * (_velocityL - _velocityR)
                         * (_densityL + _densityR)
                         * (_cL + _cR);
    // Make sure it's not negative
    pPrim = std::max(_tol, pPrim);

    // Check to see if we should use the primitive variable approximation or not
    if ( ((pMax/pMin) <= 2.0) and (pMin <= pPrim) and (pPrim <= pMax) )
    {
        // Return pPrim and terminate this function
        return pPrim;
    }

    // If the previous statement was false then
    // Choose between 2-shock or 2-rarefaction approximations
    if (pPrim < pMin)
    {
        // 2-Rarefaction Approximation
        // NOTE: Toro disagrees with himself on this. In Toro the 2-rarefaction
        //       approximation to the pressure is different in the implementation
        //       (midway down page 157) and the equation (eqn. 4.46).  The two
        //       different equations give very different results.
        double pq, vm, ptL, ptR, p2Rare;
        pq = std::pow(_pressureL / _pressureR, (_gamma - 1.0)/(2.0 * _gamma));

        vm = ( (pq * _velocityL / _cL) + (_velocityR / _cR) + (2.0 / (_gamma - 1)) * (pq - 1.0) )
             /
             (pq/_cL + 1.0/_cR);

        ptL = 1.0 + ((_gamma - 1)/2.0) * (_velocityL - vm) / _cL;

        ptR = 1.0 + ((_gamma - 1)/2.0) * (vm - _velocityR) / _cR;

        p2Rare = 0.5 * (_pressureL * std::pow(ptL, 2.0 * _gamma / (_gamma - 1))
                      + _pressureR * std::pow(ptR, 2.0 * _gamma / (_gamma - 1)));

        return std::max(_tol, p2Rare);
    }
    else
    {
        // 2-Shock Approximation
        double gL = std::sqrt( (2 / ( (_gamma + 1) * _densityL)) /
                    (pPrim + _pressureL * ( (_gamma-1)/(_gamma+1) ) ) );
        double gR = std::sqrt( (2 / ( (_gamma + 1) * _densityR)) /
                    (pPrim + _pressureR * ( (_gamma-1)/(_gamma+1) ) ) );


        double p2Shock = // Equation on next lines for readability
        // Numerator
        (gL * _pressureL + gR * _pressureR - (_velocityR - _velocityL))
        /
        // Denominator
        (gL + gR);

        return std::max(_tol, p2Shock);
    }



}
// =============================================================================

// =============================================================================
void ExactRiemannSolver::_pressureFunctions(double const &pGuess,
                                       double const &pSide,
                                       double const &dSide,
                                       double const &cSide,
                                       double &f,
                                       double &df)
{
    if (pGuess > pSide)
    {
        // Shock
        double const aSide = (2.0 / (_gamma + 1.0)) / dSide;
        double const bSide = pSide * ((_gamma - 1.0) / (_gamma + 1.0));
        f = (pGuess - pSide) * std::sqrt(aSide / (pGuess + bSide));

        df = std::sqrt(aSide / (pGuess + bSide))
             * ( 1.0 - 0.5 * ( (pGuess - pSide) / (bSide + pGuess) ) );
    }
    else
    {
        // Rarefaction
        f = (2.0 * cSide / (_gamma - 1)) *
            (std::pow(pGuess/pSide, (_gamma - 1)/(2 * _gamma)) - 1.0);

        df = (1.0 / (dSide * cSide)) *
             std::pow( pGuess/pSide, (-1.0 - _gamma) / (2.0 * _gamma) );
    }

}
// =============================================================================

// =============================================================================
double ExactRiemannSolver::_shockSpeed(std::string const &side)
{
    // Figure out which variables to use
    double velocitySide, cSide, pressureSide;
    if (side == "left")
    {
        velocitySide = _velocityL;
        cSide        = -_cL;  // note the extra negative sign
        pressureSide = _pressureL;
    }
    else if (side == "right")
    {
        velocitySide = _velocityR;
        cSide        = _cR;
        pressureSide = _pressureR;
    }
    else
    {
        throw std::invalid_argument("Incorrect input for side into ExactRiemannSolver::_shockSpeed");
    }

    // Compute and return the shock speed
    return velocitySide + cSide * std::sqrt(
           (_pressureStar/pressureSide) * ((_gamma+1)/(2*_gamma))
           + ((_gamma-1)/(2*_gamma))
           );
}
// =============================================================================

// =============================================================================
double ExactRiemannSolver::_densityShock(std::string const &side)
{
    // Figure out which variables to use
    double densitySide, pressureSide;
    if (side == "left")
    {
        densitySide  = _densityL;
        pressureSide = _pressureL;
    }
    else if (side == "right")
    {
        densitySide  = _densityR;
        pressureSide = _pressureR;
    }
    else
    {
        throw std::invalid_argument("Incorrect input for side into ExactRiemannSolver::_densityShock");
    }

    // Compute and return the shock density
    return densitySide * (
           ( (_pressureStar/pressureSide) + ((_gamma-1)/(_gamma+1)) )
           /
           (((_gamma-1)/(_gamma+1)) * (_pressureStar/pressureSide) + 1) );
}
// =============================================================================

// =============================================================================
double ExactRiemannSolver::_densityRare(std::string const &side)
{
    // Figure out which variables to use
    double densitySide, pressureSide;
    if (side == "left")
    {
        densitySide  = _densityL;
        pressureSide = _pressureL;
    }
    else if (side == "right")
    {
        densitySide  = _densityR;
        pressureSide = _pressureR;
    }
    else
    {
        throw std::invalid_argument("Incorrect input for side into ExactRiemannSolver::_densityShock");
    }

    // Compute and return the rarefaction density
    return densitySide * std::pow(_pressureStar/pressureSide, 1/_gamma);
}
// =============================================================================