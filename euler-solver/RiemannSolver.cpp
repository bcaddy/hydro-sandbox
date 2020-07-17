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
#include <iostream>

#include "RiemannSolver.h"

// =============================================================================
void RiemannSolver::riemannMain(double const &densityR,
                                double const &velocityR,
                                double const &pressureR,
                                double const &densityL,
                                double const &velocityL,
                                double const &pressureL,
                                double const &posOverT,
                                double &energyFlux,
                                double &momentumFlux,
                                double &massFlux)
{
    // Copy arguments to member variables
    _densityR     = densityR;
    _velocityR    = velocityR;
    _pressureR    = pressureR;
    _densityL     = densityL;
    _velocityL    = velocityL;
    _pressureL    = pressureL;
    _posOverT     = posOverT;

    // Compute the sound speeds
    _cR = std::sqrt(_gamma * pressureR / densityR);
    _cL = std::sqrt(_gamma * pressureL / densityL);

    // Check for a vacuum
    if ((2 / (_gamma - 1)) * (_cL + _cR) < (_velocityR - _velocityL))
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
            if (shockSpeed > 0.0)
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

            if ((rareSpeedHead > 0.0) && (rareSpeedTail > 0.0))
            {
                // We're in the L state
                _pressureState = _pressureL;
                _velocityState = _velocityL;
                _densityState  = _densityL;
            }
            else if ((rareSpeedHead < 0.0) && (rareSpeedTail < 0.0))
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

                double coef = (2 / (_gamma + 1))
                              + (_gamma - 1) / ((_gamma + 1) * _cL)
                              * (_velocityL - _posOverT);
                _pressureState = _pressureL * std::pow(coef, 2 * _gamma / (_gamma - 1));
                _densityState  = _densityL * std::pow(coef, 2 / (_gamma - 1));;
            }
        }

    }
    else
    {
        // We're in the R or R_* state
        if (_pressureStar > _pressureL)
        {
            // The Right non-linear wave is a shockwave
            double shockSpeed = _shockSpeed("right");

            // Decide between R and R_* state
            if (shockSpeed <= 0.0)
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
            rareSpeedHead = _velocityL + _cR;

            if ((rareSpeedHead < 0.0) && (rareSpeedTail < 0.0))
            {
                // We're in the R state
                _pressureState = _pressureR;
                _velocityState = _velocityR;
                _densityState  = _densityR;
            }
            else if ((rareSpeedHead > 0.0) && (rareSpeedTail > 0.0))
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
                                 + _velocityR * ((_gamma - 1) / 2)
                                 + _posOverT);

                double coef = (2 / (_gamma + 1))
                              - (_gamma - 1) / ((_gamma + 1) * _cR)
                              * (_velocityR - _posOverT);
                _pressureState = _pressureL * std::pow(coef, 2 * _gamma / (_gamma - 1));
                _densityState  = _densityL * std::pow(coef, 2 / (_gamma - 1));;
            }
        }
    }

    // Compute and return the fluxes
    massFlux = _densityState * _velocityState;
    momentumFlux = _densityState * std::pow(_velocityState, 2) + _pressureState;
    energyFlux = _velocityState * _pressureState / (_gamma - 1)
                 + 0.5 * _densityState * std::pow(_velocityState, 3)
                 + _velocityState * _pressureState;
}
// =============================================================================

// =============================================================================
double RiemannSolver::_computePressureStar()
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
        pTemp = pStar - ( (fR + fL - (_velocityR - _velocityL)) / (dfL + dfR) );

        // Check for positivity and for convergence
        if (pTemp < 0.0)
        {
            pStar = _tol;
        }
        else if ( (std::abs(pTemp - pStar) / (0.5 * (pTemp + pStar))) < _tol)
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
    if ( ((pMax/pMin) <= 2.0) && (pMin < pPrim) && (pPrim < pMax) )
    {
        // Return pPrim and terminate this function
        return pPrim;
    }

    // If the previous statement was false then
    // Choose between 2-shock or 2-rarefaction approximations
    if (pPrim < pMin)
    {
        // 2-Rarefaction Approximation
        double p2Rare = // Equation on next lines for readability
        std::pow(
                // Numerator
                ( _cL + _cR - 0.5 * (_gamma - 1) * (_velocityR - _velocityL) )
                /
                // Denominator
                ( std::pow(_cL/_pressureL, (_gamma - 1)/(2*_gamma) )
                + std::pow(_cR/_pressureR, (_gamma - 1)/(2*_gamma) ) )
                // Exponent
                , 2*_gamma/(_gamma-1));

        return std::max(_tol, p2Rare);
    }
    else
    {
        // 2-Shock Approximation
        double gL = std::sqrt( 2 / ( (_gamma + 1) * _densityL *
                    (pPrim + _pressureL * ( (_gamma-1)/(_gamma+1) )) ) );
        double gR = std::sqrt( 2 / ( (_gamma + 1) * _densityR *
                    (pPrim + _pressureR * ( (_gamma-1)/(_gamma+1) )) ) );

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
void RiemannSolver::_pressureFunctions(double const &pGuess,
                                       double const &pSide,
                                       double const &dSide,
                                       double const &cSide,
                                       double &f,
                                       double &df)
{
    if (pGuess > pSide)
    {
        // Shock
        f = (pGuess - pSide) * std::sqrt(2 /
            (dSide * ( pGuess * (_gamma + 1) + pSide * (_gamma - 1))));

        df = (1 -
        ( (pGuess - pSide) / (2 * (pSide * ((_gamma-1)/(_gamma+1)) + pGuess ))))
        * std::sqrt(2 / (dSide * (pGuess * (_gamma + 1) + pSide * (_gamma - 1))));
    }
    else
    {
        // Rarefaction
        f = (2 * cSide / (_gamma - 1)) *
            (std::pow(pGuess/pSide , (_gamma - 1)/(2 * _gamma)) - 1);

        df = (1 / (pSide * cSide)) *
             std::pow( pGuess/pSide, (1 - _gamma) / (2 * _gamma) );
    }

}
// =============================================================================

// =============================================================================
double RiemannSolver::_shockSpeed(std::string const &side)
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
        throw std::invalid_argument("Incorrect input for side into RiemannSolver::_shockSpeed");
    }

    // Compute and return the shock speed
    return velocitySide + cSide * std::sqrt(
           (_pressureStar/pressureSide) * ((_gamma+1)/(2*_gamma))
           + ((_gamma-1)/(2*_gamma))
           );
}
// =============================================================================

// =============================================================================
double RiemannSolver::_densityShock(std::string const &side)
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
        throw std::invalid_argument("Incorrect input for side into RiemannSolver::_densityShock");
    }

    // Compute and return the shock density
    return densitySide * (
           ( (_pressureStar/pressureSide) + ((_gamma-1)/(_gamma+1)) / (((_gamma-1)
           /
           (_gamma+1)) * (_pressureStar/pressureSide)-1) )
           );
}
// =============================================================================

// =============================================================================
double RiemannSolver::_densityRare(std::string const &side)
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
        throw std::invalid_argument("Incorrect input for side into RiemannSolver::_densityShock");
    }

    // Compute and return the rarefaction density
    return densitySide * std::pow(_pressureStar/pressureSide, 1/_gamma);
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