/*!
 * \file RiemannSolver.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains the RiemannSolver class for solving the Riemann Problem.
 * \version 0.1
 * \date 2020-07-14
 *
 * \copyright Copyright (c) 2020
 *
 */

#pragma once

#include <string>

/*!
 * \brief Solves the Riemann problem exactly using the same exact Riemann solver
 *        as in Toro "Riemann Solver and Numerical Methods for Fluid Dynamics 3ed"
 */
class RiemannSolver
{
private:
    // =========================================================================
    // Declare all the private variables
    // =========================================================================

    /// The ratio of specific heats gamma
    double const _gamma;

    /// The tolerance for the Newton-Raphson iterations used to compute
    /// _pressureStar
    double const _tol = 1.0E-6;

    /// The density on the right side of the interface
    double _densityR;
    /// The velocity on the right side of the interface
    double _velocityR;
    /// The pressure on the right side of the interface
    double _pressureR;

    /// The density on the left side of the interface
    double _densityL;
    /// The velocity on the left side of the interface
    double _velocityL;
    /// The pressure on the left side of the interface
    double _pressureL;

    /// The density in the star region
    double _densityStar;
    /// The velocity in the star region
    double _velocityStar;
    /// The pressure in the star region
    double _pressureStar;

    /// The density of the found state
    double _densityState;
    /// The velocity of the found state
    double _velocityState;
    /// The pressure of the found state
    double _pressureState;

    /// The sound speeds on the right side of the interface
    double _cR;

    /// The sound speeds on the left side of the interface
    double _cL;

    /// The position divided by the time.
    double _posOverT;

    /// The energy flux
    double _energyFlux;
    /// The momentum flux
    double _momentumFlux;
    /// The mass flux
    double _massFlux;

    // =========================================================================
    // End declaring all the private variables
    // Start declargin all the private methods
    // =========================================================================
    /*!
     * \brief Compute the pressure in the star region using the
     * RiemannSolver::_guessPressureStar function and Newton-Raphson iterations
     *
     * \return double The pressure in the star region
     */
    double _computePressureStar();

    /*!
     * \brief Provide an estimate of the pressure in the star region
     *
     * \return double An estimate of the pressure in the star region
     */
    double _guessPressureStar();

    /*!
     * \brief Compute the pressure functions. See Toro section 4.3.2 for more
     *
     * \param[in] pGuess The current guess for the pressure in the star region
     * \param[in] pSide The pressure on whichever side we're computing
     * \param[in] dSide The density on whichver side we're computing
     * \param[in] cSide The sound speed on whichver side we're computing
     * \param[out] f The output of the pressure function
     * \param[out] df The output of the derivative of the pressure function
     */
    void _pressureFunctions(double const &pGuess,
                            double const &pSide,
                            double const &dSide,
                            double const &cSide,
                            double &f,
                            double &df);

    /*!
     * \brief Compute the speed of a shockwave
     *
     * \param side Which side to compute the shock on
     * \return double The speed of the shock wave
     */
    double _shockSpeed(std::string const &side);

    /*!
     * \brief Compute the density of the star region next to a shock
     *
     * \param side Which side to compute the density on
     * \return double The density next to the shock
     */
    double _densityShock(std::string const &side);
    // =========================================================================
    // End declaring all the private methods
    // Start declargin all the public members
    // =========================================================================
public:
    /*!
     * \brief Solves the Riemann problem exactly by calling and organizing member
     * functions. Uses the same exact Riemann solver as in Toro "Riemann Solvers
     * and Numerical Methods for Fluid Dynamics 3ed"
     *
     * \param[in] densityR  The density on the right side of the interface
     * \param[in] velocityR The velocity on the right side of the interface
     * \param[in] pressureR The pressure on the right side of the interface
     * \param[in] densityL The density on the left side of the interface
     * \param[in] velocityL The velocity on the left side of the interface
     * \param[in] pressureL The pressure on the left side of the interface
     * \param[in] posOverT The value of the position divided by the current time.
     * Alway equal to zero for numerical solutions
     * \param[out] energyFlux The energy flux that is being solved for
     * \param[out] momentumFlux The momentum flux that is being solved for
     * \param[out] massFlux The mass flux that is being solved for
     */
    void riemannMain(double const &densityR,
                     double const &velocityR,
                     double const &pressureR,
                     double const &densityL,
                     double const &velocityL,
                     double const &pressureL,
                     double const &posOverT,
                     double const &energyFlux,
                     double const &momentumFlux,
                     double const &massFlux);

    /*!
     * \brief Construct a new Riemann Solver object
     *
     * \param[in] gamma The ratio of specific heats
     */
    RiemannSolver(double const &gamma);
    /*!
     * \brief Destroy the Riemann Solver object. Uses default destructors
     *
     */
    ~RiemannSolver() = default;
};