/*!
 * \file ExactRiemannSolver.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains the ExactRiemannSolver class for solving the Riemann Problem.
 * \version 0.1
 * \date 2020-07-14
 *
 * \copyright Copyright (c) 2020
 *
 */

#pragma once

#include <string>

#include "RiemannSolver.h"

/*!
 * \brief Solves the Riemann problem exactly using the same exact Riemann solver
 *        as in Toro "Riemann Solver and Numerical Methods for Fluid Dynamics 3ed"
 */
class ExactRiemannSolver : public RiemannSolver
{
private:
    // =========================================================================
    // Declare all the private variables
    // =========================================================================
    /// The tolerance for the Newton-Raphson iterations used to compute
    /// _pressureStar
    double const _tol = 1.0E-6;

    /// The velocity in the star region
    double _velocityStar;
    /// The pressure in the star region
    double _pressureStar;

    /// The sound speeds on the left side of the interface
    double _cL;
    /// The sound speeds on the right side of the interface
    double _cR;

    /// The position divided by the time.
    double _posOverT;

    // =========================================================================
    // End declaring all the private variables
    // Start declargin all the private methods
    // =========================================================================
    /*!
     * \brief Compute the pressure in the star region using the
     * ExactRiemannSolver::_guessPressureStar function and Newton-Raphson iterations
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

    /*!
     * \brief Compute the density of the star region next to a rarefaction
     *
     * \param side Which side to compute the density on
     * \return double The density next to the rarefaction
     */
    double _densityRare(std::string const &side);
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
     * \param[in] densityL The density on the left side of the interface
     * \param[in] velocityL The velocity on the left side of the interface
     * \param[in] pressureL The pressure on the left side of the interface
     * \param[in] densityR  The density on the right side of the interface
     * \param[in] velocityR The velocity on the right side of the interface
     * \param[in] pressureR The pressure on the right side of the interface
     * \param[out] densityFlux The density flux that is being solved for
     * \param[out] momentumFlux The momentum flux that is being solved for
     * \param[out] energyFlux The energy flux that is being solved for
     * \param[in] posOverT OPTIONAL: The value of the position divided by the current time.
     * Alway equal to zero for numerical solutions and as such defaults to it
     */
    void riemannMain(double const &densityL,
                     double const &velocityL,
                     double const &pressureL,
                     double const &densityR,
                     double const &velocityR,
                     double const &pressureR,
                     double &densityFlux,
                     double &momentumFlux,
                     double &energyFlux,
                     double const &posOverT = 0.0);

    /*!
     * \brief Construct a new Riemann Solver object
     *
     * \param[in] gamma The ratio of specific heats
     */
    ExactRiemannSolver(double const &gamma) : RiemannSolver(gamma) {};
    /*!
     * \brief Destroy the Riemann Solver object. Uses default destructors
     *
     */
    ~ExactRiemannSolver() = default;
};