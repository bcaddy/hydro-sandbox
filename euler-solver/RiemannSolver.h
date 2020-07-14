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

    /// The sound speeds on the right side of the interface
    double _cR;

    /// The sound speeds on the left side of the interface
    double _cL;

    /// The position
    double _position;

    /// The time
    double _time;

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
     * \param[in] position The position, always equal to zero for numerical
     *            algorithms
     * \param[in] time The time, only used in the context of position/time so 
     *            usually irrelevant. Just make sure to set it to some non-zero 
     *            value
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
                     double const &position,
                     double const &time,
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