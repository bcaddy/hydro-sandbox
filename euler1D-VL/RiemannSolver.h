/*!
 * \file RiemannSolver.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains the base virtual class for all Riemann Solvers
 * \version 0.1
 * \date 2020-10-15
 *
 * \copyright Copyright (c) 2020
 *
 */

#pragma once

/*!
 * \brief A (mostly) virtual base class for all Riemann solvers. Contains all
 * the member variables that are common to all Riemann solver along with a pure
 * virtual declaration of the riemannMain method
 *
 */
class RiemannSolver
{
protected:
    /// The ratio of specific heats gamma
    double const _gamma;

    /// The density on the left side of the interface
    double _densityL;
    /// The velocity on the left side of the interface
    double _velocityL;
    /// The pressure on the left side of the interface
    double _pressureL;

    /// The density on the right side of the interface
    double _densityR;
    /// The velocity on the right side of the interface
    double _velocityR;
    /// The pressure on the right side of the interface
    double _pressureR;

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
     * \param[out] energyFlux The energy flux that is being solved for
     * \param[out] momentumFlux The momentum flux that is being solved for
     * \param[out] densityFlux The density flux that is being solved for
     * \param[in] posOverT OPTIONAL: and only present on the Exact Riemann
     * Solver The value of the position divided by the current time. Alway equal
     * to zero for numerical solutions and as such defaults to it
     */
    virtual void riemannMain(double const &densityL,
                             double const &velocityL,
                             double const &pressureL,
                             double const &densityR,
                             double const &velocityR,
                             double const &pressureR,
                             double &energyFlux,
                             double &momentumFlux,
                             double &densityFlux,
                             double const &posOverT = 0.0) = 0;

    RiemannSolver(double const &gamma) : _gamma(gamma){};
    virtual ~RiemannSolver() = default;
};