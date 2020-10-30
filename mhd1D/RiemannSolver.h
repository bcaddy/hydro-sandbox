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
 * \brief A (mostly) pure virtual base class for all Riemann solvers. Contains
 * all the member variables that are common to all Riemann solver along with a
 * pure virtual declaration of the riemannMain method
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

    /// The density of the found state
    double _densityState;
    /// The velocity of the found state
    double _velocityState;
    /// The pressure of the found state
    double _pressureState;

public:
    /*!
     * \brief Solves the Riemann problem by calling and organizing member
     * functions.
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
     * \param[in] posOverT OPTIONAL: The value of the position divided by the
     * current time. Alway equal to zero for numerical solutions and as such
     * defaults to it
     */
    virtual void riemannMain(double const &densityL,
                             double const &velocityL,
                             double const &pressureL,
                             double const &densityR,
                             double const &velocityR,
                             double const &pressureR,
                             double &densityFlux,
                             double &momentumFlux,
                             double &energyFlux,
                             double const &posOverT = 0.0) = 0;

    /*!
     * \brief Get the Density State object
     *
     * \return double The density of the current state
     */
    double getDensityState(){return _densityState;};

    /*!
     * \brief Get the Velocity State object
     *
     * \return double The velocity of the current state
     */
    double getVelocityState(){return _velocityState;};

    /*!
     * \brief Get the Pressure State object
     *
     * \return double The pressure of the current state
     */
    double getPressureState(){return _pressureState;};



    RiemannSolver(double const &gamma) : _gamma(gamma){};
    virtual ~RiemannSolver() = default;
};