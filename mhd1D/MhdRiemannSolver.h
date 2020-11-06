/*!
 * \file MhdRiemannSolver.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains the base virtual class for all MHD Riemann Solvers
 * \version 0.1
 * \date 2020-10-15
 *
 * \copyright Copyright (c) 2020
 *
 */

#pragma once

#include <vector>

/*!
 * \brief A (mostly) pure virtual base class for all MHD Riemann solvers. Contains
 * all the member variables that are common to all Riemann solver along with a
 * pure virtual declaration of the riemannMain method
 *
 */
class MhdRiemannSolver
{
protected:
    /// The ratio of specific heats gamma
    double const _gamma;

    /// The density on the left side of the interface
    double _densityL;
    /// The velocity on the left side of the interface
    std::vector<double> _velocityL;
    /// The pressure on the left side of the interface
    double _pressureL;
    /// The magnetic field
    std::vector<double> _magneticL;

    /// The density on the right side of the interface
    double _densityR;
    /// The velocity on the right side of the interface
    std::vector<double> _velocityR;
    /// The pressure on the right side of the interface
    double _pressureR;
    /// The magnetic field
    std::vector<double> _magneticR;

    /// The density of the found state
    double _densityState;
    /// The velocity of the found state
    std::vector<double> _velocityState;
    /// The pressure of the found state
    double _pressureState;
    /// The magnetic field
    std::vector<double> _magneticState;

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
                             std::vector<double> const &velocityL,
                             double const &pressureL,
                             std::vector<double> const &magneticL,
                             double const &densityR,
                             std::vector<double> const &velocityR,
                             double const &pressureR,
                             std::vector<double> magneticR,
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
     * \return std::vector<double> The velocity of the current state
     */
    std::vector<double> getVelocityState(){return _velocityState;};

    /*!
     * \brief Get the Pressure State object
     *
     * \return double The pressure of the current state
     */
    double getPressureState(){return _pressureState;};

    /*!
     * \brief Get the magnetic field State object
     *
     * \return std::vector<double> The magnetic field of the current state
     */
    std::vector<double> getMagneticState(){return _magneticState;};




    /*!
     * \brief Construct a new Riemann Solver object.
     *
     * \param gamma The ratio of the specific heats
     */
    MhdRiemannSolver(double const &gamma)
        : _gamma(gamma),
          _magneticL(3),
          _magneticR(3),
          _magneticState(3)
        {};
    virtual ~MhdRiemannSolver() = default;
};