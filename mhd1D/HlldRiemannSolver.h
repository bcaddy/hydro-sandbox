/*!
 * \file HlldRiemannSolver.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains the HlldRiemannSolver class for solving the MHD Riemann
 * Problem.
 * \version 0.1
 * \date 2020-10-09
 *
 * \copyright Copyright (c) 2020
 *
 */

#pragma once

#include "MhdRiemannSolver.h"

/*!
 * \brief Solves the Riemann problem approximately using the HLLD Riemann solver
 *        detailed in Miyoshi & Kusano 2005
 */
class HlldRiemannSolver : public MhdRiemannSolver
{
private:
    /// The energy on the left side of the interface
    double _energyL;

    /// The energy on the right side of the interface
    double _energyR;

    /// The approximate wave speed of the left magnetosonic wave
    double _sL;
    /// The approximate wave speed of the right magnetosonic wave
    double _sR;

    /// The approximate wave speed of the left Alfven wave
    double _sStarL;
    /// The approximate wave speed of the right Alfven wave
    double _sStarR;

    /// The approximate wave speed of the middel contact wave
    double _sM;

    /*!
     * \brief Compute the S_M, S_L, and S_R wave speeds for the HLLD solver
     *
     */
    void _computeWaveSpeeds();

    /*!
     * \brief Compute the F_L or F_R flux directly from the interface states
     *
     * \param[in] density The density of the state
     * \param[in] velocity The velocity of the state
     * \param[in] pressure The pressure of the state
     * \param[in] energy The energy of the state
     * \param[out] densityFlux The density flux
     * \param[out] momentumFlux The momentum flux
     * \param[out] energyFlux The energy flux
     */
    void _computeStandardFluxes(double const &density,
                                double const &velocity,
                                double const &pressure,
                                double const &energy,
                                double &densityFlux,
                                double &momentumFlux,
                                double &energyFlux);

    /*!
     * \brief Compute the F_L^* or F_R^* flux by computing approximate values
     * for the star state and using them to compute the star state HLLD flux
     *
     * \param[in] density The density of the state
     * \param[in] velocity The velocity of the state
     * \param[in] pressure The pressure of the state
     * \param[in] energy The energy of the state
     * \param[in] sSide The wave speed estimate for that side
     * \param[out] densityFlux The density flux
     * \param[out] momentumFlux The momentum flux
     * \param[out] energyFlux The energy flux
     */
    void _computeStarFluxes(double const &density,
                            double const &velocity,
                            double const &pressure,
                            double const &energy,
                            double const &sSide,
                            double &densityFlux,
                            double &momentumFlux,
                            double &energyFlux);

public:
        /*!
     * \brief Solves the Riemann problem approximately by calling and organizing
     * member functions. Uses the same HLLD Riemann solver as in Miyoshi &
     * Kusano 2005
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
    void riemannMain(double const &densityL,
                     double const &velocityL,
                     double const &pressureL,
                     double const &densityR,
                     double const &velocityR,
                     double const &pressureR,
                     double &densityFlux,
                     double &momentumFlux,
                     double &energyFlux,
                     double const &posOverT = 0);

    /*!
     * \brief Construct a new Riemann Solver object
     *
     * \param[in] gamma The ratio of specific heats
     */
    HlldRiemannSolver(double const &gamma): MhdRiemannSolver(gamma) {};

    /*!
     * \brief Destroy the Hlld Riemann Solver object using default destructors
     *
     */
    ~HlldRiemannSolver() = default;
};
