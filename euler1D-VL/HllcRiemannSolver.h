/*!
 * \file HllcRiemannSolver.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains the HllcRiemannSolver class for solving the Riemann Problem.
 * \version 0.1
 * \date 2020-10-09
 *
 * \copyright Copyright (c) 2020
 *
 */

#pragma once

#include "RiemannSolver.h"

/*!
 * \brief Solves the Riemann problem approximately using the HLLC Riemann solver
 *        detailed in Batten et al. 1997 <em>"On the Choice of Wavespeeds for the
 *        HLLC Riemann Solver"</em>
 */
class HllcRiemannSolver : public RiemannSolver
{
private:
    /// The energy on the left side of the interface
    double _energyL;

    /// The energy on the right side of the interface
    double _energyR;

    /// The approximate wave speed of the left acoustic wave
    double _sL;

    /// The approximate wave speed of the middel contact wave
    double _sM;

    /// The approximate wave speed of the right acoustic wave
    double _sR;

    /*!
     * \brief Compute the S_M, S_L, and S_R wave speeds for the HLLC solver
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
     * for the star state and using them to compute the star state HLLC flux
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
     * member functions. Uses the same HLLC Riemann solver as in Batten et al.
     * 1997 <em>"On the Choice of Wavespeeds for the HLLC Riemann Solver"</em>
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
     */
    void riemannMain(double const &densityL,
                     double const &velocityL,
                     double const &pressureL,
                     double const &densityR,
                     double const &velocityR,
                     double const &pressureR,
                     double &energyFlux,
                     double &momentumFlux,
                     double &densityFlux,
                     double const &);

    /*!
     * \brief Construct a new Riemann Solver object
     *
     * \param[in] gamma The ratio of specific heats
     */
    HllcRiemannSolver(double const &gamma): RiemannSolver(gamma) {};

    /*!
     * \brief Destroy the Hllc Riemann Solver object using default destructors
     *
     */
    ~HllcRiemannSolver() = default;
};
