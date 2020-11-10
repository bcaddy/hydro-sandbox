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

    /// The density in the left star state
    double _densityStarL;
    /// The density in the right star state
    double _densityStarR;

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

    /// Indicates if we'll run into a divide by zero error. True(1) for yes,
    /// False(0) for no
    bool _divZero;

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
     * \param[in] pressureTot The total pressure of the state
     * \param[in] magnetic the magnetic field of the state
     * \param[in] energy The energy of the state
     * \param[out] densityFlux The density flux
     * \param[out] momentumFlux The momentum flux
     * \param[out] magneticFlux The magnetic field flux
     * \param[out] energyFlux The energy flux
     */
    void _computeStandardFluxes(double const &density,
                                std::vector<double> const &velocity,
                                double const &pressureTot,
                                std::vector<double> const &magnetic,
                                double const &energy,
                                double &densityFlux,
                                std::vector<double> &momentumFlux,
                                std::vector<double> &magneticFlux,
                                double &energyFlux);

    /*!
     * \brief Compute the velocity and magnetic field in the star region. Also,
     * determine the _divZero member variable
     *
     * \param[in] density The density of that side of the interface
     * \param[in] velocity The velocity of that side of the interface
     * \param[in] pressure The pressure of that side of the interface
     * \param[in] magnetic The magnetic field of that side of the interface
     * \param[in] sSide The speed of the magnetosonic wave for the side being computed
     * \param[out] velocityStar The velocity in the star state
     * \param[out] magneticStar The magnetic field in the star state
     */
    void _computeVandBStar(double const &density,
                           std::vector<double> const &velocity,
                           double const &pressure,
                           std::vector<double> const &magnetic,
                           double const &sSide,
                           std::vector<double> &velocityStar,
                           std::vector<double> &magneticStar);

    /*!
     * \brief Compute the F_L^* or F_R^* flux by computing approximate values
     * for the star state and using them to compute the star state HLLD flux
     *
     * \param[in] density The density of the state
     * \param[in] velocity The velocity of the state
     * \param[in] pressure The pressure of the state
     * \param[in] pressureTot The total pressure of the state
     * \param[in] magnetic the magnetic field of the state
     * \param[in] energy The energy of the state
     * \param[in] sSide The magnetosonic wave speed estimate for that side
     * \param[in] densityStarSide The density on that side of the star region
     * \param[out] densityFlux The density flux
     * \param[out] momentumFlux The momentum flux
     * \param[out] magneticFlux The magnetic field flux
     * \param[out] energyFlux The energy flux
     */
    void _computeStarFluxes(double const &density,
                            std::vector<double> const &velocity,
                            double const &pressure,
                            double const &pressureTot,
                            std::vector<double> const &magnetic,
                            double const &energy,
                            double const &sSide,
                            double const &densityStarSide,
                            double &densityFlux,
                            std::vector<double> &momentumFlux,
                            std::vector<double> &magneticFlux,
                            double &energyFlux);

    /*!
     * \brief Compute the F_L^** or F_R^** flux by computing approximate values
     * for the double star state and using them to compute the double star state
     * HLLD flux
     *
     * \param[in] magnetic the magnetic field of the state
     * \param[in] sStarSide The Alfven wave speed estimate for that side
     * \param[in] densityStarSide The density on that side of the star region
     * \param[in] sideSign Which side to compute. -1.0 is left and 1.0 is right.
     * This is used in the computation of the double star state energy and so it
     * MUST be one of those two values
     * \param[out] densityFlux The density flux
     * \param[out] momentumFlux The momentum flux
     * \param[out] magneticFlux The magnetic field flux
     * \param[out] energyFlux The energy flux
     */
    void _computeDblStarFluxes(std::vector<double> const &magnetic,
                               double const &sStarSide,
                               double const &densityStarSide,
                               double const &sideSign,
                               double &densityFlux,
                               std::vector<double> &momentumFlux,
                               std::vector<double> &magneticFlux,
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
     * \param[in] magneticL The magnetic field on the left side of the interface
     * \param[in] densityR  The density on the right side of the interface
     * \param[in] velocityR The velocity on the right side of the interface
     * \param[in] pressureR The pressure on the right side of the interface
     * \param[in] magneticR The magnetic field on the right side of the interface
     * \param[out] densityFlux The density flux that is being solved for
     * \param[out] momentumFlux The momentum flux that is being solved for
     * \param[out] magneticFlux The magnetic field flux that is being solved for
     * \param[out] energyFlux The energy flux that is being solved for
     * \param[in] posOverT OPTIONAL: The value of the position divided by the
     * current time. Alway equal to zero for numerical solutions and as such
     * defaults to it
     */
    void riemannMain(double const &densityL,
                     std::vector<double> const &velocityL,
                     double const &pressureL,
                     std::vector<double> const &magneticL,
                     double const &densityR,
                     std::vector<double> const &velocityR,
                     double const &pressureR,
                     std::vector<double> const &magneticR,
                     double &densityFlux,
                     std::vector<double> &momentumFlux,
                     std::vector<double> &magneticFlux,
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
