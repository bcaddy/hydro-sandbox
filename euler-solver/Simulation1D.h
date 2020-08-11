/*!
 * \file Simulation1D.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains the Simulation1D class for solving the 1D Euler equations.
 * \date 2020-06-25
 *
 * \copyright Copyright (c) 2020
 *
 */
#pragma once

#include <string>
#include <vector>
#include "Grid1D.h"
#include "ExactRiemannSolver.h"

/*!
 * \brief Solves the Euler equations
 * \details This class handles all the heavy lifting of running a 1D hydro
 * simulation. It uses the Grid1D class for the grid and then computes time steps,
 * interface steps, solves the Riemann problem using the ExactRiemannSolver class,
 * does the conservative update, keeps track of the time, etc.
 *
 */
class Simulation1D
{
private:
    /// The physical length of the simulation in meters
    double const _physLen;
    /// Courant–Friedrichs–Lewy (CFL) Number
    double const _cflNum;
    /// The physical size of each cell in meters
    double const _deltaX;
    /// The ratio of specific heats $\gamma$
    double const _gamma;
    /// The time step for each interation
    double _timeStep;

    /// The vector to store the density at the left side of the interface.
    /// _densityInterfaceL is the density state on the left side of the i-1/2
    /// interface
    std::vector<double> _densityInterfaceL;
    /// The vector to store the velocity at the left side of the interface.
    /// _velocityInterfaceL is the velocity state on the left side of the i-1/2
    /// interface
    std::vector<double> _velocityInterfaceL;
    /// The vector to store the pressure at the left side of the interface.
    /// _pressureInterfaceL is the pressure state on the left side of the i-1/2
    /// interface
    std::vector<double> _pressureInterfaceL;

    /// The vector to store the density at the right side of the interface.
    /// _densityInterfaceR is the density state on the right side of the i-1/2
    /// interface
    std::vector<double> _densityInterfaceR;
    /// The vector to store the velocity at the right side of the interface.
    /// _velocityInterfaceR is the velocity state on the right side of the i-1/2
    /// interface
    std::vector<double> _velocityInterfaceR;
    /// The vector to store the pressure at the right side of the interface.
    /// _pressureInterfaceR is the pressure state on the right side of the i-1/2
    /// interface
    std::vector<double> _pressureInterfaceR;


    /// The vector to store the density fluxes. _densityFlux[i] is the density
    /// flux through the i-1/2 interface
    std::vector<double> _densityFlux;
    /// The vector to store the momentum fluxes. _momentumFlux[i] is the momentum
    /// flux through the i-1/2 interface
    std::vector<double> _momentumFlux;
    /// The vector to store the energy fluxes. _energyFlux[i] is the energy
    /// flux through the i-1/2 interface
    std::vector<double> _energyFlux;

    /// The object used to solve the Riemann Problem. See ExactRiemannSolver for the
    /// full documentation.
    ExactRiemannSolver _riemannSolver;

    /*!
     * \brief Set the initial conditions. Currently only supports a Sod shock
     * tube.
     *
     * \param[in] initialConditionsKind The type of initial conditions to use.
     *            Options are "sod", "indexCheck", "advectionStep",  and
     *            "advectionGauss".
     *            | Keyword        | Initial Conditions                                                       |
     *            |----------------|--------------------------------------------------------------------------|
     *            | sod            | A sod shock tube                                                         |
     *            | indexCheck     | Set each primitive variable to the value of the grid index at that point |
     *            | advectionStep  | Step function advection                                                  |
     *            | advectionGauss | Gaussion function advection                                              |
     */
    void _setInitialConditions(std::string const &initialConditionsKind);

    /*!
     * \brief Slope of a given primitive variable
     *
     * \details Computes the slope. Uses Simulation1D::limiterKind to choose
     *          which limiter to use.
     *          | limiterKind | Which Limiter/Slope is used?                |
     *          |-------------|---------------------------------------------|
     *          | zeroSlope   | Always return a slope of zero               |
     *          | centerDiff  | A centered difference, no limiting          |
     *          | minMod      | Minmod limiter                              |
     *          | MC          | Monotonized Central Difference (MC) Limiter |
     *
     * \param[in] primitive The primitive variable to compute the slope of.
     * \param[in] idx The cell in which to compute the slope.
     * \return double The limited slope.
     */
    double _slope(std::vector<double> const &primitive,
                  size_t const &idx);

    /*!
     * \brief Computes the interface states using the Piecewise Linear Method
     *        (PLM) from "Introduction to Computational Astrophysical
     *        Hydrodynamics" by Zingale, git version: 4de1fef51af5 Section 8.2.2
     */
    void _piecewiseLinearReconstruction();

    /*!
     * \brief Computes the interface states using the Piecewise Constant Method (PCM)
     *
     */
    void _piecewiseConstantReconstruction();

public:
    /// The primary grid
    Grid1D grid;

    /// The current time in the simulation
    double currentTime;

    /// The reconstruction scheme to use
    std::string const reconstructionKind;

    /// The kind of slope limiter to use
    std::string const limiterKind;

    /*!
     * \brief Compute the time step using the CFL condition
     */
    void computeTimeStep();

    /*!
     * \brief Get the Time Step object
     *
     * \return double The value of the time step
     */
    double getTimeStep() { return _timeStep; };

    /*!
     * \brief Increases the current time by the value of Simulation1D::_timeStep
     *
     */
    void updateCurrentTime() { currentTime += _timeStep; };

    /*!
     * \brief Compute all the interface states.
     *
     * \details Calls either Simulation1D::_piecewiseConstant or
     *          Simulation1D::_piecewiseLinear depending on if
     *          Simulation1D::reconstructionKind is "PCM" or "PCL" respectively
     */
    void interfaceStates();

    /*!
     * \brief Solves the Riemann problem exactly by calling the main function of
     * the ExactRiemannSolver class
     */
    void solveRiemann();


    /*!
     * \brief Performs the conservative update
     *
     * \param[in] idxInput Which cell we're computing the conservative update for
     * \param[in] densityFluxL The density flux on the left side
     * \param[in] momentumFluxL The momentum flux on the left side
     * \param[in] energyFluxL The energy flux on the left side
     * \param[in] densityFluxR The density flux on the left side
     * \param[in] momentumFluxR The momentum flux on the left side
     * \param[in] energyFluxR The energy flux on the left side
     */
    void conservativeUpdate();

    /*!
     * \brief Construct a new Simulation1D object.
     *
     * \details Sets the private const variables (_limiterKind, _physLen, _cflNum,
     * _deltaX). Initializes the grid and the temporary grid and then sets the
     * initial conditions of the grid.
     *
     * \param[in] physicalLength Physical length of the simulation in meters.
     * \param[in] gamma The ratio of specific heats
     * \param[in] CFL The CFL Number
     * \param[in] reals The number of real grid cells
     * \param[in] ghosts The number of ghost cells
     * \param[in] initialConditionsKind Which initial conditions to use
     * \param[in] reconstructionKind Which kind of interface reconstruction to
     *            use. Option are "PCM" for Piecewise Constant Method and "PLM"
     *            for Piecewise Linear Method
     * \param[in] limiterKind What kind of limiter to use. Options are
     *            "centerDiff", "minMod", and "MC"
     * \param[in] boundaryConditions Which kind of boundary conditions to use
     * \param[in] saveDir The directory to save the grid to
     */
    Simulation1D(double const &physicalLength,
                 double const &gamma,
                 double const &CFL,
                 size_t const &reals,
                 size_t const &ghosts,
                 std::string const &initialConditionsKind,
                 std::string const &reconstructionKind,
                 std::string const &limiterKind,
                 std::string const &boundaryConditions,
                 std::string const &saveDir);
    /*!
     * \brief Destroy the Simulation1D object. Uses the default constructor
     */
    ~Simulation1D()=default;
};