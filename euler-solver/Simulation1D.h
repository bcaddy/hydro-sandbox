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
#include "Grid1D.h"

/*!
 * \brief Solves the Euler equations
 * \details This class handles all the heavy lifting of running a 1D hydro 
 * simulation. It uses the Grid1D class for the grid and then computes time steps,
 * interface steps, solves the Riemann problem, etc. 
 * 
 */
class Simulation1D
{
private:
    /// The type of limiter used
    std::string const _limiterKind;
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

    /// The temporary grid used for storing the next time step while it's being
    /// computed
    Grid1D _tempGrid;

    /*!
     * \brief Set the initial conditions. Currently only supports a Sod shock 
     * tube.
     * 
     * \param[in] initialConditionsKind The type of initial conditions to use.
     * currently unused.
     * 
     * \todo Add more initial conditions. Make it a template?
     */
    void _setInitialConditions(std::string const &initialConditionsKind);

    /*!
     * \brief Compute the MC limited slope of a given primitive variable
     * 
     * \todo Add other kinds of limiters. Templated?
     * 
     * \param[in] primitive The primitive variable to compute the slope of.
     * \param[in] i The cell in which to compute the slope.
     * \return double The limited slope.
     */
    double _slope(std::vector<double> const &primitive,
                 size_t const &i);

    /*!
     * \brief Compute the eigenvalues and vectors of the Euler equations for a
     * given cell.
     * 
     * \details The eigenvalues and vectors compute are given in *Introduction
     * to Computational Astrophysical Hydrodynamics* by Zingale pages 96-97, 
     * exercise 7.2 and 7.3. The eigenvector matrices are defined by equations
     * 7.28 and 7.29.
     * 
     * \param[in] idx The cell in which to compute the eigenvalues and
     *                eigenvectors
     * \param[out] eigVal A std::vector containing the eigenvalues in the order
     *                    \f$ \lambda^-, \lambda^0, \lambda^+ \f$.
     * \param[out] rEigVec A 2D std::vector which contains the right eigenvectors.
     *                     They are all column vectors and are arranged side by
     *                     side in the same -, 0, + ordering as the eigenvalues. 
     *                     Indices are [row][column].
     * \param[out] lEigVec A 2D std::vector which contains the left eigenvectors.
     *                     They are all row vectors and are arranged stacked in
     *                     the same -, 0, + ordering as the eigenvalues.Indices 
     *                     are [row][column].
     */
    void _computeEigens(size_t const &idx,
                        std::vector<double> &eigVal,
                        std::vector<std::vector<double>> &rEigVec,
                        std::vector<std::vector<double>> &lEigVec);

public:
    /// The primary grid
    Grid1D grid;

    /// The current time in the simulation
    double currentTime;

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
     * \brief Compute the states on either side of the interface to the left or
     *        right of a given cell.
     * 
     * \param[in] idxInput Which cell to find the possible interface states for
     * \param[in] side Which side of the cell to find the interface states.
     *                 Only options are "left" or "right"
     * \param[out] leftSideOfInterface The state on the left side of the interface.
     *                                 The order within the vector is density,
     *                                 velocity, pressure.
     * \param[out] rightSideOfInterface The state on the right side of the interface.
     *                                  The order within the vector is density,
     *                                  velocity, pressure.
     */
    void interfaceStates(size_t const &idxInput,
                         std::string const &side,
                         std::vector<double> &leftSideOfInterface,
                         std::vector<double> &rightSideOfInterface);

    // TODO Documentation
    void solveRiemann();


    // TODO Documentation
    void conservativeUpdate();

    /*!
     * \brief Update Simulation1D::grid to the new values that have been 
     * computed. Basically this just copies the values of _tempGrid into grid.
     */
    void updateGrid();

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
     * \param[in] limiterKindConstructor Which limiter to use
     * \param[in] saveDir The directory to save the grid to
     */
    Simulation1D(double const &physicalLength,
                 double const &gamma,
                 double const &CFL,
                 size_t const &reals,
                 size_t const &ghosts,
                 std::string const &initialConditionsKind,
                 std::string const &limiterKindConstructor,
                 std::string const &saveDir);
    /*!
     * \brief Destroy the Simulation1D object. Uses the default constructor
     */
    ~Simulation1D()=default;
};