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
     * \param initialConditionsKind The type of initial conditions to use. 
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
     * \param primitive The primitive variable to compute the slope of.
     * \param i The cell in which to compute the slope.
     * \return double The limited slope.
     */
    double _slope(std::vector<double> const &primitive,
                 size_t const &i);

    // TODO Documentation
    void _computeEigens(size_t const &idx,
                        std::vector<double> eigVal,
                        std::vector<std::vector<double>> rEigVec,
                        std::vector<std::vector<double>> lEigVec);

public:
    /// The primary grid
    Grid1D grid;

    /*!
     * \brief Compute the time step using the CFL condition
     */
    void computeTimeStep();

    // TODO Documentation
    void interfaceStates(size_t const &idxInput,
                         std::string const &side,
                         std::vector<double> &leftSideOfInterface,
                         std::vector<double> &rightSideOfInterface);

    // TODO Documentation
    void solveRiemann();


    // TODO Documentation
    void computeFluxes();


    // TODO Documentation
    void conservativeUpdate();

    /*!
     * \brief Update Simulation1D::grid to the new values that have been 
     * computed. Basically this just copies the values of _tempGrid into grid.
     */
    void updateGrid();

    /*!
     * \brief Get the Time Step object
     * 
     * \return double The value of the time step
     */
    double getTimeStep() {return _timeStep;};

    /*!
     * \brief Construct a new Simulation1D object.
     * 
     * \details Sets the private const variables (_limiterKind, _physLen, _cflNum,
     * _deltaX). Initializes the grid and the temporary grid and then sets the 
     * initial conditions of the grid.
     * 
     * \param physicalLength Physical length of the simulation in meters.
     * \param gamma The ratio of specific heats
     * \param CFL The CFL Number
     * \param reals The number of real grid cells
     * \param ghosts The number of ghost cells
     * \param initialConditionsKind Which initial conditions to use
     * \param limiterKindConstructor Which limiter to use
     * \param saveDir The directory to save the grid to
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