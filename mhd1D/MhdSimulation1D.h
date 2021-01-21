/*!
 * \file MhdSimulation1D.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains the MhdSimulation1D class for solving the 1D Euler equations.
 * \date 2020-06-25
 *
 * \copyright Copyright (c) 2020
 *
 */
#pragma once

#include <string>
#include <vector>
#include <memory>

#include "Grid1D.h"
#include "PrimitiveGrid1D.h"
#include "MhdRiemannSolver.h"

/*!
 * \brief Solves the Euler equations
 * \details This class handles all the heavy lifting of running a 1D hydro
 * simulation. It uses the Grid1D class for the grid and then computes time steps,
 * interface steps, solves the Riemann problem using the RiemannSolver class,
 * does the conservative update, keeps track of the time, etc.
 *
 */
class MhdSimulation1D
{
private:
    /// The physical length of the simulation in meters
    double const _physLen;
    /// Courant–Friedrichs–Lewy (CFL) Number
    double const _cflNum;
    /// The physical size of each cell in meters
    double const _deltaX;
    /// The physical size of each cell in meters
    double const _deltaY;
    /// The physical size of each cell in meters
    double const _deltaZ;
    /// The ratio of specific heats $\gamma$
    double const _gamma;
    /// The time step for each interation
    double _timeStep;

    /// The number of ghost cells
    size_t _numGhosts = 3;

    /// The cell centered magnetic fields
    std::vector<std::vector<double>> _magCentered;

    /// The vector to store the primitive variables at the left side of the
    /// interface. _interfaceL.member[i] is the state on the left side
    /// of the i-1/2 interface
    PrimitiveGrid1D _interfaceL;

    /// The vector to store the primitive variables at the right side of the
    /// interface. _interfaceR.member[i] is the state on the right side
    /// of the i-1/2 interface
    PrimitiveGrid1D _interfaceR;

    /// The vector to store the fluxes. _flux.variable[i] is the flux through
    /// the i-1/2 interface
    Grid1D _flux;

    /// The grid to store the half time step conserved variables
    Grid1D _gridHalfTime;

    /// The 4D vector of shape Nx3x3x3 to store the 3D edge fields in each
    /// cell. These are the i-1/2 edges pointing in the x, y, and z directions
    /// respectively
    std::vector<std::vector<std::vector<std::vector<double>>>> _edgeFields;

    /// The 5D vector used to hold the velocities from the Riemann solve.
    /// Indices are as follows
    /// [x position]
    /// [y position]
    /// [z position]
    /// [which face, i-1/2, j-1/2, or z-1/2 in that order]
    /// [the velocity in each direction, x, y, and z respectively]
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> _ctVelocities;

    /// The object used to solve the Riemann Problem. Uses the RiemannSolver
    /// virtual base class and chooses between different Riemann Solver
    /// algorithms
    std::unique_ptr<MhdRiemannSolver> _riemannSolver;

    /*!
     * \brief Set the initial conditions. Currently only supports a Sod shock
     * tube.
     *
     * \param[in] initialConditionsKind The type of initial conditions to use.
     *            Options are:
     *            | Keyword        | Initial Conditions                                                       |
     *            |----------------|--------------------------------------------------------------------------|
     *            | dwShockTube    | MHD Shock tube from Dai & Woodward 1994                                  |
     *            | bwShockTube    | MHD Shock tube from Brio & Wu 1988                                       |
     *            | indexCheck     | Set each primitive variable to the value of the grid index at that point |
     *            | advectionStep  | Step function advection                                                  |
     *            | advectionGauss | Gaussion function advection                                              |
     */
    void _setInitialConditions(std::string const &initialConditionsKind);

    /*!
     * \brief Compute the cell centered magnetic field
     *
     * \param[in] activeGrid The grid to compute centered fields for
     */
    void _centeredMagneticField(Grid1D const &activeGrid);

    /*!
     * \brief Slope of a given primitive variable
     *
     * \details Computes the slope. Uses MhdSimulation1D::limiterKind to choose
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
     * \brief Compute the Constrained Transport electric fields using the
     * algorithm from Stone & Gardiner 2009
     *
     * \param[in] activeGrid The grid to compute the CT fields for
     */
    void _ctElectricFields(Grid1D const &activeGrid);

    /*!
     * \brief Return the mathematical modulo of x mod 3. Used extensively in CT
     * implementations
     *
     * \param[in] x The number to modulo with 3
     * \return int x modulo 3
     */
    int _mod3(int const &x)
        {return (x % 3 + 3) % 3;};

    // TODO add documentation once this is figured out
    double _ctSlope();

    /*!
     * \brief Computes the interface states using the Piecewise Linear Method
     *        (PLM) from "Introduction to Computational Astrophysical
     *        Hydrodynamics" by Zingale, git version: 4de1fef51af5 Section 8.2.2
     *
     * \param[in] workingGrid The grid object that will be used for the
     * reconstruction
     */
    void _piecewiseLinearReconstruction(Grid1D const &workingGrid);

    /*!
     * \brief Computes the interface states using the Piecewise Constant Method
     * (PCM)
     *
     * \param[in] workingGrid The grid object that will be used for the
     * reconstruction
     */
    void _piecewiseConstantReconstruction(Grid1D const &workingGrid);

public:
    /// The primary grid
    Grid1D grid;

    /// The current time in the simulation
    double currentTime;

    /// The reconstruction scheme to use for the half time step reconstruction
    std::string const reconstructionKind;

    /// The kind of slope limiter to use
    std::string const limiterKind;

    /// The Riemann solver to use
    std::string const riemannSolverKind;

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
     * \brief Increases the current time by the value of MhdSimulation1D::_timeStep
     *
     */
    void updateCurrentTime() { currentTime += _timeStep; };

    /*!
     * \brief Update the ghost cells in both grids
     *
     */
    void updateBoundaries()
    {
        grid.updateBoundaries(_gamma);
        _gridHalfTime.updateBoundaries(_gamma);
    };

    /*!
     * \brief Compute all the interface states.
     *
     * \details Calls either MhdSimulation1D::_piecewiseConstant or
     *          MhdSimulation1D::_piecewiseLinear depending on if
     *          MhdSimulation1D::reconstructionKind is "PCM" or "PCL" respectively
     *          or just calls MhdSimulation1D::_piecewiseConstant for the first
     *          reconstruction in the Van Leer algorithm
     *
     * \param[in] algoStep What part of the Van Leer algorithm to compute the
     * interface states for. Options are "first reconstruction" for the zero
     * order first reconstruction or "second reconstruction" for the higher
     * order, half time step reconstruction
     */
    void interfaceStates(std::string const &algoStep);

    /*!
     * \brief Solves the Riemann problem by calling the main function of
     * the RiemannSolver class
     */
    void solveRiemann();


    /*!
     * \brief Performs the conservative update. Through the input parameter the
     * user can decide whether to update the half time step grid using values
     * from the real grid or to update the real grid a full time step using
     * values from the real grid.
     *
     * \param[in] timeChoice Choose whether to update the main grid for the full
     * time step or the temporary grid at the half time step. Options are "full
     * time update" or "half time update" respectively.
     */
    void conservativeUpdate(std::string const &timeChoice);

    /*!
     * \brief Construct a new MhdSimulation1D object.
     *
     * \details Sets the private const variables (_limiterKind, _physLen, _cflNum,
     * _deltaX). Initializes the grid and the temporary grid and then sets the
     * initial conditions of the grid.
     *
     * \param[in] physicalLength Physical length of the simulation in meters.
     * \param[in] gamma The ratio of specific heats
     * \param[in] CFL The CFL Number
     * \param[in] reals The number of real grid cells
     * \param[in] initialConditionsKind Which initial conditions to use
     * \param[in] reconstructionKind Which kind of interface reconstruction to
     *            use. Option are "PCM" for Piecewise Constant Method and "PLM"
     *            for Piecewise Linear Method
     * \param[in] limiterKind What kind of limiter to use. Options are
     *            "centerDiff", "minMod", and "MC"
     * \param[in] riemannSolverKind What kind of MHD Riemann solver to use.
     *            The only option currently is "HLLD"
     * \param[in] boundaryConditions Which kind of boundary conditions to use
     * \param[in] saveDir The directory to save the grid to
     */
    MhdSimulation1D(double const &physicalLength,
                    double const &gamma,
                    double const &CFL,
                    size_t const &reals,
                    std::string const &initialConditionsKind,
                    std::string const &reconstructionKind,
                    std::string const &limiterKind,
                    std::string const &riemannSolverKind,
                    std::string const &boundaryConditions,
                    std::string const &saveDir);
    /*!
     * \brief Destroy the MhdSimulation1D object. Uses the default constructor
     */
    ~MhdSimulation1D()=default;
};