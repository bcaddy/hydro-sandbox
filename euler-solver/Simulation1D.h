/*!
 * \file Simulation1D.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains the Simulation1D class for solving the 1D Euler equations.
 * \date 2020-06-25
 * 
 * \copyright Copyright (c) 2020
 * 
 */


#include <string>
#include "Grid1D.h"

class Simulation1D
{
private:
    std::string const _limiterKind;

    double const _physLen;
    double const _cflNum;
    double const _deltaX;
    double _timeStep;

    Grid1D _tempGrid;

    void _setInitialConditions(std::string const &initialConditionsKind);

    double _slope(std::vector<double> const &primitive,
                 size_t const &i);

public:
    Grid1D grid;

    void computeTimeStep();

    void interfaceStates();

    void solveRiemann();

    void computeFluxes();

    void updatedGrid();

    double getTimeStep() {return _timeStep;};


    Simulation1D(double const &physicalLength,
                 double const &CFL,
                 size_t const &reals,
                 size_t const &ghosts,
                 std::string const &initialConditionsKindConstructor,
                 std::string const &limiterKindConstructor,
                 std::string const &saveDirConstructor);
    ~Simulation1D()=default;
};

Simulation1D::Simulation1D(double const &physicalLength,
                           double const &CFL,
                           size_t const &reals,
                           size_t const &ghosts,
                           std::string const &initialConditionsKind,
                           std::string const &limiterKindConstructor,
                           std::string const &saveDir)

    // Start by initializing all the const member variables
    : _limiterKind(limiterKindConstructor),
      _physLen(physicalLength),
      _cflNum(CFL),
      _deltaX(_physLen / static_cast<double>(reals))
{
    // Now we initialize the grids.
    grid.init(reals, ghosts, saveDir);
    _tempGrid.init(reals, ghosts, "no saving");

    // And set the initial conditions
    _setInitialConditions(initialConditionsKind);
}