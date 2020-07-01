#include <cmath>
#include <stdexcept>

#include "Simulation1D.h"

// =============================================================================
void Simulation1D::_setInitialConditions(std::string const &initialConditionsKind)
{
    // Set up Sod initial conditions
    // TODO Add other kinds of initial conditions

    size_t const half  = grid.numTotCells / 2;

    // Iterate over just the real cells on the left side
    for (size_t i = grid.numGhostCells; 
         i < half; 
         i++)
    {
        grid.velocity[i] = 0.;
        grid.density[i]  = 1.;
        grid.pressure[i] = 1.;
    }

    // Iterate over the real cells on the right side
    for (size_t i = half;
         i < (grid.numTotCells - grid.numGhostCells);
         i++)
    {
        grid.velocity[i] = 0;
        grid.density[i]  = 0.125;  // 1/8th
        grid.pressure[i] = 0.1;
    }
}
// =============================================================================

// =============================================================================
double Simulation1D::_slope(std::vector<double> const &primitive,
              size_t const &i)
{
    // MC limiter
    if (_limiterKind == "MC")
    {
        // Declare variables
        double outValue;
        double leftDerive, rightDerive, leftRightProduct;

        // Compute the derivatives
        leftDerive =  (primitive[i] - primitive[i-1]) / _deltaX;
        rightDerive = (primitive[i + 1] - primitive[i]) / _deltaX;
        leftRightProduct = leftDerive * rightDerive;

        // Choose what value to output
        if ((std::abs(leftDerive) < std::abs(rightDerive)) && (leftRightProduct > 0.))
        {
            outValue = leftDerive;
        }
        else if ((std::abs(leftDerive) > std::abs(rightDerive)) && (leftRightProduct > 0.))
        {
            outValue = rightDerive;
        }
        else
        {
            outValue = 0.;
        }

        return outValue;
    }
    else
    {
        throw std::invalid_argument("Nonexistant limiter kind used.");
    }
    
}
// =============================================================================

// =============================================================================
void Simulation1D::computeTimeStep()
{
    // I don't want to compute this in every iteration
    double const cflTimesDeltaX = _cflNum * _deltaX;

    // Set the timestep to the value determined by the first real cell
    _timeStep = cflTimesDeltaX / std::abs(grid.velocity[grid.numGhostCells]);

    // Go through the entire grid, compute the time step for each cell, and
    // choose the smallest one by setting _timeStep equal to it.
    for (size_t i = grid.numGhostCells + 1;
         i < (grid.numTotCells - grid.numGhostCells); 
         i++)
    {
        double deltatTemp = cflTimesDeltaX / std::abs(grid.velocity[i]);
        if (_timeStep > deltatTemp)
        {
            _timeStep = deltatTemp;
        }
    }
}
// =============================================================================

// =============================================================================
void Simulation1D::updateGrid()
{
    // Copy every real element in Simulation1D::_tempGrid to Simulation1D::grid
    grid.velocity = _tempGrid.velocity;
    grid.density  = _tempGrid.density;
    grid.pressure = _tempGrid.pressure;
    grid.siEnergy = _tempGrid.siEnergy;
}
// =============================================================================

// =============================================================================
// Constructor
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
// =============================================================================