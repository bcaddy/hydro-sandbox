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
    if (limiterKind == "MC")
    {
        // Declare variables
        double outValue;
        double leftDerive, rightDerive, leftRightProduct;

        // Compute the derivatives
        leftDerive =  (primitive[i] - primitive[i-1]) / deltaX;
        rightDerive = (primitive[i + 1] - primitive[i]) / deltaX;
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
    timeStep = std::abs(grid.velocity[grid.numGhostCells]);

    for (int i = grid.numGhostCells + 1;
         i < (grid.numTotCells - grid.numGhostCells); 
         i++)
    {
        double deltatTemp = CFLnum * deltaX / std::abs(grid.velocity[i]);
        if (timeStep > deltatTemp)
        {
            timeStep = deltatTemp;
        }
    }
}
// =============================================================================