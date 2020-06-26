#include <cmath>

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