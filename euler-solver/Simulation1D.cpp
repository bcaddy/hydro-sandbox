#include <cmath>

#include "Simulation1D.h"

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