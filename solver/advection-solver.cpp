/*
================================================================================
 Second order linear advection solver
 Written by Robert Caddy.  Created on May 21, 2020

 Description: 
     This is a basic Second order linear advection solver. See exercises 5.1-5.4 
     in "Introduction to Computational Astrophysical Hydrodynamics" by Michael 
     Zingale.

 Dependencies:
     TODO
================================================================================
*/

#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <cmath>
using namespace std::chrono;
using std::cout;
using std::endl;

#include "helper.h"
#include "initial-conditions.h"
#include "interface.h"
#include "limiters.h"

int main()
{
    // Start timer and open file
    auto start = high_resolution_clock::now();
    std::ofstream outFile("../data/results.csv");

    // Settings
    const double length       = 1.;                                      // Length of problem in meters
    const int    PhysSize     = 1000;                                    //Size of the physical part of the array
    const double deltax       = length / static_cast<double>(PhysSize);  // Length of each cell in meters
    const int numGhosts       = 2;                                       // Number of ghost cells
    const int size            = PhysSize + 2 * numGhosts;                // total size of the array
    const double CFLNum       = 0.4;                                     // CFL Number
    const double vel          = 1.;                                      // Velocity in meters/second
    const double period       = length / std::abs(vel);                  // Time for one period
    const double maxTime      = 1.*period;                               // Time to simlate to
    const std::string LimType = "MC";                                    // Limiter Type

    // Conserved quantity
    std::vector<double> a(size); // Actual array
    std::vector<double> aTemp(size); // Array to store the updates in

    // Set initial conditions
    setInitialConditions(a, size, "top-hat");
    saveArray(a, outFile, numGhosts);

    //=== Begin the main evolution loop ========================================
    int step = 0;
    double time = 0;
    while (time <= maxTime)
    {
        // Compute the time step using the CFL condition
        double deltat = CFLNum * deltax / std::abs(vel);

        for (int i = numGhosts; i < (size-numGhosts); i++)
        {
            // Set boundary conditions (periodic)
            for (int j = 0; j < numGhosts; j++)
            {
                a[j] = a.end()[-(2*numGhosts-j)];
                a.end()[-(numGhosts-j)] = a[j+numGhosts];
            }

            // Computer interface states and solve Riemann problem
            double LeftInterface;
            double RightInterface;
            if (vel >= 0.)
            {
                double derivA = SlopeLimiter(a[i - 1], a[i], a[i + 1], deltax, LimType);
                RightInterface = a[i] + (deltax / 2) * (1 - (deltat / deltax) * vel) * derivA;

                derivA = SlopeLimiter(a[i - 2], a[i - 1], a[i], deltax, LimType);
                LeftInterface = a[i - 1] + (deltax / 2) * (1 - (deltat / deltax) * vel) * derivA;
            }
            else// i.e. (vel < 0.)
            {
                double derivA = SlopeLimiter(a[i], a[i + 1], a[i + 2], deltax, LimType);
                RightInterface = a[i+1] - (deltax/2) * (1+(deltat/deltax)*vel) * derivA;

                derivA = SlopeLimiter(a[i - 1], a[i], a[i + 1], deltax, LimType);
                LeftInterface = a[i] - (deltax / 2) * (1 + (deltat / deltax) * vel) * derivA;
            };

            // Compute conservative update
            double LeftFlux  = LeftInterface * vel;
            double RightFlux = RightInterface * vel;

            double FluxDerivative = (LeftFlux - RightFlux)/deltax;
            aTemp[i] = (FluxDerivative * deltat) + a[i];
        }; // End of loop to interate through array

        // Copy values from aTemp to a
        for (int i = numGhosts; i < (size - numGhosts); i++)
        {
            a[i] = aTemp[i];
        };

        // Save
        saveArray(a, outFile, numGhosts);

        // Message
        cout << "Completeted step: " << step << endl;

        // Update time and step number
        time += deltat;
        step++;

    }; // End of evolution loop


    // Close output file
    outFile.close();

    // Stop timer and print execution time. Time options are nanoseconds, 
    // microseconds, milliseconds, seconds, minutes, hours. To pick one just
    // change `using FpTime` and the cout statement suitably.
    auto stop = high_resolution_clock::now();
    using FpTime = duration<float, seconds::period>;
    static_assert(treat_as_floating_point<FpTime::rep>::value, "Rep required to be floating point");
    auto duration = FpTime(stop - start);
    cout << "Time to execute: " << duration.count() << " seconds";
    return 0;
}
