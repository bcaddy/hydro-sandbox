/*
================================================================================
 Second order linear solver for Burgers' Equation with minmod limiting
 Written by Robert Caddy.  Created on May 21, 2020

 Description: 
     This is a basic Second order Burgers' Equation solver. See exercises 6.2
     in "Introduction to Computational Astrophysical Hydrodynamics" by Michael 
     Zingale.

 Dependencies:
     helper.cpp
     stl
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

int main()
{
    // Start timer and open file
    auto start = high_resolution_clock::now();
    std::ofstream outFile("../data/results.csv");

    // Setup Initial conditions
    const double length   = 1.;                                      // Length of problem in meters
    const int    PhysSize = 256;                                     // Size of the physical part of the array
    const double deltax   = length / static_cast<double>(PhysSize);  // Length of each cell in meters
    const int numGhosts   = 2;                                       // Number of ghost cells
    const int size        = PhysSize + 2 * numGhosts;                // total size of the array
    const double CFLNum   = 0.8;                                     // CFL Number
    const double maxTime  = .2;                                      // Time to simlate to

    // Conserved quantity
    std::vector<double> uVel(size);      // Actual array
    std::vector<double> uVelTemp(size);  // Array to store the updates in

    // Set initial conditions
    setInitialConditions(uVel, size, "vel-step");
    saveArray(uVel, outFile, numGhosts);

    //=== Begin the main evolution loop ========================================
    int step = 0;
    double time = 0;
    while (time <= maxTime)
    {
        // Compute the time step using the CFL condition
        double deltat = std::abs(uVel[numGhosts]);
        for (int i = numGhosts+1; i < (size - numGhosts); i++)
        {
            double deltatTemp = CFLNum * deltax / std::abs(uVel[i]);
            if (deltat > deltatTemp)
            {
                deltat = deltatTemp;
            }
        }
        cout << deltat << endl;
        for (int i = numGhosts; i < (size-numGhosts); i++)
        {
            // Set boundary conditions (periodic)
            for (int j = 0; j < numGhosts; j++)
            {
                uVel[j] = uVel.end()[-(2*numGhosts-j)];
                uVel.end()[-(numGhosts-j)] = uVel[j+numGhosts];
            }

            // Computer interface states and solve Riemann problem
            // TODO begin
            double LeftInterface  = VelInterface(uVel[i-2], 
                                                 uVel[i-1], 
                                                 uVel[i],
                                                 uVel[i+1], 
                                                 deltax,
                                                 deltat);
            double RightInterface = VelInterface(uVel[i-1],
                                                 uVel[i], 
                                                 uVel[i+1],
                                                 uVel[i+2], 
                                                 deltax,
                                                 deltat);

            // Compute conservative update
            double LeftFlux  = 0.5 * std::pow(LeftInterface, 2.) ;
            double RightFlux = 0.5 * std::pow(RightInterface, 2.) ;

            double FluxDerivative = (LeftFlux - RightFlux)/deltax;
            uVelTemp[i] = (FluxDerivative * deltat) + uVel[i];
        }; // End of loop to interate through array

        // Copy values from uVelTemp to a
        for (int i = numGhosts; i < (size - numGhosts); i++)
        {
            uVel[i] = uVelTemp[i];
        };
        
        // Save output
        saveArray(uVel, outFile, numGhosts);

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
