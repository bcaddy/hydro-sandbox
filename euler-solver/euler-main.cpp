/*******************************************************************************
 * \file euler-main.cpp
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Solves the Euler equations in 1D using a second order Gudonov Scheme
 *
 * \date 2020-06-23
 *
 * \copyright Copyright (c) 2020
 *
 * \details This program serves as a place for me to learn computational hydrodynamics
 * and as a testbed for my future additions to [Cholla](https://github.com/cholla-hydro/cholla).
 *
 ******************************************************************************/

#include <iostream>
#include <string>
#include <chrono>
using std::cout;
using std::cin;
using std::endl;

#include "Grid1D.h"
#include "Simulation1D.h"

/*!
 * \brief Main function that invokes class methods and provides user output.
 *
 * The main function provides all, or at least most, of the calls to class
 * methods. The classes used here don't do much calling between member functions
 * and so this function handles most if this. That is primarily for readability.
 *
 * Main also handles standard output and timing to inform the user about what is
 * going on in the code.
 *
 * \return 0 on a successful run
 */
int main()
{
    // Start clock
    auto start = std::chrono::high_resolution_clock::now();

    // ===== Settings ==========================================================
    double const physicalLength        = 100.;
    double const gamma                 = 5./3.;
    double const cfl                   = 0.8;
    double const maxTime               = 10.;
    size_t const numRealCells          = 1000;
    size_t const numGhostCells         = 2;
    std::string  initialConditionsKind = "sod";
    std::string  boundaryConditions    = "sod";
    std::string  saveDir               = "../data/";
    // ===== End Settings ======================================================

    // ===== Initialize Simulation Class =======================================
    Simulation1D sim(physicalLength,
                     gamma,
                     cfl,
                     numRealCells,
                     numGhostCells,
                     initialConditionsKind,
                     boundaryConditions,
                     saveDir);

    //=== Begin the main evolution loop ========================================
    size_t step = 0;

    while (sim.currentTime <= maxTime)
    {
        // Compute the time step using the CFL condition
        sim.computeTimeStep();

        for (size_t i = sim.grid.numGhostCells;
             i < (sim.grid.numTotCells-sim.grid.numGhostCells);
             i++)
        {
            // Set boundary conditions (periodic)
            sim.grid.updateBoundaries();

            // Compute interface states on the left side.
            //   note that the order within vectors is density, velocity, pressure
            std::vector<double> leftSideOfInterface, rightSideOfInterface;
            sim.interfaceStates(i,
                                "left",
                                leftSideOfInterface,
                                rightSideOfInterface);

            // Solve Riemann problem on the left side
            double leftEnergyFlux, leftMomentumFlux, leftMassFlux;
            sim.solveRiemann(rightSideOfInterface[0],
                             rightSideOfInterface[1],
                             rightSideOfInterface[2],
                             leftSideOfInterface[0],
                             leftSideOfInterface[1],
                             leftSideOfInterface[2],
                             0.0, // position over t
                             leftEnergyFlux,
                             leftMomentumFlux,
                             leftMassFlux);

            // Compute interface states on the right side.
            //   note that the order within vectors is density, velocity, pressure
            sim.interfaceStates(i,
                                "right",
                                leftSideOfInterface,
                                rightSideOfInterface);

            // Solve Riemann problem on the left side
            double rightEnergyFlux, rightMomentumFlux, rightMassFlux;
            sim.solveRiemann(rightSideOfInterface[0],
                             rightSideOfInterface[1],
                             rightSideOfInterface[2],
                             leftSideOfInterface[0],
                             leftSideOfInterface[1],
                             leftSideOfInterface[2],
                             0.0, // position over t
                             rightEnergyFlux,
                             rightMomentumFlux,
                             rightMassFlux);

            // Compute conservative update
            sim.conservativeUpdate(i,
                                   leftMassFlux,
                                   leftMomentumFlux,
                                   leftEnergyFlux,
                                   rightMassFlux,
                                   rightMomentumFlux,
                                   rightEnergyFlux);
        }; // End of loop to interate through array

        // Copy values from uVelTemp to a
        sim.updateGrid();

        // Save output
        sim.grid.saveState();

        // Message
        cout << "Completeted step: " << step
        << ":   Time step = " << sim.getTimeStep()
        << ":   Time = " << sim.currentTime << endl;

        // Update time and step number
        sim.updateCurrentTime();
        step++;

    };
    // ===== End of evolution loop =============================================


    // Stop timer and print execution time. Time options are nanoseconds,
    // microseconds, milliseconds, seconds, minutes, hours. To pick one just
    // change `using FpTime` and the cout statement suitably.
    auto stop = std::chrono::high_resolution_clock::now();
    using FpTime = std::chrono::duration<float, std::chrono::seconds::period>;
    static_assert(std::chrono::treat_as_floating_point<FpTime::rep>::value,"Rep required to be floating point");
    auto duration = FpTime(stop - start);
    cout << "Time to execute: " << duration.count() << " seconds";
    return 0;
}