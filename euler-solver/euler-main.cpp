/*******************************************************************************
 * \file euler-main.cpp
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Solves the Euler equations in 1D using a second order Gudonov Scheme
 *
 * \date 2020-06-23
 *
 * \copyright Copyright (c) 2020
 *
 * \details This program serves as a place for me to learn computational
 * hydrodynamics and as a testbed for my future additions to
 * [Cholla](https://github.com/cholla-hydro/cholla).
 *
 ******************************************************************************/

#include <iostream>
#include <string>
using std::cout;
using std::cin;
using std::endl;

#include "Grid1D.h"
#include "Simulation1D.h"
#include "PerfTimer.h"

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
    // Initialize timers and start overall timer
    PerfTimer overallTimer("Overall Timer");
    PerfTimer timeStepTimer("Time Step Timer");
    PerfTimer interfaceTimer("Interface Reconstruction Timer");
    PerfTimer riemannTimer("Riemann Solver Timer");
    overallTimer.startTimer();

    // ===== Settings ==========================================================
    double const physicalLength        = 1.;
    double const gamma                 = 1.4;
    double const cfl                   = 0.8;
    double const maxTime               = 0.2;
    size_t const numRealCells          = 100;
    size_t const numGhostCells         = 2;
    std::string  initialConditionsKind = "sod";
    std::string  boundaryConditions    = "sod";
    std::string  reconstructionKind    = "PLM";
    std::string  limiterKind           = "MC"; // Options: zeroSlope, centerDiff, minMod, or MC
    std::string  saveDir               = "../data/";
    // ===== End Settings ======================================================

    // ===== Initialize Simulation Class =======================================
    Simulation1D sim(physicalLength,
                     gamma,
                     cfl,
                     numRealCells,
                     numGhostCells,
                     initialConditionsKind,
                     reconstructionKind,
                     limiterKind,
                     boundaryConditions,
                     saveDir);
    // ===== End initializing Simulation Class =================================

    // Save the initial state
    sim.grid.saveState();

    //=== Begin the main evolution loop ========================================
    size_t step = 1;
    while (sim.currentTime <= maxTime)
    {
        // Compute the time step using the CFL condition
        timeStepTimer.startTimer();
        sim.computeTimeStep();
        timeStepTimer.stopTimer();

        // Set boundary conditions (sod)
        sim.grid.updateBoundaries(gamma);

        // Compute interface states.
        // note that the order within vectors is density, velocity, pressure
        interfaceTimer.startTimer();
        sim.interfaceStates();
        interfaceTimer.stopTimer();

        // Solve Riemann problem on all the interfaces
        riemannTimer.startTimer();
        sim.solveRiemann();
        riemannTimer.stopTimer();

        // Compute conservative update
        sim.conservativeUpdate();

        // Save output
        sim.grid.saveState();

        // Update current time
        sim.updateCurrentTime();

        // Message
        cout << "Completeted step: " << step
        << ",   Time step = " << sim.getTimeStep()
        << ",   Simulation Time = " << sim.currentTime << endl;

        // Increment the step number
        step++;
    };
    // ===== End of evolution loop =============================================


    // Print timing results
    overallTimer.stopTimer();
    overallTimer.reportStats();
    timeStepTimer.reportStats();
    interfaceTimer.reportStats();
    riemannTimer.reportStats();
}