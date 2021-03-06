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
    PerfTimer firstInterfaceTimer("First Interface Reconstruction Timer");
    PerfTimer secondInterfaceTimer("Second Interface Reconstruction Timer");
    PerfTimer riemannTimer("Riemann Solver Timer");
    overallTimer.startTimer();

    // ===== Settings ==========================================================
    double const physicalLength        = 1.;
    double const gamma                 = 1.4;
    double const cfl                   = 0.4;
    double const maxTime               = 0.2;
    size_t const numRealCells          = 1000;
    std::string  initialConditionsKind = "sod";
    std::string  boundaryConditions    = "sod";
    std::string  reconstructionKind    = "PLM";
    std::string  limiterKind           = "MC";  // Options: zeroSlope, centerDiff, minMod, or MC
    std::string  riemannSolverKind     = "HLLC";  // Options: "HLLC" & "exact"
    std::string  saveDir               = "../data/";
    // ===== End Settings ======================================================

    // ===== Initialize Simulation Class =======================================
    Simulation1D sim(physicalLength,
                     gamma,
                     cfl,
                     numRealCells,
                     initialConditionsKind,
                     reconstructionKind,
                     limiterKind,
                     riemannSolverKind,
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

        // Set boundary conditions
        sim.updateBoundaries();

        // Compute first order interface states
        firstInterfaceTimer.startTimer();
        sim.interfaceStates("first reconstruction");
        firstInterfaceTimer.stopTimer();

        // Solve the first Riemann problem
        riemannTimer.startTimer();
        sim.solveRiemann();
        riemannTimer.stopTimer();

        // Update the half time step grid
        sim.conservativeUpdate("half time update");

        // Compute interface states from the half time step grid
        secondInterfaceTimer.startTimer();
        sim.interfaceStates("second reconstruction");
        secondInterfaceTimer.stopTimer();

        // Solve Riemann problem on all the interfaces
        riemannTimer.startTimer();
        sim.solveRiemann();
        riemannTimer.stopTimer();

        // Compute conservative update
        sim.conservativeUpdate("full time update");

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
    firstInterfaceTimer.reportStats();
    secondInterfaceTimer.reportStats();
    riemannTimer.reportStats();
}
