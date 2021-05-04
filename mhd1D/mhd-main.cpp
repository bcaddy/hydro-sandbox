/*******************************************************************************
 * \file mhd-main.cpp
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Solves the ideal MHD in their Eulerian form in 1D using the VL+CT
 * algorithm from Stone & Gardiner 2009 and the HLLD Riemann solver
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

#include "MhdSimulation1D.h"
#include "PerfTimer.h"

/*!
 * \brief Solves the ideal MHD in their Eulerian form in 1D using the VL+CT
 * algorithm from Stone & Gardiner 2009 and the HLLD Riemann solver. This main
 * function that invokes class methods and provides user output.
 *
 * \details The main function provides all, or at least most, of the calls to class
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
    PerfTimer ctFieldsTimer("CT Fields Timer");
    PerfTimer firstInterfaceTimer("First Interface Reconstruction Timer");
    PerfTimer secondInterfaceTimer("Second Interface Reconstruction Timer");
    PerfTimer riemannTimer("Riemann Solver Timer");
    overallTimer.startTimer();

    // ===== Settings ==========================================================
    double const physicalLength        = 1.;
    double const gamma                 = 5./3.;  // 5./3. for most things, chollaSod uses 1.4
    double const cfl                   = 0.4;
    double const maxTime               = 1.0; //0.1 for B&W shock tube, 0.2 for D&W and sod
    size_t const numRealCells          = 10;
    std::string  initialConditionsKind = "singleWaveCR";
    std::string  boundaryConditions    = "periodic";
    std::string  reconstructionKind    = "PLM";
    std::string  limiterKind           = "MC";  // Options: zeroSlope, centerDiff, minMod, or MC
    std::string  riemannSolverKind     = "HLLD";  // Options: "HLLD"
    std::string  saveDir               = "../data/";
    // ===== End Settings ======================================================

    // ===== Initialize Simulation Class =======================================
    MhdSimulation1D sim(physicalLength,
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
        // Set boundary conditions
        sim.updateBoundaries();

        // Compute the time step using the CFL condition
        timeStepTimer.startTimer();
        sim.computeTimeStep();
        timeStepTimer.stopTimer();

        // Compute first order interface states
        firstInterfaceTimer.startTimer();
        sim.interfaceStates("first reconstruction");
        firstInterfaceTimer.stopTimer();

        // Solve the first Riemann problem
        riemannTimer.startTimer();
        sim.solveRiemann();
        riemannTimer.stopTimer();

        // Compute the CT Fields
        ctFieldsTimer.startTimer();
        sim.ctElectricFields("Half time");
        ctFieldsTimer.stopTimer();

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

        // Compute the CT Fields
        ctFieldsTimer.startTimer();
        sim.ctElectricFields("Full time");
        ctFieldsTimer.stopTimer();

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
    ctFieldsTimer.reportStats();
}
