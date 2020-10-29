/*!
* \file test-Riemann.cpp
* \author Robert 'Bob' Caddy (rvc@pitt.edu)
* \brief Test the exact Riemann solver from Toro using a Sod Shock tube
* \version 0.1
* \date 2020-07-23
*
* \copyright Copyright (c) 2020
*
*/

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <vector>
#include <memory>
#include "RiemannSolver.h"
#include "HllcRiemannSolver.h"
#include "ExactRiemannSolver.h"
#include "Grid1D.h"
#include "HydroUtilities.h"

using std::cout;
using std::cin;
using std::endl;
using namespace HydroUtilities;

int main()
{
    // Start timer
    auto start = std::chrono::high_resolution_clock::now();


    // ===== Settings for Sod Shock Tube =======================================
    double const densityR   = 1.0;
    double const velocityR  = 0.0;
    double const pressureR  = 1.0;
    double const densityL   = 0.1;
    double const velocityL  = 0.0;
    double const pressureL  = 0.1;
    double const gamma      = 1.4;
    double const xSize      = 1.0;
    double const tMax       = 0.2;
    size_t const numSamples = 100;
    double const diaphPos   = 0.5 * xSize;
    double const deltaX     = xSize / static_cast<double>(numSamples);
    std::string  solverKind = "exact";
    // ===== End Settings for Sod Shock Tube  ==================================

    // ===== Grid & Solver ===================================================
    Grid1D grid(numSamples, 0, "pass", "../data/");

    std::unique_ptr<RiemannSolver> riemannSolver = std::unique_ptr<RiemannSolver>(new ExactRiemannSolver(gamma));
    // std::unique_ptr<RiemannSolver> riemannSolver = std::unique_ptr<RiemannSolver>(new HllcRiemannSolver(gamma));
    // ===== End Arrays & Solver ===============================================

    // ===== Compute the Solution ==============================================
    for (size_t i = 0; i < numSamples; i++)
    {
        // Compute posOverT
        double position  = (static_cast<double>(i) - 0.5) * deltaX;
        double posOverT  = (position - diaphPos) / tMax;

        // Run solver
        double denFlux, momFlux, eneFlux;  // Variables for fluxes that I'm not using
        riemannSolver->riemannMain(densityR,
                                   velocityR,
                                   pressureR,
                                   densityL,
                                   velocityL,
                                   pressureL,
                                   denFlux,
                                   momFlux,
                                   eneFlux,
                                   posOverT);

        // Get the state values
        grid.density[i]  = riemannSolver->getDensityState();
        grid.momentum[i] = computeMomentum(riemannSolver->getVelocityState(),
                                           riemannSolver->getDensityState());
        grid.energy[i]   = computeEnergy(riemannSolver->getPressureState(),
                                         riemannSolver->getDensityState(),
                                         riemannSolver->getVelocityState(), gamma);
    }
    // ===== Done Computing the Solution =======================================

    // ===== Save the arrays ===================================================
    // save a bunch so visualizer works properly
    for (size_t i = 0; i < 10; i++)
    {
        grid.saveState();
    }

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