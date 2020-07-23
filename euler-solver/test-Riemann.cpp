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
#include "RiemannSolver.h"
using std::cout;
using std::cin;
using std::endl;

int main()
{
    // Start timer
    auto start = std::chrono::high_resolution_clock::now();


    // ===== Settings for Sod Shock Tube =======================================
    double const densityL   = 1.0;
    double const velocityL  = 0.0;
    double const pressureL  = 1.0;
    double const densityR   = 0.125;
    double const velocityR  = 0.0;
    double const pressureR  = 0.1;
    double const gamma      = 1.4;
    double const xSize      = 1.;
    double const tMax       = 0.2;
    size_t const numSamples = 1000;
    double const diaphPos   = 0.5 * xSize;
    double const deltaX     = xSize / static_cast<double>(numSamples);
    // ===== End Settings for Sod Shock Tube  ===================================

    // // ===== Settings for Sod Shock Tube Swapped ===============================
    // double const densityR   = 1.0;
    // double const velocityR  = 0.0;
    // double const pressureR  = 1.0;
    // double const densityL   = 0.125;
    // double const velocityL  = 0.0;
    // double const pressureL  = 0.1;
    // double const gamma      = 1.4;
    // double const xSize      = 1.;
    // double const tMax       = 0.2;
    // size_t const numSamples = 1000;
    // double const diaphPos   = 0.5 * xSize;
    // double const deltaX     = xSize / static_cast<double>(numSamples);
    // // ===== End Settings for Sod Shock Tube Swapped ===========================

    // ===== Arrays & Solver ===================================================
    std::vector<double> density(numSamples);
    std::vector<double> pressure(numSamples);
    std::vector<double> velocity(numSamples);

    RiemannSolver solver(gamma);
    // ===== End Arrays & Solver ===============================================

    // ===== Compute the Solution ==============================================
    for (size_t i = 0; i < numSamples; i++)
    {
        // Compute posOverT
        double position  = (static_cast<double>(i) - 0.5) * deltaX;
        double posOverT  = (position - diaphPos) / tMax;

        // Run solver
        double junk1, junk2, junk3;  // Variables for fluxes that I'm not using
        solver.riemannMain(densityR,
                           velocityR,
                           pressureR,
                           densityL,
                           velocityL,
                           pressureL,
                           posOverT,
                           junk1,
                           junk2,
                           junk3);

        // Get the state values
        density[i]  = solver.getDensityState();
        velocity[i] = solver.getVelocityState();
        pressure[i] = solver.getPressureState();
    }
    // ===== Done Computing the Solution =======================================

    // ===== Save the arrays ===================================================
    // save twice so visualizer works properly
    std::string const savePath = "../data/";

    std::ofstream densitySaveFile ;
    std::ofstream velocitySaveFile;
    std::ofstream pressureSaveFile;

    densitySaveFile.open(savePath + "Density.csv");
    velocitySaveFile.open(savePath + "Velocity.csv");
    pressureSaveFile.open(savePath + "Pressure.csv");

    velocitySaveFile << velocity[0];
    densitySaveFile  << density[0];
    pressureSaveFile << pressure[0];

    for (size_t i = 1; i < numSamples; i++)
    {
        velocitySaveFile << "," << velocity[i];
        densitySaveFile  << "," << density[i];
        pressureSaveFile << "," << pressure[i];
    }
    velocitySaveFile << std::endl;
    densitySaveFile  << std::endl;
    pressureSaveFile << std::endl;

    velocitySaveFile << velocity[0];
    densitySaveFile  << density[0];
    pressureSaveFile << pressure[0];

    for (size_t i = 1; i < numSamples; i++)
    {
        velocitySaveFile << "," << velocity[i];
        densitySaveFile  << "," << density[i];
        pressureSaveFile << "," << pressure[i];
    }
    velocitySaveFile << std::endl;
    densitySaveFile  << std::endl;
    pressureSaveFile << std::endl;


    densitySaveFile.close();
    velocitySaveFile.close();
    pressureSaveFile.close();



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