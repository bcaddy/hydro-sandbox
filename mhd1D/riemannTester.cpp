/*!
* \file riemannTester.cpp
* \author Robert 'Bob' Caddy (rvc@pitt.edu)
* \brief Test a Riemann solver
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
#include "MhdRiemannSolver.h"
#include "HlldRiemannSolver.h"
#include "Grid1D.h"
#include "mhdUtilities.h"

using std::cout;
using std::cin;
using std::endl;
using namespace mhdUtilities;

/*!
 * \brief Test a Riemann solver by giving it interface states and a
 * position/time argument
 *
 * \return int
 */
int main()
{
    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    // ===== Settings for Sod Shock Tube =======================================
    // Dai & Woodward Shock tube
    // double const coef                   = 1. / std::sqrt(4. * M_PI);
    // double const densityL               = 1.08;
    // std::vector<double> const velocityL = {1.2, 0.01, 0.5};
    // std::vector<double> const magneticL = {4. * coef, 3.6 * coef, 2.0 * coef};
    // double const pressureL              = 0.95;
    // double const densityR               = 1.0;
    // std::vector<double> const velocityR = {0.0, 0.0, 0.0};
    // std::vector<double> const magneticR = {4.0*coef, 4.0*coef, 2.0*coef};
    // double const pressureR              = 1.0;
    // double const tMax                   = 0.2;

    // Brio & Wu Shock tube
    double const densityL               = 1.0;
    double const pressureL              = 1.0;
    std::vector<double> const velocityL = {0.0, 0.0, 0.0};
    std::vector<double> const magneticL = {0.75, 1.0, 0.0};
    double const densityR               = 0.125;
    double const pressureR              = 0.1;
    std::vector<double> const velocityR = {0.0, 0.0, 0.0};
    std::vector<double> const magneticR = {0.75, -1.0, 0.0};
    double const tMax                   = 0.1;

    double const gamma                  = 5./3.;
    double const xSize                  = 1.0;
    size_t const numSamples             = 800;
    double const diaphPos               = 0.5 * xSize;
    double const deltaX                 = xSize / static_cast<double>(numSamples);
    // ===== End Settings for Sod Shock Tube  ==================================

    // ===== Grid & Solver ===================================================
    Grid1D grid(numSamples, 0, "pass", "../data/");

    std::unique_ptr<MhdRiemannSolver> riemannSolver = std::unique_ptr<MhdRiemannSolver>(new HlldRiemannSolver(gamma));
    // ===== End Arrays & Solver ===============================================

    // ===== Compute the Solution ==============================================
    for (size_t i = 0; i < numSamples; i++)
    {
        // Compute posOverT
        double position  = (static_cast<double>(i) - 0.5) * deltaX;
        double posOverT  = (position - diaphPos) / tMax;

        // Run solver
        double denFlux, eneFlux;  // Variables for fluxes that I'm not using
        std::vector<double> magFlux(3), momFlux(3);
        riemannSolver->riemannMain(densityL,
                                   velocityL,
                                   pressureL,
                                   magneticL,
                                   densityR,
                                   velocityR,
                                   pressureR,
                                   magneticR,
                                   denFlux,
                                   momFlux,
                                   magFlux,
                                   eneFlux,
                                   posOverT);

        // Get the state values
        grid.density[i]  = riemannSolver->getDensityState();
        for (size_t j = 0; j < 3; j++)
        {
            grid.momentum[i][j] = computeMomentum(riemannSolver->getVelocityState()[j],
                                                  riemannSolver->getDensityState());
        }
        grid.magnetic[i] = riemannSolver->getMagneticState();
        grid.energy[i]   = computeEnergy(riemannSolver->getPressureState(),
                                         riemannSolver->getDensityState(),
                                         riemannSolver->getVelocityState(),
                                         riemannSolver->getMagneticState(),
                                         gamma);
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