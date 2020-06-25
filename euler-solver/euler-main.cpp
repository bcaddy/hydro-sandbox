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
using namespace std::chrono;
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
    auto start = high_resolution_clock::now();

    Simulation1D sim(3.,
                     .4,
                     10,
                     2,
                     "Sod",
                     "MC",
                     "here");

    cout << sim.physLen << endl;
    cout << sim.CFLnum << endl;
    cout << sim.deltaX << endl;
    cout << sim.limiterKind << endl;
    cout << sim.initialConditionsKind << endl;
    cout << sim.saveDir << endl;

    // Stop timer and print execution time. Time options are nanoseconds,
    // microseconds, milliseconds, seconds, minutes, hours. To pick one just
    // change `using FpTime` and the cout statement suitably.
    auto stop = high_resolution_clock::now();
    using FpTime = duration<float, seconds::period>;
    static_assert(treat_as_floating_point<FpTime::rep>::value,"Rep required to be floating point");
    auto duration = FpTime(stop - start);
    cout << "Time to execute: " << duration.count() << " seconds";
    return 0;
}