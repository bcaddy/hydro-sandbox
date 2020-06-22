/*
================================================================================
 Euler Equation Solver
 Written by Robert Caddy.  Created on June 22, 2020

 Description:
     This program solves the Euler equations in 1D primarily through using the
     member functions of the Grid1D class. The algorithm used is a second order
     Godunov Method.

 Dependencies:
     Grid1D class
     primitives class
================================================================================
*/

#include <iostream>
#include <string>
#include <chrono>
using namespace std::chrono;
using std::cout;
using std::cin;
using std::endl;

#include "primitives.h"

int main()
{
    // Start timer
    auto start = high_resolution_clock::now();


    primitives temp(10);

    temp.velocity[3] = 13.;

    for (size_t i = 0; i < 10; i++)
    {
        cout << temp.velocity[i] << endl;
    }


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