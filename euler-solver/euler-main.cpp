

#include <iostream>
#include <string>
#include <chrono>
using namespace std::chrono;
using std::cout;
using std::cin;
using std::endl;

#include "Grid1D.h"
// #include "Simulation1D"

int main()
{
    // Start clock
    auto start = high_resolution_clock::now();

    Grid1D temp(10, 2, "/Users/Bob/Desktop/PhD-Research/hydro-sandbox/data/");
    temp.velocity[6] = 3;
    temp.density[6] = 4;
    temp.siEnergy[6] = 2;
    cout << temp.ComputeTotalSpecEnergyElement(6);



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