#include <iostream>
#include <cmath>

#include "initial-conditions.h"

using std::cout;
using std::endl;

namespace
{
    // This namespace effectively makes the listed function private to this
    // file and they cannot be accessed elsewhere.
    void TopHat(std::vector<double> &arr, const int &size);
    void VelStep(std::vector<double> &arr, const int &size, const double &sign);
    void VelSine(std::vector<double> &arr, const int &size, const double &sign);
}

// =============================================================================
int setInitialConditions(std::vector<double> &arr,
                         const int &size,
                         const std::string kind)
/*
Set the initial conditions. Options are 'top-hat', 'vel-step', and 'vel-sine'
*/
{
    if (kind == "top-hat") {TopHat(arr, size);}

    else if (kind == "vel-step-pos") {VelStep(arr, size, 1.);}
    else if (kind == "vel-step-neg") {VelStep(arr, size, -1.);}
    else if (kind == "vel-sine-pos") {VelSine(arr, size, 1.);}
    else if (kind == "vel-sine-neg") {VelSine(arr, size, -1.);}
    else
    {
        cout << "The initial condition you chose is unavailble. Please check"
             << "for typos or choose an available initial condition" << endl;
        return 1;
    }

    return 1;
}
// =============================================================================

namespace
{ // Same namespace from earlier
// =============================================================================
void TopHat(std::vector<double> &arr,
            const int &size)
{
    // Top hat
    for (int i = 0; i < size; i++)
    {
        if (i <= (size / 10))
        {
            arr[i] = 0.;
        }
        else if (((size / 10) < i) && (i <= (2 * size / 10)))
        {
            arr[i] = 1.;
        }
        else
        {
            arr[i] = 0.;
        }
    }
}
// =============================================================================

// =============================================================================
void VelStep(std::vector<double> &arr,
             const int &size,
             const double &sign)
{
    // Velocity Step function. u=1 for left half and u=2 for right half
    for (int i = 0; i < size; i++)
    {
        if (i <= (size / 2))
        {
            arr[i] = sign;
        }
        else
        {
            arr[i] = 2. * sign;
        }
    }
}
// =============================================================================

// =============================================================================
void VelSine(std::vector<double> &arr,
             const int &size,
             const double &sign)
{
    // Velocity sine function. One full period in the middle half of the
    // simulation
    std::vector<double> sine(size / 2);
    const double twoPi = 2 * 3.14159265;
    const double sizeSine = static_cast<double>(size / 2);

    for (int i = 0; i < (size / 2); i++)
    {
        sine[i] = 1. + 0.5 * std::sin(twoPi * static_cast<double>(i) / sizeSine);
    }

    int sineIndex = 0;
    for (int i = 0; i < size; i++)
    {
        if ((i > (size / 4)) && (i < (3 * size / 4)))
        {
            arr[i] = sign * sine[sineIndex];
            sineIndex++;
        }
        else
        {
            arr[i] = sign;
        }
    }
}
// =============================================================================
}