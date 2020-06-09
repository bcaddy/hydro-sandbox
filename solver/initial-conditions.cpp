#include <iostream>
#include <cmath>

#include "initial-conditions.h"

using std::cout;
using std::endl;

int setInitialConditions(std::vector<double> &arr,
                         const int &size,
                         const std::string kind)
/*
Set the initial conditions. Options are 'top-hat', 'vel-step', and 'vel-sine'
*/
{
    if (kind == "top-hat")
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
    else if (kind == "vel-step")
    {
        // Velocity Step function. u=1 for left half and u=2 for right half
        for (int i = 0; i < size; i++)
        {
            if (i <= (size / 2))
            {
                arr[i] = 1.;
            }
            else
            {
                arr[i] = 2.;
            }
        }
    }
    else if (kind == "vel-sine")
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
                arr[i] = sine[sineIndex];
                sineIndex++;
            }
            else
            {
                arr[i] = 1.;
            }
        }
    }
    else
    {
        cout << "The initial condition you chose is unavailble. Please check"
             << "for typos or choose an available initial condition" << endl;
        return 1;
    }

    return 1;
}