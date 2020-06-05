#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using std::cout;
using std::endl;

#include "helper.h"

void saveArray(const std::vector<double> &arr,
               std::ofstream &fileName,
               const int &numGhosts)
{
    int size = arr.size();

    if (fileName.is_open())
    {
        fileName << arr[numGhosts];
        for (int i = numGhosts+1; i < (size - numGhosts); i++)
        {
            fileName << "," << arr[i];
        }
        fileName << endl;
    }
    else
    {
        cout << "File not open";
    }
}

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
        std::vector<double> sine(size/2);
        const double twoPi=2*3.14159265;
        const double sizeSine = static_cast<double>(size / 2);

        for (int i = 0; i < (size / 2); i++)
        {
            sine[i] = 1. + 0.5*std::sin(twoPi * static_cast<double>(i) / sizeSine);
        }

        int sineIndex = 0;
        for (int i = 0; i < size; i++)
        {
            if ( (i > (size/4)) && (i < (3*size/4)) )
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
        cout << "The initial condition you chose is unavailble. Please check" <<
                 "for typos or choose an available initial condition" << endl;
        return 1;
    }
    
    return 1;
}

double minModLimiter(const double &a0,
                     const double &a1,
                     const double &a2,
                     const double &deltax)
{
    // Declare variables
    double outValue;
    double leftDerive, rightDerive, leftRight;

    // Compute the derivatives
    leftDerive = (a1 - a0) / deltax;
    rightDerive = (a2 - a1) / deltax;
    leftRight = leftDerive * rightDerive;

    // Choose what value to output
    if ((std::abs(leftDerive) < std::abs(rightDerive)) && (leftRight > 0.))
    {
        outValue = leftDerive;
    }
    else if ((std::abs(leftDerive) > std::abs(rightDerive)) && (leftRight > 0.))
    {
        outValue = rightDerive;
    }
    else
    {
        outValue = 0.;
    }
    
    return outValue;
}

double VelInterface(const double &a,
                    const double &b,
                    const double &c,
                    const double &d,
                    const double &deltax,
                    const double &deltat)
/* 
a, b, c, and d are the four elements of the array that are needed to compute
the interface where a has the lowest index and d has the higest. 

For the i-1/2 interface  | For the i+1/2 interface
a = uVel[i-2];           | a = uVel[i-1]; 
b = uVel[i-1];           | b = uVel[i];
c = uVel[i];             | c = uVel[i+1];  
d = uVel[i+1];           | d = uVel[i+2];
*/
{
    double uL, uR, Interface;

    uL = b +
         0.5 * deltax * (1 - b * deltat / deltax) *
         minModLimiter(a, b, c, deltax);

    uR = c -
         0.5 * deltax * (1 + c * deltat / deltax) *
         minModLimiter(b, c, d, deltax);

    if (uL > uR)
    {
        double S = 0.5 * (uL + uR);
        if (S > 0.)
        {
            Interface = uL;
        }
        else if (S < 0.)
        {
            Interface = uR;
        }
        else // (S == 0)
        {
            Interface = 0.;
        }
    }
    // else if (uL < uR)
    // {
    //     double S = 0.5 * (uL + uR);
    //     if (S > 0.)
    //     {
    //         Interface = uR;
    //     }
    //     else if (S < 0.)
    //     {
    //         Interface = uL;
    //     }
    //     else // (S == 0)
    //     {
    //         Interface = 0.;
    //     }
    // }
    else // ie the two velocities are the same
    {
        if (uL > 0.)
        {
            Interface = uL;
        }
        else if (uR > 0.)
        {
            Interface = uR;
        }
        else // (S == 0)
        {
            Interface = 0.;
        }
    }

    return Interface;
}
