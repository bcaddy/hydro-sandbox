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
        for (int i = numGhosts; i < (size - numGhosts); i++)
        {
            fileName << arr[i] << ",";
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
            sine[i] = std::sin(twoPi * static_cast<double>(i) / sizeSine);
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