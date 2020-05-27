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