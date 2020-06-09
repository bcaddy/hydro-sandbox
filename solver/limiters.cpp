#include <cmath>
#include <string>
#include <iostream>
using std::cout;
using std::endl;

namespace
{
    // This namespace effectively makes the listed function private to this
    // file and they cannot be accessed elsewhere.
    double minModLimiter(const double &a0,
                         const double &a1,
                         const double &a2,
                         const double &deltax);
}

// =============================================================================
double SlopeLimiter(const double &a0,
                    const double &a1,
                    const double &a2,
                    const double &deltax,
                    const std::string &kind)
{
    if (kind == "minMod")
    {
        return minModLimiter(a0, a1, a2, deltax);
    }
    // else if (kind == "MC")
    // {
        
    // }
    else
    {
        cout << "The limiter you chose is unavailble. Please check"
             << "for typos or choose an available initial condition" << endl;
        return 0.0;
    }
}
// =============================================================================

namespace
{
// =============================================================================
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
// =============================================================================
}