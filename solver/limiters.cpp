#include <cmath>

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