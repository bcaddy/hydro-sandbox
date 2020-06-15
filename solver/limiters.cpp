#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>
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
    double MCLimiter(const double &a0,
                     const double &a1,
                     const double &a2,
                     const double &deltax);
    double centeredDif(const double &a0,
                       const double &a2,
                       const double &deltax);
    double forwardDif(const double &a1,
                       const double &a2,
                       const double &deltax);
    double backwardDif(const double &a0,
                       const double &a1,
                       const double &deltax);
}

// =============================================================================
double SlopeLimiter(const double &a0,
                    const double &a1,
                    const double &a2,
                    const double &deltax,
                    const std::string &kind)
{
    if (kind == "minMod") {return minModLimiter(a0, a1, a2, deltax);}
    else if (kind == "MC") {return MCLimiter(a0, a1, a2, deltax);}
    else if (kind == "none-centered") {return centeredDif(a0, a2, deltax);}
    else if (kind == "none-forward") {return forwardDif(a1, a2, deltax);}
    else if (kind == "none-backward") {return backwardDif(a0, a1, deltax);}
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

// =============================================================================
double MCLimiter(const double &a0,
                     const double &a1,
                     const double &a2,
                     const double &deltax)
    /*
    Implementation of the Monotonized Central Difference (MC) Limiter
    */
{
    double xi = (a2 - a1) * (a1 - a0);

    if (xi > 0)
    {
        double centerDif, forwardDif, backwardDif, sign;
        centerDif   =     std::abs(a2 - a0) / (2 * deltax);
        forwardDif  = 2 * std::abs(a2 - a1) / (deltax);
        backwardDif = 2 * std::abs(a1 - a0) / (deltax);

        sign = ((a2-a0) < 0)? -1:1; // equivalent to sign(a2-a0)

        return sign * std::min({centerDif, forwardDif, backwardDif});
    }
    else
    {
        return 0.;
    }

}
// =============================================================================

// =============================================================================
// Forward, backward, and centered diffs
double centeredDif(const double &a0,
                   const double &a2,
                   const double &deltax)
{
    return (a2-a0)/deltax;
}
double forwardDif(const double &a1,
                   const double &a2,
                   const double &deltax)
{
    return (a2 - a1) / deltax;
}
double backwardDif(const double &a0,
                   const double &a1,
                   const double &deltax)
{
    return (a1 - a0) / deltax;
}
// =============================================================================
}