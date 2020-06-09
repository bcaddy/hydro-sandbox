#include "limiters.h"

// =============================================================================
double VelInterface(const double &a,
                    const double &b,
                    const double &c,
                    const double &d,
                    const double &deltax,
                    const double &deltat,
                    const std::string &LimType)
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
             SlopeLimiter(a, b, c, deltax, LimType);

    uR = c -
         0.5 * deltax * (1 + c * deltat / deltax) *
             SlopeLimiter(b, c, d, deltax, LimType);

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
// =============================================================================

// =============================================================================
double AdvectionInterface(const double &a,
                          const double &b,
                          const double &c,
                          const double &d,
                          const double &vel,
                          const double &deltax,
                          const double &deltat,
                          const std::string &LimType)
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
    if (vel >= 0.)
    {
        // Compute the left side of the interface
        return b + 0.5 * deltax * (1 - b * deltat / deltax) *
               SlopeLimiter(a, b, c, deltax, LimType);
    }
    else // i.e. (vel < 0.)
    {
        // Compute the right side of the interface
        return c - 0.5 * deltax * (1 + c * deltat / deltax) *
               SlopeLimiter(b, c, d, deltax, LimType);
    };
}
// =============================================================================