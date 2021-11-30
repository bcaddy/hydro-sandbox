#include <limits>
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <tuple>

#pragma once

// =============================================================================
#define SMALL_NUMBER 1e-8
#define TINY_NUMBER 1e-20
double SQR(double const &A){return A*A;};
// =============================================================================

// =============================================================================
struct Cons1DS
{
    double d;
    double E;
    double Mx;
    double My;
    double Mz;
    double Bx;
    double By;
    double Bz;
    Cons1DS()
        :
        d (std::numeric_limits<double>::quiet_NaN()),
        E (std::numeric_limits<double>::quiet_NaN()),
        Mx(std::numeric_limits<double>::quiet_NaN()),
        My(std::numeric_limits<double>::quiet_NaN()),
        Mz(std::numeric_limits<double>::quiet_NaN()),
        Bx(std::numeric_limits<double>::quiet_NaN()),
        By(std::numeric_limits<double>::quiet_NaN()),
        Bz(std::numeric_limits<double>::quiet_NaN())
    {};
    Cons1DS(double const &dInput,
            double const &EInput,
            double const &MxInput,
            double const &MyInput,
            double const &MzInput,
            double const &BxInput,
            double const &ByInput,
            double const &BzInput)
        :
        d (dInput),
        E (EInput),
        Mx(MxInput),
        My(MyInput),
        Mz(MzInput),
        Bx(BxInput),
        By(ByInput),
        Bz(BzInput)
    {};
};
// =============================================================================

// =============================================================================
struct Prim1DS
{
    double d;
    double P;
    double Vx;
    double Vy;
    double Vz;
    double Bx;
    double By;
    double Bz;
    Prim1DS()
        :
        d (std::numeric_limits<double>::quiet_NaN()),
        P (std::numeric_limits<double>::quiet_NaN()),
        Vx(std::numeric_limits<double>::quiet_NaN()),
        Vy(std::numeric_limits<double>::quiet_NaN()),
        Vz(std::numeric_limits<double>::quiet_NaN()),
        Bx(std::numeric_limits<double>::quiet_NaN()),
        By(std::numeric_limits<double>::quiet_NaN()),
        Bz(std::numeric_limits<double>::quiet_NaN())
    {};
    Prim1DS(double const &dInput,
            double const &PInput,
            double const &VxInput,
            double const &VyInput,
            double const &VzInput,
            double const &BxInput,
            double const &ByInput,
            double const &BzInput)
        :
        d (dInput),
        P (PInput),
        Vx(VxInput),
        Vy(VyInput),
        Vz(VzInput),
        Bx(BxInput),
        By(ByInput),
        Bz(BzInput)
    {};
};
// =============================================================================

// =============================================================================
Prim1DS Cons1D_to_Prim1D(const Cons1DS &pU,
                         const double &pBx,
                         double const &Gamma)
{
    Prim1DS Prim1D;

    double di = 1.0/pU.d;

    Prim1D.d  = pU.d;
    Prim1D.Vx = pU.Mx*di;
    Prim1D.Vy = pU.My*di;
    Prim1D.Vz = pU.Mz*di;

    Prim1D.P = pU.E - 0.5*(SQR(pU.Mx)+SQR(pU.My)+SQR(pU.Mz))*di;
    Prim1D.P -= 0.5*(SQR(pBx) + SQR(pU.By) + SQR(pU.Bz));
    Prim1D.P *= Gamma;
    Prim1D.P = std::max(Prim1D.P,TINY_NUMBER);

    Prim1D.By = pU.By;
    Prim1D.Bz = pU.Bz;

    return Prim1D;
}
// =============================================================================

// =============================================================================
Cons1DS Prim1D_to_Cons1D(const Prim1DS pW,
                         double const &Gamma)
{
  Cons1DS Cons1D;

  Cons1D.d  = pW.d;
  Cons1D.Mx = pW.d*pW.Vx;
  Cons1D.My = pW.d*pW.Vy;
  Cons1D.Mz = pW.d*pW.Vz;

  Cons1D.E = pW.P/Gamma + 0.5 * pW.d * (SQR(pW.Vx) + SQR(pW.Vy) + SQR(pW.Vz));
  Cons1D.E += 0.5*(SQR(pW.Bx) + SQR(pW.By) + SQR(pW.Bz));

  Cons1D.By = pW.By;
  Cons1D.Bz = pW.Bz;

  return Cons1D;
}
// =============================================================================

// =============================================================================
std::tuple<std::vector<std::string>,
           std::vector<Prim1DS>,
           std::vector<Prim1DS>,
           std::vector<double>>  generateTestCases()
{
    std::vector<std::string> names;
    std::vector<Prim1DS> leftPrimitive, rightPrimitive;
    std::vector<double> gammaVector;

    // =========================================================================
    // Brio & Wu Shock Tube
    // =========================================================================
    {
        // Setup the constants
        double const gamma = 2.;
        double const Vz = 0.0;
        double const Bx = 0.75;
        double const Bz = 0.0;
        std::string rootName = "Brio & Wu, ";
        Prim1DS leftICs               (1.0,      1.0     ,  0.0,       0.0,      Vz, Bx,  1.0     , Bz);
        Prim1DS leftFastRareLeftSide  (0.978576, 0.957621,  0.038603, -0.011074, Vz, Bx,  0.970288, Bz);
        Prim1DS leftFastRareRightSide (0.671655, 0.451115,  0.647082, -0.238291, Vz, Bx,  0.578240, Bz);
        Prim1DS compoundLeftSide      (0.814306, 0.706578,  0.506792, -0.911794, Vz, Bx, -0.108819, Bz);
        Prim1DS compoundPeak          (0.765841, 0.624742,  0.523701, -1.383720, Vz, Bx, -0.400787, Bz);
        Prim1DS compoundRightSide     (0.695211, 0.515237,  0.601089, -1.583720, Vz, Bx, -0.537027, Bz);
        Prim1DS contactLeftSide       (0.680453, 0.515856,  0.598922, -1.584490, Vz, Bx, -0.533616, Bz);
        Prim1DS contactRightSide      (0.231160, 0.516212,  0.599261, -1.584820, Vz, Bx, -0.533327, Bz);
        Prim1DS slowShockLeftSide     (0.153125, 0.191168,  0.086170, -0.683303, Vz, Bx, -0.850815, Bz);
        Prim1DS slowShockRightSide    (0.117046, 0.087684, -0.238196, -0.165561, Vz, Bx, -0.903407, Bz);
        Prim1DS rightFastRareLeftSide (0.117358, 0.088148, -0.228756, -0.158845, Vz, Bx, -0.908335, Bz);
        Prim1DS rightFastRareRightSide(0.124894, 0.099830, -0.003132, -0.002074, Vz, Bx, -0.999018, Bz);
        Prim1DS rightICs              (0.128,    0.1,       0.0,       0.0,      Vz, Bx, -1.0,      Bz);

        // Compare ICs
        names.push_back(rootName + "initial conditions interface");
        leftPrimitive.push_back (leftICs);
        rightPrimitive.push_back(rightICs);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "initial conditions interface, reversed");
        leftPrimitive.push_back (rightICs);
        rightPrimitive.push_back(leftICs);
        gammaVector.push_back(gamma);

        // Results across each wave from left to right
        names.push_back(rootName + "Left Fast Rarefaction");
        leftPrimitive.push_back (leftFastRareLeftSide);
        rightPrimitive.push_back(leftFastRareRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "Compound, left to right side");
        leftPrimitive.push_back (compoundLeftSide);
        rightPrimitive.push_back(compoundRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "Compound, left to peak");
        leftPrimitive.push_back (compoundLeftSide);
        rightPrimitive.push_back(compoundPeak);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "Compound, peak to right");
        leftPrimitive.push_back (compoundPeak);
        rightPrimitive.push_back(compoundRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "Slow shock");
        leftPrimitive.push_back (slowShockLeftSide);
        rightPrimitive.push_back(slowShockRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "Right Fast Rarefaction");
        leftPrimitive.push_back (rightFastRareLeftSide);
        rightPrimitive.push_back(rightFastRareRightSide);
        gammaVector.push_back(gamma);

        // Results across each wave from right to left
        rootName = rootName + "Reversed, ";

        names.push_back(rootName + "Left Fast Rarefaction");
        leftPrimitive.push_back(leftFastRareRightSide);
        rightPrimitive.push_back (leftFastRareLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "Compound, left to right side");
        leftPrimitive.push_back(compoundRightSide);
        rightPrimitive.push_back (compoundLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "Compound, left to peak");
        leftPrimitive.push_back(compoundPeak);
        rightPrimitive.push_back (compoundLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "Compound, peak to right");
        leftPrimitive.push_back(compoundRightSide);
        rightPrimitive.push_back (compoundPeak);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "Slow shock");
        leftPrimitive.push_back(slowShockRightSide);
        rightPrimitive.push_back (slowShockLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "Right Fast Rarefaction");
        leftPrimitive.push_back(rightFastRareRightSide);
        rightPrimitive.push_back (rightFastRareLeftSide);
        gammaVector.push_back(gamma);
    }




        // names.push_back(rootName + "");
        // leftPrimitive.push_back (LeftSide);
        // rightPrimitive.push_back(RightSide);
        // gammaVector.push_back(gamma);










    // Check that everything is the same length
    if ( not ((names.size() == leftPrimitive.size())
               and (names.size() == rightPrimitive.size())
               and (names.size() == gammaVector.size())
             )
       )
    {
        throw std::invalid_argument("Not all vectors are the same size");
    }

    return {names, leftPrimitive, rightPrimitive, gammaVector};
}
// =============================================================================

// =============================================================================
void printResults(Cons1DS const &conservedLeft,
                  Cons1DS const &conservedRight,
                  Prim1DS const &primLeft,
                  Prim1DS const &primRight,
                  Cons1DS const &fluxes,
                  std::string const &name,
                  std::string const &state,
                  double const &gamma)
{
    int maxWidth = std::numeric_limits<double>::max_digits10;
    std::cout.precision(maxWidth);
    auto spacer = std::setw(maxWidth+4);

    std::cout << std::endl;
    std::cout << "State = " << state << std::endl;
    std::cout << "Gamma = " << gamma << std::endl;
    std::cout
        << "Test Name: " << name << std::endl
        << " -------------------------------------------------------------------------------------------------------------------------------------------" << std::endl
        << " | " << std::setw(15) << "Field"             << " | " << spacer << "Conserved Left" << " | " << spacer << "Conserved Right" << " | " << spacer << "Primitive Left" << " | " << spacer << "Primitive Right" << " | " << spacer << "Conserved Fluxes"   << " | " << std::endl
        << " |-----------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|" << std::endl
        << " | " << std::setw(15) << "Density"           << " | " << spacer << conservedLeft.d  << " | " << spacer << conservedRight.d  << " | " << spacer << primLeft.d  << " | " << spacer << primRight.d  << " | "  << spacer << fluxes.d  << " | " << std::endl
        << " | " << std::setw(15) << "Energy/Pressure"   << " | " << spacer << conservedLeft.E  << " | " << spacer << conservedRight.E  << " | " << spacer << primLeft.P  << " | " << spacer << primRight.P  << " | "  << spacer << fluxes.E  << " | " << std::endl
        << " | " << std::setw(15) << "Momentum"          << " | " << spacer << conservedLeft.Mx << " | " << spacer << conservedRight.Mx << " | " << spacer << primLeft.Vx << " | " << spacer << primRight.Vx << " | "  << spacer << fluxes.Mx << " | " << std::endl
        << " | " << std::setw(15) << "Momentum"          << " | " << spacer << conservedLeft.My << " | " << spacer << conservedRight.My << " | " << spacer << primLeft.Vy << " | " << spacer << primRight.Vy << " | "  << spacer << fluxes.My << " | " << std::endl
        << " | " << std::setw(15) << "Momentum"          << " | " << spacer << conservedLeft.Mz << " | " << spacer << conservedRight.Mz << " | " << spacer << primLeft.Vz << " | " << spacer << primRight.Vz << " | "  << spacer << fluxes.Mz << " | " << std::endl
        << " | " << std::setw(15) << "Magnetic"          << " | " << spacer << conservedLeft.Bx << " | " << spacer << conservedRight.Bx << " | " << spacer << primLeft.Bx << " | " << spacer << primRight.Bx << " | "  << spacer << fluxes.Bx << " | " << std::endl
        << " | " << std::setw(15) << "Magnetic"          << " | " << spacer << conservedLeft.By << " | " << spacer << conservedRight.By << " | " << spacer << primLeft.By << " | " << spacer << primRight.By << " | "  << spacer << fluxes.By << " | " << std::endl
        << " | " << std::setw(15) << "Magnetic"          << " | " << spacer << conservedLeft.Bz << " | " << spacer << conservedRight.Bz << " | " << spacer << primLeft.Bz << " | " << spacer << primRight.Bz << " | "  << spacer << fluxes.Bz << " | " << std::endl
        << " -------------------------------------------------------------------------------------------------------------------------------------------" << std::endl
    ;
}
// =============================================================================
