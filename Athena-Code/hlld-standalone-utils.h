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
    if (false)
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

    // =========================================================================
    // Dai & Woodward Shock Tube
    // =========================================================================
    if (false)
    {
        double const gamma = 5./3.;
        double const coef = 1. / (std::sqrt(4 * M_PI));
        double const Bx = 4. * coef;
        std::string rootName = "Dai & Woodward, ";
        Prim1DS leftICs                (1.08,    1.0,      0.0,       0.0,       0.0,      Bx,      3.6*coef, 2*coef);
        Prim1DS leftFastShockLeftSide  (1.09406, 0.970815, 1.176560,  0.021003,  0.506113, 1.12838, 1.105355, 0.614087);
        Prim1DS leftFastShockRightSide (1.40577, 1.494290, 0.693255,  0.210562,  0.611423, 1.12838, 1.457700, 0.809831);
        Prim1DS leftRotationLeftSide   (1.40086, 1.485660, 0.687774,  0.215124,  0.609161, 1.12838, 1.458735, 0.789960);
        Prim1DS leftRotationRightSide  (1.40119, 1.486570, 0.687504,  0.330268,  0.334140, 1.12838, 1.588975, 0.475782);
        Prim1DS leftSlowShockLeftSide  (1.40519, 1.493710, 0.685492,  0.326265,  0.333664, 1.12838, 1.575785, 0.472390);
        Prim1DS leftSlowShockRightSide (1.66488, 1.984720, 0.578545,  0.050746,  0.250260, 1.12838, 1.344490, 0.402407);
        Prim1DS contactLeftSide        (1.65220, 1.981250, 0.578296,  0.049683,  0.249962, 1.12838, 1.346155, 0.402868);
        Prim1DS contactRightSide       (1.49279, 1.981160, 0.578276,  0.049650,  0.249924, 1.12838, 1.346180, 0.402897);
        Prim1DS rightSlowShockLeftSide (1.48581, 1.956320, 0.573195,  0.035338,  0.245592, 1.12838, 1.370395, 0.410220);
        Prim1DS rightSlowShockRightSide(1.23813, 1.439000, 0.450361, -0.275532,  0.151746, 1.12838, 1.609775, 0.482762);
        Prim1DS rightRotationLeftSide  (1.23762, 1.437950, 0.450102, -0.274410,  0.145585, 1.12838, 1.606945, 0.493879);
        Prim1DS rightRotationRightSide (1.23747, 1.437350, 0.449993, -0.180766, -0.090238, 1.12838, 1.503855, 0.752090);
        Prim1DS rightFastShockLeftSide (1.22305, 1.409660, 0.424403, -0.171402, -0.085701, 1.12838, 1.447730, 0.723864);
        Prim1DS rightFastShockRightSide(1.00006, 1.000100, 0.000121, -0.000057, -0.000028, 1.12838, 1.128435, 0.564217);
        Prim1DS rightICs               (1.0,     0.2,      0.0,       0.0,       1.0,      Bx,      4*coef,   2*coef);

        // Compare ICs
        names.push_back(rootName + "initial conditions interface");
        leftPrimitive.push_back (leftICs);
        rightPrimitive.push_back(rightICs);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "initial conditions interface, reversed");
        leftPrimitive.push_back (rightICs);
        rightPrimitive.push_back(leftICs);
        gammaVector.push_back(gamma);

        // Results across each wave
        names.push_back(rootName + "Left Fast Shock");
        leftPrimitive.push_back (leftFastShockLeftSide);
        rightPrimitive.push_back(leftFastShockRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "left rotation wave");
        leftPrimitive.push_back (leftRotationLeftSide);
        rightPrimitive.push_back(leftRotationRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "left slow shock");
        leftPrimitive.push_back (leftSlowShockLeftSide);
        rightPrimitive.push_back(leftSlowShockRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "contact wave");
        leftPrimitive.push_back (contactLeftSide);
        rightPrimitive.push_back(contactRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "right slow shock");
        leftPrimitive.push_back (rightSlowShockLeftSide);
        rightPrimitive.push_back(rightSlowShockRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "right rotation wave");
        leftPrimitive.push_back (rightRotationLeftSide);
        rightPrimitive.push_back(rightRotationRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "right fast shock");
        leftPrimitive.push_back (rightFastShockLeftSide);
        rightPrimitive.push_back(rightFastShockRightSide);
        gammaVector.push_back(gamma);

        // Results across each wave from right to left
        rootName = rootName + "Reversed, ";

        names.push_back(rootName + "Left Fast Shock");
        leftPrimitive.push_back(leftFastShockRightSide);
        rightPrimitive.push_back (leftFastShockLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "left rotation wave");
        leftPrimitive.push_back(leftRotationRightSide);
        rightPrimitive.push_back (leftRotationLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "left slow shock");
        leftPrimitive.push_back(leftSlowShockRightSide);
        rightPrimitive.push_back (leftSlowShockLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "contact wave");
        leftPrimitive.push_back(contactRightSide);
        rightPrimitive.push_back (contactLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "right slow shock");
        leftPrimitive.push_back(rightSlowShockRightSide);
        rightPrimitive.push_back (rightSlowShockLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "right rotation wave");
        leftPrimitive.push_back(rightRotationRightSide);
        rightPrimitive.push_back (rightRotationLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "right fast shock");
        leftPrimitive.push_back(rightFastShockRightSide);
        rightPrimitive.push_back (rightFastShockLeftSide);
        gammaVector.push_back(gamma);

    }

    // =========================================================================
    // Ryu & Jones 4d Shock Tube
    // =========================================================================
    if (false)
    {
        double const gamma = 5./3.;
        double const Bx = 0.7;
        std::string rootName = "Ryu & Jones 4d, ";
        Prim1DS leftICs                   (1.0,      1.0,       0.0,      0.0,          0.0,          Bx, 0.0,           0.0);
        Prim1DS hydroRareLeftSide         (0.990414, 0.984076,  0.012415, 1.458910e-58, 6.294360e-59, Bx, 1.252355e-57,  5.366795e-58);
        Prim1DS hydroRareRightSide        (0.939477, 0.901182,  0.079800, 1.557120e-41, 7.505190e-42, Bx, 1.823624e-40,  8.712177e-41);
        Prim1DS switchOnSlowShockLeftSide (0.939863, 0.901820,  0.079142, 1.415730e-02, 7.134030e-03, Bx, 2.519650e-02,  1.290082e-02);
        Prim1DS switchOnSlowShockRightSide(0.651753, 0.490103,  0.322362, 8.070540e-01, 4.425110e-01, Bx, 6.598380e-01,  3.618000e-01);
        Prim1DS contactLeftSide           (0.648553, 0.489951,  0.322525, 8.072970e-01, 4.426950e-01, Bx, 6.599295e-01,  3.618910e-01);
        Prim1DS contactRightSide          (0.489933, 0.489980,  0.322518, 8.073090e-01, 4.426960e-01, Bx, 6.599195e-01,  3.618850e-01);
        Prim1DS slowShockLeftSide         (0.496478, 0.489823,  0.308418, 8.060830e-01, 4.420150e-01, Bx, 6.686695e-01,  3.666915e-01);
        Prim1DS slowShockRightSide        (0.298260, 0.198864, -0.016740, 2.372870e-01, 1.287780e-01, Bx, 8.662095e-01,  4.757390e-01);
        Prim1DS rotationLeftSide          (0.298001, 0.198448, -0.017358, 2.364790e-01, 1.278540e-01, Bx, 8.669425e-01,  4.750845e-01);
        Prim1DS rotationRightSide         (0.297673, 0.197421, -0.018657, 1.059540e-02, 9.996860e-01, Bx, 9.891580e-01,  1.024949e-04);
        Prim1DS fastRareLeftSide          (0.297504, 0.197234, -0.020018, 1.137420e-02, 1.000000e+00, Bx, 9.883860e-01, -4.981931e-17);
        Prim1DS fastRareRightSide         (0.299996, 0.199995, -0.000033, 1.855120e-05, 1.000000e+00, Bx, 9.999865e-01,  1.737190e-16);
        Prim1DS rightICs                  (0.3,      0.2,       0.0,      0.0,          1.0,          Bx, 1.0,           0.0);

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
        names.push_back(rootName + "hydro rarefaction");
        leftPrimitive.push_back (hydroRareLeftSide);
        rightPrimitive.push_back(hydroRareRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "switch on slow shock");
        leftPrimitive.push_back (switchOnSlowShockLeftSide);
        rightPrimitive.push_back(switchOnSlowShockRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "contact wave");
        leftPrimitive.push_back (contactLeftSide);
        rightPrimitive.push_back(contactRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "slow shock");
        leftPrimitive.push_back (slowShockLeftSide);
        rightPrimitive.push_back(slowShockRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "rotation wave");
        leftPrimitive.push_back (rotationLeftSide);
        rightPrimitive.push_back(rotationRightSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "fast rarefaction");
        leftPrimitive.push_back (fastRareLeftSide);
        rightPrimitive.push_back(fastRareRightSide);
        gammaVector.push_back(gamma);

        // Results across each wave from right to left
        rootName = rootName + "Reversed, ";

        names.push_back(rootName + "hydro rarefaction");
        leftPrimitive.push_back(hydroRareRightSide);
        rightPrimitive.push_back (hydroRareLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "switch on slow shock");
        leftPrimitive.push_back(switchOnSlowShockRightSide);
        rightPrimitive.push_back (switchOnSlowShockLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "contact wave");
        leftPrimitive.push_back(contactRightSide);
        rightPrimitive.push_back (contactLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "slow shock");
        leftPrimitive.push_back(slowShockRightSide);
        rightPrimitive.push_back (slowShockLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "rotation wave");
        leftPrimitive.push_back(rotationRightSide);
        rightPrimitive.push_back (rotationLeftSide);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "fast rarefaction");
        leftPrimitive.push_back(fastRareRightSide);
        rightPrimitive.push_back (fastRareLeftSide);
        gammaVector.push_back(gamma);
    }

    // =========================================================================
    // EFR
    // =========================================================================
    if (false)
    {
        double const gamma = 5./3.;
        double const V0 = 2.;
        double const Vy = 0.0;
        double const Vz = 0.0;
        double const Bx = 0.0;
        double const Bz = 0.0;
        std::string rootName = "EFR, , ";
        Prim1DS leftICs               (1.0,      0.45,     -V0,       Vy, Vz, Bx, 0.5,      Bz);
        Prim1DS leftRarefactionCenter (0.368580, 0.111253, -1.180830, Vy, Vz, Bx, 0.183044, Bz);
        Prim1DS leftVxTurnOver        (0.058814, 0.008819, -0.125475, Vy, Vz, Bx, 0.029215, Bz);
        Prim1DS midPoint              (0.034658, 0.006776,  0.000778, Vy, Vz, Bx, 0.017333, Bz);
        Prim1DS rightVxTurnOver       (0.062587, 0.009521,  0.152160, Vy, Vz, Bx, 0.031576, Bz);
        Prim1DS rightRarefactionCenter(0.316485, 0.089875,  1.073560, Vy, Vz, Bx, 0.159366, Bz);
        Prim1DS rightICs              (1.0,      0.45,      V0,       Vy, Vz, Bx, 0.5,      Bz);

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
        names.push_back(rootName + "left IC to left rarefaction center");
        leftPrimitive.push_back (leftICs);
        rightPrimitive.push_back(leftRarefactionCenter);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "left rarefaction center to left turnover");
        leftPrimitive.push_back (leftRarefactionCenter);
        rightPrimitive.push_back(leftVxTurnOver);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "left turnover to center");
        leftPrimitive.push_back (leftVxTurnOver);
        rightPrimitive.push_back(midPoint);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "center to right turnover");
        leftPrimitive.push_back (midPoint);
        rightPrimitive.push_back(rightVxTurnOver);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "right turnover to right rarefaction center");
        leftPrimitive.push_back (rightVxTurnOver);
        rightPrimitive.push_back(rightRarefactionCenter);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "right rarefaction center to right IC");
        leftPrimitive.push_back (rightRarefactionCenter);
        rightPrimitive.push_back(rightICs);
        gammaVector.push_back(gamma);

        // Results across each wave from right to left
        rootName = rootName + "Reversed, ";

        names.push_back(rootName + "left IC to left rarefaction center");
        leftPrimitive.push_back(leftRarefactionCenter);
        rightPrimitive.push_back (leftICs);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "left rarefaction center to left turnover");
        leftPrimitive.push_back(leftVxTurnOver);
        rightPrimitive.push_back (leftRarefactionCenter);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "left turnover to center");
        leftPrimitive.push_back(midPoint);
        rightPrimitive.push_back (leftVxTurnOver);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "center to right turnover");
        leftPrimitive.push_back(rightVxTurnOver);
        rightPrimitive.push_back (midPoint);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "right turnover to right rarefaction center");
        leftPrimitive.push_back(rightRarefactionCenter);
        rightPrimitive.push_back (rightVxTurnOver);
        gammaVector.push_back(gamma);

        names.push_back(rootName + "right rarefaction center to right IC");
        leftPrimitive.push_back(rightICs);
        rightPrimitive.push_back (rightRarefactionCenter);
        gammaVector.push_back(gamma);

    }

    // =========================================================================
    // Potential Divide by Zero Error Shock Tube
    // =========================================================================
    if (true)
    {
        // This is the case where:
        // Sm     = Vx_side
        // S_side = Vx_side +/- c_f_side
        // Bx^2 >= gamma * pressure_side
        // By_side = Bz_side = 0
        // Which simplifies (according to Athena solver) to
        // (rho_side * S_side * Sm - Bx^2) < (small_number * total_pressure_side * density_side * S_side * (S_side - S^star_side))
        double const gamma = 5./3.;

        double const rhoL      = 0.0;
        double const pressureL = 1.0;
        double const vxL       = 0.0;
        double const vyL       = 0.0;
        double const vzL       = 0.0;
        double const bxL       = 0.0;
        double const byL       = 0.0;
        double const bzL       = 0.0;

        double const rhoR      = 0.0;
        double const pressureR = 1.0;
        double const vxR       = 0.0;
        double const vyR       = 0.0;
        double const vzR       = 0.0;
        double const bxR       = 0.0;
        double const byR       = 0.0;
        double const bzR       = 0.0;

        std::string rootName = "divide by zero error, ";

        Prim1DS leftICs (rhoL, pressureL, vxL, vyL, vzL, bxL, byL, bzL);
        Prim1DS rightICs(rhoR, pressureR, vxR, vyR, vzR, bxR, byR, bzR);

        // Compare ICs
        names.push_back(rootName);
        leftPrimitive.push_back (leftICs);
        rightPrimitive.push_back(rightICs);
        gammaVector.push_back(gamma);
    }


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
        << " | " << std::setw(15) << "Vel/Momentum x"    << " | " << spacer << conservedLeft.Mx << " | " << spacer << conservedRight.Mx << " | " << spacer << primLeft.Vx << " | " << spacer << primRight.Vx << " | "  << spacer << fluxes.Mx << " | " << std::endl
        << " | " << std::setw(15) << "Vel/Momentum y"    << " | " << spacer << conservedLeft.My << " | " << spacer << conservedRight.My << " | " << spacer << primLeft.Vy << " | " << spacer << primRight.Vy << " | "  << spacer << fluxes.My << " | " << std::endl
        << " | " << std::setw(15) << "Vel/Momentum z"    << " | " << spacer << conservedLeft.Mz << " | " << spacer << conservedRight.Mz << " | " << spacer << primLeft.Vz << " | " << spacer << primRight.Vz << " | "  << spacer << fluxes.Mz << " | " << std::endl
        << " | " << std::setw(15) << "Magnetic x"        << " | " << spacer << conservedLeft.Bx << " | " << spacer << conservedRight.Bx << " | " << spacer << primLeft.Bx << " | " << spacer << primRight.Bx << " | "  << spacer << fluxes.Bx << " | " << std::endl
        << " | " << std::setw(15) << "Magnetic y"        << " | " << spacer << conservedLeft.By << " | " << spacer << conservedRight.By << " | " << spacer << primLeft.By << " | " << spacer << primRight.By << " | "  << spacer << fluxes.By << " | " << std::endl
        << " | " << std::setw(15) << "Magnetic x"        << " | " << spacer << conservedLeft.Bz << " | " << spacer << conservedRight.Bz << " | " << spacer << primLeft.Bz << " | " << spacer << primRight.Bz << " | "  << spacer << fluxes.Bz << " | " << std::endl
        << " -------------------------------------------------------------------------------------------------------------------------------------------" << std::endl
    ;
}
// =============================================================================
