#include <cmath>
// #include <stdio.h>
// #include <stdlib.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <stdexcept>

// #include "../defs.h"
// #include "../athena.h"
// #include "../globals.h"
// #include "prototypes.h"
// #include "../prototypes.h"

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
                         const double  pBx,
                         double const &Gamma)
{
  Cons1DS Cons1D;

  Cons1D.d  = pW.d;
  Cons1D.Mx = pW.d*pW.Vx;
  Cons1D.My = pW.d*pW.Vy;
  Cons1D.Mz = pW.d*pW.Vz;

  Cons1D.E = pW.P/Gamma + 0.5 * pW.d * (SQR(pW.Vx) + SQR(pW.Vy) + SQR(pW.Vz));
  Cons1D.E += 0.5*(SQR(pBx) + SQR(pW.By) + SQR(pW.Bz));

  Cons1D.By = pW.By;
  Cons1D.Bz = pW.Bz;

  return Cons1D;
}
// =============================================================================

// =============================================================================
void printState(std::string const &state)
{
    std::cout << std::endl << "State = " << state << std::endl;
}
// =============================================================================

// =============================================================================
void printResults(Cons1DS const &conservedLeft,
                  Cons1DS const &conservedRight,
                  Cons1DS const &fluxes,
                  std::string const &name)
{
    int maxWidth = std::numeric_limits<double>::max_digits10;
    std::cout.precision(maxWidth);
    auto spacer = std::setw(maxWidth+4);

    std::cout
        << "Test Name: " << name << std::endl
        << " ------------------------------------------------------------------------------------" << std::endl
        << " | " << std::setw(8) << "Field"    << " | " << spacer << "Conserved Left" << " | " << spacer << "Conserved Right" << " | " << spacer << "Fluxes"   << " | " << std::endl
        << " |----------|-----------------------|-----------------------|-----------------------|" << std::endl
        << " | " << std::setw(8) << "Density"  << " | " << spacer << conservedLeft.d  << " | " << spacer << conservedRight.d  << " | "  << spacer << fluxes.d  << " | " << std::endl
        << " | " << std::setw(8) << "Energy"   << " | " << spacer << conservedLeft.E  << " | " << spacer << conservedRight.E  << " | "  << spacer << fluxes.E  << " | " << std::endl
        << " | " << std::setw(8) << "Momentum" << " | " << spacer << conservedLeft.Mx << " | " << spacer << conservedRight.Mx << " | "  << spacer << fluxes.Mx << " | " << std::endl
        << " | " << std::setw(8) << "Momentum" << " | " << spacer << conservedLeft.My << " | " << spacer << conservedRight.My << " | "  << spacer << fluxes.My << " | " << std::endl
        << " | " << std::setw(8) << "Momentum" << " | " << spacer << conservedLeft.Mz << " | " << spacer << conservedRight.Mz << " | "  << spacer << fluxes.Mz << " | " << std::endl
        << " | " << std::setw(8) << "Magnetic" << " | " << spacer << conservedLeft.Bx << " | " << spacer << conservedRight.Bx << " | "  << spacer << fluxes.Bx << " | " << std::endl
        << " | " << std::setw(8) << "Magnetic" << " | " << spacer << conservedLeft.By << " | " << spacer << conservedRight.By << " | "  << spacer << fluxes.By << " | " << std::endl
        << " | " << std::setw(8) << "Magnetic" << " | " << spacer << conservedLeft.Bz << " | " << spacer << conservedRight.Bz << " | "  << spacer << fluxes.Bz << " | " << std::endl
        << " ------------------------------------------------------------------------------------" << std::endl
    ;
}
// =============================================================================

// =============================================================================
void fluxes(const Cons1DS Ul,    // Left conserved state
            const Cons1DS Ur,    // Right conserved state
            const Prim1DS Wl,    // Left primitive state
            const Prim1DS Wr,    // Right primitive state
            const double Bxi,    // Mag field normal to interface
            Cons1DS &pFlux,      // Output flux
            double const &Gamma) // Gamma
{
    Cons1DS Ulst,Uldst,Urdst,Urst;       /* Conserved variable for all states */
    Prim1DS Wlst,Wrst;                   /* Primitive variables for all states */
    Cons1DS Fl,Fr;                       /* Fluxes for left & right states */
    double spd[5];                        /* signal speeds, left to right */
/*  double maxspd; */
    double sdl,sdr,sdml,sdmr;             /* S_i-u_i, S_i-S_M (i=L or R) */
    double pbl,pbr;                       /* Magnetic pressures */
    double cfl,cfr,cfmax;                 /* Cf (left & right), max(cfl,cfr) */
    double gpl,gpr,gpbl,gpbr;             /* gamma*P, gamma*P + B */
    double sqrtdl,sqrtdr;                 /* sqrt of the L* & R* densities */
    double invsumd;                       /* 1/(sqrtdl + sqrtdr) */
    double ptl,ptr,ptst;                  /* total pressures */
    double vbstl,vbstr;                   /* v_i* dot B_i* for i=L or R */
    double Bxsig;                         /* sign(Bx) = 1 for Bx>0, -1 for Bx<0 */
    double Bxsq;                          /* Bx^2 */
    double tmp;                      /* Temp variable for repeated calculations */
#if (NSCALARS > 0)
    int n;
#endif



/*--- Step 1. ------------------------------------------------------------------
 * Convert left- and right- states in conserved to primitive variables.
 */

/*
    pbl = Cons1D_to_Prim1D(&Ul,&Wl,&Bxi);
    pbr = Cons1D_to_Prim1D(&Ur,&Wr,&Bxi);
*/

/*--- Step 2. ------------------------------------------------------------------
 * Compute left & right wave speeds according to Miyoshi & Kusano, eqn. (67)
 */

    Bxsq = Bxi*Bxi;
    pbl = 0.5*(Bxsq + SQR(Wl.By) + SQR(Wl.Bz));
    pbr = 0.5*(Bxsq + SQR(Wr.By) + SQR(Wr.Bz));
    gpl  = Gamma * Wl.P;
    gpr  = Gamma * Wr.P;
    gpbl = gpl + 2.0*pbl;
    gpbr = gpr + 2.0*pbr;

    cfl = sqrt((gpbl + sqrt(SQR(gpbl)-4.0*gpl*Bxsq))/(2.0*Wl.d));
    cfr = sqrt((gpbr + sqrt(SQR(gpbr)-4.0*gpr*Bxsq))/(2.0*Wr.d));
    cfmax = std::max(cfl,cfr);

    if(Wl.Vx <= Wr.Vx) {
        spd[0] = Wl.Vx - cfmax;
        spd[4] = Wr.Vx + cfmax;
    }
    else {
        spd[0] = Wr.Vx - cfmax;
        spd[4] = Wl.Vx + cfmax;
    }

/*  maxspd = MAX(fabs(spd[0]),fabs(spd[4])); */

/*--- Step 3. ------------------------------------------------------------------
 * Compute L/R fluxes
 */

    /* total pressure */
    ptl = Wl.P + pbl;
    ptr = Wr.P + pbr;

    Fl.d  = Ul.Mx;
    Fl.Mx = Ul.Mx*Wl.Vx + ptl - Bxsq;
    Fl.My = Ul.d*Wl.Vx*Wl.Vy - Bxi*Ul.By;
    Fl.Mz = Ul.d*Wl.Vx*Wl.Vz - Bxi*Ul.Bz;
    Fl.E  = Wl.Vx*(Ul.E + ptl - Bxsq) - Bxi*(Wl.Vy*Ul.By + Wl.Vz*Ul.Bz);
    Fl.By = Ul.By*Wl.Vx - Bxi*Wl.Vy;
    Fl.Bz = Ul.Bz*Wl.Vx - Bxi*Wl.Vz;

    Fr.d  = Ur.Mx;
    Fr.Mx = Ur.Mx*Wr.Vx + ptr - Bxsq;
    Fr.My = Ur.d*Wr.Vx*Wr.Vy - Bxi*Ur.By;
    Fr.Mz = Ur.d*Wr.Vx*Wr.Vz - Bxi*Ur.Bz;
    Fr.E  = Wr.Vx*(Ur.E + ptr - Bxsq) - Bxi*(Wr.Vy*Ur.By + Wr.Vz*Ur.Bz);
    Fr.By = Ur.By*Wr.Vx - Bxi*Wr.Vy;
    Fr.Bz = Ur.Bz*Wr.Vx - Bxi*Wr.Vz;

#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) {
        Fl.s[n] = Fl.d*Wl.r[n];
        Fr.s[n] = Fr.d*Wr.r[n];
    }
#endif

/*--- Step 4. ------------------------------------------------------------------
 * Return upwind flux if flow is supersonic
 */

    if(spd[0] >= 0.0){
        pFlux = Fl;
        printState("L");
        return;
    }

    if(spd[4] <= 0.0){
        pFlux = Fr;
        printState("R");
        return;
    }

/*--- Step 5. ------------------------------------------------------------------
 * Compute middle and Alfven wave speeds
 */

    sdl = spd[0] - Wl.Vx;
    sdr = spd[4] - Wr.Vx;

    /* S_M: eqn (38) of Miyoshi & Kusano */
    spd[2] = (sdr*Wr.d*Wr.Vx - sdl*Wl.d*Wl.Vx - ptr + ptl) /
                     (sdr*Wr.d-sdl*Wl.d);

    sdml   = spd[0] - spd[2];
    sdmr   = spd[4] - spd[2];
    /* eqn (43) of Miyoshi & Kusano */
    Ulst.d = Ul.d * sdl/sdml;
    Urst.d = Ur.d * sdr/sdmr;
    sqrtdl = sqrt(Ulst.d);
    sqrtdr = sqrt(Urst.d);

    /* eqn (51) of Miyoshi & Kusano */
    spd[1] = spd[2] - fabs(Bxi)/sqrtdl;
    spd[3] = spd[2] + fabs(Bxi)/sqrtdr;

/*--- Step 6. ------------------------------------------------------------------
 * Compute intermediate states
 */

    ptst = ptl + Ul.d*sdl*(sdl-sdml);

/* Ul* */
    /* eqn (39) of M&K */
    Ulst.Mx = Ulst.d * spd[2];
//   if((fabs(spd[2]/Wl.Vx-1.0)<SMALL_NUMBER) ||
//      (fabs(spd[2])/fabs(spd[0]) <= SMALL_NUMBER &&
//       fabs(Wl.Vx)/fabs(spd[0]) <= SMALL_NUMBER)) {
//     Ulst.My = Ulst.d * Wl.Vy;
//     Ulst.Mz = Ulst.d * Wl.Vz;
//
//     Ulst.By = Ul.By;
//     Ulst.Bz = Ul.Bz;
//   }
    if (fabs(Ul.d*sdl*sdml-Bxsq) < SMALL_NUMBER*ptst) {
        /* Degenerate case */
        Ulst.My = Ulst.d * Wl.Vy;
        Ulst.Mz = Ulst.d * Wl.Vz;

        Ulst.By = Ul.By;
        Ulst.Bz = Ul.Bz;
    }
    else {
        /* eqns (44) and (46) of M&K */
        tmp = Bxi*(sdl-sdml)/(Ul.d*sdl*sdml-Bxsq);
        Ulst.My = Ulst.d * (Wl.Vy - Ul.By*tmp);
        Ulst.Mz = Ulst.d * (Wl.Vz - Ul.Bz*tmp);
//     if(Ul.By == 0.0 && Ul.Bz == 0.0) {
//       Ulst.By = 0.0;
//       Ulst.Bz = 0.0;
//     }
//     else {
//       /* eqns (45) and (47) of M&K */
//       tmp = (Ul.d*SQR(sdl)-Bxsq)/(Ul.d*sdl*sdml - Bxsq);
//       Ulst.By = Ul.By * tmp;
//       Ulst.Bz = Ul.Bz * tmp;
//     }

        /* eqns (45) and (47) of M&K */
        tmp = (Ul.d*SQR(sdl)-Bxsq)/(Ul.d*sdl*sdml - Bxsq);
        Ulst.By = Ul.By * tmp;
        Ulst.Bz = Ul.Bz * tmp;
    }
    vbstl = (Ulst.Mx*Bxi+Ulst.My*Ulst.By+Ulst.Mz*Ulst.Bz)/Ulst.d;
    /* eqn (48) of M&K */
    Ulst.E = (sdl*Ul.E - ptl*Wl.Vx + ptst*spd[2] +
                        Bxi*(Wl.Vx*Bxi+Wl.Vy*Ul.By+Wl.Vz*Ul.Bz - vbstl))/sdml;
    Wlst = Cons1D_to_Prim1D(Ulst,Bxi,Gamma);


/* Ur* */
    /* eqn (39) of M&K */
    Urst.Mx = Urst.d * spd[2];
//   if((fabs(spd[2]/Wr.Vx-1.0)<SMALL_NUMBER) ||
//      (fabs(spd[2])/fabs(spd[4]) <= SMALL_NUMBER &&
//       fabs(Wr.Vx)/fabs(spd[4]) <= SMALL_NUMBER)) {
//     Urst.My = Urst.d * Wr.Vy;
//     Urst.Mz = Urst.d * Wr.Vz;
//
//     Urst.By = Ur.By;
//     Urst.Bz = Ur.Bz;
//   }
    if (fabs(Ur.d*sdr*sdmr-Bxsq) < SMALL_NUMBER*ptst) {
        /* Degenerate case */
        Urst.My = Urst.d * Wr.Vy;
        Urst.Mz = Urst.d * Wr.Vz;

        Urst.By = Ur.By;
        Urst.Bz = Ur.Bz;
    }
    else {
        /* eqns (44) and (46) of M&K */
        tmp = Bxi*(sdr-sdmr)/(Ur.d*sdr*sdmr-Bxsq);
        Urst.My = Urst.d * (Wr.Vy - Ur.By*tmp);
        Urst.Mz = Urst.d * (Wr.Vz - Ur.Bz*tmp);

//     if(Ur.By == 0.0 && Ur.Bz == 0.0) {
//       Urst.By = 0.0;
//       Urst.Bz = 0.0;
//     }
//     else {
//       /* eqns (45) and (47) of M&K */
//       tmp = (Ur.d*SQR(sdr)-Bxsq)/(Ur.d*sdr*sdmr - Bxsq);
//       Urst.By = Ur.By * tmp;
//       Urst.Bz = Ur.Bz * tmp;
//     }

        /* eqns (45) and (47) of M&K */
        tmp = (Ur.d*SQR(sdr)-Bxsq)/(Ur.d*sdr*sdmr - Bxsq);
        Urst.By = Ur.By * tmp;
        Urst.Bz = Ur.Bz * tmp;
    }
    vbstr = (Urst.Mx*Bxi+Urst.My*Urst.By+Urst.Mz*Urst.Bz)/Urst.d;
    /* eqn (48) of M&K */
    Urst.E = (sdr*Ur.E - ptr*Wr.Vx + ptst*spd[2] +
                        Bxi*(Wr.Vx*Bxi+Wr.Vy*Ur.By+Wr.Vz*Ur.Bz - vbstr))/sdmr;
    Wrst = Cons1D_to_Prim1D(Urst,Bxi,Gamma);


/* Ul** and Ur** - if Bx is zero, same as *-states */
//   if(Bxi == 0.0) {
    if(0.5*Bxsq < SMALL_NUMBER*ptst) {
        Uldst = Ulst;
        Urdst = Urst;
    }
    else {
        invsumd = 1.0/(sqrtdl + sqrtdr);
        if(Bxi > 0.0) Bxsig =  1.0;
        else          Bxsig = -1.0;

        Uldst.d = Ulst.d;
        Urdst.d = Urst.d;

        Uldst.Mx = Ulst.Mx;
        Urdst.Mx = Urst.Mx;

        /* eqn (59) of M&K */
        tmp = invsumd*(sqrtdl*Wlst.Vy + sqrtdr*Wrst.Vy + Bxsig*(Urst.By-Ulst.By));
        Uldst.My = Uldst.d * tmp;
        Urdst.My = Urdst.d * tmp;

        /* eqn (60) of M&K */
        tmp = invsumd*(sqrtdl*Wlst.Vz + sqrtdr*Wrst.Vz + Bxsig*(Urst.Bz-Ulst.Bz));
        Uldst.Mz = Uldst.d * tmp;
        Urdst.Mz = Urdst.d * tmp;

        /* eqn (61) of M&K */
        tmp = invsumd*(sqrtdl*Urst.By + sqrtdr*Ulst.By +
                                     Bxsig*sqrtdl*sqrtdr*(Wrst.Vy-Wlst.Vy));
        Uldst.By = Urdst.By = tmp;

        /* eqn (62) of M&K */
        tmp = invsumd*(sqrtdl*Urst.Bz + sqrtdr*Ulst.Bz +
                                     Bxsig*sqrtdl*sqrtdr*(Wrst.Vz-Wlst.Vz));
        Uldst.Bz = Urdst.Bz = tmp;

        /* eqn (63) of M&K */
        tmp = spd[2]*Bxi + (Uldst.My*Uldst.By + Uldst.Mz*Uldst.Bz)/Uldst.d;
        Uldst.E = Ulst.E - sqrtdl*Bxsig*(vbstl - tmp);
        Urdst.E = Urst.E + sqrtdr*Bxsig*(vbstr - tmp);
    }

/*--- Step 7. ------------------------------------------------------------------
 * Compute flux
 */

    if(spd[1] >= 0.0) {
/* return Fl* */
        printState("L*");
        pFlux.d  = Fl.d  + spd[0]*(Ulst.d  - Ul.d);
        pFlux.Mx = Fl.Mx + spd[0]*(Ulst.Mx - Ul.Mx);
        pFlux.My = Fl.My + spd[0]*(Ulst.My - Ul.My);
        pFlux.Mz = Fl.Mz + spd[0]*(Ulst.Mz - Ul.Mz);
        pFlux.E  = Fl.E  + spd[0]*(Ulst.E  - Ul.E);
        pFlux.By = Fl.By + spd[0]*(Ulst.By - Ul.By);
        pFlux.Bz = Fl.Bz + spd[0]*(Ulst.Bz - Ul.Bz);
    }
    else if(spd[2] >= 0.0) {
/* return Fl** */
        printState("L**");
        tmp = spd[1] - spd[0];
        pFlux.d  = Fl.d  - spd[0]*Ul.d  - tmp*Ulst.d  + spd[1]*Uldst.d;
        pFlux.Mx = Fl.Mx - spd[0]*Ul.Mx - tmp*Ulst.Mx + spd[1]*Uldst.Mx;
        pFlux.My = Fl.My - spd[0]*Ul.My - tmp*Ulst.My + spd[1]*Uldst.My;
        pFlux.Mz = Fl.Mz - spd[0]*Ul.Mz - tmp*Ulst.Mz + spd[1]*Uldst.Mz;
        pFlux.E  = Fl.E  - spd[0]*Ul.E  - tmp*Ulst.E  + spd[1]*Uldst.E;
        pFlux.By = Fl.By - spd[0]*Ul.By - tmp*Ulst.By + spd[1]*Uldst.By;
        pFlux.Bz = Fl.Bz - spd[0]*Ul.Bz - tmp*Ulst.Bz + spd[1]*Uldst.Bz;
    }
    else if(spd[3] > 0.0) {
/* return Fr** */
        printState("R**");
        tmp = spd[3] - spd[4];
        pFlux.d  = Fr.d  - spd[4]*Ur.d  - tmp*Urst.d  + spd[3]*Urdst.d;
        pFlux.Mx = Fr.Mx - spd[4]*Ur.Mx - tmp*Urst.Mx + spd[3]*Urdst.Mx;
        pFlux.My = Fr.My - spd[4]*Ur.My - tmp*Urst.My + spd[3]*Urdst.My;
        pFlux.Mz = Fr.Mz - spd[4]*Ur.Mz - tmp*Urst.Mz + spd[3]*Urdst.Mz;
        pFlux.E  = Fr.E  - spd[4]*Ur.E  - tmp*Urst.E  + spd[3]*Urdst.E;
        pFlux.By = Fr.By - spd[4]*Ur.By - tmp*Urst.By + spd[3]*Urdst.By;
        pFlux.Bz = Fr.Bz - spd[4]*Ur.Bz - tmp*Urst.Bz + spd[3]*Urdst.Bz;
    }
    else {
/* return Fr* */
        printState("R*");
        pFlux.d  = Fr.d  + spd[4]*(Urst.d  - Ur.d);
        pFlux.Mx = Fr.Mx + spd[4]*(Urst.Mx - Ur.Mx);
        pFlux.My = Fr.My + spd[4]*(Urst.My - Ur.My);
        pFlux.Mz = Fr.Mz + spd[4]*(Urst.Mz - Ur.Mz);
        pFlux.E  = Fr.E  + spd[4]*(Urst.E  - Ur.E);
        pFlux.By = Fr.By + spd[4]*(Urst.By - Ur.By);
        pFlux.Bz = Fr.Bz + spd[4]*(Urst.Bz - Ur.Bz);
    }

/* Fluxes of passively advected scalars, computed from density flux */
#if (NSCALARS > 0)
    if (pFlux.d >= 0.0) {
        for (n=0; n<NSCALARS; n++) pFlux.s[n] = pFlux.d*Wl.r[n];
    } else {
        for (n=0; n<NSCALARS; n++) pFlux.s[n] = pFlux.d*Wr.r[n];
    }
#endif
    return;
}
// =============================================================================

// =============================================================================
int main()
{
    // Vectors to store input and output for each test
    std::vector<std::string> names;
    std::vector<Cons1DS> leftConserved, rightConserved, outFlux;
    std::vector<double> gamma;
    std::vector<double> Bx;

    // =========================================================================
    // Brio & Wu
    // Initial Conditions
    // | Field    | Left | Right  |
    // | Density  | 1.0  |  0.125 |
    // | Pressure | 1.0  |  0.1   |
    // | VelX     | 0.0  |  0.0   |
    // | VelY     | 0.0  |  0.0   |
    // | VelZ     | 0.0  |  0.0   |
    // | MagX     | 0.75 |  0.75  |
    // | MagY     | 1.0  | -1.0   |
    // | MagZ     | 0.0  |  0.0   |
    // =========================================================================
    // All ones
    names.push_back("all ones");
    Bx.push_back(1);
    leftConserved.push_back(Cons1DS( 1., 1., 1., 1., 1., Bx.back(), 1., 1.));
    rightConserved.push_back(Cons1DS(1., 1., 1., 1., 1., Bx.back(), 1., 1.));
    gamma.push_back(1.4);

    // =========================================================================
    // End Brio & Wu
    // =========================================================================

    // Check that everything is the same length
    if ( not ((names.size() == leftConserved.size())
               and (names.size() == rightConserved.size())
               and (names.size() == gamma.size())
               and (names.size() == Bx.size())
             )
       )
    {
        throw std::invalid_argument("Not all vectors are the same size");
    }

    outFlux.resize(names.size());
    for (size_t i = 0; i < names.size(); i++)
    {
        // Generate the conserved variables
        Prim1DS primitiveLeft  = Cons1D_to_Prim1D(leftConserved.at(i),  Bx.at(i), gamma.at(i));
        Prim1DS primitiveRight = Cons1D_to_Prim1D(rightConserved.at(i), Bx.at(i), gamma.at(i));

        // Compute fluxes
        fluxes(leftConserved.at(i),
               rightConserved.at(i),
               primitiveLeft,
               primitiveRight,
               Bx.at(i),
               outFlux.at(i),
               gamma.at(i));

        // Return Values
        printResults(leftConserved.at(i), rightConserved.at(i), outFlux.at(i), names.at(i));
    }

    return 0;
}
// =============================================================================
