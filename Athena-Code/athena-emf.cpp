// From athena4.2/src/integrators/integrate_3d_vl.c lines 2295-2503

static void integrate_emf1_corner(const GridS *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku;
  Real de1_l2, de1_r2, de1_l3, de1_r3;

  il = pG->is-(nghost-1);   iu = pG->ie+(nghost-1);
  jl = pG->js-(nghost-1);   ju = pG->je+(nghost-1);
  kl = pG->ks-(nghost-1);   ku = pG->ke+(nghost-1);

  for (k=kl; k<=ku+1; k++) {
    for (j=jl; j<=ju+1; j++) {
      for (i=il; i<=iu; i++) {
/* NOTE: The x2-Flux of By is -E1. */
/*       The x3-Flux of Bz is +E1. */
	if (x2Flux[k-1][j][i].d > 0.0)
	  de1_l3 = x3Flux[k][j-1][i].Bz - emf1_cc[k-1][j-1][i];
	else if (x2Flux[k-1][j][i].d < 0.0)
	  de1_l3 = x3Flux[k][j][i].Bz - emf1_cc[k-1][j][i];
	else {
	  de1_l3 = 0.5*(x3Flux[k][j-1][i].Bz - emf1_cc[k-1][j-1][i] +
			x3Flux[k][j  ][i].Bz - emf1_cc[k-1][j  ][i] );
	}

	if (x2Flux[k][j][i].d > 0.0)
	  de1_r3 = x3Flux[k][j-1][i].Bz - emf1_cc[k][j-1][i];
	else if (x2Flux[k][j][i].d < 0.0)
	  de1_r3 = x3Flux[k][j][i].Bz - emf1_cc[k][j][i];
	else {
	  de1_r3 = 0.5*(x3Flux[k][j-1][i].Bz - emf1_cc[k][j-1][i] +
			x3Flux[k][j  ][i].Bz - emf1_cc[k][j  ][i] );
	}

	if (x3Flux[k][j-1][i].d > 0.0)
	  de1_l2 = -x2Flux[k-1][j][i].By - emf1_cc[k-1][j-1][i];
	else if (x3Flux[k][j-1][i].d < 0.0)
	  de1_l2 = -x2Flux[k][j][i].By - emf1_cc[k][j-1][i];
	else {
	  de1_l2 = 0.5*(-x2Flux[k-1][j][i].By - emf1_cc[k-1][j-1][i]
			-x2Flux[k  ][j][i].By - emf1_cc[k  ][j-1][i] );
	}

	if (x3Flux[k][j][i].d > 0.0)
	  de1_r2 = -x2Flux[k-1][j][i].By - emf1_cc[k-1][j][i];
	else if (x3Flux[k][j][i].d < 0.0)
	  de1_r2 = -x2Flux[k][j][i].By - emf1_cc[k][j][i];
	else {
	  de1_r2 = 0.5*(-x2Flux[k-1][j][i].By - emf1_cc[k-1][j][i]
			-x2Flux[k  ][j][i].By - emf1_cc[k  ][j][i] );
	}

        emf1[k][j][i] = 0.25*(  x3Flux[k][j][i].Bz + x3Flux[k][j-1][i].Bz
                              - x2Flux[k][j][i].By - x2Flux[k-1][j][i].By
			      + de1_l2 + de1_r2 + de1_l3 + de1_r3);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void integrate_emf2_corner(const GridS *pG)
 *  \brief Integrates face centered B-fluxes to compute corner EMFs.
 *
 *   Note:
 *-  x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 *-  x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 *-  x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 *-  x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 *-  x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 *-  x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX
 *- first_order_correction() - Added by N. Lemaster to run supersonic turbulence
 */
static void integrate_emf2_corner(const GridS *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku;
  Real de2_l1, de2_r1, de2_l3, de2_r3;

  il = pG->is-(nghost-1);   iu = pG->ie+(nghost-1);
  jl = pG->js-(nghost-1);   ju = pG->je+(nghost-1);
  kl = pG->ks-(nghost-1);   ku = pG->ke+(nghost-1);

  for (k=kl; k<=ku+1; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu+1; i++) {
/* NOTE: The x1-Flux of Bz is +E2. */
/*       The x3-Flux of By is -E2. */
	if (x1Flux[k-1][j][i].d > 0.0)
	  de2_l3 = -x3Flux[k][j][i-1].By - emf2_cc[k-1][j][i-1];
	else if (x1Flux[k-1][j][i].d < 0.0)
	  de2_l3 = -x3Flux[k][j][i].By - emf2_cc[k-1][j][i];
	else {
	  de2_l3 = 0.5*(-x3Flux[k][j][i-1].By - emf2_cc[k-1][j][i-1]
			-x3Flux[k][j][i  ].By - emf2_cc[k-1][j][i  ] );
	}

	if (x1Flux[k][j][i].d > 0.0)
	  de2_r3 = -x3Flux[k][j][i-1].By - emf2_cc[k][j][i-1];
	else if (x1Flux[k][j][i].d < 0.0)
	  de2_r3 = -x3Flux[k][j][i].By - emf2_cc[k][j][i];
	else {
	  de2_r3 = 0.5*(-x3Flux[k][j][i-1].By - emf2_cc[k][j][i-1]
			-x3Flux[k][j][i  ].By - emf2_cc[k][j][i  ] );
	}

	if (x3Flux[k][j][i-1].d > 0.0)
	  de2_l1 = x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i-1];
	else if (x3Flux[k][j][i-1].d < 0.0)
	  de2_l1 = x1Flux[k][j][i].Bz - emf2_cc[k][j][i-1];
	else {
	  de2_l1 = 0.5*(x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i-1] +
			x1Flux[k  ][j][i].Bz - emf2_cc[k  ][j][i-1] );
	}

	if (x3Flux[k][j][i].d > 0.0)
	  de2_r1 = x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i];
	else if (x3Flux[k][j][i].d < 0.0)
	  de2_r1 = x1Flux[k][j][i].Bz - emf2_cc[k][j][i];
	else {
	  de2_r1 = 0.5*(x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i] +
			x1Flux[k  ][j][i].Bz - emf2_cc[k  ][j][i] );
	}

	emf2[k][j][i] = 0.25*(  x1Flux[k][j][i].Bz + x1Flux[k-1][j][i  ].Bz
                              - x3Flux[k][j][i].By - x3Flux[k  ][j][i-1].By
			      + de2_l1 + de2_r1 + de2_l3 + de2_r3);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void integrate_emf3_corner(const GridS *pG)
 *  \brief Integrates face centered B-fluxes to compute corner EMFs.
 *
 *   Note:
 *-  x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 *-  x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 *-  x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 *-  x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 *-  x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 *-  x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX
 *- first_order_correction() - Added by N. Lemaster to run supersonic turbulence
 */
static void integrate_emf3_corner(const GridS *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku;
  Real de3_l1, de3_r1, de3_l2, de3_r2;
  Real rsf=1.0,lsf=1.0;

  il = pG->is-(nghost-1);   iu = pG->ie+(nghost-1);
  jl = pG->js-(nghost-1);   ju = pG->je+(nghost-1);
  kl = pG->ks-(nghost-1);   ku = pG->ke+(nghost-1);

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju+1; j++) {
      for (i=il; i<=iu+1; i++) {
#ifdef CYLINDRICAL
        rsf = pG->ri[i]/pG->r[i];  lsf = pG->ri[i]/pG->r[i-1];
#endif
/* NOTE: The x1-Flux of By is -E3. */
/*       The x2-Flux of Bx is +E3. */
	if (x1Flux[k][j-1][i].d > 0.0)
	  de3_l2 = (x2Flux[k][j][i-1].Bz - emf3_cc[k][j-1][i-1])*lsf;
	else if (x1Flux[k][j-1][i].d < 0.0)
	  de3_l2 = (x2Flux[k][j][i].Bz - emf3_cc[k][j-1][i])*rsf;
	else {
	  de3_l2 = 0.5*((x2Flux[k][j][i-1].Bz - emf3_cc[k][j-1][i-1])*lsf +
			(x2Flux[k][j][i  ].Bz - emf3_cc[k][j-1][i  ])*rsf );
	}

	if (x1Flux[k][j][i].d > 0.0)
	  de3_r2 = (x2Flux[k][j][i-1].Bz - emf3_cc[k][j][i-1])*lsf;
	else if (x1Flux[k][j][i].d < 0.0)
	  de3_r2 = (x2Flux[k][j][i].Bz - emf3_cc[k][j][i])*rsf;
	else {
	  de3_r2 = 0.5*((x2Flux[k][j][i-1].Bz - emf3_cc[k][j][i-1])*lsf +
			(x2Flux[k][j][i  ].Bz - emf3_cc[k][j][i  ])*rsf );
	}

	if (x2Flux[k][j][i-1].d > 0.0)
	  de3_l1 = -x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i-1];
	else if (x2Flux[k][j][i-1].d < 0.0)
	  de3_l1 = -x1Flux[k][j][i].By - emf3_cc[k][j][i-1];
	else {
	  de3_l1 = 0.5*(-x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i-1]
			-x1Flux[k][j  ][i].By - emf3_cc[k][j  ][i-1] );
	}

	if (x2Flux[k][j][i].d > 0.0)
	  de3_r1 = -x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i];
	else if (x2Flux[k][j][i].d < 0.0)
	  de3_r1 = -x1Flux[k][j][i].By - emf3_cc[k][j][i];
	else {
	  de3_r1 = 0.5*(-x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i]
			-x1Flux[k][j  ][i].By - emf3_cc[k][j  ][i] );
	}

	emf3[k][j][i] = 0.25*(  x2Flux[k][j  ][i-1].Bz + x2Flux[k][j][i].Bz
			      - x1Flux[k][j-1][i  ].By - x1Flux[k][j][i].By
			      + de3_l1 + de3_r1 + de3_l2 + de3_r2);
      }
    }
  }

  return;
}