/* FILE: demagcoef.h            -*-Mode: c++-*-
 *
 * Demag coefficients.
 *
 * Constant magnetization demag routines, based on formulae presented in
 * "A Generalization of the Demagnetizing Tensor for Nonuniform
 * Magnetization," by Andrew J. Newell, Wyn Williams, and David
 * J. Dunlop, Journal of Geophysical Research, vol 98, p 9551-9555, June
 * 1993.  This formulae clearly satisfy necessary symmetry and scaling
 * properties, which is not true of the formulae presented in
 * "Magnetostatic Energy Calculations in Two- and Three-Dimensional
 * Arrays of Ferromagnetic Prisms," M. Maicas, E. Lopez, M. CC. Sanchez,
 * C. Aroca and P. Sanchez, IEEE Trans Mag, vol 34, May 1998, p601-607.
 * (Side note: The units in the latter paper are apparently cgs.)  It
 * appears likely that there is an error in the latter paper (attempts
 * to implement the formulae presented there did not produce the proper
 * symmetries), as well as in the older paper, "Magnetostatic
 * Interaction Fields for a Three-Dimensional Array of Ferromagnetic
 * Cubes," Manfred E. Schabes and Amikam Aharoni, IEEE Trans Mag, vol
 * 23, November 1987, p3882-3888.  (Note: The Newell paper deals with
 * uniformly sized rectangular prisms, the Maicas paper allows
 * non-uniformly sized rectangular prisms, and the Schabes paper only
 * considers cubes.)
 *
 *   The kernel here is based on an analytically derived energy, and the
 * effective (discrete) demag field is calculated from the (discrete)
 * energy.
 *
 */

#include "nb.h"

/* End includes */

// Maybe these routines should be wrapped up in a class?

OC_REALWIDE Newell_f(OC_REALWIDE x,OC_REALWIDE y,OC_REALWIDE z);
OC_REALWIDE Newell_g(OC_REALWIDE x,OC_REALWIDE y,OC_REALWIDE z);

OC_REALWIDE
CalculateSDA00(OC_REALWIDE x, OC_REALWIDE y, OC_REALWIDE z,
	       OC_REALWIDE dx,OC_REALWIDE dy,OC_REALWIDE dz);

inline OC_REALWIDE
CalculateSDA11(OC_REALWIDE x, OC_REALWIDE y, OC_REALWIDE z,
	       OC_REALWIDE dx,OC_REALWIDE dy,OC_REALWIDE dz)
{ return CalculateSDA00(y,x,z,dy,dx,dz); }

inline OC_REALWIDE
CalculateSDA22(OC_REALWIDE x, OC_REALWIDE y, OC_REALWIDE z,
	       OC_REALWIDE dx,OC_REALWIDE dy,OC_REALWIDE dz)
{ return CalculateSDA00(z,y,x,dz,dy,dx); }


OC_REALWIDE
CalculateSDA01(OC_REALWIDE x, OC_REALWIDE y, OC_REALWIDE z,
	       OC_REALWIDE dx,OC_REALWIDE dy,OC_REALWIDE dz);

inline OC_REALWIDE
CalculateSDA02(OC_REALWIDE x, OC_REALWIDE y, OC_REALWIDE z,
	       OC_REALWIDE dx,OC_REALWIDE dy,OC_REALWIDE dz)
{ return CalculateSDA01(x,z,y,dx,dz,dy); }

inline OC_REALWIDE
CalculateSDA12(OC_REALWIDE x, OC_REALWIDE y, OC_REALWIDE z,
	       OC_REALWIDE dx,OC_REALWIDE dy,OC_REALWIDE dz)
{ return CalculateSDA01(y,z,x,dy,dz,dx); }


OC_REALWIDE SelfDemagNx(OC_REALWIDE xsize,OC_REALWIDE ysize,OC_REALWIDE zsize);
OC_REALWIDE SelfDemagNy(OC_REALWIDE xsize,OC_REALWIDE ysize,OC_REALWIDE zsize);
OC_REALWIDE SelfDemagNz(OC_REALWIDE xsize,OC_REALWIDE ysize,OC_REALWIDE zsize);

OC_REALWIDE DemagNxxAsymptotic(OC_REALWIDE x, OC_REALWIDE y, OC_REALWIDE z,
                            OC_REALWIDE dx,OC_REALWIDE dy,OC_REALWIDE dz);

OC_REALWIDE DemagNxyAsymptotic(OC_REALWIDE x, OC_REALWIDE y, OC_REALWIDE z,
                            OC_REALWIDE dx,OC_REALWIDE dy,OC_REALWIDE dz);

OC_REALWIDE DemagNyyAsymptotic(OC_REALWIDE x, OC_REALWIDE y, OC_REALWIDE z,
                            OC_REALWIDE dx,OC_REALWIDE dy,OC_REALWIDE dz);

OC_REALWIDE DemagNzzAsymptotic(OC_REALWIDE x, OC_REALWIDE y, OC_REALWIDE z,
                            OC_REALWIDE dx,OC_REALWIDE dy,OC_REALWIDE dz);

OC_REALWIDE DemagNxzAsymptotic(OC_REALWIDE x, OC_REALWIDE y, OC_REALWIDE z,
                            OC_REALWIDE dx,OC_REALWIDE dy,OC_REALWIDE dz);

OC_REALWIDE DemagNyzAsymptotic(OC_REALWIDE x, OC_REALWIDE y, OC_REALWIDE z,
                            OC_REALWIDE dx,OC_REALWIDE dy,OC_REALWIDE dz);
