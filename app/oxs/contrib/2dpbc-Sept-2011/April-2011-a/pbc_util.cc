#include <math.h>
#include <stdlib.h>   
#include "nb.h"
#include "pbc_util.h"
#include "demagcoef.h"

#define ADD_INF 1
////////////////////////////////////////////////////////////////////////////
// Routines to do accurate summation
extern "C" {
static int AS_Compare(const void* px,const void* py)
{
  // Comparison based on absolute values
  double x=fabs(*((const double *)px));
  double y=fabs(*((const double *)py));
  if(x<y) return 1;
  if(x>y) return -1;
  return 0;
}
}

double
AccurateSum(int n,double *arr)
{
  // Order by decreasing magnitude
  qsort(arr,n,sizeof(double),AS_Compare);

  double sum,corr,y,u,t,v,z,x,tmp;

  sum=arr[0]; corr=0;
  for(int i=1;i<n;i++) {
    x=arr[i];
    y=corr+x;
    tmp=y-corr;
    u=x-tmp;
    t=y+sum;
    tmp=t-sum;
    v=y-tmp;
    z=u+v;
    sum=t+z;
    tmp=sum-t;
    corr=z-tmp;
  }
  return sum;
}	


double
DemagTensorNormal(enum TensorComponent comp,double x,double y,double z,double a,double b,double c)
{
	switch(comp){
		case xx:return CalculateSDA00(x,y,z,a,b,c);
		case yy:return CalculateSDA00(y,x,z,b,a,c);
		case zz:return CalculateSDA00(z,y,x,c,b,a);
		case xy:return CalculateSDA01(x,y,z,a,b,c);
		case xz:return CalculateSDA01(x,z,y,a,c,b);
		case yz:return CalculateSDA01(y,z,x,b,c,a);
	}
}


double
DemagTensorAsymptotic(enum TensorComponent comp,double x,double y,double z,double a,double b,double c)
{
	switch(comp){
		case xx:return DemagNxxAsymptotic(x,y,z,a,b,c);
		case yy:return DemagNxxAsymptotic(y,x,z,b,a,c);
		case zz:return DemagNxxAsymptotic(z,y,x,c,b,a);
		case xy:return DemagNxyAsymptotic(x,y,z,a,b,c);
		case xz:return DemagNxyAsymptotic(x,z,y,a,c,b);
		case yz:return DemagNxyAsymptotic(y,z,x,b,c,a);
	}
}

double 
DemagTensorDipolar(enum TensorComponent comp,double x,double y,double z){
switch(comp)
{
	case xx:return DemagNxxDipolar(x,y,z);
  case yy:return DemagNxxDipolar(y,x,z);
	case zz:return DemagNxxDipolar(z,y,x);
	case xy:return DemagNxyDipolar(x,y,z);
	case xz:return DemagNxyDipolar(x,z,y);
	case yz:return DemagNxyDipolar(y,z,x);
}
}


double 
DemagTensorInfinite(enum TensorComponent comp,double x,double y,double z,double X0,double Y0)
{
	switch(comp){
		case xx:return Nxxinf(x,y,z,X0,Y0);
		case yy:return Nyyinf(x,y,z,X0,Y0);
		case zz:return Nzzinf(x,y,z,X0,Y0);
		case xy:return Nxyinf(x,y,z,X0,Y0);
		case xz:return Nxzinf(x,y,z,X0,Y0);
		case yz:return Nyzinf(x,y,z,X0,Y0);
	}

}

REALWIDE DemagNxxAsymptotic(REALWIDE x, REALWIDE y, REALWIDE z,
                            REALWIDE dx,REALWIDE dy,REALWIDE dz)
{ // Asymptotic approximation to Nxx term.
  REALWIDE R2 = x*x + y*y + z*z;

  if(R2<=0.0) {
    // Asymptotic expansion doesn't apply for R==0.  Fall back
    // to self-demag calculation.
    return SelfDemagNx(dx,dy,dz);
  }

  REALWIDE R  = sqrt(R2);

  REALWIDE sx2 = x*x/R2;
  REALWIDE sy2 = y*y/R2;
  REALWIDE sz2 = z*z/R2;

  REALWIDE sx4 = sx2*sx2;
  REALWIDE sy4 = sy2*sy2; 
  REALWIDE sz4 = sz2*sz2;

  REALWIDE sx6 = sx4*sx2;
  REALWIDE sy6 = sy4*sy2; 
  REALWIDE sz6 = sz4*sz2;

  REALWIDE dx2 = dx*dx;
  REALWIDE dy2 = dy*dy;
  REALWIDE dz2 = dz*dz;

  REALWIDE dx4 = dx2*dx2;
  REALWIDE dy4 = dy2*dy2;
  REALWIDE dz4 = dz2*dz2;

  REALWIDE term3 = 2*sx2 - sy2 - sz2;

  REALWIDE term5 = 0.0;
  if(dx2!=dy2 || dx2!=dz2 || dy2!=dz2) { // Non-cube case
    REALWIDE a1 =   8*dx2  -  4*dy2  -  4*dz2;
    REALWIDE a2 = -24*dx2  + 27*dy2  -  3*dz2;
    REALWIDE a3 = -24*dx2  -  3*dy2  + 27*dz2;
    REALWIDE a4 =   3*dx2  -  4*dy2  +  1*dz2;
    REALWIDE a5 =   6*dx2  -  3*dy2  -  3*dz2;
    REALWIDE a6 =   3*dx2  +  1*dy2  -  4*dz2;
    term5 = a1*sx4 + a2*sx2*sy2 + a3*sx2*sz2
      + a4*sy4 + a5*sy2*sz2 + a6*sz4;
    term5 *= 0.25;
  }

  REALWIDE term7 = 0.0;
  {
    REALWIDE b1  =   32*dx4  -  40*dx2*dy2  -  40*dx2*dz2   +  12*dy4   +  10*dy2*dz2   +  12*dz4;
    REALWIDE b2  = -240*dx4  + 580*dx2*dy2  +  20*dx2*dz2   - 202*dy4   -  75*dy2*dz2   +  22*dz4;
    REALWIDE b3  = -240*dx4  +  20*dx2*dy2  + 580*dx2*dz2   +  22*dy4   -  75*dy2*dz2   - 202*dz4;
    REALWIDE b4  =  180*dx4  - 505*dx2*dy2  +  55*dx2*dz2   + 232*dy4   -  75*dy2*dz2   +   8*dz4;
    REALWIDE b5  =  360*dx4  - 450*dx2*dy2  - 450*dx2*dz2   - 180*dy4   + 900*dy2*dz2   - 180*dz4;
    REALWIDE b6  =  180*dx4  +  55*dx2*dy2  - 505*dx2*dz2   +   8*dy4   -  75*dy2*dz2   + 232*dz4;
    REALWIDE b7  =  -10*dx4  +  30*dx2*dy2  -   5*dx2*dz2   -  16*dy4   +  10*dy2*dz2   -   2*dz4;
    REALWIDE b8  =  -30*dx4  +  55*dx2*dy2  +  20*dx2*dz2   +   8*dy4   -  75*dy2*dz2   +  22*dz4;
    REALWIDE b9  =  -30*dx4  +  20*dx2*dy2  +  55*dx2*dz2   +  22*dy4   -  75*dy2*dz2   +   8*dz4;
    REALWIDE b10 =  -10*dx4  -   5*dx2*dy2  +  30*dx2*dz2   -   2*dy4   +  10*dy2*dz2   -  16*dz4;

    term7 = b1*sx6 + b2*sx4*sy2 + b3*sx4*sz2 + b4*sx2*sy4 + b5*sx2*sy2*sz2
      + b6*sx2*sz4 + b7*sy6 + b8*sy4*sz2 + b9*sy2*sz4 + b10*sz6;
    term7 *= (REALWIDE(1.)/REALWIDE(16.));
  }

  REALWIDE Nxx = (-dx*dy*dz/(4*WIDE_PI)) * (((term7/R2 + term5)/R2 + term3)/(R2*R));
  // Error should be of order 1/R^9

  return Nxx;
}

REALWIDE DemagNxyAsymptotic(REALWIDE x, REALWIDE y, REALWIDE z,
                            REALWIDE dx,REALWIDE dy,REALWIDE dz)
{ // Asympotic approximation to Nxy term.
  REALWIDE R2 = x*x + y*y + z*z;

  if(R2<=0.0) {
    // Asymptotic expansion doesn't apply for R==0.  Fall back
    // to self-demag calculation.
    return static_cast<REALWIDE>(0.0);
  }

  REALWIDE R  = sqrt(R2);

  REALWIDE sx2 = x*x/R2;
  REALWIDE sy2 = y*y/R2;
  REALWIDE sz2 = z*z/R2;

  REALWIDE sx4 = sx2*sx2;
  REALWIDE sy4 = sy2*sy2; 
  REALWIDE sz4 = sz2*sz2;

  REALWIDE dx2 = dx*dx;
  REALWIDE dy2 = dy*dy;
  REALWIDE dz2 = dz*dz;

  REALWIDE dx4 = dx2*dx2;
  REALWIDE dy4 = dy2*dy2;
  REALWIDE dz4 = dz2*dz2;

  REALWIDE term3 = 3;

  REALWIDE term5 = 0.0;
  if(dx2!=dy2 || dx2!=dz2 || dy2!=dz2) { // Non-cube case
    REALWIDE a1 =   4*dx2  -  3*dy2  -  1*dz2;
    REALWIDE a2 =  -3*dx2  +  4*dy2  -  1*dz2;
    REALWIDE a3 =  -3*dx2  -  3*dy2  +  6*dz2;
    term5 = a1*sx2 + a2*sy2 + a3*sz2;
    term5 *= REALWIDE(5.)/REALWIDE(4.);
  }

  REALWIDE term7 = 0.0;
  {
    REALWIDE b1  =   16*dx4  -  30*dx2*dy2  -  10*dx2*dz2   +  10*dy4   +   5*dy2*dz2   +   2*dz4;
    REALWIDE b2  =  -40*dx4  + 105*dx2*dy2  -   5*dx2*dz2   -  40*dy4   -   5*dy2*dz2   +   4*dz4;
    REALWIDE b3  =  -40*dx4  -  15*dx2*dy2  + 115*dx2*dz2   +  20*dy4   -  35*dy2*dz2   -  32*dz4;
    REALWIDE b4  =   10*dx4  -  30*dx2*dy2  +   5*dx2*dz2   +  16*dy4   -  10*dy2*dz2   +   2*dz4;
    REALWIDE b5  =   20*dx4  -  15*dx2*dy2  -  35*dx2*dz2   -  40*dy4   + 115*dy2*dz2   -  32*dz4;
    REALWIDE b6  =   10*dx4  +  15*dx2*dy2  -  40*dx2*dz2   +  10*dy4   -  40*dy2*dz2   +  32*dz4;
    term7 = b1*sx4 + b2*sx2*sy2 + b3*sx2*sz2
      + b4*sy4 + b5*sy2*sz2 + b6*sz4;
    term7 *= (REALWIDE(7.)/REALWIDE(16.));
  }

  REALWIDE Nxy = (-dx*dy*dz*x*y/(4*WIDE_PI*R2)) * (((term7/R2 + term5)/R2 + term3)/(R2*R));
  // Error should be of order 1/R^9

  return Nxy;
}
