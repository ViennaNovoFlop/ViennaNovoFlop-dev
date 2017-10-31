/* FILE: xpfloat.cc                    -*-Mode: c++-*-
 *
 * Class for extended floating point precision using
 * compensated summation.
 * 
 * NOTICE: Please see the file ../../LICENSE
 *
 * Last modified on: $Date: 2010-07-16 22:33:58 $
 * Last modified by: $Author: donahue $
 */

#include "xpfloat.h"

/* End includes */     

void Nb_Xpfloat::Accum(const Nb_Xpfloat &y)
{
  if((x>0 && y.x>0) || (x<0 && y.x<0)) {
    Accum(y.corr);  Accum(y.x);
  } else {
    Accum(y.x);  Accum(y.corr);
  }
}

Nb_Xpfloat& Nb_Xpfloat::operator-=(const Nb_Xpfloat &y)
{
  if((x>0 && y.x<0) || (x<0 && y.x>0)) {
    Accum(-y.corr);  Accum(-y.x);
  } else {
    Accum(-y.x);  Accum(-y.corr);
  }
  return *this;
}

Nb_Xpfloat operator+(const Nb_Xpfloat& x,const Nb_Xpfloat& y)
{
  Nb_Xpfloat z(x);
  z.Accum(y);
  return z;
}

Nb_Xpfloat operator-(const Nb_Xpfloat& x,const Nb_Xpfloat& y)
{
  Nb_Xpfloat z(y);
  z *= -1.;
  z.Accum(x);
  return z;
}

Nb_Xpfloat operator*(const Nb_Xpfloat& x,OC_REALWIDE y)
{
  Nb_Xpfloat z(x);
  z*=y;
  return z;
}

Nb_Xpfloat& Nb_Xpfloat::operator*=(const Nb_Xpfloat& y)
{
  Nb_Xpfloat z(corr*y.corr);
  z.Accum(x*y.corr);
  z.Accum(y.x*corr);
  z.Accum(x*y.x);

  x = z.x;
  corr = z.corr;

  return *this;
}

Nb_Xpfloat operator*(const Nb_Xpfloat& x,const Nb_Xpfloat& y)
{
  Nb_Xpfloat z(x.corr*y.corr);
  z.Accum(x.x*y.corr);
  z.Accum(y.x*x.corr);
  z.Accum(x.x*y.x);
  return z;
}
