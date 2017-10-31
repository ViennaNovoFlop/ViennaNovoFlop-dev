/* FILE: xpfloat.h                    -*-Mode: c++-*-
 *
 * Class for extended floating point precision using
 * compensated summation.
 * 
 * NOTICE: Please see the file ../../LICENSE
 *
 * Last modified on: $Date: 2010-07-16 22:33:59 $
 * Last modified by: $Author: donahue $
 */

#include <assert.h>

#include "oc.h"

#ifndef _NB_XPFLOAT_H
#define _NB_XPFLOAT_H

#if OC_USE_SSE
# include <emmintrin.h>
#endif

class Nb_Xpfloat
{
  friend void Nb_XpfloatDualAccum(Nb_Xpfloat&,OC_REALWIDE,
                                  Nb_Xpfloat&,OC_REALWIDE);
#if OC_USE_SSE
  friend void Nb_XpfloatDualAccum(Nb_Xpfloat&,Nb_Xpfloat&,__m128d);
#endif
private:
  OC_REALWIDE x,corr;
public:
  Nb_Xpfloat(OC_REALWIDE newx=0.0) { x=newx; corr=0.; }
  Nb_Xpfloat(const Nb_Xpfloat& y) { x=y.x; corr=y.corr; }

  void Set(OC_REALWIDE newx) { x=newx; corr=0.; }
  Nb_Xpfloat& operator=(OC_REALWIDE y) { Set(y); return *this; }
  Nb_Xpfloat& operator=(const Nb_Xpfloat& y) {
    x=y.x; corr=y.corr; return *this;
  }

  OC_REALWIDE GetValue() const { return x; }
  operator OC_REALWIDE() { return GetValue(); }

  // Defining Accum(OC_REALWIDE) here, as opposed to in the xpfloat.cc
  // file, allows significant speedups with the g++ compiler.  It makes
  // small difference with the Intel C++ compiler (icpc), presumably
  // because of icpc's interprocedural optimization ability?
  void Accum(OC_REALWIDE y) {
#if !OC_USE_SSE
    volatile OC_REALWIDE sum1; // Mark as volatile to protect against
    volatile OC_REALWIDE sum2; // problems arising from extra precision.
    volatile OC_REALWIDE corrtemp;

    // Calculate sum
    sum1=y+corr;
    sum2=x+sum1;

    // Determine new correction term.  Do in two steps, as opposed to
    // "corrtemp = (x-sum2) + sum1", because a) some compilers ignore
    // parentheses (with regards to ordering) at high optimization levels,
    // and b) we want to be sure to drop any extra precision.
    corrtemp = x - sum2;
    corrtemp += sum1;
    
    // Store results
    corr = corrtemp;
    x=sum2;
#else // SSE
    // SSE implementation of the above.  Since SSE double precision
    // doesn't carry any extra precision, we don't need to use "volatile"
    __m128d wx = _mm_set_sd(x);
    __m128d sum1 = _mm_add_sd(_mm_set_sd(y),_mm_set_sd(corr));
    __m128d sum2 = _mm_add_sd(wx,sum1);
    __m128d corrtemp = _mm_add_sd(_mm_sub_sd(wx,sum2),sum1);
    corr = _mm_cvtsd_f64(corrtemp);
    x = _mm_cvtsd_f64(sum2);
#endif // SSE
  }

  Nb_Xpfloat& operator+=(OC_REALWIDE y) { Accum(y); return *this; }
  Nb_Xpfloat& operator-=(OC_REALWIDE y) { Accum(-y); return *this; }
  Nb_Xpfloat& operator*=(OC_REALWIDE y) { x*=y; corr*=y; return *this; }


  void Accum(const Nb_Xpfloat& y);
  Nb_Xpfloat& operator+=(const Nb_Xpfloat& y) { Accum(y); return *this; }
  Nb_Xpfloat& operator-=(const Nb_Xpfloat& y);
  Nb_Xpfloat& operator*=(const Nb_Xpfloat& y);

  friend Nb_Xpfloat operator*(const Nb_Xpfloat& x,const Nb_Xpfloat& y);
};

Nb_Xpfloat operator+(const Nb_Xpfloat& x,const Nb_Xpfloat& y);
Nb_Xpfloat operator-(const Nb_Xpfloat& x,const Nb_Xpfloat& y);
Nb_Xpfloat operator*(const Nb_Xpfloat& x,OC_REALWIDE y);
inline Nb_Xpfloat operator*(OC_REALWIDE y,const Nb_Xpfloat& x) { return x*y; }
Nb_Xpfloat operator*(const Nb_Xpfloat& x,const Nb_Xpfloat& y);

#if OC_USE_SSE
inline void
Nb_XpfloatDualAccum(Nb_Xpfloat& xpA,Nb_Xpfloat& xpB,__m128d y)
{ // Import y holds two packed double precision values.  The
  // lower half of y holds the value to be accumulated into xpA,
  // and the upper half holds the value to be accumulated into xpB.

  // IMPORTANT NOTE : This code assumes that xpA and xpB are
  // ***DIFFERENT*** objects.  If they are the same object, then one of
  // the accums will be lost!!!!!!!

  // For implementation notes, see the Nb_Xpfloat::Accum routine above.
  assert(&xpA != &xpB);  // See IMPORTANT NOTE above.

  __m128d wx = _mm_set_pd(xpB.x,xpA.x);
  __m128d sum1 = _mm_add_pd(y,_mm_set_pd(xpB.corr,xpA.corr));
  __m128d sum2 = _mm_add_pd(wx,sum1);
  __m128d corrtemp = _mm_add_pd(_mm_sub_pd(wx,sum2),sum1);

  _mm_storel_pd(&xpA.x,sum2);
  _mm_storeh_pd(&xpB.x,sum2);
  _mm_storel_pd(&xpA.corr,corrtemp);
  _mm_storeh_pd(&xpB.corr,corrtemp);

}
#endif

inline void
Nb_XpfloatDualAccum(Nb_Xpfloat& xpA,OC_REALWIDE yA,
                    Nb_Xpfloat& xpB,OC_REALWIDE yB)
{ // Updates two Nb_Xpfloat objects at the same time.  The main
  // reason for the existence of this function is double-scalar
  // processing offered in the SSE implementation.

  // IMPORTANT NOTE : This code assumes that xpA and xpB are
  // ***DIFFERENT*** objects.  If they are the same object, then one of
  // the accums will be lost!!!!!!!

  // For implementation notes, see the Nb_Xpfloat::Accum routine above.


  assert(&xpA != &xpB);  // See IMPORTANT NOTE above.

#if !OC_USE_SSE
  volatile OC_REALWIDE sumA1,     sumB1;
  volatile OC_REALWIDE sumA2,     sumB2;
  volatile OC_REALWIDE corrtempA, corrtempB;

  sumA1=yA+xpA.corr;           sumB1=yB+xpB.corr;
  sumA2=xpA.x+sumA1;           sumB2=xpB.x+sumB1;
  corrtempA = xpA.x - sumA2;   corrtempB = xpB.x - sumB2;
  corrtempA += sumA1;          corrtempB += sumB1;

  xpA.corr = corrtempA;
  xpA.x=sumA2;

  xpB.corr = corrtempB;
  xpB.x=sumB2;

#else // SSE

  Nb_XpfloatDualAccum(xpA,xpB,_mm_set_pd(yB,yA));

#endif // SSE
}

#endif // _NB_XPFLOAT_H
