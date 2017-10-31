/* FILE: threevector.h        -*-Mode: c++-*-
 *
 * Simple 3D vector class.  *All* client code should use the
 * ThreeVector typedef in case this class gets moved into an
 * OOMMF extension library.
 *
 */

#ifndef _OXS_THREEVECTOR
#define _OXS_THREEVECTOR

#if OC_USE_SSE
# include <emmintrin.h>
#endif

#include "oc.h"

OC_USE_STD_NAMESPACE;  // For std::fabs()

/* End includes */

class Oxs_ThreeVector;  // Forward declaration for typedef
typedef Oxs_ThreeVector ThreeVector;

class Oxs_ThreeVector {
public:
  OC_REAL8m x,y,z;

  // Constructors
  Oxs_ThreeVector(): x(0.), y(0.), z(0.) {}
  Oxs_ThreeVector(OC_REAL8m xi,OC_REAL8m yi,OC_REAL8m zi)
    : x(xi), y(yi), z(zi) {}
  Oxs_ThreeVector(OC_REAL8m *arr) : x(arr[0]), y(arr[1]), z(arr[2]) {}
  Oxs_ThreeVector(const Oxs_ThreeVector &v) : x(v.x), y(v.y), z(v.z) {}

  // Assignment operators
  Oxs_ThreeVector& Set(OC_REAL8m xi,OC_REAL8m yi,OC_REAL8m zi)
    { x=xi; y=yi; z=zi; return *this; }
  Oxs_ThreeVector& operator=(const Oxs_ThreeVector &v) {
    x=v.x; y=v.y; z=v.z;
    return *this;
  }
  Oxs_ThreeVector& operator+=(const Oxs_ThreeVector& v) { 
    x+=v.x; y+=v.y; z+=v.z;
    return *this;
  }
  Oxs_ThreeVector& operator-=(const Oxs_ThreeVector& v) {
    x-=v.x; y-=v.y; z-=v.z;
    return *this;
  }

  Oxs_ThreeVector& operator^=(const Oxs_ThreeVector& w) { // Cross product
#if 0
    OC_REAL8m tx = y * w.z  -  z * w.y;
    OC_REAL8m ty = z * w.x  -  x * w.z;
    OC_REAL8m tz = x * w.y  -  y * w.x;
    x = tx; y = ty; z = tz;
#else // In principle, the ordering below reduces register pressure
    OC_REAL8m p12 = x*w.y;
    OC_REAL8m p13 = x*w.z;
    OC_REAL8m p23 = y*w.z;
    p12 -= y*w.x;
    y = z*w.x - p13;
    x = p23 - z*w.y;
    z = p12;
#endif

    return *this;
  }

  Oxs_ThreeVector& operator*=(OC_REAL8m a) { x*=a; y*=a; z*=a; return *this; }

  Oxs_ThreeVector Accum(OC_REAL8m a,const Oxs_ThreeVector& v) {
    x+=a*v.x;  y+=a*v.y;  z+=a*v.z;  return *this;
  }


  Oxs_ThreeVector& wxvxw(const Oxs_ThreeVector& w) {
    // Performs w x *this x w

    OC_REAL8m tx = y * w.z  -  z * w.y;
    OC_REAL8m ty = z * w.x  -  x * w.z;
    z = ty * w.x  - tx * w.y ;

    OC_REAL8m tz = x * w.y  -  y * w.x;
    x = tz * w.y  - ty * w.z ;
    y = tx * w.z  - tz * w.x ;

    return *this;
  }

  // Misc
  OC_REAL8m MagSq() const { return x*x+y*y+z*z; }
  OC_REAL8m TaxicabNorm() const { return fabs(x)+fabs(y)+fabs(z); }  // aka L1-norm
  void SetMag(OC_REAL8m mag);  // Adjusts size to "mag"
  void Random(OC_REAL8m mag);  // Makes a random vector of size "mag"
  OC_REAL8m MakeUnit(); // High-precision version of SetMag(1.0).
  /// Return value is original MagSq(), i.e., on entry.
};

// Test operators on Oxs_ThreeVector's
inline OC_BOOL operator==(const Oxs_ThreeVector& lhs,
		const Oxs_ThreeVector& rhs)
{ return ((lhs.x==rhs.x) && (lhs.y==rhs.y) && (lhs.z==rhs.z)); }

inline OC_BOOL operator!=(const Oxs_ThreeVector& lhs,
		const Oxs_ThreeVector& rhs)
{ return ((lhs.x!=rhs.x) || (lhs.y!=rhs.y) || (lhs.z!=rhs.z)); }

// The next two operators are defined so MSVC++ 5.0 will accept
// vector<ThreeVector>, but are left undefined because there
// is no way to define them that makes sense in all contexts.
OC_BOOL operator<(const Oxs_ThreeVector&,const Oxs_ThreeVector&);
OC_BOOL operator>(const Oxs_ThreeVector&,const Oxs_ThreeVector&);


// Binary operators on Oxs_ThreeVector's
inline const Oxs_ThreeVector
operator+(const Oxs_ThreeVector& lhs,
	  const Oxs_ThreeVector& rhs)
{ Oxs_ThreeVector result(lhs); return result+=rhs; }

inline const Oxs_ThreeVector
operator-(const Oxs_ThreeVector& lhs,
	  const Oxs_ThreeVector& rhs)
{ Oxs_ThreeVector result(lhs); return result-=rhs; }

// Cross product
inline const Oxs_ThreeVector
operator^(const Oxs_ThreeVector& v,
	  const Oxs_ThreeVector& w)
{
#if 0
  OC_REAL8m tx = v.y * w.z  -  w.y * v.z;
  OC_REAL8m ty = w.x * v.z  -  v.x * w.z;
  OC_REAL8m tz = v.x * w.y  -  w.x * v.y;
#else // In principle, the ordering below reduces register pressure
  OC_REAL8m tz = v.x*w.y;
  OC_REAL8m ty = v.x*w.z;
  OC_REAL8m tx = v.y*w.z;
  tz -= v.y*w.x;
  ty  = v.z*w.x - ty;
  tx -= v.z*w.y;
#endif
  return Oxs_ThreeVector(tx,ty,tz);
}


// Product against scalar
inline const Oxs_ThreeVector
operator*(OC_REAL8m scalar,
	  const Oxs_ThreeVector& vec)
{ Oxs_ThreeVector result(vec); return result*=scalar; }

inline const Oxs_ThreeVector
operator*(const Oxs_ThreeVector& vec,
	  OC_REAL8m scalar)
{ Oxs_ThreeVector result(vec); return result*=scalar; }


// Dot product
inline OC_REAL8m
operator*(const Oxs_ThreeVector& lhs,const Oxs_ThreeVector& rhs)
{ return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z; }

////////////////////////////////////////////////////////////////////////
// "Pair" routines, that process two vectors at a time.  If SSE is
// available, then these versions can be significantly faster than the
// serial routines.
#if OC_USE_SSE
inline __m128d Oxs_ThreeVectorPairMagSq
(const __m128d vx,const __m128d vy,const __m128d vz)
{ // Returns |vec0|^2 in in lower half, |vec1|^2 in upper half.
  __m128d txsq = _mm_mul_pd(vx,vx);
  __m128d tysq = _mm_mul_pd(vy,vy);
  __m128d tzsq = _mm_mul_pd(vz,vz);
  return _mm_add_pd(_mm_add_pd(txsq,tysq),tzsq);
}

inline __m128d Oxs_ThreeVectorPairMagSq
(const Oxs_ThreeVector& vec0,
 const Oxs_ThreeVector& vec1)
{ // Returns |vec0|^2 in in lower half, |vec1|^2 in upper half.
  return Oxs_ThreeVectorPairMagSq(_mm_set_pd(vec1.x,vec0.x),
                                  _mm_set_pd(vec1.y,vec0.y),
                                  _mm_set_pd(vec1.z,vec0.z));
}
#endif // OC_USE_SSE

inline void Oxs_ThreeVectorPairMagSq
(const Oxs_ThreeVector& vec0,
 const Oxs_ThreeVector& vec1,
 OC_REAL8m& magsq0, OC_REAL8m& magsq1)
{ // Returns |vec0|^2 in mag0, and |vec1|^2 in mag1
#if !OC_USE_SSE
  magsq0 = vec0.x*vec0.x + vec0.y*vec0.y + vec0.z*vec0.z;
  magsq1 = vec1.x*vec1.x + vec1.y*vec1.y + vec1.z*vec1.z;
#else
  __m128d tmagsq = Oxs_ThreeVectorPairMagSq(vec0,vec1);
  _mm_storel_pd(&magsq0,tmagsq);
  _mm_storeh_pd(&magsq1,tmagsq);
#endif
}

#if OC_USE_SSE
void Oxs_ThreeVectorPairMakeUnit
(__m128d& tx,__m128d& ty,__m128d& tz);
// Here imports/exports tx,ty,tz are packed double-precision
// components of two three vectors,
//  tx = (t1.x,t0.x), ty = (t1.y,t0.y), tz = (t1.z,t0.z)
// where t0 is placed in the lower-half of the paced words.
// Code for this routine is in the file threevector.cc.

void Oxs_ThreeVectorPairMakeUnit
(__m128d tx,__m128d ty,__m128d tz,
 Oxs_ThreeVector& vec0,
 Oxs_ThreeVector& vec1,
 OC_REAL8m* magsq0=0, OC_REAL8m* magsq1=0);
// Here imports tx,ty,tz are packed double-precision components of two
// three vectors,
//  tx = (t1.x,t0.x), ty = (t1.y,t0.y), tz = (t1.z,t0.z)
// where t0 is placed in the lower-half of the paced words.
// Code for this routine is in the file threevector.cc.
#endif

void Oxs_ThreeVectorPairMakeUnit
(Oxs_ThreeVector& vec0,
 Oxs_ThreeVector& vec1,
 OC_REAL8m* magsq0=0, OC_REAL8m* magsq1=0);
// Parameters vec0 and vec1 are both imports and exports.
// Code for this routine is in the file threevector.cc.

#endif // _OXS_THREEVECTOR
