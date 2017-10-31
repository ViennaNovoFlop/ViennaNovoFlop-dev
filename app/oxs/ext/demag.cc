/* FILE: demag.cc            -*-Mode: c++-*-
 *
 * Average H demag field across rectangular cells.  This is a modified
 * version of the simpledemag class, which uses symmetries in the
 * interaction coefficients to reduce memory usage.
 *
 * The formulae used are reduced forms of equations in A. J. Newell,
 * W. Williams, and D. J. Dunlop, "A Generalization of the Demagnetizing
 * Tensor for Nonuniform Magnetization," Journal of Geophysical Research
 * - Solid Earth 98, 9551-9555 (1993).
 *
 * This code uses the Oxs_FFT3v classes to perform direct FFTs of the
 * import magnetization ThreeVectors.  This Oxs_Demag class is a
 * drop-in replacement for an older Oxs_Demag class that used the
 * scalar Oxs_FFT class.  That older class has been renamed
 * Oxs_DemagOld, and is contained in the demagold.* files.
 *
 */

#include "demag.h"  // Includes definition of OOMMF_THREADS macro

////////////////// SINGLE-THREADED IMPLEMENTATION  ///////////////
#if !OOMMF_THREADS

#include <assert.h>
#include <string>

#include "nb.h"
#include "director.h"
#include "key.h"
#include "mesh.h"
#include "meshvalue.h"
#include "simstate.h"
#include "threevector.h"
#include "energy.h"             // Needed to make MSVC++ 5 happy

#include "rectangularmesh.h"
#include "demagcoef.h"

OC_USE_STRING;

/* End includes */

// Oxs_Ext registration support
OXS_EXT_REGISTER(Oxs_Demag);

#ifndef VERBOSE_DEBUG
# define VERBOSE_DEBUG 0
#endif

// Size of threevector.  This macro is defined for code legibility
// and maintenance; it should always be "3".
#define ODTV_VECSIZE 3

// Size of complex value, in real units.  This macro is defined for code
// legibility and maintenance; it should always be "2".
#define ODTV_COMPLEXSIZE 2

// Constructor
Oxs_Demag::Oxs_Demag(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Oxs_Energy(name,newdtr,argstr),
    rdimx(0),rdimy(0),rdimz(0),cdimx(0),cdimy(0),cdimz(0),
    adimx(0),adimy(0),adimz(0),
    mesh_id(0),
    A(0),Hxfrm(0),asymptotic_radius(-1),Mtemp(0),
    embed_convolution(0),embed_block_size(0)
{
  asymptotic_radius = GetRealInitValue("asymptotic_radius",32.0);
  /// Value of -1 disables use of asymptotic approximation.

  cache_size = 1024*GetIntInitValue("cache_size_KB",1024);
  /// Cache size in KB.  Default is 1 MB.  Code wants bytes, so multiply
  /// user input by 1024.  cache_size is used to set embed_block_size in
  /// FillCoefficientArrays member function.

  zero_self_demag = GetIntInitValue("zero_self_demag",0);
  /// If true, then diag(1/3,1/3,1/3) is subtracted from the self-demag
  /// term.  In particular, for cubic cells this makes the self-demag
  /// field zero.  This will change the value computed for the demag
  /// energy by a constant amount, but since the demag field is changed
  /// by a multiple of m, the torque and therefore the magnetization
  /// dynamics are unaffected.

  VerifyAllInitArgsUsed();
}

Oxs_Demag::~Oxs_Demag() {
#if REPORT_TIME
  Oc_TimeVal cpu,wall;

  inittime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ...   init%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  fftforwardtime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ...  f-fft%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }
  fftinversetime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ...  i-fft%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  fftxforwardtime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ... f-fftx%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }
  fftxinversetime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ... i-fftx%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  fftyforwardtime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ... f-ffty%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }
  fftyinversetime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ... i-ffty%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  convtime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ...   conv%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  dottime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ...    dot%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

#endif // REPORT_TIME
  ReleaseMemory();
}

OC_BOOL Oxs_Demag::Init()
{
#if REPORT_TIME
  Oc_TimeVal cpu,wall;

  inittime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ...   init%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  fftforwardtime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ...  f-fft%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }
  fftinversetime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ...  i-fft%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  fftxforwardtime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ... f-fftx%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }
  fftxinversetime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ... i-fftx%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  fftyforwardtime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ... f-ffty%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }
  fftyinversetime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ... i-ffty%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  convtime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ...   conv%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  dottime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      subtime ...    dot%7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  inittime.Reset();
  fftforwardtime.Reset();
  fftinversetime.Reset();
  fftxforwardtime.Reset();
  fftxinversetime.Reset();
  fftyforwardtime.Reset();
  fftyinversetime.Reset();
  convtime.Reset();
  dottime.Reset();
#endif // REPORT_TIME
  mesh_id = 0;
  ReleaseMemory();
  return Oxs_Energy::Init();
}

void Oxs_Demag::ReleaseMemory() const
{ // Conceptually const
  if(A!=0) { delete[] A; A=0; }
  if(Hxfrm!=0)       { delete[] Hxfrm;       Hxfrm=0;       }
  if(Mtemp!=0)       { delete[] Mtemp;       Mtemp=0;       }
  rdimx=rdimy=rdimz=0;
  cdimx=cdimy=cdimz=0;
  adimx=adimy=adimz=0;
}

void Oxs_Demag::FillCoefficientArrays(const Oxs_Mesh* genmesh) const
{ // This routine is conceptually const.

  const Oxs_RectangularMesh* mesh
    = dynamic_cast<const Oxs_RectangularMesh*>(genmesh);
  if(mesh==NULL) {
    String msg=String("Object ")
      + String(genmesh->InstanceName())
      + String(" is not a rectangular mesh.");
    throw Oxs_ExtError(this,msg);
  }

  // Clean-up from previous allocation, if any.
  ReleaseMemory();

#if REPORT_TIME
    inittime.Start();
#endif // REPORT_TIME
  // Fill dimension variables
  rdimx = mesh->DimX();
  rdimy = mesh->DimY();
  rdimz = mesh->DimZ();
  if(rdimx==0 || rdimy==0 || rdimz==0) return; // Empty mesh!

  // Initialize fft object.  If a dimension equals 1, then zero
  // padding is not required.  Otherwise, zero pad to at least
  // twice the dimension.
  Oxs_FFT3DThreeVector::RecommendDimensions((rdimx==1 ? 1 : 2*rdimx),
                                            (rdimy==1 ? 1 : 2*rdimy),
                                            (rdimz==1 ? 1 : 2*rdimz),
                                            cdimx,cdimy,cdimz);
  OC_INDEX xfrm_size = ODTV_VECSIZE * 2 * cdimx * cdimy * cdimz;
  // "ODTV_VECSIZE" here is because we work with arrays if ThreeVectors,
  // and "2" because these are complex (as opposed to real)
  // quantities.
  if(xfrm_size<cdimx || xfrm_size<cdimy || xfrm_size<cdimz ||
     long(xfrm_size) != 
     long(2*ODTV_VECSIZE)*long(cdimx)*long(cdimy)*long(cdimz)) {
    // Partial overflow check
    char msgbuf[1024];
    Oc_Snprintf(msgbuf,sizeof(msgbuf),
                ": Product 2*ODTV_VECSIZE*cdimx*cdimy*cdimz = "
                "2*%d*%d*%d*%d too big to fit in a OC_INDEX variable",
                ODTV_VECSIZE,cdimx,cdimy,cdimz);
    String msg =
      String("OC_INDEX overflow in ")
      + String(InstanceName())
      + String(msgbuf);
    throw Oxs_ExtError(this,msg);
  }

  Mtemp = new OXS_FFT_REAL_TYPE[ODTV_VECSIZE*rdimx*rdimy*rdimz];
  /// Temporary space to hold Ms[]*m[].  The plan is to make this space
  /// unnecessary by introducing FFT functions that can take Ms as input
  /// and do the multiplication on the fly.


  // The following 3 statements are cribbed from
  // Oxs_FFT3DThreeVector::SetDimensions().  The corresponding
  // code using that class is
  //
  //  Oxs_FFT3DThreeVector fft;
  //  fft.SetDimensions(rdimx,rdimy,rdimz,cdimx,cdimy,cdimz);
  //  fft.GetLogicalDimensions(ldimx,ldimy,ldimz);
  //
  fftx.SetDimensions(rdimx, (cdimx==1 ? 1 : 2*(cdimx-1)), rdimy);
  ffty.SetDimensions(rdimy, cdimy,
                     ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx,
                     ODTV_VECSIZE*cdimx);
  fftz.SetDimensions(rdimz, cdimz,
                     ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx*cdimy,
                     ODTV_VECSIZE*cdimx*cdimy);
  OC_INDEX ldimx,ldimy,ldimz; // Logical dimensions
  // The following 3 statements are cribbed from
  // Oxs_FFT3DThreeVector::GetLogicalDimensions()
  ldimx = fftx.GetLogicalDimension();
  ldimy = ffty.GetLogicalDimension();
  ldimz = fftz.GetLogicalDimension();

  adimx = (ldimx/2)+1;
  adimy = (ldimy/2)+1;
  adimz = (ldimz/2)+1;

#if VERBOSE_DEBUG && !defined(NDEBUG)
  fprintf(stderr,"RDIMS: (%d,%d,%d)\n",rdimx,rdimy,rdimz); /**/
  fprintf(stderr,"CDIMS: (%d,%d,%d)\n",cdimx,cdimy,cdimz); /**/
  fprintf(stderr,"LDIMS: (%d,%d,%d)\n",ldimx,ldimy,ldimz); /**/
  fprintf(stderr,"ADIMS: (%d,%d,%d)\n",adimx,adimy,adimz); /**/
#endif // NDEBUG

  // Dimension of array necessary to hold 3 sets of full interaction
  // coefficients in real space:
  OC_INDEX scratch_size = ODTV_VECSIZE * ldimx * ldimy * ldimz;
  if(scratch_size<ldimx || scratch_size<ldimy || scratch_size<ldimz) {
    // Partial overflow check
    String msg =
      String("OC_INDEX overflow in ")
      + String(InstanceName())
      + String(": Product 3*8*rdimx*rdimy*rdimz too big"
               " to fit in a OC_INDEX variable");
    throw Oxs_ExtError(this,msg);
  }

  // Allocate memory for FFT xfrm target H, and scratch space
  // for computing interaction coefficients
  Hxfrm = new OXS_FFT_REAL_TYPE[xfrm_size];
  OXS_FFT_REAL_TYPE* scratch = new OXS_FFT_REAL_TYPE[scratch_size];

  if(Hxfrm==NULL || scratch==NULL) {
    // Safety check for those machines on which new[] doesn't throw
    // BadAlloc.
    String msg = String("Insufficient memory in Demag setup.");
    throw Oxs_ExtError(this,msg);
  }

  // According (16) in Newell's paper, the demag field is given by
  //                        H = -N*M
  // where N is the "demagnetizing tensor," with components Nxx, Nxy,
  // etc.  With the '-1' in 'scale' we store '-N' instead of 'N',
  // so we don't have to multiply the output from the FFT + iFFT
  // by -1 GetEnergy() below.

  // Fill interaction matrices with demag coefs from Newell's paper.
  // Note that A00, A11 and A22 are even in x,y and z.
  // A01 is odd in x and y, even in z.
  // A02 is odd in x and z, even in y.
  // A12 is odd in y and z, even in x.
  // We use these symmetries to reduce storage requirements.  If
  // f is real and even, then f^ is also real and even.  If f
  // is real and odd, then f^ is (pure) imaginary and odd.
  // As a result, the transform of each of the A## interaction
  // matrices will be real, with the same even/odd properties.
  //
  // Notation:  A00:=fs*Nxx, A01:=fs*Nxy, A02:=fs*Nxz,
  //                         A11:=fs*Nyy, A12:=fs*Nyz
  //                                      A22:=fs*Nzz
  //  where fs = -1/((ldimx/2)*ldimy*ldimz)

  OC_REALWIDE dx = mesh->EdgeLengthX();
  OC_REALWIDE dy = mesh->EdgeLengthY();
  OC_REALWIDE dz = mesh->EdgeLengthZ();
  // For demag calculation, all that matters is the relative
  // size of dx, dy and dz.  To help insure we don't run
  // outside floating point range, rescale these values so
  // largest is 1.0
  OC_REALWIDE maxedge=dx;
  if(dy>maxedge) maxedge=dy;
  if(dz>maxedge) maxedge=dz;
  dx/=maxedge; dy/=maxedge; dz/=maxedge;

  OC_REALWIDE scale = 1.0/(4*WIDE_PI*dx*dy*dz);

  // Also throw in FFT scaling.  This allows the "NoScale" FFT routines
  // to be used.  NB: There is effectively a "-1" built into the
  // differencing sections below, because we compute d^6/dx^2 dy^2 dz^2
  // instead of -d^6/dx^2 dy^2 dz^2 as required.
  // Note: Using an Oxs_FFT3DThreeVector fft object, this would be just
  //    scale *= fft.GetScaling();
  scale *= fftx.GetScaling() * ffty.GetScaling() * fftz.GetScaling();

  // Calculate Nxx, Nxy and Nxz in first octant.
  // Step 1: Evaluate f & g at each cell site.  Offset by (-dx,-dy,-dz)
  //  so we can do 2nd derivative operations "in-place".
  OC_INDEX i,j,k;
  OC_INDEX kstop=1; if(rdimz>1) kstop=rdimz+2;
  OC_INDEX jstop=1; if(rdimy>1) jstop=rdimy+2;
  OC_INDEX istop=1; if(rdimx>1) istop=rdimx+2;
  OC_INDEX sstridey = ODTV_VECSIZE*ldimx;
  OC_INDEX sstridez = sstridey*ldimy;
  for(k=0;k<kstop;k++) {
    OC_INDEX kindex = k*sstridez;
    OC_REALWIDE z = dz*(k-1);
    for(j=0;j<jstop;j++) {
      OC_INDEX jkindex = kindex + j*sstridey;
      OC_REALWIDE y = dy*(j-1);
      for(i=0;i<istop;i++) {
        OC_INDEX index = ODTV_VECSIZE*i+jkindex;
#ifndef NDEBUG
        if(index>=scratch_size) {
          String msg = String("Programming error:"
                              " array index out-of-bounds.");
          throw Oxs_ExtError(this,msg);
        }
#endif // NDEBUG
        OC_REALWIDE x = dx*(i-1);
        // Nyy(x,y,z) = Nxx(y,x,z);  Nzz(x,y,z) = Nxx(z,y,x);
        // Nxz(x,y,z) = Nxy(x,z,y);  Nyz(x,y,z) = Nxy(y,z,x);
        scratch[index]   = scale*Newell_f(x,y,z);  // For Nxx
        scratch[index+1] = scale*Newell_g(x,y,z);  // For Nxy
        scratch[index+2] = scale*Newell_g(x,z,y);  // For Nxz
      }
    }
  }

  // Step 2a: Do d^2/dz^2
  if(kstop==1) {
    // Only 1 layer in z-direction of f/g stored in scratch array.
    for(j=0;j<jstop;j++) {
      OC_INDEX jkindex = j*sstridey;
      OC_REALWIDE y = dy*(j-1);
      for(i=0;i<istop;i++) {
        OC_INDEX index = ODTV_VECSIZE*i+jkindex;
        OC_REALWIDE x = dx*(i-1);
        // Function f is even in each variable, so for example
        //    f(x,y,-dz) - 2f(x,y,0) + f(x,y,dz)
        //        =  2( f(x,y,-dz) - f(x,y,0) )
        // Function g(x,y,z) is even in z and odd in x and y,
        // so for example
        //    g(x,-dz,y) - 2g(x,0,y) + g(x,dz,y)
        //        =  2g(x,0,y) = 0.
        // Nyy(x,y,z) = Nxx(y,x,z);  Nzz(x,y,z) = Nxx(z,y,x);
        // Nxz(x,y,z) = Nxy(x,z,y);  Nyz(x,y,z) = Nxy(y,z,x);
        scratch[index]   -= scale*Newell_f(x,y,0);
        scratch[index]   *= 2;
        scratch[index+1] -= scale*Newell_g(x,y,0);
        scratch[index+1] *= 2;
        scratch[index+2] = 0;
      }
    }
  } else {
    for(k=0;k<rdimz;k++) {
      OC_INDEX kindex = k*sstridez;
      for(j=0;j<jstop;j++) {
        OC_INDEX jkindex = kindex + j*sstridey;
        for(i=0;i<istop;i++) {
          OC_INDEX index = ODTV_VECSIZE*i+jkindex;
          scratch[index]   += -2*scratch[index+sstridez]
                               + scratch[index+2*sstridez];
          scratch[index+1] += -2*scratch[index+sstridez+1]
                               + scratch[index+2*sstridez+1];
          scratch[index+2] += -2*scratch[index+sstridez+2]
                               + scratch[index+2*sstridez+2];
        }
      }
    }
    for(k=rdimz;k<ldimz;k++) { // Zero-fill overhang.  This fills
      /// into region that will be filled with reflected data,
      /// and so the span could be reduced.
      OC_INDEX kindex = k*sstridez;
      for(j=0;j<ldimy;j++) {
        OC_INDEX jkindex = kindex + j*sstridey;
        for(i=0;i<ODTV_VECSIZE*ldimx;i++) {
          scratch[jkindex+i] = 0.0;
        }
      }
    }
  }
  // Step 2b: Do d^2/dy^2
  if(jstop==1) {
    // Only 1 layer in y-direction of f/g stored in scratch array.
    for(k=0;k<rdimz;k++) {
      OC_INDEX kindex = k*sstridez;
      OC_REALWIDE z = dz*k;
      for(i=0;i<istop;i++) {
        OC_INDEX index = ODTV_VECSIZE*i+kindex;
        OC_REALWIDE x = dx*(i-1);
        // Function f is even in each variable, so for example
        //    f(x,y,-dz) - 2f(x,y,0) + f(x,y,dz)
        //        =  2( f(x,y,-dz) - f(x,y,0) )
        // Function g(x,y,z) is even in z and odd in x and y,
        // so for example
        //    g(x,-dz,y) - 2g(x,0,y) + g(x,dz,y)
        //        =  2g(x,0,y) = 0.
        // Nyy(x,y,z) = Nxx(y,x,z);  Nzz(x,y,z) = Nxx(z,y,x);
        // Nxz(x,y,z) = Nxy(x,z,y);  Nyz(x,y,z) = Nxy(y,z,x);
        scratch[index]   -= scale * 
          ((Newell_f(x,0,z-dz)+Newell_f(x,0,z+dz))-2*Newell_f(x,0,z));
        scratch[index]   *= 2;
        scratch[index+1]  = 0.0;
        scratch[index+2] -= scale * 
          ((Newell_g(x,z-dz,0)+Newell_g(x,z+dz,0))-2*Newell_g(x,z,0));
        scratch[index+2] *= 2;
      }
    }
  } else {
    for(k=0;k<rdimz;k++) {
      OC_INDEX kindex = k*sstridez;
      for(j=0;j<rdimy;j++) {
        OC_INDEX jkindex = kindex + j*sstridey;
        for(i=0;i<istop;i++) {
          OC_INDEX index = ODTV_VECSIZE*i+jkindex;
          scratch[index]   += -2*scratch[index+sstridey]
                               + scratch[index+2*sstridey];
          scratch[index+1] += -2*scratch[index+sstridey+1]
                               + scratch[index+2*sstridey+1];
          scratch[index+2] += -2*scratch[index+sstridey+2]
                               + scratch[index+2*sstridey+2];
        }
      }
    }
    for(k=0;k<rdimz;k++) { // Zero-fill overhang.  This fills
      /// into region that will be filled with reflected data,
      /// and so the span could be reduced.
      OC_INDEX kindex = k*sstridez;
      for(j=rdimy;j<ldimy;j++) {
        OC_INDEX jkindex = kindex + j*sstridey;
        for(i=0;i<ODTV_VECSIZE*ldimx;i++) {
          scratch[jkindex+i] = 0.0;
        }
      }
    }
  }

  // Step 2c: Do d^2/dx^2
  if(istop==1) {
    // Only 1 layer in x-direction of f/g stored in scratch array.
    for(k=0;k<rdimz;k++) {
      OC_INDEX kindex = k*sstridez;
      OC_REALWIDE z = dz*k;
      for(j=0;j<rdimy;j++) {
        OC_INDEX index = kindex + j*sstridey;
        OC_REALWIDE y = dy*j;
        // Function f is even in each variable, so for example
        //    f(x,y,-dz) - 2f(x,y,0) + f(x,y,dz)
        //        =  2( f(x,y,-dz) - f(x,y,0) )
        // Function g(x,y,z) is even in z and odd in x and y,
        // so for example
        //    g(x,-dz,y) - 2g(x,0,y) + g(x,dz,y)
        //        =  2g(x,0,y) = 0.
        // Nyy(x,y,z) = Nxx(y,x,z);  Nzz(x,y,z) = Nxx(z,y,x);
        // Nxz(x,y,z) = Nxy(x,z,y);  Nyz(x,y,z) = Nxy(y,z,x);
        scratch[index]   -= scale * 
          ((4*Newell_f(0,y,z)
            +Newell_f(0,y+dy,z+dz)+Newell_f(0,y-dy,z+dz)
            +Newell_f(0,y+dy,z-dz)+Newell_f(0,y-dy,z-dz))
           -2*(Newell_f(0,y+dy,z)+Newell_f(0,y-dy,z)
               +Newell_f(0,y,z+dz)+Newell_f(0,y,z-dz)));
        scratch[index]   *= 2;                       // For Nxx
        scratch[index+2]  = scratch[index+1] = 0.0;  // For Nxy & Nxz
      }
    }
  } else {
    for(k=0;k<rdimz;k++) {
      OC_INDEX kindex = k*sstridez;
      for(j=0;j<rdimy;j++) {
        OC_INDEX jkindex = kindex + j*sstridey;
        for(i=0;i<rdimx;i++) {
          OC_INDEX index = ODTV_VECSIZE*i+jkindex;
          scratch[index]   += -2*scratch[index+  ODTV_VECSIZE]
                               + scratch[index+2*ODTV_VECSIZE];
          scratch[index+1] += -2*scratch[index+  ODTV_VECSIZE+1]
                               + scratch[index+2*ODTV_VECSIZE+1];
          scratch[index+2] += -2*scratch[index+  ODTV_VECSIZE+2]
                               + scratch[index+2*ODTV_VECSIZE+2];
        }
      }
    }
    for(k=0;k<rdimz;k++) { // Zero-fill overhang.  This fills
      /// into region that will be filled with reflected data,
      /// and so the span could be reduced.
      OC_INDEX kindex = k*sstridez;
      for(j=0;j<rdimy;j++) {
        OC_INDEX jkindex = kindex + j*sstridey;
        for(i=rdimx;i<ldimx;i++) {
          OC_INDEX index = ODTV_VECSIZE*i+jkindex;
          scratch[index+2] = scratch[index+1] = scratch[index] = 0.0;
        }
      }
    }
  }

  // Step 2.5: Use asymptotic (dipolar + higher) approximation for far field
  /*   Dipole approximation:
   *
   *                        / 3x^2-R^2   3xy       3xz    \
   *             dx.dy.dz   |                             |
   *  H_demag = ----------- |   3xy   3y^2-R^2     3yz    |
   *             4.pi.R^5   |                             |
   *                        \   3xz      3yz     3z^2-R^2 /
   */
  // See Notes IV, 26-Feb-2007, p102.
  if(asymptotic_radius>=0.0) {
    // Note that all distances here are in "reduced" units,
    // scaled so that the largest of dx, dy, and dz = 1.0.
    OC_REALWIDE asymptotic_radius_sq = asymptotic_radius*asymptotic_radius;
    OC_REALWIDE fft_scaling = -1 *
      fftx.GetScaling() * ffty.GetScaling() * fftz.GetScaling();
    /// Note: Since H = -N*M, and by convention with the rest of this
    /// code, we store "-N" instead of "N" so we don't have to multiply
    /// the output from the FFT + iFFT by -1 in GetEnergy() below.

    OC_REALWIDE xtest = static_cast<OC_REALWIDE>(rdimx)*dx; xtest *= xtest;

    for(k=0;k<rdimz;++k) {
      OC_INDEX kindex = k*sstridez;
      OC_REALWIDE z = dz*k;
      OC_REALWIDE zsq = z*z;
      for(j=0;j<rdimy;++j) {
        OC_INDEX jkindex = kindex + j*sstridey;
        OC_REALWIDE y = dy*j;
        OC_REALWIDE ysq = y*y;

        OC_INDEX istart = 0;
        OC_REALWIDE test = asymptotic_radius_sq-ysq-zsq;
        if(test>0) {
          if(test>xtest) istart = rdimx+1;
          else           istart = static_cast<OC_INDEX>(ceil(sqrt(test)/dx));
        }

        for(i=istart;i<rdimx;++i) {
          OC_INDEX index = ODTV_VECSIZE*i+jkindex;
          OC_REALWIDE x = dx*i;
          scratch[index]   = fft_scaling*DemagNxxAsymptotic(x,y,z,dx,dy,dz);
          scratch[index+1] = fft_scaling*DemagNxyAsymptotic(x,y,z,dx,dy,dz);
          scratch[index+2] = fft_scaling*DemagNxzAsymptotic(x,y,z,dx,dy,dz);
        }
      }
    }
  }

  // Step 3: Use symmetries to reflect into other octants.
  //     Also, at each coordinate plane, set to 0.0 any term
  //     which is odd across that boundary.  It should already
  //     be close to 0, but will likely be slightly off due to
  //     rounding errors.
  // Symmetries: A00, A11, A22 are even in each coordinate
  //             A01 is odd in x and y, even in z.
  //             A02 is odd in x and z, even in y.
  //             A12 is odd in y and z, even in x.
  for(k=0;k<rdimz;k++) {
    OC_INDEX kindex = k*sstridez;
    for(j=0;j<rdimy;j++) {
      OC_INDEX jkindex = kindex + j*sstridey;
      for(i=0;i<rdimx;i++) {
        OC_INDEX index = ODTV_VECSIZE*i+jkindex;

        if(i==0 || j==0) scratch[index+1] = 0.0;  // A01
        if(i==0 || k==0) scratch[index+2] = 0.0;  // A02

        OC_REALWIDE tmpA00 = scratch[index];
        OC_REALWIDE tmpA01 = scratch[index+1];
        OC_REALWIDE tmpA02 = scratch[index+2];
        if(i>0) {
          OC_INDEX tindex = ODTV_VECSIZE*(ldimx-i)+j*sstridey+k*sstridez;
          scratch[tindex]   =     tmpA00;
          scratch[tindex+1] =  -1*tmpA01;
          scratch[tindex+2] =  -1*tmpA02;
        }
        if(j>0) {
          OC_INDEX tindex = ODTV_VECSIZE*i+(ldimy-j)*sstridey+k*sstridez;
          scratch[tindex]   =     tmpA00;
          scratch[tindex+1] =  -1*tmpA01;
          scratch[tindex+2] =     tmpA02;
        }
        if(k>0) {
          OC_INDEX tindex = ODTV_VECSIZE*i+j*sstridey+(ldimz-k)*sstridez;
          scratch[tindex]   =     tmpA00;
          scratch[tindex+1] =     tmpA01;
          scratch[tindex+2] =  -1*tmpA02;
        }
        if(i>0 && j>0) {
          OC_INDEX tindex
            = ODTV_VECSIZE*(ldimx-i) + (ldimy-j)*sstridey + k*sstridez;
          scratch[tindex]   =     tmpA00;
          scratch[tindex+1] =     tmpA01;
          scratch[tindex+2] =  -1*tmpA02;
        }
        if(i>0 && k>0) {
          OC_INDEX tindex
            = ODTV_VECSIZE*(ldimx-i) + j*sstridey + (ldimz-k)*sstridez;
          scratch[tindex]   =     tmpA00;
          scratch[tindex+1] =  -1*tmpA01;
          scratch[tindex+2] =     tmpA02;
        }
        if(j>0 && k>0) {
          OC_INDEX tindex
            = ODTV_VECSIZE*i + (ldimy-j)*sstridey + (ldimz-k)*sstridez;
          scratch[tindex]   =     tmpA00;
          scratch[tindex+1] =  -1*tmpA01;
          scratch[tindex+2] =  -1*tmpA02;
        }
        if(i>0 && j>0 && k>0) {
          OC_INDEX tindex
            = ODTV_VECSIZE*(ldimx-i) + (ldimy-j)*sstridey + (ldimz-k)*sstridez;
          scratch[tindex]   =     tmpA00;
          scratch[tindex+1] =     tmpA01;
          scratch[tindex+2] =     tmpA02;
        }
      }
    }
  }

  // Special "SelfDemag" code may be more accurate at index 0,0,0.
  // Note: Using an Oxs_FFT3DThreeVector fft object, this would be
  //    scale *= fft.GetScaling();
  OC_REALWIDE selfscale
    = -1 * fftx.GetScaling() * ffty.GetScaling() * fftz.GetScaling();
  scratch[0] = SelfDemagNx(dx,dy,dz);
  if(zero_self_demag) scratch[0] -= 1./3.;
  scratch[0] *= selfscale;

  scratch[1] = 0.0;  // Nxy[0] = 0.

  scratch[2] = 0.0;  // Nxz[0] = 0.

#if VERBOSE_DEBUG && !defined(NDEBUG)
  for(k=0;k<ldimz;++k) {
    for(j=0;j<ldimy;++j) {
      for(i=0;i<ldimx;++i) {
        OC_INDEX index = ODTV_VECSIZE*((k*ldimy+j)*ldimx+i);
        printf("A00[%02d][%02d][%02d] = %#25.12f\n",
               i,j,k,0.5*scratch[index]);
        printf("A01[%02d][%02d][%02d] = %#25.12f\n",
               i,j,k,0.5*scratch[index+1]);
        printf("A02[%02d][%02d][%02d] = %#25.12f\n",
               i,j,k,0.5*scratch[index+2]);
      }
    }
  }
  fflush(stdout);
#endif // NDEBUG

  // Step 4: Transform into frequency domain.  These lines are cribbed
  // from the corresponding code in Oxs_FFT3DThreeVector.
  // Note: Using an Oxs_FFT3DThreeVector fft object, this would be just
  //    fft.AdjustInputDimensions(ldimx,ldimy,ldimz);
  //    fft.ForwardRealToComplexFFT(scratch,Hxfrm);
  //    fft.AdjustInputDimensions(rdimx,rdimy,rdimz); // Safety
  {
    fftx.AdjustInputDimensions(ldimx,ldimy);
    ffty.AdjustInputDimensions(ldimy,
                               ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx,
                               ODTV_VECSIZE*cdimx);
    fftz.AdjustInputDimensions(ldimz,
                               ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx*cdimy,
                               ODTV_VECSIZE*cdimx*cdimy);

    OC_INDEX rxydim = ODTV_VECSIZE*ldimx*ldimy;
    OC_INDEX cxydim = ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx*cdimy;

    for(OC_INDEX m=0;m<ldimz;++m) {
      // x-direction transforms in plane "m"
      fftx.ForwardRealToComplexFFT(scratch+m*rxydim,Hxfrm+m*cxydim);
      
      // y-direction transforms in plane "m"
      ffty.ForwardFFT(Hxfrm+m*cxydim);
    }
    fftz.ForwardFFT(Hxfrm); // z-direction transforms

    fftx.AdjustInputDimensions(rdimx,rdimy);   // Safety
    ffty.AdjustInputDimensions(rdimy,
                               ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx,
                               ODTV_VECSIZE*cdimx);
    fftz.AdjustInputDimensions(rdimz,
                               ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx*cdimy,
                               ODTV_VECSIZE*cdimx*cdimy);

  }

  // Copy results from scratch into A00, A01, and A02.  We only need
  // store 1/8th of the results because of symmetries.
  OC_INDEX astridey = adimx;
  OC_INDEX astridez = astridey*adimy;
  OC_INDEX a_size = astridez*adimz;
  A = new A_coefs[a_size];

  OC_INDEX cstridey = 2*ODTV_VECSIZE*cdimx; // "2" for complex data
  OC_INDEX cstridez = cstridey*cdimy;
  for(k=0;k<adimz;k++) for(j=0;j<adimy;j++) for(i=0;i<adimx;i++) {
    OC_INDEX aindex = i+j*astridey+k*astridez;
    OC_INDEX hindex = 2*ODTV_VECSIZE*i+j*cstridey+k*cstridez;
    A[aindex].A00 = Hxfrm[hindex];   // A00
    A[aindex].A01 = Hxfrm[hindex+2]; // A01
    A[aindex].A02 = Hxfrm[hindex+4]; // A02
    // The A## values are all real-valued, so we only need to pull the
    // real parts out of Hxfrm, which are stored in the even offsets.
  }

  // Repeat for Nyy, Nyz and Nzz. //////////////////////////////////////
  // Step 1: Evaluate f & g at each cell site.  Offset by (-dx,-dy,-dz)
  //  so we can do 2nd derivative operations "in-place".
  for(k=0;k<kstop;k++) {
    OC_INDEX kindex = k*sstridez;
    OC_REALWIDE z = dz*(k-1);
    for(j=0;j<jstop;j++) {
      OC_INDEX jkindex = kindex + j*sstridey;
      OC_REALWIDE y = dy*(j-1);
      for(i=0;i<istop;i++) {
        OC_INDEX index = ODTV_VECSIZE*i+jkindex;
        OC_REALWIDE x = dx*(i-1);
        // Nyy(x,y,z) = Nxx(y,x,z);  Nzz(x,y,z) = Nxx(z,y,x);
        // Nxz(x,y,z) = Nxy(x,z,y);  Nyz(x,y,z) = Nxy(y,z,x);
        scratch[index]   = scale*Newell_f(y,x,z);  // For Nyy
        scratch[index+1] = scale*Newell_g(y,z,x);  // For Nyz
        scratch[index+2] = scale*Newell_f(z,y,x);  // For Nzz
      }
    }
  }

  // Step 2a: Do d^2/dz^2
  if(kstop==1) {
    // Only 1 layer in z-direction of f/g stored in scratch array.
    for(j=0;j<jstop;j++) {
      OC_INDEX jkindex = j*sstridey;
      OC_REALWIDE y = dy*(j-1);
      for(i=0;i<istop;i++) {
        OC_INDEX index = ODTV_VECSIZE*i+jkindex;
        OC_REALWIDE x = dx*(i-1);
        // Function f is even in each variable, so for example
        //    f(x,y,-dz) - 2f(x,y,0) + f(x,y,dz)
        //        =  2( f(x,y,-dz) - f(x,y,0) )
        // Function g(x,y,z) is even in z and odd in x and y,
        // so for example
        //    g(x,-dz,y) - 2g(x,0,y) + g(x,dz,y)
        //        =  2g(x,0,y) = 0.
        // Nyy(x,y,z) = Nxx(y,x,z);  Nzz(x,y,z) = Nxx(z,y,x);
        // Nxz(x,y,z) = Nxy(x,z,y);  Nyz(x,y,z) = Nxy(y,z,x);
        scratch[index]   -= scale*Newell_f(y,x,0);  // For Nyy
        scratch[index]   *= 2;
        scratch[index+1]  = 0.0;                    // For Nyz
        scratch[index+2] -= scale*Newell_f(0,y,x);  // For Nzz
        scratch[index+2] *= 2;
      }
    }
  } else {
    for(k=0;k<rdimz;k++) {
      OC_INDEX kindex = k*sstridez;
      for(j=0;j<jstop;j++) {
        OC_INDEX jkindex = kindex + j*sstridey;
        for(i=0;i<istop;i++) {
          OC_INDEX index = ODTV_VECSIZE*i+jkindex;
          scratch[index]   += -2*scratch[index+sstridez]
                               + scratch[index+2*sstridez];
          scratch[index+1] += -2*scratch[index+sstridez+1]
                               + scratch[index+2*sstridez+1];
          scratch[index+2] += -2*scratch[index+sstridez+2]
                               + scratch[index+2*sstridez+2];
        }
      }
    }
    for(k=rdimz;k<ldimz;k++) { // Zero-fill overhang.  This fills
      /// into region that will be filled with reflected data,
      /// and so the span could be reduced.
      OC_INDEX kindex = k*sstridez;
      for(j=0;j<ldimy;j++) {
        OC_INDEX jkindex = kindex + j*sstridey;
        for(i=0;i<ODTV_VECSIZE*ldimx;i++) {
          scratch[jkindex+i] = 0.0;
        }
      }
    }
  }
  // Step 2b: Do d^2/dy^2
  if(jstop==1) {
    // Only 1 layer in y-direction of f/g stored in scratch array.
    for(k=0;k<rdimz;k++) {
      OC_INDEX kindex = k*sstridez;
      OC_REALWIDE z = dz*k;
      for(i=0;i<istop;i++) {
        OC_INDEX index = ODTV_VECSIZE*i+kindex;
        OC_REALWIDE x = dx*(i-1);
        // Function f is even in each variable, so for example
        //    f(x,y,-dz) - 2f(x,y,0) + f(x,y,dz)
        //        =  2( f(x,y,-dz) - f(x,y,0) )
        // Function g(x,y,z) is even in z and odd in x and y,
        // so for example
        //    g(x,-dz,y) - 2g(x,0,y) + g(x,dz,y)
        //        =  2g(x,0,y) = 0.
        // Nyy(x,y,z) = Nxx(y,x,z);  Nzz(x,y,z) = Nxx(z,y,x);
        // Nxz(x,y,z) = Nxy(x,z,y);  Nyz(x,y,z) = Nxy(y,z,x);
        scratch[index]   -= scale * 
          ((Newell_f(0,x,z-dz)+Newell_f(0,x,z+dz))-2*Newell_f(0,x,z));
        scratch[index]   *= 2;   // For Nyy
        scratch[index+1] = 0.0;  // For Nyz
        scratch[index+2] -= scale * 
          ((Newell_f(z-dz,0,x)+Newell_f(z+dz,0,x))-2*Newell_f(z,0,x));
        scratch[index+2] *= 2;   // For Nzz
      }
    }
  } else {
    for(k=0;k<rdimz;k++) {
      OC_INDEX kindex = k*sstridez;
      for(j=0;j<rdimy;j++) {
        OC_INDEX jkindex = kindex + j*sstridey;
        for(i=0;i<istop;i++) {
          OC_INDEX index = ODTV_VECSIZE*i+jkindex;
          scratch[index]   += -2*scratch[index+sstridey]
                               + scratch[index+2*sstridey];
          scratch[index+1] += -2*scratch[index+sstridey+1]
                               + scratch[index+2*sstridey+1];
          scratch[index+2] += -2*scratch[index+sstridey+2]
                               + scratch[index+2*sstridey+2];
        }
      }
    }
    for(k=0;k<rdimz;k++) { // Zero-fill overhang.  This fills
      /// into region that will be filled with reflected data,
      /// and so the span could be reduced.
      OC_INDEX kindex = k*sstridez;
      for(j=rdimy;j<ldimy;j++) {
        OC_INDEX jkindex = kindex + j*sstridey;
        for(i=0;i<ODTV_VECSIZE*ldimx;i++) {
          scratch[jkindex+i] = 0.0;
        }
      }
    }
  }
  // Step 2c: Do d^2/dx^2
  if(istop==1) {
    // Only 1 layer in x-direction of f/g stored in scratch array.
    for(k=0;k<rdimz;k++) {
      OC_INDEX kindex = k*sstridez;
      OC_REALWIDE z = dz*k;
      for(j=0;j<rdimy;j++) {
        OC_INDEX index = kindex + j*sstridey;
        OC_REALWIDE y = dy*j;
        // Function f is even in each variable, so for example
        //    f(x,y,-dz) - 2f(x,y,0) + f(x,y,dz)
        //        =  2( f(x,y,-dz) - f(x,y,0) )
        // Function g(x,y,z) is even in z and odd in x and y,
        // so for example
        //    g(x,-dz,y) - 2g(x,0,y) + g(x,dz,y)
        //        =  2g(x,0,y) = 0.
        // Nyy(x,y,z) = Nxx(y,x,z);  Nzz(x,y,z) = Nxx(z,y,x);
        // Nxz(x,y,z) = Nxy(x,z,y);  Nyz(x,y,z) = Nxy(y,z,x);
        scratch[index]   -= scale * 
          ((4*Newell_f(y,0,z)
            +Newell_f(y+dy,0,z+dz)+Newell_f(y-dy,0,z+dz)
            +Newell_f(y+dy,0,z-dz)+Newell_f(y-dy,0,z-dz))
           -2*(Newell_f(y+dy,0,z)+Newell_f(y-dy,0,z)
              +Newell_f(y,0,z+dz)+Newell_f(y,0,z-dz)));
        scratch[index]   *= 2;  // For Nyy
        scratch[index+1] -= scale * 
          ((4*Newell_g(y,z,0)
            +Newell_g(y+dy,z+dz,0)+Newell_g(y-dy,z+dz,0)
            +Newell_g(y+dy,z-dz,0)+Newell_g(y-dy,z-dz,0))
           -2*(Newell_g(y+dy,z,0)+Newell_g(y-dy,z,0)
               +Newell_g(y,z+dz,0)+Newell_g(y,z-dz,0)));
        scratch[index+1] *= 2;  // For Nyz
        scratch[index+2] -= scale * 
          ((4*Newell_f(z,y,0)
            +Newell_f(z+dz,y+dy,0)+Newell_f(z+dz,y-dy,0)
            +Newell_f(z-dz,y+dy,0)+Newell_f(z-dz,y-dy,0))
           -2*(Newell_f(z,y+dy,0)+Newell_f(z,y-dy,0)
              +Newell_f(z+dz,y,0)+Newell_f(z-dz,y,0)));
        scratch[index+2] *= 2;  // For Nzz
      }
    }
  } else {
    for(k=0;k<rdimz;k++) {
      OC_INDEX kindex = k*sstridez;
      for(j=0;j<rdimy;j++) {
        OC_INDEX jkindex = kindex + j*sstridey;
        for(i=0;i<rdimx;i++) {
          OC_INDEX index = ODTV_VECSIZE*i+jkindex;
          scratch[index]   += -2*scratch[index+  ODTV_VECSIZE]
                               + scratch[index+2*ODTV_VECSIZE];
          scratch[index+1] += -2*scratch[index+  ODTV_VECSIZE+1]
                               + scratch[index+2*ODTV_VECSIZE+1];
          scratch[index+2] += -2*scratch[index+  ODTV_VECSIZE+2]
                               + scratch[index+2*ODTV_VECSIZE+2];
        }
      }
    }
    for(k=0;k<rdimz;k++) { // Zero-fill overhang.  This fills
      /// into region that will be filled with reflected data,
      /// and so the span could be reduced.
      OC_INDEX kindex = k*sstridez;
      for(j=0;j<rdimy;j++) {
        OC_INDEX jkindex = kindex + j*sstridey;
        for(i=rdimx;i<ldimx;i++) {
          OC_INDEX index = ODTV_VECSIZE*i+jkindex;
          scratch[index+2] = scratch[index+1] = scratch[index] = 0.0;
        }
      }
    }
  }

  // Step 2.5: Use asymptotic (dipolar + higher) approximation for far field
  /*   Dipole approximation:
   *
   *                        / 3x^2-R^2   3xy       3xz    \
   *             dx.dy.dz   |                             |
   *  H_demag = ----------- |   3xy   3y^2-R^2     3yz    |
   *             4.pi.R^5   |                             |
   *                        \   3xz      3yz     3z^2-R^2 /
   */
  // See Notes IV, 26-Feb-2007, p102.
  if(asymptotic_radius>=0.0) {
    // Note that all distances here are in "reduced" units,
    // scaled so that the largest of dx, dy, and dz = 1.0.
    OC_REALWIDE asymptotic_radius_sq = asymptotic_radius*asymptotic_radius;
    OC_REALWIDE fft_scaling = -1 *
      fftx.GetScaling() * ffty.GetScaling() * fftz.GetScaling();
    /// Note: Since H = -N*M, and by convention with the rest of this
    /// code, we store "-N" instead of "N" so we don't have to multiply
    /// the output from the FFT + iFFT by -1 in GetEnergy() below.

    OC_REALWIDE xtest = static_cast<OC_REALWIDE>(rdimx)*dx; xtest *= xtest;

    for(k=0;k<rdimz;++k) {
      OC_INDEX kindex = k*sstridez;
      OC_REALWIDE z = dz*k;
      OC_REALWIDE zsq = z*z;
      for(j=0;j<rdimy;++j) {
        OC_INDEX jkindex = kindex + j*sstridey;
        OC_REALWIDE y = dy*j;
        OC_REALWIDE ysq = y*y;

        OC_INDEX istart = 0;
        OC_REALWIDE test = asymptotic_radius_sq-ysq-zsq;
        if(test>0) {
          if(test>xtest) istart = rdimx+1;
          else           istart = static_cast<OC_INDEX>(ceil(sqrt(test)/dx));
        }

        for(i=istart;i<rdimx;++i) {
          OC_INDEX index = ODTV_VECSIZE*i+jkindex;
          OC_REALWIDE x = dx*i;
          scratch[index]   = fft_scaling*DemagNyyAsymptotic(x,y,z,dx,dy,dz);
          scratch[index+1] = fft_scaling*DemagNyzAsymptotic(x,y,z,dx,dy,dz);
          scratch[index+2] = fft_scaling*DemagNzzAsymptotic(x,y,z,dx,dy,dz);
        }
      }
    }
  }

  // Step 3: Use symmetries to reflect into other octants.
  //     Also, at each coordinate plane, set to 0.0 any term
  //     which is odd across that boundary.  It should already
  //     be close to 0, but will likely be slightly off due to
  //     rounding errors.
  // Symmetries: A00, A11, A22 are even in each coordinate
  //             A01 is odd in x and y, even in z.
  //             A02 is odd in x and z, even in y.
  //             A12 is odd in y and z, even in x.
  for(k=0;k<rdimz;k++) {
    OC_INDEX kindex = k*sstridez;
    for(j=0;j<rdimy;j++) {
      OC_INDEX jkindex = kindex + j*sstridey;
      for(i=0;i<rdimx;i++) {
        OC_INDEX index = ODTV_VECSIZE*i+jkindex;

        if(j==0 || k==0) scratch[index+1] = 0.0;  // A12

        OC_REALWIDE tmpA11 = scratch[index];
        OC_REALWIDE tmpA12 = scratch[index+1];
        OC_REALWIDE tmpA22 = scratch[index+2];
        if(i>0) {
          OC_INDEX tindex = ODTV_VECSIZE*(ldimx-i)+j*sstridey+k*sstridez;
          scratch[tindex]   =     tmpA11;
          scratch[tindex+1] =     tmpA12;
          scratch[tindex+2] =     tmpA22;
        }
        if(j>0) {
          OC_INDEX tindex = ODTV_VECSIZE*i+(ldimy-j)*sstridey+k*sstridez;
          scratch[tindex]   =     tmpA11;
          scratch[tindex+1] =  -1*tmpA12;
          scratch[tindex+2] =     tmpA22;
        }
        if(k>0) {
          OC_INDEX tindex = ODTV_VECSIZE*i+j*sstridey+(ldimz-k)*sstridez;
          scratch[tindex]   =     tmpA11;
          scratch[tindex+1] =  -1*tmpA12;
          scratch[tindex+2] =     tmpA22;
        }
        if(i>0 && j>0) {
          OC_INDEX tindex
            = ODTV_VECSIZE*(ldimx-i) + (ldimy-j)*sstridey + k*sstridez;
          scratch[tindex]   =     tmpA11;
          scratch[tindex+1] =  -1*tmpA12;
          scratch[tindex+2] =     tmpA22;
        }
        if(i>0 && k>0) {
          OC_INDEX tindex
            = ODTV_VECSIZE*(ldimx-i) + j*sstridey + (ldimz-k)*sstridez;
          scratch[tindex]   =     tmpA11;
          scratch[tindex+1] =  -1*tmpA12;
          scratch[tindex+2] =     tmpA22;
        }
        if(j>0 && k>0) {
          OC_INDEX tindex
            = ODTV_VECSIZE*i + (ldimy-j)*sstridey + (ldimz-k)*sstridez;
          scratch[tindex]   =     tmpA11;
          scratch[tindex+1] =     tmpA12;
          scratch[tindex+2] =     tmpA22;
        }
        if(i>0 && j>0 && k>0) {
          OC_INDEX tindex
            = ODTV_VECSIZE*(ldimx-i) + (ldimy-j)*sstridey + (ldimz-k)*sstridez;
          scratch[tindex]   =     tmpA11;
          scratch[tindex+1] =     tmpA12;
          scratch[tindex+2] =     tmpA22;
        }
      }
    }
  }

  // Special "SelfDemag" code may be more accurate at index 0,0,0.
  scratch[0] = SelfDemagNy(dx,dy,dz);
  if(zero_self_demag) scratch[0] -= 1./3.;
  scratch[0] *= selfscale;

  scratch[1] = 0.0;  // Nyz[0] = 0.

  scratch[2] = SelfDemagNz(dx,dy,dz);
  if(zero_self_demag) scratch[2] -= 1./3.;
  scratch[2] *= selfscale;

#if VERBOSE_DEBUG && !defined(NDEBUG)
  for(k=0;k<ldimz;++k) {
    for(j=0;j<ldimy;++j) {
      for(i=0;i<ldimx;++i) {
        OC_INDEX index = ODTV_VECSIZE*((k*ldimy+j)*ldimx+i);
        printf("A11[%02d][%02d][%02d] = %#25.12f\n",
               i,j,k,0.5*scratch[index]);
        printf("A12[%02d][%02d][%02d] = %#25.12f\n",
               i,j,k,0.5*scratch[index+1]);
        printf("A22[%02d][%02d][%02d] = %#25.12f\n",
               i,j,k,0.5*scratch[index+2]);
      }
    }
  }
  fflush(stdout);
#endif // NDEBUG

  // Step 4: Transform into frequency domain.  These lines are cribbed
  // from the corresponding code in Oxs_FFT3DThreeVector.
  // Note: Using an Oxs_FFT3DThreeVector fft object, this would be just
  //    fft.AdjustInputDimensions(ldimx,ldimy,ldimz);
  //    fft.ForwardRealToComplexFFT(scratch,Hxfrm);
  //    fft.AdjustInputDimensions(rdimx,rdimy,rdimz); // Safety
  {
    fftx.AdjustInputDimensions(ldimx,ldimy);
    ffty.AdjustInputDimensions(ldimy,
                               ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx,
                               ODTV_VECSIZE*cdimx);
    fftz.AdjustInputDimensions(ldimz,
                               ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx*cdimy,
                               ODTV_VECSIZE*cdimx*cdimy);

    OC_INDEX rxydim = ODTV_VECSIZE*ldimx*ldimy;
    OC_INDEX cxydim = ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx*cdimy;

    for(OC_INDEX m=0;m<ldimz;++m) {
      // x-direction transforms in plane "m"
      fftx.ForwardRealToComplexFFT(scratch+m*rxydim,Hxfrm+m*cxydim);
      
      // y-direction transforms in plane "m"
      ffty.ForwardFFT(Hxfrm+m*cxydim);
    }
    fftz.ForwardFFT(Hxfrm); // z-direction transforms

    fftx.AdjustInputDimensions(rdimx,rdimy);   // Safety
    ffty.AdjustInputDimensions(rdimy,
                               ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx,
                               ODTV_VECSIZE*cdimx);
    fftz.AdjustInputDimensions(rdimz,
                               ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx*cdimy,
                               ODTV_VECSIZE*cdimx*cdimy);

  }

  // At this point we no longer need the "scratch" array, so release it.
  delete[] scratch;

  // Copy results from scratch into A11, A12, and A22.  We only need
  // store 1/8th of the results because of symmetries.
  for(k=0;k<adimz;k++) for(j=0;j<adimy;j++) for(i=0;i<adimx;i++) {
    OC_INDEX aindex =   i+j*astridey+k*astridez;
    OC_INDEX hindex = 2*ODTV_VECSIZE*i+j*cstridey+k*cstridez;
    A[aindex].A11 = Hxfrm[hindex];   // A11
    A[aindex].A12 = Hxfrm[hindex+2]; // A12
    A[aindex].A22 = Hxfrm[hindex+4]; // A22
    // The A## values are all real-valued, so we only need to pull the
    // real parts out of Hxfrm, which are stored in the even offsets.
  }


  // Do we want to embed "convolution" computation inside z-axis FFTs?
  // If so, setup control variables.
  OC_INDEX footprint
    = ODTV_COMPLEXSIZE*ODTV_VECSIZE*sizeof(OXS_FFT_REAL_TYPE) // Data
    + sizeof(A_coefs)                           // Interaction matrix
    + 2*ODTV_COMPLEXSIZE*sizeof(OXS_FFT_REAL_TYPE); // Roots of unity
  footprint *= cdimz;
  OC_INDEX trialsize = cache_size/(2*footprint); // "2" is fudge factor
  if(trialsize>cdimx) trialsize=cdimx;
  if(cdimz>1 && trialsize>4) {
    // Note: If cdimz==1, then the z-axis FFT is a nop, so there is
    // nothing to embed the "convolution" with and we are better off
    // using the non-embedded code.
    embed_convolution = 1;
    embed_block_size = trialsize;
  } else {
    embed_convolution = 0;
    embed_block_size = 0;  // A cry for help...
  }

#if REPORT_TIME
    inittime.Stop();
#endif // REPORT_TIME
}

void Oxs_Demag::GetEnergy
(const Oxs_SimState& state,
 Oxs_EnergyData& oed
 ) const
{
  OC_INDEX i,j,k;

  // (Re)-initialize mesh coefficient array if mesh has changed.
  if(mesh_id != state.mesh->Id()) {
    mesh_id = 0; // Safety
    FillCoefficientArrays(state.mesh);
    mesh_id = state.mesh->Id();
  }

  const Oxs_MeshValue<ThreeVector>& spin = state.spin;
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);

  // Use supplied buffer space, and reflect that use in oed.
  oed.energy = oed.energy_buffer;
  oed.field = oed.field_buffer;
  Oxs_MeshValue<OC_REAL8m>& energy = *oed.energy_buffer;
  Oxs_MeshValue<ThreeVector>& field = *oed.field_buffer;
  energy.AdjustSize(state.mesh);
  field.AdjustSize(state.mesh);

  // Fill Mtemp with Ms[]*spin[].  The plan is to eventually
  // roll this step into the forward FFT routine.
  const OC_INDEX rsize = Ms.Size();
  assert(rdimx*rdimy*rdimz == rsize);
  for(i=0;i<rsize;++i) {
    OC_REAL8m scale = Ms[i];
    const ThreeVector& vec = spin[i];
    Mtemp[3*i]   = scale*vec.x;
    Mtemp[3*i+1] = scale*vec.y;
    Mtemp[3*i+2] = scale*vec.z;
  }

  if(!embed_convolution) {
    // Do not embed convolution inside z-axis FFTs.  Instead,
    // first compute full forward FFT, then do the convolution
    // (really matrix-vector A^*M^ multiply), and then do the
    // full inverse FFT.
    
    // Calculate FFT of Mtemp
#if REPORT_TIME
    fftforwardtime.Start();
#endif // REPORT_TIME
    // Transform into frequency domain.  These lines are cribbed from the
    // corresponding code in Oxs_FFT3DThreeVector.
    // Note: Using an Oxs_FFT3DThreeVector object, this would be just
    //    fft.ForwardRealToComplexFFT(Mtemp,Hxfrm);
    {
      OC_INDEX rxydim = ODTV_VECSIZE*rdimx*rdimy;
      OC_INDEX cxydim = ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx*cdimy;
      for(OC_INDEX m=0;m<rdimz;++m) {
        // x-direction transforms in plane "m"
        fftx.ForwardRealToComplexFFT(Mtemp+m*rxydim,Hxfrm+m*cxydim);
        // y-direction transforms in plane "m"
        ffty.ForwardFFT(Hxfrm+m*cxydim);
      }
      fftz.ForwardFFT(Hxfrm); // z-direction transforms
    }
#if REPORT_TIME
    fftforwardtime.Stop();
#endif // REPORT_TIME

    // Calculate field components in frequency domain.  Make use of
    // realness and even/odd properties of interaction matrices Axx.
    // Note that in transform space only the x>=0 half-space is
    // stored.
    // Symmetries: A00, A11, A22 are even in each coordinate
    //             A01 is odd in x and y, even in z.
    //             A02 is odd in x and z, even in y.
    //             A12 is odd in y and z, even in x.
    assert(adimx>=cdimx);
    assert(cdimy-adimy<adimy);
    assert(cdimz-adimz<adimz);
#if REPORT_TIME
    convtime.Start();
#endif // REPORT_TIME
    const OC_INDEX  jstride = ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx;
    const OC_INDEX  kstride = jstride*cdimy;
    const OC_INDEX ajstride = adimx;
    const OC_INDEX akstride = ajstride*adimy;
    for(k=0;k<adimz;++k) {
      // k>=0
      OC_INDEX  kindex = k*kstride;
      OC_INDEX akindex = k*akstride;
      for(j=0;j<adimy;++j) {
        // j>=0, k>=0
        OC_INDEX  jindex =  kindex + j*jstride;
        OC_INDEX ajindex = akindex + j*ajstride;
        for(i=0;i<cdimx;++i) {
          OC_INDEX  index = ODTV_COMPLEXSIZE*ODTV_VECSIZE*i+jindex;
          OXS_FFT_REAL_TYPE Hx_re = Hxfrm[index];
          OXS_FFT_REAL_TYPE Hx_im = Hxfrm[index+1];
          OXS_FFT_REAL_TYPE Hy_re = Hxfrm[index+2];
          OXS_FFT_REAL_TYPE Hy_im = Hxfrm[index+3];
          OXS_FFT_REAL_TYPE Hz_re = Hxfrm[index+4];
          OXS_FFT_REAL_TYPE Hz_im = Hxfrm[index+5];

          const A_coefs& Aref = A[ajindex+i]; 

          Hxfrm[index]   = Aref.A00*Hx_re + Aref.A01*Hy_re + Aref.A02*Hz_re;
          Hxfrm[index+1] = Aref.A00*Hx_im + Aref.A01*Hy_im + Aref.A02*Hz_im;
          Hxfrm[index+2] = Aref.A01*Hx_re + Aref.A11*Hy_re + Aref.A12*Hz_re;
          Hxfrm[index+3] = Aref.A01*Hx_im + Aref.A11*Hy_im + Aref.A12*Hz_im;
          Hxfrm[index+4] = Aref.A02*Hx_re + Aref.A12*Hy_re + Aref.A22*Hz_re;
          Hxfrm[index+5] = Aref.A02*Hx_im + Aref.A12*Hy_im + Aref.A22*Hz_im;
        }
      }
      for(j=adimy;j<cdimy;++j) {
        // j<0, k>=0
        OC_INDEX  jindex =  kindex + j*jstride;
        OC_INDEX ajindex = akindex + (cdimy-j)*ajstride;
        for(i=0;i<cdimx;++i) {
          OC_INDEX  index = ODTV_COMPLEXSIZE*ODTV_VECSIZE*i+jindex;
          OXS_FFT_REAL_TYPE Hx_re = Hxfrm[index];
          OXS_FFT_REAL_TYPE Hx_im = Hxfrm[index+1];
          OXS_FFT_REAL_TYPE Hy_re = Hxfrm[index+2];
          OXS_FFT_REAL_TYPE Hy_im = Hxfrm[index+3];
          OXS_FFT_REAL_TYPE Hz_re = Hxfrm[index+4];
          OXS_FFT_REAL_TYPE Hz_im = Hxfrm[index+5];

          const A_coefs& Aref = A[ajindex+i]; 

          // Flip signs on a01 and a12 as compared to the j>=0
          // case because a01 and a12 are odd in y.
          Hxfrm[index]   =  Aref.A00*Hx_re - Aref.A01*Hy_re + Aref.A02*Hz_re;
          Hxfrm[index+1] =  Aref.A00*Hx_im - Aref.A01*Hy_im + Aref.A02*Hz_im;
          Hxfrm[index+2] = -Aref.A01*Hx_re + Aref.A11*Hy_re - Aref.A12*Hz_re;
          Hxfrm[index+3] = -Aref.A01*Hx_im + Aref.A11*Hy_im - Aref.A12*Hz_im;
          Hxfrm[index+4] =  Aref.A02*Hx_re - Aref.A12*Hy_re + Aref.A22*Hz_re;
          Hxfrm[index+5] =  Aref.A02*Hx_im - Aref.A12*Hy_im + Aref.A22*Hz_im;
        }
      }
    }
    for(k=adimz;k<cdimz;++k) {
      // k<0
      OC_INDEX  kindex = k*kstride;
      OC_INDEX akindex = (cdimz-k)*akstride;
      for(j=0;j<adimy;++j) {
        // j>=0, k<0
        OC_INDEX  jindex =  kindex + j*jstride;
        OC_INDEX ajindex = akindex + j*ajstride;
        for(i=0;i<cdimx;++i) {
          OC_INDEX  index = ODTV_COMPLEXSIZE*ODTV_VECSIZE*i+jindex;
          OXS_FFT_REAL_TYPE Hx_re = Hxfrm[index];
          OXS_FFT_REAL_TYPE Hx_im = Hxfrm[index+1];
          OXS_FFT_REAL_TYPE Hy_re = Hxfrm[index+2];
          OXS_FFT_REAL_TYPE Hy_im = Hxfrm[index+3];
          OXS_FFT_REAL_TYPE Hz_re = Hxfrm[index+4];
          OXS_FFT_REAL_TYPE Hz_im = Hxfrm[index+5];

          const A_coefs& Aref = A[ajindex+i]; 

          // Flip signs on a02 and a12 as compared to the k>=0, j>=0 case
          // because a02 and a12 are odd in z.
          Hxfrm[index]   =  Aref.A00*Hx_re + Aref.A01*Hy_re - Aref.A02*Hz_re;
          Hxfrm[index+1] =  Aref.A00*Hx_im + Aref.A01*Hy_im - Aref.A02*Hz_im;
          Hxfrm[index+2] =  Aref.A01*Hx_re + Aref.A11*Hy_re - Aref.A12*Hz_re;
          Hxfrm[index+3] =  Aref.A01*Hx_im + Aref.A11*Hy_im - Aref.A12*Hz_im;
          Hxfrm[index+4] = -Aref.A02*Hx_re - Aref.A12*Hy_re + Aref.A22*Hz_re;
          Hxfrm[index+5] = -Aref.A02*Hx_im - Aref.A12*Hy_im + Aref.A22*Hz_im;
        }
      }
      for(j=adimy;j<cdimy;++j) {
        // j<0, k<0
        OC_INDEX  jindex =  kindex + j*jstride;
        OC_INDEX ajindex = akindex + (cdimy-j)*ajstride;
        for(i=0;i<cdimx;++i) {
          OC_INDEX  index = ODTV_COMPLEXSIZE*ODTV_VECSIZE*i+jindex;
          OXS_FFT_REAL_TYPE Hx_re = Hxfrm[index];
          OXS_FFT_REAL_TYPE Hx_im = Hxfrm[index+1];
          OXS_FFT_REAL_TYPE Hy_re = Hxfrm[index+2];
          OXS_FFT_REAL_TYPE Hy_im = Hxfrm[index+3];
          OXS_FFT_REAL_TYPE Hz_re = Hxfrm[index+4];
          OXS_FFT_REAL_TYPE Hz_im = Hxfrm[index+5];

          const A_coefs& Aref = A[ajindex+i]; 

          // Flip signs on a01 and a02 as compared to the k>=0, j>=0 case
          // because a01 is odd in y and even in z,
          //     and a02 is odd in z and even in y.
          // No change to a12 because it is odd in both y and z.
          Hxfrm[index]   =  Aref.A00*Hx_re - Aref.A01*Hy_re - Aref.A02*Hz_re;
          Hxfrm[index+1] =  Aref.A00*Hx_im - Aref.A01*Hy_im - Aref.A02*Hz_im;
          Hxfrm[index+2] = -Aref.A01*Hx_re + Aref.A11*Hy_re + Aref.A12*Hz_re;
          Hxfrm[index+3] = -Aref.A01*Hx_im + Aref.A11*Hy_im + Aref.A12*Hz_im;
          Hxfrm[index+4] = -Aref.A02*Hx_re + Aref.A12*Hy_re + Aref.A22*Hz_re;
          Hxfrm[index+5] = -Aref.A02*Hx_im + Aref.A12*Hy_im + Aref.A22*Hz_im;
        }
      }
    }
#if REPORT_TIME
    convtime.Stop();
#endif // REPORT_TIME

#if REPORT_TIME
    fftinversetime.Start();
#endif // REPORT_TIME
    // Transform back into space domain.  These lines are cribbed from the
    // corresponding code in Oxs_FFT3DThreeVector.
    // Note: Using an Oxs_FFT3DThreeVector object, this would be
    //     assert(3*sizeof(OXS_FFT_REAL_TYPE)==sizeof(ThreeVector));
    //     void* fooptr = static_cast<void*>(&(field[0]));
    //     fft.InverseComplexToRealFFT(Hxfrm,
    //                static_cast<OXS_FFT_REAL_TYPE*>(fooptr));
    {
      OC_INDEX m;
      OC_INDEX rxydim = ODTV_VECSIZE*rdimx*rdimy;
      OC_INDEX cxydim = ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx*cdimy;
      assert(3*sizeof(OXS_FFT_REAL_TYPE)<=sizeof(ThreeVector));
      OXS_FFT_REAL_TYPE* fptr
        = static_cast<OXS_FFT_REAL_TYPE*>(static_cast<void*>(&field[OC_INDEX(0)]));
      fftz.InverseFFT(Hxfrm); // z-direction transforms
      for(m=0;m<rdimz;++m) {
        // y-direction transforms
        ffty.InverseFFT(Hxfrm+m*cxydim);
        // x-direction transforms
        fftx.InverseComplexToRealFFT(Hxfrm+m*cxydim,fptr+m*rxydim);
      }

      if(3*sizeof(OXS_FFT_REAL_TYPE)<sizeof(ThreeVector)) {
        // The fftx.InverseComplexToRealFFT calls above assume the
        // target is an array of OXS_FFT_REAL_TYPE.  If ThreeVector is
        // not tightly packed, then this assumption is false; however we
        // can correct the problem by expanding the results in-place.
        // The only setting I know of where ThreeVector doesn't tight
        // pack is under the Borland bcc32 compiler on Windows x86 with
        // OXS_FFT_REAL_TYPE equal to "long double".  In that case
        // sizeof(long double) == 10, but sizeof(ThreeVector) == 36.
        for(m = rsize - 1; m>=0 ; --m) {
          ThreeVector temp(fptr[ODTV_VECSIZE*m],fptr[ODTV_VECSIZE*m+1],fptr[ODTV_VECSIZE*m+2]);
          field[m] = temp;
        }
      }

    }
#if REPORT_TIME
    fftinversetime.Stop();
#endif // REPORT_TIME
  } else { // if(!embed_convolution)
    // Embed "convolution" (really matrix-vector multiply A^*M^) inside
    // z-axis FFTs.  First compute full forward x- and y-axis FFTs.
    // Then, do a small number of z-axis forward FFTs, followed by the
    // the convolution for the corresponding elements, and after that
    // the corresponding number of inverse FFTs.  The number of z-axis
    // forward and inverse FFTs to do in each sandwich is given by the
    // class member variable embed_block_size.
    //    NB: In this branch, the fftforwardtime and fftinversetime timer
    // variables measure the time for the x- and y-axis transforms only.
    // The convtime timer variable includes not only the "convolution"
    // time, but also the wrapping z-axis FFT times.

    // Calculate x- and y-axis FFTs of Mtemp.
    {
      OC_INDEX rxydim = ODTV_VECSIZE*rdimx*rdimy;
      OC_INDEX cxydim = ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx*cdimy;
      for(OC_INDEX m=0;m<rdimz;++m) {
        // x-direction transforms in plane "m"
#if REPORT_TIME
        fftxforwardtime.Start();
#endif // REPORT_TIME
        fftx.ForwardRealToComplexFFT(Mtemp+m*rxydim,Hxfrm+m*cxydim);
#if REPORT_TIME
        fftxforwardtime.Stop();
        fftyforwardtime.Start();
#endif // REPORT_TIME
        // y-direction transforms in plane "m"
        ffty.ForwardFFT(Hxfrm+m*cxydim);
#if REPORT_TIME
        fftyforwardtime.Stop();
#endif // REPORT_TIME
      }
    }

    // Do z-axis FFTs with embedded "convolution" operations.

    // Calculate field components in frequency domain.  Make use of
    // realness and even/odd properties of interaction matrices Axx.
    // Note that in transform space only the x>=0 half-space is
    // stored.
    // Symmetries: A00, A11, A22 are even in each coordinate
    //             A01 is odd in x and y, even in z.
    //             A02 is odd in x and z, even in y.
    //             A12 is odd in y and z, even in x.
    assert(adimx>=cdimx);
    assert(cdimy-adimy<adimy);
    assert(cdimz-adimz<adimz);
#if REPORT_TIME
    convtime.Start();
#endif // REPORT_TIME
    const OC_INDEX  jstride = ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx;
    const OC_INDEX  kstride = jstride*cdimy;
    const OC_INDEX ajstride = adimx;
    const OC_INDEX akstride = ajstride*adimy;

    for(j=0;j<adimy;++j) {
      // j>=0
      OC_INDEX  jindex = j*jstride;
      OC_INDEX ajindex = j*ajstride;
      fftz.AdjustArrayCount(ODTV_VECSIZE*embed_block_size);
      for(OC_INDEX m=0;m<cdimx;m+=embed_block_size) {
        // Do one block of forward z-direction transforms
        OC_INDEX istop = m + embed_block_size;
        if(embed_block_size>cdimx-m) {
          // Partial block
          fftz.AdjustArrayCount(ODTV_VECSIZE*(cdimx-m));
          istop = cdimx;
        }
        fftz.ForwardFFT(Hxfrm+jindex+m*ODTV_COMPLEXSIZE*ODTV_VECSIZE);
        // Do matrix-vector multiply ("convolution") for block
        for(k=0;k<adimz;++k) {
          // j>=0, k>=0
          OC_INDEX  kindex =  jindex + k*kstride;
          OC_INDEX akindex = ajindex + k*akstride;
          for(i=m;i<istop;++i) {
            OC_INDEX  index = ODTV_COMPLEXSIZE*ODTV_VECSIZE*i+kindex;
            OXS_FFT_REAL_TYPE Hx_re = Hxfrm[index];
            OXS_FFT_REAL_TYPE Hx_im = Hxfrm[index+1];
            OXS_FFT_REAL_TYPE Hy_re = Hxfrm[index+2];
            OXS_FFT_REAL_TYPE Hy_im = Hxfrm[index+3];
            OXS_FFT_REAL_TYPE Hz_re = Hxfrm[index+4];
            OXS_FFT_REAL_TYPE Hz_im = Hxfrm[index+5];

            const A_coefs& Aref = A[akindex+i]; 

            Hxfrm[index]   = Aref.A00*Hx_re + Aref.A01*Hy_re + Aref.A02*Hz_re;
            Hxfrm[index+1] = Aref.A00*Hx_im + Aref.A01*Hy_im + Aref.A02*Hz_im;
            Hxfrm[index+2] = Aref.A01*Hx_re + Aref.A11*Hy_re + Aref.A12*Hz_re;
            Hxfrm[index+3] = Aref.A01*Hx_im + Aref.A11*Hy_im + Aref.A12*Hz_im;
            Hxfrm[index+4] = Aref.A02*Hx_re + Aref.A12*Hy_re + Aref.A22*Hz_re;
            Hxfrm[index+5] = Aref.A02*Hx_im + Aref.A12*Hy_im + Aref.A22*Hz_im;
          }
        }
        for(k=adimz;k<cdimz;++k) {
          // j>=0, k<0
          OC_INDEX  kindex =  jindex + k*kstride;
          OC_INDEX akindex = ajindex + (cdimz-k)*akstride;
          for(i=m;i<istop;++i) {
            OC_INDEX  index = ODTV_COMPLEXSIZE*ODTV_VECSIZE*i+kindex;
            OXS_FFT_REAL_TYPE Hx_re = Hxfrm[index];
            OXS_FFT_REAL_TYPE Hx_im = Hxfrm[index+1];
            OXS_FFT_REAL_TYPE Hy_re = Hxfrm[index+2];
            OXS_FFT_REAL_TYPE Hy_im = Hxfrm[index+3];
            OXS_FFT_REAL_TYPE Hz_re = Hxfrm[index+4];
            OXS_FFT_REAL_TYPE Hz_im = Hxfrm[index+5];

            const A_coefs& Aref = A[akindex+i]; 

            // Flip signs on a02 and a12 as compared to the k>=0, j>=0 case
            // because a02 and a12 are odd in z.
            Hxfrm[index]   =  Aref.A00*Hx_re + Aref.A01*Hy_re - Aref.A02*Hz_re;
            Hxfrm[index+1] =  Aref.A00*Hx_im + Aref.A01*Hy_im - Aref.A02*Hz_im;
            Hxfrm[index+2] =  Aref.A01*Hx_re + Aref.A11*Hy_re - Aref.A12*Hz_re;
            Hxfrm[index+3] =  Aref.A01*Hx_im + Aref.A11*Hy_im - Aref.A12*Hz_im;
            Hxfrm[index+4] = -Aref.A02*Hx_re - Aref.A12*Hy_re + Aref.A22*Hz_re;
            Hxfrm[index+5] = -Aref.A02*Hx_im - Aref.A12*Hy_im + Aref.A22*Hz_im;
          }
        }
        // Do inverse z-direction transforms for block
        fftz.InverseFFT(Hxfrm+jindex+m*ODTV_COMPLEXSIZE*ODTV_VECSIZE);
      }
    }
    for(j=adimy;j<cdimy;++j) {
      // j<0
      OC_INDEX  jindex = j*jstride;
      OC_INDEX ajindex = (cdimy-j)*ajstride;
      fftz.AdjustArrayCount(ODTV_VECSIZE*embed_block_size);
      for(OC_INDEX m=0;m<cdimx;m+=embed_block_size) {
        // Do one block of forward z-direction transforms
        OC_INDEX istop = m + embed_block_size;
        if(embed_block_size>cdimx-m) {
          // Partial block
          fftz.AdjustArrayCount(ODTV_VECSIZE*(cdimx-m));
          istop = cdimx;
        }
        fftz.ForwardFFT(Hxfrm+jindex+m*ODTV_COMPLEXSIZE*ODTV_VECSIZE);
        // Do matrix-vector multiply ("convolution") for block
        for(k=0;k<adimz;++k) {
          // j<0, k>=0
          OC_INDEX  kindex =  jindex + k*kstride;
          OC_INDEX akindex = ajindex + k*akstride;
          for(i=m;i<istop;++i) {
            OC_INDEX  index = ODTV_COMPLEXSIZE*ODTV_VECSIZE*i+kindex;
            OXS_FFT_REAL_TYPE Hx_re = Hxfrm[index];
            OXS_FFT_REAL_TYPE Hx_im = Hxfrm[index+1];
            OXS_FFT_REAL_TYPE Hy_re = Hxfrm[index+2];
            OXS_FFT_REAL_TYPE Hy_im = Hxfrm[index+3];
            OXS_FFT_REAL_TYPE Hz_re = Hxfrm[index+4];
            OXS_FFT_REAL_TYPE Hz_im = Hxfrm[index+5];

            const A_coefs& Aref = A[akindex+i];

            // Flip signs on a01 and a12 as compared to the j>=0
            // case because a01 and a12 are odd in y.
            Hxfrm[index]   =  Aref.A00*Hx_re - Aref.A01*Hy_re + Aref.A02*Hz_re;
            Hxfrm[index+1] =  Aref.A00*Hx_im - Aref.A01*Hy_im + Aref.A02*Hz_im;
            Hxfrm[index+2] = -Aref.A01*Hx_re + Aref.A11*Hy_re - Aref.A12*Hz_re;
            Hxfrm[index+3] = -Aref.A01*Hx_im + Aref.A11*Hy_im - Aref.A12*Hz_im;
            Hxfrm[index+4] =  Aref.A02*Hx_re - Aref.A12*Hy_re + Aref.A22*Hz_re;
            Hxfrm[index+5] =  Aref.A02*Hx_im - Aref.A12*Hy_im + Aref.A22*Hz_im;
          }
        }
        for(k=adimz;k<cdimz;++k) {
          // j<0, k<0
          OC_INDEX  kindex =  jindex + k*kstride;
          OC_INDEX akindex = ajindex + (cdimz-k)*akstride;
          for(i=m;i<istop;++i) {
            OC_INDEX  index = ODTV_COMPLEXSIZE*ODTV_VECSIZE*i+kindex;
            OXS_FFT_REAL_TYPE Hx_re = Hxfrm[index];
            OXS_FFT_REAL_TYPE Hx_im = Hxfrm[index+1];
            OXS_FFT_REAL_TYPE Hy_re = Hxfrm[index+2];
            OXS_FFT_REAL_TYPE Hy_im = Hxfrm[index+3];
            OXS_FFT_REAL_TYPE Hz_re = Hxfrm[index+4];
            OXS_FFT_REAL_TYPE Hz_im = Hxfrm[index+5];

            const A_coefs& Aref = A[akindex+i]; 

            // Flip signs on a01 and a02 as compared to the k>=0, j>=0 case
            // because a01 is odd in y and even in z,
            //     and a02 is odd in z and even in y.
            // No change to a12 because it is odd in both y and z.
            Hxfrm[index]   =  Aref.A00*Hx_re - Aref.A01*Hy_re - Aref.A02*Hz_re;
            Hxfrm[index+1] =  Aref.A00*Hx_im - Aref.A01*Hy_im - Aref.A02*Hz_im;
            Hxfrm[index+2] = -Aref.A01*Hx_re + Aref.A11*Hy_re + Aref.A12*Hz_re;
            Hxfrm[index+3] = -Aref.A01*Hx_im + Aref.A11*Hy_im + Aref.A12*Hz_im;
            Hxfrm[index+4] = -Aref.A02*Hx_re + Aref.A12*Hy_re + Aref.A22*Hz_re;
            Hxfrm[index+5] = -Aref.A02*Hx_im + Aref.A12*Hy_im + Aref.A22*Hz_im;
          }
        }
        // Do inverse z-direction transforms for block
        fftz.InverseFFT(Hxfrm+jindex+m*ODTV_COMPLEXSIZE*ODTV_VECSIZE);
      }
    }
#if REPORT_TIME
    convtime.Stop();
#endif // REPORT_TIME

    // Do inverse y- and x-axis FFTs, to complete transform back into
    // space domain.
    {
      OC_INDEX m;
      OC_INDEX rxydim = ODTV_VECSIZE*rdimx*rdimy;
      OC_INDEX cxydim = ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx*cdimy;
      assert(3*sizeof(OXS_FFT_REAL_TYPE)<=sizeof(ThreeVector));
      OXS_FFT_REAL_TYPE* fptr
        = static_cast<OXS_FFT_REAL_TYPE*>(static_cast<void*>(&field[OC_INDEX(0)]));
      for(m=0;m<rdimz;++m) {
        // y-direction transforms
#if REPORT_TIME
        fftyinversetime.Start();
#endif // REPORT_TIME
        ffty.InverseFFT(Hxfrm+m*cxydim);
#if REPORT_TIME
        fftyinversetime.Stop();
        fftxinversetime.Start();
#endif // REPORT_TIME
        // x-direction transforms
        fftx.InverseComplexToRealFFT(Hxfrm+m*cxydim,fptr+m*rxydim);
#if REPORT_TIME
        fftxinversetime.Stop();
#endif // REPORT_TIME
      }

      if(3*sizeof(OXS_FFT_REAL_TYPE)<sizeof(ThreeVector)) {
        // The fftx.InverseComplexToRealFFT calls above assume the
        // target is an array of OXS_FFT_REAL_TYPE.  If ThreeVector is
        // not tightly packed, then this assumption is false; however we
        // can correct the problem by expanding the results in-place.
        // The only setting I know of where ThreeVector doesn't tight
        // pack is under the Borland bcc32 compiler on Windows x86 with
        // OXS_FFT_REAL_TYPE equal to "long double".  In that case
        // sizeof(long double) == 10, but sizeof(ThreeVector) == 36.
#if REPORT_TIME
        fftxinversetime.Start();
#endif // REPORT_TIME
        for(m = rsize - 1; m>=0 ; --m) {
          ThreeVector temp(fptr[ODTV_VECSIZE*m],fptr[ODTV_VECSIZE*m+1],fptr[ODTV_VECSIZE*m+2]);
          field[m] = temp;
        }
#if REPORT_TIME
        fftxinversetime.Stop();
#endif // REPORT_TIME
      }

    }

  } // if(!embed_convolution)

#if REPORT_TIME
  dottime.Start();
#endif // REPORT_TIME
  // Calculate pointwise energy density: -0.5*MU0*<M,H>
  const OXS_FFT_REAL_TYPE emult =  -0.5 * MU0;
  for(i=0;i<rsize;++i) {
    OXS_FFT_REAL_TYPE dot = spin[i]*field[i];
    energy[i] = emult * dot * Ms[i];
  }
#if REPORT_TIME
  dottime.Stop();
#endif // REPORT_TIME
}

#endif // OOMMF_THREADS
