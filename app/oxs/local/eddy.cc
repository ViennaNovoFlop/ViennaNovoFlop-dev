/* FILE: demag-threaded.cc            -*-Mode: c++-*-
 *
 * Threaded version of demag.cc (q.v.).
 */

#include "eddy.h"  // Includes definition of OOMMF_THREADS macro

////////////////// MULTI-THREADED IMPLEMENTATION  ///////////////
#if OOMMF_THREADS

#include <assert.h>
#include <string>
#include <vector>

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
OXS_EXT_REGISTER(Oxs_Eddy);

#ifndef VERBOSE_DEBUG
# define VERBOSE_DEBUG 0
#endif

// Size of threevector.  This macro is defined for code legibility
// and maintenance; it should always be "3".
#define ODTV_VECSIZE 3

// Size of complex value, in real units.  This macro is defined for code
// legibility and maintenance; it should always be "2".
#define ODTV_COMPLEXSIZE 2

class _Oxs_EddyJobControl {
public:

  _Oxs_EddyJobControl() : imax(0),big_block_limit(0),big_blocksize(0),
                            small_blocksize(0),next_job_start(0) {}

  // Note: The initializer is not thread-safe.  It should be called
  // only from the control thread.
  void Init(OC_INDEX i_imax,OC_INT4m threadcount,OC_INDEX minjobsize);

  // ClaimJob is thread safe.  It is called from inside worker
  // threads to get job assignments.  The return value istop
  // is guaranteed to be less than imax.
  void ClaimJob(OC_INDEX& istart,OC_INDEX& istop);

private:
  // Variable next_job_start marks the next available job.  Jobs are
  // incremented by big_blocksize until big_block_limit is reached,
  // after which jobs are handed out of size small_blocksize.
  OC_INDEX imax;
  OC_INDEX big_block_limit;
  OC_INDEX big_blocksize;
  OC_INDEX small_blocksize;
  OC_INDEX next_job_start;
  Oxs_Mutex job_mutex;
};

void _Oxs_EddyJobControl::Init
(OC_INDEX i_imax,
 OC_INT4m threadcount,
 OC_INDEX minjobsize)
{
  imax = i_imax;
  next_job_start = 0;
  assert(imax>=1 && threadcount>=1 && minjobsize>=1);
  if(threadcount<1) threadcount = 1; // Safety
  if(minjobsize<1) minjobsize = 1;

  // big_proportion + small_proportion should be <= 1.0
#if 1
  const OC_REAL8m small_proportion = 0.05;
  const OC_REAL8m big_proportion   = 0.95;
#else
  // On some systems (8 threads, non-numa), slightly more
  // and slightly larger small blocks works slightly better.
  const OC_REAL8m small_proportion = 0.10;
  const OC_REAL8m big_proportion   = 0.80;
#endif

  small_blocksize = OC_INDEX(floor(0.5+(small_proportion*imax)/threadcount));
  // Round small_blocksize to an integral multiple of minjobsize
  if(small_blocksize<=minjobsize) {
    small_blocksize = minjobsize;
  } else {
    small_blocksize = minjobsize*((small_blocksize+minjobsize/2)/minjobsize);
  }

  // Select big_block_limit close to big_proportion of imax, but
  // round so that remainder is an integral multiple of small_blocksize.
  OC_INDEX remainder = imax - OC_INDEX(floor(0.5+big_proportion*imax));
  if(0<remainder && remainder<=small_blocksize) {
    remainder = small_blocksize;
  } else {
    remainder
      = small_blocksize*(OC_INDEX(floor((0.5+remainder)/small_blocksize)));
  }
  big_block_limit = imax - remainder;

  // Pick big_blocksize close to (but not smaller than)
  // big_block_limit/threadcount, and make it a multiple of minjobsize.
  // The size of the last (short) big block will be big_block_limit -
  // big_blocksize*(threadcount-1) (see ::ClaimJob); if imax is a
  // multiple of minjobsize, then big_block_limit will also be a
  // multiple of minjobsize, and so in this case the last (short) big
  // block will also be a multiple of minjobsize.
  big_blocksize = minjobsize*
    ((big_block_limit+minjobsize*threadcount-1)/(minjobsize*threadcount));

#if 0
  // On some systems (8 threads, non-numa), breaking up the
  // big blocks seems to help.
  if(big_blocksize>=4*minjobsize) {
    big_blocksize = (big_blocksize+3)/4;
  }
#endif

static int foocount=0; /**/ // asdf
if(foocount<5) {
  ++foocount;
  fprintf(stderr,"Total limit=%ld, big limit=%ld,"
          " big size=%ld, small size=%ld\n",
          long(imax),long(big_block_limit),
          long(big_blocksize),long(small_blocksize));
 } /**/ // asdf

}

void _Oxs_EddyJobControl::ClaimJob(OC_INDEX& istart,OC_INDEX& istop)
{
  OC_INDEX tmp_start,tmp_stop;
  job_mutex.Lock();
  if((tmp_start = next_job_start)<big_block_limit) {
    if(next_job_start + big_blocksize>big_block_limit) {
      tmp_stop = next_job_start = big_block_limit;
    } else {
      tmp_stop = (next_job_start += big_blocksize);
    }
  } else {
    tmp_stop = (next_job_start += small_blocksize);
  }
  job_mutex.Unlock();
  if(tmp_stop>imax) tmp_stop=imax;  // Guarantee that istop is in-range.
  istart = tmp_start;
  istop  = tmp_stop;
}

////////////////////////////////////////////////////////////////////////
// Oxs_Eddy::Oxs_FFTLocker provides thread-specific instance of FFT
// objects.

Oxs_Eddy::Oxs_FFTLocker::Oxs_FFTLocker
(const Oxs_Eddy::Oxs_FFTLocker_Info& info)
  : ifftx_scratch(0), fftz_Hwork(0), fftyconvolve_Hwork(0),
    A_copy(0),
    ifftx_scratch_size(0), fftz_Hwork_size(0), fftyconvolve_Hwork_size(0),
    A_copy_size(0)
{
  // Check import data
  assert(info.rdimx>0 && info.rdimy>0 && info.rdimz>0 &&
         info.cdimx>0 && info.cdimy>0 && info.cdimz>0 &&
         info.embed_block_size>0);
         
  // FFT objects
  fftx.SetDimensions(info.rdimx,
                     (info.cdimx==1 ? 1 : 2*(info.cdimx-1)), info.rdimy);
  ffty.SetDimensions(info.rdimy, info.cdimy,
                     ODTV_COMPLEXSIZE*ODTV_VECSIZE*info.cdimx,
                     ODTV_VECSIZE*info.cdimx);
  fftz.SetDimensions(info.rdimz, info.cdimz,
                     ODTV_COMPLEXSIZE*ODTV_VECSIZE*info.cdimx*info.cdimy,
                     ODTV_VECSIZE*info.cdimx*info.cdimy);

  // Workspace
  ifftx_scratch_size = info.rdimx*ODTV_VECSIZE;
  ifftx_scratch = static_cast<OXS_FFT_REAL_TYPE*>
    (Oc_AllocThreadLocal(ifftx_scratch_size*sizeof(OXS_FFT_REAL_TYPE)));
  if(info.cdimz>1) {
    fftz_Hwork_size
      = 2*ODTV_VECSIZE*ODTV_COMPLEXSIZE*info.embed_block_size*info.cdimz;
    fftz_Hwork = static_cast<OXS_FFT_REAL_TYPE*>
      (Oc_AllocThreadLocal(fftz_Hwork_size*sizeof(OXS_FFT_REAL_TYPE)));
  } else {
    fftyconvolve_Hwork_size
      = ODTV_VECSIZE*ODTV_COMPLEXSIZE*info.embed_block_size*info.cdimy;
    fftyconvolve_Hwork = static_cast<OXS_FFT_REAL_TYPE*>
      (Oc_AllocThreadLocal(fftyconvolve_Hwork_size*sizeof(OXS_FFT_REAL_TYPE)));
  }
}

Oxs_Eddy::Oxs_FFTLocker::~Oxs_FFTLocker()
{
  if(ifftx_scratch) Oc_FreeThreadLocal(ifftx_scratch,
                       ifftx_scratch_size*sizeof(OXS_FFT_REAL_TYPE));
  if(fftz_Hwork)    Oc_FreeThreadLocal(fftz_Hwork,
                       fftz_Hwork_size*sizeof(OXS_FFT_REAL_TYPE));
  if(fftyconvolve_Hwork)
                    Oc_FreeThreadLocal(fftyconvolve_Hwork,
                       fftyconvolve_Hwork_size*sizeof(OXS_FFT_REAL_TYPE));
  if(A_copy)        Oc_FreeThreadLocal(A_copy,
                       A_copy_size*sizeof(Oxs_Eddy::A_coefs));
}

////////////////////////////////////////////////////////////////////////
// Oxs_Eddy Constructor
Oxs_Eddy::Oxs_Eddy(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Oxs_Energy(name,newdtr,argstr),
    rdimx(0),rdimy(0),rdimz(0),cdimx(0),cdimy(0),cdimz(0),
    adimx(0),adimy(0),adimz(0),
    mesh_id(0),
    A(0),Hxfrm(0),asymptotic_radius(-1),Jtemp(0),
    MaxThreadCount(Oc_GetMaxThreadCount()),
    embed_block_size(0)
{
  asymptotic_radius = GetRealInitValue("asymptotic_radius",32.0);
  /// Default value of -1 disables use of asymptotic approximation.

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
// Base magnitude of the current flow.  Positive values mean
// current flow is in the direction specified by J_direction.
// Negative values mean actual current flow is in direction
// opposite J_direction.
OXS_GET_INIT_EXT_OBJECT("J",Oxs_ScalarField,J_init);    // Required

// Stage-varying multiplier for J.  Note sign convention, as
// outlined in previous stanza.
if(HasInitValue("J_profile")) {
  has_J_profile = 1;
  String cmdoptreq = GetStringInitValue("J_profile_args",
                                        "stage stage_time total_time");
  J_profile_opts.push_back(Oxs_TclCommandLineOption("stage",1));
  J_profile_opts.push_back(Oxs_TclCommandLineOption("stage_time",1));
  J_profile_opts.push_back(Oxs_TclCommandLineOption("total_time",1));
  J_profile_cmd.SetBaseCommand(InstanceName(),
                 director->GetMifInterp(),
                 GetStringInitValue("J_profile"),
                 Oxs_ParseTclCommandLineRequest(InstanceName(),
                                                J_profile_opts,
                                                cmdoptreq));
}

// "mp" is required input
OXS_GET_INIT_EXT_OBJECT("mp",Oxs_VectorField,mp_init);

  VerifyAllInitArgsUsed();
}

Oxs_Eddy::~Oxs_Eddy() {
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

OC_BOOL Oxs_Eddy::Init()
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
  mp.Release();
  J.Release();
  J_scratch.Release();
  
  return Oxs_Energy::Init();
}

void Oxs_Eddy::ReleaseMemory() const
{ // Conceptually const
  if(A!=0) { delete[] A; A=0; }
  if(Hxfrm!=0)       { delete[] Hxfrm;       Hxfrm=0;       }
  if(Jtemp!=0)       { delete[] Jtemp;       Jtemp=0;       }
  rdimx=rdimy=rdimz=0;
  cdimx=cdimy=cdimz=0;
  adimx=adimy=adimz=0;
}

void Oxs_Eddy::FillCoefficientArrays(const Oxs_Mesh* genmesh) const
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

  Jtemp = new OXS_FFT_REAL_TYPE[ODTV_VECSIZE*rdimx*rdimy*rdimz];
  /// Temporary space to hold Ms[]*m[].  The plan is to make this space
  /// unnecessary by introducing FFT functions that can take Ms as input
  /// and do the multiplication on the fly.

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

  // Compute block size for "convolution" embedded with inner FFT's.
  OC_INDEX footprint
    = ODTV_COMPLEXSIZE*ODTV_VECSIZE*sizeof(OXS_FFT_REAL_TYPE) // Data
    + sizeof(A_coefs)                           // Interaction matrix
    + 2*ODTV_COMPLEXSIZE*sizeof(OXS_FFT_REAL_TYPE); // Roots of unity
  if(cdimz<2) {
    footprint *= cdimy;  // Embed convolution with y-axis FFT's
  } else {
    footprint *= cdimz;  // Embed convolution with z-axis FFT's
  }
  OC_INDEX trialsize = cache_size/(3*footprint); // "3" is fudge factor
  // Tweak trialsize to a good divisor of cdimx
  if(trialsize<=1) {
    trialsize = 1;
  } else if(trialsize<cdimx) {
    OC_INDEX bc = (cdimx + trialsize/2)/trialsize;
    trialsize = (cdimx + bc - 1)/bc;
    if(trialsize>32) {
      // Round to next multiple of eight
      if(trialsize%8) { trialsize += 8 - (trialsize%8); }
    }
  }
  if(trialsize>cdimx) trialsize = cdimx;
  embed_block_size = trialsize;
  fprintf(stderr,"Embed block size=%ld\n",long(embed_block_size)); /**/


  // The following 3 statements are cribbed from
  // Oxs_FFT3DThreeVector::SetDimensions().  The corresponding
  // code using that class is
  //
  //  Oxs_FFT3DThreeVector fft;
  //  fft.SetDimensions(rdimx,rdimy,rdimz,cdimx,cdimy,cdimz);
  //  fft.GetLogicalDimensions(ldimx,ldimy,ldimz);
  //
  // Set up non-threaded instances of FFT objects.  Threaded code sets
  // up it own copies using keyword "Oxs_EddyFFTObject" to
  // Oxs_ThreadRunObj::local_locker
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

  // Dimensions for the demag tensor A.
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
               " to fit in an OC_INDEX variable");
    throw Oxs_ExtError(this,msg);
  }

  // Allocate memory for FFT xfrm target H, and scratch space
  // for computing interaction coefficients
  Hxfrm = new OXS_FFT_REAL_TYPE[xfrm_size];
  OXS_FFT_REAL_TYPE* scratch = new OXS_FFT_REAL_TYPE[scratch_size];

  if(Hxfrm==NULL || scratch==NULL) {
    // Safety check for those machines on which new[] doesn't throw
    // BadAlloc.
    String msg = String("Insufficient memory in Eddy setup.");
    throw Oxs_ExtError(this,msg);
  }

  // According (16) in Newell's paper, the demag field is given by
  //                        H = -N*M
  // where N is the "demagnetizing tensor," with components Nxx, Nxy,
  // etc.  With the '-1' in 'scale' we store '-N' instead of 'N',
  // so we don't have to multiply the output from the FFT + iFFT
  // by -1 in GetEnergy() below.

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

#if REPORT_TIME
    inittime.Stop();
#endif // REPORT_TIME
}


class _Oxs_EddyFFTxThread : public Oxs_ThreadRunObj {
public:
  const Oxs_MeshValue<ThreeVector>* mp;
  const Oxs_MeshValue<OC_REAL8m>* J;

  OXS_FFT_REAL_TYPE* rarr;
  OXS_FFT_REAL_TYPE* carr;

  Oxs_Eddy::Oxs_FFTLocker_Info locker_info;
  Oxs_FFT1DThreeVector* fftx;

  static _Oxs_EddyJobControl job_control;

  OC_INDEX mp_xdim;
  OC_INDEX mp_xydim;

  OC_INDEX j_dim;
  OC_INDEX j_rstride;
  OC_INDEX j_cstride;

  OC_INDEX k_rstride;
  OC_INDEX k_cstride;

  OC_INDEX jk_max;

  enum { INVALID, FORWARD, INVERSE } direction;
  _Oxs_EddyFFTxThread()
    : rarr(0),carr(0),
      fftx(0),
      mp_xdim(0),mp_xydim(0),
      j_dim(0),j_rstride(0),j_cstride(0),
      k_rstride(0),k_cstride(0),
      jk_max(0),
      direction(INVALID) {}
  void Cmd(int threadnumber, void* data);
};

_Oxs_EddyJobControl _Oxs_EddyFFTxThread::job_control;

void _Oxs_EddyFFTxThread::Cmd(int /* threadnumber */, void* /* data */)
{
  // Thread local storage
  if(!fftx) {
    Oxs_ThreadMapDataObject* foo = local_locker.GetItem(locker_info.name);
    if(!foo) {
      // Oxs_FFTLocker object not constructed
      foo = new Oxs_Eddy::Oxs_FFTLocker(locker_info);
      local_locker.AddItem(locker_info.name,foo);
    }
    fftx = &(dynamic_cast<Oxs_Eddy::Oxs_FFTLocker*>(foo)->fftx);
    if(!fftx) {
      Oxs_ThreadError::SetError(String("Error in"
         "_Oxs_EddyFFTxThread::Cmd(): fftx downcast failed."));
      return;
    }
  }

  while(1) {
    OC_INDEX jkstart,jkstop;
    job_control.ClaimJob(jkstart,jkstop);
    if(jkstart>=jk_max) break;

    OC_INDEX kstart = jkstart/j_dim;
    OC_INDEX jstart = jkstart - kstart*j_dim;

    OC_INDEX kstop = jkstop/j_dim;
    OC_INDEX jstop = jkstop - kstop*j_dim;

    for(OC_INDEX k=kstart;k<=kstop;++k) {
      OC_INDEX j_line_stop = j_dim;
      if(k==kstop) j_line_stop = jstop;
      if(jstart<j_line_stop) {
        fftx->AdjustArrayCount(j_line_stop - jstart);
        if(direction == FORWARD) {
          const OC_INDEX istart = jstart*mp_xdim + k*mp_xydim;
          fftx->ForwardRealToComplexFFT(static_cast<const OC_REAL8m*>(&((*mp)[istart].x)), // CHEAT
                                        carr+jstart*j_cstride+k*k_cstride,
                                        static_cast<const OC_REAL8m*>(&((*J)[istart]))); // CHEAT
        } else { // direction == INVERSE
          fftx->InverseComplexToRealFFT(carr+jstart*j_cstride+k*k_cstride,
                                        rarr+jstart*j_rstride+k*k_rstride);
        }
      }
      jstart=0;
    }
  }
}

class _Oxs_EddyiFFTxDotThread : public Oxs_ThreadRunObj {
public:
  OXS_FFT_REAL_TYPE* carr;
  const Oxs_MeshValue<ThreeVector> *mp_ptr;
  const Oxs_MeshValue<OC_REAL8m> *J_ptr;
  Oxs_ComputeEnergyData* oced_ptr;

  Oxs_Eddy::Oxs_FFTLocker_Info locker_info;
  Oxs_Eddy::Oxs_FFTLocker* locker;

  static _Oxs_EddyJobControl job_control;

  Nb_Xpfloat energy_sum;

  OC_INDEX rdimx;

  OC_INDEX j_dim;
  OC_INDEX j_rstride;
  OC_INDEX j_cstride;

  OC_INDEX k_rstride;
  OC_INDEX k_cstride;

  OC_INDEX jk_max;

  _Oxs_EddyiFFTxDotThread()
    : carr(0),
      mp_ptr(0), J_ptr(0), oced_ptr(0), locker(0),
      energy_sum(0.0),
      rdimx(0), 
      j_dim(0), j_rstride(0), j_cstride(0),
      k_rstride(0), k_cstride(0),
      jk_max(0) {}
  void Cmd(int threadnumber, void* data);
};

_Oxs_EddyJobControl _Oxs_EddyiFFTxDotThread::job_control;

void _Oxs_EddyiFFTxDotThread::Cmd(int /* threadnumber */, void* /* data */)
{
  // Thread local storage
  if(!locker) {
    Oxs_ThreadMapDataObject* foo = local_locker.GetItem(locker_info.name);
    if(!foo) {
      // Oxs_FFTLocker object not constructed
      foo = new Oxs_Eddy::Oxs_FFTLocker(locker_info);
      local_locker.AddItem(locker_info.name,foo);
    }
    locker = dynamic_cast<Oxs_Eddy::Oxs_FFTLocker*>(foo);
    if(!locker) {
      Oxs_ThreadError::SetError(String("Error in"
         "_Oxs_EddyFFTxDotThread::Cmd(): locker downcast failed."));
      return;
    }
  }
  Oxs_FFT1DThreeVector* const fftx = &(locker->fftx);
  OXS_FFT_REAL_TYPE* scratch = locker->ifftx_scratch;

  const OXS_FFT_REAL_TYPE emult =  -0.5 * MU0;
  energy_sum = 0.0;

  const Oxs_MeshValue<ThreeVector>& mp = *mp_ptr;
  const Oxs_MeshValue<OC_REAL8m>& J = *J_ptr;
  Oxs_ComputeEnergyData& oced = *oced_ptr;

  const OC_INDEX ijstride = j_rstride/ODTV_VECSIZE;
  const OC_INDEX ikstride = k_rstride/ODTV_VECSIZE;

  fftx->AdjustArrayCount(1);

  while(1) {
    OC_INDEX jkstart,jkstop;
    job_control.ClaimJob(jkstart,jkstop);
    if(jkstart>=jk_max) break;

    OC_INDEX kstart = jkstart/j_dim;
    OC_INDEX jstart = jkstart - kstart*j_dim;

    OC_INDEX kstop = jkstop/j_dim;
    OC_INDEX jstop = jkstop - kstop*j_dim;

    for(OC_INDEX k=kstart;k<=kstop;++k,jstart=0) {
      OC_INDEX j_line_stop = j_dim;
      if(k==kstop) j_line_stop = jstop;
      if(jstart>=j_line_stop) continue;

      for(OC_INDEX j=jstart;j<j_line_stop;++j) {
        fftx->InverseComplexToRealFFT(carr+j*j_cstride+k*k_cstride,scratch);

        const OC_INDEX ioffset = (j*ijstride+k*ikstride);

        if(oced.H) {
          for(OC_INDEX i=0;i<rdimx;++i) {
            (*oced.H)[ioffset + i].Set(scratch[3*i],scratch[3*i+1],scratch[3*i+2]);
          }
        }

        if(oced.H_accum) {
          for(OC_INDEX i=0;i<rdimx;++i) {
            (*oced.H_accum)[ioffset + i]
              += ThreeVector(scratch[3*i],scratch[3*i+1],scratch[3*i+2]);
          }
        }

        OC_INDEX i;
        for(i=0;i<rdimx;++i) {
          const ThreeVector& m = mp[ioffset+i];
          const OXS_FFT_REAL_TYPE& tHx = scratch[3*i];
          const OXS_FFT_REAL_TYPE& tHy = scratch[3*i+1];
          const OXS_FFT_REAL_TYPE& tHz = scratch[3*i+2];

          OC_REAL8m ei = emult * J[ioffset+i] * (m.x * tHx + m.y * tHy + m.z * tHz);

          OC_REAL8m tx = m.y * tHz - m.z * tHy; // mxH
          OC_REAL8m ty = m.z * tHx - m.x * tHz;
          OC_REAL8m tz = m.x * tHy - m.y * tHx;

          energy_sum += ei;

          if(oced.energy)       (*oced.energy)[ioffset + i] = ei;
          if(oced.energy_accum) (*oced.energy_accum)[ioffset + i] += ei;

          if(oced.mxH)          (*oced.mxH)[ioffset + i] = ThreeVector(tx,ty,tz);
          if(oced.mxH_accum)    (*oced.mxH_accum)[ioffset + i] += ThreeVector(tx,ty,tz);

#if defined(__GNUC__) && __GNUC__ == 4 \
  && defined(__GNUC_MINOR__) && __GNUC_MINOR__ <= 1
          Oc_Nop(i);  // gcc 4.1.2 20080704 (Red Hat 4.1.2-44) (others?)
          /// has an optimization bug that segfaults the above code block
          /// if -fprefetch-loop-arrays is enabled.  This nop call is a
          /// workaround.
          /// i686-apple-darwin9-gcc-4.0.1 (OS X Leopard) (others?) has an
          /// optimization bug that segfaults the above code block if
          /// -fstrict-aliasing (or -fast, which includes -fstrict-aliasing)
          /// is enabled.  This nop call is a workaround.
#endif

        } // for-i

      } // for-j
    }   // for-k
  }     // while(1)
}

class _Oxs_EddyFFTyThread : public Oxs_ThreadRunObj {
public:
  OXS_FFT_REAL_TYPE* carr;

  Oxs_Eddy::Oxs_FFTLocker_Info locker_info;
  Oxs_FFTStrided* ffty;

  static _Oxs_EddyJobControl job_control;

  OC_INDEX k_stride;
  OC_INDEX k_dim;
  OC_INDEX i_dim;

  enum { INVALID, FORWARD, INVERSE } direction;
  _Oxs_EddyFFTyThread()
    : carr(0), ffty(0),
      k_stride(0), k_dim(0),
      i_dim(0), direction(INVALID) {}
  void Cmd(int threadnumber, void* data);
};

_Oxs_EddyJobControl _Oxs_EddyFFTyThread::job_control;

void _Oxs_EddyFFTyThread::Cmd(int /* threadnumber */, void* /* data */)
{
  // Thread local storage
  if(!ffty) {
    Oxs_ThreadMapDataObject* foo = local_locker.GetItem(locker_info.name);
    if(!foo) {
      // Oxs_FFTLocker object not constructed
      foo = new Oxs_Eddy::Oxs_FFTLocker(locker_info);
      local_locker.AddItem(locker_info.name,foo);
    }
    ffty = &(dynamic_cast<Oxs_Eddy::Oxs_FFTLocker*>(foo)->ffty);
    if(!ffty) {
      Oxs_ThreadError::SetError(String("Error in"
         "_Oxs_EddyFFTyThread::Cmd(): ffty downcast failed."));
      return;
    }
  }

  const OC_INDEX istride = ODTV_COMPLEXSIZE;

  while(1) {
    OC_INDEX ikstart,ikstop;
    job_control.ClaimJob(ikstart,ikstop);

    OC_INDEX kstart = ikstart/i_dim;
    OC_INDEX istart = ikstart - kstart*i_dim;

    OC_INDEX kstop = ikstop/i_dim;
    OC_INDEX istop = ikstop - kstop*i_dim;

    if(kstart>=k_dim) break;
    for(OC_INDEX k=kstart;k<=kstop;++k,istart=0) {
      OC_INDEX i_line_stop = i_dim;
      if(k==kstop) i_line_stop = istop;
      if(istart>=i_line_stop) continue;
      ffty->AdjustArrayCount(i_line_stop - istart);
      if(direction == FORWARD) {
        ffty->ForwardFFT(carr+istart*istride+k*k_stride);
      } else { // direction == INVERSE
        ffty->InverseFFT(carr+istart*istride+k*k_stride);
      }
    }

  }
}

class _Oxs_EddyFFTyConvolveThread : public Oxs_ThreadRunObj {
  // *ONLY* for use when cdimz==1
public:
  OXS_FFT_REAL_TYPE* Hxfrm;
  const Oxs_Eddy::A_coefs* A;

  Oxs_Eddy::Oxs_FFTLocker_Info locker_info;
  Oxs_Eddy::Oxs_FFTLocker* locker;

  static _Oxs_EddyJobControl job_control;

  OC_INDEX embed_block_size;
  OC_INDEX jstride,ajstride;
  OC_INDEX i_dim;

  OC_INDEX rdimy,adimy,cdimy;

  _Oxs_EddyFFTyConvolveThread()
    : Hxfrm(0), A(0), locker(0),
      embed_block_size(0),
      jstride(0),ajstride(0),
      i_dim(0),
      rdimy(0), adimy(0), cdimy(0) {}
  void Cmd(int threadnumber, void* data);
};

_Oxs_EddyJobControl _Oxs_EddyFFTyConvolveThread::job_control;

void _Oxs_EddyFFTyConvolveThread::Cmd(int /* threadnumber */, void* /* data */)
{
  // Thread local storage
  if(!locker) {
    Oxs_ThreadMapDataObject* foo = local_locker.GetItem(locker_info.name);
    if(!foo) {
      // Oxs_FFTLocker object not constructed
      foo = new Oxs_Eddy::Oxs_FFTLocker(locker_info);
      local_locker.AddItem(locker_info.name,foo);
    }
    locker = dynamic_cast<Oxs_Eddy::Oxs_FFTLocker*>(foo);
    if(!locker) {
      Oxs_ThreadError::SetError(String("Error in"
         "_Oxs_EddyFFTyConvolveThread::Cmd(): locker downcast failed."));
      return;
    }
  }
  Oxs_FFTStrided* const ffty = &(locker->ffty);

  // Hwork:  Data is copied from Hxfrm into and out of this space
  // on each m increment.
  const OC_INDEX Hwstride = ODTV_VECSIZE*ODTV_COMPLEXSIZE*embed_block_size;
  OXS_FFT_REAL_TYPE* const Hwork = locker->fftyconvolve_Hwork;

  // Adjust ffty to use Hwork
  ffty->AdjustInputDimensions(rdimy,Hwstride,
                              ODTV_VECSIZE*embed_block_size);

  const OC_INDEX istride = ODTV_COMPLEXSIZE*ODTV_VECSIZE;
  while(1) {
    OC_INDEX istart,istop;
    job_control.ClaimJob(istart,istop);

    if(istart>=i_dim) break;

    for(OC_INDEX ix=istart;ix<istop;ix+=embed_block_size) {
      OC_INDEX j;

      OC_INDEX ix_end = ix + embed_block_size;
      if(ix_end>istop) ix_end = istop;

      // Copy data from Hxfrm into Hwork
      const size_t Hcopy_line_size
        = static_cast<size_t>(istride*(ix_end-ix))*sizeof(OXS_FFT_REAL_TYPE);
      for(j=0;j<rdimy;++j) {
        const OC_INDEX windex = j*Hwstride;
        const OC_INDEX hindex = j*jstride + ix*istride;
        memcpy(Hwork+windex,Hxfrm+hindex,Hcopy_line_size);
      }

      ffty->AdjustArrayCount((ix_end-ix)*ODTV_VECSIZE);
      ffty->ForwardFFT(Hwork);
        
      { // j==0
        for(OC_INDEX i=ix;i<ix_end;++i) {
          const Oxs_Eddy::A_coefs& Aref = A[i];
          {
            OC_INDEX  index = istride*(i-ix);
            OXS_FFT_REAL_TYPE Hx_re = Hwork[index];
            OXS_FFT_REAL_TYPE Hx_im = Hwork[index+1];
            OXS_FFT_REAL_TYPE Hy_re = Hwork[index+2];
            OXS_FFT_REAL_TYPE Hy_im = Hwork[index+3];
            OXS_FFT_REAL_TYPE Hz_re = Hwork[index+4];
            OXS_FFT_REAL_TYPE Hz_im = Hwork[index+5];

            Hwork[index]   = Aref.A00*Hx_re + Aref.A01*Hy_re + Aref.A02*Hz_re;
            Hwork[index+1] = Aref.A00*Hx_im + Aref.A01*Hy_im + Aref.A02*Hz_im;
            Hwork[index+2] = Aref.A01*Hx_re + Aref.A11*Hy_re + Aref.A12*Hz_re;
            Hwork[index+3] = Aref.A01*Hx_im + Aref.A11*Hy_im + Aref.A12*Hz_im;
            Hwork[index+4] = Aref.A02*Hx_re + Aref.A12*Hy_re + Aref.A22*Hz_re;
            Hwork[index+5] = Aref.A02*Hx_im + Aref.A12*Hy_im + Aref.A22*Hz_im;
          }
        }
      }

      for(j=1;j<cdimy/2;++j) {
        OC_INDEX ajindex = j*ajstride;
        OC_INDEX  jindex = j*Hwstride;
        OC_INDEX  j2index = (cdimy-j)*Hwstride;

        for(OC_INDEX i=ix;i<ix_end;++i) {
          const Oxs_Eddy::A_coefs& Aref = A[ajindex+i];
          { // j>0
            OC_INDEX  index = jindex + istride*(i-ix);
            OXS_FFT_REAL_TYPE Hx_re = Hwork[index];
            OXS_FFT_REAL_TYPE Hx_im = Hwork[index+1];
            OXS_FFT_REAL_TYPE Hy_re = Hwork[index+2];
            OXS_FFT_REAL_TYPE Hy_im = Hwork[index+3];
            OXS_FFT_REAL_TYPE Hz_re = Hwork[index+4];
            OXS_FFT_REAL_TYPE Hz_im = Hwork[index+5];

            Hwork[index]   = Aref.A00*Hx_re + Aref.A01*Hy_re + Aref.A02*Hz_re;
            Hwork[index+1] = Aref.A00*Hx_im + Aref.A01*Hy_im + Aref.A02*Hz_im;
            Hwork[index+2] = Aref.A01*Hx_re + Aref.A11*Hy_re + Aref.A12*Hz_re;
            Hwork[index+3] = Aref.A01*Hx_im + Aref.A11*Hy_im + Aref.A12*Hz_im;
            Hwork[index+4] = Aref.A02*Hx_re + Aref.A12*Hy_re + Aref.A22*Hz_re;
            Hwork[index+5] = Aref.A02*Hx_im + Aref.A12*Hy_im + Aref.A22*Hz_im;
          }
          { // j2<0
            OC_INDEX  index2 = j2index + istride*(i-ix);
            OXS_FFT_REAL_TYPE Hx_re = Hwork[index2];
            OXS_FFT_REAL_TYPE Hx_im = Hwork[index2+1];
            OXS_FFT_REAL_TYPE Hy_re = Hwork[index2+2];
            OXS_FFT_REAL_TYPE Hy_im = Hwork[index2+3];
            OXS_FFT_REAL_TYPE Hz_re = Hwork[index2+4];
            OXS_FFT_REAL_TYPE Hz_im = Hwork[index2+5];

            // Flip signs on a01 and a12 as compared to the j>=0
            // case because a01 and a12 are odd in y.
            Hwork[index2]   =  Aref.A00*Hx_re - Aref.A01*Hy_re + Aref.A02*Hz_re;
            Hwork[index2+1] =  Aref.A00*Hx_im - Aref.A01*Hy_im + Aref.A02*Hz_im;
            Hwork[index2+2] = -Aref.A01*Hx_re + Aref.A11*Hy_re - Aref.A12*Hz_re;
            Hwork[index2+3] = -Aref.A01*Hx_im + Aref.A11*Hy_im - Aref.A12*Hz_im;
            Hwork[index2+4] =  Aref.A02*Hx_re - Aref.A12*Hy_re + Aref.A22*Hz_re;
            Hwork[index2+5] =  Aref.A02*Hx_im - Aref.A12*Hy_im + Aref.A22*Hz_im;
          }
        }
      }

      // Note special case: If cdimy==adimy==1, then cdimy/2=0, and so j
      // coming out from the previous loop is 1, which is different from
      // cdimy/2.  In this case we want to run *only* the j=0 loop farther
      // above, and not either of the others.
      for(;j<adimy;++j) {
        OC_INDEX ajindex = j*ajstride;
        OC_INDEX  jindex = j*Hwstride;
        for(OC_INDEX i=ix;i<ix_end;++i) {
          const Oxs_Eddy::A_coefs& Aref = A[ajindex+i];
          { // j>0
            OC_INDEX  index = jindex + istride*(i-ix);
            OXS_FFT_REAL_TYPE Hx_re = Hwork[index];
            OXS_FFT_REAL_TYPE Hx_im = Hwork[index+1];
            OXS_FFT_REAL_TYPE Hy_re = Hwork[index+2];
            OXS_FFT_REAL_TYPE Hy_im = Hwork[index+3];
            OXS_FFT_REAL_TYPE Hz_re = Hwork[index+4];
            OXS_FFT_REAL_TYPE Hz_im = Hwork[index+5];

            Hwork[index]   = Aref.A00*Hx_re + Aref.A01*Hy_re + Aref.A02*Hz_re;
            Hwork[index+1] = Aref.A00*Hx_im + Aref.A01*Hy_im + Aref.A02*Hz_im;
            Hwork[index+2] = Aref.A01*Hx_re + Aref.A11*Hy_re + Aref.A12*Hz_re;
            Hwork[index+3] = Aref.A01*Hx_im + Aref.A11*Hy_im + Aref.A12*Hz_im;
            Hwork[index+4] = Aref.A02*Hx_re + Aref.A12*Hy_re + Aref.A22*Hz_re;
            Hwork[index+5] = Aref.A02*Hx_im + Aref.A12*Hy_im + Aref.A22*Hz_im;
          }
        }
      }

      ffty->InverseFFT(Hwork);

      // Copy data from Hwork back into Hxfrm
      for(j=0;j<rdimy;++j) {
        const OC_INDEX windex = j*Hwstride;
        const OC_INDEX hindex = j*jstride + ix*istride;
        memcpy(Hxfrm+hindex,Hwork+windex,Hcopy_line_size);
      }

    }
  }
}

class _Oxs_EddyFFTzConvolveThread : public Oxs_ThreadRunObj {
public:
  OXS_FFT_REAL_TYPE* Hxfrm;
  const Oxs_Eddy::A_coefs* A;

  Oxs_Eddy::Oxs_FFTLocker_Info locker_info;
  Oxs_Eddy::Oxs_FFTLocker* locker;

  static _Oxs_EddyJobControl job_control;

  OC_INDEX thread_count;
  OC_INDEX cdimx,cdimy,cdimz;
  OC_INDEX adimx,adimy,adimz;
  OC_INDEX rdimz;
  OC_INDEX embed_block_size;
  OC_INDEX jstride,ajstride;
  OC_INDEX kstride,akstride;
  _Oxs_EddyFFTzConvolveThread()
    : Hxfrm(0), A(0), locker(0),
      thread_count(0),
      cdimx(0),cdimy(0),cdimz(0),
      adimx(0),adimy(0),adimz(0),rdimz(0),
      embed_block_size(0),
      jstride(0),ajstride(0),kstride(0),akstride(0) {}
  void Cmd(int threadnumber, void* data);
};

_Oxs_EddyJobControl _Oxs_EddyFFTzConvolveThread::job_control;

void _Oxs_EddyFFTzConvolveThread::Cmd(int /* threadnumber */, void* /* data */)
{
  // Thread local storage
  if(!locker) {
    Oxs_ThreadMapDataObject* foo = local_locker.GetItem(locker_info.name);
    if(!foo) {
      // Oxs_FFTLocker object not constructed
      foo = new Oxs_Eddy::Oxs_FFTLocker(locker_info);
      local_locker.AddItem(locker_info.name,foo);
    }
    locker = dynamic_cast<Oxs_Eddy::Oxs_FFTLocker*>(foo);
    if(!locker) {
      Oxs_ThreadError::SetError(String("Error in"
         "_Oxs_EddyFFTzConvolveThread::Cmd(): locker downcast failed."));
      return;
    }
  }
  Oxs_FFTStrided* const fftz = &(locker->fftz);

  // Hwork:  Data is copied from Hxfrm into and out of this space
  // on each m increment.  Hwork1 shadows the active j>=0 block of
  // Hxfrm, Hwork2 the j<0 block.
  const OC_INDEX Hwstride = ODTV_VECSIZE*ODTV_COMPLEXSIZE*embed_block_size;
  OXS_FFT_REAL_TYPE* const Hwork1 = locker->fftz_Hwork;
  OXS_FFT_REAL_TYPE* const Hwork2 = Hwork1 + Hwstride * cdimz;

  if(locker->A_copy == NULL) {
    size_t asize = static_cast<size_t>(adimx*adimy*adimz);
    locker->A_copy_size = asize;
    locker->A_copy = static_cast<Oxs_Eddy::A_coefs*>
      (Oc_AllocThreadLocal(asize*sizeof(Oxs_Eddy::A_coefs)));
    memcpy(locker->A_copy,A,asize*sizeof(Oxs_Eddy::A_coefs));
  }
  const Oxs_Eddy::A_coefs* const Acopy = locker->A_copy;

  // Adjust fftz to use Hwork
  fftz->AdjustInputDimensions(rdimz,Hwstride,
                              ODTV_VECSIZE*embed_block_size);

  while(1) {
    OC_INDEX i,j,k;

    OC_INDEX jstart,jstop;
    job_control.ClaimJob(jstart,jstop);
    if(jstart>=jstop) break;

    for(j=jstart;j<jstop;++j) {
      // j>=0
      const OC_INDEX  jindex = j*jstride;
      const OC_INDEX ajindex = j*ajstride;

      const OC_INDEX j2 = cdimy - j;
      const OC_INDEX j2index = j2*jstride;

      fftz->AdjustArrayCount(ODTV_VECSIZE*embed_block_size);
      for(OC_INDEX m=0;m<cdimx;m+=embed_block_size) {

        // Do one block of forward z-direction transforms
        OC_INDEX istop_tmp = m + embed_block_size;
        if(embed_block_size>cdimx-m) {
          // Partial block
          fftz->AdjustArrayCount(ODTV_VECSIZE*(cdimx-m));
          istop_tmp = cdimx;
        }
        const OC_INDEX istop = istop_tmp;

        // Copy data into Hwork
        const size_t Hcopy_line_size
          = static_cast<size_t>(ODTV_COMPLEXSIZE*ODTV_VECSIZE*(istop-m))
          *sizeof(OXS_FFT_REAL_TYPE);
        for(k=0;k<rdimz;++k) {
          const OC_INDEX windex = k*Hwstride;
          const OC_INDEX h1index = k*kstride + m*ODTV_COMPLEXSIZE*ODTV_VECSIZE
            + jindex;
          memcpy(Hwork1+windex,Hxfrm+h1index,Hcopy_line_size);
        }
        if(adimy<=j2 && j2<cdimy) {
          for(k=0;k<rdimz;++k) {
            const OC_INDEX windex = k*Hwstride;
            const OC_INDEX h2index = k*kstride + m*ODTV_COMPLEXSIZE*ODTV_VECSIZE
              + j2index;
            memcpy(Hwork2+windex,Hxfrm+h2index,Hcopy_line_size);
          }
        }

        fftz->ForwardFFT(Hwork1);
        if(adimy<=j2 && j2<cdimy) {
          fftz->ForwardFFT(Hwork2);
        }

        // Do matrix-vector multiply ("convolution") for block
        for(k=0;k<adimz;++k) {
          // k>=0
          const OC_INDEX windex = k*Hwstride;
          const OC_INDEX akindex = ajindex + k*akstride;
          for(i=m;i<istop;++i) {
            const Oxs_Eddy::A_coefs& Aref = Acopy[akindex+i];
            const OC_INDEX index = ODTV_COMPLEXSIZE*ODTV_VECSIZE*(i-m)+windex;
            {
              OXS_FFT_REAL_TYPE Hx_re = Hwork1[index];
              OXS_FFT_REAL_TYPE Hx_im = Hwork1[index+1];
              OXS_FFT_REAL_TYPE Hy_re = Hwork1[index+2];
              OXS_FFT_REAL_TYPE Hy_im = Hwork1[index+3];
              OXS_FFT_REAL_TYPE Hz_re = Hwork1[index+4];
              OXS_FFT_REAL_TYPE Hz_im = Hwork1[index+5];

              Hwork1[index]   = Aref.A00*Hx_re + Aref.A01*Hy_re + Aref.A02*Hz_re;
              Hwork1[index+1] = Aref.A00*Hx_im + Aref.A01*Hy_im + Aref.A02*Hz_im;
              Hwork1[index+2] = Aref.A01*Hx_re + Aref.A11*Hy_re + Aref.A12*Hz_re;
              Hwork1[index+3] = Aref.A01*Hx_im + Aref.A11*Hy_im + Aref.A12*Hz_im;
              Hwork1[index+4] = Aref.A02*Hx_re + Aref.A12*Hy_re + Aref.A22*Hz_re;
              Hwork1[index+5] = Aref.A02*Hx_im + Aref.A12*Hy_im + Aref.A22*Hz_im;
            }
            if(adimy<=j2 && j2<cdimy) {
              OXS_FFT_REAL_TYPE Hx_re = Hwork2[index];
              OXS_FFT_REAL_TYPE Hx_im = Hwork2[index+1];
              OXS_FFT_REAL_TYPE Hy_re = Hwork2[index+2];
              OXS_FFT_REAL_TYPE Hy_im = Hwork2[index+3];
              OXS_FFT_REAL_TYPE Hz_re = Hwork2[index+4];
              OXS_FFT_REAL_TYPE Hz_im = Hwork2[index+5];

              // Flip signs on a01 and a12 as compared to the j>=0
              // case because a01 and a12 are odd in y.
              Hwork2[index]   =  Aref.A00*Hx_re - Aref.A01*Hy_re + Aref.A02*Hz_re;
              Hwork2[index+1] =  Aref.A00*Hx_im - Aref.A01*Hy_im + Aref.A02*Hz_im;
              Hwork2[index+2] = -Aref.A01*Hx_re + Aref.A11*Hy_re - Aref.A12*Hz_re;
              Hwork2[index+3] = -Aref.A01*Hx_im + Aref.A11*Hy_im - Aref.A12*Hz_im;
              Hwork2[index+4] =  Aref.A02*Hx_re - Aref.A12*Hy_re + Aref.A22*Hz_re;
              Hwork2[index+5] =  Aref.A02*Hx_im - Aref.A12*Hy_im + Aref.A22*Hz_im;
            }
          }
        }
        for(k=adimz;k<cdimz;++k) {
          // k<0
          const OC_INDEX windex = k*Hwstride;
          const OC_INDEX akindex = ajindex + (cdimz-k)*akstride;
          for(i=m;i<istop;++i) {
            const Oxs_Eddy::A_coefs& Aref = Acopy[akindex+i];
            const OC_INDEX index = ODTV_COMPLEXSIZE*ODTV_VECSIZE*(i-m)+windex;
            {
              OXS_FFT_REAL_TYPE Hx_re = Hwork1[index];
              OXS_FFT_REAL_TYPE Hx_im = Hwork1[index+1];
              OXS_FFT_REAL_TYPE Hy_re = Hwork1[index+2];
              OXS_FFT_REAL_TYPE Hy_im = Hwork1[index+3];
              OXS_FFT_REAL_TYPE Hz_re = Hwork1[index+4];
              OXS_FFT_REAL_TYPE Hz_im = Hwork1[index+5];

              // Flip signs on a02 and a12 as compared to the k>=0, j>=0 case
              // because a02 and a12 are odd in z.
              Hwork1[index]   =  Aref.A00*Hx_re + Aref.A01*Hy_re - Aref.A02*Hz_re;
              Hwork1[index+1] =  Aref.A00*Hx_im + Aref.A01*Hy_im - Aref.A02*Hz_im;
              Hwork1[index+2] =  Aref.A01*Hx_re + Aref.A11*Hy_re - Aref.A12*Hz_re;
              Hwork1[index+3] =  Aref.A01*Hx_im + Aref.A11*Hy_im - Aref.A12*Hz_im;
              Hwork1[index+4] = -Aref.A02*Hx_re - Aref.A12*Hy_re + Aref.A22*Hz_re;
              Hwork1[index+5] = -Aref.A02*Hx_im - Aref.A12*Hy_im + Aref.A22*Hz_im;
            }
            if(adimy<=j2 && j2<cdimy) {
              OXS_FFT_REAL_TYPE Hx_re = Hwork2[index];
              OXS_FFT_REAL_TYPE Hx_im = Hwork2[index+1];
              OXS_FFT_REAL_TYPE Hy_re = Hwork2[index+2];
              OXS_FFT_REAL_TYPE Hy_im = Hwork2[index+3];
              OXS_FFT_REAL_TYPE Hz_re = Hwork2[index+4];
              OXS_FFT_REAL_TYPE Hz_im = Hwork2[index+5];

              // Flip signs on a01 and a02 as compared to the k>=0, j>=0 case
              // because a01 is odd in y and even in z,
              //     and a02 is odd in z and even in y.
              // No change to a12 because it is odd in both y and z.
              Hwork2[index]   =  Aref.A00*Hx_re - Aref.A01*Hy_re - Aref.A02*Hz_re;
              Hwork2[index+1] =  Aref.A00*Hx_im - Aref.A01*Hy_im - Aref.A02*Hz_im;
              Hwork2[index+2] = -Aref.A01*Hx_re + Aref.A11*Hy_re + Aref.A12*Hz_re;
              Hwork2[index+3] = -Aref.A01*Hx_im + Aref.A11*Hy_im + Aref.A12*Hz_im;
              Hwork2[index+4] = -Aref.A02*Hx_re + Aref.A12*Hy_re + Aref.A22*Hz_re;
              Hwork2[index+5] = -Aref.A02*Hx_im + Aref.A12*Hy_im + Aref.A22*Hz_im;
            }
          }
        }
        // Do inverse z-direction transforms for block
        fftz->InverseFFT(Hwork1);
        if(adimy<=j2 && j2<cdimy) {
          fftz->InverseFFT(Hwork2);
        }

        // Copy data out of Hwork
        for(k=0;k<rdimz;++k) {
          const OC_INDEX windex = k*Hwstride;
          const OC_INDEX h1index = k*kstride + m*ODTV_COMPLEXSIZE*ODTV_VECSIZE
            + jindex;
          memcpy(Hxfrm+h1index,Hwork1+windex,Hcopy_line_size);
        }
        if(adimy<=j2 && j2<cdimy) {
          for(k=0;k<rdimz;++k) {
            const OC_INDEX windex = k*Hwstride;
            const OC_INDEX h2index = k*kstride + m*ODTV_COMPLEXSIZE*ODTV_VECSIZE
              + j2index;
            memcpy(Hxfrm+h2index,Hwork2+windex,Hcopy_line_size);
          }
        }
      }
    }
  }
}


void Oxs_Eddy::ComputeEnergy
(const Oxs_SimState& state,
 Oxs_ComputeEnergyData& oced
 ) const
{
    const OC_INDEX size = state.mesh->Size();
  // (Re)-initialize mesh coefficient array if mesh has changed.
  if(mesh_id != state.mesh->Id()) {
    mesh_id = 0; // Safety
    const Oxs_Mesh* mesh = state.mesh;
    FillCoefficientArrays(mesh);
    mp_init->FillMeshValue(mesh,mp);
    J_init->FillMeshValue(mesh,J);
    mesh_id = state.mesh->Id();

    for(OC_INDEX i=1 ; i < size; i++) {
         J[i] *= -1.;//*MU0;
       }
    }

    // Take care of time dependence -> multiply Jmult wit J indexwise
    OC_REAL8m Jmult = EvaluateJProfileScript(state.stage_number,
                          state.stage_elapsed_time,
                          state.stage_start_time+state.stage_elapsed_time);

    J_scratch.AdjustSize(state.mesh);

    for(OC_INDEX i=1; i< size; i++) {
         J_scratch[i]  = J[i];
         J_scratch[i] *= Jmult; 
    }

// For debugging purposes -> time profile
//    if(state.stage_number == 10) {
//       OC_CHAR buff[1042];
//       Oc_Snprintf(buff,sizeof(buff),"Jmult: %f",Jmult) ;
//       throw Oxs_ExtError(this,buff);
//      }


  // Fill Jtemp with mp[]*J[].  The plan is to eventually
  // roll this step into the forward FFT routine.
  assert(rdimx*rdimy*rdimz == mp.Size());

  const OC_INDEX rxdim = ODTV_VECSIZE*rdimx;
  const OC_INDEX cxdim = ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx;
  const OC_INDEX rxydim = rxdim*rdimy;
  const OC_INDEX cxydim = cxdim*cdimy;

  // Calculate x- and y-axis FFTs of Jtemp.
  OC_INT4m ithread;
#if REPORT_TIME
  fftxforwardtime.Start();
#endif // REPORT_TIME
  {
    vector<_Oxs_EddyFFTxThread> fftx_thread;
    fftx_thread.resize(MaxThreadCount);

    _Oxs_EddyFFTxThread::job_control.Init(rdimy*rdimz,MaxThreadCount,8);

    for(ithread=0;ithread<MaxThreadCount;++ithread) {
      fftx_thread[ithread].mp = &mp;
      fftx_thread[ithread].J = &J_scratch;
      fftx_thread[ithread].rarr = Jtemp;
      fftx_thread[ithread].carr = Hxfrm;
      fftx_thread[ithread].locker_info.Set(rdimx,rdimy,rdimz,
                                           cdimx,cdimy,cdimz,
                                           embed_block_size,
                                           MakeLockerName());
      fftx_thread[ithread].mp_xdim  = rdimx;
      fftx_thread[ithread].mp_xydim = rdimx*rdimy;

      fftx_thread[ithread].j_dim = rdimy;
      fftx_thread[ithread].j_rstride = rxdim;
      fftx_thread[ithread].j_cstride = cxdim;

      fftx_thread[ithread].k_rstride = rxydim;
      fftx_thread[ithread].k_cstride = cxydim;

      fftx_thread[ithread].jk_max = rdimy*rdimz;

      fftx_thread[ithread].direction = _Oxs_EddyFFTxThread::FORWARD;
      if(ithread>0) threadtree.Launch(fftx_thread[ithread],0);
    }
    threadtree.LaunchRoot(fftx_thread[0],0);
  }
#if REPORT_TIME
  fftxforwardtime.Stop();
#endif // REPORT_TIME


  if(cdimz<2) {
#if REPORT_TIME
    convtime.Start();
#endif // REPORT_TIME
    {
      vector<_Oxs_EddyFFTyConvolveThread> ffty_thread;
      ffty_thread.resize(MaxThreadCount);

      _Oxs_EddyFFTyConvolveThread::job_control.Init(cdimx,MaxThreadCount,1);

      for(ithread=0;ithread<MaxThreadCount;++ithread) {
        ffty_thread[ithread].Hxfrm = Hxfrm;
        ffty_thread[ithread].A = A;
        ffty_thread[ithread].locker_info.Set(rdimx,rdimy,rdimz,
                                             cdimx,cdimy,cdimz,
                                             embed_block_size,
                                             MakeLockerName());
        ffty_thread[ithread].embed_block_size = embed_block_size;
        ffty_thread[ithread].jstride = ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx;
        ffty_thread[ithread].ajstride = adimx;
        ffty_thread[ithread].i_dim = cdimx;
        ffty_thread[ithread].rdimy = rdimy;
        ffty_thread[ithread].adimy = adimy;
        ffty_thread[ithread].cdimy = cdimy;

        if(ithread>0) threadtree.Launch(ffty_thread[ithread],0);
      }
      threadtree.LaunchRoot(ffty_thread[0],0);
    }
#if REPORT_TIME
    convtime.Stop();
#endif // REPORT_TIME

  } else { // cdimz<2

#if REPORT_TIME
    fftyforwardtime.Start();
#endif // REPORT_TIME
    {
      vector<_Oxs_EddyFFTyThread> ffty_thread;
      ffty_thread.resize(MaxThreadCount);

      _Oxs_EddyFFTyThread::job_control.Init(cdimx*ODTV_VECSIZE*rdimz,
                                             MaxThreadCount,16);

      for(ithread=0;ithread<MaxThreadCount;++ithread) {
        ffty_thread[ithread].carr = Hxfrm;
        ffty_thread[ithread].locker_info.Set(rdimx,rdimy,rdimz,
                                             cdimx,cdimy,cdimz,
                                             embed_block_size,
                                             MakeLockerName());
        ffty_thread[ithread].k_stride = cxydim;
        ffty_thread[ithread].k_dim = rdimz;
        ffty_thread[ithread].i_dim = cdimx*ODTV_VECSIZE;
        ffty_thread[ithread].direction = _Oxs_EddyFFTyThread::FORWARD;
        if(ithread>0) threadtree.Launch(ffty_thread[ithread],0);
      }
      threadtree.LaunchRoot(ffty_thread[0],0);
    }
#if REPORT_TIME
    fftyforwardtime.Stop();
#endif // REPORT_TIME

    // Do z-axis FFTs with embedded "convolution" operations.
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
    {
      const OC_INDEX  jstride = ODTV_COMPLEXSIZE*ODTV_VECSIZE*cdimx;
      const OC_INDEX  kstride = jstride*cdimy;
      const OC_INDEX ajstride = adimx;
      const OC_INDEX akstride = ajstride*adimy;

      // Multi-thread
      vector<_Oxs_EddyFFTzConvolveThread> fftzconv;
      fftzconv.resize(MaxThreadCount);

      _Oxs_EddyFFTzConvolveThread::job_control.Init(adimy,MaxThreadCount,1);

      for(ithread=0;ithread<MaxThreadCount;++ithread) {
        fftzconv[ithread].Hxfrm = Hxfrm;
        fftzconv[ithread].A = A;
        fftzconv[ithread].locker_info.Set(rdimx,rdimy,rdimz,
                                          cdimx,cdimy,cdimz,
                                          embed_block_size,
                                          MakeLockerName());
        fftzconv[ithread].thread_count = MaxThreadCount;

        fftzconv[ithread].cdimx = cdimx;
        fftzconv[ithread].cdimy = cdimy;
        fftzconv[ithread].cdimz = cdimz;
        fftzconv[ithread].adimx = adimx;
        fftzconv[ithread].adimy = adimy;
        fftzconv[ithread].adimz = adimz;
        fftzconv[ithread].rdimz = rdimz;

        fftzconv[ithread].embed_block_size = embed_block_size;
        fftzconv[ithread].jstride = jstride;
        fftzconv[ithread].ajstride = ajstride;
        fftzconv[ithread].kstride = kstride;
        fftzconv[ithread].akstride = akstride;
        if(ithread>0) threadtree.Launch(fftzconv[ithread],0);
      }
      threadtree.LaunchRoot(fftzconv[0],0);
    }
#if REPORT_TIME
    convtime.Stop();
#endif // REPORT_TIME

    // Do inverse y- and x-axis FFTs, to complete transform back into
    // space domain.
#if REPORT_TIME
    fftyinversetime.Start();
#endif // REPORT_TIME
    {
      vector<_Oxs_EddyFFTyThread> ffty_thread;
      ffty_thread.resize(MaxThreadCount);

      _Oxs_EddyFFTyThread::job_control.Init(cdimx*ODTV_VECSIZE*rdimz,
                                             MaxThreadCount,16);

      for(ithread=0;ithread<MaxThreadCount;++ithread) {
        ffty_thread[ithread].carr = Hxfrm;
        ffty_thread[ithread].locker_info.Set(rdimx,rdimy,rdimz,
                                             cdimx,cdimy,cdimz,
                                             embed_block_size,
                                             MakeLockerName());
        ffty_thread[ithread].k_stride = cxydim;
        ffty_thread[ithread].k_dim = rdimz;
        ffty_thread[ithread].i_dim = cdimx*ODTV_VECSIZE;
        ffty_thread[ithread].direction = _Oxs_EddyFFTyThread::INVERSE;
        if(ithread>0) threadtree.Launch(ffty_thread[ithread],0);
      }
      threadtree.LaunchRoot(ffty_thread[0],0);
    }
#if REPORT_TIME
    fftyinversetime.Stop();
#endif // REPORT_TIME
  }  // cdimz<2
#if REPORT_TIME
  fftxinversetime.Start();
#endif // REPORT_TIME

  vector<_Oxs_EddyiFFTxDotThread> fftx_thread;
  fftx_thread.resize(MaxThreadCount);

  _Oxs_EddyiFFTxDotThread::job_control.Init(rdimy*rdimz,MaxThreadCount,8);

  for(ithread=0;ithread<MaxThreadCount;++ithread) {
    fftx_thread[ithread].carr = Hxfrm;

    fftx_thread[ithread].mp_ptr = &mp;
    fftx_thread[ithread].J_ptr   = &J_scratch;
    fftx_thread[ithread].oced_ptr = &oced;
    fftx_thread[ithread].locker_info.Set(rdimx,rdimy,rdimz,
                                         cdimx,cdimy,cdimz,
                                         embed_block_size,
                                         MakeLockerName());
    fftx_thread[ithread].rdimx = rdimx;
    fftx_thread[ithread].j_dim = rdimy;
    fftx_thread[ithread].j_rstride = rxdim;
    fftx_thread[ithread].j_cstride = cxdim;

    fftx_thread[ithread].k_rstride = rxydim;
    fftx_thread[ithread].k_cstride = cxydim;

    fftx_thread[ithread].jk_max = rdimy*rdimz;

    if(ithread>0) threadtree.Launch(fftx_thread[ithread],0);
  }
  threadtree.LaunchRoot(fftx_thread[0],0);

  Nb_Xpfloat tempsum = fftx_thread[0].energy_sum;
  for(ithread=1;ithread<MaxThreadCount;++ithread) {
    tempsum += fftx_thread[ithread].energy_sum;
  }
  oced.energy_sum = tempsum.GetValue() * state.mesh->Volume(0); // All
  /// cells have same volume in an Oxs_RectangularMesh.

#if REPORT_TIME
  fftxinversetime.Stop();
#endif // REPORT_TIME
}

OC_REAL8m Oxs_Eddy::EvaluateJProfileScript
(OC_UINT4m stage,
 OC_REAL8m stage_time,
 OC_REAL8m total_time) const
{
  if(!has_J_profile) return 1.0;

  int index;
  if((index = J_profile_opts[0].position)>=0) { // stage
    J_profile_cmd.SetCommandArg(index,stage);
  }
  if((index = J_profile_opts[1].position)>=0) { // stage_time
    J_profile_cmd.SetCommandArg(index,stage_time);
  }
  if((index = J_profile_opts[2].position)>=0) { // total_time
    J_profile_cmd.SetCommandArg(index,total_time);
  }

  J_profile_cmd.SaveInterpResult();
  J_profile_cmd.Eval();
  if(J_profile_cmd.GetResultListSize()!=1) {
    String msg
      = String("Return script value is not a single scalar: ")
      + J_profile_cmd.GetWholeResult();
    J_profile_cmd.RestoreInterpResult();
    throw Oxs_ExtError(this,msg.c_str());
  }
  OC_REAL8m result;
  J_profile_cmd.GetResultListItem(0,result);
  J_profile_cmd.RestoreInterpResult();

  return result;
}

#endif // OOMMF_THREADS
