/* FILE: evoc.cc                   -*-Mode: c++-*-
 *
 * Stuff that should be moved eventually into the oc extension.
 *
 * NOTICE: Please see the file ../../LICENSE
 *
 * Last modified on: $Date: 2010-07-16 22:33:55 $
 * Last modified by: $Author: donahue $
 */

// Note: We need to include evoc.h before testing the
// OC_SYSTEM_TYPE macro, because OC_SYSTEM_TYPE is declared
// in oommf/pkg/oc/<platform>/ocport.h, which is included
// in a round-about fashion by evoc.h
#include "evoc.h"

#include <string.h>

#if (OC_SYSTEM_TYPE == OC_UNIX)
#include <sys/time.h>  // Some of these may be OS dependent...
#include <sys/times.h>
#endif // OC_SYSTEM_TYPE

#include <time.h>

#include "errhandlers.h"

/* End includes */     

int Verbosity=10; // Default value for Verbosity

ClassDoc::ClassDoc(const OC_CHAR *new_classname,const OC_CHAR *new_maintainer,
		   const OC_CHAR *new_revision,const OC_CHAR *new_revdate)
{
  classname  = new_classname;
  maintainer = new_maintainer;
  revision   = new_revision;
  revdate    = new_revdate;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////// Tcl/Tk Interface ////////////////////////////
//////////////////////////////////////////////////////////////////////////

// C-pointer <--> ascii string conversion class.
//   Currently this routine maps things through %lx, but this
// may change in the future to use a hash table.  (I'm assuming
// that void * can be converted to an unsigned long without
// loss of data.  There is a check on this in PtrToAscii for safety.)

const ClassDoc Omf_AsciiPtr::class_doc("Omf_AsciiPtr",
	      "Michael J. Donahue (michael.donahue@nist.gov)",
	      "1.0.0","19-Sep-1997");

const size_t Omf_AsciiPtr::ascii_string_width(2*sizeof(unsigned long)+8);
/// Each byte goes to 2 hex digits (hence 2*), +1 for the trailing
/// null, +2 in case a "0x" is prepended, +5 for safety.

void Omf_AsciiPtr::PtrToAscii(const void* ptr,char* buf)
{
#define MEMBERNAME "PtrToAscii(const void* ptr,char* buf)"
  ptr_int_hack cvt;
  cvt.lval = 0;
  cvt.ptr = const_cast<void*>(ptr);
  sprintf(buf,"%lx",cvt.lval);
  if(strlen(buf)+1>GetAsciiSize()) {
    FatalError(-1,STDDOC,"Buffer overflow; "
	       "Omf_AsciiPtr::ascii_string_width too small");
  }

  // Check readback
  if(AsciiToPtr(buf)!=ptr) {
    char errmsg[1024];
    void* checkval = AsciiToPtr(buf);
    Oc_Snprintf(errmsg,sizeof(errmsg),
                "Ptr -> lx-format -> ptr conversion error;"
                " Input %lX, export %lX, diff=%ld",
                (unsigned long)ptr,(unsigned long)checkval,
                (long)((char*)checkval-(char*)ptr));
    FatalError(-1,STDDOC,errmsg);
  }

  return;
#undef MEMBERNAME
}

void* Omf_AsciiPtr::AsciiToPtr(const char* buf)
{ // NOTE: ONLY strings produced by Omf_AsciiPtr::PtrToAscii
  //  should ever be sent to this routine.
#define MEMBERNAME "AsciiToPtr(const char* buf,void* &ptr)"
  ptr_int_hack cvt;
  cvt.ptr = 0;
  sscanf(buf,"%lx",&cvt.lval);
  Oc_Nop(&cvt);  // Some compilers (cough, Solaris CC, cough) at
  /// high optimization swallow the union.  Call Oc_Nop to force
  /// memory write.
  return cvt.ptr;
#undef MEMBERNAME
}

//////////////////////////////////////////////////////////////////////////
// Oc_TimeVal class

const ClassDoc Oc_TimeVal::class_doc("Oc_TimeVal",
	      "Michael J. Donahue (michael.donahue@nist.gov)",
	      "1.0.0","May-1998");

void Oc_TimeVal::Print(FILE* fptr)
{ // For debugging
#if (TCL_MAJOR_VERSION < 8) || \
     (TCL_MAJOR_VERSION == 8 && TCL_MINOR_VERSION < 4)
  // Tcl_WideUInt type and friends first appear in Tcl 8.4.
  // For older versions, print as double, which provides 52
  // bits in the mantissa.
  fprintf(fptr," ticks_per_second=%.0f\n",
          OC_TIMEVAL_TO_DOUBLE(ticks_per_second));
  fprintf(fptr,"        max_ticks=%.0f\n",
          OC_TIMEVAL_TO_DOUBLE(max_ticks));
  fprintf(fptr,"            ticks=%.0f\n",
          OC_TIMEVAL_TO_DOUBLE(ticks));
  fprintf(fptr,"         overflow=%.0f\n",
          OC_TIMEVAL_TO_DOUBLE(overflow));
  fflush(fptr);
#else // Tcl 8.4 or later
  // Note: Tcl_WideUInt type is unsigned integer at least 64-bits
  //  wide.  The TCL_LL_MODIFIER macro is magic string to get proper
  //  format modifier.  These work across a variety of C compilers
  //  supported by Tcl.  (See tcl.h for details.)  In principle C++
  //  may be different in this regard (for example, "long long" types
  //  are part of C99 spec but not officially part of C++ as of this
  //  date (May 2008)), but in practice C++ compilers are largely built
  //  on top of C compilers, and share much of C99 as extensions.
  fprintf(fptr," ticks_per_second=%" TCL_LL_MODIFIER "u\n",
          static_cast<Tcl_WideUInt>(ticks_per_second));
  fprintf(fptr,"        max_ticks=%" TCL_LL_MODIFIER "u\n",
          static_cast<Tcl_WideUInt>(max_ticks));
  fprintf(fptr,"            ticks=%" TCL_LL_MODIFIER "u\n",
          static_cast<Tcl_WideUInt>(ticks));
  fprintf(fptr,"         overflow=%" TCL_LL_MODIFIER "u\n",
          static_cast<Tcl_WideUInt>(overflow));
  fflush(fptr);
#endif
}

OC_BOOL Oc_TimeVal::IsValid() const
{
  if(ticks>max_ticks) return 0;
  return 1;
}

OC_BOOL AreCompatible(const Oc_TimeVal& time1,const Oc_TimeVal& time2)
{
  if(time1.ticks_per_second==time2.ticks_per_second &&
     time1.max_ticks==time2.max_ticks) return 1;
  return 0;
}

void Oc_TimeVal::Reset()
{
  ticks=0;
  overflow=0;
}

void Oc_TimeVal::Reset(OC_TIMEVAL_TICK_TYPE _ticks_per_second,
		       OC_TIMEVAL_TICK_TYPE _max_ticks)
{
  ticks_per_second=_ticks_per_second;
  max_ticks=_max_ticks;
  Reset();
}

Oc_TimeVal::Oc_TimeVal()
{
  Reset(1,(OC_TIMEVAL_TICK_TYPE)-1);
}

Oc_TimeVal::Oc_TimeVal(const Oc_TimeVal &time)
{
  Reset(time.ticks_per_second,time.max_ticks);
  ticks=time.ticks;
  overflow=time.overflow;
}

Oc_TimeVal::Oc_TimeVal(OC_TIMEVAL_TICK_TYPE _ticks_per_second,
		       OC_TIMEVAL_TICK_TYPE _max_ticks)
{
  Reset(_ticks_per_second,_max_ticks);
}

void Oc_TimeVal::SetTicks(OC_TIMEVAL_TICK_TYPE _ticks)
{
#define MEMBERNAME "SetTicks"
  if(_ticks>max_ticks)
    FatalError(-1,STDDOC,"Import _ticks bigger than max_ticks");
  ticks=_ticks;
#undef MEMBERNAME
}

Oc_TimeVal::operator double() const
{ // Returns time in seconds, in floating point.
  // Guard carefully against overflow on integer types
  double bigticks= OC_TIMEVAL_TO_DOUBLE(overflow)
    *(OC_TIMEVAL_TO_DOUBLE(max_ticks)+1.0);
  double smallticks=OC_TIMEVAL_TO_DOUBLE(ticks);
  return (bigticks+smallticks)/OC_TIMEVAL_TO_DOUBLE(ticks_per_second);
}

Oc_TimeVal& Oc_TimeVal::operator=(const Oc_TimeVal& time)
{
  Reset(time.ticks_per_second,time.max_ticks);
  ticks=time.ticks;
  overflow=time.overflow;
  return *this;
}

Oc_TimeVal& Oc_TimeVal::operator+=(const Oc_TimeVal& time)
{
#define MEMBERNAME "operator+="
  // Insure times use compatible bases
  if(!AreCompatible(*this,time))
    FatalError(-1,STDDOC,"Attempt to add incompatible Oc_TimeVal's\n"
	       "  *this is (%lu,%lu), time is (%lu,%lu)",
	       ticks_per_second,max_ticks,
	       time.ticks_per_second,time.max_ticks);
  overflow+=time.overflow;
  if( max_ticks-time.ticks >= ticks ) {
    ticks+=time.ticks;
  } else {
    overflow++;
    // Be careful to protect against overflow
    OC_TIMEVAL_TICK_TYPE tmpticks = max_ticks - time.ticks;
    ticks -= (tmpticks + 1);
  }
  return *this;
#undef MEMBERNAME
}

Oc_TimeVal& Oc_TimeVal::operator-=(const Oc_TimeVal& time)
{ // NOTE: Truncates to zero if time>*this.
#define MEMBERNAME "operator-="
  // Insure times use compatible bases
  if(!AreCompatible(*this,time))
    FatalError(-1,STDDOC,"Attempt to subtract incompatible Oc_TimeVal's");
  if(overflow<time.overflow) {
    Reset();
  } else {
    overflow-=time.overflow;
    if(ticks>=time.ticks) {
      ticks-=time.ticks;
    } else {
      if(overflow<1) Reset();
      else {
	overflow--;
        // Be careful to protect against overflow
        OC_TIMEVAL_TICK_TYPE tmpticks = max_ticks - time.ticks;
        ticks += tmpticks + 1;
      }
    }
  }
  return *this;
#undef MEMBERNAME
}


Oc_TimeVal operator+(const Oc_TimeVal& time1,const Oc_TimeVal& time2)
{ // This routine assumes time1 and time2 are valid.

  // Insure times use compatible bases
  if(!AreCompatible(time1,time2))
    PlainError(1,"Attempt to add incompatible Oc_TimeVal's");

  Oc_TimeVal time3(time1);
  time3.overflow+=time2.overflow;
  if( time3.max_ticks-time2.ticks >= time3.ticks)
    time3.ticks+=time2.ticks;
  else {
    time3.overflow++;
    // Be careful to protect against overflow
    OC_TIMEVAL_TICK_TYPE tmpticks = time3.max_ticks - time2.ticks;
    time3.ticks -= (tmpticks + 1);
  }

  return time3;
}

Oc_TimeVal operator-(const Oc_TimeVal& time1,const Oc_TimeVal& time2)
{ // This routine assumes time1 and time2 are valid.
  // If time1<time2, then returns zero time.

  // Insure times use compatible bases
  if(!AreCompatible(time1,time2))
    PlainError(1,"Attempt to subtract incompatible Oc_TimeVal's");

  Oc_TimeVal zero(time1);
  zero.ticks=0; zero.overflow=0;

  Oc_TimeVal time3(time1);
  if(time3.overflow<time2.overflow) return zero;
  time3.overflow-=time2.overflow;
  
  if(time3.ticks>=time2.ticks)
    time3.ticks-=time2.ticks;
  else {
    if(time3.overflow<1) return zero;
    time3.overflow--;
    // Be careful to protect against overflow
    OC_TIMEVAL_TICK_TYPE tmpticks = time3.max_ticks - time2.ticks;
    time3.ticks += tmpticks + 1;
  }

  return time3;
}

Oc_TimeVal GetBigger(const Oc_TimeVal& time1,const Oc_TimeVal& time2)
{ // This routine assumes time1 and time2 are valid.

  // Insure times use compatible bases
  if(!AreCompatible(time1,time2))
    PlainError(1,"Attempt to compare incompatible Oc_TimeVal's");

  if(time1.overflow > time2.overflow) return time1;
  if(time1.overflow < time2.overflow) return time2;
  if(time1.ticks > time2.ticks)       return time1;
  return time2;
}

Oc_TimeVal GetSmaller(const Oc_TimeVal& time1,const Oc_TimeVal& time2)
{ // This routine assumes time1 and time2 are valid.

  // Insure times use compatible bases
  if(!AreCompatible(time1,time2))
    PlainError(1,"Attempt to compare incompatible Oc_TimeVal's");

  if(time1.overflow < time2.overflow) return time1;
  if(time1.overflow > time2.overflow) return time2;
  if(time1.ticks < time2.ticks)       return time1;
  return time2;
}

//////////////////////////////////////////////////////////////////////////
// Oc_Times function
//   System independent replacement for the Unix times(2)
// function. Returns the cpu and elapsed times for the current
// process, relative to the first time this routine is called.  If
// this routine is called early in the process initialization, then
// the returned times will be effectively process times.
//   Resolution is system dependent.  This routine tries to correct for
// counter overflow and wrap-around by keeping track of the system
// tick count from the last call, and if the new return value is
// smaller than the last then a wrap-around of 1 period is assumed.
// This will not work properly if the time between calls to this
// function is larger than the wrap-around period.  For an 4-byte wide
// clock_t with CLOCKS_PER_SEC at 1024, the time to overflow is
// about 48.5 days (Windows NT).  With CLOCKS_PER_SEC at 1000000
// overflow time is just over 71 minutes (Linux/x86).  OTOH, an 8-byte
// wide clock_t with a nanosecond tick rate takes over 584 years to
// overflow.  NOTE 1: Even though CLOCKS_PER_SEC is 1e6 on
// Linux/x86, the granularity of all user-accessible clocks appears to
// be only 1/100 seconds (=1/CLK_TCK).  On Linux/AXP, CLK_TCK is 1000
// (and clock_t is 8 bytes wide). NOTE 2: The wrap around periods are
// generally determined as clock_t(-1).  This ASSUMES clock_t is an
// unsigned type (and 2's complement integer arithmetic).  We should
// probably put a check for this in varinfo.cc/ocport.h.
//   If suitable system timing call(s) can't be determined, then this
// routine will return 0('s).
//
// For Unix platforms, use the following macros:
//   HAS_TIMES
//   HAS_CLOCK
//   HAS_GETTIMEOFDAY
// If HAS_TIMES is defined, then Oc_Times will use the system times()
// command.  Otherwise, if HAS_CLOCK is defined, then clock() will be
// used to determine the cpu time; if HAS_GETTIMEOFDAY is defined, then
// gettimeofday() will be used to set the wall (elapsed) time.
// On Windows platforms, clock() and GetTickCount() are used.
// Alternately, define
//   NO_CLOCKS
// to have these routines always return 0.  This always overrides the
// other macros.

///// TEMPORARY MACROS, TO BE FIXED UP BY DGP ////////////////////
#if !defined(NO_CLOCKS) && (OC_SYSTEM_TYPE != OC_WINDOWS)
#  ifndef CLK_TCK
#    ifdef _SC_CLK_TCK
#      define CLK_TCK sysconf(_SC_CLK_TCK)
#    elif defined(HZ)
#      define CLK_TCK HZ
#    endif
#  endif
#  ifdef CLK_TCK
#    define HAS_TIMES
#  endif
#  ifdef CLOCKS_PER_SEC
#    define HAS_CLOCK
#  endif
#  define HAS_GETTIMEOFDAY
#endif
//////////////////////////////////////////////////////////////////

#ifdef NO_CLOCKS
void Oc_Times(Oc_TimeVal& cpu_time,Oc_TimeVal& wall_time)
{
  cpu_time.Reset();  wall_time.Reset();  // Return zeros
}
#else // NO_CLOCKS
void Oc_Times(Oc_TimeVal& cpu_time,Oc_TimeVal& wall_time,OC_BOOL reset)
{
  TCL_DECLARE_MUTEX(time_mutex);  // Thread-safe hack.  Ideally,
  /// this code should be re-written to better support thread safety.
  Tcl_MutexLock(&time_mutex);
  try {
    static OC_BOOL first_time=1;
    if(reset) {
      first_time = 1;
    }

#if (OC_SYSTEM_TYPE == OC_WINDOWS)
#ifdef HAS_GETPROCESSTIMES
    // Use GetProcessTimes if available, because the clock() function
    // appears to return wall time on Win2K.  Note: The tick interval
    // for GetProcessTimes is 100 ns.
    static Oc_TimeVal cpu_accum(10000000,(OC_TIMEVAL_TICK_TYPE)(-1));
    static Oc_TimeVal cpu_last(10000000,(OC_TIMEVAL_TICK_TYPE)(-1));
    static Oc_TimeVal cpu_now(10000000,(OC_TIMEVAL_TICK_TYPE)(-1));
    static Oc_TimeVal wall_accum(1000,(OC_TIMEVAL_TICK_TYPE)DWORD(-1));
    static Oc_TimeVal wall_last(1000,(OC_TIMEVAL_TICK_TYPE)DWORD(-1));
    static Oc_TimeVal wall_now(1000,(OC_TIMEVAL_TICK_TYPE)DWORD(-1));
    if(reset) {
      cpu_accum.Reset();   cpu_last.Reset();   cpu_now.Reset();
      wall_accum.Reset();  wall_last.Reset();  wall_now.Reset();
    }
    HANDLE current_process = GetCurrentProcess();
    FILETIME create_time,exit_time,kernel_time,user_time;
    if(GetProcessTimes(current_process,&create_time,&exit_time,
                       &kernel_time,&user_time)) {
      OC_TIMEVAL_TICK_TYPE kticks = 
        (static_cast<OC_TIMEVAL_TICK_TYPE>(kernel_time.dwHighDateTime)<<32)
        | static_cast<OC_TIMEVAL_TICK_TYPE>(kernel_time.dwLowDateTime);
      OC_TIMEVAL_TICK_TYPE uticks =
        (static_cast<OC_TIMEVAL_TICK_TYPE>(user_time.dwHighDateTime)<<32)
        | static_cast<OC_TIMEVAL_TICK_TYPE>(user_time.dwLowDateTime);
      cpu_now.ticks = kticks + uticks;
      cpu_now.overflow=cpu_last.overflow;
      if(cpu_now.ticks<cpu_last.ticks) cpu_now.overflow++;
    } else {
      cpu_now.ticks = 0;
      cpu_now.overflow = 0; // Punt
    }
    wall_now.ticks=GetTickCount();
    wall_now.overflow=wall_last.overflow;
    if(wall_now.ticks<wall_last.ticks) wall_now.overflow++;
#else // !HAS_GETPROCESSTIMES
    // Use clock() to get cpu time, GetTickCount() to get wall time.
    // GetTickCount() is in the Windows API; it returns number of ms
    // since Windows was started
    static Oc_TimeVal cpu_accum(CLOCKS_PER_SEC,
                                (OC_TIMEVAL_TICK_TYPE)clock_t(-1));
    static Oc_TimeVal cpu_last(CLOCKS_PER_SEC,
                               (OC_TIMEVAL_TICK_TYPE)clock_t(-1));
    static Oc_TimeVal cpu_now(CLOCKS_PER_SEC,
                              (OC_TIMEVAL_TICK_TYPE)clock_t(-1));
    static Oc_TimeVal wall_accum(1000,(OC_TIMEVAL_TICK_TYPE)DWORD(-1));
    static Oc_TimeVal wall_last(1000,(OC_TIMEVAL_TICK_TYPE)DWORD(-1));
    static Oc_TimeVal wall_now(1000,(OC_TIMEVAL_TICK_TYPE)DWORD(-1));
    if(reset) {
      cpu_accum.Reset();   cpu_last.Reset();   cpu_now.Reset();
      wall_accum.Reset();  wall_last.Reset();  wall_now.Reset();
    }
    cpu_now.ticks=clock();
    wall_now.ticks=GetTickCount();
    cpu_now.overflow=cpu_last.overflow;
    if(cpu_now.ticks<cpu_last.ticks) cpu_now.overflow++;
    wall_now.overflow=wall_last.overflow;
    if(wall_now.ticks<wall_last.ticks) wall_now.overflow++;
#endif // HAS_GETPROCESSTIMES
#else // OC_SYSTEM_TYPE != OC_WINDOWS
    // gettimeofday tends to have better resolution than times()
    // so use gettimeofday if possible.

#if defined(HAS_TIMES) && (!defined(HAS_GETTIMEOFDAY) || !defined(HAS_CLOCK))
    // We are going to use these below:
    static const long clktck = CLK_TCK;
    static struct tms buf;
#endif

# if defined(HAS_GETTIMEOFDAY)
    // Store microseconds in .tick, seconds in .overflow
    static Oc_TimeVal wall_accum(1000000,999999);
    static Oc_TimeVal wall_last(1000000,999999);
    static Oc_TimeVal wall_now(1000000,999999);
    if(reset) {
      wall_accum.Reset();  wall_last.Reset();  wall_now.Reset();
    }
    static struct timeval tv;
    gettimeofday(&tv,NULL);
    wall_now.overflow=tv.tv_sec;  // Seconds
    wall_now.ticks=tv.tv_usec;    // Microseconds
# elif defined(HAS_TIMES)
    // Use the Unix times(2) call
    static Oc_TimeVal wall_accum(clktck,(OC_TIMEVAL_TICK_TYPE)clock_t(-1));
    static Oc_TimeVal wall_last(clktck,(OC_TIMEVAL_TICK_TYPE)clock_t(-1));
    static Oc_TimeVal wall_now(clktck,(OC_TIMEVAL_TICK_TYPE)clock_t(-1));
    if(reset) {
      wall_accum.Reset();  wall_last.Reset();  wall_now.Reset();
    }
    wall_now.ticks = (OC_TIMEVAL_TICK_TYPE)times(&buf);
    wall_now.overflow=wall_last.overflow;
    if(wall_now.ticks<wall_last.ticks) wall_now.overflow++;
# else // !HAS_GETTIMEOFDAY
    static Oc_TimeVal wall_accum,wall_last,wall_now;
# endif // HAS_GETTIMEOFDAY

# ifdef HAS_CLOCK
    static Oc_TimeVal cpu_accum(CLOCKS_PER_SEC,
                                (OC_TIMEVAL_TICK_TYPE)clock_t(-1));
    static Oc_TimeVal cpu_last(CLOCKS_PER_SEC,
                               (OC_TIMEVAL_TICK_TYPE)clock_t(-1));
    static Oc_TimeVal cpu_now(CLOCKS_PER_SEC,
                              (OC_TIMEVAL_TICK_TYPE)clock_t(-1));
    if(reset) {
      cpu_accum.Reset();   cpu_last.Reset();   cpu_now.Reset();
    }
    cpu_now.ticks=clock();
    cpu_now.overflow=cpu_last.overflow;
    if(cpu_now.ticks<cpu_last.ticks) cpu_now.overflow++;
# elif defined(HAS_TIMES)
    // Use the Unix times(2) call
    static Oc_TimeVal cpu_accum(clktck,(OC_TIMEVAL_TICK_TYPE)clock_t(-1));
    static Oc_TimeVal cpu_last(clktck,(OC_TIMEVAL_TICK_TYPE)clock_t(-1));
    static Oc_TimeVal cpu_now(clktck,(OC_TIMEVAL_TICK_TYPE)clock_t(-1));
    if(reset) {
      cpu_accum.Reset();   cpu_last.Reset();   cpu_now.Reset();
    }
    cpu_now.ticks=(OC_TIMEVAL_TICK_TYPE)(buf.tms_utime+buf.tms_stime);
    cpu_now.overflow=cpu_last.overflow;
    if(cpu_now.ticks<cpu_last.ticks) cpu_now.overflow++;
# else // !HAS_CLOCK
    static Oc_TimeVal cpu_accum,cpu_last,cpu_now;
# endif // HAS_CLOCK

#endif // ...OC_SYSTEM_TYPE == OC_WINDOWS
    if(!first_time) {
      cpu_accum+=(cpu_now-cpu_last);
      wall_accum+=(wall_now-wall_last);
    } else {
      first_time=0;
    }
    cpu_last=cpu_now;     wall_last=wall_now;
    cpu_time=cpu_accum;   wall_time=wall_accum;
  } catch(...) {
    Tcl_MutexUnlock(&time_mutex);
    throw;
  }
  Tcl_MutexUnlock(&time_mutex);
}
#endif // NO_CLOCKS

// Tcl wrapper for Oc_Times
int OcTimes(ClientData,Tcl_Interp *interp,int argc,CONST84 char** argv)
{
  static char buf[256];
  Tcl_ResetResult(interp);
  if(argc<1 || argc>2) {
    Oc_Snprintf(buf,sizeof(buf),
		"wrong # args: should be \"%.100s ?reset?\"",argv[0]);
    Tcl_AppendResult(interp,buf,(char *)NULL);
    return TCL_ERROR;
  }

  OC_BOOL reset = 0;
  if(argc==2) {
    reset = atoi(argv[1]);
  }

  Oc_TimeVal cpu_time,wall_time;
  Oc_Times(cpu_time,wall_time,reset);

  Oc_Snprintf(buf,sizeof(buf),"%.17g %.17g",
	      double(cpu_time),double(wall_time));
  Tcl_AppendResult(interp,buf,(char *)NULL);

  return TCL_OK;
}
