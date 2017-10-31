/* FILE: imports.h                    -*-Mode: c++-*-
 *
 *	Directives to include header files from outside the extension.
 *
 * NOTICE: Plase see the file ../../LICENSE
 *
 * Last modified on: $Date: 2010-07-16 22:34:00 $
 * Last modified by: $Author: donahue $
 *
 * NOTE: Code which requires tcl.h or tk.h should *not* #include
 * tcl.h or tk.h directly, but should rather #include either this
 * file or oc.h (which #includes this file).  This is necessary
 * to insure that macros such as TCL_THREADS are properly defined
 * during the #include processing of tcl.h.
 */

#ifndef _OC_IMPORTS
#define _OC_IMPORTS

#include <stdio.h>
#include "ocport.h"

/*
 * Need 'extern "C"' to guarantee that typedefs of different function
 * types established in tcl.h, tk.h, or here are always declared with
 * C linkage.
 */
#ifdef __cplusplus
extern "C" {
/*
 * Also, when compiling with a C++ compiler, create macro definitions
 * in tcl.h and tk.h as would be suitable for a standard C compiler
 * using stdarg.h.
 */
#  ifndef HAS_STDARG
#    define HAS_STDARG
#  endif
#endif

/***********************************************************************
 * Defensive programming against varying contents of different versions 
 * of tcl.h
 **********************************************************************/
#if OOMMF_THREADS
# ifndef TCL_THREADS
#  define TCL_THREADS /* Enable Tcl threads support API */
# endif
#endif

/*
 * BIG KLUDGE: If we are building under Cygwin, then a straight
 * include of tcl.h will evaluate macros assuming a unix build.  But
 * if we are using a standard Windows build of Tcl under Cygwin, then
 * the Tcl libraries use Windows-style interfaces, not unix.  This can
 * lead to runtime errors; one example: the thread mainline return
 * type, "Tcl_ThreadCreateType", is void in unix but unsigned int in
 * Windows.  Detect this combination, and pre-define the __WIN32__
 * macro so that the declarations in tcl.h match the run-time library.
 */
#if OC_SYSTEM_TYPE==OC_UNIX && OC_SYSTEM_SUBTYPE==OC_CYGWIN \
    && OC_TCL_TYPE==OC_WINDOWS && !defined(__WIN32__)
# define __WIN32__
# include <tcl.h>
# undef __WIN32__
#else
# include <tcl.h>
#endif

/*
 * Verify version of tcl.h matches the Tcl version for which ocport.h was
 * configured.  Otherwise libraries selected for linking may not be
 * compatible with headers used for compiling, leading to a run-time error.
 */

#if ((TCL_MAJOR_VERSION != CONFIG_TCL_MAJOR_VERSION) \
    || (TCL_MINOR_VERSION != CONFIG_TCL_MINOR_VERSION))
#  error "tcl.h version mismatch"
#endif

/*
 * Safety that all macros have a definition (even in Tcl 7.5).
 */

#ifndef TCL_PATCH_LEVEL
#define TCL_PATCH_LEVEL CONFIG_TCL_PATCH_LEVEL
#define TCL_RELEASE_SERIAL CONFIG_TCL_RELEASE_SERIAL
#define TCL_RELEASE_LEVEL CONFIG_TCL_RELEASE_LEVEL
#endif

/*
 * Tcl 7.5p0 had no prototype for Tcl_Free()
 */
#if (TCL_MAJOR_VERSION == 7) && (TCL_MINOR_VERSION == 5) \
	&& (TCL_RELEASE_LEVEL == 2) && (TCL_RELEASE_LEVEL == 0)
#define Tcl_Free Tcl_Ckfree
#endif

/*
 * Provide a fake implementation of Tcl_GetStringResult for pre-8
 * Tcl libraries.
 */
#if (TCL_MAJOR_VERSION < 8)
#define Tcl_GetStringResult(interp) (interp->result)
#endif

/*
 * Provide Tcl_GetChannelHandle() to pre-8.0 interpreters
 * The TCL_STORAGE_CLASS stuff makes sure we don't make
 * the MSVC compiler think that the Tcl library will be
 * supplying the function.
 *
 * NB: 8.1aX, 8.1b1  not supported here.
 */

#if (TCL_MAJOR_VERSION < 8)
#undef TCL_STORAGE_CLASS
#define TCL_STORAGE_CLASS
EXTERN int Tcl_GetChannelHandle _ANSI_ARGS_((Tcl_Channel, int, ClientData*));
#undef TCL_STORAGE_CLASS
#define TCL_STORAGE_CLASS DLLIMPORT
#define Tcl_GetChannelOption(i,c,o,d) Tcl_GetChannelOption(c,o,d)
#endif

/*
 * In the tcl.h file released with Tcl 8.0.3 and Tcl 8.0.4, the
 * prototype for Tcl_AppInit is mistakenly declared with linkage
 * type __declspec(dllimport) on Windows platforms when tcl.h
 * is included in other apps, causing warnings from the MSVC
 * compiler.  The following declaration and macro work around
 * that bug.
 */

#if (TCL_MAJOR_VERSION == 8) && (TCL_MINOR_VERSION == 0) \
    && (TCL_RELEASE_LEVEL == 2) && (TCL_RELEASE_SERIAL > 2) \
    && (TCL_RELEASE_SERIAL < 5) && (OC_SYSTEM_TYPE == OC_WINDOWS)
Tcl_AppInitProc Oc_AppInit;
#define Tcl_AppInit Oc_AppInit
#endif

/*
 * Provide Tcl_PkgPresent() to pre-8.0.6 interpreters
 * The TCL_STORAGE_CLASS stuff makes sure we don't make
 * the MSVC compiler think that the Tcl library will be
 * supplying the function.
 *
 * NB: 8.1aX, 8.1b1  not supported here.
 */

#if ((TCL_MAJOR_VERSION < 8) || ((TCL_MAJOR_VERSION == 8) \
    && (TCL_MINOR_VERSION == 0) && (TCL_RELEASE_SERIAL < 6)))
#undef TCL_STORAGE_CLASS
#define TCL_STORAGE_CLASS
EXTERN char * Tcl_PkgPresent _ANSI_ARGS_((Tcl_Interp *, char *, char *, int));
#undef TCL_STORAGE_CLASS
#define TCL_STORAGE_CLASS DLLIMPORT
#endif

/*
 * A "declaration" for panic first appears in the 8.1b2 tcl.h .
 * (Really it's a macro definition to Tcl_Panic. "panic" is deprecated.)
 * For earlier Tcl versions, supply the missing prototype, but use the
 * panic() supplied by the Tcl library.
 *
 * NB: 8.1aX, 8.1b1  not supported here.
 */

#if ((TCL_MAJOR_VERSION < 8) \
    || ((TCL_MAJOR_VERSION == 8) && (TCL_MINOR_VERSION == 0)))
typedef void (Tcl_PanicProc) _ANSI_ARGS_(TCL_VARARGS(char *, format));
EXTERN void panic _ANSI_ARGS_(TCL_VARARGS(char *,format));
#endif

/*
 * The Tcl_SaveResult(), Tcl_RestoreResult(), and Tcl_DiscardResult()
 * routines and the Tcl_SavedResult struct were introduced in Tcl 8.1.
 * Provide a simplified form of them for pre-8.1 interpreters.
 *
 * Likewise for Tcl_UtfToExternalDString().
 *
 * NB: 8.1aX, 8.1b1  not supported here.
 */
#if ((TCL_MAJOR_VERSION < 8) \
    || ((TCL_MAJOR_VERSION == 8) && (TCL_MINOR_VERSION == 0)))
#undef TCL_STORAGE_CLASS
#define TCL_STORAGE_CLASS
typedef struct Tcl_SavedResult {
    char *result;
    Tcl_FreeProc *freeProc;
} Tcl_SavedResult;
EXTERN void Tcl_DiscardResult _ANSI_ARGS_((Tcl_SavedResult *));
EXTERN void Tcl_RestoreResult _ANSI_ARGS_((Tcl_Interp *, Tcl_SavedResult *));
EXTERN void Tcl_SaveResult _ANSI_ARGS_((Tcl_Interp *, Tcl_SavedResult *));
#undef TCL_STORAGE_CLASS
#define TCL_STORAGE_CLASS DLLIMPORT
#define Tcl_UtfToExternalDString(e, s, l, d) \
        (Tcl_DStringInit(d), Tcl_DStringAppend(d, s, l))
#endif

#ifndef CONST84
#    define CONST84
#    define OC_USING_CONST84 0
#else
#    define OC_USING_CONST84 1
#endif

#ifndef TCL_INTEGER_SPACE
#define TCL_INTEGER_SPACE 24
#endif

/*
 * Tcl_Tell() and Tcl_Seek() changed signatures in Tcl 8.4.  Define
 * a new type for their return type to bridge the difference.
 */
#ifdef TCL_WIDE_INT_TYPE
  typedef Tcl_WideInt Oc_SeekPos;
#else
  typedef int Oc_SeekPos;
#endif

/***********************************************************************
 * Defensive programming against varying contents of different versions 
 * of tk.h
 **********************************************************************/
#include <tk.h>

/*
 * Verify version of tk.h matches the Tk version for which ocport.h was
 * configured.  Otherwise libraries selected for linking may not be
 * compatible with headers used for compiling, leading to a run-time error.
 */

#if ((TK_MAJOR_VERSION != CONFIG_TK_MAJOR_VERSION) \
    || (TK_MINOR_VERSION != CONFIG_TK_MINOR_VERSION))
#  error "tk.h version mismatch"
#endif

/*
 * Safety that all macros have a definition (even in Tk 4.1).
 */

#ifndef TK_PATCH_LEVEL
#define TK_PATCH_LEVEL CONFIG_TK_PATCH_LEVEL
#define TK_RELEASE_SERIAL CONFIG_TK_RELEASE_SERIAL
#define TK_RELEASE_LEVEL CONFIG_TK_RELEASE_LEVEL
#endif

/*
 * Tk_SafeInit() was introduced with Tk 8.0
 */

#if (TK_MAJOR_VERSION < 8)
#define Tk_SafeInit ((Tcl_PackageInitProc *) NULL)
#endif

/*
 * Tk_InitConsoleChannels() and Tk_CreateConsoleWindow() first worked
 * in Tk 8.2.0.
 */
#if (TK_MAJOR_VERSION < 8) \
	|| ((TK_MAJOR_VERSION == 8) && (TK_MINOR_VERSION < 2)) \
	|| ((OC_SYSTEM_TYPE == OC_UNIX) && (OC_SYSTEM_SUBTYPE != OC_DARWIN))
#undef TCL_STORAGE_CLASS
#define TCL_STORAGE_CLASS
EXTERN void Tk_InitConsoleChannels _ANSI_ARGS_((Tcl_Interp *));
#define Tk_CreateConsoleWindow(interp) (0)
#undef TCL_STORAGE_CLASS
#define TCL_STORAGE_CLASS DLLIMPORT
#endif

/*
 * A pretty ugly hack to get around exit crashes in Tk 8.1.x on MS Windows.
 * The Tk Dll cleanup routines try to call Tcl_GetThreadData, but calls
 * to Tcl_* crash, unless the Tk DLL has called Tcl_InitStubs.  This is
 * normally done in Tk_Init, but since we're "statically" linking Tk,
 * we sometimes omit Tk_Init.  The routine Tk_Main also calls Tcl_InitStubs
 * from within the Tk DLL, so if we always call Tk_Main, we're safe.
 */
#if (OC_SYSTEM_TYPE == OC_WINDOWS) && (TK_MAJOR_VERSION == 8) \
        && (TK_MINOR_VERSION == 1) && (TK_RELEASE_LEVEL == 2) \
        && (TK_RELEASE_SERIAL < 2)
#define Tcl_Main Tk_Main
#endif

/* End includes */     /* Optional directive to pimake */

#ifdef __cplusplus
}	/* end of extern "C" */
#endif

/*
 * Utility routines
 */


//////////////////////////////////////////////////////////////////////////
// Random number generator routines.  The Oc_Random class (with
// Oc_RandomState helper) is an implementation of the GLIBC random()
// function.  This class provides the default functionality of the
// OMF_RANDOM macro declared in ocport.h
//    The Oc_UnifRand() function returns a random value in [0,1] with
// unif. distrib.  It is a wrapper for OMF_RANDOM.  It may be
// (re)initialized by calling Oc_Srand() or Oc_Srand(seed).  In the
// first case the seed is determined by sampling the system clock.

class Oc_RandomState {
  // Support class for Oc_Random
public:
  Oc_RandomState() : ihead(-1) {
    Init(1); // Default seed
  }
private:
  friend class Oc_Random;

  enum { SIZE = 32, SEP  =  3, WARMUP = 310 };

  int ihead;
  OC_UINT4m arr[SIZE+SEP];

  void Init(OC_UINT4m seed);
  OC_UINT4m Step();
};

class Oc_Random {
  // An implementation of glibc random(), based on notes
  // by Peter Selinger, "The GLIBC random number generator,"
  // 4-Jan-2007.
  // NOTE: The generator state is stored in a static variable.
  //       This means that the state is shared program-wide,
  //       across threads.  Mutexes protect against re-entrancy
  //       problems, but if multiple threads access Oc_Random
  //       then the results can vary from run-to-run, depending
  //       on the relative access order between threads.
public:
#if OOMMF_THREADS
  static Tcl_Mutex random_state_mutex;  // Thread-safe hack.
#endif
  static void Srandom(OC_UINT4m seed) {
    Tcl_MutexLock(&random_state_mutex);
    state.Init(seed);
    Tcl_MutexUnlock(&random_state_mutex);
  }
  static OC_INT4m Random() {
    Tcl_MutexLock(&random_state_mutex);
    OC_UINT4m step_result = state.Step();
    Tcl_MutexUnlock(&random_state_mutex);
    return static_cast<OC_INT4m>(step_result>>1);
  }
  static OC_INT4m MaxValue() { return 0x7fffffff; }
  /// NB: This function is referenced by procs.tcl to build
  ///     ocport.h.  Changes to MaxValue() may need reflection
  ///     there as well.
private:
  static Oc_RandomState state;
};

void Oc_Srand();
void Oc_Srand(unsigned int seed);
double Oc_UnifRand();
/*
 * Oc_UnifRand() random value in [0,1] with unif. distrib.
 * The random number generator may be (re)initialized by
 * calling Oc_Srand() or Oc_Srand(seed).  In the first case
 * the seed is determined by sampling the system clock.
 */

/* Utility math commands. */
double Oc_Polynomial(double x, double *coef, int N);
double Oc_Erf(double x);
double Oc_Hypot(double x,double y);

/*
 * Extra goodies for Tcl's expr command.  The Tcl_CmdProc version
 * is intended for use at the script level when setting up slave
 * interpreters.  The main interpreter has this done for it during
 * Oc_Init().  Slave interpreters get these by redefinition of
 * "interp" in custom.tcl.
 */
double Oc_Atan2(double y,double x);
void Oc_AddTclExprExtensions(Tcl_Interp* interp);

#ifdef __cplusplus
extern "C" {
#endif

/* Wrappers for the above. */
Tcl_CmdProc OcSrand;
Tcl_CmdProc OcUnifRand;
#if OC_TCL_TYPE==OC_WINDOWS && defined(__CYGWIN__)
Tcl_CmdProc OcCygwinChDir;
#endif //  OC_WINDOWS && __CYGWIN__
Tcl_MathProc Oc_TclWrappedAtan2;
Tcl_MathProc Oc_Exp;
Tcl_MathProc Oc_Pow;

Tcl_CmdProc OcAddTclExprExtensions;

#ifdef __cplusplus
}	/* end of extern "C" */
#endif

#endif /* _OC_IMPORTS */
