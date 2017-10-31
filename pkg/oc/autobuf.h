/* FILE: autobuf.h                    -*-Mode: c++-*-
 *
 * This class was initially created to convert string literals into
 * the writable strings needed as arguments to old-versions of many
 * Tcl library calls.  The Oc_AutoBuf class is now used in many areas
 * of the code where a simple dynamic buffer is needed with automatic
 * and exception-safe cleanup.
 * 
 * NOTICE: Please see the file ../../LICENSE
 *
 * Last modified on: $Date: 2009-10-29 05:11:53 $
 * Last modified by: $Author: donahue $
 */

#ifndef _OC_AUTOBUF
#define _OC_AUTOBUF

#include <stddef.h>  // size_t definition

#include "imports.h" // OC_USING_CONST84 definition

/* End includes */     /* Optional directive to pimake */

class Oc_AutoBuf
{
private:
  static const size_t reclaim_bound;
  static const size_t reclaim_proportion;
  size_t bufsize; // Length of buffer pointed to by str. ( >= strlen(str)+1.)
  char *str;
public:
  ~Oc_AutoBuf() { delete[] str; }

  Oc_AutoBuf() : bufsize(0), str(NULL) { Dup(""); }
  Oc_AutoBuf(const char* const_str)
    : bufsize(0), str(NULL) { Dup(const_str); }

  Oc_AutoBuf(const Oc_AutoBuf& other);
  // Note: Proper copy constructor is needed for use of Oc_AutoBuf as
  // arguments or as members of arguments to C++ exception "throw"
  // commands.  It doesn't hurt to have a proper assignment operator
  // either.

#if (OC_SYSTEM_TYPE == OC_WINDOWS)
  // Specialized conversion function for Windows.  Converts from wide
  // char to char.  This is mainly needed in the case UNICODE is
  // defined, because then TCHAR is a wide char (TCHAR is a normal
  // char otherwise.)
  Oc_AutoBuf(const WCHAR* const_str);
#endif // WINDOWS

  Oc_AutoBuf& operator=(const Oc_AutoBuf& other);

  operator char*() const { return str; }
  char* GetStr() const { return str; }
  char* operator()(const char* const_str);
  void Dup(const char* const_str);
  void Append(const char* const_str);
  void SetLength(size_t new_length); // bufsize = length + 1
  size_t GetLength() { return bufsize-1; }
  const char& operator[](int i) const {
#ifndef NDEBUG
    if(i<0 || static_cast<size_t>(i) > bufsize) {
      char errmsg[]
	= "Out-of-bounds access in Oc_AutoBuf::operator[] const";
      panic(errmsg); // errmsg buffer used to support Tcl 8.3 and earlier,
      /// for which panic takes char* (as opposed to const char*)
    }
#endif
    return str[i];
  }
  char& operator[](int i) {
#ifndef NDEBUG
    if(i<0 || static_cast<size_t>(i) > bufsize) {
      char errmsg[]
	= "Out-of-bounds access in Oc_AutoBuf::operator[]";
      panic(errmsg); // errmsg buffer used to support Tcl 8.3 and earlier,
      /// for which panic takes char* (as opposed to const char*)
    }
#endif
    return str[i];
  }

}; 

// NOTE: OC_USING_CONST84 is #define'd in oc/imports.h
// The OC_CONST84_CHAR macro can be used to wrap "const char*"
// arguments to the Tcl API that changed from char* to const char*
// between Tcl 8.3 and 8.4.  On Tcl 8.4 and later OC_CONST84_CHAR
// translates to a nop.  On Tcl 8.3 and earlier OC_CONST84_CHAR
// copies the const char* into a temporary Oc_AutoBuf, and forwards
// a pointer to the char buf in Oc_AutoBuf.
#if OC_USING_CONST84
# define OC_CONST84_CHAR(x) (x)
#else
# define OC_CONST84_CHAR(x) (Oc_AutoBuf(x).GetStr())
#endif



  
#if (OC_SYSTEM_TYPE == OC_WINDOWS)
  // Specialized conversion function for Windows.  Converts from char
  // to wide char.  This is mainly needed in the case UNICODE is
  // defined, because then TCHAR is a wide char (TCHAR is a normal
  // char otherwise.)
class Oc_AutoWideBuf
{
private:
  size_t bufsize; // Length of buffer pointed to by str. ( >= strlen(str)+1.)
  WCHAR *wstr;
public:
  ~Oc_AutoWideBuf() { if(bufsize>0) { delete[] wstr; } }
  Oc_AutoWideBuf()
    : bufsize(0), wstr(NULL) { Set(""); }
  Oc_AutoWideBuf(const char* const_str)
    : bufsize(0), wstr(NULL) { Set(const_str); }

  void Set(const char* const_str);

  operator WCHAR*() const { return wstr; }
  WCHAR* GetStr() const { return wstr; }
  WCHAR* operator()(const char* const_str);

  // The following are (to date) undefined
  Oc_AutoWideBuf(const Oc_AutoWideBuf& other);
  Oc_AutoWideBuf& operator=(const Oc_AutoWideBuf& other);

};

// Code to convert from char* to TCHAR*
// WARNING: If UNICODE is not defined (the usual case for us), then the
//    Oc_AutoTBuf typedef resolves to Oc_AutoBuf.  This is problematic
//    in that the other case (less common for us) the Oc_AutoTBuf
//    typedef resolves to Oc_AutoWideBuf, which has *less* total
//    functionality than the former.  In principle this should not be
//    a problem, because Oc_AutoWideBuf provides all the functionality
//    intended for Oc_AutoTBuf, but the risk is that an unwary
//    programmer may inadvertently use Oc_AutoBuf functionality on an
//    Oc_AutoTBuf variable, and this error won't get caught until
//    someone (perhaps in the field) tries to do a build with UNICODE
//    defined.

#ifdef UNICODE
typedef Oc_AutoWideBuf Oc_AutoTBuf;
#else
typedef Oc_AutoBuf Oc_AutoTBuf;
#endif

#endif // WINDOWS

#endif /* _OC_AUTOBUF */

