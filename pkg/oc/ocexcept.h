/* FILE: except.h                    -*-Mode: c++-*-
 *
 *   Simple exception object for C++ exceptions.
 * 
 * NOTICE: Please see the file ../../LICENSE
 *
 * Last modified on: $Date: 2008-09-09 23:52:31 $
 * Last modified by: $Author: donahue $
 */

#ifndef _OC_EXCEPT
#define _OC_EXCEPT

#include "autobuf.h"

/* End includes */     /* Optional directive to pimake */

////////////////////////////////////////////////////////////////////////
// Simple exception object.  Sample usage is:
//
//   OC_THROW(Oc_Exception(__FILE__,__LINE__,"My_Class","My_Function",
//            4096,"Too many fubars (%d) for input string %.3500s",
//            fubar_count,fubar_input);
//
// The fifth argument to the Oc_Exception constructor, errmsg_size (here
// 4096), must be large enough to hold the result of evaluating the
// format string errfmt with the appended arguments.  To is strongly
// recommended that any %s format directives include a precision
// setting (".3500" in the above example) to protect against overflow
// in filling the extended error message buffer.
//
class Oc_Exception {
 public:
  Oc_Exception
  (const char* file,      // File from which exception is thrown
   int lineno,            // Line number from which exception is thrown
   const char* classname, // Name of class throwing exception
   const char* funcname,  // Name of function throwing exception
   int errmsg_size,       // Buffer size needed to hold extd. error msg
   const char* errfmt,    // Format string for extended error message
   ...);                  // arguments for errfmt
  // Any of the char* may be set to NULL, in which case an empty string
  // will be substituted.

  // ConstructMessage is intended to build a generic error message
  // out of the data stored in this object.  We may add interfaces
  // to the private data if there arises a need for specially
  // tailored error messages.  As a convenience, the return value is
  // a pointer to the import msg buffer.
  const char* ConstructMessage(Oc_AutoBuf& msg) const;

  // Routines to modify errmsg.  These are intended for use
  // by routines that catch and rethrow the exception, in
  // order to add additional info.
  void PrependMessage(const char* prefix);
  void PostpendMessage(const char* suffix);

 private:
  Oc_AutoBuf file;      // File from which exception is thrown
  int        lineno;    // Line number from which exception is thrown
  Oc_AutoBuf classname; // Name of class throwing exception, or
  /// empty if exception is being thrown from outside a class.
  Oc_AutoBuf funcname;  // Name of function throwing exception
  Oc_AutoBuf errmsg;    // Extended error information
};

#endif // _OC_EXCEPT
