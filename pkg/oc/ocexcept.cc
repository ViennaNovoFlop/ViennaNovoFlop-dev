/* FILE: except.cc                    -*-Mode: c++-*-
 *
 *   Simple exception object for C++ exceptions.
 * 
 * NOTICE: Please see the file ../../LICENSE
 *
 * Last modified on: $Date: 2008-09-09 23:52:31 $
 * Last modified by: $Author: donahue $
 */

#include <string.h>

#include "ocexcept.h"
#include "messages.h"

/* End includes */     /* Optional directive to pimake */

#define INT_DISPLAY_LENGTH_BOUND 256   // Should be safe.

// Base constructor
Oc_Exception::Oc_Exception
(const char* file_in,      // File from which exception is thrown
 int lineno_in,            // Line number from which exception is thrown
 const char* classname_in, // Name of class throwing exception
 const char* funcname_in,  // Name of function throwing exception
 int errmsg_size_in,       // Buffer size needed to hold extd. error msg
 const char* errfmt_in,   // Format string for extended error message
 ...)                      // arguments for errfmt
  : file(file_in), lineno(lineno_in),
    classname(classname_in), funcname(funcname_in)
{
  errmsg.SetLength(errmsg_size_in);
  
  va_list arg_ptr;
  va_start(arg_ptr,errfmt_in);
  Oc_Vsnprintf(errmsg.GetStr(),errmsg_size_in,errfmt_in,arg_ptr);
  va_end(arg_ptr);
}


void Oc_Exception::PrependMessage(const char* prefix)
{
  size_t prefix_len = strlen(prefix);
  Oc_AutoBuf tmp;
  tmp.SetLength(prefix_len+errmsg.GetLength());
  strcpy(tmp.GetStr(),prefix);
  strcpy(tmp.GetStr()+prefix_len,errmsg.GetStr());
  errmsg = tmp;
}

void Oc_Exception::PostpendMessage(const char* suffix)
{
  size_t suffix_len = strlen(suffix);
  errmsg.SetLength(errmsg.GetLength()+suffix_len);
  strcpy(errmsg.GetStr()+suffix_len,suffix);
}


// ConstructMessage is intended to build a generic error essage
// out of the data stored in this object.  We may add interfaces
// to the private data if there arises a need for specially
// tailored error messages.
const char* Oc_Exception::ConstructMessage(Oc_AutoBuf& msg) const {
  // File info
  Oc_AutoBuf fileinfo;
  const char* fileinfo_fmt
    = "Error detected; refer to program source file %s, line %d, ";
  int fileinfo_bufsize = strlen(fileinfo)
    + strlen(file.GetStr()) + INT_DISPLAY_LENGTH_BOUND;
  fileinfo.SetLength(fileinfo_bufsize);
  Oc_Snprintf(fileinfo.GetStr(),fileinfo_bufsize,
              fileinfo_fmt,file.GetStr(),lineno);

  // Class info
  Oc_AutoBuf proginfo;
  if(strlen(classname.GetStr())>0) {
    const char* proginfo_fmt="in class %s, member function %s.";
    int proginfo_bufsize = strlen(proginfo_fmt)
      + strlen(classname.GetStr()) + strlen(funcname.GetStr());
    proginfo.SetLength(proginfo_bufsize);
    Oc_Snprintf(proginfo.GetStr(),proginfo_bufsize,
                proginfo_fmt,classname.GetStr(),funcname.GetStr());
  } else {
    const char* proginfo_fmt="in function %s.";
    int proginfo_bufsize = strlen(proginfo_fmt)
      + strlen(funcname.GetStr());
    proginfo.SetLength(proginfo_bufsize);
    Oc_Snprintf(proginfo.GetStr(),proginfo_bufsize,
                proginfo_fmt,funcname.GetStr());
  }

  // Build message.
  const char* combine_fmt = "%s%s\nERROR DETAILS: %s";
  msg.SetLength(strlen(combine_fmt)
                +strlen(fileinfo.GetStr())
                +strlen(proginfo.GetStr())
                +strlen(errmsg.GetStr()));
  Oc_Snprintf(msg.GetStr(),msg.GetLength(),combine_fmt,
              fileinfo.GetStr(),proginfo.GetStr(),errmsg.GetStr());

  return msg.GetStr();
}
