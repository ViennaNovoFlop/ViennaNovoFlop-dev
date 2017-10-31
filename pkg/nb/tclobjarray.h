/* FILE: tclobjarray.h          -*-Mode: c++-*-
 *
 * Wrapper for Tcl_Obj objects.  Inert if building against Tcl 8.0
 * or earlier.
 *
 * NOTICE: Please see the file ../../LICENSE
 *
 * Last modified on: $Date: 2007-03-19 19:12:26 $
 * Last modified by: $Author: donahue $
 */

#ifndef _NB_TCLOBJARRAY
#define _NB_TCLOBJARRAY

#include <string>
#include "oc.h"

OC_USE_STRING;         // Map String --> std::string

/* End includes */

#if (TCL_MAJOR_VERSION > 8 || (TCL_MAJOR_VERSION==8 && TCL_MINOR_VERSION>0))
////////////////////////////////////////////////////////////////////////
// Wrapper for Tcl_Obj**
typedef Tcl_Obj* NbTclObjPtr;
class Nb_TclObjArray
{
private:
  int size;
  Tcl_Obj** arr;
  void Alloc(int arrsize);
  void Free();
public:
  Nb_TclObjArray() : size(0), arr(NULL) {}
  Nb_TclObjArray(int arrsize) : size(0), arr(NULL) { Alloc(arrsize); }
  ~Nb_TclObjArray() { Free(); }

  void Resize(int newsize) { Free(); Alloc(newsize); }

  NbTclObjPtr& operator[](int index);
  const NbTclObjPtr& operator[](int index) const;

  int Size() const { return size; }
  Tcl_Obj** Array() { return arr; }
  Tcl_Obj* const * Array() const { return arr; }

  // Write interface routines.  These handle the 
  // copy-on-write semantics of the Tcl referencing
  // system.
  void WriteString(int index,const char*);
  void WriteString(int index,const String& str);
  void WriteDouble(int index,double val);
  void WriteInt(int index,int val);
  void WriteLong(int index,long val);
  /// Add more as needed

  // Read interface routines.
  String GetString(int index) const;

};
#endif // Tcl version check


#endif // _NB_TCLOBJARRAY
