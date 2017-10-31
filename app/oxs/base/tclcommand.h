/* FILE: tclcommand.h                 -*-Mode: c++-*-
 *
 * Wrapper for Tcl_Eval().
 *
 */

#ifndef _OXS_TCLCOMMAND
#define _OXS_TCLCOMMAND

#include <string>
#include <vector>

#include "oc.h"
#include "util.h"

OC_USE_STRING;

/* End includes */

// Use the TclObj interface with Tcl later than 8.0; earlier than that
// we must use the string interface.  However, one can set the macro
// OXS_TCL_COMMAND_USE_STRING_INTERFACE to 1 to force the use of the
// string interface.  This may be useful for debugging purposes.
#ifndef OXS_TCL_COMMAND_USE_STRING_INTERFACE
# if TCL_MAJOR_VERSION < 8 || (TCL_MAJOR_VERSION==8 && TCL_MINOR_VERSION==0)
#  define OXS_TCL_COMMAND_USE_STRING_INTERFACE 1
# else
#  define OXS_TCL_COMMAND_USE_STRING_INTERFACE 0
# endif
#endif

// Oxs_TclCommand is a wrapper for Tcl_Eval.  The script to be evaluated
// is divided into two parts: the base command and the args.  Setting
// the args and evaluating the command are viewed as a function call on
// the base command, and as such are conceptually const.  For this
// reason, the args and result are stored in mutable storage, so that
// they may be changed by a const Oxs_TclCommand object.
//   For similar reasons, saving and restoring the interp result is
// also considered conceptually const.
class Oxs_TclCommand
{
private:
  String exception_prefix;
  Tcl_Interp* interp;
  mutable vector<Tcl_SavedResult> result_stack;
  mutable OC_BOOL no_too_few_restores_warning; // If this is set,
  /// then no warning is raised on cleanup if there are unpopped
  /// results.  This is intended for exception handling.  It
  /// resets to 0 after each cleanup.
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  String command_base;  // Non-mutable!
  mutable vector<String> command_args;
  mutable String eval_result;
  mutable OC_BOOL eval_list_set;
  mutable Oxs_SplitList eval_list; // Cached list version of eval_result
  void FillEvalList() const;
#else
  int base_command_size;
  mutable Nb_TclObjArray objcmd;
  mutable Tcl_Obj* eval_result;
#endif // Tcl version check
  int extra_args_count;
public:
  Oxs_TclCommand();
  ~Oxs_TclCommand();

  void Dump(String& contents) const; // For debugging

  void SetBaseCommand(const char* exception_prefix_,
		      Tcl_Interp* interp_,
		      const String& cmdbase,
		      int extra_args_count_);
  void ReleaseBaseCommand();

  void SaveInterpResult() const;
  void RestoreInterpResult() const;
  void DiscardInterpResult() const;

  void SetCommandArg(int index,const char* arg) const;
  void SetCommandArg(int index,const String& arg) const
    { SetCommandArg(index,arg.c_str()); }
  void SetCommandArg(int index,OC_INT4m arg) const;
  void SetCommandArg(int index,OC_UINT4m arg) const;
  void SetCommandArg(int index,OC_REAL8m arg) const;

  void Eval() const;

  String GetWholeResult() const;

  int GetResultListSize() const;
  void GetResultListItem(int index,String& result) const;
  void GetResultListItem(int index,OC_INT4m& result) const;
  void GetResultListItem(int index,OC_UINT4m& result) const;
  void GetResultListItem(int index,OC_REAL8m& result) const;
#if !OC_REAL8m_IS_REAL8
  void GetResultListItem(int index,OC_REAL8& result) const;
#elif !OC_REAL8m_IS_DOUBLE
  void GetResultListItem(int index,double& result) const;
#endif
  void GetResultList(vector<String>& result_list) const;
};


/////////////////////////////////////////////////////////////////////
//
// Oxs_ParseTclCommandLineRequest, intended to be used in support of
// Oxs_TclCommand. The label and param_count fields of each
// Oxs_TclCommandLineOptions object in the import/export "options"
// parameter should be filled on entry. Import request_string is
// a Tcl list, where each list item matches the label portion of
// some Oxs_TclCommandLineOptions object in "options".
// Oxs_ParseTclCommandLineRequest fills in the position field of
// each Oxs_TclCommandLineOptions object corresponding to the
// position of that option in the request_string.  If the option
// is not used, the position field is set to -1.  The position
// is padded to account for the number of parameters in each
// preceding option in request_string.  The return value is the
// number of total parameters associated with the options selected
// in request_string.  Exceptions are thrown if request_string
// contains an option not included in "options", or if the same
// options is selected more than once.
//   NOTE: If the param_count value for an Oxs_TclCommandLineOption
// is <1, then that option is not available in the current instance.
// This is a convenience for clients of Oxs_TclCommand, because it
// allows an option to hold a place in the Oxs_TclCommandLineOption
// "options" vector returned by Oxs_ParseTclCommandLineRequest, w/o
// allowing it to be used.  (Intended for options which are
// conditionally available, where <1 indicates that the necessary
// prerequisites are missing.)
//
class Oxs_TclCommandLineOption
{
private:
  String label; // Identifying string for option
  int param_count; // Number of parameters associated with option
  /// I'd like to make label and param_count const, but that would
  /// disallow operator=, which is needed if the class is used
  /// inside STL containers such as vector.
public:
  const String& Label() { return label; }
  int ParamCount() { return param_count; }
  int position; // Command line position; -1 indicates option
  /// not used.
  Oxs_TclCommandLineOption() {};
  Oxs_TclCommandLineOption(String label_,int param_count_)
    : label(label_), param_count(param_count_), position(-1) {}
};
OC_BOOL operator<(const Oxs_TclCommandLineOption&, const Oxs_TclCommandLineOption&);
OC_BOOL operator>(const Oxs_TclCommandLineOption&, const Oxs_TclCommandLineOption&);
OC_BOOL operator==(const Oxs_TclCommandLineOption&, const Oxs_TclCommandLineOption&);
OC_BOOL operator<(const Tcl_SavedResult&, const Tcl_SavedResult&);
OC_BOOL operator>(const Tcl_SavedResult&, const Tcl_SavedResult&);
OC_BOOL operator==(const Tcl_SavedResult&, const Tcl_SavedResult&);

int Oxs_ParseTclCommandLineRequest
(const char* exception_prefix,
 vector<Oxs_TclCommandLineOption>& options,
 const String& request_string);

#endif // _OXS_TCLCOMMAND
