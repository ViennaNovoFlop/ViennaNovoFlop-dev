/* FILE: tclcommand.cc                 -*-Mode: c++-*-
 *
 * Wrapper for Tcl_Eval().
 *
 */

#include <exception>

#include "tclcommand.h"
#include "util.h"

OC_USE_STD_NAMESPACE;

/* End includes */

void Oxs_TclCommand::Dump(String& contents) const
{
  char buf[4096];
  contents = "exception_prefix: \"";
  contents += exception_prefix;
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  contents += "\"\ncommand_base: \"";
  contents += command_base;
  contents += "\"\n";
  Oc_Snprintf(buf,sizeof(buf),"extra args count: %d\n",
	      extra_args_count);
  contents += buf;
  Oc_Snprintf(buf,sizeof(buf),"argc count: %d\n",command_args.size());
  contents += buf;
  for(int i=0;i<static_cast<int>(command_args.size());i++) {
    contents += " ->";
    contents += command_args[i];
    contents += "<-\n";
  }
  contents += "eval_result: \"";
  contents += eval_result;
  contents += "\"";
#else
  Oc_Snprintf(buf,sizeof(buf),"\"\nbase command size: %d\n",
	      base_command_size);
  contents += buf;
  Oc_Snprintf(buf,sizeof(buf),"extra args count: %d\n",
	      extra_args_count);
  contents += buf;
  Oc_Snprintf(buf,sizeof(buf),"objcmd size: %d\n",objcmd.Size());
  contents += buf;
  Tcl_Obj* const * arr = objcmd.Array();
  for(int i=0;i<objcmd.Size();i++) {
    Oc_Snprintf(buf,sizeof(buf)," Addr: %p ->%s<-\n",
		arr[i],objcmd.GetString(i).c_str());
    contents += buf;
    Oc_Snprintf(buf,sizeof(buf),
		"   refcount: %d, byteptr: %p, length: %d\n",
		arr[i]->refCount,arr[i]->bytes,arr[i]->length);
    contents += buf;
  }
  if(eval_result==NULL) {
      contents += "eval_result, Addr: <nilptr>";
  } else {
    Oc_Snprintf(buf,sizeof(buf),
		"eval_result, Addr: %p, ->%s<-\n",
		eval_result,Tcl_GetString(eval_result));
    contents += buf;
    Oc_Snprintf(buf,sizeof(buf),
		"   refcount: %d, byteptr: %p, length: %d\n",
		eval_result->refCount,eval_result->bytes,eval_result->length);
    contents += buf;
    if(eval_result->typePtr==NULL) {
      contents += "   Invalid type";
    } else {
      Oc_Snprintf(buf,sizeof(buf),
		  "   Type: %s",eval_result->typePtr->name);
      contents += buf;
    }
  }
#endif
}

Oxs_TclCommand::Oxs_TclCommand()
  : interp(NULL), no_too_few_restores_warning(0),
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
    eval_list_set(0),
#else
    base_command_size(0),objcmd(0),eval_result(NULL),
#endif // Tcl version check
    extra_args_count(0)
{}

Oxs_TclCommand::~Oxs_TclCommand()
{
  ReleaseBaseCommand();
}

void Oxs_TclCommand::ReleaseBaseCommand()
{
  extra_args_count=0;
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  eval_list_set=0;
  eval_list.Release();
  eval_result.erase();
  command_base.erase();
  command_args.clear();
#else
  if(eval_result!=NULL) {
    Tcl_DecrRefCount(eval_result);
    eval_result=NULL;
  }
  objcmd.Resize(0);
  base_command_size=0;
#endif // Tcl version check

  // Empty result stack.  We discard all stacked results, on the
  // theory that if the destructor is being called with items still
  // in the stack, then we are probably in an exceptional (read,
  // error) state, and perhaps the current Tcl result holds some
  // useful error information.  The alternative is to dump all
  // the results except the one at the bottom of the stack,
  // with a command like
  //
  // if(it == result_stack.rend()-1) Tcl_RestoreResult(interp,&(*it));
  //
  // The boolean "unbalanced" is set true if the result stack is
  // not empty.  This presumably indicates a programming error.
  // After freeing resources, we print a warning if unbalanced
  // is true, provided an exception is not currently being processed,
  // and also provided an earlier exception did not explicitly
  // request no warning by setting the member variable
  // no_too_few_restores_warning true.
  vector<Tcl_SavedResult>::reverse_iterator it = result_stack.rbegin();
  OC_BOOL unbalanced = (it != result_stack.rend());
  while(it != result_stack.rend()) {
    Tcl_DiscardResult(&(*it));
    ++it;
  }
  result_stack.clear();

  interp=NULL;

  if(unbalanced && !uncaught_exception()
     && !no_too_few_restores_warning) {
    PlainWarning("%s --- Unbalanced Tcl save/restore result calls"
		 " detected in Oxs_TclCommand: too few restores.",
		 exception_prefix.c_str());
  }

  no_too_few_restores_warning = 0; // Reset in any case

  exception_prefix.erase();
}

void Oxs_TclCommand::SetBaseCommand
(const char* exception_prefix_,
 Tcl_Interp* interp_,
 const String& cmdbase,
 int extra_args_count_)
{
  // Free old command, if any
  ReleaseBaseCommand();

  // Setup new command
  exception_prefix = exception_prefix_;
  interp = interp_;
  extra_args_count = extra_args_count_;
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  // Tcl string interface.
  command_base = cmdbase;
  for(int iarg=0;iarg<extra_args_count;iarg++) {
    command_args.push_back(String("")); // Dummy fill
  }
#else
  // Tcl obj interface
  Oxs_SplitList cmdlist;
  cmdlist.Split(cmdbase);
  base_command_size = cmdlist.Count();
  objcmd.Resize(base_command_size+extra_args_count);
  for(int iarg=0;iarg<base_command_size;iarg++) {
    objcmd.WriteString(iarg,cmdlist[iarg]);
  }
#endif // Tcl version check
}

void Oxs_TclCommand::SaveInterpResult() const
{ // Conceptually const
  result_stack.resize(result_stack.size()+1);
  Tcl_SaveResult(interp,&(result_stack.back()));
}

void Oxs_TclCommand::RestoreInterpResult() const
{ // Conceptually const
  if(result_stack.empty()) {
    String msg = exception_prefix
      + String(" --- Tcl result restore call"
	       " with no interp result saved.");
    OXS_THROW(Oxs_ProgramLogicError,msg);
  }
  Tcl_RestoreResult(interp,&(result_stack.back()));
  result_stack.pop_back();
}

void Oxs_TclCommand::DiscardInterpResult() const
{ // Conceptually const
  if(result_stack.empty()) {
    String msg = exception_prefix
      + String(" --- Tcl result discard call"
	       " with no interp result saved.");
    OXS_THROW(Oxs_ProgramLogicError,msg);
  }
  Tcl_DiscardResult(&(result_stack.back()));
  result_stack.pop_back();
}

void Oxs_TclCommand::SetCommandArg(int index,const char* arg) const
{ // Conceptually const
#ifndef NDEBUG
  if(index>=extra_args_count || index<0) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[4096];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::SetCommandArg: "
		"Array out-of-bounds; index=%d, arg=%.3000s",
		exception_prefix.c_str(),index,arg);
    OXS_THROW(Oxs_BadIndex,buf);
  }
#endif
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  command_args[index] = arg;
#else
  objcmd.WriteString(base_command_size+index,arg);
#endif // Tcl version check
}

void Oxs_TclCommand::SetCommandArg(int index,OC_INT4m arg) const
{ // Conceptually const
#ifndef NDEBUG
  if(index>=extra_args_count || index<0) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[1024];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::SetCommandArg: "
		"Array out-of-bounds; index=%d, arg=%ld",
		exception_prefix.c_str(),index,(long)arg);
    OXS_THROW(Oxs_BadIndex,buf);
  }
#endif
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  char buf[64];
  Oc_Snprintf(buf,sizeof(buf),"%ld",static_cast<long>(arg));
  command_args[index] = buf;
#else
  objcmd.WriteLong(base_command_size+index,static_cast<long>(arg));
#endif // Tcl version check
}

void Oxs_TclCommand::SetCommandArg(int index,OC_UINT4m arg) const
{ // Conceptually const
#ifndef NDEBUG
  if(index>=extra_args_count || index<0) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[1024];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::SetCommandArg: "
		"Array out-of-bounds; index=%d, arg=%ld",
		exception_prefix.c_str(),index,(long)arg);
    OXS_THROW(Oxs_BadIndex,buf);
  }
#endif
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  char buf[64];
  Oc_Snprintf(buf,sizeof(buf),"%lu",static_cast<unsigned long>(arg));
  command_args[index] = buf;
#else
  objcmd.WriteLong(base_command_size+index,static_cast<long>(arg));
  /// NOTE: there is no unsigned Tcl obj.  We could overflow
  /// check, but for now we'll just assume arg will fit into a long.
#endif // Tcl version check
}

void Oxs_TclCommand::SetCommandArg(int index,OC_REAL8m arg) const
{ // Conceptually const
#ifndef NDEBUG
  if(index>=extra_args_count || index<0) {
    char buf[1024];
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::SetCommandArg: "
		"Array out-of-bounds; index=%d, arg=%.17g",
		exception_prefix.c_str(),index,
                static_cast<double>(arg));
    OXS_THROW(Oxs_BadIndex,buf);
  }
#endif
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  char buf[64];
  Oc_Snprintf(buf,sizeof(buf),"%.17g",static_cast<double>(arg));
  command_args[index] = buf;
#else
  objcmd.WriteDouble(base_command_size+index,static_cast<double>(arg));
#endif // Tcl version check
}


void Oxs_TclCommand::Eval() const
{ // Conceptually const
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  // For Tcl eval string interface
  if(eval_list_set) {
    eval_list.Release();
    eval_list_set=0;
  }
  String script = command_base
    + String(" ") + Oxs_MergeList(command_args);

  int errcode = TCL_OK;
  try {
    errcode = Tcl_Eval(interp,OC_CONST84_CHAR(script.c_str()));
  } catch(...) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// This is expected to unbalance the results stack.  There does
    /// not seem to be any real reason to print an unbalanced result
    /// stack warning too.
    throw;
  }
  eval_result = Tcl_GetStringResult(interp);
#else // ! (Tcl version < 8.1)
  // Use Tcl eval obj interface
  if(eval_result!=NULL) {
    Tcl_DecrRefCount(eval_result); // Free old result
    eval_result=NULL; // Safety
  }

  int errcode = TCL_OK;
  try {
    errcode = Tcl_EvalObjv(interp,objcmd.Size(),objcmd.Array(),
			   TCL_EVAL_GLOBAL);
  } catch(...) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// This is expected to unbalance the results stack.  There does
    /// not seem to be any real reason to print an unbalanced result
    /// stack warning too.
    throw;
  }
  eval_result = Tcl_GetObjResult(interp);
  Tcl_IncrRefCount(eval_result);
#endif // (Tcl version < 8.1)

  if(errcode!=TCL_OK) {
    String msg = exception_prefix
      + String(" --- Error evaluating Tcl script: ");
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
    msg += script + String("\n--- Error message: ") + eval_result;
#else
    msg += Oxs_MergeList(objcmd)
      + String("\n--- Error message: ")
      + String(Tcl_GetString(eval_result));
#endif
    // Extended error info
    const char* ei = Tcl_GetVar(interp,OC_CONST84_CHAR("errorInfo"),
                                TCL_GLOBAL_ONLY);
    const char* ec = Tcl_GetVar(interp,OC_CONST84_CHAR("errorCode"),
                                TCL_GLOBAL_ONLY);
    if(ei==NULL) ei = "";   if(ec==NULL) ec = "";

    // NOTE TO SELF: Change this to 0 to test Don's improved
    // DisplayError proc in oxsi.tcl (mjd, 25-June-2002)
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// This is expected to unbalance the results stack, unless the
    /// client wrapped up the Eval call in a try block.  There does
    /// not seem to be any real reason to print an unbalanced result
    /// stack warning too.
    OXS_TCLTHROW(msg,String(ei),String(ec));
  }
}

String Oxs_TclCommand::GetWholeResult() const
{
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  return eval_result;
#else // ! (Tcl version < 8.1)
  if(eval_result==NULL) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    OXS_THROW(Oxs_ProgramLogicError,"No eval result;"
	      " Eval() has not been called on Oxs_TclCommand");
  }
  return String(Tcl_GetString(eval_result));
#endif // (Tcl version < 8.1)
}

#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
void Oxs_TclCommand::FillEvalList() const
{
  if(eval_list.Split(eval_result.c_str())!=TCL_OK) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    String msg = exception_prefix
      + String(" --- Tcl script return is not a list: ")
      + eval_result;
    OXS_THROW(Oxs_TclBadReturnType,msg);
  }
  eval_list_set=1;
}
#endif

int Oxs_TclCommand::GetResultListSize() const
{
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  if(!eval_list_set) FillEvalList();
  return eval_list.Count();
#else // ! (Tcl version < 8.1)
  if(eval_result==NULL) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    OXS_THROW(Oxs_ProgramLogicError,"No eval result;"
	      " Eval() has not been called on Oxs_TclCommand");
  }
  int size;
  if(Tcl_ListObjLength(NULL,eval_result,&size)!=TCL_OK) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    String msg = exception_prefix
      + String(" --- Tcl script return is not a list: ")
      + GetWholeResult();
    OXS_THROW(Oxs_TclBadReturnType,msg);
  }
  return size;
#endif // (Tcl version < 8.1)
}

void Oxs_TclCommand::GetResultListItem(int index,String& result) const
{
#ifndef NDEBUG
  int size = GetResultListSize();
  if(index>=size || index<0) {
    char buf[1024];
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::GetResultListItem:"
		" Array out-of-bounds; index=%d"
		" (should be >0 and <%d).",
		exception_prefix.c_str(),index,size);
    OXS_THROW(Oxs_BadIndex,buf);
  }
#endif
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  if(!eval_list_set) FillEvalList();
  result = eval_list[index];
#else // ! (Tcl version < 8.1)
  if(eval_result==NULL) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    OXS_THROW(Oxs_ProgramLogicError,"No eval result;"
	      " Eval() has not been called on Oxs_TclCommand");
  }
  Tcl_Obj* obj;
  if(Tcl_ListObjIndex(NULL,eval_result,index,&obj)!=TCL_OK) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    String msg = exception_prefix
      + String(" --- Tcl script return is not a list: ")
      + GetWholeResult();
    OXS_THROW(Oxs_TclBadReturnType,msg);
  }
  if(obj==NULL) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[1024];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::GetResultListItem:"
		" Array out-of-bounds; index=%d.",
		exception_prefix.c_str(),index);
    OXS_THROW(Oxs_BadIndex,buf);
  }
  result = Tcl_GetString(obj);
#endif // (Tcl version < 8.1)
}

void Oxs_TclCommand::GetResultListItem(int index,OC_INT4m& result) const
{
#ifndef NDEBUG
  int size = GetResultListSize();
  if(index>=size || index<0) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[1024];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::GetResultListItem:"
		" Array out-of-bounds; index=%d"
		" (should be >0 and <%d).",
		exception_prefix.c_str(),index,size);
    OXS_THROW(Oxs_BadIndex,buf);
  }
#endif
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  if(!eval_list_set) FillEvalList();
  char* cptr;
  result = static_cast<OC_INT4m>(strtol(eval_list[index],&cptr,10));
  if(eval_list[index][0]=='\0' || *cptr!='\0') {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[2048];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::GetResultListItem:"
		" Tcl script return item at index=%d"
		" is not an integer: %.500s",
		exception_prefix.c_str(),index,eval_list[index]);
    OXS_THROW(Oxs_BadIndex,buf);
  }
#else // ! (Tcl version < 8.1)
  if(eval_result==NULL) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    OXS_THROW(Oxs_ProgramLogicError,"No eval result;"
	      " Eval() has not been called on Oxs_TclCommand");
  }
  Tcl_Obj* obj;
  if(Tcl_ListObjIndex(NULL,eval_result,index,&obj)!=TCL_OK) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    String msg = exception_prefix
      + String(" --- Tcl script return is not a list: ")
      + GetWholeResult();
    OXS_THROW(Oxs_TclBadReturnType,msg);
  }
  if(obj==NULL) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[1024];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::GetResultListItem:"
		" Array out-of-bounds; index=%d.",
		exception_prefix.c_str(),index);
    OXS_THROW(Oxs_BadIndex,buf);
  }
  long lresult;
  if(Tcl_GetLongFromObj(NULL,obj,&lresult)!=TCL_OK) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[2048];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::GetResultListItem:"
		" Tcl script return item at index=%d"
		" is not a long integer: %.500s",
		exception_prefix.c_str(),index,Tcl_GetString(obj));
    OXS_THROW(Oxs_BadIndex,buf);
  }
  result = static_cast<OC_INT4m>(lresult);
#endif // (Tcl version < 8.1)
}

void Oxs_TclCommand::GetResultListItem(int index,OC_UINT4m& result) const
{
#ifndef NDEBUG
  int size = GetResultListSize();
  if(index>=size || index<0) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[1024];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::GetResultListItem:"
		" Array out-of-bounds; index=%d"
		" (should be >0 and <%d).",
		exception_prefix.c_str(),index,size);
    OXS_THROW(Oxs_BadIndex,buf);
  }
#endif
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  if(!eval_list_set) FillEvalList();
  char* cptr;
  result = static_cast<OC_UINT4m>(strtoul(eval_list[index],&cptr,10));
  if(eval_list[index][0]=='\0' || *cptr!='\0') {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[2048];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::GetResultListItem:"
		" Tcl script return item at index=%d"
		" is not an integer: %.500s",
		exception_prefix.c_str(),index,eval_list[index]);
    OXS_THROW(Oxs_BadIndex,buf);
  }
#else // ! (Tcl version < 8.1)
  if(eval_result==NULL) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    OXS_THROW(Oxs_ProgramLogicError,"No eval result;"
	      " Eval() has not been called on Oxs_TclCommand");
  }
  Tcl_Obj* obj;
  if(Tcl_ListObjIndex(NULL,eval_result,index,&obj)!=TCL_OK) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    String msg = exception_prefix
      + String(" --- Tcl script return is not a list: ")
      + GetWholeResult();
    OXS_THROW(Oxs_TclBadReturnType,msg);
  }
  if(obj==NULL) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[1024];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::GetResultListItem:"
		" Array out-of-bounds; index=%d.",
		exception_prefix.c_str(),index);
    OXS_THROW(Oxs_BadIndex,buf);
  }
  long lresult;
  if(Tcl_GetLongFromObj(NULL,obj,&lresult)!=TCL_OK) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[2048];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::GetResultListItem:"
		" Tcl script return item at index=%d"
		" is not a long integer: %.500s",
		exception_prefix.c_str(),index,Tcl_GetString(obj));
    OXS_THROW(Oxs_BadIndex,buf);
  }
  result = static_cast<OC_UINT4m>(lresult);
  // Note: There is no unsigned Tcl obj.
#endif // (Tcl version < 8.1)
}

void Oxs_TclCommand::GetResultListItem(int index,OC_REAL8m& result) const
{
#ifndef NDEBUG
  int size = GetResultListSize();
  if(index>=size || index<0) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[1024];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::GetResultListItem:"
		" Array out-of-bounds; index=%d"
		" (should be >0 and <%d).",
		exception_prefix.c_str(),index,size);
    OXS_THROW(Oxs_BadIndex,buf);
  }
#endif
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  if(!eval_list_set) FillEvalList();
  OC_BOOL err;
  result = static_cast<OC_REAL8m>(Nb_Atof(eval_list[index],err));
  if(err) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[2048];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::GetResultListItem:"
		" Tcl script return item at index=%d"
		" is not a floating point number: %.500s",
		exception_prefix.c_str(),index,eval_list[index]);
    OXS_THROW(Oxs_BadIndex,buf);
  }
#else // Obj interface
  if(eval_result==NULL) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    OXS_THROW(Oxs_ProgramLogicError,"No eval result;"
	      " Eval() has not been called on Oxs_TclCommand");
  }
  Tcl_Obj* obj;
  if(Tcl_ListObjIndex(NULL,eval_result,index,&obj)!=TCL_OK) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    String msg = exception_prefix
      + String(" --- Tcl script return is not a list: ")
      + GetWholeResult();
    OXS_THROW(Oxs_TclBadReturnType,msg);
  }
  if(obj==NULL) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[1024];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::GetResultListItem:"
		" Array out-of-bounds; index=%d.",
		exception_prefix.c_str(),index);
    OXS_THROW(Oxs_BadIndex,buf);
  }
  double tmpresult;
  if(Tcl_GetDoubleFromObj(NULL,obj,&tmpresult)!=TCL_OK) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    char buf[2048];
    Oc_Snprintf(buf,sizeof(buf),
		"%.800s --- void Oxs_TclCommand::GetResultListItem:"
		" Tcl script return item at index=%d"
		" is not a floating point number: %.500s",
		exception_prefix.c_str(),index,Tcl_GetString(obj));
    OXS_THROW(Oxs_BadIndex,buf);
  }
  result = static_cast<OC_REAL8m>(tmpresult);
#endif // (Tcl version < 8.1)
}

#if !OC_REAL8m_IS_REAL8
void Oxs_TclCommand::GetResultListItem(int index,OC_REAL8& result) const
{
  OC_REAL8m result8m;
  GetResultListItem(index,result8m);
  result = result8m;
}
#elif !OC_REAL8m_IS_DOUBLE
void Oxs_TclCommand::GetResultListItem(int index,double& result) const
{
  OC_REAL8m result8m;
  GetResultListItem(index,result8m);
  result = result8m;
}
#endif

void Oxs_TclCommand::GetResultList(vector<String>& result_list) const
{
#if OXS_TCL_COMMAND_USE_STRING_INTERFACE
  if(!eval_list_set) FillEvalList();
  eval_list.FillParams(result_list);
#else // ! (Tcl version < 8.1)
  if(eval_result==NULL) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    OXS_THROW(Oxs_ProgramLogicError,"No eval result;"
	      " Eval() has not been called on Oxs_TclCommand");
  }
  int obj_count;
  Tcl_Obj** obj_arr;
  if(Tcl_ListObjGetElements(NULL,eval_result,&obj_count,&obj_arr)!=TCL_OK) {
    no_too_few_restores_warning = 1; // We're throwing an exception.
    /// See notes in ::Eval() about result stack and exceptions
    String msg = exception_prefix
      + String(" --- Tcl script return is not a list: ")
      + GetWholeResult();
    OXS_THROW(Oxs_TclBadReturnType,msg);
  }
  result_list.clear();
  for(int i=0;i<obj_count;i++) {
    result_list.push_back(String(Tcl_GetString(obj_arr[i])));
  }
#endif // (Tcl version < 8.1)
}

/////////////////////////////////////////////////////////////////////
//
// Oxs_ParseTclCommandLineRequest, intended to be used in support of
// Oxs_TclCommand. The label and param_count fields of each
// Oxs_TclCommandLineOption object in the import/export "options"
// parameter should be filled on entry. Import request_string is
// a Tcl list, where each list item matches the label portion of
// some Oxs_TclCommandLineOption object in "options".
// Oxs_ParseTclCommandLineRequest fills in the position field of
// each Oxs_TclCommandLineOption object corresponding to the
// position of that option in the request_string.  If the option
// is not used, the position field is set to -1.  The position
// is padded to account for the number of parameters in each
// preceding option in request_string.  The return value is the
// number of total parameters associated with the options selected
// in request_string.  Exceptions are thrown if request_string
// contains an option not included in "options", or if the same
// option is selected more than once.
//   NOTE: If the param_count value for an Oxs_TclCommandLineOption
// is <1, then that option is not available in the current instance.
// This is a convenience for clients of Oxs_TclCommand, because it
// allows an option to hold a place in the Oxs_TclCommandLineOption
// "options" vector returned by Oxs_ParseTclCommandLineRequest, w/o
// allowing it to be used.  (Intended for options which are
// conditionally available, where <1 indicates that the necessary
// prerequisites are missing.))
//
int Oxs_ParseTclCommandLineRequest
(const char* exception_prefix,
 vector<Oxs_TclCommandLineOption>& options,
 const String& request_string)
{
  // Initialize options vector
  vector<Oxs_TclCommandLineOption>::iterator opit = options.begin();
  while(opit != options.end()) {
    opit->position = -1;
    ++opit;
  }

  // Split request string
  Oxs_SplitList request;
  request.Split(request_string);

  // Process request string
  int total_param_count=0;
  for(int i=0;i<request.Count();i++) {
    opit = options.begin();
    while(opit != options.end()) {
      if(opit->Label().compare(request[i])==0) {
	// Match
	if(opit->position >= 0) {
	  // Option already selected
	  String msg = String(exception_prefix)
	    + String(" --- Duplicate option request: \"")
	    + request[i] + String("\"");
	  OXS_THROW(Oxs_BadParameter,msg);
	}
	if(opit->ParamCount()<1) {
	  // Option not supported in this instance
	  String msg = String(exception_prefix)
	    + String(" --- Invalid option request: \"")
	    + request[i]
	    + String("\"; prerequisite(s) missing.");
	  OXS_THROW(Oxs_BadParameter,msg);
	}
	opit->position = total_param_count;
	total_param_count += opit->ParamCount();
	break;
      }
      ++opit;
    }
    if(opit == options.end()) {
      // No match
      String msg = String(exception_prefix)
	+ String("Unrecognized option request: \"")
	+ request[i] + String("\"");
      OXS_THROW(Oxs_BadParameter,msg);
    }
  }

  return total_param_count;
}

