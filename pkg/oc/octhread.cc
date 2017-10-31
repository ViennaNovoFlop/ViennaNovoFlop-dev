/* FILE: octhread.cc                   -*-Mode: c++-*-
 *
 *      Support for C-level multi-threaded (parallel) code
 *
 *
 * NOTICE: Please see the file ../../LICENSE
 *
 * Last modified on: $Date: 2010-06-18 21:49:34 $
 * Last modified by: $Author: donahue $
 */

#include <errno.h>

#include "ocexcept.h"
#include "octhread.h"
#include "messages.h"

#if OC_USE_NUMA
// Note: The numa-devel package is needed to get numaif.h.
# include <numaif.h>
#endif // OC_USE_NUMA

/* End includes */     /* Optional directive to pimake */

/*
 * Number of threads to run
 */
#if OOMMF_THREADS
static int oommf_thread_count = 1;
int Oc_GetMaxThreadCount() { return oommf_thread_count; }

void Oc_SetMaxThreadCount(int threads) {
  if(threads<1) {
    OC_THROW(Oc_Exception(__FILE__,__LINE__,"",
                          "Oc_SetMaxThreadCount",
                          1024,
                          "Import parameter error:"
                          " \"threads\" value %d must be >=1.",threads));
  }
  oommf_thread_count = threads;
}
#endif // OOMMF_THREADS

////////////////////////////////////////////////////////////////////////
// Tcl wrappers for thread code
int OcHaveThreads(ClientData, Tcl_Interp *interp, int argc,CONST84 char **argv) 
{
  Tcl_ResetResult(interp);
  if (argc != 1) {
    Tcl_AppendResult(interp, argv[0], " takes no arguments", (char *) NULL);
    return TCL_ERROR;
  }
  char buf[256];
  Oc_Snprintf(buf,sizeof(buf),"%d",Oc_HaveThreads());
  Tcl_AppendResult(interp,buf,(char *)NULL);
  return TCL_OK;
}

int
OcGetMaxThreadCount(ClientData, Tcl_Interp *interp, int argc,CONST84 char **argv) 
{
  Tcl_ResetResult(interp);
  if (argc != 1) {
    Tcl_AppendResult(interp, argv[0], " takes no arguments", (char *) NULL);
    return TCL_ERROR;
  }
  char buf[256];
  Oc_Snprintf(buf,sizeof(buf),"%d",Oc_GetMaxThreadCount());
  Tcl_AppendResult(interp,buf,(char *)NULL);
  return TCL_OK;
}

int
OcSetMaxThreadCount(ClientData,Tcl_Interp *interp,int argc,CONST84 char** argv)
{
  char buf[1024];
  Tcl_ResetResult(interp);
  if(argc!=2) {
    Oc_Snprintf(buf,sizeof(buf),"Oc_SetMaxThreadCount must be called with 1 argument,"
	    " thread_count (%d arguments passed)",argc-1);
    Tcl_AppendResult(interp,buf,(char *)NULL);
    return TCL_ERROR;
  }
  int thread_count=atoi(argv[1]);

  if(thread_count<1) {
    Oc_Snprintf(buf,sizeof(buf),"Oc_SetMaxThreadCount input must be a positive integer"
	    " (got %d)", thread_count);
    Tcl_AppendResult(interp,buf,(char *)NULL);
    return TCL_ERROR;
  }

  try {
    Oc_SetMaxThreadCount(thread_count);
  } catch(char* errmsg) {
    Tcl_AppendResult(interp,errmsg,(char *)NULL);
    return TCL_ERROR;
  }

  return TCL_OK;
}


#if OC_USE_NUMA
/* The client selects which nodes to interleave memory
 * across by setting the import nodes;  If this list
 * is empty, then the default node#select arrays are
 * used instead.
 *   If import max_threads is >0, then the nodes
 * array is truncated at max_threads entries.
 *   These tables could no doubt be improved, especially
 * for the larger arrays.  It would also be nice if the
 * auto selection code could find out which nodes are
 * loaded (at least on problem initialization time) and
 * assign to less busy nodes.
 *   Addendum: These tables are based on core orderings
 * from AMD Opteron systems in the 2008-2010 time frame.
 * Orderings on Intel systems appear to be significantly
 * different: On AMD, cores on one processor or node are
 * numbered successively, which means for best memory
 * bandwidth one should skip node numbers.  On Intel,
 * successive numbers skip first across processors and
 * then across nodes, so that successive core numbers give
 * good memory through-put and good cache locality, but
 * the mapping is harder to anticipate so automated core
 * selection is harder.
 */
static int node32select[] = {  8,  9, 10, 11,
                              12, 13, 14, 15,
                              16, 17, 18, 19,
                              20, 21, 22, 23,
                              28, 29, 30, 31,
                              24, 25, 26, 27,
                               4,  5,  6,  7,
                               1,  2,  3,  0 } ;
static int node12select[] = { 2, 3, 4, 5, 10, 11, 6, 7, 8, 9, 1, 0} ;
static int node8select[] = { 2, 3, 4, 7, 5, 6, 1, 0 } ;
static int node4select[] = { 2, 3, 1, 0 } ;
static int node2select[] = { 1, 0 } ;
static int node1select[] = { 0 } ;
static vector<int> nodeselect;

static int numa_ready = 0;

void Oc_NumaInit(int max_threads,vector<int>& nodes) {
  // Init numa library
  numa_ready = 0;
  if(numa_available() == -1) {
    OC_THROW("ProgramLogicError: NUMA (non-uniform memory access)"
             " library not available.");
  }
#if 0
  numa_set_bind_policy(1);  /**/ // TEST
  numa_set_strict(1);
#endif
  // Set up nodes select list
  nodeselect.clear();
  if(nodes.size() == 0) {
    // Use default node list
    int* ns=NULL;
    int ns_size=0;
    switch(numa_max_node()+1) {
    case 1:
      ns=node1select;
      ns_size=1;
      break;
    case 2:
      ns=node2select;
      ns_size = 2;
      break;
    case 4:
      ns=node4select;
      ns_size = 4;
      break;
    case 8:
      ns=node8select;
      ns_size = 8;
      break;
    case 12:
      ns=node12select;
      ns_size = 12;
      break;
    case 32:
      ns=node32select;
      ns_size = 32;
      break;
    default:
      OC_THROW(Oc_Exception(__FILE__,__LINE__,"","Oc_NumaInit",1024,
                  "NUMA (non-uniform memory access) init error;"
                  " auto-config not supported for %d node machines.",
                  numa_max_node()+1));
    }
    if(max_threads>0 && max_threads<ns_size) ns_size=max_threads;
    for(int i=0;i<ns_size;++i) nodeselect.push_back(ns[i]);
  } else {
    const int maxnode = numa_max_node();
    int ns_size = nodes.size();
    if(max_threads>0 && max_threads<ns_size) ns_size=max_threads;
    for(int i=0;i<ns_size;++i) {
      if(nodes[i]>maxnode) {
        OC_THROW(Oc_Exception(__FILE__,__LINE__,"","Oc_NumaInit",1024,
           "NUMA (non-uniform memory access) init error;"
           " requested node (%d) is outside machine node range [0-%d].",
           nodes[i],maxnode));
      }
      nodeselect.push_back(nodes[i]);
    }
  }

  // Memory interleaving
#if defined(LIBNUMA_API_VERSION) && LIBNUMA_API_VERSION>=2
  struct bitmask* mask = numa_allocate_nodemask();
  numa_bitmask_clearall(mask);
  int ibitmax = static_cast<int>(nodeselect.size());
  for(int ibit=0;ibit<ibitmax;++ibit) {
    numa_bitmask_setbit(mask,static_cast<unsigned int>(nodeselect[ibit]));
  }
  numa_set_interleave_mask(mask);
  numa_free_nodemask(mask);
#else // LIBNUMA_API_VERSION
  nodemask_t mask;
  nodemask_zero(&mask);
  int ibitmax = static_cast<int>(nodeselect.size());
  for(int ibit=0;ibit<ibitmax;++ibit) {
    nodemask_set(&mask,nodeselect[ibit]);
  }
  numa_set_interleave_mask(&mask);
#endif // LIBNUMA_API_VERSION
  numa_ready = 1;
}

int Oc_NumaReady()
{
  return numa_ready;
}

void Oc_NumaRunOnNode(int thread)
{
  if(nodeselect.size()>0) {
    numa_run_on_node(nodeselect[thread % nodeselect.size()]);
  }
}

// Oc_SetMemoryPolicyFirstTouch sets the node allocation policy for
//  the given memory address range to "First Touch", meaning each page
//  is allocated to the node running the thread that first accesses
//  a location on that page.  This overrides the Oc_Numa default of
//  interleaved across the nodes specified in Oc_NumaInit.  Note that
//  Oc_SetMemoryPolicyFirstTouch must be called after memory allocation
//  (so an address is held) but before memory use.
void Oc_SetMemoryPolicyFirstTouch(char* start,OC_INDEX len)
{
  // What to do if the numa library is not initialized?  Options are
  // throw an exception or else return without doing anything.  The
  // advantage of the latter is that it may allow NUMA-compiled code
  // to run on systems w/o a working NUMA installation (the libraries
  // can be installed but won't be functional if the procfs isn't
  // setup for numa).  The advantage of the former is that it may
  // catch some coding errors.  For now, try the latter.
  if(!numa_ready) return;

  // For some reason, mbind reports an EINVAL (illegal parameter) error
  // if called as below if the system has only one node.  This routine
  // should be a nop in that case, so just handle this case specially.
  if(numa_max_node()==0) return;

  // The documentation is not as clear as it could be on this point, but
  // calling mbind with MPOL_PREFERRED and an empty node mask implements
  // a page-by-page first touch policy.
  unsigned long first_touch_node_mask = 0x0;
  int errcode = mbind(start,len,MPOL_PREFERRED,&first_touch_node_mask,
                      numa_max_node()+1,0);
  if(errcode != 0) {
    int errsav = errno;  // Save code from global errno
    OC_THROW(Oc_Exception(__FILE__,__LINE__,"",
                          "Oc_SetMemoryPolicyFirstTouch",
                          256,"%.100s",strerror(errsav)));
  }
}

int Oc_NumaGetRunNode(int thread)
{
  if(nodeselect.size()>0) {
    return nodeselect[thread % nodeselect.size()];
  }
  return -1;
}

int Oc_NumaGetLocalNode()
{ // Returns lowest-numbered node in run mask for current thread.
#if defined(LIBNUMA_API_VERSION) && LIBNUMA_API_VERSION>=2
  struct bitmask* mask = numa_get_run_node_mask();
  const int top_node = numa_max_node();
  int node = 0;
  while(node <= top_node) {
    if(numa_bitmask_isbitset(mask,(unsigned int)node)) break;
    ++node;
  }
  numa_free_nodemask(mask); // Docs don't say who is responsible for
  /// freeing mask, but some sample code in the numa compatibility
  /// library does this, so we will too.
  if(node>top_node) return -1;
  return node;
#else
  nodemask_t mask = numa_get_run_node_mask();
  int node = 0;
  const int top_node = numa_max_node();
  while(node <= top_node) {
    if(nodemask_isset(&mask,node)) return node;
    ++node;
  }
  return -1;
#endif
}

int Oc_NumaNodeFirstThread(int node)
{ /* Returns the smallest thread number running on node.  This can be
   * used in conjuction with Oc_NumaGetRunNode to build and access
   * node-based memory structures.  For example, if you place a copy of
   * "A" on each node, with an array of pointers, say type* A[n] where n
   * runs over available nodes, then A[0] points to the copy used by
   * thread 0 (which resides on node Oc_NumaGetRunNode(0)), and in
   * general thread k should use the copy stored at A[j] where j =
   * Oc_NumaNodeFirstThread(Oc_NumaGetRunNode(k)).
   *   Returns -1 if node is not used by any thread.
   */
  const int node_count = nodeselect.size();
  int i = 0;
  while(i<node_count && nodeselect[i] != node) ++i;
  if(i<node_count) return i;
  return -1;
}

#endif // OC_USE_NUMA

////////////////////////////////////////////////////////////////////////
// Memory allocation
#if OC_USE_NUMA
static int using_numa_alloc = -1;
void* Oc_AllocThreadLocal(size_t size)
{
  if(using_numa_alloc<0) {
    // Initialize
    using_numa_alloc = Oc_NumaAvailable();
  }
  if(!using_numa_alloc) {
    // Handle case where code was built for NUMA, but
    // NUMA not presently available.
    void* foo = malloc(size);
    if(!foo) { OC_THROW("Out of memory."); }
    return foo;
  }
#if 0
  void* foo = numa_alloc_local(size);  // Broken???
  if(!foo) { OC_THROW("Out of memory."); }
  return foo;
#else
  int node = Oc_NumaGetLocalNode();
  if(node<0) {
    OC_THROW("Unable to determine local run node.");
  }
  void* foo = numa_alloc_onnode(size,node);
  if(!foo) { OC_THROW("Out of memory."); }
  return foo;
#endif
}
void Oc_FreeThreadLocal(void* start,size_t size)
{
  if(!start) return; // Free of NULL pointer always allowed.
  if(using_numa_alloc<0) {
    OC_THROW("Oc_FreeThreadLocal() called before Oc_AllocThreadLocal()");
  }
  if(using_numa_alloc) {
    numa_free(start,size);
  } else {
    free(start);
  }
}
#else
void* Oc_AllocThreadLocal(size_t size)
{
  void* foo = malloc(size);
  if(!foo) { OC_THROW("Out of memory."); }
  return foo;
}
void Oc_FreeThreadLocal(void* start,size_t /* size */)
{
  if(!start) return; // Free of NULL pointer always allowed.
  free(start);
}
#endif

#if OC_USE_NUMA
void Oc_NumaNodemaskStringRep(nodemask_t imask,Oc_AutoBuf& ab)
{ // Converts a NUMA nodemask_t to a string representation.
  const int maxnode = numa_max_node();
  unsigned int buf = 0;
  unsigned int bit = 1;

  Oc_AutoBuf reverse_ab;
  for(int i=0;i<=maxnode;++i) {
    if(nodemask_isset(&imask,i)) {
      buf |= bit;
    }
    if( (bit <<= 1) == 0x10 ) {
      char cbuf[2];
      Oc_Snprintf(cbuf,sizeof(cbuf),"%X",buf);
      reverse_ab.Append(cbuf); // "Append" adds to the right, so
      buf = 0;                /// this stores the digits backwards.
      bit = 1;
    }
  }
  if(bit>1) { // Dump leftovers
    char cbuf[2];
    Oc_Snprintf(cbuf,sizeof(cbuf),"%X",buf);
    reverse_ab.Append(cbuf);
  }

  // Reverse and return
  const int strsize = static_cast<int>(strlen(reverse_ab.GetStr()));
  if(strsize==0) {
    ab = "0";
  } else {
    ab.SetLength(strsize);
    for(int j=0;j<strsize;++j) {
      ab[j] = reverse_ab[strsize-1-j];
    }
    ab[strsize] = '\0';
  }
}
void Oc_NumaGetInterleaveMask(Oc_AutoBuf& ab)
{ // Wrapper around nodemask_t numa_get_interleave_mask(void),
  // which fills export "ab" with a hexadecimal rendering of
  // the node mask.
  if(!numa_ready) {
    ab = "numa library not initialized";
  } else {
    nodemask_t imask = numa_get_interleave_mask();
    Oc_NumaNodemaskStringRep(imask,ab);
  }
}
#endif // OC_USE_NUMA

////////////////////////////////////////////////////////////////////////
// Tcl wrappers for NUMA code.

// Note 1: The numa_available() function from the NUMA library returns a
//   negative value on error, a 0 or positive value on success.  This
//   is rather non-intuitive; the Oc_NumaAvailable routine and the
//   associated Tcl wrapper return 1 on success, 0 on error (i.e., if
//   NUMA routines are not available).
// Note 2: Unlike Oc_NumaInit, the OcNumaInit wrapper also calls
//   Oc_NumaRunOnNode to tie the current thread to the first node.
int
OcNumaAvailable(ClientData, Tcl_Interp *interp, int argc,CONST84 char **argv) 
{
  Tcl_ResetResult(interp);
  if (argc != 1) {
    Tcl_AppendResult(interp, argv[0], " takes no arguments", (char *) NULL);
    return TCL_ERROR;
  }
  char buf[256];
  Oc_Snprintf(buf,sizeof(buf),"%d", Oc_NumaAvailable() );
  Tcl_AppendResult(interp,buf,(char *)NULL);
  return TCL_OK;
}

int
OcNumaInit(ClientData,Tcl_Interp *interp,int argc,CONST84 char** argv)
{
  char buf[1024];
  Tcl_ResetResult(interp);
  if(argc!=3) {
    Oc_Snprintf(buf,sizeof(buf),"Oc_NumaInit must be called with 2 arguments,"
	    " max_thread_count and node_list (%d arguments passed)",argc-1);
    Tcl_AppendResult(interp,buf,(char *)NULL);
    return TCL_ERROR;
  }
  int max_thread_count=atoi(argv[1]);

  if(max_thread_count<1) {
    Oc_Snprintf(buf,sizeof(buf),"Import max_thread_count to Oc_NumaInit"
                " must be a positive integer (got %d)", max_thread_count);
    Tcl_AppendResult(interp,buf,(char *)NULL);
    return TCL_ERROR;
  }

  vector<int> nodes;
  int node_argc;
  CONST84 char** node_argv;
  if(TCL_OK != Tcl_SplitList(interp, argv[2], &node_argc, &node_argv)) {
    Oc_Snprintf(buf,sizeof(buf),"Import node_list to Oc_NumaInit"
                " is not a proper list (got \"%.500s\")", argv[2]);
    Tcl_AppendResult(interp,buf,(char *)NULL);
    return TCL_ERROR;
  }
  for(int i=0;i<node_argc;++i) {
    char* endptr;
    long int number = strtol(node_argv[i],&endptr,10);
    if(node_argv[i] == endptr || *endptr != '\0') {
      // Invalid element
      Oc_Snprintf(buf,sizeof(buf),"Import node_list to Oc_NumaInit"
                  " is not a list of numbers (element %d is \"%.500s\")",
                  i,node_argv[i]);
      Tcl_AppendResult(interp,buf,(char *)NULL);
      Tcl_Free((char *)node_argv);
      return TCL_ERROR;
    }
    if(number<0) {
      // Invalid element
      Oc_Snprintf(buf,sizeof(buf),"Import node_list to Oc_NumaInit"
                  " is not a list of non-negative numbers (element %d is %d)",
                  i,number);
      Tcl_AppendResult(interp,buf,(char *)NULL);
      Tcl_Free((char *)node_argv);
      return TCL_ERROR;
    }
    nodes.push_back(number);
  }
  Tcl_Free((char *)node_argv);

  try {
    Oc_NumaInit(max_thread_count,nodes);
    Oc_NumaRunOnNode(0); // Tie main thread (i.e., this one) to first node
  } catch(char* errmsg) {
    Tcl_AppendResult(interp,errmsg,(char *)NULL);
    return TCL_ERROR;
  }

  return TCL_OK;
}

