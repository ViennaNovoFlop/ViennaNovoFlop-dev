/* FILE: oxsthread.h              -*-Mode: c++-*-
 *
 * OXS thread support.
 *
 * Threaded/Non-Threaded Note: oxsthread.h and oxsthread.cc are designed
 * so that a threaded program (#if OOMMF_THREADS) can compile and run
 * properly in a non-threaded (#if !OOMMF_THREADS) context, provided that
 * one uses only the following functions and classes:
 *
 *           Oxs_ThreadPrintf
 *           Oxs_ThreadRunObj
 *           Oxs_Mutex
 *           Oxs_ThreadControl
 *           Oxs_ThreadTree
 *           Oxs_ThreadMapDataObject
 *           Oxs_ThreadLocalMap
 *
 *           Oxs_ThreadBush
 *           Oxs_ThreadTwig
 */

#ifndef _OXS_THREAD
#define _OXS_THREAD

#include <assert.h>
#include <map>    // Used for thread local storage
#include <vector>

#include <string.h>

#include "oc.h"  /* includes tcl.h */

/* End includes */

// TODO: Check mutex and thread handling in the case where an exception
// is thrown inside one thread while other threads are running.

// Thread-safe printf
void Oxs_ThreadPrintf(FILE* fptr,const char *format, ...);

class Oxs_ThreadControl; // Forward reference

// C++ wrapper around Tcl_Mutex
#if OOMMF_THREADS
class Oxs_Mutex {
  friend class Oxs_ThreadControl;
private:
  Tcl_Mutex mutex;
  volatile int lock_held;
public:
  Oxs_Mutex() : mutex(0), lock_held(0) {}
  ~Oxs_Mutex() { Tcl_MutexFinalize(&mutex); }

  int IsLockProbablyFree() {
    // This routine returns true is the mutex is *probably* free, by
    // looking at the value of lock_held.  This is intended as a
    // probabilistic aid to cooperative locking.  The idea is that a
    // thread can call here to check if the mutex is free, and if not
    // then do some other useful work and try again later.  However, it
    // can happen that IsLockProbablyFree returns true, but that an
    // immediately following Lock call will block because some other
    // thread grabbed the mutex between the IsLockProbablyFree call and
    // the Lock call (or already had it and just hadn't set lock_held
    // yet, or had reset lock_held but not yet released the mutex).  So,
    // no guarantees, but most of the time one hopes that a Lock call
    // will following a true return from IsLockProbablyFree will return
    // without blocking.
    //    BTW, in some circumstances on some shared memory architectures,
    // a read check (of say lock_held) may be possible without causing
    // any processor-to-processor communication.
    //    Note: Presumably lock_held is written and read atomically (since
    // it is an "int"), but I don't know if that really matters.
    return !lock_held;
  }

  void Lock()   {
    Tcl_MutexLock(&mutex);
    const_cast<int&>(lock_held)=1;  // lock_held not volatile inside mutex
  }

  void Unlock() {
    const_cast<int&>(lock_held)=0;
    Tcl_MutexUnlock(&mutex);
  }

  void Reset() {
    if(lock_held) {
      Oxs_ThreadPrintf(stderr,
                "\n*** Error? Oxs_Mutex finalized while locked. ***\n");
    }
    Tcl_MutexFinalize(&mutex);
    lock_held = 0;
    mutex = 0;
  }

};

class Oxs_ThreadError {
private:
  static Oxs_Mutex mutex;
  static int error;
  static String msg;
public:
  static void SetError(String errmsg) {
    mutex.Lock();
    if(error) msg += String("\n---\n");
    error=1;
    msg += errmsg;
    mutex.Unlock();
  }
  static void ClearError() {
    mutex.Lock();
    error=0;
    msg.clear();
    mutex.Unlock();
  }
  static int IsError(String* errmsg=0)  {
    int check;
    mutex.Lock();
    check=error;
    if(check && errmsg!=0) *errmsg = msg;
    mutex.Unlock();
    return check;
  }
};

#else  // OOMMF_THREADS
class Oxs_Mutex {
public:
  Oxs_Mutex() {}
  ~Oxs_Mutex() {}
  void Lock() {}
  void Unlock() {}
  int IsLockProbablyFree() { return 1; }
};

class Oxs_ThreadError {
public:
  static void SetError(String) {}
  static void ClearError() {}
  static int IsError(String* =0)  { return 0; }
};

#endif // OOMMF_THREADS


// Thread control trio: mutex, condition, counter
#if OOMMF_THREADS
class Oxs_ThreadControl {
private:
  Oxs_Mutex mutex;    // In general, this mutex should be acquired
  Tcl_Condition cond; // before changing either cond or count.
public:
  int count;    // In typical usage, count is used to protect against
  /// spurious notifies; if count>0 then a thread is waiting on cond.
  /// Code should Lock mutex before changing count.

  Oxs_ThreadControl() : cond(0), count(0) {}
  ~Oxs_ThreadControl() {
    Tcl_ConditionFinalize(&cond);    cond=0; 
    mutex.Lock();
    count=0;  // Safety
    mutex.Unlock();
    // Oxs_Mutex "mutex" is automatically finalized on destruction.
  }

  void Lock() { mutex.Lock(); }
  void Unlock() { mutex.Unlock(); }

  void GetStatus(int& locked, int& count_value) {
    // For debugging only.  Results are not reliable, in "thread" sense.
    locked = mutex.IsLockProbablyFree();
    count_value = count;
  }

  // Code should acquire mutex (through Lock()) before calling
  // Wait() or Notify().  Mutex is automatically freed while
  // inside Tcl_ConditionWait(), but automatically re-acquired
  // upon exit of Tcl_ConditionWait().
  void Wait(OC_UINT4m wait_interval=0) {
    // Import wait_interval is time in microseconds.  Set to
    // 0 for no timeout.
    if(wait_interval==0) {
      mutex.lock_held = 0;
      Tcl_ConditionWait(&cond,&mutex.mutex,NULL);
      mutex.lock_held = 1;
    } else {
      Tcl_Time ttinterval;
      ttinterval.sec = (long)(wait_interval/1000000);
      ttinterval.usec = (long)(wait_interval - 1000000*ttinterval.sec);
      mutex.lock_held = 0;
      Tcl_ConditionWait(&cond,&mutex.mutex,&ttinterval);
      mutex.lock_held = 1;
    }
  }
  void Notify() {
    Tcl_ConditionNotify(&cond);
  }
};
#else
// Thread control trio: mutex, condition, counter
class Oxs_ThreadControl {
public:
  int count;
  Oxs_ThreadControl() : count(0) {}
  ~Oxs_ThreadControl() { count=0; }
  void Lock() {}
  void Unlock() {}
  void GetStatus(int& locked,int& count_value) {
    locked = 0; count_value = count;
  }
  void Wait(OC_UINT4m /* wait_interval */=0) {}
  void Notify() {}
};
#endif

#if OOMMF_THREADS
// IMPLEMENTATION FOR THREADED BUILDS //////////////////////////////////

// NOTE: Oxs_Thread is for oxsthread internal use only!  See
// Threaded/Non-Threaded Note above

class Oxs_ThreadRunObj; // forward declaration

// Object for representation/control of actual thread
class Oxs_Thread {
 public:
  Oxs_Thread(int thread_number);
  /// Note: (Child) thread numbers start at 1, not 0.  Thread_number 0 is
  /// reserved for the main (parent) thread.

  ~Oxs_Thread();
  void RunCmd(Oxs_ThreadControl& stop,Oxs_ThreadRunObj& runobj,void* data);

  // When start.count>0, control OS thread is paused, waiting for a
  // new command to run.  start.count==0 means thread is active.
  int IsRunning() { return (0==start.count); }

  void GetStatus(String& results); // For debugging.

 private:
  friend Tcl_ThreadCreateProc _Oxs_Thread_threadmain;
  friend class Oxs_ThreadRunObj;
  Tcl_ThreadId thread_id;   // Assigned by Tcl and/or OS
  int thread_number;        // Assigned by Oxs_ThreadTree static control.
  /// These numbers are 1-based, because thread_number 0 is reserved
  /// for the main (parent) thread.

  Oxs_ThreadControl start;  // Unique to thread object
  Oxs_ThreadControl* stop;  // Assigned by and unique to tree object
  Oxs_ThreadRunObj* runobj; // Assigned by user code through Oxs_ThreadTree
  void* data;               // Ditto

  // Disable copy constructor and operator=() by declaring
  // without providing implementations.
  Oxs_Thread(const Oxs_Thread&);
  Oxs_Thread& operator=(const Oxs_Thread&);

};

// Manager for Oxs_Thread objects.  Call Oxs_ThreadTree to
// start and stop threads.
class Oxs_ThreadTree {
public:
  Oxs_ThreadTree() : threads_unjoined(0) {}
  ~Oxs_ThreadTree();

  static void EndThreads();

  void Launch(Oxs_ThreadRunObj& runobj,void* data); // Launch
  /// threads, one at a time.
  void Join(); // Joins all running threads.

  void LaunchRoot(Oxs_ThreadRunObj& runobj,void* data); // Runs
  /// command on root (i.e., 0) thread, and calls Join() to
  /// wait for all children to finish, with error handling.

  // Run on a range of threads.  Creates new threads as needed.
  // Set first_thread to 0 to include main (host) thread.  Set
  // last_thread to -1 to include uppermost existing thread.
  // This routine is self-joining.  Do not use this routine
  // in parallel with Launch.  RunOnThreadRange will block
  // until there are no running threads.
  void RunOnThreadRange(int first_thread,int last_thread,
                        Oxs_ThreadRunObj& runobj,void* data);

  void DeleteLockerItem(String name); // Deletes local locker item
  /// named "name" in all existing threads (including thread 0).
  /// No error is raised if name does not exist in some or all
  /// of the thread lockers.

  static void GetStatus(String& results); // For debugging.  Returns
  /// string detailing "probable" status for all mutexes and thread
  /// controls that Oxs_ThreadTree knows about, in human readable
  /// format.

private:
  Oxs_ThreadControl stop; // One per tree
  int threads_unjoined;   // Count of launched, unjoined threads.
  /// Mirrors stop.count, except modified only by *this, not
  /// client threads.

  static Oxs_Mutex launch_mutex;
  static std::vector<Oxs_Thread*> threads;
  /// Note: Child thread_number's are one larger than their index into
  /// the "threads" vector, because thread_number 0 is reserved for the
  /// main (parent) thread, which is not referenced in "threads".

  // Disable copy constructor and operator=() by declaring
  // without providing implementations.
  Oxs_ThreadTree(const Oxs_ThreadTree&);
  Oxs_ThreadTree& operator=(const Oxs_ThreadTree&);
};

#endif // OOMMF_THREADS

////////////////////////////////////////////////////////////////////////
// Base class that ThreadMap objects inherit from.  It is just a wrapper
// around a virtual destructor, so the map object can do its own
// clean-ups.
class Oxs_ThreadMapDataObject {
public:
  virtual ~Oxs_ThreadMapDataObject() {}
};

typedef map<String,Oxs_ThreadMapDataObject*> OXS_THREADMAP;

// Auto-initializing map object, used for thread local storage.
class Oxs_ThreadLocalMap {
private:
  Tcl_ThreadDataKey* mapkey;

#if OOMMF_THREADS
  OXS_THREADMAP* GetLockerPointer();
  /// Automatically initialized on first access in each thread.
#else
  static OXS_THREADMAP nothreads_locker;
  OXS_THREADMAP* GetLockerPointer() { return &nothreads_locker; }
#endif

  // Routine to remove all items from locker map, and call destructor on
  // each.  This routine should be called as part of problem release.
  // To help prevent abuse, this is a private function, accessible only
  // by friends --- it is not called by any Oxs_ThreadLocalMap member,
  // including in particular ~Oxs_ThreadLocalMap.
  static void DeleteLocker(Tcl_ThreadDataKey* mapkey);

  // Friend access for DeleteLocker
  friend Tcl_ThreadCreateType _Oxs_Thread_threadmain(ClientData);
#if OOMMF_THREADS
  friend void Oxs_ThreadTree::EndThreads();
#else
  friend class Oxs_ThreadTree;
#endif

public:
  Oxs_ThreadLocalMap(Tcl_ThreadDataKey* key) : mapkey(key) {}
  ~Oxs_ThreadLocalMap() {} // Note: Destructor does *not* delete
  /// any objects in locker map.  Instead, friend functions should
  /// call DeleteLocker during problem release.

  void AddItem(String name,Oxs_ThreadMapDataObject* obj); // Throws
  /// an exception if "name" is already in map

  Oxs_ThreadMapDataObject* GetItem(String name);  // Returns
  /// null pointer if "name" is not in map.

  Oxs_ThreadMapDataObject* UnmapItem(String name); // Removes
  /// item from map, and returns pointer to item.  Destructor
  /// for item is not called.  Throws an exception if "name"
  /// is not in map.

  void DeleteItem(String name); // Removes item from map and
  /// calls item's destructor.  Throws an exception if "name"
  /// is not in map.

};

// Base class that Oxs_Thread client should inherit from.
// This class is used in both threaded and non-threaded
// builds.
class Oxs_ThreadRunObj {
private:
  static Tcl_ThreadDataKey thread_data_map; // Key for thread-local
  // storage.  This storage is a map<string,void*> object.

  // Friend access for thread_data_map key.
  friend Tcl_ThreadCreateType _Oxs_Thread_threadmain(ClientData);
#if OOMMF_THREADS
  friend void Oxs_ThreadTree::EndThreads();
#else
  friend class Oxs_ThreadTree;
#endif

protected:
  // Allow children of Oxs_ThreadRunObj access to thread local storage.
  Oxs_ThreadLocalMap local_locker;

public:
  Oxs_ThreadRunObj() : local_locker(&thread_data_map) {}
  virtual void Cmd(int threadnumber,void* data) =0;
  virtual ~Oxs_ThreadRunObj() {}
};


#if !OOMMF_THREADS
// IMPLEMENTATION FOR NON-THREADED BUILDS //////////////////////////////

// Manager for Oxs_Thread objects.  Call Oxs_ThreadTree to
// start and stop threads.
class Oxs_ThreadTree {
public:
  Oxs_ThreadTree() {}
  ~Oxs_ThreadTree() {}
  static void EndThreads() {
    Oxs_ThreadLocalMap::DeleteLocker(&Oxs_ThreadRunObj::thread_data_map);
    // Note: thread_data_map is not used by no-threads version of
    // DeleteLocker, so in this case we may alternatively pass in a NULL
    // pointer.
    Oxs_ThreadError::ClearError();
  }

  void Launch(Oxs_ThreadRunObj& runobj,void* data) {
    // Launch "thread" and block until completion.
    runobj.Cmd(0,data);
  }

  void LaunchRoot(Oxs_ThreadRunObj& runobj,void* data) {
    // In non-threaded case, same as Launch
    Launch(runobj,data);
  }

  void RunOnThreadRange(int first_thread,int /* last_thread */,
                        Oxs_ThreadRunObj& runobj,void* data) {
    if(first_thread == 0) runobj.Cmd(0,data);
  }

  void DeleteLockerItem(String name) {
    Oxs_ThreadLocalMap locker(&Oxs_ThreadRunObj::thread_data_map);
    /// Actually, There is only one locker in no-thread case, so
    /// key value doesn't matter.
    if(locker.GetItem(name)) {
      locker.DeleteItem(name);
    }
  }

  static void GetStatus(String& results) { // For debugging.
    results = String("No threads.");
  }

  void Join() {}
};
#endif // !OOMMF_THREADS


////////////////////////////////////////////////////////////////////////
// THREAD BUSH AND TWIGS ///////////////////////////////////////////////
//
// This is an threading interface alternative to ThreadTree and friends.
// Unlike ThreadTree, which create semi-permanent threads that are
// reused thoughout problem lifetime, ThreadBush creates one-off threads
// for immediate use, which are terminated at the completion of one job.
//
// To use this interface, first create a child class off Oxs_ThreadTwig,
// which defines the job to be done through the virtual Task() member
// function call.  You can pack whatever data is needed into the child
// class, but only one instance of this class is used in all threads, so
// any differences between the threads has to be selected via the one
// parameter to Oxs_ThreadTwig::Task, namely "oc_thread_id".  When the
// child class is defined, then create one instance of the child class
// and one instance of the Oxs_ThreadBush class.  Call the RunAllThreads
// member of Oxs_ThreadBush with a pointer to the child instance.  Note
// that RunAllThreads is a regular member function, as opposed to a
// static member.  This is because each Oxs_ThreadBush instance has its
// own mutex and condition variables; this allows multiple
// Oxs_ThreadBush objects to be run concurrently.  Of course, in that
// case it is up to the clients to insure against deadlocks.

class Oxs_ThreadTwig {
  // This base class encapsulates a thread job for Oxs_ThreadBush
public:
  virtual void Task(int oc_thread_id) =0;
  virtual ~Oxs_ThreadTwig() {}
};

#if OOMMF_THREADS
Tcl_ThreadCreateProc _Oxs_ThreadBush_main; // Forward declaration
/// of helper function for Oxs_ThreadBush::RunAllThreads().
#endif // OOMMF_THREADS

class Oxs_ThreadBush
{
public:
  void RunAllThreads(Oxs_ThreadTwig* twig);
private:
  Oxs_ThreadControl task_stop;
  struct TwigBundle {
  public:
    Oxs_ThreadTwig* twig;
    Oxs_ThreadControl* thread_control;
    int oc_thread_id;
    TwigBundle() : twig(0), thread_control(0), oc_thread_id(-1) {}
  };
#if OOMMF_THREADS
  friend Tcl_ThreadCreateProc _Oxs_ThreadBush_main;
#endif // OOMMF_THREADS
};

#if !OOMMF_THREADS
inline void Oxs_ThreadBush::RunAllThreads(Oxs_ThreadTwig* twig)
{
  twig->Task(0);
  // All done!
}
#endif // OOMMF_THREADS


////////////////////////////////////////////////////////////////////////
// Oxs_StripedArray is an array wrapper that stripes the memory across
// memory nodes.  This is only meaningful on a NUMA (non-uniform memory
// access) system.  In other cases there is effectively only one node,
// so there is no striping.  However, in all cases the start of the
// array is aligned: If the requested memory size is larger than the memory
// page size, then the array start is aligned to a page boundary.  If
// the requested memory size is smaller than the memory page size, then
// the array start is aligned to the cache line size.
//   Note: This class is defined here, inside oxsthread.h, because
// it is accessed by the Oxs_JobControl class.

// Oxs_StripedArray helper thread function
class _Oxs_StripedArray_StripeThread : public Oxs_ThreadTwig {
  // This class represents the threads that do the actual initializing
  // of the allocated memory.  It is the act of initializing the memory
  // that causes it to be assigned to a specific memory node.  Provided
  // that the memory was allocated in 
private:
  char* const bufbase;
  const OC_INDEX full_size;
  const OC_INDEX strip_size;
  const int strip_count;
  const int augment_count;

   // Declare but don't define the following members
  _Oxs_StripedArray_StripeThread();
  _Oxs_StripedArray_StripeThread(const _Oxs_StripedArray_StripeThread&);
  _Oxs_StripedArray_StripeThread&
                       operator=(const _Oxs_StripedArray_StripeThread&);

public:
  _Oxs_StripedArray_StripeThread(char* new_bufbase,
                                 OC_INDEX new_full_size,
                                 OC_INDEX new_strip_size,
                                 int new_strip_count,
                                 int new_augment_count)
    : bufbase(new_bufbase), full_size(new_full_size),
      strip_size(new_strip_size),
      strip_count(new_strip_count),
      augment_count(new_augment_count) {}

  static void SplitAlgorithm(int threadid,
                             OC_INDEX fullsize,
                             OC_INDEX stripsize,
                             int stripcount,
                             int augmentcount,
                             OC_INDEX& mystart,
                             OC_INDEX& mystop) {
    assert(threadid<stripcount);
    mystart = threadid*stripsize;
    OC_INDEX mysize = stripsize;
    if(threadid<augmentcount) {
      mystart += threadid * OC_PAGESIZE;
      mysize += OC_PAGESIZE;
    } else {
      mystart += augmentcount * OC_PAGESIZE;
    }
    mystop = mystart + mysize;
    if(threadid < stripcount - 1) {
      if(mystop > fullsize) {
        // This shouldn't happen, unless full_size
        // is comparable to OC_PAGESIZE.
        mystop = fullsize;
      }
    } else {
      mystop = fullsize; // Last strip always gets all that's left
    }
  }

  void Task(int oc_thread_id) {
    // NB: The algorithm here needs to agree with the code
    // in Oxs_StripedArray<T>::GetStripPosition().
    OC_INDEX mystart,mystop;
    SplitAlgorithm(oc_thread_id,
                   full_size,strip_size,
                   strip_count,augment_count,
                   mystart,mystop);
    OC_INDEX mysize = mystop - mystart;
    if(mysize>0) {
      memset(bufbase+mystart,0,mysize);
    }
  }
};

template<class T> class Oxs_StripedArray
{
private:
  char* datablock;
  T* const arr;
  OC_INDEX arr_size;
  OC_INDEX strip_size;  // Nominal strip size, in bytes
  int strip_count;      // Number of strips
  int augment_count;    // Number of strips augmented with an extra page
  
  // Disable all copy operators by declaring them but not providing
  // a definition.
  Oxs_StripedArray(const Oxs_StripedArray<T>&);
  Oxs_StripedArray<T>& operator=(const Oxs_StripedArray<T>&);
public:
  void Free() {
    if(datablock) {
      // Implement "placement delete"
      while(arr_size>0) arr[--arr_size].~T(); // Explicit destructor call
      const_cast<T*&>(arr)=0;
      delete[] datablock;
      datablock=0;
    }
    const_cast<T*&>(arr) = 0;      // Safety
    arr_size = 0;
    strip_size = 0;
    strip_count = 0;
    augment_count = 0;
  }

  Oxs_StripedArray() : datablock(0), arr(0),
                       arr_size(0), strip_size(0),
                       strip_count(0), augment_count(0) {}
  ~Oxs_StripedArray() { Free(); }

  void Swap(Oxs_StripedArray<T>& other) {
    char* tdb = datablock;
    datablock = other.datablock;
    other.datablock=tdb;

    T* tarr = arr;
    const_cast<T*&>(arr) = other.arr;
    const_cast<T*&>(other.arr) = tarr;

    OC_INDEX tsize = arr_size;
    arr_size = other.arr_size;
    other.arr_size = tsize;

    OC_INDEX ssize = strip_size;
    strip_size = other.strip_size;
    other.strip_size = ssize;

    int tsc = strip_count;
    strip_count = other.strip_count;
    other.strip_count = tsc;

    int tac = augment_count;
    augment_count = other.augment_count;
    other.augment_count = tac;
  }

  void SetSize(OC_INDEX newsize);

  T* GetArrBase() const { return arr; }

  // Get offsets for a given strip, where
  //           0 <= strip_number < Oc_GetMaxThreadCount(),
  // in units of sizeof(T), relative to arr.
  void GetStripPosition(int number,
                        OC_INDEX& start_offset,
                        OC_INDEX& stop_offset) const;
};

template<class T> void Oxs_StripedArray<T>::GetStripPosition
(int number,
 OC_INDEX& start_offset,
 OC_INDEX& stop_offset) const
{ // NB: The algorithm here needs to agree with the code
  // in _Oxs_StripedArray_StripeThread::Task().

  assert(0<=number && number<strip_count);

  // First, compute strip position in bytes.
  OC_INDEX mystart, mystop;
  _Oxs_StripedArray_StripeThread::SplitAlgorithm
    (number,arr_size*sizeof(T),strip_size,strip_count,augment_count,
     mystart,mystop);

  // Convert units to sizeof(T).  If sizeof(T) does not divide
  // strip_size, then the convention in terms of T is that the strip
  // starts with the first full T element entirely inside the strip, and
  // ends with the last T element with at least one byte in the strip.
  start_offset = (mystart + sizeof(T) - 1)/sizeof(T);
  OC_INDEX endpt = arr_size;
  if(number<strip_count-1) { // Last thread gets (entire) remainder
    OC_INDEX testpt = (mystop + sizeof(T) - 1)/sizeof(T);
    if(testpt<endpt) endpt = testpt;
  }
  stop_offset = endpt;
}

template<class T> void Oxs_StripedArray<T>::SetSize
(OC_INDEX newsize)
{
  if(newsize == arr_size) return; // Nothing to do.
  Free();  // Release existing memory, if any.
  if(newsize<=0) return;

  // Determine strip sizes
  strip_count = Oc_GetMaxThreadCount();
  if(strip_count>1 && newsize*sizeof(T)>OC_PAGESIZE) {
    strip_size = (newsize/strip_count)*sizeof(T);
    strip_size -= strip_size % OC_PAGESIZE; // Reduce to page boundary

    // strip_size * strip_count is guaranteed not larger than newsize,
    // but may fall short by up to strip_count*(OC_PAGESIZE+sizeof(T))
    // We don't want this to bunch up
    // on last strip, so instead distribute excess across the strips,
    // starting at the beginning, one page each as needed.
    OC_INDEX leftovers = newsize*sizeof(T) - strip_size*strip_count;
    augment_count = static_cast<int>(leftovers/OC_PAGESIZE);
  } else {
    // Put everything into the first strip.
    strip_size = newsize;
    augment_count = 0;
  }

  // Allocate and initialize memory using "placement new".
  OC_INDEX alignment = 0;
  OC_INDEX blocksize = newsize*sizeof(T);
  if(blocksize>=OC_PAGESIZE) alignment = OC_PAGESIZE;
  else                       alignment = OC_CACHE_LINESIZE;
  blocksize += alignment - 1;

  datablock = new char[blocksize];
  if(!datablock) {
    // Handling for case where "new" is old broken kind that
    // doesn't throw an exception on error.
    char tbuf[1024];
    Oc_Snprintf(tbuf,sizeof(tbuf),
                "Oxs_StripedArray: Failure to allocate memory block"
                " of size %ld bytes",(long)blocksize);
    OXS_THROW(Oxs_NoMem,tbuf);
  }
  OC_INDEX alignoff = (OC_INDEX)datablock % alignment;
  if(alignoff>0) {
    alignoff = alignment - alignoff;
  }

  // Set memory policy (memory node selection) for datablock to
  // "first touch"
  Oc_SetMemoryPolicyFirstTouch(datablock+alignoff,newsize*sizeof(T));

  // Stripe data
  _Oxs_StripedArray_StripeThread stripe(datablock+alignoff,
                                        newsize*sizeof(T),
                                        strip_size,
                                        strip_count,augment_count);
  Oxs_ThreadBush thread_bush;
  thread_bush.RunAllThreads(&stripe);

  const_cast<T*&>(arr)
    = new(datablock + alignoff)T[newsize];  // Placement new
  /// Might want to put a fence here, to insure that arr gets set
  /// before arr_size, for any threads that are relying on arr_size to
  /// decide whether or not arr is usable.  In any event, it is
  /// probably safer to set "arr" before "arr_size".
  arr_size = newsize;
}


////////////////////////////////////////////////////////////////////////
// THREAD JOB CONTROL //////////////////////////////////////////////////

#if OOMMF_THREADS //////////////////////////////////////////////////////
template<class T> class Oxs_JobControl
{
private:
  struct JobData {
    OC_INDEX start;
    OC_INDEX stop;
    Oxs_Mutex mutex;
    int all_done;
    JobData() : start(-1), stop(-1), all_done(0) {}
  };

  JobData* threadjob;

  void ReassignJobs(OC_INT4m thread_id, OC_INDEX& start, OC_INDEX& stop);
  Oxs_Mutex rejuggle; // Used for reassigning unfinished jobs.

  OC_INDEX job_size;

  OC_INT4m thread_count;

  int all_done;

  enum { BASE_PIECE_COUNT=32 };

public:
  void Free() {
    if(threadjob != 0) delete[] threadjob;
    threadjob = 0;
    job_size = -1;
    thread_count = -1;
    all_done = 0;
  }

  Oxs_JobControl() : threadjob(0),
                     job_size(-1), thread_count(-1), all_done(0) {}

  ~Oxs_JobControl() { Free(); }

  void Init(OC_INT4m number_of_threads,
            const Oxs_StripedArray<T>* arrblock);

  void GetJob(OC_INT4m thread_id,
              OC_INDEX &start, OC_INDEX& stop);

};

template<class T> void Oxs_JobControl<T>::Init
(OC_INT4m number_of_threads,
 const Oxs_StripedArray<T>* arrblock)
{
  if(thread_count != number_of_threads ||
     threadjob==0) {
    // (re)Alloc.  Otherwise, reuse old space.
    Free();
    thread_count = number_of_threads;
    threadjob = new JobData[thread_count];
  }

  for(OC_INT4m i=0;i<thread_count;++i) {
    arrblock->GetStripPosition(i,threadjob[i].start,threadjob[i].stop);
    threadjob[i].all_done = 0;
  }

  all_done = 0;
}

template<class T> inline void Oxs_JobControl<T>::GetJob
(OC_INT4m thread_id,
 OC_INDEX &start,
 OC_INDEX& stop)
{
  JobData& job = threadjob[thread_id];

  OC_INDEX test_start,test_stop;

  job.mutex.Lock();
  test_start = job.start;
  job.start = test_stop  = job.stop;
  int fini = job.all_done;
  job.mutex.Unlock();

  if(test_start<test_stop) {
    // Job found in bucket for thread_id
    start = test_start;
    stop  = test_stop;
    return;
  }

  // Otherwise, all jobs for thread_id have been completed.
  // See if we can find some jobs in a different thread bucket.
  if(fini) {
    // No; all jobs assigned
    start = stop = -1;
  } else {
    ReassignJobs(thread_id,start,stop);
  }

  return;
}

template<class T> void Oxs_JobControl<T>::ReassignJobs
(OC_INT4m thread_id,
 OC_INDEX &start,
 OC_INDEX& stop)
{
  // Reassignments not currently supported
  start = stop = -1;

  JobData& job = threadjob[thread_id];
  job.mutex.Lock();
  job.all_done = 1;
  job.mutex.Unlock();

  return;
}

#else  // !OOMMF_THREADS ///////////////////////////////////////////////

template<class T> class Oxs_JobControl
{ // Unthreaded version
private:
  struct JobRange {
    OC_INDEX start;
    OC_INDEX stop;
    JobRange() : start(0), stop(0) {}
  };
  JobRange threadjobs;
public:
  void Free() {
    threadjobs.start = threadjobs.stop = 0;
  }

  Oxs_JobControl() {}
  ~Oxs_JobControl() { Free(); }

  void Init(OC_INT4m
#ifndef NDEBUG
            number_of_threads
#endif
            , const Oxs_StripedArray<T>* arrblock) {
    assert(number_of_threads == 1);
    arrblock->GetStripPosition(0,threadjobs.start,threadjobs.stop);
  }

  void GetJob(OC_INT4m /* thread_id */,
              OC_INDEX &start, OC_INDEX& stop) {
    start = threadjobs.start;
    stop  = threadjobs.stop;
    threadjobs.start = threadjobs.stop = -1;
  }
};
#endif // OOMMF_THREADS


#endif // _OXS_THREAD
