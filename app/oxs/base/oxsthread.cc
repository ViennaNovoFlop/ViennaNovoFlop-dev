/* FILE: oxsthread.cc              -*-Mode: c++-*-
 *
 * OXS thread support.
 *
 */

#include "oxsexcept.h"
#include "oxsthread.h"

/* End includes */

// Thread-safe fprintf
void Oxs_ThreadPrintf(FILE* fptr,const char *format, ...) {
  static Oxs_Mutex oxs_threadprintf_mutex;
  oxs_threadprintf_mutex.Lock();

  va_list arg_ptr;
  va_start(arg_ptr,format);
  vfprintf(fptr,format,arg_ptr);
  va_end(arg_ptr);

  fflush(fptr);

  oxs_threadprintf_mutex.Unlock();
}

#if OOMMF_THREADS
Oxs_Mutex Oxs_ThreadError::mutex;
int Oxs_ThreadError::error(0);
String Oxs_ThreadError::msg;
#endif

// Notes on Tcl_MutexLock, Tcl_ConditionNotify, and Tcl_ConditionWait:
//  1) A lock should be held on a mutex before issuing a
//     notify or performing a wait.
//  2) When a wait is issued, the mutex lock is released.
//     When the thread awakens, it automatically reclaims
//     the mutex lock.
//  3) After issuing a notify, the code needs to explicitly
//     release the mutex lock.  Threads waiting on the
//     notify condition will not begin executing until
//     they can grab the mutex (as explained in 2).
//  4) After a notify is issued, there is no guarantee that
//     upon release of the mutex that the next mutex lock
//     will go to a thread waiting at a condition wait.  In
//     particular, consider this scenario:
//
//        Initially, count=0, thread 1 has a lock on mutex.
//
//        Thread 1:
//           Tcl_ConditionWait(&cond,&mutex,NULL);
//           cout << count;
//           Tcl_MutexUnlock(mutex);
//
//        Threads 2 & 3:
//           Tcl_Lock(&mutex);
//           ++count;
//           Tcl_ConditionNotify(&cond);
//           Tcl_Unlock(&mutex);
//
//     The value for count output by thread 1 may be either 1 or 2.
//     For example, suppose thread 2 succeeds in acquiring mutex first.
//     Then thread 2 increments count to 1, and notifies on cond.  Then
//     if thread 1 acquires mutex next, count will be 1.  However, it
//     is just as likely that thread 3 will acquire mutex, in which case
//     thread 3 will increment count to 2.  Assuming that no other
//     threads are waiting for mutex, then thread 1 will finally acquire
//     mutex, and the value for count at this time will be 2.
//        Similarly, even if there is only one thread waiting on a
//     particular condition, the number of times that condition is
//     notified may be greater than the number of times the waiting
//     thread is awaken.  For example, if Thread 1 above is changed to
//
//        Thread 1: (BAD BAD BAD!)
//           Tcl_ConditionWait(&cond,&mutex,NULL);
//           cout << count;
//           Tcl_ConditionWait(&cond,&mutex,NULL);
//           cout << count;
//           Tcl_MutexUnlock(mutex);
//
//     then this may hang at the second wait.  The correct way to
//     wait on multiple notifies is like so:
//
//        Thread 1: (GOOD)
//           while(count<2) { // Assume mutex held coming in
//              Tcl_ConditionWait(&cond,&mutex,NULL);
//           }
//           count << count;
//           Tcl_MutexUnlock(mutex);
//


////////////////////////////////////////////////////////////////////////
// Thread local storage.  Each thread gets its own instance.

#if OOMMF_THREADS //    --- threads version ----------------------------

void Oxs_ThreadLocalMap::DeleteLocker(Tcl_ThreadDataKey* mapkey)
{ // Destroy locker items, version for threaded builds
  void** vtmp = (void**)Tcl_GetThreadData(mapkey,sizeof(void*));
  assert(vtmp != NULL);
  if(*vtmp != NULL) {
    OXS_THREADMAP* lockerptr = static_cast<OXS_THREADMAP*>(*vtmp);
    OXS_THREADMAP::iterator it;
    for(it = lockerptr->begin();it!=lockerptr->end();++it) {
      // Delete object pointed to by Oxs_ThreadMapDataObject*
      delete it->second;
    }
    delete lockerptr;
    *vtmp = NULL; // Safety
  }
}

OXS_THREADMAP*
Oxs_ThreadLocalMap::GetLockerPointer()
{ // Returns pointer to thread-specific instance of map.
  // If map does not exist, then automatically creates it.
  // (Threaded version.)

  void** vtmp = (void**)Tcl_GetThreadData(mapkey,sizeof(void*));
  assert(vtmp != NULL);

  OXS_THREADMAP* lockerptr = NULL;

  if(*vtmp == NULL) {
    // Map doesn't exist.  Create it.
    lockerptr = new OXS_THREADMAP;
    *vtmp = lockerptr;
  } else {
    // Map already exists.  Return pointer.
    lockerptr = static_cast<OXS_THREADMAP*>(*vtmp);
  }

  return lockerptr;
}

#else // OOMMF_THREADS   --- no-threads version ------------------------

OXS_THREADMAP Oxs_ThreadLocalMap::nothreads_locker;

void Oxs_ThreadLocalMap::DeleteLocker(Tcl_ThreadDataKey*)
{ // Destroy locker items, version for non-threaded builds
  OXS_THREADMAP::iterator it;
  for(it = nothreads_locker.begin();it!=nothreads_locker.end();++it) {
    // Delete object pointed to by Oxs_ThreadMapDataObject*
    delete it->second;
  }
  nothreads_locker.clear();
}

// No-threads version of GetLockerPointer() is inlined in oxsthread.h

#endif // OOMMF_THREADS

void
Oxs_ThreadLocalMap::AddItem
(String name,
 Oxs_ThreadMapDataObject* obj)
{ // Throws an exception if "name" is already in map
  OXS_THREADMAP& locker = *GetLockerPointer();

  OXS_THREADMAP::iterator it = locker.find(name);
  if(it != locker.end()) {
    OXS_THROW(Oxs_ProgramLogicError, String("Map item \"") + name
              + String("\" already in map."));
  }
  // New value
  locker[name] = obj;
}

Oxs_ThreadMapDataObject*
Oxs_ThreadLocalMap::GetItem
(String name)
{ // Returns null pointer if "name" is not in map.
  OXS_THREADMAP& locker = *GetLockerPointer();
  OXS_THREADMAP::iterator it = locker.find(name);
  if(it == locker.end()) {
    return NULL;
  }
  return it->second;
}

Oxs_ThreadMapDataObject*
Oxs_ThreadLocalMap::UnmapItem
(String name)
{ // Removes item from map, and returns pointer to item.  Destructor
  // for item is not called.  Throws an exception if "name" is not in
  // map.
  OXS_THREADMAP& locker = *GetLockerPointer();
  OXS_THREADMAP::iterator it = locker.find(name);
  if(it == locker.end()) {
    OXS_THROW(Oxs_ProgramLogicError, String("Map item \"") + name
              + String("\" not in map; can't be unmapped."));
  }
  Oxs_ThreadMapDataObject* objptr = it->second;
  locker.erase(it);
  return objptr;
}

void
Oxs_ThreadLocalMap::DeleteItem
(String name)
{ // Removes item from map and calls item's destructor.  Throws an
  // exception if "name" is not in map.
  Oxs_ThreadMapDataObject* objptr = UnmapItem(name);
  delete objptr;
}

Tcl_ThreadDataKey Oxs_ThreadRunObj::thread_data_map;

#if OOMMF_THREADS
////////////////////////////////////////////////////////////////////////
// Mainline for actual (OS) threads
// Note: The return type differs by OS; on unix it is void, on
//       Windows it is unsigned.

Tcl_ThreadCreateProc _Oxs_Thread_threadmain;
Tcl_ThreadCreateType _Oxs_Thread_threadmain(ClientData clientdata)
{ // NB: Be careful, this routine and any routine called from inside
  //     must be re-entrant!

  // Cast and dereference clientdata
  // NB: It is assumed that the referenced object remains valid
  // throughout the life of the thread.
  Oxs_Thread& othread = *static_cast<Oxs_Thread *>(clientdata);

  // Break out references to object data.  Copy thread_number on the
  // assumption that it is fixed throughout object lifetime.
  // The stop, runobj, and data objects may change between calls
  // (i.e., wakeups); we hold a references to these pointers to keep
  // this up-to-date.
  Oxs_ThreadControl& start  = othread.start;
  int thread_number         = othread.thread_number;
  Oxs_ThreadControl* &stop  = othread.stop;
  Oxs_ThreadRunObj* &runobj = othread.runobj;
  void* &data               = othread.data;

  Oc_NumaRunOnNode(thread_number);

  // Notify parent Oxs_Thread that we are ready and waiting
  start.Lock();
  start.count = 1; // Ready to wait
  start.Notify();

  // Error handling
  char errbuf[256];
  String threadstr;
  Oc_Snprintf(errbuf,sizeof(errbuf),"\nException thrown in thread %d",
              thread_number);
  threadstr = errbuf;

  // Action loop
  while(1) {
    while(start.count>0) {
      start.Wait(0);
    }
    if(0==runobj) break;
    try {
      runobj->Cmd(thread_number,data);
    } catch(String& smsg) {
      Oxs_ThreadError::SetError(smsg + threadstr);
    } catch(const char* cmsg) {
      Oxs_ThreadError::SetError(String(cmsg) + threadstr);
    } catch(...) {
      Oc_Snprintf(errbuf,sizeof(errbuf),
                  "Unrecognized exception thrown in thread %d\n",
                  thread_number);
      Oxs_ThreadError::SetError(String(errbuf));
    }
    stop->Lock();
    --stop->count; // NB: stop->mutex must be acquired before
                  ///     changing stop->count.
    stop->Notify();
    start.count = 1; // Ready to wait.  Set this before freeing stop
    /// mutex, so that the thread signals as ready if tested from
    /// inside Oxs_ThreadTree::Launch.
    stop->Unlock();
  }
  start.count = 1;
  start.Notify();
  start.Unlock();

  // Clean-up thread-local storage
  Oxs_ThreadLocalMap::DeleteLocker(&Oxs_ThreadRunObj::thread_data_map);

  Tcl_ExitThread(TCL_OK);
  TCL_THREAD_CREATE_RETURN;
}

////////////////////////////////////////////////////////////////////////
// Oxs_Thread class, which is control object for actual thread.

Oxs_Thread::Oxs_Thread
(int thread_number_x)
  : thread_number(thread_number_x),
    stop(0), runobj(0), data(0)
{
  start.Lock();
  if(TCL_OK
     != Tcl_CreateThread(&thread_id, _Oxs_Thread_threadmain, this,
                         TCL_THREAD_STACK_DEFAULT, TCL_THREAD_NOFLAGS)) {
    start.Unlock();
    OXS_THROW(Oxs_BadResourceAlloc,"Thread creation failed."
              "(Is the Tcl library thread enabled?)");
  }
  start.count = 0;
  while(0 == start.count) {
    // I don't know that this loop is really necessary, because
    // this is the only thread waiting on start.cond and only
    // one thread is notifying on start.cond, but all the
    // references say to do this, and it shouldn't hurt.
    start.Wait(0);
  }
  start.Unlock();
}

Oxs_Thread::~Oxs_Thread()
{
  start.Lock();
  stop = 0;
  runobj = 0;
  data = 0;
  start.count = 0;
  start.Notify();
  while(0 == start.count) {
    // Is this loop necessary?  Cf. notes in check_start
    // loop in Oxs_Thread constructor.
    start.Wait(0);
  }
  start.Unlock();
  // Note: start.cond and start.mutex are automatically finalized
  // inside the Oxs_ThreadControl destructor
}

void Oxs_Thread::RunCmd
(Oxs_ThreadControl& stop_x,
 Oxs_ThreadRunObj& runobj_x,
 void* data_x)
{
  // runobj == NULL is used as a message to threadmain to
  // terminate.  Don't expose this implementation detail
  // to the user.  Instead, raise an error.
  if(0 == &runobj_x) {
    OXS_THROW(Oxs_BadParameter,"NULL thread runobj reference.");
  }
  start.Lock();
  stop = &stop_x;
  runobj = &runobj_x;
  data = data_x;
  start.count = 0; // Run signal
  start.Notify();
  start.Unlock();
}

void Oxs_Thread::GetStatus(String& results)
{ // For debugging.
  char buf[1024];
  int start_locked,start_count;
  start.GetStatus(start_locked,start_count);
  Oc_Snprintf(buf,sizeof(buf),"start (%p) mutex %s/count %d",
              &start,(start_locked ? "locked" : "free"),start_count);
  results += String(buf);
  
  if(stop) {
    int stop_locked,stop_count;
    stop->GetStatus(stop_locked,stop_count);
    Oc_Snprintf(buf,sizeof(buf),", stop (%p) mutex %s/count %d",
                stop,(stop_locked ? "locked" : "free"),stop_count);
    results += String(buf);
  }
}

////////////////////////////////////////////////////////////////////////
// ThreadTree class, manager for Oxs_Thread objects
// Note: stop.mutex is grabbed by Oxs_ThreadTree upon the
//       first Launch command, and is held until Join
//       is called.  Each tree has its own stop control,
//       but the thread array is shared by all trees.
Oxs_Mutex Oxs_ThreadTree::launch_mutex;
std::vector<Oxs_Thread*> Oxs_ThreadTree::threads;

void Oxs_ThreadTree::Launch(Oxs_ThreadRunObj& runobj,void* data)
{ // Launched threads, one at a time.
  if(0==threads_unjoined) stop.Lock();

  launch_mutex.Lock(); // Make sure no new threads get
  /// run while we're looking for a free thread.
  int free_thread=0;
  while(size_t(free_thread)<threads.size()) {
    if(!threads[free_thread]->IsRunning()) break;
    ++free_thread;
  }
  if(size_t(free_thread)==threads.size()) {
    // All threads in use; create a new thread.
    // Note: The thread_number is one larger than the index into
    // the "threads" vector, because thread_number 0 is reserved
    // for the main (parent) thread, which is not referenced in
    // the "threads" vector.  (The "threads" vector only holds
    // child threads.)
    threads.push_back(new Oxs_Thread(free_thread+1));
Oxs_ThreadPrintf(stderr,"Starting child thread #%d\n",free_thread+1); /**/
  }
  ++threads_unjoined;
  ++stop.count;
  threads[free_thread]->RunCmd(stop,runobj,data);
  launch_mutex.Unlock();
}

void Oxs_ThreadTree::Join()
{ // Joins all running threads, and releases mutex_stop
  String errmsg;

  if(0==threads_unjoined) {
    if(Oxs_ThreadError::IsError(&errmsg)) OXS_THROW(Oxs_BadThread,errmsg);
    return;
  }

  // Note: If threads_unjoined>0, then we should be
  //       holding stop.mutex
  while(stop.count>0) {
    stop.Wait(0);
  }
  threads_unjoined = 0;
  stop.Unlock();
  if(Oxs_ThreadError::IsError(&errmsg)) OXS_THROW(Oxs_BadThread,errmsg);
}

void Oxs_ThreadTree::LaunchRoot(Oxs_ThreadRunObj& runobj,void* data)
{ // Runs command on root (i.e., 0) thread, and calls Join() to wait for
  // all children to finish, with error handling.
  try {
    runobj.Cmd(0,data);
  } catch(...) {
    if(threads_unjoined>0) {
      try {
        Join();
      } catch(...) {
        // Ignore for now.  More generally, should probably
        // catch Oxs_BadThread exceptions and append them
        // to exception from root thread.
      }
    }
    throw;
  }

  // Wait for children to finish
  if(threads_unjoined>0) Join();

  String errmsg;
  if(Oxs_ThreadError::IsError(&errmsg)) OXS_THROW(Oxs_BadThread,errmsg);
}

void
Oxs_ThreadTree::RunOnThreadRange
(int first_thread,
 int last_thread,
 Oxs_ThreadRunObj& runobj,
 void* data)
{ // Run runobj on threads first_thread through last_thread, inclusive.
  // Threads are numbered from 0, where 0 is the main thread (i.e., the
  // thread that calls into here).  New threads will be created as
  // needed to satisfy the thread range.  If last_thread is -1, then the
  // top thread is taken to be the highest existing thread --- the
  // first_thread request is satisfied first, so that threads may be
  // created even if last_thread == -1.  To run on all existing threads
  // (including thread 0), set first_thread = 0 and last_thread = -1.
  // NOTE: If there are any threads running at entry, then this routine
  // will block until all are joined.

  String errmsg;
  if(Oxs_ThreadError::IsError(&errmsg)) OXS_THROW(Oxs_BadThread,errmsg);

  launch_mutex.Lock(); // Make sure nothing new gets launched
  /// before this routine is complete.

  // Join any running threads, and grab stop.mutex
  Join();
  stop.Lock();

  // Create new threads as needed.
  int i = threads.size();
  while(i<first_thread) {
    threads.push_back(new Oxs_Thread(++i));
  }
  while(i<last_thread) {
    threads.push_back(new Oxs_Thread(++i));
  }
  if(last_thread == -1) last_thread = threads.size();
  if(first_thread<0) first_thread = 0;

  // Run on child threads
  for(i=first_thread;i<=last_thread;++i) {
    if(i<1) continue;
    ++threads_unjoined;
    ++stop.count;
    threads[i-1]->RunCmd(stop,runobj,data);
  }

  if(first_thread==0) {
    // Run on main thread
    try {
      runobj.Cmd(0,data);
    } catch(...) {
      if(threads_unjoined>0) {
        Join();
      } else {
        stop.Unlock();  // Make sure this gets released.
      }
      launch_mutex.Unlock();
      throw;
    }
  }

  // Wait for children to finish
  if(threads_unjoined>0) {
    Join();
  } else {
    stop.Unlock();  // Make sure this gets released.
  }

  launch_mutex.Unlock();
  if(Oxs_ThreadError::IsError(&errmsg)) OXS_THROW(Oxs_BadThread,errmsg);
}

// Support class for Oxs_ThreadTree::DeleteLockerItem
class _Oxs_ThreadTree_DLI : public Oxs_ThreadRunObj {
private:
  String item_name;
public:
  _Oxs_ThreadTree_DLI(String name) : item_name(name) {}
  void Cmd(int /* threadnumber */, void* /* data */) {
    if(local_locker.GetItem(item_name)) {
      local_locker.DeleteItem(item_name);
    }
  }
};

void Oxs_ThreadTree::DeleteLockerItem(String name)
{ // Deletes local locker item
  _Oxs_ThreadTree_DLI foo(name);
  RunOnThreadRange(0,-1,foo,0);
}


void
Oxs_ThreadTree::GetStatus(String& results)
{ // This routine is intended fo debugging purposes only.  It returns
  // string detailing "probable" status for all mutexes Oxs_ThreadTree
  // knows about, in human readable format.
  results.clear();
  results += String("launch_mutex ");
  if(launch_mutex.IsLockProbablyFree()) {
    results += String("free");
  } else {
    results += String("locked");
  }

  char buf[1024];
  for(size_t i=0;i<threads.size();++i) {
    Oc_Snprintf(buf,sizeof(buf),"\nThread %2d: ",i);
    results += String(buf);
    threads[i]->GetStatus(results);
  }

  return;
}

Oxs_ThreadTree::~Oxs_ThreadTree()
{ // Note: The mutex and condition variables in stop are automatically
  // finalized.

  // If we get here during error processing, it is possible that there
  // are unjoined threads running.  Try to kill these.
  // NB: This code is not well-tested and is probably holey.
  if(threads_unjoined>0) {
    // Note: If threads_unjoined>0, then we should be
    //       holding stop.mutex
    while(threads_unjoined>0 && stop.count>0) {
      stop.Wait(10000000);  // Wait time in microseconds
      --threads_unjoined;
    }
    threads_unjoined = 0;
    stop.Unlock();
  }
  if(stop.count>0) {
    Oxs_ThreadPrintf(stderr,
                     "WARNING: Oxs_ThreadTree destruction with %d"
                     " unjoined threads\n",stop.count);
  }
}

void Oxs_ThreadTree::EndThreads()
{ // Note: This is a static member function of Oxs_ThreadTree
  for(size_t i=0;i<threads.size();++i) {
Oxs_ThreadPrintf(stderr,"Ending child thread #%d\n",i+1); /**/
    delete threads[i];
  }
  threads.clear();
  launch_mutex.Reset();

  // Destroy local locker data in thread 0 (mainline)
  Oxs_ThreadLocalMap::DeleteLocker(&Oxs_ThreadRunObj::thread_data_map);
  Oxs_ThreadError::ClearError();
}

#endif // OOMMF_THREADS


////////////////////////////////////////////////////////////////////////
// THREAD BUSH AND TWIGS ///////////////////////////////////////////////

#if OOMMF_THREADS

Tcl_ThreadCreateType _Oxs_ThreadBush_main(ClientData clientdata)
{ // Helper function for Oxs_ThreadBush::RunAllThreads().  There
  // is a declaration of this function in oxsthread.h

  // Cast and dereference clientdata
  // NB: It is assumed that the referenced object remains valid
  // throughout the life of the thread.
  Oxs_ThreadBush::TwigBundle& bundle
    = *static_cast<Oxs_ThreadBush::TwigBundle *>(clientdata);

  const int oc_thread_id = bundle.oc_thread_id;
  Oxs_ThreadTwig& twig = *(bundle.twig);
  Oxs_ThreadControl& thread_control = *(bundle.thread_control);

  Oc_NumaRunOnNode(oc_thread_id);

  // Error handling
  char errbuf[256];
  String threadstr;
  Oc_Snprintf(errbuf,sizeof(errbuf),
              "\nException thrown in thread %d",
              oc_thread_id);
  threadstr = errbuf;

  // Action loop
  try {
    twig.Task(oc_thread_id);
  } catch(String& smsg) {
    Oxs_ThreadError::SetError(smsg + threadstr);
  } catch(const char* cmsg) {
    Oxs_ThreadError::SetError(String(cmsg) + threadstr);
  } catch(...) {
    Oc_Snprintf(errbuf,sizeof(errbuf),
                "Unrecognized exception thrown in thread %d\n",
                oc_thread_id);
    Oxs_ThreadError::SetError(String(errbuf));
  }

  thread_control.Lock();
  --thread_control.count;
  thread_control.Notify();
  thread_control.Unlock();
  // NOTE: After this point, the bundle ptr may go invalid.

  // All done.
  Tcl_ExitThread(TCL_OK);

  TCL_THREAD_CREATE_RETURN; // pro forma
}

void Oxs_ThreadBush::RunAllThreads(Oxs_ThreadTwig* twig)
{
  const int thread_count = Oc_GetMaxThreadCount();
  task_stop.count = thread_count;  // Do we need to Lock() here???

  TwigBundle* bundle = new TwigBundle[thread_count];

  // Launch all threads (including thread 0)
  for(int i=0;i<thread_count;++i) {
    Tcl_ThreadId tcl_id;  // Not used
    bundle[i].twig = twig;
    bundle[i].thread_control = &task_stop;
    bundle[i].oc_thread_id = i;
    if(TCL_OK
       != Tcl_CreateThread(&tcl_id, _Oxs_ThreadBush_main, &bundle[i],
                           TCL_THREAD_STACK_DEFAULT, TCL_THREAD_NOFLAGS)) {
      OXS_THROW(Oxs_BadResourceAlloc,"Thread creation failed."
                "(Is the Tcl library thread enabled?)");
    }
  }

  // Wait for all threads to finish
  task_stop.Lock();
  while(task_stop.count>0) task_stop.Wait();
  task_stop.Unlock();

  delete[] bundle; // Child threads must not access bundle array
  /// after they decrement task_stop.count.

  // Error check
  String errmsg;
  if(Oxs_ThreadError::IsError(&errmsg)) OXS_THROW(Oxs_BadThread,errmsg);

  // All done!
}

#endif
