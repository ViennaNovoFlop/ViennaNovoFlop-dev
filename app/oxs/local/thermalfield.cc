/* FILE: thermalfield.cc            -*-Mode: c++-*-
 * 24.1.2013 Thomas Windbacher (TW)
 * Thermal Field derived from Uniaxial Anisotropy and Oxs_Energy class.
 * See: Finocchio, J. Appl. Phys. 99 doi:10.1063/1.2177049
 */

#include "oc.h"
#include "nb.h"
#include "threevector.h"
#include "director.h"
#include "simstate.h"
#include "ext.h"
#include "key.h"
#include "mesh.h"
#include "meshvalue.h"
#include "uniformscalarfield.h"
#include "uniformvectorfield.h"
#include "thermalfield.h"
#include "energy.h"		// Needed to make MSVC++ 5 happy

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(Oxs_ThermalField);

/* End includes */


// Constructor
Oxs_ThermalField::Oxs_ThermalField(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Oxs_ChunkEnergy(name,newdtr,argstr), mesh_id(0)
{
  // Process arguments
  if(HasInitValue("T")){
    OXS_GET_INIT_EXT_OBJECT("T",Oxs_ScalarField,T_init);
  } else {
      // Set default temperature scalar field to 300 K
      // Oxs_Ext* foo=MakeNew("Oxs_UniformScalarField",director,"value 300.0");
      // T_init.SetAsOwner(dynamic_cast<Oxs_ScalarField*>(foo));
      throw Oxs_ExtError(this,"Error Oxs_ThermalField: No temperature T defined.\n");
  }
 if(HasInitValue("alpha")){
  OXS_GET_INIT_EXT_OBJECT("alpha",Oxs_ScalarField,alpha_init);
  } else {
      // Set default alpha to 0.01
      // Oxs_Ext* foo=MakeNew("Oxs_UniformScalarField",director,"value 0.01");
      // alpha_init.SetAsOwner(dynamic_cast<Oxs_ScalarField*>(foo));
      throw Oxs_ExtError(this,"Error Oxs_ThermalField: No damping coefficient alpha defined.\n");
  } 

  D = GetRealInitValue("D",1.0);
  //in_time_step = GetRealInitValue("min_time_step",1.e-14);
  time_step = GetRealInitValue("time_step",1.e-14);
  fixed_time_step = GetIntInitValue("fixed_time_step",0);
  VerifyAllInitArgsUsed();
 // start the clock used for seeding
  beginning = myclock::now();
 // start random number generators and normal distributions for x, y, and z direction
 // obtain a seed from the timer and initilaize random number setup
  dx = myclock::now() - beginning;
  unsigned seedx = dx.count();

  generatorx.seed(seedx);
  
  distributionx = std::normal_distribution<OC_REAL8m> (0.0,D);
  distributionx.reset();


}

OC_BOOL Oxs_ThermalField::Init()
{
  mesh_id = 0;
  T.Release();
  alpha.Release();
 // start the clock used for seeding
  beginning = myclock::now();
  dx = myclock::now() - beginning;
  unsigned seedx = dx.count();

  distributionx.reset();
  generatorx.seed(seedx);

  return Oxs_ChunkEnergy::Init();
}


void Oxs_ThermalField::IntegEnergy
(const Oxs_SimState& state,
 Oxs_ComputeEnergyDataThreaded& ocedt,
 Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
 OC_INDEX node_start,OC_INDEX node_stop
 ) const
{

 // start the clock used for seeding
 //beginning = myclock::now();

  const Oxs_Mesh* mesh = state.mesh;
  const Oxs_MeshValue<OC_REAL8m>& Ms_inverse = *(state.Ms_inverse);
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);
  const Oxs_MeshValue<ThreeVector>& spin = state.spin;
  Nb_Xpfloat energy_sum  = 0.0;
  const OC_REAL8m KB     = 1.3806488e-23;  // J/K
  const OC_REAL8m GAMMAE = 1.760859708e11; // 1/(s*T)
  const OC_REAL8m GAMMA  = 2.211e5; // 1/(s*T)
  const OC_REAL8m sign   = -1.;
        OC_REAL8m dt     = -1.; 
  alpha_init->FillMeshValue(mesh,alpha);
  T_init->FillMeshValue(mesh,T);

  if(fixed_time_step != 0) {
     dt = time_step;
    } else {
      dt =  state.last_timestep;
      if(state.last_timestep < time_step) { dt = time_step; }
    }
// if(mesh_id != state.mesh->Id()) {                                                                                                
// This is either the first pass through, or else mesh                                                                         
// has changed.                                                                                                                
//  mesh_id = 0;

// These values do not change wihtin the subsequent loop so it is more efficient to calculate only once
    OC_REAL8m pre_factor  = 2.0;
    pre_factor           *= KB;
    //pre_factor           /= GAMMAE;
    pre_factor           /= GAMMA;
    pre_factor           /= dt;
    pre_factor            = sqrt(pre_factor);

 for(OC_INDEX i=node_start;i<node_stop;++i) {


    	OC_REAL8m cell_alpha  = alpha[i];
    	OC_REAL8m cell_T      = T[i]; 
    	OC_REAL8m cell_Volume = mesh->Volume(i); 
    	ThreeVector         m = spin[i]; 
    	ThreeVector       mxH = m;

    	OC_REAL8m cell_factor    = cell_alpha;
        cell_factor   *= cell_T;
       	cell_factor   /= cell_Volume;
       	cell_factor   *= Ms_inverse[i];
       	cell_factor   /= (1.+cell_alpha*cell_alpha);
       	cell_factor    = sqrt(cell_factor);
       	cell_factor   *= pre_factor;

    	//if(cell_factor==0.0) {
      	//if(ocedt.energy) (*ocedt.energy)[i] = 0.0;
      	//if(ocedt.H)      (*ocedt.H)[i].Set(0.,0.,0.);
      	//if(ocedt.mxH)    (*ocedt.mxH)[i].Set(0.,0.,0.);
    	//if(ocedt.energy)       (*ocedt.energy)[i] = 0.;
    	//if(ocedt.energy_accum) (*ocedt.energy_accum)[i] += 0.;
      	//continue;
    	//}

    	//ThreeVector xi = ThreeVector(distributionx(generatorx),distributiony(generatory),distributionz(generatorz)); 
    	ThreeVector xi = ThreeVector(distributionx(generatorx),distributionx(generatorx),distributionx(generatorx)); 

    	ThreeVector     H  =  xi*cell_factor;
                      mxH ^=  H;
    	OC_REAL8m   mdotH  = (m*H); 
    	OC_REAL8m   ei     = sign;
                    ei    *= MU0;
                    ei    *= Ms[i];
                    ei    *= mdotH;
                

  
    	if(ocedt.energy)       (*ocedt.energy)[i] = ei;
    	if(ocedt.energy_accum) (*ocedt.energy_accum)[i] += ei;
    	if(ocedt.H)       (*ocedt.H)[i] = H;
    	if(ocedt.H_accum) (*ocedt.H_accum)[i] += H;
    	if(ocedt.mxH)       (*ocedt.mxH)[i] = mxH;
    	if(ocedt.mxH_accum) (*ocedt.mxH_accum)[i] += mxH;

    	energy_sum += ei*cell_Volume;
  }

//  mesh_id = state.mesh->Id();
  ocedtaux.energy_total_accum += energy_sum.GetValue();
//}
}


void Oxs_ThermalField::ComputeEnergyChunk
(const Oxs_SimState& state,
 Oxs_ComputeEnergyDataThreaded& ocedt,
 Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
 OC_INDEX node_start,OC_INDEX node_stop,
 int threadnumber
 ) const
{
  if(node_stop>state.mesh->Size() || node_start>node_stop) {
    throw Oxs_ExtError(this,"Programming error:"
                       " Invalid node_start/node_stop values");
  }

  if(mesh_id !=  state.mesh->Id()) {
    // This is either the first pass through, or else mesh
    // has changed.  Initialize/update data fields.
    // NB: At a lower level, this may potentially involve calls back
    // into the Tcl interpreter.  Per Tcl spec, only the thread
    // originating the interpreter is allowed to make calls into it, so
    // only threadnumber == 0 can do this processing.  Any other thread
    // must block until that processing is complete.
    thread_control.Lock();
    if(Oxs_ThreadError::IsError()) {
      if(thread_control.count>0) {
        // Release a blocked thread
        thread_control.Notify();
      }
      thread_control.Unlock();
      return; // What else?
    }
    if(threadnumber != 0) {
      if(mesh_id != state.mesh->Id()) {
        // If above condition is false, then the main thread came
        // though and initialized everything between the time of
        // the previous check and this thread's acquiring of the
        // thread_control mutex; in which case, "never mind".
        // Otherwise:
        ++thread_control.count; // Multiple threads may progress to this
        /// point before the main thread (threadnumber == 0) grabs the
        /// thread_control mutex.  Keep track of how many, so that
        /// afterward they may be released, one by one.  (The main
        /// thread will Notify control_wait.cond once; after that
        /// as each waiting thread is released, the newly released
        /// thread sends a Notify to wake up the next one.
        thread_control.Wait(0);
        --thread_control.count;
        int condcheckerror=0;
        if(mesh_id !=  state.mesh->Id()) {
          // Error?
          condcheckerror=1;
          Oxs_ThreadPrintf(stderr,"Invalid condition in"
                           " Oxs_ThermalField::ComputeEnergyChunk(),"
                           " thread number %d\n",threadnumber);
        }
        if(thread_control.count>0) {
          // Free a waiting thread.
          thread_control.Notify();
        }
        thread_control.Unlock();
        if(condcheckerror || Oxs_ThreadError::IsError()) {
          return; // What else?
        }
      } else {
        if(thread_control.count>0) {
          // Free a waiting thread.  (Actually, it can occur that the
          // thread_control will be grabbed by another thread that is
          // blocked at the first thread_control mutex Lock() call above
          // rather than on the ConditionWait, in which case this
          // ConditionNotify will be effectively lost.  But that is
          // okay, because then *that* thread will Notify when it
          // releases the mutex.)
          thread_control.Notify();
        }
        thread_control.Unlock();
      }
    } else {
      // Main thread (threadnumber == 0)
      try {
        T_init->FillMeshValue(state.mesh,T);
        alpha_init->FillMeshValue(state.mesh,alpha);
        const OC_INDEX size = state.mesh->Size();
        mesh_id = state.mesh->Id();
      } catch(Oxs_ExtError& err) {
        // Leave unmatched mesh_id as a flag to check
        // Oxs_ThreadError for an error.
        Oxs_ThreadError::SetError(String(err));
        if(thread_control.count>0) {
          thread_control.Notify();
        }
        thread_control.Unlock();
        throw;
      } catch(String& serr) {
        // Leave unmatched mesh_id as a flag to check
        // Oxs_ThreadError for an error.
        Oxs_ThreadError::SetError(serr);
        if(thread_control.count>0) {
          thread_control.Notify();
        }
        thread_control.Unlock();
        throw;
      } catch(const char* cerr) {
        // Leave unmatched mesh_id as a flag to check
        // Oxs_ThreadError for an error.
        Oxs_ThreadError::SetError(String(cerr));
        if(thread_control.count>0) {
          thread_control.Notify();
        }
        thread_control.Unlock();
        throw;
      } catch(...) {
        // Leave unmatched mesh_id as a flag to check
        // Oxs_ThreadError for an error.
        Oxs_ThreadError::SetError(String("Error in "
            "Oxs_ThermalField::ComputeEnergyChunk"));
        if(thread_control.count>0) {
          thread_control.Notify();
        }
        thread_control.Unlock();
        throw;
      }
      if(thread_control.count>0) {
        // Free a waiting thread.  (Actually, it can occur that the
        // thread_control will be grabbed by another thread that is
        // blocked at the first thread_control mutex Lock() call above
        // rather than on the ConditionWait, in which case this
        // ConditionNotify will be effectively lost.  But that is
        // okay, because then *that* thread will Notify when it
        // releases the mutex.)
        thread_control.Notify();
      }
      thread_control.Unlock();
    }
  }

    IntegEnergy(state,ocedt,ocedtaux,node_start,node_stop);
    return;
  }

 void Oxs_ThermalField::ComputeEnergyChunkInitialize(const Oxs_SimState& state,
					    Oxs_ComputeEnergyDataThreaded& ocedt,
					    Oxs_ComputeEnergyDataThreadedAux& ocedtaux,
					    int threadnumber) const 
{

  mesh_id = 0;
  T.Release();
  alpha.Release();
 // start the clock used for seeding
  beginning = myclock::now();
  dx = myclock::now() - beginning;
  unsigned seedx = dx.count();

  generatorx.seed(seedx);
  distributionx = std::normal_distribution<OC_REAL8m> (0.0,D);
  distributionx.reset();

}
