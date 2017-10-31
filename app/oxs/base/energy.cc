/* FILE: energy.cc                 -*-Mode: c++-*-
 *
 * Abstract energy class, derived from Oxs_Ext class.
 *
 */

#include <assert.h>
#include <string>

#include "energy.h"
#include "mesh.h"

OC_USE_STRING;

/* End includes */

#ifdef EXPORT_CALC_COUNT
void Oxs_Energy::FillCalcCountOutput(const Oxs_SimState& state)
{
  calc_count_output.cache.state_id = state.Id();
  calc_count_output.cache.value    = static_cast<OC_REAL8m>(calc_count);
}
#endif // EXPORT_CALC_COUNT



void Oxs_Energy::SetupOutputs()
{
  energy_sum_output.Setup(this,InstanceName(),"Energy","J",1,
			  &Oxs_Energy::UpdateStandardOutputs);
  field_output.Setup(this,InstanceName(),"Field","A/m",1,
		     &Oxs_Energy::UpdateStandardOutputs);
  energy_density_output.Setup(this,InstanceName(),"Energy density","J/m^3",1,
		     &Oxs_Energy::UpdateStandardOutputs);
#ifdef EXPORT_CALC_COUNT
  calc_count_output.Setup(this,InstanceName(),"Calc count","",0,
			  &Oxs_Energy::FillCalcCountOutput);
#endif // EXPORT_CALC_COUNT
  // Note: MS VC++ 6.0 requires fully qualified member names

  // Question: Do we want to add a mxH output?

  // Register outputs
  energy_sum_output.Register(director,0);
  field_output.Register(director,0);
  energy_density_output.Register(director,0);
#ifdef EXPORT_CALC_COUNT
  calc_count_output.Register(director,0);
#endif // EXPORT_CALC_COUNT

  // Eventually, caching should be handled by controlling Tcl script.
  // Until then, request caching of scalar energy output by default.
  energy_sum_output.CacheRequestIncrement(1);
}

// Constructors
Oxs_Energy::Oxs_Energy
( const char* name,     // Child instance id
  Oxs_Director* newdtr  // App director
  ) : Oxs_Ext(name,newdtr),calc_count(0)
{ SetupOutputs(); }

Oxs_Energy::Oxs_Energy
( const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr    // MIF block argument string
  ) : Oxs_Ext(name,newdtr,argstr),calc_count(0)
{ SetupOutputs(); }

//Destructor
Oxs_Energy::~Oxs_Energy() {
#if REPORT_TIME
  Oc_TimeVal cpu,wall;
  energytime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"GetEnergy time (secs)%7.2f cpu /%7.2f wall,"
            " module %s (%u evals)\n",double(cpu),double(wall),
            InstanceName(),GetEnergyEvalCount());
  }
#endif // REPORT_TIME
}

// Default problem initializer routine.  Any child that redefines
// this function should embed a call to this Init() inside
// the child specific version.
OC_BOOL Oxs_Energy::Init()
{
  if(!Oxs_Ext::Init()) return 0;

#if REPORT_TIME
  Oc_TimeVal cpu,wall;
  energytime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"GetEnergy time (secs)%7.2f cpu /%7.2f wall,"
            " module %s (%u evals)\n",double(cpu),double(wall),
            InstanceName(),GetEnergyEvalCount());
  }
  energytime.Reset();
#endif // REPORT_TIME

  calc_count=0;
  return 1;
}


// Standard output object update interface
void Oxs_Energy::UpdateStandardOutputs(const Oxs_SimState& state)
{
  if(state.Id()==0) { // Safety
    return;
  }

#if REPORT_TIME
  energytime.Start();
#endif // REPORT_TIME

  Oxs_ComputeEnergyData oced(state);

  // Dummy buffer space.
  Oxs_MeshValue<OC_REAL8m> dummy_energy;
  Oxs_MeshValue<ThreeVector> dummy_field;
  oced.scratch_energy = &dummy_energy;
  oced.scratch_H      = &dummy_field;

  if(energy_density_output.GetCacheRequestCount()>0) {
    energy_density_output.cache.state_id=0;
    oced.scratch_energy = oced.energy = &energy_density_output.cache.value;
    oced.energy->AdjustSize(state.mesh);
  }

  if(field_output.GetCacheRequestCount()>0) {
    field_output.cache.state_id=0;
    oced.scratch_H = oced.H = &field_output.cache.value;
    oced.H->AdjustSize(state.mesh);
  }

  ++calc_count;
  ComputeEnergy(state,oced);

  if(energy_density_output.GetCacheRequestCount()>0) {
    energy_density_output.cache.state_id=state.Id();
  }
  if(field_output.GetCacheRequestCount()>0) {
    field_output.cache.state_id=state.Id();
  }
  if(energy_sum_output.GetCacheRequestCount()>0) {
    energy_sum_output.cache.value=oced.energy_sum;
    energy_sum_output.cache.state_id=state.Id();
  }

#if REPORT_TIME
  energytime.Stop();
#endif // REPORT_TIME

}

///////////////////////////////////////////////////////////////////
// Energy and energy derivatives calculation function.  This
// should be "conceptually const."  The energy return is average
// energy density for the corresponding cell, in J/m^3.  The
// field is in A/m, pE_pt is in J/s.
//   There are two GetEnergy functions.  The first is a
// private, virtual member function of Oxs_Energy.  It takes as
// imports Oxs_MeshValue references into which to store the
// energy and H results.  The second is a non-virtual public
// member function, GetEnergyData, that passes in Oxs_MeshValue
// references that *may* be used to store the results, and
// separate parameters that return pointers to where the
// results actually were stored.  If output caching in enabled,
// the Oxs_Energy base class calls the virtual GetEnergy
// function with references to the output cache Oxs_MeshValue
// objects.  Otherwise the import buffer space is used.
//   NOTE: The pE_pt export is the partial derivative of energy
// with respect to time.  For most energy terms this will be
// 0.  It will only be non-zero if there is an explicit dependence
// on time, as for example with a time-varying applied field.
void Oxs_Energy::GetEnergyData
(const Oxs_SimState& state,
 Oxs_EnergyData& oed
 )
{
  if(oed.energy_buffer==NULL || oed.field_buffer==NULL) {
    // Bad input
    String msg = String("Oxs_EnergyData object in function"
			" Oxs_Energy::GetEnergyData"
			" contains NULL buffer pointers.");
    throw Oxs_ExtError(this,msg.c_str());
  }

  if(state.Id()==0) {
    String msg = String("Programming error:"
			" Invalid (unlocked) state detected"
			" in Oxs_Energy::GetEnergyData.");
    throw Oxs_ExtError(this,msg.c_str());
  }

#if REPORT_TIME
  energytime.Start();
#endif // REPORT_TIME

  if(field_output.GetCacheRequestCount()>0) {
    field_output.cache.state_id=0;
    oed.field_buffer = &field_output.cache.value;
  }

  if(energy_density_output.GetCacheRequestCount()>0) {
    energy_density_output.cache.state_id=0;
    oed.energy_buffer = &energy_density_output.cache.value;
  }

  oed.energy_buffer->AdjustSize(state.mesh);
  oed.field_buffer->AdjustSize(state.mesh);

  ++calc_count;
  GetEnergy(state,oed);

  if(field_output.GetCacheRequestCount()>0) {
    field_output.cache.state_id=state.Id();
  }
  if(energy_density_output.GetCacheRequestCount()>0) {
    energy_density_output.cache.state_id=state.Id();
  }
  if(energy_sum_output.GetCacheRequestCount()>0) {
    if(oed.energy_sum.IsSet()) {
      energy_sum_output.cache.value=oed.energy_sum;
    } else {
      OC_INDEX size = state.mesh->Size();
      Nb_Xpfloat energy_sum=0;
      const Oxs_MeshValue<OC_REAL8m>& ebuf = *oed.energy;
      const Oxs_Mesh& mesh = *state.mesh;
      for(OC_INDEX i=0; i<size; ++i) {
	energy_sum += ebuf[i] * mesh.Volume(i);
      }
      energy_sum_output.cache.value=energy_sum.GetValue();
      oed.energy_sum.Set(energy_sum.GetValue()); // Might as well set
      /// this field if we are doing the calculation anyway.
    }
    energy_sum_output.cache.state_id=state.Id();
  }
#if REPORT_TIME
  energytime.Stop();
#endif // REPORT_TIME
}


////////////////////////////////////////////////////////////////////////
// The ComputeEnergy interface replaces the older GetEnergy interface.
// The parameter list is similar, but ComputeEnergy uses the
// Oxs_ComputeEnergyData data struct in place Oxs_EnergyData.  The
// state_id, scratch_energy and scratch_H members of
// Oxs_ComputeEnergyData must be set on entry to ComputeEnergy.  The
// scratch_* members must be non-NULL, but the underlying
// Oxs_MeshValue objects will be size adjusted as (and if) needed.
// The scratch_* members are need for backward compatibility with
// older (pre Oct 2008) Oxs_Energy child classes, but also for
// Oxs_Energy classes like Oxs_Demag that always require space for
// field output.  Member "scratch_energy" is expressly allowed to be
// the same as member "energy", and likewise for "scratch_H" and "H".
//
// The remaining Oxs_MeshValue pointers are output requests.  They can
// be NULL, in which case the output is not requested, or non-NULL, in
// which case output is requested.  If output is requested, then the
// corresponding Oxs_MeshValue object will be filled.  (Note that the
// usual ComputeEnergy caller, AccumEnergyAndTorque, may adjust some
// of these pointers to point into Oxs_Energy standard output cache
// space, but the ComputeEnergy function itself plays no such games.)
// Any of these members that are non-NULL must be pre-sized
// appropriately for the given mesh.  This sizing is done automatically
// by AccumEnergyAndTorque for the "energy", "H", and "mxH" members,
// but not for the "accum" members.
//
// The remaining members, energy_sum and pE_pt are exports that are
// always filled by ComputeEnergy.
//
// The main practical advantage of ComputeEnergy over GetEnergy
// is that use of the "accum" fields can allow significant reduction
// in memory bandwidth use in evolvers.  This can be especially
// important in multi-threaded settings.
//
// The following is a default implementation of the virtual
// ComputeEnergy function.  It is effectively a wrapper/adapter
// to the deprecated GetEnergy function.  New Oxs_Energy child
// classes should override this function with a child-specific
// version, and define their GetEnergy function as a simple wrapper
// to GetEnergyAlt (q.v.).
//
void Oxs_Energy::ComputeEnergy
(const Oxs_SimState& state,
 Oxs_ComputeEnergyData& oced) const
{
  const Oxs_Mesh* mesh = state.mesh;

  if(oced.scratch_energy==NULL || oced.scratch_H==NULL) {
    // Bad input
    String msg = String("Oxs_ComputeEnergyData object in function"
			" Oxs_Energy::ComputeEnergy"
			" contains NULL scratch pointers.");
    throw Oxs_ExtError(this,msg.c_str());
  }

  if((oced.energy_accum!=0 && !oced.energy_accum->CheckMesh(mesh))
     || (oced.H_accum!=0   && !oced.H_accum->CheckMesh(mesh))
     || (oced.mxH_accum!=0 && !oced.mxH_accum->CheckMesh(mesh))
     || (oced.energy!=0    && !oced.energy->CheckMesh(mesh))
     || (oced.H!=0         && !oced.H->CheckMesh(mesh))
     || (oced.mxH!=0       && !oced.mxH->CheckMesh(mesh))) {
    // Bad input
    String msg = String("Oxs_ComputeEnergyData object in function"
			" Oxs_Energy::ComputeEnergy"
			" contains ill-sized buffers.");
    throw Oxs_ExtError(this,msg.c_str());
  }

  Oxs_EnergyData oed(state);
  if(oced.energy) oed.energy_buffer = oced.energy;
  else            oed.energy_buffer = oced.scratch_energy;
  if(oced.H)      oed.field_buffer  = oced.H;
  else            oed.field_buffer  = oced.scratch_H;

  // Although not stated in the interface docs, some Oxs_Energy children
  // assume that the oed energy and field buffers are pre-sized on
  // entry.  For backwards compatibility, make this so.
  oed.energy_buffer->AdjustSize(mesh);
  oed.field_buffer->AdjustSize(mesh);

  GetEnergy(state,oed);

  // Accum as requested
  OC_INDEX i;
  const OC_INDEX size = mesh->Size();

  const OC_BOOL have_energy_sum = oed.energy_sum.IsSet();
  if(have_energy_sum) oced.energy_sum = oed.energy_sum;
  else                oced.energy_sum = 0.0;

  if(oced.energy_accum) {
    Oxs_MeshValue<OC_REAL8m>& energy_accum = *oced.energy_accum;
    const Oxs_MeshValue<OC_REAL8m>& energy = *(oed.energy.Get());
    for(i=0; i<size; ++i) {
      energy_accum[i] += energy[i];
      if(!have_energy_sum) {
        oced.energy_sum += energy[i] * mesh->Volume(i);
      }
    }
  } else if(!have_energy_sum) {
    const Oxs_MeshValue<OC_REAL8m>& energy = *(oed.energy.Get());
    for(i=0; i<size; ++i) {
      oced.energy_sum += energy[i] * mesh->Volume(i);
    }
  }

  const Oxs_MeshValue<ThreeVector>& spin = state.spin;
  const Oxs_MeshValue<ThreeVector>& H = *(oed.field.Get());
  if(oced.mxH_accum && oced.H_accum) {
    for(i=0; i<size; ++i) {
      (*oced.H_accum)[i] += H[i];
      ThreeVector temp = spin[i];
      temp ^= H[i];
      (*oced.mxH_accum)[i] += temp;
    }
  } else if(oced.mxH_accum) {
    for(i=0; i<size; ++i) {
      ThreeVector temp = spin[i];
      temp ^= H[i];
      (*oced.mxH_accum)[i] += temp;
    }
  } else if(oced.H_accum) {
    (*oced.H_accum) += H;
  }

  // Copy energy and field results, as needed
  if(oced.energy && oced.energy != oed.energy.Get()) (*oced.energy) = (*oed.energy);
  if(oced.H      && oced.H      != oed.field.Get())  (*oced.H)      = (*oed.field);

  // mxH requested?
  if(oced.mxH) {
    for(i=0; i<size; ++i) {
      ThreeVector temp = spin[i];
      temp ^= H[i];
      (*oced.mxH)[i] = temp;
    }
  }

  // pE_pt
 if(oed.pE_pt.IsSet()) oced.pE_pt = oed.pE_pt;
 else                  oced.pE_pt = 0.0;
}

////////////////////////////////////////////////////////////////////////
// GetEnergyAlt is an adapter that allows a child-defined ComputeEnergy
// class to provide GetEnergy services.  Such children can define their
// GetEnergy to be a simple wrapper to GetEnergyAlt.  NOTE: Children
// must NOT wrap GetEnergyAlt with GetEnergy without also overriding the
// default ComputeEnergy implementation.  Otherwise an infinite circular
// call sequence will result:
//  GetEnergy -> GetEnergyAlt -> ComputeEnergy -> GetEnergy -> ...
//
void Oxs_Energy::GetEnergyAlt
(const Oxs_SimState& state,
 Oxs_EnergyData& oed) const
{
  if(oed.energy_buffer==NULL || oed.field_buffer==NULL) {
    // Bad input
    String msg = String("Oxs_EnergyData object in function"
			" Oxs_Energy::GetEnergyAlt"
			" contains NULL buffer pointers.");
    throw Oxs_ExtError(this,msg.c_str());
  }

  Oxs_ComputeEnergyData oced(state);

  oed.energy = oced.energy = oced.scratch_energy = oed.energy_buffer;
  oed.field  = oced.H      = oced.scratch_H      = oed.field_buffer;

  oced.energy->AdjustSize(state.mesh);
  oced.H->AdjustSize(state.mesh);

  ComputeEnergy(state,oced);

  oed.energy_sum = oced.energy_sum;
  oed.pE_pt      = oced.pE_pt;
}


////////////////////////////////////////////////////////////////////////
// The AccumEnergyAndTorque method is similar to GetEnergyData, except
// that the former uses the newer ComputeEnergy interface rather than
// the deprecated GetEnergy interface.
//
// The "oced.scratch_*" members must be filled before entry.  They
// point to scratch space that may or may not be used.  This space
// will be resized as needed.  If desired, "scratch_energy" may point
// to the same place as "energy", and likewise for "scratch_H" and
// "H".
//
// The remaining Oxs_MeshValue pointers are output requests.  If one
// of these pointers is NULL on entry, then that output is requested.
// The output will to into the pointed-to space, *unless* that output
// is associated with one of the Oxs_Energy standard outputs
// (energy_density_output or field_output) and caching of that output
// is enabled, in which case the output will go into the output cache
// and the corresponding oced pointer (energy or H) will be changed to
// point to the cache.  Pay ATTENTION to this point: the arrays sent
// in *may* not be the ones used, so clients should always check and
// use the oced pointers directly, rather than the arrays values sent
// in.  Also, if caching is enabled, then on return the energy and/or
// H pointers in oced will be set to the cache, even if energy and/or
// H were set to NULL on entry.
//
// The oced.*_accum members are accumulated (added) into rather than
// set.
//
void Oxs_Energy::AccumEnergyAndTorque
(const Oxs_SimState& state,
 Oxs_ComputeEnergyData& oced)
{
  if(state.Id()==0) {
    String msg = String("Programming error:"
			" Invalid (unlocked) state detected"
			" in Oxs_Energy::AccumEnergyAndTorque");
    throw Oxs_ExtError(this,msg.c_str());
  }

  if(oced.scratch_energy==NULL || oced.scratch_H==NULL) {
    // Bad input
    String msg = String("Oxs_ComputeEnergyData object in function"
			" Oxs_Energy::AccumEnergyAndTorque"
			" contains NULL scratch pointers.");
    throw Oxs_ExtError(this,msg.c_str());
  }

  if((oced.energy_accum!=0 && !oced.energy_accum->CheckMesh(state.mesh))
     || (oced.H_accum!=0   && !oced.H_accum->CheckMesh(state.mesh))
     || (oced.mxH_accum!=0 && !oced.mxH_accum->CheckMesh(state.mesh))) {
    // Bad input
    String msg = String("Oxs_ComputeEnergyData object in function"
			" Oxs_Energy::AccumEnergyAndTorque"
			" contains ill-sized accum buffers.");
    throw Oxs_ExtError(this,msg.c_str());
  }

#if REPORT_TIME
  energytime.Start();
#endif // REPORT_TIME

  if(energy_density_output.GetCacheRequestCount()>0) {
    energy_density_output.cache.state_id=0;
    oced.energy = &energy_density_output.cache.value;
  }

  if(field_output.GetCacheRequestCount()>0) {
    field_output.cache.state_id=0;
    oced.H = &field_output.cache.value;
  }

  if(oced.energy) oced.energy->AdjustSize(state.mesh);
  if(oced.H)      oced.H->AdjustSize(state.mesh);
  if(oced.mxH)    oced.mxH->AdjustSize(state.mesh);

  ++calc_count;
  ComputeEnergy(state,oced);

  if(field_output.GetCacheRequestCount()>0) {
    field_output.cache.state_id=state.Id();
  }
  if(energy_density_output.GetCacheRequestCount()>0) {
    energy_density_output.cache.state_id=state.Id();
  }
  if(energy_sum_output.GetCacheRequestCount()>0) {
    energy_sum_output.cache.value=oced.energy_sum;
    energy_sum_output.cache.state_id=state.Id();
  }
#if REPORT_TIME
  energytime.Stop();
#endif // REPORT_TIME
}

///////////////////////// CHUNK ENERGY /////////////////////////////////

struct Oxs_ComputeEnergies_ChunkStruct {
public:
  Oxs_ChunkEnergy* energy;
  Oxs_ComputeEnergyDataThreaded ocedt;
  Oxs_ComputeEnergyDataThreadedAux ocedtaux;
  Oxs_ComputeEnergies_ChunkStruct()
    : energy(0) {}
};

class Oxs_ComputeEnergiesChunkThread : public Oxs_ThreadRunObj {
public:
  static Oxs_JobControl<ThreeVector> job_basket;
  /// job_basket is static, so only one "set" of this class is allowed.

  const Oxs_SimState* state;
  vector<Oxs_ComputeEnergies_ChunkStruct> energy_terms;

  Oxs_MeshValue<ThreeVector>* mxH;
  Oxs_MeshValue<ThreeVector>* mxH_accum;
  Oxs_MeshValue<ThreeVector>* mxHxm;
  const vector<OC_INDEX>* fixed_spins;
  OC_REAL8m max_mxH;

  OC_INDEX cache_blocksize;

  OC_BOOL accums_initialized;

  Oxs_ComputeEnergiesChunkThread()
    : state(0),
      mxH(0),mxH_accum(0),
      mxHxm(0), fixed_spins(0),
      max_mxH(0.0),
      cache_blocksize(0), accums_initialized(0) {}

  void Cmd(int threadnumber, void* data);

  static void Init(int thread_count,
                   const Oxs_StripedArray<ThreeVector>* arrblock) {
    job_basket.Init(thread_count,arrblock);
  }

  // Note: Default copy constructor and assignment operator,
  // and destructor.
};

Oxs_JobControl<ThreeVector> Oxs_ComputeEnergiesChunkThread::job_basket;

void
Oxs_ComputeEnergiesChunkThread::Cmd
(int threadnumber,
 void* /* data */)
{
  OC_REAL8m max_mxH_sq = 0.0;
  const Oxs_MeshValue<ThreeVector>& spin = state->spin;
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state->Ms);

  // In chunk post-processing segment, the torque at fixed spins
  // is forced to zero.  Variable i_fixed holds the latest working
  // position in the fixed_spins array between chunks.
  OC_INDEX i_fixed = 0;
  OC_INDEX i_fixed_total = 0;
  if(fixed_spins) i_fixed_total = fixed_spins->size();

  while(1) {
    // Claim a chunk
    OC_INDEX index_start,index_stop;
    job_basket.GetJob(threadnumber,index_start,index_stop);

    if(index_start>=index_stop) break;

    // We claim by blocksize, but work by cache_blocksize (which is
    // presumably smaller).  We want blocksize big enough to reduce
    // mutex collisions and other overhead, but cache_blocksize small
    // enough that the spin data can reside in cache acrosss all the
    // energy terms.

    for(OC_INDEX icache_start=index_start;
        icache_start<index_stop; icache_start+=cache_blocksize) {
      OC_INDEX icache_stop = icache_start + cache_blocksize;
      if(icache_stop>index_stop) icache_stop = index_stop;

      // Process chunk
      OC_UINT4m energy_item = 0;
      for(vector<Oxs_ComputeEnergies_ChunkStruct>::iterator eit
            = energy_terms.begin();
          eit != energy_terms.end() ; ++eit, ++energy_item) {

        // Set up some refs for convenience
        Oxs_ChunkEnergy& eterm = *(eit->energy);
        Oxs_ComputeEnergyDataThreaded& ocedt = eit->ocedt;
        Oxs_ComputeEnergyDataThreadedAux& ocedtaux = eit->ocedtaux;
#if REPORT_TIME
# if 0  // Individual chunk times currently meaningless,
        //and may slow code due to mutex blocks.
        ocedtaux.energytime.Start();
# endif 
#endif // REPORT_TIME
        if(!accums_initialized && energy_item==0) {
          // Note: Each thread has its own copy of the ocedt and
          // ocedtaux data, so we can tweak these as desired without
          // stepping on other threads

          // Move each accum pointer to corresponding non-accum member
          // for initialization.
          assert(ocedt.mxH == 0);
          Oxs_MeshValue<OC_REAL8m>* energy_accum_save = ocedt.energy_accum;
          if(ocedt.energy == 0) ocedt.energy = ocedt.energy_accum;
          ocedt.energy_accum = 0;

          Oxs_MeshValue<ThreeVector>* H_accum_save = ocedt.H_accum;
          if(ocedt.H == 0)      ocedt.H      = ocedt.H_accum;
          ocedt.H_accum      = 0;

          // Note: ocedt.mxH should always be zero, but check here anyway
          // for easier code maintenance.
          Oxs_MeshValue<ThreeVector>* mxH_accum_save = ocedt.mxH_accum;
          if(ocedt.mxH == 0)    ocedt.mxH    = ocedt.mxH_accum;
          ocedt.mxH_accum    = 0;

          eterm.ComputeEnergyChunk(*state,ocedt,ocedtaux,
                                   icache_start,icache_stop,
                                   threadnumber);

          // Copy data as necessary
          if(energy_accum_save) {
            if(ocedt.energy != energy_accum_save) {
              for(OC_INDEX i=icache_start;i<icache_stop;++i) {
                (*energy_accum_save)[i] = (*(ocedt.energy))[i];
              }
            } else {
              ocedt.energy = 0;
            }
          }
          ocedt.energy_accum = energy_accum_save;

          if(H_accum_save) {
            if(ocedt.H != H_accum_save) {
              for(OC_INDEX i=icache_start;i<icache_stop;++i) {
                (*H_accum_save)[i] = (*(ocedt.H))[i];
              }
            } else {
              ocedt.H = 0;
            }
          }
          ocedt.H_accum = H_accum_save;

          if(mxH_accum_save) {
            if(ocedt.mxH != mxH_accum_save) {
              // This branch should never run
              abort();
            } else {
              ocedt.mxH = 0;
            }
          }
          ocedt.mxH_accum = mxH_accum_save;

        } else {
          // Standard processing: accum elements already initialized.
          eterm.ComputeEnergyChunk(*state,ocedt,ocedtaux,
                                   icache_start,icache_stop,
                                   threadnumber);
        }

#if REPORT_TIME
# if 0  // Individual chunk times currently meaningless,
        //and may slow code due to mutex blocks.
        ocedtaux.energytime.Stop();
# endif
#endif // REPORT_TIME
      }

      // Post-processing, for this energy term and chunk.

      // Zero torque on fixed spins.  This code assumes that, 1) the
      // fixed_spins list is sorted in increasing order, and 2) the
      // chunk indices come in strictly monotonically increasing
      // order.
      // NB: Outside this loop, "i_fixed" stores the search start
      // location for the next chunk.
      while(i_fixed < i_fixed_total) {
        OC_INDEX index = (*fixed_spins)[i_fixed];
        if(index <  icache_start) { ++i_fixed; continue; }
        if(index >= icache_stop) break;
        if(mxH)       (*mxH)[index].Set(0.,0.,0.);
        if(mxH_accum) (*mxH_accum)[index].Set(0.,0.,0.);
        ++i_fixed;
      }

      // Note: The caller must pre-size mxHxm as appropriate.
      if(mxHxm) {
        // Compute mxHxm and max_mxH
        if(mxH_accum == mxHxm) {
          for(OC_INDEX i=icache_start;i<icache_stop;++i) {
            if(Ms[i]==0.0) { // Ignore zero-moment spins
              (*mxHxm)[i].Set(0.,0.,0.);
              continue;
            }
            OC_REAL8m tx = (*mxHxm)[i].x;
            OC_REAL8m ty = (*mxHxm)[i].y;
            OC_REAL8m tz = (*mxHxm)[i].z;
            OC_REAL8m mx = spin[i].x;  
            OC_REAL8m my = spin[i].y;  
            OC_REAL8m mz = spin[i].z;  

            OC_REAL8m magsq = tx*tx + ty*ty + tz*tz;
            if(magsq > max_mxH_sq) max_mxH_sq = magsq;

            (*mxHxm)[i].x = ty*mz - tz*my;
            (*mxHxm)[i].y = tz*mx - tx*mz;
            (*mxHxm)[i].z = tx*my - ty*mx;
          }
        } else {
          // It is the responsibility of the caller to make certain
          // that if mxHxm is non-zero, then so is mxH_accum.
          assert(mxH_accum != 0);
          for(OC_INDEX i=icache_start;i<icache_stop;++i) {
            if(Ms[i]==0.0) { // Ignore zero-moment spins
              (*mxH_accum)[i].Set(0.,0.,0.);
              (*mxHxm)[i].Set(0.,0.,0.);
              continue;
            }
            OC_REAL8m tx = (*mxH_accum)[i].x;
            OC_REAL8m ty = (*mxH_accum)[i].y;
            OC_REAL8m tz = (*mxH_accum)[i].z;
            OC_REAL8m mx = spin[i].x;
            OC_REAL8m my = spin[i].y;
            OC_REAL8m mz = spin[i].z;

            OC_REAL8m magsq = tx*tx + ty*ty + tz*tz;
            if(magsq > max_mxH_sq) max_mxH_sq = magsq;

            (*mxHxm)[i].x = ty*mz - tz*my;
            (*mxHxm)[i].y = tz*mx - tx*mz;
            (*mxHxm)[i].z = tx*my - ty*mx;
          }
        }
      } else if(mxH_accum) {
          for(OC_INDEX i=icache_start;i<icache_stop;++i) {
            if(Ms[i]==0.0) { // Ignore zero-moment spins
              (*mxH_accum)[i].Set(0.,0.,0.);
              continue;
            }
            OC_REAL8m tx = (*mxH_accum)[i].x;
            OC_REAL8m ty = (*mxH_accum)[i].y;
            OC_REAL8m tz = (*mxH_accum)[i].z;
            OC_REAL8m magsq = tx*tx + ty*ty + tz*tz;
            if(magsq > max_mxH_sq) max_mxH_sq = magsq;
          }
      }
      // Otherwise, don't compute max_mxH
    }
  }

  max_mxH = sqrt(max_mxH_sq);
}

void Oxs_ComputeEnergies
(const Oxs_SimState& state,
 Oxs_ComputeEnergyData& oced,
 const vector<Oxs_Energy*>& energies,
 Oxs_ComputeEnergyExtraData& oceed)
{ // Compute sums of energies, fields, and/or torques for all energies
  // in "energies" import.  On entry, oced.energy_accum, oced.H_accum,
  // and oced.mxH_accum should be set or null as desired.
  // oced.scratch_energy and oced.scratch_H must be non-null.
  // oced.energy, oced.H and oced.mxH *must* be *null* on entry.  This
  // routine does not fill these fields, but rather the accumulated
  // values are collected as necessary in oced.*_accum entries.
  // (However, pointers to energy_accum, H_accum, and mxH_accum may
  // be temporarily swapped to energy, H, and mxH for initialization
  // purposes.  This is transparent to the Oxs_ComputeEnergies caller,
  // but will exercise the non-accum portions of callees.)
  //   This routine handles outputs, energy calculation counts, and
  // timers appropriately.
  //   Those "energies" members that are actually Oxs_ChunkEnergies will
  // use the ComputeEnergyChunk interface in a collated fashion to help
  // minimize memory bandwidth usage.  On threaded OOMMF builds, these
  // calls will be run in parallel and load balanced.  Also, the number
  // of threads launched will not exceed the number of chunks.  This is
  // to insure that the main thread (threadnumber == 0) has an
  // opportunity to run for initialization purposes in the
  // Oxs_ChunkEnergy::ComputeEnergyChunk() function.  (Oxs_ChunkEnergy
  // classes that make (or may make) call into the Tcl interpreter must
  // use threadnumber == 0 for those calls, as per Tcl specs.  So if
  // all threads with threadnumber != 0 block on ComputeEnergyChunk()
  // entry, then the threadnumber == 0 is guaranteed at least one call
  // into ComputeEnergyChunk().
  //
  // Update May-2009: The now preferred initialization method is to
  // use ComputeEnergyChunkInitialize.  The guarantee that threadnumber
  // 0 will always run is honored for backward compatibility, but new
  // code should use ComputeEnergyChunkInitialize instead.
  //
  //    Data in Oxs_ComputeEnergyExtraData are filled in on the backside
  // of the chunk compute code.  These results could be computed by the
  // client, but doing it here gives improved cache locality.

  if(state.Id()==0) {
    String msg = String("Programming error:"
			" Invalid (unlocked) state detected"
			" in Oxs_ComputeEnergies");
    throw Oxs_ExtError(msg);
  }

  if(oced.scratch_energy==NULL || oced.scratch_H==NULL) {
    // Bad input
    String msg = String("Oxs_ComputeEnergyData object in function"
			" Oxs_ComputeEnergies"
			" contains NULL scratch pointers.");
    throw Oxs_ExtError(msg);
  }

  if(oced.energy != NULL || oced.H != NULL || oced.mxH != NULL) {
    String msg = String("Programming error in function"
			" Oxs_ComputeEnergies:"
			" non-NULL energy, H, and/or mxH imports.");
    throw Oxs_ExtError(msg);
  }

  const int thread_count = Oc_GetMaxThreadCount();

  if(oced.energy_accum) {
    oced.energy_accum->AdjustSize(state.mesh);
  }
  if(oced.H_accum) {
    oced.H_accum->AdjustSize(state.mesh);
  }
  if(oced.mxH_accum) {
    oced.mxH_accum->AdjustSize(state.mesh);
  }
  if(oceed.mxHxm) {
    oceed.mxHxm->AdjustSize(state.mesh);
  }
  oced.energy_sum = 0.0;
  oced.pE_pt = 0.0;
  oceed.max_mxH = 0.0;

  if(energies.size() == 0) {
    // No energies.  Zero requested outputs and return.
    OC_INDEX size = state.mesh->Size();
    if(oced.energy_accum) {
      for(OC_INDEX i=0; i<size; ++i) {
        (*(oced.energy_accum))[i] = 0.0;
      }
    }
    if(oced.H_accum) {
      for(OC_INDEX i=0; i<size; ++i) {
        (*(oced.H_accum))[i] = ThreeVector(0.0,0.0,0.0);
      }
    }
    if(oced.mxH_accum) {
      for(OC_INDEX i=0; i<size; ++i) {
        (*(oced.mxH_accum))[i] = ThreeVector(0.0,0.0,0.0);
      }
    }
    if(oceed.mxHxm) {
      for(OC_INDEX i=0; i<size; ++i) {
        (*(oceed.mxHxm))[i] = ThreeVector(0.0,0.0,0.0);
      }
    }
    return;
  }

  if(oced.mxH_accum==0 && oceed.mxHxm!=0) {
    // Hack mxHxm into mxH_accum.  We can identify this situation
    // by checking mxH_accum == mxHxm, and undo at the end.  Also
    // The Oxs_ComputeEnergiesChunkThread objects know about this
    // and respond appropriately.
    oced.mxH_accum = oceed.mxHxm;
  }

  vector<Oxs_ComputeEnergies_ChunkStruct> chunk;
  vector<Oxs_Energy*> nonchunk;

  // Initialize those parts of ChunkStruct that are independent
  // of any particular energy term.
  Oxs_ComputeEnergies_ChunkStruct foo;
  foo.ocedt.state_id = state.Id();
  foo.ocedt.scratch_energy = oced.scratch_energy;
  foo.ocedt.scratch_H      = oced.scratch_H;
  foo.ocedt.energy_accum   = oced.energy_accum;
  foo.ocedt.H_accum        = oced.H_accum;
  foo.ocedt.mxH_accum      = oced.mxH_accum;
  for(vector<Oxs_Energy*>::const_iterator it = energies.begin();
      it != energies.end() ; ++it ) {
    Oxs_ChunkEnergy* ceptr =
      dynamic_cast<Oxs_ChunkEnergy*>(*it);
    if(ceptr != NULL) {
      // Set up and initialize chunk energy structures
      foo.energy = ceptr;
      if(ceptr->energy_density_output.GetCacheRequestCount()>0) {
        ceptr->energy_density_output.cache.state_id=0;
        foo.ocedt.energy = &(ceptr->energy_density_output.cache.value);
        foo.ocedt.energy->AdjustSize(state.mesh);
      }
      if(ceptr->field_output.GetCacheRequestCount()>0) {
        ceptr->field_output.cache.state_id=0;
        foo.ocedt.H = &(ceptr->field_output.cache.value);
        foo.ocedt.H->AdjustSize(state.mesh);
      }
      chunk.push_back(foo);
    } else {
      nonchunk.push_back(*it);
    }
  }

  // The "accum" elements are initialized on the first pass by
  // moving each accum pointer to the corresponding non-accum member.
  // After filling by the first energy term, the pointers are moved
  // back to the accum member.  This way we avoid a pass through
  // memory storing zeros, and a pass through memory loading zeros.
  // Zero load/stores are cheap in the chunk memory case, because
  // in that case the load/stores are just to and from cache, but
  // we prefer here to run non-chunk energies first so that we
  // can compute mxHxm and max |mxH| on the backsize of the chunk
  // energy runs (where m and mxH are in cache and so don't have
  // to be loaded).  Create a boolean to track initialization.
  OC_BOOL accums_initialized = 0;

  // Non-chunk energies //////////////////////////////////////
  for(vector<Oxs_Energy*>::const_iterator ncit = nonchunk.begin();
      ncit != nonchunk.end() ; ++ncit ) {
    Oxs_Energy& eterm = *(*ncit);  // Convenience

#if REPORT_TIME
    eterm.energytime.Start();
#endif // REPORT_TIME

    Oxs_ComputeEnergyData term_oced(state);
    term_oced.scratch_energy = oced.scratch_energy;
    term_oced.scratch_H      = oced.scratch_H;
    term_oced.energy_accum = oced.energy_accum;
    term_oced.H_accum      = oced.H_accum;
    term_oced.mxH_accum    = oced.mxH_accum;

    if(eterm.energy_density_output.GetCacheRequestCount()>0) {
      eterm.energy_density_output.cache.state_id=0;
      term_oced.energy = &(eterm.energy_density_output.cache.value);
      term_oced.energy->AdjustSize(state.mesh);
    }

    if(eterm.field_output.GetCacheRequestCount()>0) {
      eterm.field_output.cache.state_id=0;
      term_oced.H = &(eterm.field_output.cache.value);
      term_oced.H->AdjustSize(state.mesh);
    }

    if(!accums_initialized) {
      // Initialize by filling
      term_oced.energy_accum = 0;
      term_oced.H_accum = 0;
      term_oced.mxH_accum = 0;
      if(term_oced.energy == 0) term_oced.energy = oced.energy_accum;
      if(term_oced.H == 0)      term_oced.H      = oced.H_accum;
      if(term_oced.mxH == 0)    term_oced.mxH    = oced.mxH_accum;
    }

    ++(eterm.calc_count);
    eterm.ComputeEnergy(state,term_oced);

    if(eterm.field_output.GetCacheRequestCount()>0) {
      eterm.field_output.cache.state_id=state.Id();
    }
    if(eterm.energy_density_output.GetCacheRequestCount()>0) {
      eterm.energy_density_output.cache.state_id=state.Id();
    }
    if(eterm.energy_sum_output.GetCacheRequestCount()>0) {
      eterm.energy_sum_output.cache.value=term_oced.energy_sum;
      eterm.energy_sum_output.cache.state_id=state.Id();
    }

    if(!accums_initialized) {
      // If output buffer spaced was used instead of accum space, then
      // copy from output buffer to accum space.  This hurts from a
      // memory bandwidth perspective, but is rather hard to avoid.
      // (Options: Do accum initialization in chunk-energy branch,
      // but that hurts with respect to mxHxm and max |mxH| computations.
      // Or one could have the ComputeEnergy class fill more than one
      // array with the non-accum output (say, via a parameter that
      // says to set to accum rather than add to accum), but that is
      // rather awkward.  Instead, we assume that if the user wants
      // high speed then he won't enable term energy or H outputs.)
      if(oced.energy_accum && term_oced.energy != oced.energy_accum) {
        *(oced.energy_accum) = *(term_oced.energy);
      }
      if(oced.H_accum      && term_oced.H      != oced.H_accum) {
        *(oced.H_accum) = *(term_oced.H);
      }
      if(oced.mxH_accum    && term_oced.mxH    != oced.mxH_accum) {
        *(oced.mxH_accum) = *(term_oced.mxH);
      }
      accums_initialized = 1;
    }

    oced.energy_sum += term_oced.energy_sum;
    oced.pE_pt += term_oced.pE_pt;

#if REPORT_TIME
    eterm.energytime.Stop();
#endif // REPORT_TIME
  }


  // Chunk energies ///////////////////////////////////////////
#if REPORT_TIME
  Oxs_ChunkEnergy::chunktime.Start();
#endif

  // Compute cache_blocksize
  const OC_INDEX meshsize = state.mesh->Size();
  const OC_INDEX cache_size = 1024 * 1024;  // Should come from
  /// platform file or perhaps sysconf().

  const OC_INDEX recsize = sizeof(ThreeVector) + sizeof(OC_REAL8m);
  /// May want to query individual energies for this.

#define FUDGE 8
  OC_INDEX tcblocksize = (cache_size>FUDGE*recsize ?
			cache_size/(FUDGE*recsize) : 1);
  if(thread_count*tcblocksize>meshsize) {
    tcblocksize = meshsize/thread_count;
  }
  if(0 == tcblocksize) {
    tcblocksize = 1;    // Safety
  } else if(0 != tcblocksize%16) {
    tcblocksize += 16 - (tcblocksize%16);  // Make multiple of 16
  }
  const OC_INDEX cache_blocksize = tcblocksize;

  // Thread control
  static Oxs_ThreadTree threadtree;

  Oxs_ComputeEnergiesChunkThread::Init(thread_count,
                                       state.spin.GetArrayBlock());

  vector<Oxs_ComputeEnergiesChunkThread> chunk_thread;
  chunk_thread.resize(thread_count);
  chunk_thread[0].state     = &state;
  chunk_thread[0].energy_terms = chunk; // Make copies.
  chunk_thread[0].mxH       = oced.mxH;
  chunk_thread[0].mxH_accum = oced.mxH_accum;
  chunk_thread[0].mxHxm     = oceed.mxHxm;
  chunk_thread[0].fixed_spins = oceed.fixed_spin_list;
  chunk_thread[0].cache_blocksize = cache_blocksize;
  chunk_thread[0].accums_initialized = accums_initialized;

  // Initialize chunk energy computations
  for(vector<Oxs_ComputeEnergies_ChunkStruct>::iterator it
        = chunk.begin(); it != chunk.end() ; ++it ) {
    Oxs_ChunkEnergy& eterm = *(it->energy);  // For code clarity
    Oxs_ComputeEnergyDataThreaded& ocedt = it->ocedt;
    Oxs_ComputeEnergyDataThreadedAux& ocedtaux = it->ocedtaux;
    eterm.ComputeEnergyChunkInitialize(state,ocedt,ocedtaux,
                                       thread_count);
  }

  for(int ithread=1;ithread<thread_count;++ithread) {
    chunk_thread[ithread] = chunk_thread[0];
    threadtree.Launch(chunk_thread[ithread],0);
  }
  threadtree.LaunchRoot(chunk_thread[0],0);

  // Note: If chunk.size()>0, then we are guaranteed that accums are
  // initialized.  If accums_initialized is ever needed someplace
  // downstream, then uncomment the following line:
  // if(chunk.size()>0) accums_initialized = 1;

  // Finalize chunk energy computations
  for(OC_INDEX ei=0;static_cast<size_t>(ei)<chunk.size();++ei) {

    Oxs_ChunkEnergy& eterm = *(chunk[ei].energy);  // Convenience
    const Oxs_ComputeEnergyDataThreaded& ocedt = chunk[ei].ocedt;
    const Oxs_ComputeEnergyDataThreadedAux& ocedtaux = chunk[ei].ocedtaux;

    eterm.ComputeEnergyChunkFinalize(state,ocedt,ocedtaux,
                                     thread_count);

    ++(eterm.calc_count);

    // For each energy term, loop though all threads and sum
    // energy and pE_pt contributions.
    OC_REAL8m pE_pt_term = chunk[ei].ocedtaux.pE_pt_accum;
    for(int ithread=0;ithread<thread_count;++ithread) {
      pE_pt_term
        += chunk_thread[ithread].energy_terms[ei].ocedtaux.pE_pt_accum;
    }
    oced.pE_pt += pE_pt_term;

    OC_REAL8m energy_term = chunk[ei].ocedtaux.energy_total_accum;
    for(int ithread=0;ithread<thread_count;++ithread) {
      energy_term
        += chunk_thread[ithread].energy_terms[ei].ocedtaux.energy_total_accum;
    }
    oced.energy_sum += energy_term;

    if(eterm.energy_sum_output.GetCacheRequestCount()>0) {
      eterm.energy_sum_output.cache.value=energy_term;
      eterm.energy_sum_output.cache.state_id=state.Id();
    }

    if(eterm.field_output.GetCacheRequestCount()>0) {
      eterm.field_output.cache.state_id=state.Id();
    }

    if(eterm.energy_density_output.GetCacheRequestCount()>0) {
      eterm.energy_density_output.cache.state_id=state.Id();
    }

#if REPORT_TIME
    Nb_StopWatch bar;
    bar.ThreadAccum(chunk[ei].ocedtaux.energytime);
    for(int ithread=0;ithread<thread_count;++ithread) {
      bar.ThreadAccum
        (chunk_thread[ithread].energy_terms[ei].ocedtaux.energytime);

    }
    eterm.energytime.Accum(bar);
#endif // REPORT_TIME
  }

  if(oceed.mxHxm!=0 && oced.mxH_accum == oceed.mxHxm) {
    // Undo mxHxm hack
    oced.mxH_accum = 0;
  }

  oceed.max_mxH = 0.0;
  for(vector<Oxs_ComputeEnergiesChunkThread>::const_iterator cect
        = chunk_thread.begin(); cect != chunk_thread.end() ; ++cect ) {
    if(cect->max_mxH > oceed.max_mxH) oceed.max_mxH = cect->max_mxH;
  }

#if REPORT_TIME
  Oxs_ChunkEnergy::chunktime.Stop();
#endif

}

#if REPORT_TIME
Nb_StopWatch Oxs_ChunkEnergy::chunktime;

void Oxs_ChunkEnergy::ReportTime()
{
  Oc_TimeVal cpu,wall;
  chunktime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"GetEnergy time (secs)%7.2f cpu /%7.2f wall,"
            " ChunkEnergies total (%u evals)\n",
            double(cpu),double(wall),GetEnergyEvalCount());
    chunktime.Reset();  // Only print once (per run).
  }
}
#endif // REPORT_TIME
