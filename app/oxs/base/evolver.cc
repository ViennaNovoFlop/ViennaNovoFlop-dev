/* FILE: evolver.cc                 -*-Mode: c++-*-
 *
 * Abstract evolver class, derived from Oxs_Ext class.
 *
 */

#include <algorithm>   // "sort" algorithm used on fixed_spin_list

#include "director.h"
#include "evolver.h"

/* End includes */

// Fixed spin helper functions
void Oxs_Evolver::InitFixedRegions(const vector<String> &names)
{ // The first element of the names vector should be the name
  // of the atlas to use; the remaining fields are regions in
  // which the spins are fixed.  This function should be called
  // from the Oxs_Evolver (parent or child) constructor.

  fixed_region_ids.clear();  // Empty fixed regions list

  OC_UINT4m name_count = static_cast<OC_UINT4m>(names.size());
  if(name_count<1) { // Empty name list
    atlas_key.Release();
    return;
  }

  vector<String> atlas_list;
  Oxs_SplitList atlas_splitlist;
  atlas_splitlist.Split(names[0].c_str());
  atlas_splitlist.FillParams(atlas_list);
  OXS_GET_EXT_OBJECT(atlas_list,Oxs_Atlas,atlas);

  atlas_key.Set(atlas.GetPtr());
  atlas_key.GetReadReference();
  // Hold read lock until *this is deleted.

  for(OC_UINT4m i=1;i<name_count;i++) {
    OC_INT4m id = atlas->GetRegionId(names[i]);
    if(id<0) {
      char buf[1024];
      Oc_Snprintf(buf,sizeof(buf),
		  "Input error for evolver object \"%s\";"
		  " Section \"%s\""
		  " is not a known region in atlas \"%s\".  "
		  "Known regions:",
		  InstanceName(),names[i].c_str(),
		  atlas->InstanceName());
      String msg = String(buf);
      vector<String> regions;
      atlas->GetRegionList(regions);
      for(unsigned int j=0;j<regions.size();++j) {
	msg += String("\n ");
	msg += regions[j];
      }
      throw Oxs_ExtError(msg);
    }
    fixed_region_ids.push_back(id);
  }

}

void Oxs_Evolver::UpdateFixedSpinList(const Oxs_Mesh* mesh)
{
  if(mesh_id == mesh->Id()) return; // Mesh already set
  mesh_id = mesh->Id();
  fixed_spin_list.clear();
  OC_UINT4m fixed_region_count
    = static_cast<OC_UINT4m>(fixed_region_ids.size());
  OC_INDEX mesh_size = mesh->Size();
  if(atlas_key.GetPtr()==NULL || fixed_region_count==0 || mesh_size==0)
    return; // No fixed spins
  for(OC_INDEX i=0;i<mesh_size;i++) {
    // Find center point associated with spin i
    ThreeVector location;
    mesh->Center(i,location);
    // Find region containing spin i center point
    OC_INT4m region_id = atlas->GetRegionId(location);
    // Check to see if region_id is a fixed spin region
    if(region_id>=0) {
      for(OC_UINT4m j=0;j<fixed_region_count;j++) {
	if(region_id==fixed_region_ids[j]) {
	  fixed_spin_list.push_back(i); // Record spin i in
	  break;                        // fixed spin list.
	}
      }
    }
  }

  // Sort into increasing order.  This may improve cache locality
  // for accesses by solvers, but in particular is assumed by
  // the threaded portion of the energy/field evaluation code,
  // Oxs_ComputeEnergies (see oxs/base/energy.cc).
  sort(fixed_spin_list.begin(),fixed_spin_list.end());

}


// Constructors

// Barebones constructor for children that manually interpret
// their argstr.  If using this, then SetFixedRegions should
// be called in the child constructor, mimicking the calls in
// the other Oxs_Evolver constructor, which automatically
// parses argstr.
Oxs_Evolver::Oxs_Evolver
( const char* name,     // Child instance id
  Oxs_Director* newdtr  // App director
  ) : Oxs_Ext(name,newdtr), mesh_id(0)
{}

// Constructor with automatic argstr parsing.  Also automatically
// interprets fixed regions list, if any.
Oxs_Evolver::Oxs_Evolver
( const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr    // MIF block argument string
  ) : Oxs_Ext(name,newdtr,argstr), mesh_id(0)
{
  // Initialized fixed spin data structures.
  vector<String> fspins;
  if(FindInitValue("fixed_spins",fspins)) {
    // Fixed spin request
    InitFixedRegions(fspins);
    DeleteInitValue("fixed_spins");
  }
}

//Destructor
Oxs_Evolver::~Oxs_Evolver() {}

// Parent init function
OC_BOOL Oxs_Evolver::Init()
{
  if(!Oxs_Ext::Init()) return 0;

  mesh_id = 0;
  fixed_spin_list.clear();
  return 1;
}
