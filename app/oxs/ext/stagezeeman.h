/* FILE: stagezeeman.h         -*-Mode: c++-*-
 *
 * Zeeman (applied field) energy, derived from Oxs_Energy class.
 * The field may be spatially varying, and may change from stage
 * to stage, but is fixed within each stage.  The fields may be
 * specified either by a list of filenames (one per stage), or
 * else by a script that generates a vector field initialization
 * string.
 */

#ifndef _OXS_STAGEZEEMAN
#define _OXS_STAGEZEEMAN

#include <string>
#include <vector>

#include "energy.h"
#include "output.h"
#include "tclcommand.h"
#include "threevector.h"
#include "vectorfield.h"

OC_USE_STD_NAMESPACE;
OC_USE_STRING;

/* End includes */

class Oxs_StageZeeman:public Oxs_Energy {
private:
  OC_REAL8m hmult;
  OC_UINT4m number_of_stages;

  mutable OC_UINT4m mesh_id;
  mutable OC_BOOL stage_valid;
  mutable OC_UINT4m working_stage;     // Stage index
  mutable Oxs_OwnedPointer<Oxs_VectorField> stagefield_init;
  mutable Oxs_MeshValue<ThreeVector> stagefield;
  mutable ThreeVector max_field;
  /// stagefield is a cached value filled by
  /// stagefield_init when a change in mesh or
  /// stage is detected.  stagefield_init is updated
  /// when a change in stage is detected.  The stage_valid
  /// boolean is set false in the constructor and Init();
  /// It's purpose is to insure proper setup of stagefield_init
  /// on the first pass.
  ///  The max_field cache value is the field at the point
  /// of maximum magnitude in stagefield.  It is cached
  /// for Bapp output below.  Units are A/m.

  // Vector field may be specified by *either* a list of
  // files, or else a Tcl command that returns a vector
  // field spec.  This is set up in the Oxs_StageZeeman
  // constructor.  The choice is determined elsewhere by
  // examining the length of filelist; if it is >0, then
  // the list of files method is used, otherwise cmd is
  // called with the stage number as the single appended
  // argument.
  vector<String> filelist;
  Oxs_TclCommand cmd;

  // Cache update routines
  void ChangeFieldInitializer(OC_UINT4m stage,const Oxs_Mesh* mesh) const;
  void FillStageFieldCache(const Oxs_Mesh* mesh) const;
  void UpdateCache(const Oxs_SimState& state) const;

  // Supplied outputs, in addition to those provided by Oxs_Energy.
  // These return the applied field at the point of maximum applied
  // field magnitude.
  Oxs_ScalarOutput<Oxs_StageZeeman> Bapp_output;
  Oxs_ScalarOutput<Oxs_StageZeeman> Bappx_output;
  Oxs_ScalarOutput<Oxs_StageZeeman> Bappy_output;
  Oxs_ScalarOutput<Oxs_StageZeeman> Bappz_output;
  void Fill__Bapp_output(const Oxs_SimState& state);

protected:
  virtual void GetEnergy(const Oxs_SimState& state,
			 Oxs_EnergyData& oed) const;

public:
  virtual const char* ClassName() const; // ClassName() is
  /// automatically generated by the OXS_EXT_REGISTER macro.
  Oxs_StageZeeman(const char* name,     // Child instance id
		  Oxs_Director* newdtr, // App director
		  const char* argstr);  // MIF input block parameters
  virtual ~Oxs_StageZeeman();
  virtual OC_BOOL Init();
  virtual void StageRequestCount(unsigned int& min,
				 unsigned int& max) const;
};


#endif // _OXS_STAGEZEEMAN