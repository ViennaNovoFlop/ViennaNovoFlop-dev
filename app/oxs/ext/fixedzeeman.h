/* FILE: fixedzeeman.h            -*-Mode: c++-*-
 *
 * Fixed (in time) Zeeman energy/field, derived from Oxs_Energy class.
 *
 */

#ifndef _OXS_FIXEDZEEMAN
#define _OXS_FIXEDZEEMAN

#include "oc.h"
#include "director.h"
#include "threevector.h"
#include "energy.h"
#include "simstate.h"
#include "mesh.h"
#include "meshvalue.h"
#include "vectorfield.h"

/* End includes */

class Oxs_FixedZeeman:public Oxs_Energy {
private:
  mutable OC_UINT4m mesh_id;
  OC_REAL8m field_mult;
  Oxs_OwnedPointer<Oxs_VectorField> fixedfield_init;
  mutable Oxs_MeshValue<ThreeVector> fixedfield;
  /// fixedfield is a cached value filled by
  /// fixedfield_init when a change in mesh is
  /// detected.

protected:
  virtual void GetEnergy(const Oxs_SimState& state,
			 Oxs_EnergyData& oed) const;

public:
  virtual const char* ClassName() const; // ClassName() is
  /// automatically generated by the OXS_EXT_REGISTER macro.
  Oxs_FixedZeeman(const char* name,  // Child instance id
		  Oxs_Director* newdtr, // App director
		  const char* argstr);  // MIF input block parameters

  virtual ~Oxs_FixedZeeman() {}
  virtual OC_BOOL Init();
};


#endif // _OXS_FIXEDZEEMAN
