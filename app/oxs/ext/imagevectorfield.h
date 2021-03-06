/* FILE: imagevectorfield.h      -*-Mode: c++-*-
 *
 * Image vector field object, derived from Oxs_VectorField class.
 *
 */

#ifndef _OXS_IMAGEVECTORFIELD
#define _OXS_IMAGEVECTORFIELD

#include "oc.h"

#include "vectorfield.h"

/* End includes */

class Oxs_ImageVectorField:public Oxs_VectorField {
public:
  virtual const char* ClassName() const; // ClassName() is
  /// automatically generated by the OXS_EXT_REGISTER macro.

  Oxs_ImageVectorField
  (const char* name,     // Child instance id
   Oxs_Director* newdtr, // App director
   const char* argstr);  // MIF input block parameters

  virtual ~Oxs_ImageVectorField();

  virtual void Value(const ThreeVector& pt,ThreeVector& value) const;

  virtual void FillMeshValue(const Oxs_Mesh* mesh,
			     Oxs_MeshValue<ThreeVector>& array) const
  { DefaultFillMeshValue(mesh,array); }

private:
  Oxs_Box bbox;

  enum ViewPlane { xy, zx, yz } view;
  // View plane selection:
  //  Default is "xy", for which
  //   x increases from left to right of image, and
  //   y increases from bottom to top of image.
  //  Selection "zx" specifies
  //   z increases from left to right of image, and
  //   x increases from bottom to top of image.
  //  Selection "yz" specifies
  //   y increases from left to right of image, and
  //   z increases from bottom to top of image.

  enum Exterior_Handling {
    EXTERIOR_INVALID, EXTERIOR_ERROR, EXTERIOR_BOUNDARY, EXTERIOR_DEFAULT
  } exterior;
  ThreeVector default_value;

  OC_INT4m image_width,image_height;

  vector<ThreeVector> value_array;

  OC_REAL8m column_offset,column_scale;
  OC_REAL8m row_offset,row_scale;
};

#endif // _OXS_IMAGEVECTORFIELD
