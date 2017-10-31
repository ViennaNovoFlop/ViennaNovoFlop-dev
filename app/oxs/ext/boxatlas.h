/* FILE: boxatlas.h                 -*-Mode: c++-*-
 *
 * Atlas class derived from Oxs_Atlas class that uses an import
 * image file for region demarcation.
 */

#ifndef _OXS_BOXATLAS
#define _OXS_BOXATLAS

#include <string>

#include "threevector.h"
#include "util.h"
#include "atlas.h"

OC_USE_STRING;

/* End includes */

class Oxs_BoxAtlas:public Oxs_Atlas {
private:
  String region_name;
  Oxs_Box world;
public:
  virtual const char* ClassName() const; // ClassName() is
  /// automatically generated by the OXS_EXT_REGISTER macro.

  Oxs_BoxAtlas(const char* name,
	       Oxs_Director* newdtr,
	       const char* argstr);   // MIF block argument string

  ~Oxs_BoxAtlas() {}

  void GetWorldExtents(Oxs_Box &mybox) const { mybox = world; }
  /// Fills mybox with bounding box for the atlas.

  OC_BOOL GetRegionExtents(OC_UINT4m id,Oxs_Box &mybox) const;
  /// If id is 0 or 1, sets mybox to world and returns 1.
  /// If id > 1, leaves mybox untouched and returns 0.

  OC_INT4m GetRegionId(const ThreeVector& point) const;
  /// Returns the id number for the region containing point.
  /// The return value is 0 if the point is not contained in
  /// the atlas, i.e., belongs to the "universe" region.

  OC_INT4m GetRegionId(const String& name) const;
  /// Given a region id string (name), returns
  /// the corresponding region id index.  If
  /// "name" is not included in the atlas, then
  /// -1 is returned.  Note: If name == "universe",
  /// then the return value will be 0.

  OC_BOOL GetRegionName(OC_UINT4m id,String& name) const;
  /// Given an id number, fills in "name" with
  /// the corresponding region id string.  Returns
  /// 1 on success, 0 if id is invalid.  If id is 0,
  /// then name is set to "universe", and the return
  /// value is 1.

  OC_UINT4m GetRegionCount() const { return 2; }
  /// Two regions: implicit "universe" pseudo-region, and
  /// that specified by name and identified with the world
  /// bbox.
};

#endif // _OXS_BOXATLAS

