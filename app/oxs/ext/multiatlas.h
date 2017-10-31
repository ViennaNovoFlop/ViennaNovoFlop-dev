/* FILE: multiatlas.h                 -*-Mode: c++-*-
 *
 * Combining atlas class, derived from Oxs_Atlas class.
 *
 */

#ifndef _OXS_MULTIATLAS
#define _OXS_MULTIATLAS

#include <string>
#include <vector>

#include "ext.h"
#include "threevector.h"
#include "util.h"
#include "atlas.h"

OC_USE_STD_NAMESPACE;
OC_USE_STRING;

/* End includes */

class Oxs_MultiAtlas:public Oxs_Atlas {
private:
  // Global regions
  vector<String> region_name; // List of region names.
  vector<Oxs_Box> region_bbox; // Region bounding boxes
  /// NB: region_name[0]=="universe", and region_bbox[0]
  ///  is the atlas bounding box.

  // Subatlases
  struct SubAtlas {
  public:
    Oxs_OwnedPointer<Oxs_Atlas> atlas;
    vector<OC_UINT4m> global_id; // Maps local region id's to
    /// multiatlas global id's.
  };
  Nb_ArrayWrapper<SubAtlas> subatlas;

  // NOTE: subatlas.size() is the number of subatlases.
  // region_name.size()==region_bbox.size() which is the
  // total number of distinct regions.  Distinctness of
  // regions is determined by name, via case sensitive
  // compare.  Usually the number of distinct regions will
  // be at least as large as the number of subatlases,
  // unless single regions are uses by multiple atlases.

public:
  virtual const char* ClassName() const; // ClassName() is
  /// automatically generated by the OXS_EXT_REGISTER macro.

  Oxs_MultiAtlas(const char* name,
		  Oxs_Director* newdtr,
		  const char* argstr);   // MIF block argument string

  ~Oxs_MultiAtlas() {};

  void GetWorldExtents(Oxs_Box &mybox) const { mybox = region_bbox[0]; }
  /// Fills mybox with bounding box for the atlas,
  /// excluding "universe," natch.

  OC_BOOL GetRegionExtents(OC_UINT4m id,Oxs_Box &mybox) const;
  /// If 0<id<GetRegionCount, then sets mybox to bounding box
  /// for the specified region, and returns 1.  If id is 0,
  /// sets mybox to atlas bounding box, and returns 1.
  /// Otherwise, leaves mybox untouched and returns 0.

  OC_INT4m GetRegionId(const ThreeVector& point) const;
  /// Returns the id number for the region containing point.
  /// The return value is 0 if the point is not contained in
  /// the atlas, i.e., belongs to the "universe" region.
  /// The return value is guaranteed to be in the range
  /// [0,GetRegionCount()-1].

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

  OC_UINT4m GetRegionCount() const {
    return static_cast<OC_UINT4m>(region_name.size());
  }
  /// Valid RegionId numbers range from 0 to GetRegionCount() - 1,
  /// inclusive, where 0 is the special "universe" region.
};

#endif // _OXS_MULTIATLAS
