/* FILE: rectangularmesh.h            -*-Mode: c++-*-
 *
 * Cell rectangular mesh, derived from Oxs_Mesh class.
 *
 */

#ifndef _OXS_RECTANGULARMESH
#define _OXS_RECTANGULARMESH

#include <string>
#include <vector>

#include "atlas.h"
#include "scalarfield.h"
#include "mesh.h"
#include "vf.h"

OC_USE_STD_NAMESPACE;
OC_USE_STRING;

/* End includes */


class Oxs_RectangularMesh : public Oxs_Mesh {
private:
  ThreeVector base,cellsize;
  OC_REAL8m cellvolume; // Convenience value
  OC_INDEX xdim,ydim,zdim;
  OC_INDEX xydim; // Convenience value: xdim*ydim
  OC_INDEX elementcount;  // Convenience value: xdim*ydim*zdim

  OC_BOOL GetNeighborPoint(const ThreeVector& base_pt,
			OC_INDEX ngbr_index,
			ThreeVector& ngbr_pt) const;
  // Fills ngbr_pt with location of neighbor element indexed
  // by "ngbr_index", relative to base_pt.  Returns 1 if
  // ngbr_pt < number of neighbors (currently 6); otherwise
  // 0 is returned, in which case ngbr_pt is unchanged.
  // NB: ngbr_index is 0-based.

  // Helper function for constructors
  void InitScaling(const Oxs_Box& box);

  // Constructor for private use by MakeRefinedMesh member function.
  Oxs_RectangularMesh(const char* name,
		      Oxs_Director* newdtr,
		      const ThreeVector& in_base,
		      const ThreeVector& in_cellsize,
		      OC_INDEX in_xdim,OC_INDEX in_ydim,OC_INDEX in_zdim);

  // Dummy class for throwing exceptions; this is intended for
  // internal use only
  class Oxs_RectangularMesh_WBError {};

public:
  // Oxs_Ext interface
  Oxs_RectangularMesh
  (const char* name,     // Child instance id
   Oxs_Director* newdtr, // App director
   const char* argstr);  // MIF input block parameters
  virtual ~Oxs_RectangularMesh() {}

  virtual const char* ClassName() const; // ClassName() is
  /// automatically generated by the OXS_EXT_REGISTER macro.

  // Mesh interface
  virtual void GetBoundingBox(Oxs_Box& bbox) const;
  virtual OC_INDEX Size() const { return elementcount; }
  virtual OC_REAL8m TotalVolume() const { return cellvolume*elementcount; }
  virtual OC_REAL8m Volume(OC_INDEX) const { return cellvolume; }
  /// NOTE: Volume() is thread-safe, provided that the underlying mesh
  /// object is not changed in another thread.
  virtual void Center(OC_INDEX index,ThreeVector& location) const;
  virtual OC_INDEX FindNearestIndex(const ThreeVector& location) const;
  virtual OC_BOOL HasUniformCellVolumes() const { return 1; }
  virtual OC_BOOL HasUniformCellVolumes(OC_REAL8m& vol) const {
    vol = cellvolume;
    return 1;
 }

  // Vf_Ovf20_MeshNodes interface
  virtual OC_BOOL SupportsOutputType(Vf_Ovf20_MeshType type) const {
    OC_BOOL i_do = 0;
    switch(type) {
    case vf_ovf20mesh_irregular:
    case vf_ovf20mesh_rectangular:
      i_do = 1;
      break;
    default:
      i_do = 0;  // Safety
      break;
    }
    return i_do;
  }

  virtual Vf_Ovf20_MeshType NaturalType() const {
    return vf_ovf20mesh_rectangular;
  }

  virtual void DumpGeometry(Vf_Ovf20FileHeader& header,
                            Vf_Ovf20_MeshType type) const;


  // File (channel) output routines.  These throw an exception on error.
  virtual void WriteOvf
  (Tcl_Channel channel,   // Output channel
   OC_BOOL headers,          // If false, then output only raw data
   const char* title,     // Long filename or title
   const char* desc,      // Description to embed in output file
   const char* valueunit, // Field units, such as "A/m".
   const char* meshtype,  // Either "rectangular" or "irregular"
   const char* datatype,  // Either "binary" or "text"
   const char* precision, // For binary, "4" or "8";
                          ///  for text, a printf-style format
   const Oxs_MeshValue<ThreeVector>* vec,  // Vector array
   const Oxs_MeshValue<OC_REAL8m>* scale=NULL // Optional scaling for vec
   /// Set scale to NULL to use vec values directly.
   ) const;

  // Conversion routines from Vf_Mesh to Oxs_MeshValue<ThreeVector>.
  // IsCompatible returns true iff vfmesh is a Vf_GridVec3f with
  //   dimensions identical to those of *this.
  // NB: IsCompatible only compares the mesh dimensions, not the
  //   underlying physical scale, or the cell aspect ratios.  Do
  //   we want to include such a check?, or is it more flexible
  //   to leave it out?
  // FillMeshValueExact copies the vector field held in vfmesh to the
  //   export Oxs_MeshValue<ThreeVector> vec.  This routine throws
  //   an exception on error, the primary cause of which is that
  //   vfmesh is not compatible with *this.  In other words, if you
  //   don't want to catch the exception, call IsCompatible first.
  //   The "Exact" in the name refers to the requirement that the
  //   dimensions on vfmesh exactly match thos of *this.
  OC_BOOL IsCompatible(const Vf_Mesh* vfmesh) const;
  void FillMeshValueExact(const Vf_Mesh* vfmesh,
                          Oxs_MeshValue<ThreeVector>& vec) const;

  virtual const char* NaturalOutputType() const { return "rectangular"; }

  // Summing routines.  The terms are weighted by the basis element volumes.
  virtual OC_REAL8m VolumeSum(const Oxs_MeshValue<OC_REAL8m>& scalar) const;
  virtual ThreeVector VolumeSum(const Oxs_MeshValue<ThreeVector>& vec) const;
  virtual ThreeVector VolumeSum(const Oxs_MeshValue<ThreeVector>& vec,
				const Oxs_MeshValue<OC_REAL8m>& scale) const;
  /// "scale" is elementwise scaling.

  // Extra precision summing routines.  These are the same as above,
  // except that intermediate sums are accumulated into Nb_Xpfloat
  // variables.
  virtual OC_REAL8m VolumeSumXp(const Oxs_MeshValue<OC_REAL8m>& scalar) const;
  virtual ThreeVector VolumeSumXp(const Oxs_MeshValue<ThreeVector>& vec) const;
  virtual ThreeVector VolumeSumXp(const Oxs_MeshValue<ThreeVector>& vec,
				  const Oxs_MeshValue<OC_REAL8m>& scale) const;

  // Max angle routine.  Every mesh should have some concept of a
  // "neighborhood."  This routine should return the minimum dot
  // product between neighboring vectors across the mesh.  If the
  // input Oxs_MeshValue<ThreeVector> array are unit vectors, then
  // the max angle can be construed from acos(mindot).
  // NB: This routine is boundary condition sensitive.
  virtual OC_REAL8m MinNeighborDot(const Oxs_MeshValue<ThreeVector>& vec,
				const Oxs_MeshValue<OC_REAL8m>& zero_check)
    const;

  // Rectangular mesh specific functions.
  // IMPORTANT NOTE: Client routines can assume that the storage
  // order is first on x, second on y, third on z (Fortran order),
  // i.e., (0,0,0) (1,0,0) ... (xdim-1,0,0) (0,1,0) (1,1,0) ...
  // (xdim-1,ydim-1,0) (0,0,1) (1,0,1) ... (xdim-1,ydim-1,zdim-1)
  // This locks down the representation in this function, but
  // allows improved access in the clients.  If we need to change
  // representations, we can always introduce a new mesh class. ;^)
  OC_REAL8m EdgeLengthX() const { return cellsize.x; }
  OC_REAL8m EdgeLengthY() const { return cellsize.y; }
  OC_REAL8m EdgeLengthZ() const { return cellsize.z; }
  OC_INDEX DimX() const { return xdim; }
  OC_INDEX DimY() const { return ydim; }
  OC_INDEX DimZ() const { return zdim; }

  // Secondary constructor; provides a function-level API
  // for use by other Oxs_Ext objects.
  Oxs_RectangularMesh
  (const char* name,     // Child instance id
   Oxs_Director* newdtr, // App director
   const char* argstr,   // MIF input block parameters, for parent use
   const ThreeVector& in_cellsize,
   const Oxs_Box& range_box);

  // Convenience function for clients that use multiple translations
  // of a single mesh.
  void TranslateBasePt(const ThreeVector& newbasept) { base = newbasept; }

  // Mesh refinement, via generation of a new mesh.
  void MakeRefinedMesh(const char* name,
	       OC_INDEX refine_x,OC_INDEX refine_y,OC_INDEX refine_z,
	       Oxs_OwnedPointer<Oxs_RectangularMesh>& out_mesh) const;

  // Boundary extraction.  Returns list (technically, an STL vector) of
  // indices for those elements inside base_region that have a neighbor
  // (in the 6 nearest ngbr sense) such that the first element lies on
  // one side (the "inside") of the surface specified by the
  // Oxs_ScalarField bdry_surface + bdry_value, and the neighbor lies on
  // the other (the "outside").  If the bdry_side argument is "-", then
  // the "inside" of the surface is the set of those points x for which
  // bdry_surface.Value(x)<bdry_value.  If bdry_side is "+", then the
  // "inside" of the surface is the set of those points x for which
  // bdry_surface.Value(x)>bdry_value.
  // Return value is the number of entries in the export list.
  // NB: The tested neighbors may lie outside the mesh proper, allowing
  // elements on the edge of the mesh to be specified.
  OC_INDEX BoundaryList(const Oxs_Atlas& atlas,
		      const String& region_name,
		      const Oxs_ScalarField& bdry_surface,
		      OC_REAL8m bdry_value,
		      const String& bdry_side,
		      vector<OC_INDEX> &BoundaryIndexList) const;

  // Use the following function to convert a triplet index (x,y,z)
  // to the corresponding address used to access linear MeshValue
  // arrays.
  inline OC_INDEX Index(OC_INDEX x,OC_INDEX y,OC_INDEX z) const {
    return x+y*xdim+z*xydim;
  }

  inline void GetCoords(OC_INDEX index,
                        OC_INDEX& x,OC_INDEX& y,OC_INDEX& z) const {
    z = index/xydim;
    y = (index - z*xydim)/xdim;
    x = index - z*xydim - y*xdim;
  }
  void GetCenter(OC_INDEX index, ThreeVector& vec) const;


};

#endif // _OXS_RECTANGULARMESH
