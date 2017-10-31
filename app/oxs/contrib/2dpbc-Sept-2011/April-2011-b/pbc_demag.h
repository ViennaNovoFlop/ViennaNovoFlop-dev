/* FILE: demag.h            -*-Mode: c++-*-
 *
 * Average H demag field across rectangular cells.  This is a modified
 * version of the simpledemag class, which uses symmetries in the
 * interaction coefficients to reduce memory usage.
 *
 */

#ifndef _PBC_DEMAG_2D
#define _PBC_DEMAG_2D

#include "energy.h"
#include "vf.h"
#include "fft.h"
#include "key.h"
#include "mesh.h"
#include "output.h"
#include "meshvalue.h"
#include "simstate.h"
#include "threevector.h"
#include "rectangularmesh.h"

#include "pbc_util.h"

/* End includes */

class PBC_Demag_2D : public Oxs_Energy {
private:
#if REPORT_TIME
    mutable Nb_StopWatch ffttime;
    mutable Nb_StopWatch convtime;
    mutable Nb_StopWatch dottime;
#endif // REPORT_TIME

    mutable UINT4m rdimx; // Natural size of real data
    mutable UINT4m rdimy; // Digital Mars compiler wants these as separate
    mutable UINT4m rdimz; //    statements, because of "mutable" keyword.
    mutable UINT4m cdimx; // Full size of complex data
    mutable UINT4m cdimy;
    mutable UINT4m cdimz;
    // 2*cdimx>=rdimx, cdimy>=rdimy, cdimz>=rdimz
    // cdim[xyz] should be powers of 2.
    mutable UINT4m cstridey; // Strides across complex data
    mutable UINT4m cstridez;
    // cstridey>=cdimx, cstridez>=cdimy*cstridey
    // The stride sizes for the real arrays are just double the
    // complex strides, except cstride1 and rstride1 are assumed
    // to be 1.  Total matrix size is effectively cdimz*cstridez
    // Oxs_Complex elements, or twice that many "double" elements.

    mutable UINT4m mesh_id;

    // The A## arrays hold demag coefficients, transformed into
    // frequency domain.  These are held long term.  xcomp,
    // ycomp, and zcomp are used as temporary space, first to hold
    // the transforms of Mx, My, and Mz, then to store Hx, Hy, and
    // Hz.
    mutable UINT4m adimx;
    mutable UINT4m adimy;
    mutable UINT4m adimz;
    mutable UINT4m astridey;
    mutable UINT4m astridez;
    mutable OXS_COMPLEX_REAL_TYPE *A00;
    mutable OXS_COMPLEX_REAL_TYPE *A01;
    mutable OXS_COMPLEX_REAL_TYPE *A02;
    mutable OXS_COMPLEX_REAL_TYPE *A11;
    mutable OXS_COMPLEX_REAL_TYPE *A12;
    mutable OXS_COMPLEX_REAL_TYPE *A22;
    mutable Oxs_Complex *xcomp;
    mutable Oxs_Complex *ycomp;
    mutable Oxs_Complex *zcomp;
    mutable Oxs_FFT3D fft; // All transforms are same size, so we need
    /// only one Oxs_FFT3D object.

    String tensor_file_name;
    mutable BOOL load_from_file_success;
    REAL8m pbc_2d_error;
    mutable UINT4m xdim;
    mutable UINT4m ydim;
    mutable UINT4m zdim;
    mutable UINT4m sample_repeat_n;
    mutable REAL8m asymptotic_radius;
    mutable REAL8m dipolar_radius;
    mutable REAL8m asymptotic_radius_sq;
    mutable REAL8m dipolar_radius_sq;
    mutable Oxs_MeshValue<ThreeVector> Npbc_diag;
    mutable Oxs_MeshValue<ThreeVector> Npbc_offdiag;

    void (PBC_Demag_2D::*fillcoefs)(const Oxs_Mesh*) const;
    //  void FillCoefficientArraysFast(const Oxs_Mesh* mesh) const;
    void FillCoefficientArraysStandard(const Oxs_Mesh* mesh) const;
    /// The "standard" variant is simpler but slower, and is retained
    /// mainly for testing and development purposes.

    UINT4m NextPowerOfTwo(UINT4m n) const; // Helper function
    void ReleaseMemory() const;
    double CalculateSingleTensor(enum TensorComponent comp,
            int g, double x, double y, double z, double a, double b, double c) const;
    int FindG(enum TensorComponent comp, double v, double Tx, double Ty) const;
    void CalculateDemagTensors(const Oxs_RectangularMesh* mesh) const;
    void SavePbcDemagTensor(const Oxs_Mesh *mesh) const;
    void LoadPbcDemagTensor(const Oxs_RectangularMesh* mesh) const;
    REAL8m GetTensorFromBuffer(enum TensorComponent comp, int i, int j, int k) const;
    

protected:
    virtual void GetEnergy(const Oxs_SimState& state,
            Oxs_EnergyData& oed) const;

public:
    virtual const char* ClassName() const; // ClassName() is
    /// automatically generated by the OXS_EXT_REGISTER macro.
    PBC_Demag_2D(const char* name, // Child instance id
            Oxs_Director* newdtr, // App director
            const char* argstr); // MIF input block parameters
    virtual ~PBC_Demag_2D();
    virtual BOOL Init();
};


#endif // _PBC_Demag_2D
