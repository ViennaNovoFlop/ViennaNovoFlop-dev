/* FILE: demag.cc            -*-Mode: c++-*-
 *
 * Average H demag field across rectangular cells.  This is a modified
 * version of the simpledemag class, which uses symmetries in the
 * interaction coefficients to reduce memory usage.
 *
 * The formulae used are reduced forms of equations in A. J. Newell,
 * W. Williams, and D. J. Dunlop, "A Generalization of the Demagnetizing
 * Tensor for Nonuniform Magnetization," Journal of Geophysical Research
 * - Solid Earth 98, 9551-9555 (1993).
 *
 */

#include <string>

#include "nb.h"
#include "director.h"
#include "key.h"
#include "oxsexcept.h"
#include "mesh.h"
#include "meshvalue.h"
#include "simstate.h"
#include "threevector.h"
#include "energy.h"		// Needed to make MSVC++ 5 happy

#include "rectangularmesh.h"

#include "demagcoef.h"
#include "fft.h"

#include "pbc_demag.h"
#include "pbc_util.h"

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(PBC_Demag_2D);

/* End includes */

// Helper function

OC_UINT4m PBC_Demag_2D::NextPowerOfTwo(OC_UINT4m n) const { // Returns first power of two >= n
    OC_UINT4m m = 1;
    while (m < n) m *= 2;
    return m;
}

// Constructor

PBC_Demag_2D::PBC_Demag_2D(
        const char* name, // Child instance id
        Oxs_Director* newdtr, // App director
        const char* argstr) // MIF input block parameters
: Oxs_Energy(name, newdtr, argstr),
rdimx(0), rdimy(0), rdimz(0),
cdimx(0), cdimy(0), cdimz(0), cstridey(0), cstridez(0),
mesh_id(0),
adimx(0), adimy(0), adimz(0), astridey(0), astridez(0),
A00(NULL), A01(NULL), A02(NULL), A11(NULL), A12(NULL), A22(NULL),
xcomp(NULL), ycomp(NULL), zcomp(NULL), fillcoefs(0),
tensor_file_name(""), pbc_2d_error(0),
asymptotic_radius(0), dipolar_radius(0),
asymptotic_radius_sq(0), dipolar_radius_sq(0),
xdim(0), ydim(0), zdim(0), sample_repeat_n(0),
Npbc_diag(NULL), Npbc_offdiag(NULL) {

    fillcoefs = &PBC_Demag_2D::FillCoefficientArraysStandard;

    tensor_file_name = GetStringInitValue("tensor_file_name", "");
    pbc_2d_error = GetRealInitValue("tensor_error", 1e-10);
    asymptotic_radius = GetRealInitValue("asymptotic_radius", 32.0);
    dipolar_radius = GetRealInitValue("dipolar_radius", 10000.0);
    sample_repeat_n = GetIntInitValue("sample_repeat_number", 0);

    asymptotic_radius_sq = asymptotic_radius*asymptotic_radius;
    dipolar_radius_sq = dipolar_radius*dipolar_radius;

    VerifyAllInitArgsUsed();
}

PBC_Demag_2D::~PBC_Demag_2D() {
#if REPORT_TIME
    Oc_TimeVal cpu, wall;
    ffttime.GetTimes(cpu, wall);
    if (OC_REALWIDE(wall) > 0.0) {
        fprintf(stderr, "      subtime ...   fft%7.2f cpu /%7.2f wall,"
                " (%s)\n",
                OC_REALWIDE(cpu), OC_REALWIDE(wall), InstanceName());
        convtime.GetTimes(cpu, wall);
        fprintf(stderr, "      subtime ...  conv%7.2f cpu /%7.2f wall,"
                " (%s)\n",
                OC_REALWIDE(cpu), OC_REALWIDE(wall), InstanceName());
        dottime.GetTimes(cpu, wall);
        fprintf(stderr, "      subtime ...   dot%7.2f cpu /%7.2f wall,"
                " (%s)\n",
                OC_REALWIDE(cpu), OC_REALWIDE(wall), InstanceName());
    }
#endif // REPORT_TIME
    ReleaseMemory();
}

OC_BOOL PBC_Demag_2D::Init() {
#if REPORT_TIME
    Oc_TimeVal cpu, wall;
    ffttime.GetTimes(cpu, wall);
    if (OC_REALWIDE(wall) > 0.0) {
        fprintf(stderr, "      subtime ...   fft%7.2f cpu /%7.2f wall,"
                " subtime module %s\n",
                OC_REALWIDE(cpu), OC_REALWIDE(wall), InstanceName());
        convtime.GetTimes(cpu, wall);
        fprintf(stderr, "              ...  conv%7.2f cpu /%7.2f wall,"
                " subtime module %s\n",
                OC_REALWIDE(cpu), OC_REALWIDE(wall), InstanceName());
        dottime.GetTimes(cpu, wall);
        fprintf(stderr, "      subtime ...   dot%7.2f cpu /%7.2f wall,"
                " (%s)\n",
                OC_REALWIDE(cpu), OC_REALWIDE(wall), InstanceName());
    }
    ffttime.Reset();
    convtime.Reset();
    dottime.Reset();
#endif // REPORT_TIME
    mesh_id = 0;
    ReleaseMemory();
    return Oxs_Energy::Init();
}

void PBC_Demag_2D::ReleaseMemory() const { // Conceptually const
    if (xcomp != NULL) {
        delete[] xcomp;
        xcomp = NULL;
    }
    if (ycomp != NULL) {
        delete[] ycomp;
        ycomp = NULL;
    }
    if (zcomp != NULL) {
        delete[] zcomp;
        zcomp = NULL;
    }
    if (A00 != NULL) {
        delete[] A00;
        A00 = NULL;
    }
    if (A01 != NULL) {
        delete[] A01;
        A01 = NULL;
    }
    if (A02 != NULL) {
        delete[] A02;
        A02 = NULL;
    }
    if (A11 != NULL) {
        delete[] A11;
        A11 = NULL;
    }
    if (A12 != NULL) {
        delete[] A12;
        A12 = NULL;
    }
    if (A22 != NULL) {
        delete[] A22;
        A22 = NULL;
    }
    rdimx = rdimy = rdimz = 0;
    cdimx = cdimy = cdimz = 0;
    cstridey = cstridez = 0;
    adimx = adimy = adimz = 0;
    astridey = astridez = 0;
}

void PBC_Demag_2D::FillCoefficientArraysStandard(const Oxs_Mesh* genmesh) const { // This routine is conceptually const.
    // This "OLD" variant is simpler but slower, and is retained
    // solely for testing and development purposes.

    const Oxs_RectangularMesh* mesh
            = dynamic_cast<const Oxs_RectangularMesh*> (genmesh);
    if (mesh == NULL) {
        String msg = String("Object ")
                + String(genmesh->InstanceName())
                + String(" is not a rectangular mesh.");
        throw Oxs_Ext::Error(this, msg);
    }

    // Clean-up from previous allocation, if any.
    ReleaseMemory();

    // Fill dimension variables
    rdimx = mesh->DimX();
    rdimy = mesh->DimY();
    rdimz = mesh->DimZ();
    if (rdimx == 0 || rdimy == 0 || rdimz == 0) return; // Empty mesh!

    cdimx = NextPowerOfTwo(2 * rdimx) / 2;
    if (rdimy > 1) cdimy = NextPowerOfTwo(2 * rdimy);
    else cdimy = 1;
    if (rdimz > 1) cdimz = NextPowerOfTwo(2 * rdimz);
    else cdimz = 1;
    cstridey = cdimx + 1; // Pad by one to avoid cache line entanglement
    cstridez = cstridey*cdimy;

    OC_UINT4m ctotalsize = cstridez*cdimz;
    OC_UINT4m rtotalsize = 2 * ctotalsize;
    if (rtotalsize < 2 * cdimx || rtotalsize < cdimy || rtotalsize < cdimz) {
        // Partial overflow check
        String msg = String("OC_UINT4m overflow in ") + String(InstanceName())
                + String(": Product cdimx*cdimy*cdimz too big"
                " to fit in a OC_UINT4m variable");
        throw Oxs_Ext::Error(this, msg);
    }

    // Allocate memory for interaction matrices and magnetization components
    adimx = 1 + cdimx;
    adimy = 1 + cdimy / 2;
    adimz = 1 + cdimz / 2;
    astridey = adimx;
    astridez = adimy*astridey;
    OC_UINT4m atotalsize = adimz*astridez;
    A00 = new OXS_COMPLEX_REAL_TYPE[atotalsize];
    A01 = new OXS_COMPLEX_REAL_TYPE[atotalsize];
    A02 = new OXS_COMPLEX_REAL_TYPE[atotalsize];
    A11 = new OXS_COMPLEX_REAL_TYPE[atotalsize];
    A12 = new OXS_COMPLEX_REAL_TYPE[atotalsize];
    A22 = new OXS_COMPLEX_REAL_TYPE[atotalsize];
    xcomp = new Oxs_Complex[ctotalsize];
    ycomp = new Oxs_Complex[ctotalsize];
    zcomp = new Oxs_Complex[ctotalsize];
    if (xcomp == NULL || ycomp == NULL || zcomp == NULL
            || A00 == NULL || A01 == NULL || A02 == NULL
            || A11 == NULL || A12 == NULL || A22 == NULL) {
        // Safety check for those machines on which new[] doesn't throw
        // BadAlloc.
        String msg = String("Insufficient memory in Demag constructor.");
        throw Oxs_Ext::Error(this, msg);
    }

    OC_REALWIDE *rxcomp = static_cast<OC_REALWIDE*> (static_cast<void*> (xcomp));
    OC_REALWIDE *rycomp = static_cast<OC_REALWIDE*> (static_cast<void*> (ycomp));
    OC_REALWIDE *rzcomp = static_cast<OC_REALWIDE*> (static_cast<void*> (zcomp));
    OC_UINT4m rstridey = 2 * cstridey;
    OC_UINT4m rstridez = 2 * cstridez;

    // According (16) in Newell's paper, the demag field is given by
    //                        H = -N*M
    // where N is the "demagnetizing tensor," with components Nxx, Nxy,
    // etc.  With the '-1' in 'scale' we store '-N' instead of 'N',
    // so we don't have to multiply the output from the FFT + iFFT
    // by -1 in ConstMagField() below.

    // Fill interaction matrices with demag coefs from Newell's paper.
    // Note that A00, A11 and A22 are even in x,y and z.
    // A01 is odd in x and y, even in z.
    // A02 is odd in x and z, even in y.
    // A12 is odd in y and z, even in x.
    // We use these symmetries to reduce storage requirements.  If
    // f is real and even, then f^ is also real and even.  If f
    // is real and odd, then f^ is (pure) imaginary and odd.
    // As a result, the transform of each of the Axx interaction
    // matrices will be real, with the same even/odd properties.

    OC_UINT4m index, i, j, k;

    OC_REALWIDE dx = mesh->EdgeLengthX();
    OC_REALWIDE dy = mesh->EdgeLengthY();
    OC_REALWIDE dz = mesh->EdgeLengthZ();
    // For demag calculation, all that matters is the relative
    // sizes of dx, dy and dz.  To help insure we don't run
    // outside floating point range, rescale these values so
    // largest is 1.0
    OC_REALWIDE maxedge = dx;
    if (dy > maxedge) maxedge = dy;
    if (dz > maxedge) maxedge = dz;
    dx /= maxedge;
    dy /= maxedge;
    dz /= maxedge;

    // OC_REALWIDE scale = -1./(4*PI*dx*dy*dz);
    OC_REALWIDE scale = -1.0;


    // Also throw in FFT scaling.  This allows the "NoScale" FFT routines
    // to be used.
    scale *= 1.0 / (2 * cdimx * cdimy * cdimz);

    for (index = 0; index < ctotalsize; index++) {
        // Initialize temp space to 0
        xcomp[index].Set(0., 0.);
        ycomp[index].Set(0., 0.);
        zcomp[index].Set(0., 0.);
    }


#ifdef DUMP_COEF_TEST
    fprintf(stderr, "Nxy(1,2,3,1,2,3)=%.17g   Nxy(10,1,1,1,2,3)=%.17g\n",
            (OC_REALWIDE) CalculateSDA01(1., 2., 3., 1., 2., 3.) / (4 * PI * 1. * 2. * 3.),
            (OC_REALWIDE) CalculateSDA01(10., 1., 1., 1., 2., 3.) / (4 * PI * 1. * 2. * 3.));
    fprintf(stderr, "Nxy(-1,2,3,1,2,3)=%.17g   Nxy(10,1,-1,1,2,3)=%.17g\n",
            (OC_REALWIDE) CalculateSDA01(-1., 2., 3., 1., 2., 3.) / (4 * PI * 1. * 2. * 3.),
            (OC_REALWIDE) CalculateSDA01(10., 1., -1., 1., 2., 3.) / (4 * PI * 1. * 2. * 3.));
    fprintf(stderr, "Nxy(1,1,0,1,2,3)=%.17g   Nxy(1,1,0,2,1,3)=%.17g\n",
            (OC_REALWIDE) CalculateSDA01(1., 1., 0., 1., 2., 3.) / (4 * PI * 1. * 2. * 3.),
            (OC_REALWIDE) CalculateSDA01(1., 1., 0., 2., 1., 3.) / (4 * PI * 1. * 2. * 3.));
    fprintf(stderr, "Nxy(-1,1,0,1,2,3)=%.17g   Nxy(1,-1,0,2,1,3)=%.17g\n",
            (OC_REALWIDE) CalculateSDA01(-1., 1., 0., 1., 2., 3.) / (4 * PI * 1. * 2. * 3.),
            (OC_REALWIDE) CalculateSDA01(1., -1., 0., 2., 1., 3.) / (4 * PI * 1. * 2. * 3.));
#endif // DUMP_COEF_TEST

    for (k = 0; k < rdimz; k++) for (j = 0; j < rdimy; j++) for (i = 0; i < rdimx; i++) {
                //  OC_REALWIDE x = dx*i;
                //  OC_REALWIDE y = dy*j;
                //  OC_REALWIDE z = dz*k;
                //   OC_REALWIDE a00=scale*CalculateSDA00(x,y,z,dx,dy,dz);
                //   OC_REALWIDE a01=scale*CalculateSDA01(x,y,z,dx,dy,dz);
                //   OC_REALWIDE a02=scale*CalculateSDA02(x,y,z,dx,dy,dz);
                OC_REALWIDE a00 = scale * GetTensorFromBuffer(xx, i, j, k);
                OC_REALWIDE a01 = scale * GetTensorFromBuffer(xy, i, j, k);
                ;
                OC_REALWIDE a02 = scale * GetTensorFromBuffer(xz, i, j, k);
                ;

                index = i + j * rstridey + k*rstridez;

                rxcomp[index] = a00;
                rycomp[index] = a01;
                rzcomp[index] = a02;

                if (i > 0) {
                    OC_UINT4m tindex = (2 * cdimx - i) + j * rstridey + k*rstridez;
                    rxcomp[tindex] = a00;
                    rycomp[tindex] = -a01;
                    rzcomp[tindex] = -a02;
                }
                if (j > 0) {
                    OC_UINT4m tindex = i + (cdimy - j) * rstridey + k*rstridez;
                    rxcomp[tindex] = a00;
                    rycomp[tindex] = -a01;
                    rzcomp[tindex] = a02;
                }
                if (k > 0) {
                    OC_UINT4m tindex = i + j * rstridey + (cdimz - k) * rstridez;
                    rxcomp[tindex] = a00;
                    rycomp[tindex] = a01;
                    rzcomp[tindex] = -a02;
                }
                if (i > 0 && j > 0) {
                    OC_UINT4m tindex = (2 * cdimx - i)+(cdimy - j) * rstridey + k*rstridez;
                    rxcomp[tindex] = a00;
                    rycomp[tindex] = a01;
                    rzcomp[tindex] = -a02;
                }
                if (i > 0 && k > 0) {
                    OC_UINT4m tindex = (2 * cdimx - i) + j * rstridey + (cdimz - k) * rstridez;
                    rxcomp[tindex] = a00;
                    rycomp[tindex] = -a01;
                    rzcomp[tindex] = a02;
                }
                if (j > 0 && k > 0) {
                    OC_UINT4m tindex = i + (cdimy - j) * rstridey + (cdimz - k) * rstridez;
                    rxcomp[tindex] = a00;
                    rycomp[tindex] = -a01;
                    rzcomp[tindex] = -a02;
                }
                if (i > 0 && j > 0 && k > 0) {
                    OC_UINT4m tindex = (2 * cdimx - i)+(cdimy - j) * rstridey + (cdimz - k) * rstridez;
                    rxcomp[tindex] = a00;
                    rycomp[tindex] = a01;
                    rzcomp[tindex] = a02;
                }
            }

    // Transform into frequency domain.
    fft.ForwardRealDataNoScale(xcomp, 2 * cdimx, cdimy, cdimz,
            cdimx, cdimy, cdimz, cstridey, cstridez);
    fft.ForwardRealDataNoScale(ycomp, 2 * cdimx, cdimy, cdimz,
            cdimx, cdimy, cdimz, cstridey, cstridez);
    fft.ForwardRealDataNoScale(zcomp, 2 * cdimx, cdimy, cdimz,
            cdimx, cdimy, cdimz, cstridey, cstridez);

    // Copy results from ?comp into A??.  We only need store 1/8th
    // of the results because of symmetries.
    for (k = 0; k < adimz; k++) for (j = 0; j < adimy; j++) for (i = 0; i < adimx - 1; i++) {
                OC_UINT4m aindex = i + j * astridey + k*astridez;
                index = i + j * cstridey + k*cstridez;
                A00[aindex] = xcomp[index].real();
                A01[aindex] = ycomp[index].real();
                A02[aindex] = zcomp[index].real();
            }
    // Handle special packing for i=adimx-1 plane
    for (k = 0; k < adimz; k++) for (j = 0; j < adimy; j++) {
            Oxs_Complex temp00 = fft.RetrievePackedIndex(xcomp,
                    cdimx, cdimy, cdimz, cstridey, cstridez,
                    adimx - 1, j, k);
            Oxs_Complex temp01 = fft.RetrievePackedIndex(ycomp,
                    cdimx, cdimy, cdimz, cstridey, cstridez,
                    adimx - 1, j, k);
            Oxs_Complex temp02 = fft.RetrievePackedIndex(zcomp,
                    cdimx, cdimy, cdimz, cstridey, cstridez,
                    adimx - 1, j, k);
            OC_UINT4m aindex = (adimx - 1) + j * astridey + k*astridez;
            A00[aindex] = temp00.real();
            A01[aindex] = temp01.real();
            A02[aindex] = temp02.real();
        }

    // Repeat for A11, A12 and A22. //////////////////////////////////////
    for (index = 0; index < ctotalsize; index++) {
        xcomp[index].Set(0., 0.);
        ycomp[index].Set(0., 0.);
        zcomp[index].Set(0., 0.);
    }

    for (k = 0; k < rdimz; k++) for (j = 0; j < rdimy; j++) for (i = 0; i < rdimx; i++) {
                //  OC_REALWIDE x = dx*i;    OC_REALWIDE y = dy*j;    OC_REALWIDE z = dz*k;
                //  OC_REALWIDE a11=scale*CalculateSDA11(x,y,z,dx,dy,dz);
                //  OC_REALWIDE a12=scale*CalculateSDA12(x,y,z,dx,dy,dz);
                //  OC_REALWIDE a22=scale*CalculateSDA22(x,y,z,dx,dy,dz);
                OC_REALWIDE a11 = scale * GetTensorFromBuffer(yy, i, j, k);
                OC_REALWIDE a12 = scale * GetTensorFromBuffer(yz, i, j, k);
                OC_REALWIDE a22 = scale * GetTensorFromBuffer(zz, i, j, k);

                index = i + j * rstridey + k*rstridez;
                rxcomp[index] = a11;
                rycomp[index] = a12;
                rzcomp[index] = a22;
                if (i > 0) {
                    OC_UINT4m tindex = (2 * cdimx - i) + j * rstridey + k*rstridez;
                    rxcomp[tindex] = a11;
                    rycomp[tindex] = a12;
                    rzcomp[tindex] = a22;
                }
                if (j > 0) {
                    OC_UINT4m tindex = i + (cdimy - j) * rstridey + k*rstridez;
                    rxcomp[tindex] = a11;
                    rycomp[tindex] = -a12;
                    rzcomp[tindex] = a22;
                }
                if (k > 0) {
                    OC_UINT4m tindex = i + j * rstridey + (cdimz - k) * rstridez;
                    rxcomp[tindex] = a11;
                    rycomp[tindex] = -a12;
                    rzcomp[tindex] = a22;
                }
                if (i > 0 && j > 0) {
                    OC_UINT4m tindex = (2 * cdimx - i)+(cdimy - j) * rstridey + k*rstridez;
                    rxcomp[tindex] = a11;
                    rycomp[tindex] = -a12;
                    rzcomp[tindex] = a22;
                }
                if (i > 0 && k > 0) {
                    OC_UINT4m tindex = (2 * cdimx - i) + j * rstridey + (cdimz - k) * rstridez;
                    rxcomp[tindex] = a11;
                    rycomp[tindex] = -a12;
                    rzcomp[tindex] = a22;
                }
                if (j > 0 && k > 0) {
                    OC_UINT4m tindex = i + (cdimy - j) * rstridey + (cdimz - k) * rstridez;
                    rxcomp[tindex] = a11;
                    rycomp[tindex] = a12;
                    rzcomp[tindex] = a22;
                }
                if (i > 0 && j > 0 && k > 0) {
                    OC_UINT4m tindex = (2 * cdimx - i)+(cdimy - j) * rstridey + (cdimz - k) * rstridez;
                    rxcomp[tindex] = a11;
                    rycomp[tindex] = a12;
                    rzcomp[tindex] = a22;
                }
            }


    fft.ForwardRealDataNoScale(xcomp, 2 * cdimx, cdimy, cdimz,
            cdimx, cdimy, cdimz, cstridey, cstridez);
    fft.ForwardRealDataNoScale(ycomp, 2 * cdimx, cdimy, cdimz,
            cdimx, cdimy, cdimz, cstridey, cstridez);
    fft.ForwardRealDataNoScale(zcomp, 2 * cdimx, cdimy, cdimz,
            cdimx, cdimy, cdimz, cstridey, cstridez);

    for (k = 0; k < adimz; k++) for (j = 0; j < adimy; j++) for (i = 0; i < adimx - 1; i++) {
                OC_UINT4m aindex = i + j * astridey + k*astridez;
                index = i + j * cstridey + k*cstridez;
                A11[aindex] = xcomp[index].real();
                A12[aindex] = ycomp[index].real();
                A22[aindex] = zcomp[index].real();
            }
    for (k = 0; k < adimz; k++) for (j = 0; j < adimy; j++) { // Special packing case
            Oxs_Complex temp11 = fft.RetrievePackedIndex(xcomp,
                    cdimx, cdimy, cdimz, cstridey, cstridez,
                    adimx - 1, j, k);
            Oxs_Complex temp12 = fft.RetrievePackedIndex(ycomp,
                    cdimx, cdimy, cdimz, cstridey, cstridez,
                    adimx - 1, j, k);
            Oxs_Complex temp22 = fft.RetrievePackedIndex(zcomp,
                    cdimx, cdimy, cdimz, cstridey, cstridez,
                    adimx - 1, j, k);
            OC_UINT4m aindex = (adimx - 1) + j * astridey + k*astridez;
            A11[aindex] = temp11.real();
            A12[aindex] = temp12.real();
            A22[aindex] = temp22.real();
        }

#ifdef FUBAR
    /**/
    printf("\n Index     A00        A11        A22\n");
    for (k = 0; k < adimz; k++) for (j = 0; j < adimy; j++) for (i = 0; i < adimx; i++) {
                OC_UINT4m aindex = i + j * astridey + k*astridez;
                printf("%d %d %d  %#10.4g  %#10.4g  %#10.4g\n",
                        i, j, k, A00[aindex], A11[aindex], A22[aindex]);
            }
    printf("\n Index     A01        A02        A12\n");
    for (k = 0; k < adimz; k++) for (j = 0; j < adimy; j++) for (i = 0; i < adimx; i++) {
                OC_UINT4m aindex = i + j * astridey + k*astridez;
                printf("%d %d %d  %#10.4g  %#10.4g  %#10.4g\n",
                        i, j, k, A01[aindex], A02[aindex], A12[aindex]);
            }
    /**/
#endif // FUBAR
}

void PBC_Demag_2D::GetEnergy
(const Oxs_SimState& state,
        Oxs_EnergyData& oed
        ) const {
    const Oxs_RectangularMesh* mesh
            = dynamic_cast<const Oxs_RectangularMesh*> (state.mesh);

    // (Re)-initialize mesh coefficient array if mesh has changed.
    if (mesh_id != state.mesh->Id()) {
        mesh_id = 0; // Safety

        LoadPbcDemagTensor(mesh);
        if (!load_from_file_success) {
            CalculateDemagTensors(mesh);
            SavePbcDemagTensor(state.mesh);
        }
        (this->*fillcoefs)(state.mesh);
        mesh_id = state.mesh->Id();
    }

    const Oxs_MeshValue<ThreeVector>& spin = state.spin;
    const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);

    // Use supplied buffer space, and reflect that use in oed.
    oed.energy = oed.energy_buffer;
    oed.field = oed.field_buffer;
    Oxs_MeshValue<OC_REAL8m>& energy = *oed.energy_buffer;
    Oxs_MeshValue<ThreeVector>& field = *oed.field_buffer;

    // Calculate FFT of Mx, My and Mz
    OC_UINT4m i, j, k;
    OC_UINT4m mstridey = rdimx; // Assume import mesh is tight packed
    OC_UINT4m mstridez = rdimy*mstridey;
    OC_UINT4m rstridey = 2 * cstridey;
    OC_UINT4m rstridez = 2 * cstridez;
    OXS_COMPLEX_REAL_TYPE* rxcomp
            = static_cast<OXS_COMPLEX_REAL_TYPE*> (static_cast<void*> (xcomp));
    OXS_COMPLEX_REAL_TYPE* rycomp
            = static_cast<OXS_COMPLEX_REAL_TYPE*> (static_cast<void*> (ycomp));
    OXS_COMPLEX_REAL_TYPE* rzcomp
            = static_cast<OXS_COMPLEX_REAL_TYPE*> (static_cast<void*> (zcomp));

    for (k = 0; k < rdimz; k++) for (j = 0; j < rdimy; j++) {
            OC_UINT4m mindex = j * mstridey + k*mstridez;
            OC_UINT4m rindex = j * rstridey + k*rstridez;
            for (i = 0; i < rdimx; i++) {
                OC_REAL8m scale = Ms[mindex];
                const ThreeVector& vec = spin[mindex];
                ++mindex;
                rxcomp[rindex] = scale * vec.x;
                rycomp[rindex] = scale * vec.y;
                rzcomp[rindex] = scale * vec.z;
                ++rindex;
            }
        }

#if REPORT_TIME
    ffttime.Start();
#endif // REPORT_TIME
    fft.ForwardRealDataNoScale(xcomp, rdimx, rdimy, rdimz,
            cdimx, cdimy, cdimz, cstridey, cstridez);
    fft.ForwardRealDataNoScale(ycomp, rdimx, rdimy, rdimz,
            cdimx, cdimy, cdimz, cstridey, cstridez);
    fft.ForwardRealDataNoScale(zcomp, rdimx, rdimy, rdimz,
            cdimx, cdimy, cdimz, cstridey, cstridez);
#if REPORT_TIME
    ffttime.Stop();
#endif // REPORT_TIME

    // Calculate field components in frequency domain, then do iFFT to
    // transform back to space domain.  Make use of realness and even/odd
    // properties of interaction matrices Axx.
#if REPORT_TIME
    convtime.Start();
#endif // REPORT_TIME
    OC_UINT4m block;
    for (block = 0; block < 4; block++) {
        OC_UINT4m base_offset = 0;
        if (block % 2 == 1) {
            if (cdimy == 1) continue; // There is no second "j" block
            base_offset += cstridey * (cdimy / 2);
        }
        if (block / 2 == 1) {
            if (cdimz == 1) continue; // There is no second "k" block
            base_offset += cstridez * (cdimz / 2);
        }
        Oxs_Complex* xbase = xcomp + base_offset;
        Oxs_Complex* ybase = ycomp + base_offset;
        Oxs_Complex* zbase = zcomp + base_offset;
        int jsign = 1;
        if (block % 2 == 1) jsign = -1;
        int ksign = 1;
        if (block / 2 == 1) ksign = -1;
        Oxs_Complex Mx, My, Mz;
        OXS_COMPLEX_REAL_TYPE a00, a11, a22, a01, a02, a12;
        OC_UINT4m aindex = 0, cindex;
        k = 0;
        do {
            // j==0, i==0
            cindex = k*cstridez;
            Mx = xbase[cindex];
            My = ybase[cindex];
            Mz = zbase[cindex];
            if (k == 0) {
                aindex = 0;
                if (jsign != 1) aindex += (adimy - 1) * astridey;
                if (ksign != 1) aindex += (adimz - 1) * astridez;
                OXS_COMPLEX_REAL_TYPE a00b, a11b, a22b, a01b, a02b, a12b;
                a00 = A00[aindex];
                a00b = A00[aindex + cdimx];
                a11 = A11[aindex];
                a11b = A11[aindex + cdimx];
                a22 = A22[aindex];
                a22b = A22[aindex + cdimx];
                a01 = A01[aindex];
                a01b = A01[aindex + cdimx];
                a02 = A02[aindex];
                a02b = A02[aindex + cdimx];
                a12 = A12[aindex];
                a12b = A12[aindex + cdimx];
                xbase[cindex].re = a00 * Mx.re + a01 * My.re + a02 * Mz.re; // Hx
                xbase[cindex].im = a00b * Mx.im + a01b * My.im + a02b * Mz.im;
                ybase[cindex].re = a01 * Mx.re + a11 * My.re + a12 * Mz.re; // Hy
                ybase[cindex].im = a01b * Mx.im + a11b * My.im + a12b * Mz.im;
                zbase[cindex].re = a02 * Mx.re + a12 * My.re + a22 * Mz.re; // Hz
                zbase[cindex].im = a02b * Mx.im + a12b * My.im + a22b * Mz.im;
            } else {
                // k>0
                if (ksign == 1) {
                    aindex = k*astridez;
                } else {
                    aindex = (adimz - 1 - k) * astridez + cdimx;
                }
                if (jsign != 1) aindex += (adimy - 1) * astridey;
                a00 = A00[aindex];
                a11 = A11[aindex];
                a22 = A22[aindex];
                a01 = A01[aindex];
                a02 = ksign * A02[aindex];
                a12 = ksign * A12[aindex];
                xbase[cindex] = a00 * Mx + a01 * My + a02*Mz; // Hx
                ybase[cindex] = a01 * Mx + a11 * My + a12*Mz; // Hy
                zbase[cindex] = a02 * Mx + a12 * My + a22*Mz; // Hz
            }
            // j==0, i>0
            if (ksign == 1) aindex = k * astridez;
            else aindex = (adimz - 1 - k) * astridez;
            if (jsign != 1) aindex += (adimy - 1) * astridey;
            for (i = 1; i < cdimx; i++) {
                ++cindex;
                Mx = xbase[cindex];
                My = ybase[cindex];
                Mz = zbase[cindex];
                ++aindex;
                a00 = A00[aindex];
                a11 = A11[aindex];
                a22 = A22[aindex];
                a01 = A01[aindex];
                a02 = ksign * A02[aindex];
                a12 = ksign * A12[aindex];
                xbase[cindex] = a00 * Mx + a01 * My + a02*Mz; // Hx
                ybase[cindex] = a01 * Mx + a11 * My + a12*Mz; // Hy
                zbase[cindex] = a02 * Mx + a12 * My + a22*Mz; // Hz
            }

            for (j = 1; j < cdimy / 2; j++) {
                // i==0
                cindex = j * cstridey + k*cstridez;
                Mx = xbase[cindex];
                My = ybase[cindex];
                Mz = zbase[cindex];
                if (ksign == 1) aindex = k * astridez;
                else aindex = (adimz - 1 - k) * astridez;
                if (jsign == 1) aindex += j * astridey;
                else aindex += (adimy - 1 - j) * astridey + cdimx;
                a00 = A00[aindex];
                a11 = A11[aindex];
                a22 = A22[aindex];
                a01 = jsign * A01[aindex];
                a02 = ksign * A02[aindex];
                a12 = jsign * ksign * A12[aindex];
                xbase[cindex] = a00 * Mx + a01 * My + a02*Mz; // Hx
                ybase[cindex] = a01 * Mx + a11 * My + a12*Mz; // Hy
                zbase[cindex] = a02 * Mx + a12 * My + a22*Mz; // Hz
                if (jsign != 1) aindex -= cdimx;
                for (i = 1; i < cdimx; i++) {
                    ++cindex;
                    Mx = xbase[cindex];
                    My = ybase[cindex];
                    Mz = zbase[cindex];
                    ++aindex;
                    a00 = A00[aindex];
                    a11 = A11[aindex];
                    a22 = A22[aindex];
                    a01 = jsign * A01[aindex];
                    a02 = ksign * A02[aindex];
                    a12 = jsign * ksign * A12[aindex];
                    xbase[cindex] = a00 * Mx + a01 * My + a02*Mz; // Hx
                    ybase[cindex] = a01 * Mx + a11 * My + a12*Mz; // Hy
                    zbase[cindex] = a02 * Mx + a12 * My + a22*Mz; // Hz
                }
            }
        } while (++k < cdimz / 2); // For each block, do one pass through
        /// the above loop if cdimz is 1 or 2.  If cdimz is 1 then
        /// ksign will always be 1, because the upper half "blocks"
        /// are skipped.
    }
#if REPORT_TIME
    convtime.Stop();
#endif // REPORT_TIME

#if REPORT_TIME
    ffttime.Start();
#endif // REPORT_TIME
    fft.InverseRealDataNoScale(xcomp, rdimx, rdimy, rdimz,
            cdimx, cdimy, cdimz, cstridey, cstridez);
    fft.InverseRealDataNoScale(ycomp, rdimx, rdimy, rdimz,
            cdimx, cdimy, cdimz, cstridey, cstridez);
    fft.InverseRealDataNoScale(zcomp, rdimx, rdimy, rdimz,
            cdimx, cdimy, cdimz, cstridey, cstridez);
#if REPORT_TIME
    ffttime.Stop();
#endif // REPORT_TIME

#if REPORT_TIME
    dottime.Start();
#endif // REPORT_TIME
    OXS_COMPLEX_REAL_TYPE mult = -0.5 * MU0;
    for (k = 0; k < rdimz; k++) for (j = 0; j < rdimy; j++) {
            OC_UINT4m mindex = j * mstridey + k*mstridez;
            OC_UINT4m rindex = j * rstridey + k*rstridez;
            for (i = 0; i < rdimx; i++) {
                field[mindex].x = rxcomp[rindex];
                field[mindex].y = rycomp[rindex];
                field[mindex].z = rzcomp[rindex];
                // Calculate also pointwise energy density: -0.5*MU0*<M,H>
                OXS_COMPLEX_REAL_TYPE dot = spin[mindex] * field[mindex];
                energy[mindex] = mult * Ms[mindex] * dot;
                ++mindex;
                ++rindex;
            }
        }
#if REPORT_TIME
    dottime.Stop();
#endif // REPORT_TIME
}

void PBC_Demag_2D::LoadPbcDemagTensor(
        const Oxs_RectangularMesh* mesh
        ) const {

    xdim = mesh->DimX();
    ydim = mesh->DimY();
    zdim = mesh->DimZ();

    Npbc_diag.AdjustSize(mesh);
    Npbc_offdiag.AdjustSize(mesh);

    load_from_file_success = 0;
    if (tensor_file_name.length() > 0) {
       
        String diagname = tensor_file_name;
        diagname += String("-diag.ovf");
        String offdiagname = tensor_file_name;
        offdiagname += String("-offdiag.ovf");

        Vf_FileInput* vffi = NULL;
        Vf_Mesh* file_mesh = NULL;

        // Ncorr_diag
        Nb_DString dummy;
        vffi = Vf_FileInput::NewReader(diagname.c_str(), &dummy);
        if (vffi != NULL) {
            file_mesh = vffi->NewMesh();
            delete vffi;
            if (file_mesh == NULL || !mesh->IsCompatible(file_mesh)) {
                if (file_mesh != NULL) delete file_mesh;
            } else {
                mesh->FillMeshValueExact(file_mesh, Npbc_diag);
                delete file_mesh;

                // Ncorr_offdiag
                vffi = Vf_FileInput::NewReader(offdiagname.c_str(), &dummy);
                if (vffi != NULL) {
                    file_mesh = vffi->NewMesh();
                    delete vffi;
                    if (file_mesh == NULL || !mesh->IsCompatible(file_mesh)) {
                        if (file_mesh != NULL) delete file_mesh;
                    } else {
                        mesh->FillMeshValueExact(file_mesh, Npbc_offdiag);
                        delete file_mesh;
                        load_from_file_success = 1;
                    }
                }
            }
        }
    }

}


void PBC_Demag_2D::SavePbcDemagTensor(
        const Oxs_Mesh * mesh
        ) const {

    if (tensor_file_name.length() > 0 && !load_from_file_success) {
        // Oxs_MeshValue<ThreeVector> Ncorr_diag;
        // Oxs_MeshValue<ThreeVector> Ncorr_offdiag;
        String diagname = tensor_file_name;
        diagname += String("-diag.ovf");
        String offdiagname = tensor_file_name;
        offdiagname += String("-offdiag.ovf");


        mesh->WriteOvf(diagname.c_str(), 1,
                "N-diag",
                "PBC_Demag_2D::DemagTensors:"
                " Nxx, Nyy, Nzz",
                "1", "rectangular", "binary", "8", &Npbc_diag, NULL);
        mesh->WriteOvf(offdiagname.c_str(), 1,
                "N-offdiag",
                "PBC_Demag_2D::DemagTensors:"
                " Nxy, Nxz, Nyz",
                "1", "rectangular", "binary", "8", &Npbc_offdiag, NULL);
    }

}

void PBC_Demag_2D::CalculateDemagTensors(
        const Oxs_RectangularMesh* mesh
        ) const {

    ReleaseMemory();

    Npbc_diag.AdjustSize(mesh);
    Npbc_offdiag.AdjustSize(mesh);

    OC_REALWIDE dx = mesh->EdgeLengthX();
    OC_REALWIDE dy = mesh->EdgeLengthY();
    OC_REALWIDE dz = mesh->EdgeLengthZ();

    xdim = mesh->DimX();
    ydim = mesh->DimY();
    zdim = mesh->DimZ();
    OC_UINT4m xydim = xdim*ydim;


    OC_REALWIDE maxedge = dx;
    if (dy > maxedge) maxedge = dy;
    if (dz > maxedge) maxedge = dz;
    dx /= maxedge;
    dy /= maxedge;
    dz /= maxedge;


    OC_REALWIDE x, y, z;

    int gxx, gyy, gzz;


    if (sample_repeat_n > 0) {
        gxx = sample_repeat_n;
        gyy = sample_repeat_n;
        gzz = sample_repeat_n;
    } else {
        gxx = FindG(xx, dx * dy*dz, xdim*dx, ydim * dy);
        gyy = FindG(yy, dx * dy*dz, xdim*dx, ydim * dy);
        gzz = FindG(zz, dx * dy*dz, xdim*dx, ydim * dy);
    }
    //printf("gxx=%d  gyy=%d  gzz=%d\n", gxx, gyy, gzz);

    OC_UINT4m index, i, j, k;

    for (k = 0; k < zdim; k++) {
        z = k*dz;
        for (j = 0; j < ydim; j++) {
            y = j*dy;
            for (i = 0; i < xdim; i++) {
                x = i*dx;
                index = k * xydim + j * xdim + i;
                Npbc_diag[index].x = CalculateSingleTensor(xx, gxx, x, y, z, dx, dy, dz);
                Npbc_diag[index].y = CalculateSingleTensor(yy, gyy, x, y, z, dx, dy, dz);
                //Npbc_diag[index].z=CalculateSingleTensor(zz,gzz,x,y,z,dx,dy,dz);
                Npbc_diag[index].z = -(Npbc_diag[index].x + Npbc_diag[index].y);
                Npbc_offdiag[index].x = CalculateSingleTensor(xy, gxx, x, y, z, dx, dy, dz);
                Npbc_offdiag[index].y = CalculateSingleTensor(xz, gxx, x, y, z, dx, dy, dz);
                Npbc_offdiag[index].z = CalculateSingleTensor(yz, gxx, x, y, z, dx, dy, dz);
            }
        }
    }

    Npbc_diag[0].z += 1.0;

}

int PBC_Demag_2D::FindG(
        enum TensorComponent comp,
        OC_REALWIDE v, OC_REALWIDE Tx, OC_REALWIDE Ty
        ) const {

    OC_REALWIDE tmp;

    switch (comp) {
        case xy:
        case xz:
        case yz:
        case xx:
            tmp = v / (4 * PI * (Tx * Tx) * sqrt(Tx * Tx + Ty * Ty) * pbc_2d_error);
            return (int) pow(tmp, 1 / 3.0) + 1;
        case yy:
            tmp = v / (4 * PI * (Ty * Ty) * sqrt(Tx * Tx + Ty * Ty) * pbc_2d_error);
            return (int) pow(tmp, 1 / 3.0) + 1;
        case zz:
            tmp = v * sqrt(Tx * Tx + Ty * Ty) / (4 * PI * (Tx * Ty * Tx * Ty) * pbc_2d_error);
            return (int) pow(tmp, 1 / 3.0) + 1;
    }
}

OC_REALWIDE PBC_Demag_2D::CalculateSingleTensor(
        enum TensorComponent comp, int g, OC_REALWIDE x, OC_REALWIDE y, OC_REALWIDE z,
        OC_REALWIDE a, OC_REALWIDE b, OC_REALWIDE c
        ) const {

    if ((comp == xy || comp == xz || comp == yz) && x * y == 0) return 0.0;

    OC_REALWIDE Tx = xdim*a, Ty = ydim*b, cof1 = 1 / (4 * PI * a * b * c), cof2 = a * b * c / (4 * PI);
    OC_REALWIDE* tmpx = new OC_REALWIDE[2 * g + 2];
    OC_REALWIDE* tmpy = new OC_REALWIDE[2 * g + 1];
    OC_REALWIDE tpx, tpy, radius_sq;
    for (int i = -g; i <= g; i++) {
        for (int j = -g; j <= g; j++) {
            tpx = x + i*Tx;
            tpy = y + j*Ty;
            radius_sq = tpx * tpx + tpy * tpy + z*z;
            if (radius_sq < asymptotic_radius_sq) {
                tmpy[j + g] = DemagTensorNormal(comp, tpx, tpy, z, a, b, c) * cof1;
            } else if (radius_sq < dipolar_radius_sq) {
                tmpy[j + g] = DemagTensorAsymptotic(comp, tpx, tpy, z, a, b, c);
            } else {
                printf("%f\n", radius_sq);
                tmpy[j + g] = DemagTensorDipolar(comp, tpx, tpy, z) * cof2;
            }
        }
        tmpx[i + g] = AccurateSum(2 * g + 1, tmpy);
    }

    OC_REALWIDE X0 = (g + 0.5) * Tx;
    OC_REALWIDE Y0 = (g + 0.5) * Ty;

    tmpx[2 * g + 1] = DemagTensorInfinite(comp, x, y, z, X0, Y0) * cof2 / (Tx * Ty);

    OC_REALWIDE result = AccurateSum(2 * g + 2, tmpx);

    delete[] tmpx;
    delete[] tmpy;

    return result;
}

OC_REAL8m PBC_Demag_2D::GetTensorFromBuffer(
        enum TensorComponent comp, int i, int j, int k
        ) const {

    int index = xdim * ydim * k + xdim * j + i;

    switch (comp) {
        case xx:
            return Npbc_diag[index].x;
        case yy:
            return Npbc_diag[index].y;
        case zz:
            return Npbc_diag[index].z;
        case xy:
            return Npbc_offdiag[index].x;
        case xz:
            return Npbc_offdiag[index].y;
        case yz:
            return Npbc_offdiag[index].z;
    }

}




