/* FILE: coarsedemag.cc            -*-Mode: c++-*-
 *
 * Average H demag field across rectangular cells.  Simple
 * implementation that does not take advantage of any memory
 * saving or calculation speed improvements possible using
 * symmetries in interaction tensor or realness of data in
 * spatial domain.
 *
 * Calculation is actually done on a coarser grid than that
 * defined by the mesh.
 *
 */

#define STANDARD

#include <string>

#include "oc.h"
#include "nb.h"
#include "director.h"
#include "key.h"
#include "mesh.h"
#include "meshvalue.h"
#include "simstate.h"
#include "threevector.h"
#include "energy.h"		// Needed to make MSVC++ 5 happy

#include "rectangularmesh.h"
#include "coarsedemag.h"
#include "demagcoef.h"
#include "fft.h"

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(Oxs_CoarseDemag);

/* End includes */

// Helper function
OC_INDEX Oxs_CoarseDemag::NextPowerOfTwo(OC_INDEX n) const
{ // Returns first power of two >= n
  OC_INDEX m=1;
  while(m<n) {
    m*=2;
    if(m<1) {
      char msgbuf[1024];
      Oc_Snprintf(msgbuf,sizeof(msgbuf),
                  ": Import n=%ld too big",static_cast<long int>(n));
      String msg =
        String("OC_INDEX overflow in Oxs_CoarseDemag::NextPowerOfTwo")
        + String(msgbuf);
      OXS_THROW(Oxs_BadParameter,msg);
    }
  }
  return m;
}

// Constructor
Oxs_CoarseDemag::Oxs_CoarseDemag(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Oxs_Energy(name,newdtr,argstr),
    local_correction_neighborhood(0),halfNx(0.),halfNy(0.),halfNz(0.),
    xdim(0),ydim(0),zdim(0),totalsize(0),
    pxdim(0),pydim(0),pzdim(0),ptotalsize(0),mesh_id(0),
    A00(NULL),A01(NULL),A02(NULL),A11(NULL),A12(NULL),A22(NULL),
    Mx(NULL),My(NULL),Mz(NULL),Hcomp(NULL)
{
  local_correction_neighborhood = GetUIntInitValue("corrsize",0);
  if(local_correction_neighborhood>1) {
    String msg=String("Error in Oxs_CoarseDemag initialization; "
		      "corrsize can only be 0 or 1");
    throw Oxs_ExtError(msg.c_str());
  }
  subx = GetIntInitValue("subx",1);
  suby = GetIntInitValue("suby",1);
  subz = GetIntInitValue("subz",1);
  VerifyAllInitArgsUsed();
}


OC_BOOL Oxs_CoarseDemag::Init()
{
  mesh_id = 0;
  ReleaseMemory();
  return Oxs_Energy::Init();
}

void Oxs_CoarseDemag::ReleaseMemory()
{
  if(Hcomp!=NULL) { delete[] Hcomp; Hcomp=NULL; }
  if(Mx!=NULL) { delete[] Mx; Mx=NULL; }
  if(My!=NULL) { delete[] My; My=NULL; }
  if(Mz!=NULL) { delete[] Mz; Mz=NULL; }
  if(A00!=NULL) { delete[] A00; A00=NULL; }
  if(A01!=NULL) { delete[] A01; A01=NULL; }
  if(A02!=NULL) { delete[] A02; A02=NULL; }
  if(A11!=NULL) { delete[] A11; A11=NULL; }
  if(A12!=NULL) { delete[] A12; A12=NULL; }
  if(A22!=NULL) { delete[] A22; A22=NULL; }
  xdim=ydim=zdim=totalsize=0;
  pxdim=pydim=pzdim=ptotalsize=0;
}

void Oxs_CoarseDemag::FillCoefficientArrays(const Oxs_Mesh* genmesh) const
{ // This routine is conceptually const.
  const Oxs_RectangularMesh* mesh
    = dynamic_cast<const Oxs_RectangularMesh*>(genmesh);
  if(mesh==NULL) {
    String msg=String("Object ")
      + String(genmesh->InstanceName())
      + String(" is not a rectangular mesh.");
    throw Oxs_ExtError(msg.c_str());
  }

  // Fill dimension variables
  xdim = mesh->DimX();
  ydim = mesh->DimY();
  zdim = mesh->DimZ();
  totalsize=xdim*ydim*zdim;
  if(xdim==0 || ydim==0 || zdim==0) return; // Empty mesh!
  if(totalsize < xdim || totalsize < ydim || totalsize < zdim) {
    // Partial overflow check
    String msg = String("OC_INDEX overflow in ") + String(InstanceName())
      + (": Product xdim*ydim*zdim too big to fit in a OC_INDEX variable");
    throw Oxs_ExtError(msg.c_str());
  }

  // Dimensions of coarse grid.
  OC_INDEX cxdim = 1 + (xdim - 1) / subx;
  OC_INDEX cydim = 1 + (ydim - 1) / suby;
  OC_INDEX czdim = 1 + (zdim - 1) / subz;

  pxdim = NextPowerOfTwo(2*cxdim);
  pydim = NextPowerOfTwo(2*cydim);
  pzdim = NextPowerOfTwo(2*czdim);
  ptotalsize=pxdim*pydim*pzdim;
  if(ptotalsize < pxdim || ptotalsize < pydim || ptotalsize < pzdim) {
    // Partial overflow check
    String msg = String("OC_INDEX overflow in ") + String(InstanceName())
      + String(": Product pxdim*pydim*pzdim too big"
               " to fit in a OC_INDEX variable");
    throw Oxs_ExtError(msg.c_str());
  }

  // Allocate memory for interaction matrices and magnetization components
  A00   = new Oxs_Complex[ptotalsize];
  A01   = new Oxs_Complex[ptotalsize];
  A02   = new Oxs_Complex[ptotalsize];
  A11   = new Oxs_Complex[ptotalsize];
  A12   = new Oxs_Complex[ptotalsize];
  A22   = new Oxs_Complex[ptotalsize];
  Mx    = new Oxs_Complex[ptotalsize];
  My    = new Oxs_Complex[ptotalsize];
  Mz    = new Oxs_Complex[ptotalsize];
  Hcomp = new Oxs_Complex[ptotalsize];
  if(Hcomp==NULL || Mx==NULL || My==NULL || Mz==NULL
     || A00==NULL || A01==NULL || A02==NULL
     || A11==NULL || A12==NULL || A22==NULL) {
    // Safety check for those machines on which new[] doesn't throw
    // BadAlloc.
    String msg = String("Insufficient memory in simpledemag constructor.");
    throw Oxs_ExtError(msg.c_str());
  }

  // Initialize interaction matrices to zero
  OC_INDEX pindex;
  for(pindex=0;pindex<ptotalsize;pindex++) A00[pindex].Set(0.,0.);
  for(pindex=0;pindex<ptotalsize;pindex++) A01[pindex].Set(0.,0.);
  for(pindex=0;pindex<ptotalsize;pindex++) A02[pindex].Set(0.,0.);
  for(pindex=0;pindex<ptotalsize;pindex++) A11[pindex].Set(0.,0.);
  for(pindex=0;pindex<ptotalsize;pindex++) A12[pindex].Set(0.,0.);
  for(pindex=0;pindex<ptotalsize;pindex++) A22[pindex].Set(0.,0.);

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
  OC_INDEX i,j,k;
  OC_INDEX pxydim=pxdim*pydim;

  // NOTE: dx, dy, dz of the COARSE mesh are calculated.
  OC_REALWIDE dx = mesh->EdgeLengthX();
  OC_REALWIDE dy = mesh->EdgeLengthY();
  OC_REALWIDE dz = mesh->EdgeLengthZ();
  dx *= subx;  dy *= suby;  dz *= subz;

  // For demag calculation, all that matters is the relative
  // sizes of dx, dy and dz.  To help insure we don't run
  // outside floating point range, rescale these values so
  // largest is 1.0
  OC_REALWIDE maxedge=dx;
  if(dy>maxedge) maxedge=dy;
  if(dz>maxedge) maxedge=dz;
  dx/=maxedge; dy/=maxedge; dz/=maxedge;
  OC_REALWIDE scale = -1./(4*WIDE_PI*dx*dy*dz);
  for(k=0;k<czdim;k++) for(j=0;j<cydim;j++) for(i=0;i<cxdim;i++) {
    OC_REALWIDE x = dx*i;
    OC_REALWIDE y = dy*j;
    OC_REALWIDE z = dz*k;
    OC_REALWIDE coef=scale*CalculateSDA00(x,y,z,dx,dy,dz);
    A00[i+j*pxdim+k*pxydim].Set(coef,0.);
    if(i>0) A00[(pxdim-i)+j*pxdim+k*pxydim].Set(coef,0.);
    if(j>0) A00[i+(pydim-j)*pxdim+k*pxydim].Set(coef,0.);
    if(k>0) A00[i+j*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
    if(i>0 && j>0)
      A00[(pxdim-i)+(pydim-j)*pxdim+k*pxydim].Set(coef,0.);
    if(i>0 && k>0)
      A00[(pxdim-i)+j*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
    if(j>0 && k>0)
      A00[i+(pydim-j)*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
    if(i>0 && j>0 && k>0)
      A00[(pxdim-i)+(pydim-j)*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
  }
  for(k=0;k<czdim;k++) for(j=0;j<cydim;j++) for(i=0;i<cxdim;i++) {
    OC_REALWIDE x = dx*i;
    OC_REALWIDE y = dy*j;
    OC_REALWIDE z = dz*k;
    OC_REALWIDE coef=scale*CalculateSDA11(x,y,z,dx,dy,dz);
    A11[i+j*pxdim+k*pxydim].Set(coef,0.);
    if(i>0) A11[(pxdim-i)+j*pxdim+k*pxydim].Set(coef,0.);
    if(j>0) A11[i+(pydim-j)*pxdim+k*pxydim].Set(coef,0.);
    if(k>0) A11[i+j*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
    if(i>0 && j>0)
      A11[(pxdim-i)+(pydim-j)*pxdim+k*pxydim].Set(coef,0.);
    if(i>0 && k>0)
      A11[(pxdim-i)+j*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
    if(j>0 && k>0)
      A11[i+(pydim-j)*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
    if(i>0 && j>0 && k>0)
      A11[(pxdim-i)+(pydim-j)*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
  }
  for(k=0;k<czdim;k++) for(j=0;j<cydim;j++) for(i=0;i<cxdim;i++) {
    OC_REALWIDE x = dx*i;
    OC_REALWIDE y = dy*j;
    OC_REALWIDE z = dz*k;
    OC_REALWIDE coef=scale*CalculateSDA22(x,y,z,dx,dy,dz);
    A22[i+j*pxdim+k*pxydim].Set(coef,0.);
    if(i>0) A22[(pxdim-i)+j*pxdim+k*pxydim].Set(coef,0.);
    if(j>0) A22[i+(pydim-j)*pxdim+k*pxydim].Set(coef,0.);
    if(k>0) A22[i+j*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
    if(i>0 && j>0)
      A22[(pxdim-i)+(pydim-j)*pxdim+k*pxydim].Set(coef,0.);
    if(i>0 && k>0)
      A22[(pxdim-i)+j*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
    if(j>0 && k>0)
      A22[i+(pydim-j)*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
    if(i>0 && j>0 && k>0)
      A22[(pxdim-i)+(pydim-j)*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
  }
  for(k=0;k<czdim;k++) for(j=0;j<cydim;j++) for(i=0;i<cxdim;i++) {
    OC_REALWIDE x = dx*i;
    OC_REALWIDE y = dy*j;
    OC_REALWIDE z = dz*k;
    OC_REALWIDE coef=scale*CalculateSDA01(x,y,z,dx,dy,dz);
    A01[i+j*pxdim+k*pxydim].Set(coef,0.);
    if(i>0) A01[(pxdim-i)+j*pxdim+k*pxydim].Set(-coef,0.);
    if(j>0) A01[i+(pydim-j)*pxdim+k*pxydim].Set(-coef,0.);
    if(k>0) A01[i+j*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
    if(i>0 && j>0)
      A01[(pxdim-i)+(pydim-j)*pxdim+k*pxydim].Set(coef,0.);
    if(i>0 && k>0)
      A01[(pxdim-i)+j*pxdim+(pzdim-k)*pxydim].Set(-coef,0.);
    if(j>0 && k>0)
      A01[i+(pydim-j)*pxdim+(pzdim-k)*pxydim].Set(-coef,0.);
    if(i>0 && j>0 && k>0)
      A01[(pxdim-i)+(pydim-j)*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
  }
  for(k=0;k<czdim;k++) for(j=0;j<cydim;j++) for(i=0;i<cxdim;i++) {
    OC_REALWIDE x = dx*i;
    OC_REALWIDE y = dy*j;
    OC_REALWIDE z = dz*k;
    OC_REALWIDE coef=scale*CalculateSDA02(x,y,z,dx,dy,dz);
    A02[i+j*pxdim+k*pxydim].Set(coef,0.);
    if(i>0) A02[(pxdim-i)+j*pxdim+k*pxydim].Set(-coef,0.);
    if(j>0) A02[i+(pydim-j)*pxdim+k*pxydim].Set(coef,0.);
    if(k>0) A02[i+j*pxdim+(pzdim-k)*pxydim].Set(-coef,0.);
    if(i>0 && j>0)
      A02[(pxdim-i)+(pydim-j)*pxdim+k*pxydim].Set(-coef,0.);
    if(i>0 && k>0)
      A02[(pxdim-i)+j*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
    if(j>0 && k>0)
      A02[i+(pydim-j)*pxdim+(pzdim-k)*pxydim].Set(-coef,0.);
    if(i>0 && j>0 && k>0)
      A02[(pxdim-i)+(pydim-j)*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
  }
  for(k=0;k<czdim;k++) for(j=0;j<cydim;j++) for(i=0;i<cxdim;i++) {
    OC_REALWIDE x = dx*i;
    OC_REALWIDE y = dy*j;
    OC_REALWIDE z = dz*k;
    OC_REALWIDE coef=scale*CalculateSDA12(x,y,z,dx,dy,dz);
    A12[i+j*pxdim+k*pxydim].Set(coef,0.);
    if(i>0) A12[(pxdim-i)+j*pxdim+k*pxydim].Set(coef,0.);
    if(j>0) A12[i+(pydim-j)*pxdim+k*pxydim].Set(-coef,0.);
    if(k>0) A12[i+j*pxdim+(pzdim-k)*pxydim].Set(-coef,0.);
    if(i>0 && j>0)
      A12[(pxdim-i)+(pydim-j)*pxdim+k*pxydim].Set(-coef,0.);
    if(i>0 && k>0)
      A12[(pxdim-i)+j*pxdim+(pzdim-k)*pxydim].Set(-coef,0.);
    if(j>0 && k>0)
      A12[i+(pydim-j)*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
    if(i>0 && j>0 && k>0)
      A12[(pxdim-i)+(pydim-j)*pxdim+(pzdim-k)*pxydim].Set(coef,0.);
  }

  // Setup local corrections, if requested
#ifdef STANDARD
  if(local_correction_neighborhood>0) {
    // Determine self demag coefficients for fine mesh.
    // These are used directly in GetEnergy() to
    // make local corrections at the fine mesh level.
    halfNx=SelfDemagNx(dx,dy,dz)/2.;
    halfNy=SelfDemagNy(dx,dy,dz)/2.;
    halfNz=SelfDemagNz(dx,dy,dz)/2.;
#ifndef FOO
    // Remove local contribution from coarse mesh factors.
    OC_REAL8m blockwgt;
    blockwgt = 1/OC_REAL8m(subx);
    A00[0].re += 2*halfNx*blockwgt;
    if(pxdim>1) { // What is proper way to handle cdim<2 cases?
      A00[1].re -= halfNx*blockwgt;
      A00[pxdim-1].re -= halfNx*blockwgt;
    }
    blockwgt = 1/OC_REAL8m(suby);
    A11[0].re += 2*halfNy*blockwgt;
    if(pydim>1) {
      A11[pxdim].re -= halfNy*blockwgt;
      A11[(pydim-1)*pxdim].re -= halfNy*blockwgt;
    }
    blockwgt = 1/OC_REAL8m(subz);
    A22[0].re += 2*halfNz*blockwgt;
    if(pzdim>1) {
      A22[pxydim].re -= halfNz*blockwgt;
      A22[(pzdim-1)*pxydim].re -= halfNz*blockwgt;
    }
#endif // FOO
  }
#else // STANDARD
  // Self correction term only
  if(local_correction_neighborhood>0) {
    // Determine self demag coefficients for fine mesh.
    // These are used directly in GetEnergy() to
    // make local corrections at the fine mesh level.
    halfNx=SelfDemagNx(dx,dy,dz)/2.;
    halfNy=SelfDemagNy(dx,dy,dz)/2.;
    halfNz=SelfDemagNz(dx,dy,dz)/2.;
    // Remove local contribution from coarse mesh factors.
    A00[0].re += 2*halfNx;
    A11[0].re += 2*halfNy;
    A22[0].re += 2*halfNz;
  }
#endif // STANDARD

  // Transform interaction matrices into frequency domain.
  fft.Forward(pxdim,pydim,pzdim,A00);
  fft.Forward(pxdim,pydim,pzdim,A01);
  fft.Forward(pxdim,pydim,pzdim,A02);
  fft.Forward(pxdim,pydim,pzdim,A11);
  fft.Forward(pxdim,pydim,pzdim,A12);
  fft.Forward(pxdim,pydim,pzdim,A22);
}

void Oxs_CoarseDemag::GetEnergy
(const Oxs_SimState& state,
 Oxs_EnergyData& oed
 ) const
{
  // (Re)-initialize mesh coefficient array if mesh has changed.
  if(mesh_id != state.mesh->Id()) {
    mesh_id = 0; // Safety
    FillCoefficientArrays(state.mesh);
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
  OC_INDEX xydim=xdim*ydim;
  OC_INDEX pxydim=pxdim*pydim;
  OC_INDEX index,pindex;
  OC_INDEX i,j,k,ci,cj,ck;
  OC_INDEX xlim, cxdim = 1 + (xdim - 1) / subx;
  OC_INDEX ylim, cydim = 1 + (ydim - 1)/ suby;
  OC_INDEX zlim, czdim = 1 + (zdim - 1)/ subz;
  for(pindex=0;pindex<ptotalsize;pindex++) Mx[pindex].Set(0.,0.);
  for(pindex=0;pindex<ptotalsize;pindex++) My[pindex].Set(0.,0.);
  for(pindex=0;pindex<ptotalsize;pindex++) Mz[pindex].Set(0.,0.);

  for(ck=0;ck<czdim;ck++) for(cj=0;cj<cydim;cj++) for(ci=0;ci<cxdim;ci++) {
    Oxs_ThreeVector vec;
    xlim = ci*subx + subx;
    ylim = cj*suby + suby;
    zlim = ck*subz + subz;
    for(k=ck*subz; k<zdim && k<zlim; k++) {
      for(j=cj*suby; j<ydim && j<ylim; j++) {
	for(i=ci*subx; i<xdim && i<xlim; i++) {
	  index  = i + j*xdim + k*xydim;
	  vec += Ms[index] * spin[index];
	}
      }
    }
    vec *= 1.0 / (subx * suby * subz);
    pindex = ci + cj*pxdim + ck*pxydim;
    Mx[pindex].Set(vec.x,0);
    My[pindex].Set(vec.y,0);
    Mz[pindex].Set(vec.z,0);
  }

  fft.Forward(pxdim,pydim,pzdim,Mx);
  fft.Forward(pxdim,pydim,pzdim,My);
  fft.Forward(pxdim,pydim,pzdim,Mz);
  
  // Calculate field components in frequency domain, then do iFFT
  // to transform back to space domain.

  // Hx component
  for(pindex=0;pindex<ptotalsize;pindex++) {
    Hcomp[pindex] = A00[pindex]*Mx[pindex]
      + A01[pindex]*My[pindex] + A02[pindex]*Mz[pindex];
  }
  fft.Inverse(pxdim,pydim,pzdim,Hcomp);
  for(k=0;k<zdim;k++) for(j=0;j<ydim;j++) for(i=0;i<xdim;i++) {
    index  = i+j*xdim+k*xydim;
    pindex = (i/subx)+(j/suby)*pxdim+(k/subz)*pxydim;
    field[index].x = Hcomp[pindex].real();
  }

  // Hy component
  for(pindex=0;pindex<ptotalsize;pindex++) {
    Hcomp[pindex] = A01[pindex]*Mx[pindex]
      + A11[pindex]*My[pindex] + A12[pindex]*Mz[pindex];
  }
  fft.Inverse(pxdim,pydim,pzdim,Hcomp);
  for(k=0;k<zdim;k++) for(j=0;j<ydim;j++) for(i=0;i<xdim;i++) {
    index  = i+j*xdim+k*xydim;
    pindex = (i/subx)+(j/suby)*pxdim+(k/subz)*pxydim;
    field[index].y = Hcomp[pindex].real();
  }
  
  // Hz component
  for(pindex=0;pindex<ptotalsize;pindex++) {
    Hcomp[pindex] = A02[pindex]*Mx[pindex]
      + A12[pindex]*My[pindex] + A22[pindex]*Mz[pindex];
  }
  fft.Inverse(pxdim,pydim,pzdim,Hcomp);
  for(k=0;k<zdim;k++) for(j=0;j<ydim;j++) for(i=0;i<xdim;i++) {
    index  = i+j*xdim+k*xydim;
    pindex = (i/subx)+(j/suby)*pxdim+(k/subz)*pxydim;
    field[index].z = Hcomp[pindex].real();
  }

  // Put in local field corrections, if requested
#ifdef STANDARD
  if(local_correction_neighborhood>0) {
    for(k=0;k<zdim;k++) for(j=0;j<ydim;j++) for(i=0;i<xdim;i++) {
      index  = i+j*xdim+k*xydim;
      // x-component adjustment
      OC_REAL8m charge = Ms[index]*spin[index].x;
      if(i==0) { // Add in adjustment from lefthand edge
	field[index].x -= halfNx*charge;
      }
      if(i+1<xdim) {
	OC_INDEX nindex = index+1;
	charge -= Ms[nindex]*spin[nindex].x;
	OC_REAL8m fieldadj = halfNx*charge;
	field[index].x -= fieldadj;
	field[nindex].x += fieldadj;
      } else { // Insert phantom neighbor with 0 moment
	field[index].x -= halfNx*charge;
      }
      // y-component adjustment
      charge = Ms[index]*spin[index].y;
      if(j==0) { // Add in adjustment from bottom edge
	field[index].y -= halfNy*charge;
      }
      if(j+1<ydim) {
	OC_INDEX nindex = index+xdim;
	charge -= Ms[nindex]*spin[nindex].y;
	OC_REAL8m fieldadj = halfNy*charge;
	field[index].y -= fieldadj;
	field[nindex].y += fieldadj;
      } else { // Insert phantom neighbor with 0 moment
	field[index].y -= halfNy*charge;
      }
      // z-component adjustment
      charge = Ms[index]*spin[index].z;
      if(k==0) { // Add in adjustment from back edge
	field[index].z -= halfNz*charge;
      }
      if(k+1<zdim) {
	OC_INDEX nindex = index+xydim;
	charge -= Ms[nindex]*spin[nindex].z;
	OC_REAL8m fieldadj = halfNz*charge;
	field[index].z -= fieldadj;
	field[nindex].z += fieldadj;
      } else { // Insert phantom neighbor with 0 moment
	field[index].z -= halfNz*charge;
      }
    }
  }
#else // STANDARD
  // Self-correction only
  if(local_correction_neighborhood>0) {
    for(k=0;k<zdim;k++) for(j=0;j<ydim;j++) for(i=0;i<xdim;i++) {
      index  = i+j*xdim+k*xydim;
      OC_REAL8m wgt = -2 * Ms[index];
      OC_REAL8m fieldadj;
      fieldadj = halfNx * wgt * spin[index].x;
      field[index].x += fieldadj;
      fieldadj = halfNy * wgt * spin[index].y;
      field[index].y += fieldadj;
      fieldadj = halfNz * wgt * spin[index].z;
      field[index].z += fieldadj;
    }
  }
#endif // STANDARD

  // Calculate energy from field
  const OC_REAL8m mult = -0.5 * MU0;
  for(index=0;index<totalsize;index++) {
    energy[index] = mult * Ms[index] * ( spin[index] * field[index] );
  }

}
