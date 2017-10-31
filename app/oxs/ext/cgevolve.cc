/* FILE: cgevolve.cc                 -*-Mode: c++-*-
 *
 * Concrete minimization evolver class, using conjugate gradients
 *
 */

#include <assert.h>
#include <float.h>

#include <string>

#include "nb.h"
#include "director.h"
#include "mindriver.h"
#include "simstate.h"
#include "cgevolve.h"
#include "key.h"
#include "threevector.h"
#include "oxswarn.h"
#include "energy.h"		// Needed to make MSVC++ 5 happy

#if OOMMF_THREADS
# include <vector>
# include "oxsthread.h"
#endif // OOMMF_THREADS

#if OC_USE_SSE
# include <emmintrin.h>
#endif

#if OC_USE_NUMA
# include <numaif.h>
#endif

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(Oxs_CGEvolve);

/* End includes */

// Revision information, set via CVS keyword substitution
static const Oxs_WarningMessageRevisionInfo revision_info
  (__FILE__,
   "$Revision: 1.108 $",
   "$Date: 2010-07-20 00:54:36 $",
   "$Author: donahue $",
   "Michael J. Donahue (michael.donahue@nist.gov)");

// #define INSTRUMENT

// Constructor
Oxs_CGEvolve::Oxs_CGEvolve(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Oxs_MinEvolver(name,newdtr,argstr),
    step_attempt_count(0),
    energy_calc_count(0),
    cycle_count(0),
    cycle_sub_count(0),
    bracket_count(0),
    line_minimum_count(0),
    temp_id(0)
{
  // Check code assumptions, if any
#if OC_USE_SSE
  // Some of the SSE code assumes that Oxs_ThreeVector arrays
  // are tightly packed, with a particular order on the components.
  {
    ThreeVector foo[10];
    if((char*)(foo+9) - (char*)(foo+0) != 9*3*sizeof(OC_REAL8m)) {
      throw Oxs_ExtError(this,"ThreeVectors are not tightly packed;"
                         " this will break some of the code using"
                         " SSE intrinsics.  Edit the OOMMF platform"
                         " file to disable SSE intrinsics (see config"
                         " sse_level option) and rebuild.");
    }
    if((void*)(&(foo[0].x)) != (void*)(foo+0) ||
       &(foo[0].x) + 1 != &(foo[0].y) ||
       &(foo[0].x) + 2 != &(foo[0].z)) {
      throw Oxs_ExtError(this,"ThreeVector packing order is not"
                         " x:y:z.  This will break some of the SSE"
                         " intrinsics code.  Edit the OOMMF platform"
                         " file to disable SSE intrinsics (see config"
                         " sse_level option) and rebuild.");
    }
  }
#endif // OC_USE_SSE

  // Process arguments

  // gradient_reset_angle is input in degrees, but then converted to cos
  // of that value, for ease of use inside ::SetBasePoint().
  gradient_reset_angle_cos = GetRealInitValue("gradient_reset_angle",80);
  gradient_reset_angle_cos = cos(gradient_reset_angle_cos*PI/180.);

  gradient_reset_count = GetUIntInitValue("gradient_reset_count",50);

  OC_REAL8m min_ang = GetRealInitValue("minimum_bracket_step",0.05);
  min_ang *= PI/180.; // Convert from degrees to radians
  bracket.minstep = tan(min_ang);

  OC_REAL8m max_ang = GetRealInitValue("maximum_bracket_step",10);
  max_ang *= PI/180.; // Convert from degrees to radians
  bracket.maxstep = tan(max_ang);

  if(bracket.minstep<0.0 || bracket.minstep>bracket.maxstep) {
    throw Oxs_ExtError(this,"Invalid value for minimum_bracket_step"
			 " and/or maximum_bracket_step.");
  }

  bracket.angle_precision
    = GetRealInitValue("line_minimum_angle_precision",5);
  // Convert from degrees to sin(angle).  We want sin(angle)
  // instead of cos(angle) because we are interested in differences
  // from 90 degrees (use angle-sum formula on sin(90-acos(dot))).
  bracket.angle_precision = sin(bracket.angle_precision*PI/180.);

  bracket.relative_minspan
    = GetRealInitValue("line_minimum_relwidth",1);

  bracket.energy_precision = GetRealInitValue("energy_precision",1e-10);

  String method_name = GetStringInitValue("method","Fletcher-Reeves");
  if(method_name.compare("Fletcher-Reeves")==0) {
    basept.method = Basept_Data::FLETCHER_REEVES;
  } else if(method_name.compare("Polak-Ribiere")==0) {
    basept.method = Basept_Data::POLAK_RIBIERE;
  } else {
    String msg=String("Invalid conjugate-gradient method request: ")
      + method_name
      + String("\n Should be either Fletcher-Reeves or Polak-Ribiere.");
    throw Oxs_ExtError(this,msg.c_str());
  }

  // Setup output.  Note: MSVC++ 6.0 requires fully qualified
  // member function names.
  total_H_field_output.Setup(this,InstanceName(),"H","A/m",1,
              &Oxs_CGEvolve::UpdateDerivedFieldOutputs);
  mxHxm_output.Setup(this,InstanceName(),"mxHxm","A/m",1,
	      &Oxs_CGEvolve::UpdateDerivedFieldOutputs);
  total_energy_density_output.Setup(this,InstanceName(),
              "Total energy density","J/m^3",1,
              &Oxs_CGEvolve::UpdateDerivedFieldOutputs);
  max_mxHxm_output.Setup(this,InstanceName(),"Max mxHxm","A/m",0,
	      &Oxs_CGEvolve::UpdateDerivedOutputs);
  total_energy_output.Setup(this,InstanceName(),"Total energy","J",0,
              &Oxs_CGEvolve::UpdateDerivedOutputs);
  delta_E_output.Setup(this,InstanceName(),"Delta E","J",0,
              &Oxs_CGEvolve::UpdateDerivedOutputs);
  bracket_count_output.Setup(this,InstanceName(),"Bracket count","",0,
              &Oxs_CGEvolve::UpdateDerivedOutputs);
  line_min_count_output.Setup(this,InstanceName(),"Line min count","",0,
              &Oxs_CGEvolve::UpdateDerivedOutputs);
  cycle_count_output.Setup(this,InstanceName(),
              "Cycle count","",0,
              &Oxs_CGEvolve::UpdateDerivedOutputs);
  cycle_sub_count_output.Setup(this,InstanceName(),
              "Cycle sub count","",0,
              &Oxs_CGEvolve::UpdateDerivedOutputs);
  energy_calc_count_output.Setup(this,InstanceName(),
              "Energy calc count","",0,
              &Oxs_CGEvolve::UpdateDerivedOutputs);

  total_H_field_output.Register(director,-5);
  mxHxm_output.Register(director,-5);
  total_energy_density_output.Register(director,-5);
  max_mxHxm_output.Register(director,-5);
  total_energy_output.Register(director,-5);
  delta_E_output.Register(director,-5);
  bracket_count_output.Register(director,-5);
  line_min_count_output.Register(director,-5);
  cycle_count_output.Register(director,-5);
  cycle_sub_count_output.Register(director,-5);
  energy_calc_count_output.Register(director,-5);

  VerifyAllInitArgsUsed();

#if REPORT_TIME_CGDEVEL
  timer.resize(10);
  timer_counts.resize(10);
#endif
}

OC_BOOL Oxs_CGEvolve::Init()
{
#if REPORT_TIME
  Oc_TimeVal cpu,wall;

#if REPORT_TIME_CGDEVEL
  energyobjtime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"   Step energy  .....   %7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }
#endif // REPORT_TIME_CGDEVEL

  steponlytime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"   Step-only    .....   %7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

#if REPORT_TIME_CGDEVEL
  basepttime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      Base point   .....   %7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  findbrackettime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      Find bracket .....   %7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  findlinemintime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      Find line min ....   %7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  fillbrackettime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"         Fill bracket  ....   %7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  getenergyandmxHxmtime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"         Energy + mxHxm ...   %7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  for(unsigned int ti=0;ti<timer.size();++ti) {
    timer[ti].GetTimes(cpu,wall);
    if(double(wall)>0.0) {
      fprintf(stderr,"               timer %2u ...   %7.2f cpu /%7.2f wall,"
              " (%s)\n",
              ti,double(cpu),double(wall),InstanceName());
      if(timer_counts[ti].pass_count>0) {
        fprintf(stderr,"                \\---> passes = %d,"
                " bytes=%.2f MB, %.2f GB/sec\n",
                timer_counts[ti].pass_count,
                double(timer_counts[ti].bytes)/double(1024*1024),
                double(timer_counts[ti].bytes)
                /(double(1024*1024*1024)*double(wall)));
      }
    }
    timer[ti].Reset(); 
    timer_counts[ti].Reset();
  }
#endif // REPORT_TIME_CGDEVEL


  steponlytime.Reset();
#if REPORT_TIME_CGDEVEL
  energyobjtime.Reset();
  basepttime.Reset();
  findbrackettime.Reset();
  findlinemintime.Reset();
  fillbrackettime.Reset();
  getenergyandmxHxmtime.Reset();
#endif // REPORT_TIME_CGDEVEL

#endif // REPORT_TIME

  // Initialize instance variables
  step_attempt_count=0;
  energy_calc_count=0;
  cycle_count=0;
  cycle_sub_count=0;
  bracket_count=0;
  line_minimum_count=0;
  basept.Init();
  bestpt.Init();
  bracket.Init();

  scratch_energy.Release();
  scratch_field.Release();
  temp_id=0;
  temp_energy.Release();
  temp_mxHxm.Release();

  return Oxs_MinEvolver::Init();
}


Oxs_CGEvolve::~Oxs_CGEvolve()
{
#if REPORT_TIME
  Oc_TimeVal cpu,wall;

#if REPORT_TIME_CGDEVEL
  energyobjtime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"   Step energy  .....   %7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }
#endif // REPORT_TIME_CGDEVEL

  steponlytime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"   Step-only    .....   %7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

#if REPORT_TIME_CGDEVEL
  basepttime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      Base point   .....   %7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }
  findbrackettime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      Find bracket .....   %7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }
  findlinemintime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"      Find line min ....   %7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }
  fillbrackettime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"         Fill bracket  ....   %7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }
  getenergyandmxHxmtime.GetTimes(cpu,wall);
  if(double(wall)>0.0) {
    fprintf(stderr,"         Energy + mxHxm ...   %7.2f cpu /%7.2f wall,"
            " (%s)\n",
            double(cpu),double(wall),InstanceName());
  }

  for(unsigned int ti=0;ti<timer.size();++ti) {
    timer[ti].GetTimes(cpu,wall);
    if(double(wall)>0.0) {
      fprintf(stderr,"               timer %2u ...   %7.2f cpu /%7.2f wall,"
              " (%s)\n",
              ti,double(cpu),double(wall),InstanceName());
      if(timer_counts[ti].pass_count>0) {
        fprintf(stderr,"                \\---> passes = %d,"
                " bytes=%.2f MB, %.2f GB/sec\n",
                timer_counts[ti].pass_count,
                double(timer_counts[ti].bytes)/double(1024*1024),
                double(timer_counts[ti].bytes)
                /(double(1024*1024*1024)*double(wall)));
      }
    }
  }
#endif // REPORT_TIME_CGDEVEL

#endif // REPORT_TIME
}

#if OOMMF_THREADS
class _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadA : public Oxs_ThreadRunObj {
public:
  static Oxs_Mutex job_control;
  static OC_INDEX offset;

  const Oxs_MeshValue<ThreeVector>* direction;
  const Oxs_MeshValue<ThreeVector>* best_spin;
  Oxs_MeshValue<ThreeVector>* spin;

  OC_REAL8m tsq;
  OC_REAL8m dvec_scale;
  OC_INDEX vecsize;
  OC_INDEX block_size;

  _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadA()
    : direction(0), best_spin(0), spin(0),
      tsq(0.),dvec_scale(0.),vecsize(0), block_size(0) {}

  void Cmd(int threadnumber, void* data);
};

Oxs_Mutex _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadA::job_control;
OC_INDEX _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadA::offset(0);

void _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadA::Cmd(int /* threadnumber */,
                                                  void* /* data */)
{
  const ThreeVector *scratch_direction;
  const ThreeVector *scratch_best_spin;
  ThreeVector *scratch_spin;

  while(1) {
    job_control.Lock();
    OC_INDEX istart = offset;
    OC_INDEX istop = ( offset += block_size );
    job_control.Unlock();

    if(istart>=vecsize) break;
    if(istop>vecsize) istop=vecsize;

    scratch_direction = &((*direction)[istart]);
    scratch_best_spin = &((*best_spin)[istart]);
    scratch_spin      = &((*spin)[istart]);
    const OC_INDEX jsize = istop - istart;

    // We want to process scratch_spin starting at a multiple of 16
    const unsigned long ALIGNMENT = 16;  // In bytes
    const unsigned long STRIDE = 2;      // In ThreeVectors
    const unsigned long addr
      = (unsigned long)((void*)scratch_spin) % ALIGNMENT;
    OC_INDEX j_break_A = 0;
    while(addr+(j_break_A*sizeof(ThreeVector)) % ALIGNMENT != 0) {
      if(++j_break_A>=jsize) {
        j_break_A = jsize;  // Safety
        break;
      }
    }
    OC_INDEX j_break_B = jsize;
    while((j_break_B-j_break_A) % STRIDE != 0) {
      if(--j_break_B<=j_break_A) break;
    }

    OC_INDEX j;
    for(j=0;j<j_break_A;++j) {
      const ThreeVector& dvec = scratch_direction[j];
      ThreeVector temp = scratch_best_spin[j];
      OC_REAL8m dsq = dvec.MagSq();
      temp *= sqrt(1+tsq*dsq);
      temp.Accum(dvec_scale,dvec);
      temp.MakeUnit();
      Nb_NOP(&temp);          // Safety; icpc optimizer showed problems
      scratch_spin[j] = temp; // with this stanza in some circumstances.
    }
    for(; j<j_break_B; j+=STRIDE) {
      assert(STRIDE == 2);
      const ThreeVector& dvec0 = scratch_direction[j];
      const ThreeVector& dvec1 = scratch_direction[j+1];
      const ThreeVector& bs0 = scratch_best_spin[j];
      const ThreeVector& bs1 = scratch_best_spin[j+1];

#if !OC_USE_SSE

      OC_REAL8m mult0 = sqrt(1+tsq*(dvec0.x*dvec0.x
                                 +dvec0.y*dvec0.y+dvec0.z*dvec0.z));
      OC_REAL8m mult1 = sqrt(1+tsq*(dvec1.x*dvec1.x
                                 +dvec1.y*dvec1.y+dvec1.z*dvec1.z));
      ThreeVector spin0(mult0*bs0.x + dvec_scale*dvec0.x,
                        mult0*bs0.y + dvec_scale*dvec0.y,
                        mult0*bs0.z + dvec_scale*dvec0.z);
      spin0.MakeUnit();
      ThreeVector spin1(mult1*bs1.x + dvec_scale*dvec1.x,
                        mult1*bs1.y + dvec_scale*dvec1.y,
                        mult1*bs1.z + dvec_scale*dvec1.z);
      spin1.MakeUnit();

      scratch_spin[j].x   = spin0.x;
      scratch_spin[j].y   = spin0.y;
      scratch_spin[j].z   = spin0.z;
      scratch_spin[j+1].x = spin1.x;
      scratch_spin[j+1].y = spin1.y;
      scratch_spin[j+1].z = spin1.z;

#else // OC_USE_SSE

      __m128d dvecx = _mm_set_pd(dvec1.x,dvec0.x);
      __m128d dvecy = _mm_set_pd(dvec1.y,dvec0.y);
      __m128d dvecz = _mm_set_pd(dvec1.z,dvec0.z);
      __m128d mult = Oxs_ThreeVectorPairMagSq(dvecx,dvecy,dvecz);
      mult = _mm_add_pd(_mm_mul_pd(mult,_mm_set1_pd(tsq)),_mm_set1_pd(1.0));
      mult = _mm_sqrt_pd(mult);

      __m128d spx = _mm_mul_pd(mult,_mm_set_pd(bs1.x,bs0.x));
      spx = _mm_add_pd(spx,_mm_mul_pd(_mm_set1_pd(dvec_scale),dvecx));

      __m128d spy = _mm_mul_pd(mult,_mm_set_pd(bs1.y,bs0.y));
      spy = _mm_add_pd(spy,_mm_mul_pd(_mm_set1_pd(dvec_scale),dvecy));

      __m128d spz = _mm_mul_pd(mult,_mm_set_pd(bs1.z,bs0.z));
      spz = _mm_add_pd(spz,_mm_mul_pd(_mm_set1_pd(dvec_scale),dvecz));

      Oxs_ThreeVectorPairMakeUnit(spx,spy,spz);

      // Note: BIG assumptions about spin storage layout
      // These assumptions are checked in the Oxs_CGEvolve constructor.
      double* sspinbase = (double*)(scratch_spin+j);
      _mm_stream_pd(sspinbase,   _mm_unpacklo_pd(spx,spy));
      _mm_stream_pd(sspinbase+2, _mm_shuffle_pd(spz,spx,2));
      _mm_stream_pd(sspinbase+4, _mm_unpackhi_pd(spy,spz));

#endif // !OC_USE_SSE

    }

    for(;j<jsize;++j) {
      const ThreeVector& dvec = scratch_direction[j];
      ThreeVector temp = scratch_best_spin[j];
      OC_REAL8m dsq = dvec.MagSq();
      temp *= sqrt(1+tsq*dsq);
      temp.Accum(dvec_scale,dvec);
      temp.MakeUnit();
      Nb_NOP(&temp);          // Safety; icpc optimizer showed problems
      scratch_spin[j] = temp; // with this stanza in some circumstances.
    }

  }
}

class _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB : public Oxs_ThreadRunObj {
public:
  static Oxs_Mutex job_control;
  static OC_INDEX offset;
  static Nb_Xpfloat etemp;
  static Nb_Xpfloat dtemp;
  static Nb_Xpfloat stemp;

  const Oxs_Mesh* mesh;
  const Oxs_MeshValue<OC_REAL8m>* tmpenergy;
  const Oxs_MeshValue<OC_REAL8m>* bestpt_energy;
  const Oxs_MeshValue<OC_REAL8m>* Ms;

  const Oxs_MeshValue<ThreeVector>* direction;
  const Oxs_MeshValue<ThreeVector>* mxHxm;

  OC_REAL8m offset_sq;

  OC_INDEX vecsize;
  OC_INDEX block_size;

  _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB()
    : mesh(0),tmpenergy(0),bestpt_energy(0),
      direction(0), mxHxm(0),
      offset_sq(0.), vecsize(0), block_size(0) {}

  void Cmd(int threadnumber, void* data);
};

Oxs_Mutex _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::job_control;
OC_INDEX _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::offset(0);
Nb_Xpfloat _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::etemp(0.0);
Nb_Xpfloat _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::dtemp(0.0);
Nb_Xpfloat _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::stemp(0.0);


void _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::Cmd(int /* threadnumber */,
                                                  void* /* data */)
{
  Nb_Xpfloat work_etemp = 0.0;
  Nb_Xpfloat work_dtemp = 0.0;
  Nb_Xpfloat work_stemp = 0.0;

  Nb_Xpfloat work_etemp_b = 0.0;
  Nb_Xpfloat work_dtemp_b = 0.0;
  Nb_Xpfloat work_stemp_b = 0.0;

  while(1) {
    job_control.Lock();
    OC_INDEX istart = offset;
    OC_INDEX istop = ( offset += block_size );
    job_control.Unlock();

    if(istart>=vecsize) break;
    if(istop>vecsize) istop=vecsize;

    const OC_REAL8m* stenergy    = &((*tmpenergy)[istart]);
    const OC_REAL8m* sbenergy    = &((*bestpt_energy)[istart]);
    const OC_REAL8m* sMs         = &((*Ms)[istart]);
    const ThreeVector* sdir   = &((*direction)[istart]);
    const ThreeVector* smxHxm = &((*mxHxm)[istart]);

    OC_INDEX jsize = istop - istart;
    OC_INDEX j;
    for(j=0;j<jsize%2;++j) {
      OC_REAL8m vol = mesh->Volume(istart+j);

      const ThreeVector& vtemp = sdir[j];
      OC_REAL8m scale_adj = sMs[j]*vol/sqrt(1+offset_sq * vtemp.MagSq());

      work_etemp.Accum((stenergy[j] - sbenergy[j]) * vol);
      work_dtemp.Accum((smxHxm[j]*vtemp)*scale_adj);
      work_stemp.Accum(smxHxm[j].MagSq()*scale_adj*scale_adj);
    }
    for(;j<jsize;j+=2) {
      const ThreeVector& vtemp0 = sdir[j];
      const ThreeVector& vtemp1 = sdir[j+1];
#if !OC_USE_SSE
      OC_REAL8m vol0 = mesh->Volume(istart+j);
      OC_REAL8m vol1 = mesh->Volume(istart+j+1);
      OC_REAL8m scale_adj0 = sMs[j]*vol0/sqrt(1+offset_sq * vtemp0.MagSq());
      OC_REAL8m scale_adj1 = sMs[j+1]*vol1/sqrt(1+offset_sq * vtemp1.MagSq());
      work_etemp.Accum((stenergy[j] - sbenergy[j]) * vol0);
      work_etemp_b.Accum((stenergy[j+1] - sbenergy[j+1]) * vol1);
      work_dtemp.Accum((smxHxm[j]*vtemp0)*scale_adj0);
      work_stemp.Accum(smxHxm[j].MagSq()*scale_adj0*scale_adj0);
      work_dtemp_b.Accum((smxHxm[j+1]*vtemp1)*scale_adj1);
      work_stemp_b.Accum(smxHxm[j+1].MagSq()*scale_adj1*scale_adj1);
#else // OC_USE_SSE
      // SSE2 version of the above

      // Compute scale_adj
      __m128d tx   = _mm_set_pd(vtemp1.x,vtemp0.x);
      __m128d txsq = _mm_mul_pd(tx,tx);

      __m128d ty   = _mm_set_pd(vtemp1.y,vtemp0.y);
      __m128d tysq = _mm_mul_pd(ty,ty);

      __m128d tz   = _mm_set_pd(vtemp1.z,vtemp0.z);
      __m128d tzsq = _mm_mul_pd(tz,tz);

      __m128d denom = _mm_add_pd(_mm_add_pd(txsq,tysq),tzsq);

      denom = _mm_mul_pd(denom,_mm_set1_pd(offset_sq));
      denom = _mm_add_pd(denom,_mm_set1_pd(1.0));
      denom = _mm_sqrt_pd(denom);

      __m128d vol = _mm_set_pd(mesh->Volume(istart+j+1),
                               mesh->Volume(istart+j));
      __m128d num = _mm_loadu_pd(sMs+j);
      num = _mm_mul_pd(num,vol);

      __m128d scale_adj = _mm_div_pd(num,denom);

      // work_etemp.Accum((stenergy[j] - sbenergy[j]) * vol0);
      // work_etemp_b.Accum((stenergy[j+1] - sbenergy[j+1]) * vol1);
      __m128d etemp = _mm_sub_pd(_mm_loadu_pd(stenergy+j),
                                _mm_loadu_pd(sbenergy+j));
      etemp = _mm_mul_pd(etemp,vol);
      Nb_XpfloatDualAccum(work_etemp,work_etemp_b,etemp);

      // work_dtemp.Accum((smxHxm[j]*vtemp0)*scale_adj0);
      // work_dtemp_b.Accum((smxHxm[j+1]*vtemp1)*scale_adj1);
      __m128d sx   = _mm_set_pd(smxHxm[j+1].x,smxHxm[j].x);
      __m128d sy   = _mm_set_pd(smxHxm[j+1].y,smxHxm[j].y);
      __m128d sz   = _mm_set_pd(smxHxm[j+1].z,smxHxm[j].z);
      __m128d dot = _mm_add_pd(_mm_add_pd(_mm_mul_pd(sx,tx),_mm_mul_pd(sy,ty)),
                               _mm_mul_pd(sz,tz));
      dot = _mm_mul_pd(dot,scale_adj);
      Nb_XpfloatDualAccum(work_dtemp,work_dtemp_b,dot);

      // work_stemp.Accum(smxHxm[j].MagSq()*scale_adj0*scale_adj0);
      // work_stemp_b.Accum(smxHxm[j+1].MagSq()*scale_adj1*scale_adj1);
      __m128d ssq = _mm_add_pd(_mm_add_pd(_mm_mul_pd(sx,sx),_mm_mul_pd(sy,sy)),
                               _mm_mul_pd(sz,sz));
      ssq = _mm_mul_pd(ssq,_mm_mul_pd(scale_adj,scale_adj));
      Nb_XpfloatDualAccum(work_stemp,work_stemp_b,ssq);
#endif // OC_USE_SSE
    }
  }

  // Sum results into static accumulators
  work_etemp += work_etemp_b;
  work_dtemp += work_dtemp_b;
  work_stemp += work_stemp_b;
  job_control.Lock();
  etemp += work_etemp;
  dtemp += work_dtemp;
  stemp += work_stemp;
  job_control.Unlock();
}

class _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadC : public Oxs_ThreadRunObj {
public:
  static Oxs_Mutex result_mutex;            // Note these are static, so
  static Oxs_JobControl<OC_REAL8m> job_basket; // only one "set" of this class
  static Nb_Xpfloat etemp;                  // may be run at one time.

  const Oxs_Mesh* mesh;
  const Oxs_MeshValue<OC_REAL8m>* energy;
  OC_INT4m id;

  OC_INT4m job_count; // For development/debugging

  _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadC()
    : mesh(0),energy(0), id(-1), job_count(0) {}

  static void Init(int thread_count,const Oxs_StripedArray<OC_REAL8m>* arrblock) {
    job_basket.Init(thread_count,arrblock);
    result_mutex.Lock();
    etemp = 0.0;
    result_mutex.Unlock();
  }

  void Cmd(int threadnumber, void* data);
};

Oxs_Mutex _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadC::result_mutex;
Oxs_JobControl<OC_REAL8m> _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadC::job_basket;
Nb_Xpfloat _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadC::etemp(0.0);

void _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadC::Cmd(int /* threadnumber */,
                                                  void* /* data */)
{
  Nb_Xpfloat work_etemp_a = 0.0;
  Nb_Xpfloat work_etemp_b = 0.0;
  Nb_Xpfloat work_etemp_c = 0.0;
  Nb_Xpfloat work_etemp_d = 0.0;

  const OC_REAL8m* const senergy  = &((*energy)[OC_INDEX(0)]);
  const Oxs_Mesh* const smesh = mesh;
  while(1) {

    OC_INDEX istart,istop;
    job_basket.GetJob(id,istart,istop);

    if(istart>=istop) break; // No more jobs
    ++job_count;

    const unsigned int LOOP_SIZE = 8;
    const OC_INDEX midstop = istop - (istop-istart)%LOOP_SIZE;

    OC_INDEX i;

    assert(LOOP_SIZE == 8); // The value "8" is used in innermost loop below.
    for(i=istart;i<midstop;i+=8) {
      // Code below assumes LOOP_SIZE == 8
      OC_REAL8m se0 = senergy[i]   * smesh->Volume(i); 
      OC_REAL8m se1 = senergy[i+1] * smesh->Volume(i+1);
      OC_REAL8m se2 = senergy[i+2] * smesh->Volume(i+2);
      OC_REAL8m se3 = senergy[i+3] * smesh->Volume(i+3);

      OC_REAL8m se4 = smesh->Volume(i+4);
      OC_REAL8m se5 = smesh->Volume(i+5);
      OC_REAL8m se6 = smesh->Volume(i+6);
      OC_REAL8m se7 = smesh->Volume(i+7);

      Nb_XpfloatDualAccum(work_etemp_a,se0,work_etemp_b,se1);
      Nb_XpfloatDualAccum(work_etemp_c,se2,work_etemp_d,se3);

      se4 *= senergy[i+4];
      se5 *= senergy[i+5];
      se6 *= senergy[i+6];
      se7 *= senergy[i+7];

      Nb_XpfloatDualAccum(work_etemp_a,se4,work_etemp_b,se5);
      Nb_XpfloatDualAccum(work_etemp_c,se6,work_etemp_d,se7);
    }
    for(;i<istop;++i) {
      work_etemp_a.Accum(senergy[i] * smesh->Volume(i));
    }

  }

  // Sum results into static accumulators
  work_etemp_a += work_etemp_b;
  work_etemp_c += work_etemp_d;
  work_etemp_a += work_etemp_c;
  result_mutex.Lock();
  etemp += work_etemp_a;
  result_mutex.Unlock();
}

class _Oxs_CGEvolve_SetBasePoint_ThreadA : public Oxs_ThreadRunObj {
public:
  static Oxs_Mutex job_control;
  static OC_INDEX offset;

  // Imports
  const Oxs_Mesh* mesh;
  const Oxs_MeshValue<OC_REAL8m>* Ms;
  const Oxs_MeshValue<ThreeVector>* bestpt_mxHxm;

  // basept_mxHxm is used only by Polak-Ribiere method,
  // where it is both import and export
  Oxs_MeshValue<ThreeVector>* basept_mxHxm;

  OC_REAL8m gamma_sum;    // Export --- one per thread
  OC_REAL8m g_sum_sq; // Export --- one per thread

  OC_INDEX vecsize;
  OC_INDEX block_size;

  enum Conjugate_Method { FLETCHER_REEVES, POLAK_RIBIERE } method;

  _Oxs_CGEvolve_SetBasePoint_ThreadA()
    : mesh(0), Ms(0),
       bestpt_mxHxm(0), basept_mxHxm(0),
      gamma_sum(0.0), g_sum_sq(0.0),
      vecsize(0), block_size(0), method(FLETCHER_REEVES) {}

  void Cmd(int threadnumber, void* data);
};

Oxs_Mutex _Oxs_CGEvolve_SetBasePoint_ThreadA::job_control;
OC_INDEX _Oxs_CGEvolve_SetBasePoint_ThreadA::offset(0);

void _Oxs_CGEvolve_SetBasePoint_ThreadA::Cmd(int /* threadnumber */,
                                             void* /* data */)
{
  gamma_sum = g_sum_sq = 0.0;
  while(1) {
    job_control.Lock();
    OC_INDEX istart = offset;
    OC_INDEX istop = ( offset += block_size );
    job_control.Unlock();

    if(istart>=vecsize) break;
    if(istop>vecsize) istop=vecsize;

    const OC_REAL8m* sMs         = &((*Ms)[istart]);
    const ThreeVector* sbest_mxHxm = &((*bestpt_mxHxm)[istart]);

    OC_INDEX jsize = istop - istart;

    if(method != POLAK_RIBIERE) {
      // Fletcher-Reeves method
      OC_INDEX j;
      OC_REAL8m g_sum_sq_0 = 0.0, g_sum_sq_1 = 0.0;
      for(j=0;j<jsize%2;++j) {
        OC_REAL8m cell_scale = sMs[j] * mesh->Volume(istart+j);
	g_sum_sq += sbest_mxHxm[j].MagSq() * cell_scale * cell_scale;
      }
      for(;j<jsize;j+=2) {
        OC_REAL8m cell_scale_0 = sMs[j]   * mesh->Volume(istart+j);
        OC_REAL8m cell_scale_1 = sMs[j+1] * mesh->Volume(istart+j+1);
	g_sum_sq_0 += sbest_mxHxm[j].MagSq()   * cell_scale_0 * cell_scale_0;
	g_sum_sq_1 += sbest_mxHxm[j+1].MagSq() * cell_scale_1 * cell_scale_1;
      }
      g_sum_sq += g_sum_sq_0 + g_sum_sq_1;
      gamma_sum = g_sum_sq;
    } else {
      // Polak-Ribiere method
      ThreeVector* sbase_mxHxm = &((*basept_mxHxm)[istart]);
      for(OC_INDEX j=0;j<jsize;++j) {
        OC_REAL8m cell_scale = sMs[j] * mesh->Volume(istart+j);
        OC_REAL8m cell_scale_sq = cell_scale * cell_scale;
	ThreeVector temp = sbest_mxHxm[j];
	g_sum_sq += temp.MagSq() * cell_scale_sq;
	temp -= sbase_mxHxm[j];  // Polak-Ribiere adjustment
	sbase_mxHxm[j] = sbest_mxHxm[j];
	gamma_sum += cell_scale_sq * (temp * sbest_mxHxm[j]);
      }
    }
  }
}

class _Oxs_CGEvolve_SetBasePoint_ThreadB : public Oxs_ThreadRunObj {
public:
  static Oxs_Mutex job_control;
  static OC_INDEX offset;

  // Imports
  const Oxs_Mesh* mesh;
  const Oxs_MeshValue<OC_REAL8m>* Ms;
  const Oxs_MeshValue<ThreeVector>* spin;
  const Oxs_MeshValue<ThreeVector>* bestpt_mxHxm;

  // basept_direction is both import and export
  Oxs_MeshValue<ThreeVector>* basept_direction;

  OC_REAL8m gamma;  // Import --- same for all threads
  OC_REAL8m Ep;     // Export --- one per thread
  OC_REAL8m maxmagsq;  // Export --- one per thread
  OC_REAL8m normsumsq; // Export --- one per thread

  OC_INDEX vecsize;
  OC_INDEX block_size;

  _Oxs_CGEvolve_SetBasePoint_ThreadB()
    : mesh(0), Ms(0), spin(0),
      bestpt_mxHxm(0), basept_direction(0),
      gamma(0.0), Ep(0.0), normsumsq(0.0),
      vecsize(0), block_size(0) {}

  void Cmd(int threadnumber, void* data);
};

Oxs_Mutex _Oxs_CGEvolve_SetBasePoint_ThreadB::job_control;
OC_INDEX _Oxs_CGEvolve_SetBasePoint_ThreadB::offset(0);

void _Oxs_CGEvolve_SetBasePoint_ThreadB::Cmd(int /* threadnumber */,
                                             void* /* data */)
{
  Ep = maxmagsq = normsumsq = 0.0;

  OC_REAL8m Ep0, Ep1;
  OC_REAL8m maxmagsq0, maxmagsq1, normsumsq0, normsumsq1;
  Ep0 = maxmagsq0 = normsumsq0 = 0.0;
  Ep1 = maxmagsq1 = normsumsq1 = 0.0;

  while(1) {
    job_control.Lock();
    OC_INDEX istart = offset;
    OC_INDEX istop = ( offset += block_size );
    job_control.Unlock();

    if(istart>=vecsize) break;
    if(istop>vecsize) istop=vecsize;

    const OC_REAL8m* sMs         = &((*Ms)[istart]);
    const ThreeVector* sspin  = &((*spin)[istart]);
    const ThreeVector* smxHxm = &((*bestpt_mxHxm)[istart]);
    ThreeVector* sdir   = &((*basept_direction)[istart]);

    OC_INDEX jsize = istop - istart;
    OC_INDEX j;

    for(j=0;j<jsize%2;++j) {
      OC_REAL8m cell_scale = sMs[j] * mesh->Volume(istart+j);
      ThreeVector tmxHxm = smxHxm[j];
      ThreeVector temp = gamma * sdir[j];
      temp.Accum(cell_scale,tmxHxm);
      // Make temp orthogonal to sspin[j].  If the angle between
      // temp and spin[i] is larger than about 45 degrees, then
      // it may be beneficial to divide by spin[i].MagSq(),
      // to offset effects of non-unit spin.  For small angles
      // it appears to be better to leave out this adjustment.
      // Or, better still, use temp -= (temp.m)m formulation
      // for large angles, w/o size correction.
      temp.wxvxw(sspin[j]);
      sdir[j] = temp;

      OC_REAL8m magsq = temp.MagSq();
      if(magsq>maxmagsq) maxmagsq = magsq;
      normsumsq += magsq;

      Ep += (temp * tmxHxm)*cell_scale;
      /// See mjd's NOTES II, 29-May-2002, p156.
    }
    for(;j<jsize;j+=2) {
      OC_REAL8m cell_scale0 = sMs[j] * mesh->Volume(istart+j);
      ThreeVector tmxHxm0 = smxHxm[j];
      ThreeVector temp0 = gamma * sdir[j];
      temp0.Accum(cell_scale0,tmxHxm0);
      temp0.wxvxw(sspin[j]);
      sdir[j] = temp0;

      OC_REAL8m magsq0 = temp0.MagSq();
      if(magsq0>maxmagsq0) maxmagsq0 = magsq0;
      normsumsq0 += magsq0;
      Ep0 += (temp0 * tmxHxm0)*cell_scale0;

      OC_REAL8m cell_scale1 = sMs[j+1] * mesh->Volume(istart+j+1);
      ThreeVector tmxHxm1 = smxHxm[j+1];
      ThreeVector temp1 = gamma * sdir[j+1];
      temp1.Accum(cell_scale1,tmxHxm1);
      temp1.wxvxw(sspin[j+1]);
      sdir[j+1] = temp1;

      OC_REAL8m magsq1 = temp1.MagSq();
      if(magsq1>maxmagsq1) maxmagsq1 = magsq1;
      normsumsq1 += magsq1;
      Ep1 += (temp1 * tmxHxm1)*cell_scale1;
    }
  }

  if(maxmagsq0>maxmagsq) maxmagsq = maxmagsq0;
  if(maxmagsq1>maxmagsq) maxmagsq = maxmagsq1;
  Ep += Ep0 + Ep1;
  normsumsq += normsumsq0 + normsumsq1;
}

class _Oxs_CGEvolve_SetBasePoint_ThreadC : public Oxs_ThreadRunObj { // asdf
public:
  static Oxs_Mutex job_control;
  static OC_INDEX offset;

  // Imports
  const Oxs_Mesh* mesh;
  const Oxs_MeshValue<OC_REAL8m>* Ms;
  const Oxs_MeshValue<ThreeVector>* bestpt_mxHxm;

  // basept_direction is both import and export
  Oxs_MeshValue<ThreeVector>* basept_direction;

  OC_REAL8m maxmagsq; // Export --- one per thread
  OC_REAL8m sumsq;    // Export --- one per thread

  OC_INDEX vecsize;
  OC_INDEX block_size;

  _Oxs_CGEvolve_SetBasePoint_ThreadC()
    : mesh(0), Ms(0),
      bestpt_mxHxm(0), basept_direction(0),
      maxmagsq(0.0), sumsq(0.0),
      vecsize(0), block_size(0) {}

  void Cmd(int threadnumber, void* data);
};

Oxs_Mutex _Oxs_CGEvolve_SetBasePoint_ThreadC::job_control;
OC_INDEX _Oxs_CGEvolve_SetBasePoint_ThreadC::offset(0);

void _Oxs_CGEvolve_SetBasePoint_ThreadC::Cmd(int /* threadnumber */,
                                             void* /* data */)
{
  maxmagsq = sumsq = 0.0;
  OC_REAL8m sumsq_0 = 0.0, sumsq_1 = 0.0;
  OC_REAL8m maxmagsq_0 = 0.0, maxmagsq_1 = 0.0;
  while(1) {
    job_control.Lock();
    OC_INDEX istart = offset;
    OC_INDEX istop = ( offset += block_size );
    job_control.Unlock();

    if(istart>=vecsize) break;
    if(istop>vecsize) istop=vecsize;

    const OC_REAL8m* sMs         = &((*Ms)[istart]);
    const ThreeVector* sbest_mxHxm = &((*bestpt_mxHxm)[istart]);
    ThreeVector* sbase_direction = &((*basept_direction)[istart]);

    OC_INDEX jsize = istop - istart;

    OC_INDEX j;
    for(j=0;j<jsize%2;++j) {
      OC_REAL8m cell_scale = sMs[j] * mesh->Volume(istart+j);
      sbase_direction[j] = cell_scale * sbest_mxHxm[j];
      OC_REAL8m magsq = sbase_direction[j].MagSq();
      if(magsq>maxmagsq) maxmagsq = magsq;
      sumsq += magsq;
      /// See mjd's NOTES II, 29-May-2002, p156.
    }
    for(;j<jsize;j+=2) {
      OC_REAL8m cell_scale_0 = sMs[j]   * mesh->Volume(istart+j);
      OC_REAL8m cell_scale_1 = sMs[j+1] * mesh->Volume(istart+j+1);

      ThreeVector temp0 = sbase_direction[j]
        = cell_scale_0 * sbest_mxHxm[j];
      OC_REAL8m magsq_0 = temp0.MagSq();
      if(magsq_0>maxmagsq_0) maxmagsq_0 = magsq_0;

      ThreeVector temp1 = sbase_direction[j+1]
        = cell_scale_1 * sbest_mxHxm[j+1];
      OC_REAL8m magsq_1 = temp1.MagSq();
      if(magsq_1>maxmagsq_1) maxmagsq_1 = magsq_1;

      sumsq_0 += magsq_0;
      sumsq_1 += magsq_1;
    }
  }

  sumsq += sumsq_0 + sumsq_1;
  if(maxmagsq_0>maxmagsq) maxmagsq = maxmagsq_0;
  if(maxmagsq_1>maxmagsq) maxmagsq = maxmagsq_1;
}

#endif // OOMMF_THREADS

void Oxs_CGEvolve::GetEnergyAndmxHxm
(const Oxs_SimState* state,          // Import
 Oxs_MeshValue<OC_REAL8m>& export_energy,      // Export
 Oxs_MeshValue<ThreeVector>& export_mxHxm,  // Export
 Oxs_MeshValue<ThreeVector>* export_Hptr)   // Export
{ // Fills export_energy and export_mxHxm, which must be different than
  // scratch_energy and scratch_mxHxm.  Also fills export_Hptr with
  // total field, unless export_Hptr == NULL, in which case the total
  // field is not saved.
  //   This routine also updates energy_calc_count, and fills "Total
  // energy", "Bracket_count", "Line min count", "Energy calc count",
  // and "Cycle count" derived data into state.  These class values
  // should be updated as desired before calling this routine.
  //   SIDE EFFECTS: scratch_energy & scratch_field are altered
  // Note: (mxH)xm = mx(Hxm) = -mx(mxH)

#if REPORT_TIME_CGDEVEL
  getenergyandmxHxmtime.Start();
  OC_BOOL sot_running = steponlytime.IsRunning();
  OC_BOOL bpt_running = basepttime.IsRunning(); 
  OC_BOOL fbt_running = findbrackettime.IsRunning();
  OC_BOOL flt_running = findlinemintime.IsRunning();
  OC_BOOL fill_running = fillbrackettime.IsRunning();
#endif // REPORT_TIME_CGDEVEL

  if(&export_energy == &scratch_energy) {
      throw Oxs_ExtError(this,
        "Programming error in Oxs_CGEvolve::GetEnergyAndmxHxm():"
	" export_energy is same as scratch_energy.");
  }
  if(&export_mxHxm == &scratch_field) {
      throw Oxs_ExtError(this,
        "Programming error in Oxs_CGEvolve::GetEnergyAndmxHxm():"
	" export_mxHxm is same as scratch_field.");
  }
  if(export_Hptr == &scratch_field) {
      throw Oxs_ExtError(this,
        "Programming error in Oxs_CGEvolve::GetEnergyAndmxHxm():"
	" export_Hptr is same as scratch_field.");
  }

  ++energy_calc_count;    // Update call count

  // For convenience
  const Oxs_Mesh* mesh = state->mesh;
  const OC_INDEX vecsize = mesh->Size();

  // Set up energy computation output data structure
  Oxs_ComputeEnergyData oced(*state);
  oced.scratch_energy = &scratch_energy;
  oced.scratch_H      = &scratch_field;
  oced.energy_accum   = &export_energy;
  oced.H_accum        = export_Hptr;
  oced.mxH_accum      = 0;  // Fill export_mxHxm instead
  oced.energy         = NULL;
  oced.H              = NULL;
  oced.mxH            = NULL;

  UpdateFixedSpinList(mesh);
  Oxs_ComputeEnergyExtraData oceed(GetFixedSpinList(),
                                   &export_mxHxm);

  // Compute total energy and torque
#if REPORT_TIME
    if(sot_running) {
      steponlytime.Stop();
#if REPORT_TIME_CGDEVEL
      if(bpt_running)  basepttime.Stop();
      if(fbt_running)  findbrackettime.Stop();
      if(flt_running)  findlinemintime.Stop();
      if(fill_running) fillbrackettime.Stop();
      energyobjtime.Start();
#endif // REPORT_TIME_CGDEVEL
    }
#if REPORT_TIME_CGDEVEL
    getenergyandmxHxmtime.Stop();
#endif // REPORT_TIME_CGDEVEL
#endif // REPORT_TIME
    Oxs_ComputeEnergies(*state,oced,director->GetEnergyObjects(),oceed);
#if REPORT_TIME
#if REPORT_TIME_CGDEVEL
    getenergyandmxHxmtime.Start();
#endif // REPORT_TIME_CGDEVEL
    if(sot_running) {
#if REPORT_TIME_CGDEVEL
      energyobjtime.Stop();
      if(bpt_running) basepttime.Start();
      if(fbt_running) findbrackettime.Start();
      if(flt_running) findlinemintime.Start();
      if(fill_running) fillbrackettime.Start();
#endif // REPORT_TIME_CGDEVEL
      steponlytime.Start();
    }
#endif // REPORT_TIME
    if(oced.pE_pt != 0.0) {
      String msg = 
        String("Oxs_CGEvolve::GetEnergyAndmxHxm:"
               " At least one energy object is time varying; this"
               " property is not supported by minimization evolvers.");
      throw Oxs_Ext::Error(this,msg.c_str());
    }

  // Fill supplemental derived data.
  state->AddDerivedData("Max mxHxm",oceed.max_mxH);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if 0
  { // asdfasdf
    // Flush all cpu caches
#   define FLUSH_COUNT 8
    static Oxs_MeshValue<OC_REAL8m> flusher[FLUSH_COUNT];
    if(flusher[0].Size()==0) {
      for(int i=0;i<FLUSH_COUNT;++i) {
        flusher[i].AdjustSize(mesh);
        for(OC_INDEX k=0;k<mesh->Size();++k) {
          flusher[i][k] = 1.7e-15;
        }
      }
    }

    Oxs_ThreadTree threadtree;
    const OC_INT4m MaxThreadCount = Oc_GetMaxThreadCount();
    vector<_Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadC> threadC;
    threadC.resize(MaxThreadCount);

    _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadC::Init(MaxThreadCount,
                                         export_energy.GetArrayBlock());

    for(int i=0;i<FLUSH_COUNT;++i) {
      for(OC_INT4m ithread=0;ithread<MaxThreadCount;++ithread) {
        threadC[ithread].id = ithread;
        threadC[ithread].mesh = mesh;
        // threadC[ithread].energy = &export_energy;
        threadC[ithread].energy = &flusher[i];
        if(ithread>0) threadtree.Launch(threadC[ithread],0);
      }
      threadtree.LaunchRoot(threadC[0],0);
    }
  }
#endif
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if REPORT_TIME_CGDEVEL
timer[8].Start(); /**/ // 2 secs
#endif // REPORT_TIME_CGDEVEL

#if !OOMMF_THREADS
  Nb_Xpfloat total_energy = 0.0;
  for(OC_INDEX i=0;i<vecsize;++i) {
    total_energy.Accum(export_energy[i] * mesh->Volume(i));
  }
  state->AddDerivedData("Total energy",total_energy.GetValue());
#else // OOMMF_THREADS
  {
    Oxs_ThreadTree threadtree;
    const OC_INT4m MaxThreadCount = Oc_GetMaxThreadCount();
    vector<_Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadC> threadC;
    threadC.resize(MaxThreadCount);

    _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadC::Init(MaxThreadCount,
                                         export_energy.GetArrayBlock());

    for(OC_INT4m ithread=0;ithread<MaxThreadCount;++ithread) {
      threadC[ithread].id = ithread;
      threadC[ithread].mesh = mesh;
      threadC[ithread].energy = &export_energy;
      if(ithread>0) threadtree.Launch(threadC[ithread],0);
    }
    threadtree.LaunchRoot(threadC[0],0);
    state->AddDerivedData("Total energy",
             _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadC::etemp.GetValue());
  }
#endif // OOMMF_THREADS

#if REPORT_TIME_CGDEVEL
timer[8].Stop(); /**/
++timer_counts[8].pass_count;
 timer_counts[8].bytes += vecsize*(1*sizeof(OC_REAL8m)); // Volume(i) call
 /// can be done w/o main memory access
#endif // REPORT_TIME_CGDEVEL

  state->AddDerivedData("Bracket count",
			OC_REAL8m(bracket_count));
  state->AddDerivedData("Line min count",
			OC_REAL8m(line_minimum_count));
  state->AddDerivedData("Energy calc count",
			OC_REAL8m(energy_calc_count));
  state->AddDerivedData("Cycle count",
			OC_REAL8m(cycle_count));
  state->AddDerivedData("Cycle sub count",
			OC_REAL8m(cycle_sub_count));
#if REPORT_TIME_CGDEVEL
  getenergyandmxHxmtime.Stop();
#endif // REPORT_TIME_CGDEVEL
}

void Oxs_CGEvolve::GetRelativeEnergyAndDerivative(
 const Oxs_SimState* state, // Import
 OC_REAL8m offset,             // Import
 OC_REAL8m &relenergy,         // Export
 OC_REAL8m &derivative,        // Export
 OC_REAL8m &grad_norm)         // Export
{ // Calculates total energy relative to that of best_state,
  // and derivative in base_direction (-mu0.H*base_direction).
  // The base_direction and best_energy arrays _must_ be set
  // before calling this function.

  // Uses mxHxm instead of H in calculation of derivative, because
  // energy is calculated off of normalized m, so component of H in m
  // direction doesn't actually have any effect.  Plus, experiments
  // appear to show generally faster convergence with mxHxm than with H.

  temp_id = 0;
  GetEnergyAndmxHxm(state,temp_energy,temp_mxHxm,NULL);
  temp_id = state->Id();

  // Set up some convenience variables
  const Oxs_Mesh* mesh = state->mesh; // For convenience
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state->Ms);

  const OC_INDEX vecsize = mesh->Size();

#if REPORT_TIME_CGDEVEL
timer[1].Start(); /**/ // 2secs
#endif // REPORT_TIME_CGDEVEL

#if !OOMMF_THREADS
  OC_INDEX i;
  Nb_Xpfloat etemp = 0.0;
  for(i=0;i<vecsize;++i) {
    etemp.Accum((temp_energy[i] - bestpt.energy[i]) * mesh->Volume(i));
  }

  Nb_Xpfloat dtemp = 0.0;
  Nb_Xpfloat stemp = 0.0;
  OC_REAL8m offset_sq = offset*offset;
  for(i=0;i<vecsize;++i) {
    const ThreeVector& vtemp = basept.direction[i];
    OC_REAL8m scale_adj = Ms[i]*mesh->Volume(i)
      /sqrt(1+offset_sq * vtemp.MagSq());
    dtemp.Accum((temp_mxHxm[i]*vtemp)*scale_adj);
    stemp.Accum(temp_mxHxm[i].MagSq()*scale_adj*scale_adj);
  }

  relenergy = etemp.GetValue();
  derivative = -MU0 * dtemp.GetValue();
  grad_norm = sqrt(stemp.GetValue());
  /// See mjd's NOTES II, 29-May-2002, p156, which includes
  /// the derivation of the scale_adj term above.
#else // OOMMF_THREADS
      {
        Oxs_ThreadTree threadtree;
        const OC_INT4m MaxThreadCount = Oc_GetMaxThreadCount();
        vector<_Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB> threadB;
        threadB.resize(MaxThreadCount);

        _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::job_control.Lock();
        _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::offset = 0;
        _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::etemp = 0.0;
        _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::dtemp = 0.0;
        _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::stemp = 0.0;
        _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::job_control.Unlock();

        OC_INDEX ibs = (vecsize + MaxThreadCount - 1)/MaxThreadCount;
        ibs = (ibs + 7)/8;
        if(ibs%8) { ibs += 8 - (ibs%8); }
        if(ibs> static_cast<OC_INDEX>((vecsize*3+3)/4)) ibs = vecsize;

        for(OC_INT4m ithread=0;ithread<MaxThreadCount;++ithread) {
          threadB[ithread].mesh = mesh;
          threadB[ithread].tmpenergy = &temp_energy;
          threadB[ithread].bestpt_energy = &bestpt.energy;
          threadB[ithread].Ms = &Ms;
          threadB[ithread].direction = &basept.direction;
          threadB[ithread].mxHxm = &temp_mxHxm;
          threadB[ithread].offset_sq = offset*offset;
          threadB[ithread].vecsize = vecsize;
          threadB[ithread].block_size = ibs;
          if(ithread>0) threadtree.Launch(threadB[ithread],0);
        }
        threadtree.LaunchRoot(threadB[0],0);

        relenergy
          = _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::etemp.GetValue();
        derivative
          = -MU0 * _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::dtemp.GetValue();
        grad_norm
          = sqrt(_Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadB::stemp.GetValue());
        /// See mjd's NOTES II, 29-May-2002, p156, which includes
        /// the derivation of the scale_adj term above.
      }
#endif // OOMMF_THREADS

#if REPORT_TIME_CGDEVEL
timer[1].Stop(); /**/
++timer_counts[1].pass_count;
 timer_counts[1].bytes += vecsize*(3*sizeof(OC_REAL8m)+2*sizeof(ThreeVector));
#endif // REPORT_TIME_CGDEVEL

  state->AddDerivedData("Relative energy",relenergy);
  state->AddDerivedData("Energy best state id",
			static_cast<OC_REAL8m>(bestpt.key.ObjectId()));
}

void Oxs_CGEvolve::FillBracket(OC_REAL8m offset,
			       const Oxs_SimState* oldstateptr,
			       Oxs_Key<Oxs_SimState>& statekey,
			       Bracket& endpt)
{
  // bestpt and basept must be set prior to calling this routine
#if REPORT_TIME_CGDEVEL
  fillbrackettime.Start();
#endif
  try {
    Oxs_SimState& workstate
      = statekey.GetWriteReference(); // Write lock
    Oxs_MeshValue<ThreeVector>& spin = workstate.spin;
    const Oxs_MeshValue<ThreeVector>& best_spin
      = bestpt.key.GetPtr()->spin;
    OC_INDEX vecsize = workstate.mesh->Size();
    OC_REAL8m t1sq = bestpt.offset*bestpt.offset;

#if REPORT_TIME_CGDEVEL
timer[0].Start(); /**/ // 7 secs
#endif // REPORT_TIME_CGDEVEL

#if !OOMMF_THREADS
    for(OC_INDEX i=0;i<vecsize;i++) {
      const ThreeVector& dvec = basept.direction[i];
      OC_REAL8m dsq = dvec.MagSq();
      ThreeVector temp = best_spin[i];
      temp *= sqrt(1+t1sq*dsq);
      temp.Accum((offset-bestpt.offset),dvec);
      temp.MakeUnit();
      Nb_NOP(&temp);  // Safety; icpc optimizer showed problems
      spin[i] = temp; // with this stanza in some circumstances.
    }
#else // OOMMF_THREADS
      {
        Oxs_ThreadTree threadtree;
        const OC_INT4m MaxThreadCount = Oc_GetMaxThreadCount();
        vector<_Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadA> threadA;
        threadA.resize(MaxThreadCount);

        _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadA::job_control.Lock();
        _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadA::offset = 0;
        _Oxs_CGEvolve_GetEnergyAndmxHxm_ThreadA::job_control.Unlock();

        OC_INDEX ibs = (vecsize + MaxThreadCount - 1)/MaxThreadCount;
        ibs = (ibs + 7)/8;
        if(ibs%8) { ibs += 8 - (ibs%8); }
        if(ibs> static_cast<OC_INDEX>((vecsize*3+3)/4)) ibs = vecsize;

        for(OC_INT4m ithread=0;ithread<MaxThreadCount;++ithread) {
          threadA[ithread].direction = &basept.direction;
          threadA[ithread].best_spin = &best_spin;
          threadA[ithread].spin = &spin;
          threadA[ithread].tsq = t1sq;
          threadA[ithread].vecsize = vecsize;
          threadA[ithread].dvec_scale = offset-bestpt.offset;
          threadA[ithread].block_size = ibs;
          if(ithread>0) threadtree.Launch(threadA[ithread],0);
        }
        threadtree.LaunchRoot(threadA[0],0);
      }
#endif // OOMMF_THREADS

#if REPORT_TIME_CGDEVEL
timer[0].Stop(); /**/
++timer_counts[0].pass_count;
timer_counts[0].bytes += vecsize*(3*sizeof(ThreeVector));
#endif // REPORT_TIME_CGDEVEL

    workstate.iteration_count = oldstateptr->iteration_count + 1;
    workstate.stage_iteration_count
      = oldstateptr->stage_iteration_count + 1;

    // Release write lock and update energy data
    const Oxs_SimState& newstate = statekey.GetDepReference();
#if REPORT_TIME_CGDEVEL
    fillbrackettime.Stop();
#endif
    GetRelativeEnergyAndDerivative(&newstate,offset,
                                   endpt.E,endpt.Ep,endpt.grad_norm);
#if REPORT_TIME_CGDEVEL
    fillbrackettime.Start();
#endif
    endpt.offset=offset;
  }
  catch (...) {
    statekey.GetReadReference(); // Release write lock
    throw;
  }
#if REPORT_TIME_CGDEVEL
  fillbrackettime.Stop();
#endif
}

void
Oxs_CGEvolve::UpdateBrackets
(Oxs_ConstKey<Oxs_SimState> tstate_key,
 const Bracket& tbracket,
 OC_BOOL force_bestpt)
{ // Incorporates tbracket into existing brackets, and changes
  // bestpt if appropriate.  Returns 1 if bestpt was changed,
  // in which case bracket E values and temp space data will
  // be modified.
  //  If force_bestpt is true, and temp_id agrees with
  // tstate_key, then tstate is made new "bestpt"

  // Manage brackets
  if(tbracket.offset>bracket.right.offset) {
    // Shift brackets to the right
    if(bracket.right.E<bracket.left.E
       && bracket.right.Ep<0) { // This should always be
      /// true, since otherwise (left,right] already bracketed the
      /// minimum, but check anyway for robustness.
      bracket.left=bracket.right;
    }
    bracket.right = tbracket;
  } else {
    // Collapse brackets
    OC_REAL8m energy_slack 
      = bracket.energy_precision * fabs(basept.total_energy);
    OC_REAL8m tlspan = tbracket.offset - bracket.left.offset;
    if(tbracket.Ep<0
       && (tbracket.E<=bracket.right.E || bracket.right.Ep>=0)
       && tbracket.E<bracket.left.E+energy_slack
       && fabs(tlspan*bracket.left.Ep) < energy_slack) {
      // Keep right interval.  Note: It may be that leftpt + tbracket
      // also bracket a minimum, but when these two E values are within
      // energy_slack, it is a good bet that the lefthand interval has
      // fallen below the noise floor, in which case we are better off
      // keeping the farther minimum.
      if(bracket.left.offset == bestpt.offset) {
	force_bestpt = 1; // Force shift of bestpt from
	/// old leftpt to new leftpt.  This increases
	/// bestpt E by at most energy_slack.
      }
      bracket.left=tbracket;
    } else if(tbracket.E>bracket.left.E || tbracket.Ep>=0) {
      // Keep left interval
      bracket.right=tbracket;
    } else {
      // Keep right interval
      bracket.left=tbracket;
    }
  }

  // Update bestpt
  if((force_bestpt || tbracket.E<=0)
     && temp_id==tstate_key.ObjectId()) {
    // New best point found
    temp_id=bestpt.key.ObjectId();
    bestpt.key = tstate_key;
    bestpt.key.GetReadReference(); // Set read lock
    bestpt.offset = tbracket.offset;
    bestpt.Ep = tbracket.Ep;
    bestpt.grad_norm = tbracket.grad_norm;
    bestpt.energy.Swap(temp_energy);
    bestpt.mxHxm.Swap(temp_mxHxm);
    bracket.left.E -= tbracket.E;
    bracket.right.E -= tbracket.E;
  }

}
				

void Oxs_CGEvolve::FindBracketStep(const Oxs_SimState* oldstate,
				   Oxs_Key<Oxs_SimState>& statekey)
{ // Takes one step; assumes right.offset >= left.offset and
  // minimum to be bracketed is to the right of right.offset
  assert(bracket.right.offset>=bracket.left.offset && bracket.left.Ep<0.0);
  ++bracket_count;
  OC_REAL8m offset = bracket.right.offset;
#ifndef NDEBUG
  const OC_REAL8m starting_offset = offset;
#endif
  if(offset < bracket.start_step) {
    offset = bracket.start_step;
  } else {
    if(offset<bracket.left.offset+bracket.scaled_minstep) {
      offset = bracket.left.offset*2 + bracket.scaled_minstep;
    } else {
      // Try to use derivative information at existing endpoints
      // to estimate location of minimum.
      OC_REAL8m temp1 = (bracket.left.offset-bracket.right.offset)
	*bracket.left.Ep;
      OC_REAL8m temp2 = bracket.right.Ep-bracket.left.Ep;
      if(temp2<=0) {
	offset *=4; // Aiee! We're rolling off a cliff!
      } else {
	OC_REAL8m temp3 = 8*bracket.right.offset-bracket.left.offset;
	if(2*temp1 < temp3*temp2) { // Overflow check
	  offset = bracket.left.offset + 2*temp1/temp2;
	} else {
	  offset = 8*bracket.right.offset;
	}
      }
    }
  }
  if(offset < 2*bracket.right.offset) {
    offset = 2*bracket.right.offset;
  }
  if(offset > bracket.scaled_maxstep) {
    offset = bracket.scaled_maxstep;
  }

  // Zero span check
  if(offset <= bracket.right.offset) {
    // Proposed bracket span increase is floating point 0.
    if(bracket.right.offset>0) {
      offset = bracket.right.offset*(1+16*OC_REAL8_EPSILON);
    } else {
      if(bracket.scaled_maxstep>0) {
	offset = bracket.scaled_maxstep;
      } else {
	// Punt
	offset = OC_REAL8_EPSILON;
      }
    }
  }
  if(fabs(offset-bestpt.offset)*basept.direction_max_mag
     <16*OC_REAL8_EPSILON) {
    // If LHS above is smaller that OC_REAL8_EPSILON, and we assume
    // that the bestpt spins are unit vectors, then there is a
    // good chance that the proposed offset won't actually
    // change any of the discretized spins.  For bracketing
    // purposes, it shouldn't hurt to bump this up a bit.
    if(basept.direction_max_mag >= 0.5 ||
       DBL_MAX*basept.direction_max_mag > 16*OC_REAL8_EPSILON) {
      offset += 16*OC_REAL8_EPSILON/basept.direction_max_mag;
    }
  }


  Bracket temppt;
  FillBracket(offset,oldstate,statekey,temppt);

  UpdateBrackets(statekey,temppt,0);

  // Classify situation.  We always require leftpt have Ep<0,
  // so there are two possibilities:
  //  1) temppt.E>bracket.left.E+energy_slack or temppt.Ep>=0
  //          ==> Minimum is bracketed
  //  2) temppt.E<=bracket.left.E+energy_slack and temppt.Ep<0
  //          ==> Minimum not bracketed
  OC_REAL8m energy_slack 
    = bracket.energy_precision * fabs(basept.total_energy);
  if(temppt.E<=bracket.left.E+energy_slack && temppt.Ep<0) {
    // Minimum not bracketed.
    bracket.min_bracketed=0; // Safety
    bracket.stop_span = 0.0;
    if(bracket.right.offset >= bracket.scaled_maxstep) {
      // Unable to bracket minimum inside allowed range.  Collapse
      // bracket to bestpt, and mark line minimization complete.  NOTE:
      // This may cause solver to get stuck with bestpt at basept, if
      // the surface is very flat and noisy.  May want to include some
      // logic either here or in UpdateBracket() to detect this
      // situation and move bestpt away from basept.
      bracket.min_bracketed=1;
      bracket.min_found = 1;
      bracket.right.offset = bestpt.offset;
      bracket.right.E = 0.0;
      bracket.right.Ep = bestpt.Ep;
      bracket.right.grad_norm = bestpt.grad_norm;
      bracket.left = bracket.right;
    }
  } else {
    // Minimum is bracketed.
    bracket.min_bracketed=1;
    bracket.stop_span = bracket.relative_minspan*bracket.right.offset;
    if(bracket.stop_span*basept.direction_max_mag < 4*OC_REAL8_EPSILON) {
      bracket.stop_span = 4*OC_REAL8_EPSILON/basept.direction_max_mag;
      /// This is a reasonable guess at the smallest variation in
      /// offset that leads to a detectable (i.e., discretizational)
      /// change in spins.
    }
  }

#ifndef NDEBUG
  assert(bracket.min_bracketed ||
         (bracket.right.offset>0 && bracket.right.offset>=2*starting_offset));
#endif
}


void Oxs_CGEvolve::FindLineMinimumStep(const Oxs_SimState* oldstate,
				       Oxs_Key<Oxs_SimState>& statekey)
{ // Compresses bracket span.
  // Assumes coming in and forces on exit requirement that leftpt.Ep<0,
  // and either rightpt.E>leftpt.E or rightpt.Ep>=0.
  assert(bracket.left.Ep<=0.0
	 && (bracket.right.E>bracket.left.E || bracket.right.Ep>=0));
  OC_REAL8m span = bracket.right.offset - bracket.left.offset; 
  OC_REAL8m energy_slack 
    = bracket.energy_precision * fabs(basept.total_energy);

  OC_REAL8m nudge = DBL_MAX/2;
  if(basept.direction_max_mag>=1.
     || OC_REAL8_EPSILON<nudge*basept.direction_max_mag) {
    nudge = OC_REAL8_EPSILON/basept.direction_max_mag;
    /// nudge is a reasonable guess at the smallest variation in
    /// offset that leads to a detectable (i.e., discretizational)
    /// change in spins.
    if(nudge>=0.5*span*(1-OC_REAL8_EPSILON) && nudge<=span*(1+OC_REAL8_EPSILON)) {
      // fudge nudge.  Strictly speaking, if nudge really represented
      // the numerically smallest possible non-zero spin tweak, and if
      // also span/2<nudge<span, then no tweaking could help, because
      // all points between the brackets numerically match one end or
      // the other.  In truth, however, nudge is probably not tight, so
      // if nudge is inside (span/2,span), round it down to span/2 and
      // give FindLineMinimumStep() one more shot.
      nudge = 0.5*span*(1-OC_REAL8_EPSILON);
    }
  } else {
    // Degenerate case.  basept.direction_max_mag must be exactly
    // or very nearly zero.  Punt.
    nudge = span; // This will kick out in one of the following
    /// two stanzas.
  }

  // The first test below is an orthogonality check.
  //   The second test is a rough check to see if applying the
  // conjugation procedure to the gradient at this point will yield a
  // downhill direction.  This condition is trivially obtained if
  // bestpt.Ep is 0, so it should be always possible in principle,
  // ignoring discretization effects.  This check here is not exactly
  // the same as in SetBasePoint(), because a) it uses the
  // Fletcher-Reeves conjugation procedure, regardless of the selected
  // Conjugate_Method, and b) bestpt.Ep, used here, is computed in
  // GetRelativeEnergyAndDerivative() with a "scale_adj" that differs
  // from that used in SetBasePoint().
  //   The last test group are sanity checks, and the obsoleted span
  // size control.
  if(fabs(bestpt.Ep)
     < MU0*bestpt.grad_norm*basept.direction_norm*bracket.angle_precision
     && bestpt.Ep<=MU0*basept.g_sum_sq
     && (bestpt.Ep==0 || span<=bracket.stop_span || nudge>=span)) {

    // Minimum found
    bracket.min_found=1;
    bestpt.is_line_minimum=1; // Good gradient info
    bracket.last_min_reduction_ratio=0.0;
    bracket.next_to_last_min_reduction_ratio=0.0;
    assert(bracket.left.Ep<=0.0
           && (bracket.right.E>bracket.left.E || bracket.right.Ep>=0));
    return;
  }

  if(bracket.left.Ep>=0. || nudge>=span*(1-OC_REAL8_EPSILON)
     || bracket.right.Ep==0.) {
    /// left.Ep==0 means minimum found. >0 is an error, but included
    ///    anyway for robustness.
    /// If nudge>=span, than further refinement doesn't appreciably
    ///    change any spins, so we might as well stop.
    /// right.Ep can have either sign, but regardless exact 0
    ///    indicate minimum.
    bracket.min_found=1;
    bestpt.is_line_minimum=0; // Bad gradient info
    bracket.last_min_reduction_ratio=0.0;
    bracket.next_to_last_min_reduction_ratio=0.0;
    assert(bracket.left.Ep<=0.0
	   && (bracket.right.E>bracket.left.E || bracket.right.Ep>=0));
    return;
  }

  // Pick test point.
  OC_REAL8m test_offset;
  OC_REAL8m lambda = 0.5; // Default
  Bracket testpt;
  OC_REAL8m lEp = bracket.left.Ep * span;
  OC_REAL8m rEp = bracket.right.Ep * span;
  OC_REAL8m Ediff = bracket.right.E - bracket.left.E;

  // Test for total numerics loss
  if(span<=256*bracket.stop_span
     && fabs(Ediff)<=energy_slack
     && fabs(rEp-lEp)<fabs(lEp)/16.
     && fabs(lEp) < energy_slack) {
    // Yeah, this could happen by chance, but it's
    // not very likely (I hope...) -mjd, Feb-2003
    bracket.min_found=1;
    bracket.last_min_reduction_ratio=0.0;
    bracket.next_to_last_min_reduction_ratio=0.0;
    return;
  }

  if(fabs(Ediff)>energy_slack && rEp>=0.) {
    // Cubic fit.  See mjd's NOTES II, 1-Sep-2001, p134-136.
    OC_REAL8m a = -2*Ediff +   lEp + rEp;
    OC_REAL8m b =  3*Ediff - 2*lEp - rEp;
    OC_REAL8m c = lEp;
    if(a==0.0) {
      // Quadratic
      if(b!=0.0) { // Safety check. b should be >= -c/2 > 0.
	lambda = -c/(2*b);
      }
    } else {
      OC_REAL8m disc = b*b - 3*a*c;
      if(disc<=0.0) disc=0.0;          // Safety check.  See NOTES II,
      else          disc = sqrt(disc); // 1-Sep-2001, p135.
      if(b>=0.) {
	if(fabs(c)>=b+disc) lambda = Nb_Signum(-c);
	else                lambda = -c/(b + disc);
      } else {
	if(fabs(3*a)<=(-b + disc)) lambda = Nb_Signum(a);
	else                       lambda = (-b + disc)/(3*a);
      }
    }
  } else if(rEp>=0.) {
    // E data looks iffy, so just use derivative information,
    // i.e., quadratic fit.
    OC_REAL8m Ep_diff = bracket.right.Ep - bracket.left.Ep;
    /// Since left.Ep<0,/ Ep_diff = |right.Ep| + |left.Ep|.
    lambda = -bracket.left.Ep/Ep_diff;
    // Provided Ep_diff is a valid floating point number,
    // then lambda should not overflow, and in fact must
    // lay inside [0,1].
  } else {
    // Overly large bracket, or bad rightpt.Ep.
    // Guess at minpt using leftpt.E/Ep and rightpt.E
    const OC_REAL8m reduce_limit = 1./32.; // The
    /// situation rightpt.Ep<0 is a suspicious, so limit
    /// reduction.
    OC_REAL8m numerator = -1.0 * lEp; // This is > 0.
    OC_REAL8m denominator = 2*(Ediff - lEp);
    // denominator must also be >0, because rightpt.E>=leftpt.E
    // if a minimum is bracketed with rightpt.Ep<0.
    denominator = fabs(denominator); // But play it safe, because
    /// checks below depend on denominator>0.
    if(numerator<reduce_limit*denominator) {
      lambda = reduce_limit;
    } else if(numerator>(1-reduce_limit)*denominator) {
      lambda = 1-reduce_limit;
    } else {
      lambda = numerator/denominator;
    }
  }

  // Curb reduction a bit to improve odds that minimum
  // lands in smaller interval
  const OC_REAL8m SAFETY = 1.0/(1024.0*1024.0);
  if(lambda<0.25)      lambda *= (1.0+SAFETY);
  else if(lambda>0.75) lambda *= (1.0-SAFETY);

  // Restrict reduction.
  const OC_REAL8m max_reduce_base=0.5; // Don't set this any larger than 0.707
  OC_REAL8m max_reduce=max_reduce_base;
  OC_REAL8m last_reduce = bracket.last_min_reduction_ratio;
  OC_REAL8m next_to_last_reduce = bracket.next_to_last_min_reduction_ratio;
  if(last_reduce<max_reduce) max_reduce=last_reduce;
  if(next_to_last_reduce<max_reduce) max_reduce=next_to_last_reduce;
  max_reduce *= max_reduce;
  if(span*max_reduce < OC_REAL8_EPSILON * bracket.right.offset) {
    OC_REAL8m temp = OC_REAL8_EPSILON * bracket.right.offset;
    if(temp<0.5*span) max_reduce = temp/span;
    else              max_reduce = 0.5;
  }
  if(/* span>256*nudge && */ span*max_reduce<nudge) max_reduce = nudge/span;
  /// There are some difficulties with the above nudge barrier control.
  /// In particular, even if the difference in offset values between
  /// test_offset and left/right.offset is less than OC_REAL8_EPSILON,
  /// the difference in the spin configuration at the test point may
  /// be different than the spin configuration at the closer endpt
  /// because of roundoff errors.  Some tests have appeared to indicate
  /// that this sloppiness aids in the small mxHxm situation, so we
  /// restrict the nudge barrier operation to just those intervals that
  /// are wide relative to nudge.
  assert(max_reduce<=0.5);
  if(max_reduce>0.5) max_reduce=0.5; // Safety; This shouldn't trigger.
  if(lambda>0.5) {
    // Don't reverse the order of these two branches (lambda>0.5
    // vs. lambda<=0.5), because otherwise we run afoul of a bug in g++
    // 3.2.3 with -O3 optimization.  The bug manifests by effectively
    // removing this block, allowing the solver to get stuck in a
    // sequence of barely reducing intervals.
    if(lambda > 1.0-max_reduce) lambda = 1.0-max_reduce;
  } else {
    if(lambda < max_reduce) lambda = max_reduce;
  }

  test_offset = bracket.left.offset + lambda * span;
  if(test_offset<=bracket.left.offset
     || test_offset>=bracket.right.offset) {
    // Roundoff check
    lambda = 0.5;
    test_offset = 0.5 * (bracket.left.offset + bracket.right.offset);
    if(test_offset<=bracket.left.offset 
       || test_offset>=bracket.right.offset) {
      // Interval width effectively machine 0; assume minimum reached
      bracket.min_found=1;
      bracket.last_min_reduction_ratio=0.0;
      bracket.next_to_last_min_reduction_ratio=0.0;
      return;
    }
  }

  ++line_minimum_count;
  FillBracket(test_offset,oldstate,statekey,testpt);

  // Determine which interval, [left,test] or [test,right] to keep.
  UpdateBrackets(statekey,testpt,0);
  OC_REAL8m newspan = bracket.right.offset - bracket.left.offset;
  assert(bracket.next_to_last_min_reduction_ratio
         *bracket.last_min_reduction_ratio
         *newspan/span<1-max_reduce_base*max_reduce_base);
  bracket.next_to_last_min_reduction_ratio
    = bracket.last_min_reduction_ratio;
  bracket.last_min_reduction_ratio =  newspan/span;
#ifdef INSTRUMENT
/**/ fprintf(stderr," --Span ratio:%#10.7f\n",
             static_cast<double>(newspan/span));
#endif // INSTRUMENT

  // Is minimum found?  The second test here is a rough check to
  // see if the conjugation procedure applied to mxHxm will yield
  // a downhill direction.  See note in the up-front check for
  // additional details.
  if(fabs(bestpt.Ep)
     < MU0*bestpt.grad_norm*basept.direction_norm*bracket.angle_precision
     && bestpt.Ep>MU0*basept.g_sum_sq
     && (bestpt.Ep==0 || newspan<=bracket.stop_span || nudge>=newspan)) {
    bracket.min_found=1;
    bestpt.is_line_minimum=1; // Here and in the up-front check at the
    /// top of this routine should be the only two places where
    /// is_line_minimum gets set true.
  }

  assert(bracket.left.Ep<0.0
	 && (bracket.right.E>bracket.left.E || bracket.right.Ep>=0));
}

void Oxs_CGEvolve::SetBasePoint(Oxs_ConstKey<Oxs_SimState> cstate_key)
{
  const Oxs_SimState* cstate = cstate_key.GetPtr();
  if(cstate->Id() == basept.id) {
    return;  // Already set
  }

  ++cycle_count;
  ++cycle_sub_count;

  // Some convenience variables
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(cstate->Ms);
  const Oxs_MeshValue<ThreeVector>& spin = cstate->spin;
  const Oxs_Mesh* mesh = cstate->mesh;
  const OC_INDEX size = mesh->Size();
  OC_REAL8m Ep = 0.0;
  OC_REAL8m last_scaling = basept.direction_max_mag;
  OC_REAL8m last_step = bestpt.offset;
  OC_BOOL   last_step_is_minimum = bestpt.is_line_minimum;

  // cstate and bestpt should be same
  if(!bestpt.key.SameState(cstate_key)) {
    // Fill bestpt from cstate
    bestpt.key = cstate_key;
    bestpt.key.GetReadReference(); // Set read lock
    GetEnergyAndmxHxm(cstate,bestpt.energy,bestpt.mxHxm,NULL);
    last_step_is_minimum = 0;
    /// Should last_step be set to 0 in this case???
  }
  bestpt.offset=0.0; // Position relative to basept
  bestpt.is_line_minimum=0; // Presumably not minimum wrt new direction

  // Determine new direction
  if(!basept.valid
     || cstate->stage_number!=basept.stage
     || cycle_sub_count >= gradient_reset_count
     || basept.direction.Size()!=bestpt.mxHxm.Size()
     || !last_step_is_minimum) {
    basept.valid=0;
  } else {
    // Use conjugate gradient as new direction, where
    //              n = cycle count
    //              k = cell index
    //         g_n[k] = V[k]*Ms[k]*mxHxm_n[k]
    //          gamma = g_n^2/g_(n-1)^2            (Fletcher-Reeves)
    //  or      gamma = g_n*(g_n-g_(n-1))/g_(n-1)^2  (Polak-Ribiere)
    //    Direction_n = g_n+gamma*direction_(n-1).
    // Note that the g_(n-1) (or equivalently, mxHxm[n-1])
    // array is only needed for the Polak-Ribiere method.
    // If using the Fletcher-Reeves method, only the g_(n-1)^2 
    // scalar value needs to be saved between cycles, reducing
    // memory requirements.
    //   We also perform an additional correction step to Direction_n,
    // making it orthogonal to m[k] (=m_n[k]) at each k.  The wisdom
    // of this is uncertain.
    //   See NOTES IV, 27-Dec-2004, pp. 1-3 for details, also NOTES II,
    // 29-May-2002, p. 156, and NOTES III, 30-Jan-2003, p. 35.

#if REPORT_TIME_CGDEVEL
    timer[3].Start(); /**/ // asdf
#endif // REPORT_TIME_CGDEVEL

    OC_REAL8m gamma, new_g_sum_sq;
#if !OOMMF_THREADS
    OC_INDEX i;
    gamma = new_g_sum_sq = 0.0;
    if(basept.method != Basept_Data::POLAK_RIBIERE) {
      // Fletcher-Reeves method
      for(i=0;i<size;i++) {
        OC_REAL8m cell_scale = Ms[i] * mesh->Volume(i);
	new_g_sum_sq += bestpt.mxHxm[i].MagSq() * cell_scale * cell_scale;
      }
      gamma = new_g_sum_sq;
    } else {
      // Polak-Ribiere method
      for(i=0;i<size;i++) {
        OC_REAL8m cell_scale = Ms[i] * mesh->Volume(i);
        OC_REAL8m cell_scale_sq = cell_scale * cell_scale;
	ThreeVector temp = bestpt.mxHxm[i];
	new_g_sum_sq += temp.MagSq() * cell_scale_sq;
	temp -= basept.mxHxm[i];  // Polak-Ribiere adjustment
	gamma += temp * bestpt.mxHxm[i] * cell_scale_sq;
	basept.mxHxm[i] = bestpt.mxHxm[i];
      }
    }
    gamma /= basept.g_sum_sq;
#else // !OOMMF_THREADS
    {
      Oxs_ThreadTree threadtree;
      const OC_INT4m MaxThreadCount = Oc_GetMaxThreadCount();
      vector<_Oxs_CGEvolve_SetBasePoint_ThreadA> threadA;
      threadA.resize(MaxThreadCount);

      _Oxs_CGEvolve_SetBasePoint_ThreadA::job_control.Lock();
      _Oxs_CGEvolve_SetBasePoint_ThreadA::offset = 0;
      _Oxs_CGEvolve_SetBasePoint_ThreadA::job_control.Unlock();

      OC_INDEX ibs = (size + MaxThreadCount - 1)/MaxThreadCount;
      ibs = (ibs + 7)/8;
      if(ibs%8) { ibs += 8 - (ibs%8); }
      if(ibs> static_cast<OC_INDEX>((size*3+3)/4)) ibs = size;

      OC_INT4m ithread;
      for(ithread=0;ithread<MaxThreadCount;++ithread) {
        threadA[ithread].mesh = mesh;
        threadA[ithread].Ms = &Ms;
        threadA[ithread].bestpt_mxHxm = &(bestpt.mxHxm);
        threadA[ithread].basept_mxHxm = &(basept.mxHxm);
        threadA[ithread].vecsize = size;
        threadA[ithread].block_size = ibs;
        if(basept.method == Basept_Data::POLAK_RIBIERE) {
          threadA[ithread].method
            = _Oxs_CGEvolve_SetBasePoint_ThreadA::POLAK_RIBIERE;
        } else {
          threadA[ithread].method
            = _Oxs_CGEvolve_SetBasePoint_ThreadA::FLETCHER_REEVES;
        }
        if(ithread>0) threadtree.Launch(threadA[ithread],0);
      }
      threadtree.LaunchRoot(threadA[0],0);
      // Accumulate results
      gamma = new_g_sum_sq = 0.0;
      for(ithread=0;ithread<MaxThreadCount;++ithread) {
        gamma  += threadA[ithread].gamma_sum; 
        new_g_sum_sq += threadA[ithread].g_sum_sq;
      }
      gamma /= basept.g_sum_sq;
    }
#endif // !OOMMF_THREADS

#if REPORT_TIME_CGDEVEL
timer[3].Stop(); /**/
++timer_counts[3].pass_count;
if(basept.method != Basept_Data::POLAK_RIBIERE) {
  timer_counts[3].bytes += size*(sizeof(OC_REAL8m)+sizeof(ThreeVector));
} else {
  timer_counts[3].bytes += size*(sizeof(OC_REAL8m)+3*sizeof(ThreeVector));
}
#endif // REPORT_TIME_CGDEVEL


#if REPORT_TIME_CGDEVEL
timer[5].Start(); /**/ // asdf
#endif // REPORT_TIME_CGDEVEL

#if !OOMMF_THREADS
    OC_REAL8m maxmagsq = 0.0;
    OC_REAL8m normsumsq = 0.0;
    for(i=0;i<size;++i) {
      OC_REAL8m cell_scale = Ms[i] * mesh->Volume(i);
      ThreeVector temp = gamma * basept.direction[i];
      temp += cell_scale*bestpt.mxHxm[i];
      // Make temp orthogonal to spin[i].  If the angle between
      // temp and spin[i] is larger than about 45 degrees, then
      // it may be beneficial to divide by spin[i].MagSq(),
      // to offset effects of non-unit spin.  For small angles
      // it appears to be better to leave out this adjustment.
      // Or, better still, use temp -= (temp.m)m formulation
      // for large angles, w/o size correction.
      temp.wxvxw(spin[i]);

      basept.direction[i] = temp;
      Ep += (temp*bestpt.mxHxm[i])*cell_scale;
      /// See mjd's NOTES II, 29-May-2002, p156.
      OC_REAL8m magsq = temp.MagSq();
      if(magsq>maxmagsq) maxmagsq = magsq;
      normsumsq += magsq;
    }
#else // OOMMF_THREADS
    OC_REAL8m maxmagsq, normsumsq;
    {
      Oxs_ThreadTree threadtree;
      const OC_INT4m MaxThreadCount = Oc_GetMaxThreadCount();
      vector<_Oxs_CGEvolve_SetBasePoint_ThreadB> threadB;
      threadB.resize(MaxThreadCount);

      _Oxs_CGEvolve_SetBasePoint_ThreadB::job_control.Lock();
      _Oxs_CGEvolve_SetBasePoint_ThreadB::offset = 0;
      _Oxs_CGEvolve_SetBasePoint_ThreadB::job_control.Unlock();

      OC_INDEX ibs = (size + MaxThreadCount - 1)/MaxThreadCount;
      ibs = (ibs + 7)/8;
      if(ibs%8) { ibs += 8 - (ibs%8); }
      if(ibs> static_cast<OC_INDEX>((size*3+3)/4)) ibs = size;

      OC_INT4m ithread;
      for(ithread=0;ithread<MaxThreadCount;++ithread) {
        threadB[ithread].mesh = mesh;
        threadB[ithread].Ms = &Ms;
        threadB[ithread].spin = &spin;
        threadB[ithread].bestpt_mxHxm = &(bestpt.mxHxm);
        threadB[ithread].basept_direction = &(basept.direction);
        threadB[ithread].gamma = gamma;
        threadB[ithread].vecsize = size;
        threadB[ithread].block_size = ibs;
        if(ithread>0) threadtree.Launch(threadB[ithread],0);
      }
      threadtree.LaunchRoot(threadB[0],0);
      // Accumulate results
      Ep = maxmagsq = normsumsq = 0.0;
      for(ithread=0;ithread<MaxThreadCount;++ithread) {
        if(threadB[ithread].maxmagsq > maxmagsq) {
          maxmagsq = threadB[ithread].maxmagsq;
        }
        Ep  += threadB[ithread].Ep; 
        normsumsq += threadB[ithread].normsumsq; 
      }
    }
#endif // OOMMF_THREADS


#if REPORT_TIME_CGDEVEL
timer[5].Stop(); /**/
++timer_counts[5].pass_count;
 timer_counts[5].bytes += size*(sizeof(OC_REAL8m)+4*sizeof(ThreeVector));
#endif // REPORT_TIME_CGDEVEL

    basept.direction_max_mag = sqrt(maxmagsq);
    basept.direction_norm = sqrt(normsumsq);
    basept.g_sum_sq = new_g_sum_sq;
    OC_REAL8m direction_slope = Ep;
    Ep *= -MU0; // See mjd's NOTES II, 29-May-2002, p156.
    bestpt.Ep = Ep;
    bestpt.grad_norm = sqrt(new_g_sum_sq);

    // Ep is derivative of energy in basept.direction.  If this
    // is not negative (presumably because the previous line
    // minimization did a poor job), then pitch current CG
    // data and restart using the pure gradient as the search
    // direction.  (See next code block.)
    basept.valid=0; // Use gradient unless following test is passed.
    const OC_REAL8m maxdot = basept.direction_norm*bestpt.grad_norm;
    if(direction_slope < -maxdot*sqrt(double(size))*OC_REAL8_EPSILON*8) {
      static Oxs_WarningMessage foo(3);
      foo.Send(revision_info,OC_STRINGIFY(__LINE__),
               "Bad line direction selection; poor line minimization?"
               "  If convergence is slow, try reducing the"
               " line_minimum_angle_precision setting in the"
               " Oxs_CGEvolve Specify block.");
    } else if(direction_slope > gradient_reset_angle_cos*maxdot) {
      // New direction slope is not small compared to what it would be
      // if pure gradient direction were used.  Interpret failure of
      // this test as a sign that the conjugate-gradient algorithm has
      // gone astray, in which case we leave basept.valid=0 and reset
      // direction to pure gradient (next code block).
      basept.valid=1;
    }
  }

  if(!basept.valid) {
#if REPORT_TIME_CGDEVEL
timer[6].Start(); /**/
#endif // REPORT_TIME_CGDEVEL
    // Use gradient as new direction
    basept.direction.AdjustSize(mesh);
    cycle_sub_count=0;
    if(basept.method == Basept_Data::POLAK_RIBIERE) {
      basept.mxHxm = bestpt.mxHxm;      
    }
#if !OOMMF_THREADS
    OC_REAL8m maxmagsq = 0.0;
    OC_REAL8m sumsq = 0.0;
    for(OC_INDEX i=0;i<size;i++) {
      OC_REAL8m cell_scale = Ms[i] * mesh->Volume(i);
      basept.direction[i] = bestpt.mxHxm[i];
      basept.direction[i] *= cell_scale;
      OC_REAL8m magsq = basept.direction[i].MagSq();
      if(magsq>maxmagsq) maxmagsq = magsq;
      sumsq += magsq;
      /// See mjd's NOTES II, 29-May-2002, p156.
    }
    Ep = sumsq;
#else // !OOMMF_THREADS
    OC_REAL8m maxmagsq, sumsq;
    {
      Oxs_ThreadTree threadtree;
      const OC_INT4m MaxThreadCount = Oc_GetMaxThreadCount();
      vector<_Oxs_CGEvolve_SetBasePoint_ThreadC> threadC;
      threadC.resize(MaxThreadCount);

      _Oxs_CGEvolve_SetBasePoint_ThreadC::job_control.Lock();
      _Oxs_CGEvolve_SetBasePoint_ThreadC::offset = 0;
      _Oxs_CGEvolve_SetBasePoint_ThreadC::job_control.Unlock();

      OC_INDEX ibs = (size + MaxThreadCount - 1)/MaxThreadCount;
      ibs = (ibs + 7)/8;
      if(ibs%8) { ibs += 8 - (ibs%8); }
      if(ibs> static_cast<OC_INDEX>((size*3+3)/4)) ibs = size;

      OC_INT4m ithread;
      for(ithread=0;ithread<MaxThreadCount;++ithread) {
        threadC[ithread].mesh = mesh;
        threadC[ithread].Ms = &Ms;
        threadC[ithread].bestpt_mxHxm = &(bestpt.mxHxm);
        threadC[ithread].basept_direction = &(basept.direction);
        threadC[ithread].vecsize = size;
        threadC[ithread].block_size = ibs;
        if(ithread>0) threadtree.Launch(threadC[ithread],0);
      }
      threadtree.LaunchRoot(threadC[0],0);
      // Accumulate results
      maxmagsq = sumsq = 0.0;
      for(ithread=0;ithread<MaxThreadCount;++ithread) {
        if(threadC[ithread].maxmagsq > maxmagsq) {
          maxmagsq = threadC[ithread].maxmagsq;
        }
        sumsq += threadC[ithread].sumsq;
      }
      Ep = sumsq;
    }
#endif // !OOMMF_THREADS
    basept.direction_max_mag = sqrt(maxmagsq);
    basept.g_sum_sq = sumsq;
    basept.direction_norm = sqrt(sumsq);
    Ep *= -MU0; // See mjd's NOTES II, 29-May-2002, p156.
    bestpt.Ep = Ep;
    bestpt.grad_norm = basept.direction_norm;
    basept.valid=1;
#if REPORT_TIME_CGDEVEL
timer[6].Stop(); /**/
++timer_counts[6].pass_count;
 timer_counts[6].bytes += size*(sizeof(OC_REAL8m)+2*sizeof(ThreeVector));
#endif // REPORT_TIME_CGDEVEL
  }

  // Fill remaining basept fields
  basept.id = cstate->Id();
  basept.stage = cstate->stage_number;
  if(!cstate->GetDerivedData("Total energy",basept.total_energy)) {
      throw Oxs_ExtError(this,
        "Missing \"Total energy\" data in Oxs_CGEvolve::SetBasePoint()."
	" Programming error?");
  }

  // Initialize bracket data
  bracket.min_bracketed=0;
  bracket.min_found=0;

  if(basept.direction_max_mag>=1.0
     || bracket.maxstep<basept.direction_max_mag*DBL_MAX) {
    bracket.scaled_minstep=bracket.minstep/basept.direction_max_mag;
    bracket.scaled_maxstep=bracket.maxstep/basept.direction_max_mag;
  } else {
    // Safety
    if(bracket.maxstep>0.0) {
      bracket.scaled_maxstep = 0.5*DBL_MAX;
      bracket.scaled_minstep
	= bracket.scaled_maxstep*(bracket.minstep/bracket.maxstep);
    } else {
      bracket.scaled_maxstep = 0.0;
      bracket.scaled_minstep = 0.0;
    }
  }

  // Size initial step to match that from previous
  // line minimization, if applicable
  bracket.start_step = bracket.scaled_minstep;
  if(last_step_is_minimum && last_step>0.0) {
    OC_REAL8m scaling_ratio = 1.0;
    if(basept.direction_max_mag>=1.0
       || last_scaling<basept.direction_max_mag*DBL_MAX) {
      scaling_ratio = last_scaling/basept.direction_max_mag;
    }
    if(scaling_ratio<=1.0
       || last_step<DBL_MAX/scaling_ratio) {
      bracket.start_step = last_step*scaling_ratio;
    }
    if(bracket.start_step>bracket.scaled_maxstep) {
      bracket.start_step = bracket.scaled_maxstep;
    }
  }

  // To force maximal restriction on initial bracket interval
  // reduction (see FindLineMinimumStep()), set both
  // minimum_reduction_ratio's to 1.  Moderately conservative
  // values are last_mrr=0.5, next_to_last_mrr=0.25.  Or set
  // both values to 0 to allow FindLineMinimumStep complete
  // freedom on first two steps.
  bracket.last_min_reduction_ratio=1./16.;
  bracket.next_to_last_min_reduction_ratio=1./256.;
  bracket.left.offset = bracket.right.offset = 0.;
  bracket.left.E      = bracket.right.E      = 0.;
  bracket.left.Ep     = bracket.right.Ep     = Ep;
}

void Oxs_CGEvolve::NudgeBestpt(const Oxs_SimState* oldstate,
			       Oxs_Key<Oxs_SimState>& statekey)
{ // Line minimization didn't change state, which is
  // a problem because in this case the next pass of
  // the base-direction code may re-generate the same
  // direction, in which case the algorithm gets
  // stuck.  So, we nudge.
  assert(bracket.left.Ep<0.0
	 && (bracket.right.E>bracket.left.E || bracket.right.Ep>=0));
  OC_REAL8m nudge = DBL_MAX/2;
  if(basept.direction_max_mag>1.
     || 4*OC_REAL8_EPSILON<nudge*basept.direction_max_mag) {
    nudge = 4*OC_REAL8_EPSILON/basept.direction_max_mag;
    /// nudge is a reasonable guess at the smallest variation in
    /// offset that leads to a detectable (i.e., discretizational)
    /// change in spins.
  } else {
    // Degenerate case.  basept.direction_max_mag must be exactly
    // or very nearly zero.  Punt.
    nudge = 1.0;
  }

  OC_REAL8m test_offset = 0.5;
  if(bracket.right.Ep>0.) {
    if(-bracket.left.Ep
       < test_offset*(bracket.right.Ep-bracket.left.Ep)) {
      test_offset
	= -bracket.left.Ep/(bracket.right.Ep-bracket.left.Ep);
    }
  } else {
    // Ep data is not informative
    test_offset = 0.0;  // This will get bumped up by
    /// "nudge" test below.
  }
  test_offset *= bracket.right.offset/2.0;  // Scale best guess of
  // minima to interval, and then reduce to err on the safe side.
  
  // Make sure offset can be felt.
  if(test_offset<nudge) test_offset=nudge;
  ++line_minimum_count;
  Bracket testpt;
  FillBracket(test_offset,oldstate,statekey,testpt);
  UpdateBrackets(statekey,testpt,1);
  // One may want to force testpt acceptance only if
  // testpt.E<energy_slack, but if this doesn't hold
  // then one needs a backup plan.
}


OC_BOOL Oxs_CGEvolve::InitNewStage
(const Oxs_MinDriver* /* driver */,
 Oxs_ConstKey<Oxs_SimState> state,
 Oxs_ConstKey<Oxs_SimState> /* prevstate */)
{
  SetBasePoint(state);
  return 1;
}


OC_BOOL Oxs_CGEvolve::Step(const Oxs_MinDriver* driver,
                        Oxs_ConstKey<Oxs_SimState> current_state_key,
                        const Oxs_DriverStepInfo& /* step_info */,
                        Oxs_Key<Oxs_SimState>& next_state_key)
{
  // Note: On entry, next_state_key is holding a write lock.
#ifdef INSTRUMENT
  static OC_REAL8m last_energy = 0.0;
#endif // INSTRUMENT

#if REPORT_TIME
  steponlytime.Start();
#endif // REPORT_TIME

  const Oxs_SimState* cstate = &(current_state_key.GetReadReference());

  if(++step_attempt_count == 1) {
    // Update counts from current state, if available.  This keeps
    // counts in-sync with data saved in checkpoint (restart) files.
    // Wrap-around on step_attempt_count could cause some loss here,
    // but then we would also be wrapping around on the *count's
    // below; I'm not going to lose sleep over it.  This is rather
    // hackish.  It seems like it would be better to get these
    // values set up inside the Init() function.
    OC_REAL8m temp_value;
    if(cstate->GetDerivedData("Energy calc count",temp_value)) {
      energy_calc_count = static_cast<OC_UINT4m>(temp_value);
    }
    if(cstate->GetDerivedData("Cycle count",temp_value)) {
      cycle_count = static_cast<OC_UINT4m>(temp_value);
    }
    if(cstate->GetDerivedData("Cycle sub count",temp_value)) {
      cycle_sub_count = static_cast<OC_UINT4m>(temp_value);
    }
    if(cstate->GetDerivedData("Bracket count",temp_value)) {
      bracket_count = static_cast<OC_UINT4m>(temp_value);
    }
    if(cstate->GetDerivedData("Line min count",temp_value)) {
      line_minimum_count = static_cast<OC_UINT4m>(temp_value);
    }
  }

  // Do first part of next_state structure initialization.
  driver->FillStateMemberData(*cstate,
                              next_state_key.GetWriteReference());

  if(!basept.valid
     || basept.stage != cstate->stage_number
     || bracket.min_found) {
#ifdef INSTRUMENT
/**/ OC_REAL8m new_energy=0.0;
/**/ cstate->GetDerivedData("Total energy",new_energy);
/**/ fprintf(stderr," Energy drop: %g\n",double(last_energy-new_energy));
/**/ last_energy=new_energy;
/**/ fprintf(stderr,"It=%5d  valid=%d, stages: %d/%d, min_found=%d\n",
/**/  cstate->iteration_count+1,
/**/  basept.valid,basept.stage,cstate->stage_number,bracket.min_found);
#endif // INSTRUMENT
#if REPORT_TIME_CGDEVEL
  basepttime.Start();
#endif // REPORT_TIME_CGDEVEL
    SetBasePoint(current_state_key);
#if REPORT_TIME_CGDEVEL
  basepttime.Stop();
#endif // REPORT_TIME_CGDEVEL
  }

  if(!bracket.min_bracketed) {
#if REPORT_TIME_CGDEVEL
  findbrackettime.Start();
#endif // REPORT_TIME_CGDEVEL
    FindBracketStep(cstate,next_state_key);
#if REPORT_TIME_CGDEVEL
  findbrackettime.Stop();
#endif // REPORT_TIME_CGDEVEL
  } else if(!bracket.min_found) {
#if REPORT_TIME_CGDEVEL
  findlinemintime.Start();
#endif // REPORT_TIME_CGDEVEL
    FindLineMinimumStep(cstate,next_state_key);
    if(bracket.min_found && bestpt.offset==0.0) {
      NudgeBestpt(cstate,next_state_key);
    }
#if REPORT_TIME_CGDEVEL
  findlinemintime.Stop();
#endif // REPORT_TIME_CGDEVEL
  }

  // Note: Call to GetReadReference insures that write lock is released.
  const Oxs_SimState* nstate = &(next_state_key.GetReadReference());

  // Finish filling in next_state structure 
  driver->FillStateDerivedData(*cstate,*nstate);

#if defined(TRACK_STEPSIZE)
  {
    OC_REAL8m mindot = 1.0;
    const OC_INDEX meshsize = cstate->mesh->Size();
    const Oxs_MeshValue<ThreeVector>& cspin = cstate->spin;
    const Oxs_MeshValue<ThreeVector>& nspin = nstate->spin;
    for(OC_INDEX i=0;i<meshsize;++i) {
      OC_REAL8m dot = cspin[i]*nspin[i];
      if(dot<mindot) mindot = dot;
    }
    OC_REAL8m maxang = acos(mindot)*180./PI;
    OC_REAL8m maxang_allowed = atan(bracket.maxstep)*180./PI;
    if(maxang>maxang_allowed) {
      fprintf(stderr,"MAX ANGLE VIOLATION: %g > %g\n",
              static_cast<double>(maxang),
              static_cast<double>(maxang_allowed));
      fprintf(stderr," cstate: %u/%u/%u   nstate: %u/%u/%u\n",
              cstate->Id(),cstate->iteration_count,
              cstate->stage_iteration_count,
              nstate->Id(),nstate->iteration_count,
              nstate->stage_iteration_count);
    }
  }
#endif // TRACK_STEPSIZE

#if REPORT_TIME
  steponlytime.Stop();
#endif // REPORT_TIME
  return bestpt.key.SameState(next_state_key);
}


void Oxs_CGEvolve::UpdateDerivedFieldOutputs(const Oxs_SimState& state)
{ // Fill Oxs_VectorOutput's that have CacheRequest enabled.
  if(total_H_field_output.GetCacheRequestCount()>0
     || state.Id()!=temp_id) {
    // Need to call GetEnergyAndmxHxm
    temp_id=0;
    if(total_H_field_output.GetCacheRequestCount()>0) {
      total_H_field_output.cache.state_id=0;
      GetEnergyAndmxHxm(&state,temp_energy,temp_mxHxm,
			&total_H_field_output.cache.value);
      total_H_field_output.cache.state_id=state.Id();
    } else {
      GetEnergyAndmxHxm(&state,temp_energy,temp_mxHxm,NULL);
    }
    temp_id = state.Id();
  }
  if(total_energy_density_output.GetCacheRequestCount()>0) {
    total_energy_density_output.cache.state_id = 0;
    total_energy_density_output.cache.value = temp_energy;
    total_energy_density_output.cache.state_id = state.Id();
  }
  if(mxHxm_output.GetCacheRequestCount()>0) {
    mxHxm_output.cache.state_id = 0;
    mxHxm_output.cache.value    = temp_mxHxm;
    mxHxm_output.cache.state_id = temp_id;
  }
}

void Oxs_CGEvolve::UpdateDerivedOutputs(const Oxs_SimState& state)
{ // This routine fills all the Oxs_CGEvolve Oxs_ScalarOutput's to
  // the appropriate value based on the import "state".

  energy_calc_count_output.cache.state_id
    = cycle_count_output.cache.state_id
    = cycle_sub_count_output.cache.state_id
    = max_mxHxm_output.cache.state_id
    = delta_E_output.cache.state_id
    = total_energy_output.cache.state_id
    = bracket_count_output.cache.state_id
    = line_min_count_output.cache.state_id
    = 0;  // Mark change in progress

  OC_REAL8m last_energy;
  if(!state.GetDerivedData("Energy calc count",
			   energy_calc_count_output.cache.value) ||
     !state.GetDerivedData("Cycle count",
			   cycle_count_output.cache.value) ||
     !state.GetDerivedData("Cycle sub count",
			   cycle_sub_count_output.cache.value) ||
     !state.GetDerivedData("Max mxHxm",max_mxHxm_output.cache.value) ||
     !state.GetDerivedData("Last energy",last_energy) ||
     !state.GetDerivedData("Total energy",
			   total_energy_output.cache.value) ||
     !state.GetDerivedData("Bracket count",
			   bracket_count_output.cache.value) ||
     !state.GetDerivedData("Line min count",
			   line_min_count_output.cache.value)) {
    // Missing at least some data
    temp_id=0;
    GetEnergyAndmxHxm(&state,temp_energy,temp_mxHxm,NULL);
    temp_id=state.Id();
    if(!state.GetDerivedData("Energy calc count",
			     energy_calc_count_output.cache.value) ||
       !state.GetDerivedData("Cycle count",
			     cycle_count_output.cache.value) ||
       !state.GetDerivedData("Cycle sub count",
			     cycle_sub_count_output.cache.value) ||
       !state.GetDerivedData("Max mxHxm",max_mxHxm_output.cache.value) ||
       !state.GetDerivedData("Last energy",last_energy) ||
       !state.GetDerivedData("Total energy",
			     total_energy_output.cache.value) ||
       !state.GetDerivedData("Bracket count",
			     bracket_count_output.cache.value) ||
       !state.GetDerivedData("Line min count",
			     line_min_count_output.cache.value)) {
      throw Oxs_ExtError(this,"Missing output data."
			   " Programming error?");
    }
  }
  delta_E_output.cache.value
    = total_energy_output.cache.value - last_energy;

  energy_calc_count_output.cache.state_id
    = cycle_count_output.cache.state_id
    = cycle_sub_count_output.cache.state_id
    = max_mxHxm_output.cache.state_id
    = delta_E_output.cache.state_id
    = total_energy_output.cache.state_id
    = bracket_count_output.cache.state_id
    = line_min_count_output.cache.state_id
    = state.Id();
}
