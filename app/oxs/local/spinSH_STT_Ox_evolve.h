/* FILE: spinSH_STT_Ox_evolve.cc                 -*-Mode: c++-*-
 *
 * Concrete evolver class, built on top of the Oxs_RungeKuttaEvolve
 * class, but including a spin injection field torque.
 *
 */

#ifndef _OXS_SPIN_SH_STT_OX_EVOLVE
#define _OXS_SPIN_SH_STT_OX_EVOLVE

#include "threevector.h"
#include "util.h"
#include "mesh.h"
#include "meshvalue.h"
#include "scalarfield.h"
#include "vectorfield.h"

#include "key.h"
#include "output.h"
#include "tclcommand.h"

#include "timeevolver.h"

/* End includes */
class Oxs_Spin_SH_STT_Ox_Evolve:public Oxs_TimeEvolver {
private:
  mutable OC_UINT4m mesh_id;

  // Spin-transfer polarization field data structures for Slonczewski metal barrier model
  Oxs_OwnedPointer<Oxs_ScalarField> P_init;
  Oxs_OwnedPointer<Oxs_ScalarField> J_init;
  Oxs_OwnedPointer<Oxs_VectorField> mp_init;
//Zeng DOI:10.1038/srep01426, Sankey DOI:10.1038/nphys783 --> field torque boost under voltage load

  mutable Oxs_MeshValue<OC_REAL8m> polarization;
  mutable Oxs_MeshValue<OC_REAL8m> A;
  mutable Oxs_MeshValue<OC_REAL8m> A2;
  mutable Oxs_MeshValue<OC_REAL8m> B;
  mutable Oxs_MeshValue<OC_REAL8m> B2;
  mutable Oxs_MeshValue<ThreeVector> mp;
//Zeng DOI:10.1038/srep01426, Sankey DOI:10.1038/nphys783 --> field torque boost under voltage load
  OC_REAL8m ave_J;   // Average value of J, excluding J_profile term.
  OC_REAL8m ave_ICSH; // Average value of JSH , excluding JSH_profile term.

  // Current flow direction.  Computed "free layer thickness" is
  // actually just cell size dimension parallel to current flow.
  // Current flow direction is also used to propagate mp, if enabled.
  enum JDirection { JD_INVALID,
                    JD_X_POS, JD_X_NEG,
                    JD_Y_POS, JD_Y_NEG,
                    JD_Z_POS, JD_Z_NEG };
  JDirection J_direction;

  void UpdateMeshArrays(const Oxs_Mesh*);

  void Propagate_mp(const Oxs_Mesh* mesh,
                    const Oxs_MeshValue<ThreeVector>& spin);

  // Support for time varying current:
  OC_BOOL has_J_profile;
  vector<Oxs_TclCommandLineOption> J_profile_opts;
  Oxs_TclCommand J_profile_cmd;
  OC_REAL8m EvaluateJProfileScript(OC_UINT4m stage,OC_REAL8m stage_time,
                                OC_REAL8m total_time) const;

  // Support for time varying Spin-Hall current:
  OC_BOOL has_ICSH_profile;
  vector<Oxs_TclCommandLineOption> ICSH_profile_opts;
  Oxs_TclCommand ICSH_profile_cmd;
  OC_REAL8m EvaluateICSHProfileScript(OC_UINT4m stage,OC_REAL8m stage_time,
                                OC_REAL8m total_time) const;

  // Base step size control parameters
  OC_REAL8m min_timestep;           // Seconds
  OC_REAL8m max_timestep;           // Seconds

  const OC_REAL8m max_step_decrease;        // Safety size adjusment
  const OC_REAL8m max_step_increase_limit;  // bounds.
  const OC_REAL8m max_step_increase_adj_ratio;
  OC_REAL8m max_step_increase;
  /// NOTE: These bounds do not include step_headroom, which
  /// is applied at the end.

  // Error-based step size control parameters.  Each may be disabled
  // by setting to -1.  There is an additional step size control that
  // insures that energy is monotonically non-increasing (up to
  // estimated rounding error).
  OC_REAL8m allowed_error_rate;  // Step size is adjusted so
  /// that the estimated maximum error (across all spins) divided
  /// by the step size is smaller than this value.  The units
  /// internally are radians per second, converted from the value
  /// specified in the input MIF file, which is in deg/sec.

  OC_REAL8m allowed_absolute_step_error; // Similar to allowed_error_rate,
  /// but without the step size adjustment.  Internal units are
  /// radians; MIF input units are degrees.

  OC_REAL8m allowed_relative_step_error; // Step size is adjusted so that
  /// the estimated maximum error (across all spins) divided by
  /// [maximum dm/dt (across all spins) * step size] is smaller than
  /// this value.  This value is non-dimensional, representing the
  /// allowed relative (proportional) error, presumably in (0,1).

  OC_REAL8m expected_energy_precision; // Expected relative energy
  /// precision.  Set <0 to disable energy stepsize control.

  OC_REAL8m energy_check_slack; // Allowed slack in energy stepsize
  /// control.


  OC_REAL8m reject_goal,reject_ratio;
  OC_REAL8m min_step_headroom,max_step_headroom;
  OC_REAL8m step_headroom; // Safety margin used in step size adjustment

  // Spatially variable Landau-Lifschitz-Gilbert damping coef
  Oxs_OwnedPointer<Oxs_ScalarField> alpha_init;
  mutable Oxs_MeshValue<OC_REAL8m> alpha;

  // Spatially variable Landau-Lifschitz-Gilbert gyromagnetic ratio.
  //OC_BOOL do_precess;  // If false, then do pure damping
  OC_BOOL allow_signed_gamma; // If false, then force gamma negative
  enum GammaStyle { GS_INVALID, GS_LL, GS_G }; // Landau-Lifshitz or Gilbert
  GammaStyle gamma_style;
  Oxs_OwnedPointer<Oxs_ScalarField> gamma_init;
  mutable Oxs_MeshValue<OC_REAL8m> gamma;

  Oxs_OwnedPointer<Oxs_ScalarField> g_init;
  mutable Oxs_MeshValue<OC_REAL8m> g;

  // Spin-transfer data structures due to spin Hall effect
  Oxs_OwnedPointer<Oxs_VectorField> ico_init;
  mutable Oxs_MeshValue<ThreeVector> ico;
  Oxs_OwnedPointer<Oxs_VectorField> iso_init;
  mutable Oxs_MeshValue<ThreeVector> iso;
  Oxs_OwnedPointer<Oxs_ScalarField> thetaSH_init;
  mutable Oxs_MeshValue<OC_REAL8m>  thetaSH;
  Oxs_OwnedPointer<Oxs_ScalarField> tSH_init;
  mutable Oxs_MeshValue<OC_REAL8m>  tSH;
  Oxs_OwnedPointer<Oxs_ScalarField> wSH_init;
  mutable Oxs_MeshValue<OC_REAL8m>  wSH;
  Oxs_OwnedPointer<Oxs_ScalarField> lambdaSF_init;
  mutable Oxs_MeshValue<OC_REAL8m> lambdaSF;
  Oxs_OwnedPointer<Oxs_ScalarField> Icsh_init;
//  mutable Oxs_MeshValue<OC_REAL8m> icsh;// covered by A2
  Oxs_OwnedPointer<Oxs_ScalarField> beta_init;
  mutable Oxs_MeshValue<OC_REAL8m> beta;

  // The next timestep is based on the error from the last step.  If
  // there is no last step (either because this is the first step,
  // or because the last state handled by this routine is different
  // from the incoming current_state), then timestep is calculated
  // so that max_dm_dt * timestep = start_dm.
  OC_REAL8m start_dm;

  // Stepsize control for first step of each stage after the first.
  // Choices are to use start_dm, use continuation from end of
  // previous stage, or to automatically select between the two
  // methods depending on whether or not the energy appears to
  // be continuous across the stage boundary.
  enum StageInitStepControl { SISC_INVALID, SISC_START_DM,
			      SISC_CONTINUOUS, SISC_AUTO };
  StageInitStepControl stage_init_step_control;

  // Data cached from last state
  OC_UINT4m energy_state_id;
  Oxs_MeshValue<OC_REAL8m> energy;
  OC_REAL8m next_timestep;

  // Outputs
  void UpdateDerivedOutputs(const Oxs_SimState&);
  void UpdateSpinTorqueOutputs(const Oxs_SimState&);
  Oxs_ScalarOutput<Oxs_Spin_SH_STT_Ox_Evolve> max_dm_dt_output;
  Oxs_ScalarOutput<Oxs_Spin_SH_STT_Ox_Evolve> dE_dt_output;
  Oxs_ScalarOutput<Oxs_Spin_SH_STT_Ox_Evolve> delta_E_output;
  Oxs_ScalarOutput<Oxs_Spin_SH_STT_Ox_Evolve> ave_J_output;
  Oxs_ScalarOutput<Oxs_Spin_SH_STT_Ox_Evolve> ave_ICSH_output;
  Oxs_VectorFieldOutput<Oxs_Spin_SH_STT_Ox_Evolve> dm_dt_output;
  Oxs_VectorFieldOutput<Oxs_Spin_SH_STT_Ox_Evolve> mxH_output;
  Oxs_VectorFieldOutput<Oxs_Spin_SH_STT_Ox_Evolve> spin_torque_output;
  Oxs_VectorFieldOutput<Oxs_Spin_SH_STT_Ox_Evolve> Jmp_output;

  // Scratch space
  Oxs_MeshValue<OC_REAL8m> temp_energy;
  Oxs_MeshValue<OC_REAL8m> stmpA;
  Oxs_MeshValue<ThreeVector> vtmpA;
  Oxs_MeshValue<ThreeVector> vtmpB;
  Oxs_MeshValue<ThreeVector> vtmpC;
  Oxs_MeshValue<ThreeVector> vtmpD;

  // Utility functions
  void CheckCache(const Oxs_SimState& cstate);

  void AdjustState(OC_REAL8m hstep,
		   OC_REAL8m mstep,
		   const Oxs_SimState& old_state,
		   const Oxs_MeshValue<ThreeVector>& dm_dt,
		   Oxs_SimState& new_state,
		   OC_REAL8m& norm_error) const;
  // Export new state has time index from old_state + h,
  // and spins from old state + mstep*dm_dt and re-normalized.

  void UpdateTimeFields(const Oxs_SimState& cstate,
			Oxs_SimState& nstate,
			OC_REAL8m stepsize) const;

  void NegotiateTimeStep(const Oxs_TimeDriver* driver,
			 const Oxs_SimState&  cstate,
			 Oxs_SimState& nstate,
			 OC_REAL8m stepsize,
			 OC_BOOL use_start_dm,
			 OC_BOOL& forcestep,
			 OC_BOOL& driver_set_step) const;

  OC_BOOL CheckError(OC_REAL8m global_error_order,OC_REAL8m error,
		  OC_REAL8m stepsize,OC_REAL8m reference_stepsize,
		  OC_REAL8m max_dm_dt,OC_REAL8m& new_stepsize);
  /// Returns 1 if step is good, 0 if error is too large.
  /// Export new_stepsize is set to suggested stepsize
  /// for next step.

  OC_REAL8m MaxDiff(const Oxs_MeshValue<ThreeVector>& vecA,
		 const Oxs_MeshValue<ThreeVector>& vecB);
  /// Returns maximum difference between vectors in corresponding
  /// positions in two vector fields.

  void AdjustStepHeadroom(OC_INT4m step_reject);
  /// step_reject should be 0 or 1, reflecting whether the current
  /// step was rejected or not.  This routine updates reject_ratio
  /// and adjusts step_headroom appropriately.

  // Stepper routines:  If routine needs to compute the energy
  // at the new (final) state, then it should store the final
  // energy results in temp_energy, mxH in mxH_output.cache,
  // and dm_dt into the vtmpA scratch array, fill
  // the "Timestep lower bound", "Max dm/dt", "dE/dt", and
  // "pE/pt" derived data fields in nstate, and set the export
  // value new_energy_and_dmdt_computed true.  Otherwise the export
  // value should be set false, and the client routine is responsible
  // for obtaining these values as necessary.  (If possible, it is
  // better to let the client compute these values, because the
  // client may be able to defer computation until it has decided
  // whether or not to keep the step.)

  // One would like to declare the step functions and pointer
  // to same via typedef's, but the MS VC++ 6.0 (& others?)
  // compiler doesn't handle member function typedef's properly---
  // it produces __cdecl linkage rather than instance member
  // linkage.  Typedef's on pointers to member functions work
  // okay, just not typedef's on member functions themselves.
  // So, instead we use a #define, which is ugly but portable.
#define RKStepFuncSig(NAME) \
  void NAME (                                            \
     OC_REAL8m stepsize,                                    \
     Oxs_ConstKey<Oxs_SimState> current_state,           \
     const Oxs_MeshValue<ThreeVector>& current_dm_dt,    \
     Oxs_Key<Oxs_SimState>& next_state,                  \
     OC_REAL8m& error_estimate,                             \
     OC_REAL8m& global_error_order,                         \
     OC_REAL8m& norm_error,                                 \
     OC_REAL8m& min_dE_dt, OC_REAL8m& max_dE_dt,               \
     OC_BOOL& new_energy_and_dmdt_computed)

  // Functions that calculate a single RK step
  RKStepFuncSig(TakeRungeKuttaStep2);
  RKStepFuncSig(TakeRungeKuttaStep2Heun);
  RKStepFuncSig(TakeRungeKuttaStep4);
  RKStepFuncSig(TakeRungeKuttaFehlbergStep54);
  RKStepFuncSig(TakeRungeKuttaFehlbergStep54M);
  RKStepFuncSig(TakeRungeKuttaFehlbergStep54S);

  // Pointer set at runtime during instance initialization
  // to one of the above functions single RK step functions.
  RKStepFuncSig((Oxs_Spin_SH_STT_Ox_Evolve::* rkstep_ptr));

  // Utility code used by the TakeRungeKuttaFehlbergStep54* routines.
  enum RKF_SubType { RKF_INVALID, RK547FC, RK547FM, RK547FS };
  void RungeKuttaFehlbergBase54(RKF_SubType method,
			   OC_REAL8m stepsize,
			   Oxs_ConstKey<Oxs_SimState> current_state,
			   const Oxs_MeshValue<ThreeVector>& current_dm_dt,
			   Oxs_Key<Oxs_SimState>& next_state,
			   OC_REAL8m& error_estimate,
			   OC_REAL8m& global_error_order,
			   OC_REAL8m& norm_error,
                           OC_REAL8m& min_dE_dt, OC_REAL8m& max_dE_dt,
			   OC_BOOL& new_energy_and_dmdt_computed);

  static OC_REAL8m PositiveTimestepBound(OC_REAL8m max_dm_dt);
  // Computes estimate on minimal timestep that will move at least one
  // spin an amount perceptible to the floating point representation.

  void Calculate_dm_dt
  (const Oxs_SimState& state_,
   const Oxs_MeshValue<ThreeVector>& mxH_,
   OC_REAL8m pE_pt_,
   Oxs_MeshValue<ThreeVector>& dm_dt_,
   OC_REAL8m& max_dm_dt_,OC_REAL8m& dE_dt_,OC_REAL8m& min_timestep_);
  /// Imports: mesh_, Ms_, mxH_, spin_, pE_pt
  /// Exports: dm_dt_, max_dm_dt_, dE_dt_, min_timestep_

  // Declare but leave undefined copy constructor and assignment operator
  Oxs_Spin_SH_STT_Ox_Evolve(const Oxs_Spin_SH_STT_Ox_Evolve&);
  Oxs_Spin_SH_STT_Ox_Evolve& operator=(const Oxs_Spin_SH_STT_Ox_Evolve&);

public:
  virtual const char* ClassName() const; // ClassName() is
  /// automatically generated by the OXS_EXT_REGISTER macro.
  virtual OC_BOOL Init();
  Oxs_Spin_SH_STT_Ox_Evolve(const char* name,     // Child instance id
		       Oxs_Director* newdtr, // App director
		       const char* argstr);  // MIF input block parameters
  virtual ~Oxs_Spin_SH_STT_Ox_Evolve();

  virtual  OC_BOOL
  Step(const Oxs_TimeDriver* driver,
       Oxs_ConstKey<Oxs_SimState> current_state,
       const Oxs_DriverStepInfo& step_info,
       Oxs_Key<Oxs_SimState>& next_state);
  // Returns true if step was successful, false if
  // unable to step as requested.
};

#endif // _OXS_SPIN_SH_STT_OX_EVOLVE
