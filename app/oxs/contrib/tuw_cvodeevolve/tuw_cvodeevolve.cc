/* FILE: tuw_cvodeevolve.cc                 -*-Mode: c++-*-
 *
 * Concrete evolver class, using cvode package
 *
 */

#include <float.h>

#include "nb.h"
#include "director.h"
#include "timedriver.h"
#include "simstate.h"
#include "key.h"
#include "energy.h"		// Needed to make MSVC++ 5 happy

#include "tuw_cvodeevolve.h"

// Oxs_Ext registration support
OXS_EXT_REGISTER(Tuw_CvodeEvolve);

/* End includes */

#include "cvspgmr.h"

// Constructor
Tuw_CvodeEvolve::Tuw_CvodeEvolve(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Oxs_TimeEvolver(name,newdtr,argstr),
    style(NONSTIFF), logfile(NULL),
    state_id(0), stage_number(0),
    stepsize(0.),
    cvode_initialized(0),
    cvode_reset_request(0), cvode_reset_count(0),
    cvode_mem(NULL), yout(NULL)
{
  datablock.evolveptr = this;
  datablock.mif_interp = director->GetMifInterp();

  // Process arguments
  min_timestep=GetRealInitValue("min_timestep",0.);
  max_timestep=GetRealInitValue("max_timestep",1e-9);
  if(max_timestep<=0.0) {
    char buf[4096];
    Oc_Snprintf(buf,sizeof(buf),
		"Invalid parameter value:"
		" Specified max time step is %g (should be >0.)",
		max_timestep);
    throw Oxs_Ext::Error(this,buf);
  }

  initial_timestep=GetRealInitValue("initial_timestep",1e-12);
  if(initial_timestep>max_timestep) initial_timestep=max_timestep;

  allowed_error = GetRealInitValue("error",1.0);
  if(allowed_error<=0.0) allowed_error=1.0;
  allowed_error *= PI*1e9/180.;  // Convert from deg/ns to rad/s

  reltol = GetRealInitValue("reltol",1e-5);
  abstol = GetRealInitValue("abstol",1e-8);

  renormtol = GetRealInitValue("renormtol",1e-5);

  string stylestr = GetStringInitValue("ode_style");
  if(stylestr.compare("nonstiff")==0) {
    style = NONSTIFF;
  } else if(stylestr.compare("stiff")==0) {
    style = STIFF;
  } else {
    string msg = "Tuw_CvodeEvolve::Tuw_CvodeEvolve init error;"
      " Unrecognized ode_style request: " + stylestr;
    throw Oxs_Ext::Error(msg.c_str());
  }

  datablock.do_precess = GetIntInitValue("do_precess",1);

  datablock.alpha = GetRealInitValue("alpha"); // alpha is required input

  // User may specify either gamma_G (Gilbert) or
  // gamma_LL (Landau-Lifshitz) gyromagnetic ratio.
  // Code uses "gamma" which is LL form.
  if(HasInitValue("gamma_G") && HasInitValue("gamma_LL")) {
    throw Oxs_Ext::Error(this,"Invalid Specify block; "
			 "both gamma_G and gamma_LL specified.");
  } else if(HasInitValue("gamma_G")) {
    datablock.gamma =
      GetRealInitValue("gamma_G")/(1+datablock.alpha*datablock.alpha);
  } else if(HasInitValue("gamma_LL")) {
    datablock.gamma = GetRealInitValue("gamma_LL");
  } else {
    datablock.gamma = 2.211e5/(1+datablock.alpha*datablock.alpha);
  }
  datablock.gamma = fabs(datablock.gamma); // Force positive

  exact_report_period
    = GetUIntInitValue("exact_report_period",0);

  // Setup output objects
  max_dm_dt_output.Setup(this,InstanceName(),
     "Max dm/dt","deg/ns",0,
     &Tuw_CvodeEvolve::UpdateDerivedOutputs);
  cvode_reset_count_output.Setup(this,InstanceName(),
     "CVODE reset count","",0,
     &Tuw_CvodeEvolve::UpdateDerivedOutputs);
  Delta_E_output.Setup(this,InstanceName(),
     "Delta E","J",0,
     &Tuw_CvodeEvolve::UpdateDerivedOutputs);

  VerifyAllInitArgsUsed();
}

OC_BOOL Tuw_CvodeEvolve::Init()
{
  // Register output objects
  max_dm_dt_output.Register(director,-5);
  cvode_reset_count_output.Register(director,-5);
  Delta_E_output.Register(director,-5);

  Oxs_TimeEvolver::Init();  // Initialize parent class.
  /// Do this after child output registration so that
  /// UpdateDerivedOutputs gets called before the parent
  /// total_energy_output update function.

  // Instance variables
  state_id=0;   // Mark as invalid state
  stage_number=static_cast<OC_UINT4m>(static_cast<OC_INT4m>(-1));
  stepsize=initial_timestep;
  /// Assume this is an invalid stage number.
  if(cvode_initialized) {
    N_VFree(yout);          yout=NULL;
    CVodeFree(cvode_mem);   cvode_mem=NULL;
    cvode_initialized=0;
  }
  cvode_reset_request=0;
  cvode_reset_count=0;
  string total_energy_name = InstanceName();
  total_energy_name += string(":Total energy");
  datablock.state_energy_func =
    director->FindOutputObject(total_energy_name.c_str());
  if(datablock.state_energy_func==NULL) {
    string msg = string("Tuw_CvodeEvolve::Init: ")
      + total_energy_name
      + string(" output object not found.");
    throw Oxs_Ext::Error(msg.c_str());
  }
  datablock.state_energy_func->SetOutputFormat("%.17g");
  datablock.state_energy_func->CacheRequestIncrement(1);

  if(logfile!=NULL) {
    fclose(logfile);
    logfile=NULL;
    logfilename="";
  }
  if(director->GetLogFilename(logfilename)==0) {
    logfile = fopen(logfilename.c_str(),"a");
    if(logfile==NULL) {
      string msg =
	string("Tuw_CvodeEvolve::Init: Unable to open"
	       " log file \"")
	+ logfilename
	+ string("\" for appending.");
      throw Oxs_Ext::Error(msg.c_str());
    }
  }

  return 1;
}

Tuw_CvodeEvolve::~Tuw_CvodeEvolve()
{
  if(datablock.state_energy_func!=NULL) {
    datablock.state_energy_func->CacheRequestIncrement(-1);
  }
  if(cvode_initialized) {
    N_VFree(yout);
    CVodeFree(cvode_mem);
  }
  if(logfile!=NULL) {
    fclose(logfile);
    logfile=NULL;
    logfilename="";
  }
}

void Tuw_CvodeEvolve::ComputeEnergyAndMaxDmDt
(const Oxs_SimState& state,
 OC_REAL8m& max_dm_dt,
 OC_REAL8m& total_energy)
{
  Oxs_MeshValue<OC_REAL8m>& energy = datablock.energy;
  Oxs_MeshValue<ThreeVector>& mxH = datablock.mxH;
  const Oxs_Mesh* mesh = state.mesh;
  energy.AdjustSize(mesh);
  mxH.AdjustSize(mesh);
  OC_REAL8m pE_pt=0.;
  GetEnergyDensity(state,energy,&mxH,NULL,pE_pt);

  // Zero torque at fixed spin sites
  UpdateFixedSpinList(mesh);
  const OC_UINT4m fixed_count = GetFixedSpinCount();
  for(OC_UINT4m j=0;j<fixed_count;j++) {
    mxH[GetFixedSpin(j)].Set(0.,0.,0.); // What about pE_pt???
  }

  datablock.state_energy_func->Output(&state,
				       datablock.mif_interp,0,NULL);
  OC_BOOL err;
  total_energy
    = Nb_Atof(Tcl_GetStringResult(datablock.mif_interp),err);
  if(err) {
    string msg = string("Error in ComputeEnergyAndMaxDmDt ---"
			" Energy output return value: \"")
      + string(Tcl_GetStringResult(datablock.mif_interp))
      + string("\" does not represent a floating point value.");
    throw Oxs_Ext::Error(this,msg);
  }

  // Find max |mxH|
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);
  OC_UINT4m size = mesh->Size();
  OC_REAL8m max_mxH_sq = 0.0;
  for(OC_UINT4m i=0;i<size;i++) {
    if(Ms[i]==0) continue;
    OC_REAL8m mxH_sq = mxH[i].MagSq();
    if(mxH_sq>max_mxH_sq) max_mxH_sq = mxH_sq;
  }

  // Compute max |dm_dt| from max |mxH|
  max_dm_dt = fabs(datablock.gamma)*sqrt(max_mxH_sq);
  if(!datablock.do_precess) {
    max_dm_dt *= fabs(datablock.alpha);
  } else {
    max_dm_dt *= sqrt(1 + (datablock.alpha * datablock.alpha));
  }
}

void GetDmDt(integer vecsize, real totaltime,
	     void* y_dummy, void* ydot_dummy,
	     void* f_data)
{
  Tuw_CvodeEvolveDataBlock* datablock
    = static_cast<Tuw_CvodeEvolveDataBlock*>(f_data);
  N_Vector y    = static_cast<N_Vector>(y_dummy);
  N_Vector ydot = static_cast<N_Vector>(ydot_dummy);


  Tuw_CvodeEvolve* evolver = datablock->evolveptr;
  Oxs_MeshValue<OC_REAL8m>& energy = datablock->energy;
  Oxs_MeshValue<ThreeVector>& mxH = datablock->mxH;

  // Get reference to writable, non-locked state
  Oxs_SimState& work_state = datablock->state_key->GetWriteReference();

  const Oxs_Mesh* mesh = work_state.mesh;
  OC_UINT4m size = mesh->Size();
  if(3*size != OC_UINT4m(vecsize)) {
    char buf[1024];
    Oc_Snprintf
      (buf,sizeof(buf),
       "Tuw_CvodeEvolve friend GetDmDt: 3*mesh->Size()=%lu != vecsize=%lu",
       (unsigned long)(3*size),(unsigned long)vecsize);
    OXS_THROW(Oxs_ProgramLogicError,buf);
  }

  if(vecsize != N_VLENGTH(y) || vecsize != N_VLENGTH(ydot)) {
    OXS_THROW(Oxs_ProgramLogicError,
       "Tuw_CvodeEvolve friend GetDmDt: vecsize and N_VLENGTH mismatch");
  }

  // Fill spin vector from y
  Oxs_MeshValue<ThreeVector>& spin = work_state.spin;
  spin.AdjustSize(mesh);
  OC_UINT4m i,j;
  real* y_data = N_VDATA(y);
  for(i=j=0;i<size;i++,j+=3) {
    spin[i].x = y_data[j];
    spin[i].y = y_data[j+1];
    spin[i].z = y_data[j+2];
    // Normalize?  We might want to add a correction
    // term to dm_dt, to bring y_data back towards
    // unit length.
    spin[i].MakeUnit();
  }

  // Fill work state time index
  work_state.stage_elapsed_time
    = totaltime - work_state.stage_start_time;

  // Lock state down
  const Oxs_SimState& fixed_state
    = datablock->state_key->GetReadReference();

  // Get mxH
  energy.AdjustSize(mesh);
  mxH.AdjustSize(mesh);
  OC_REAL8m pE_pt=0.;
  evolver->GetEnergyDensity(fixed_state,energy,&mxH,NULL,pE_pt);

  // Zero torque at fixed spin sites
  evolver->UpdateFixedSpinList(mesh);
  const OC_UINT4m fixed_count = evolver->GetFixedSpinCount();
  for(OC_UINT4m k=0;k<fixed_count;k++) {
    mxH[evolver->GetFixedSpin(k)].Set(0.,0.,0.);
  }

  datablock->state_energy_func->Output(&fixed_state,
				       datablock->mif_interp,0,NULL);
  OC_BOOL err;
  datablock->total_energy
    = Nb_Atof(Tcl_GetStringResult(datablock->mif_interp),err);
  if(err) {
    string msg = string("Energy output return value: \"")
      + string(Tcl_GetStringResult(datablock->mif_interp))
      + string("\" does not represent a floating point value.");
    OXS_THROW(Oxs_ProgramLogicError,msg);
  }

  // Calculate dm_dt from mxH
  OC_REAL8m coef1 = -fabs(datablock->gamma);
  OC_REAL8m coef2 = -fabs(datablock->alpha);
  OC_BOOL do_precess = datablock->do_precess;
  real* ydot_data = N_VDATA(ydot);

  const Oxs_MeshValue<OC_REAL8m>& Ms = *(fixed_state.Ms);
  
  OC_REAL8m max_dm_dt_sq=-1;
  ThreeVector scratch,dm_dt;
  for(i=j=0;i<size;i++,j+=3) {
    if(Ms[i]==0) {
      ydot_data[j]=0.0;
      ydot_data[j+1]=0.0;
      ydot_data[j+2]=0.0;
    } else {
      scratch = mxH[i];
      scratch *= coef1;
      if(do_precess) {
	dm_dt = scratch;
      } else {
	dm_dt.Set(0.0,0.0,0.0);
      }
      scratch ^= spin[i]; // -|gamma|((mxH)xm)
      scratch *= coef2;  // |alpha.gamma|((mxH)xm) = -|alpha.gamma|(mx(mxH))
      dm_dt += scratch;
      OC_REAL8m dm_dt_sq = dm_dt.MagSq();
      if(dm_dt_sq > max_dm_dt_sq) max_dm_dt_sq = dm_dt_sq;
      ydot_data[j]=dm_dt.x;
      ydot_data[j+1]=dm_dt.y;
      ydot_data[j+2]=dm_dt.z;
    }
  }
  datablock->max_dm_dt = ( max_dm_dt_sq>=0.0 ? sqrt(max_dm_dt_sq) : -1.);
}

// C-style wrapper for GetDmDt
extern "C" {
  void GetDmDt_Cstyle
  (integer vecsize, real totaltime,
   N_Vector y, N_Vector ydot,
   void* f_data)
  {
    GetDmDt(vecsize,totaltime,y,ydot,f_data);
  }
}

OC_BOOL
Tuw_CvodeEvolve::Step(const Oxs_TimeDriver* driver,
                Oxs_ConstKey<Oxs_SimState> current_state,
                const Oxs_DriverStepInfo& /* step_info */,
                Oxs_Key<Oxs_SimState>& next_state)
{
  const Oxs_SimState& cstate = current_state.GetReadReference();
  Oxs_SimState* workstate = &(next_state.GetWriteReference());
  driver->FillState(cstate,*workstate);

  if(cstate.mesh->Id() != workstate->mesh->Id()) {
    throw Oxs_Ext::Error
      (this,"Tuw_CvodeEvolve::Step: Oxs_Mesh not fixed across steps.");
  }

  if(cstate.Id() != workstate->previous_state_id) {
    throw Oxs_Ext::Error
      (this,"Tuw_CvodeEvolve::Step: State continuity break detected.");
  }

  if(cstate.stage_number != workstate->stage_number) {
    throw Oxs_Ext::Error
      (this,"Tuw_CvodeEvolve::Step: Can't \"Step\" across stages.");
  }

  if(stage_number != workstate->stage_number) {
    // Re-initialize solver across stage events
    cvode_reset_request = 1;
  }

  // Save energy from current state, so Delta E can be computed
  // below.
  datablock.state_energy_func->Output(&cstate,
				      datablock.mif_interp,0,NULL);
  OC_BOOL err;
  OC_REAL8m old_energy
    = Nb_Atof(Tcl_GetStringResult(datablock.mif_interp),err);
  if(err) {
    string msg = string("Energy output return value: \"")
      + string(Tcl_GetStringResult(datablock.mif_interp))
      + string("\" does not represent a floating point value.");
    throw Oxs_Ext::Error(this,msg);
  }

  OC_UINT4m size = cstate.mesh->Size();

  // Adjust stepsize to be within requested bounds
  if(stepsize<min_timestep) stepsize = min_timestep;
  if(stepsize>max_timestep) stepsize = max_timestep;

  // Negotiate with driver over size of next step
  workstate->last_timestep=stepsize;
  workstate->stage_start_time = cstate.stage_start_time;
  workstate->stage_elapsed_time = cstate.stage_elapsed_time + stepsize;
  /// Note: cstate and *workstate have same stage_number.
  workstate->iteration_count = cstate.iteration_count + 1;
  workstate->stage_iteration_count = cstate.stage_iteration_count + 1;
  driver->FillStateSupplemental(*workstate);
  int itask = ONE_STEP;
  if(stepsize>workstate->last_timestep) {
    // Driver cut requested stepsize.  Force CVode call to
    // integrate out to exactly stepsize.  This won't catch
    // all overstep problems, but should get most of them.
    itask = NORMAL;
  }
  stepsize = workstate->last_timestep;
  workstate->spin.AdjustSize(workstate->mesh); // Safety

  // Protect against 0.0 or effectively 0.0 stepsize.  The
  // extra 2 here is just a safety fudge factor.
  if(2*cstate.stage_elapsed_time + stepsize
     <= 2*cstate.stage_elapsed_time) {
    stepsize = 2*cstate.stage_elapsed_time*OC_REAL8_EPSILON;
  }

  real start_time = cstate.stage_start_time+cstate.stage_elapsed_time;
  // Perhaps time 0 should be stage start instead of run start?

  // Update user datablock for Cvode call
  datablock.state_key = &next_state;
  workstate = NULL; // GetDmDt locks down next_state, so we shouldn't
  /// keep a non-const pointer to it.

  // Do step
  if(!cvode_initialized || cvode_reset_request) {
    if(!cvode_initialized) {
      // Allocate memory for cvode internal spin array
      yout = N_VNew(3*size,NULL);
      if(yout==NULL) {
	throw Oxs_Ext::Error
	  ("Tuw_CvodeEvolve::Step: Insufficient memory to create yout.");
      }
    } else {
      // This is a reset request.
      CVodeFree(cvode_mem);
    }

    // Fill yout with initial spin configuration
    unsigned int i,j;
    real* yout_data = N_VDATA(yout);
    for(i=j=0;i<size;++i,j+=3) {
      const ThreeVector& tempspin = cstate.spin[i];
      yout_data[j]   = tempspin.x;
      yout_data[j+1] = tempspin.y;
      yout_data[j+2] = tempspin.z;
    }

    // Allocate and initialize CVode memory
    switch(style) {
    case NONSTIFF:
      cvode_mem = CVodeMalloc(3*size,GetDmDt_Cstyle,
			      start_time,yout,ADAMS,
			      FUNCTIONAL,SS,&reltol,&abstol,
			      &datablock,logfile,FALSE,NULL,NULL,NULL);
      break;
    case STIFF:
      cvode_mem = CVodeMalloc(3*size,GetDmDt_Cstyle,
			      start_time,yout,BDF,
			      NEWTON,SS,&reltol,&abstol,
			      &datablock,logfile,FALSE,NULL,NULL,NULL);
      break;
    default:
      throw Oxs_Ext::Error("Tuw_CvodeEvolve::Step: Invalid ODESolverStyle");
    }
    CVSpgmr(cvode_mem,NONE,MODIFIED_GS,0,0.,NULL,NULL,NULL);

    // CVode doesn't like very tiny stepsizes when starting an integration.
    if(start_time + stepsize <= start_time*(1+1024*OC_REAL8_EPSILON) ) {
      stepsize += 1024*start_time*OC_REAL8_EPSILON;
    }

    cvode_initialized = 1;
    cvode_reset_request=0;
    ++cvode_reset_count;
  }

  real end_time;
  int success
    = CVode(cvode_mem,start_time+stepsize,yout,&end_time,itask);
  if(success != SUCCESS) {
    char buf[1024];
    const char* errtype=NULL;
    switch(success) {
    case 0: // If test above should prevent this branch
      errtype = "SUCCESS";
      break;
    case -1:
      errtype = "CVODE_NO_MEM";
      break;
    case -2:
      errtype = "ILL_INPUT";
      break;
    case -3:
      errtype = "TOO_MUCH_WORK";
      break;
    case -4:
      errtype = "TOO_MUCH_ACC";
      break;
    case -5:
      errtype = "ERR_FAILURE";
      break;
    case -6:
      errtype = "CONV_FAILURE";
      break;
    case -7:
      errtype = "SETUP_FAILURE";
      break;
    case -8:
      errtype = "SOLVE_FAILURE";
      break;
    default:
      errtype = NULL; // Unknown error.
      break;
    }
    if(errtype!=NULL) {
      Oc_Snprintf(buf,sizeof(buf),
		  "Tuw_CvodeEvolve::Step: CVode error=%s",errtype);
    } else {
      Oc_Snprintf(buf,sizeof(buf),
		  "Tuw_CvodeEvolve::Step: Unrecognized CVode ERROR,"
		  " errorcode=%d",success);
    }
    if(logfile!=NULL) {
      strncat(buf,"\nFor additional details, see message preceding"
	      " this one in the Oxs log file \"",
	      sizeof(buf)-strlen(buf)-1);
      strncat(buf,logfilename.c_str(),sizeof(buf)-strlen(buf)-1);
      strncat(buf,"\".",sizeof(buf)-strlen(buf)-1);
    }
    throw Oxs_Ext::Error(buf);
  }

  // At this point, GetDmDt has locked down next_state, but the spin
  // configuration in next_state is not the same as CVode reports as the
  // final y configuration.  In fact, typically CVode does not evaluate
  // f(y) (i.e., GetDmDt(y)) at the final "corrected" y value for a
  // step.
  //   In order to save an extra and generally unneeded (and expensive)
  // DmDt evaluation at each step, the user may specify
  // exact_report_period in the Tuw_CvodeEvolve Specify block to
  // control this extra GetDmDt evaluation.  If exact_report_period
  // == 0, then this evaluation happens only at potential StageDone
  // events.  If exact_report_period>0, then it happens at each step
  // for which iteration_count%exact_report_period == 0 in addition
  // to potential StageDone events.
  //  If the extra GetDmDt evaluation is required, then the code below
  // opens a new state, thereby invalidating any cached data
  // for the old (most recently evaluated by GetDmDt) state, and
  // updates the spin and related data.
  //  Otherwise, the code below CHEATS, by casting away the const on
  // the state pointer obtained from the locked down new_state, and
  // updates the spin and related data.  Any cached data, in
  // particular max_dm_dt and energy values will not be marked invalid
  // because the state_id doesn't change.  This breaks the design
  // rules for state manipulation, so this is potentially dangerous,
  // but it appears to be the best option at present.  If problems
  // arise, the user may always override this behavior by setting
  // exact_report_period to 1.

  const Oxs_SimState* fixedstate = &(next_state.GetReadReference());
  /// Use GetReadReference() instead of GetPtr() as a safety, to make
  /// sure next_state is locked down.

  OC_BOOL state_cheat = 0;
  if(exact_report_period<1 ||
     fixedstate->iteration_count % exact_report_period != 0) {
    // Perform state cheat described above.
    workstate = const_cast<Oxs_SimState*>(fixedstate);
    state_cheat = 1;
  } else {
    // Otherwise, open new state.
    workstate = &(next_state.GetWriteReference());
  }

  // Update workstate timing information
  // NB: cstate and *workstate have same stage index.
  stepsize = end_time - start_time;
  workstate->last_timestep = stepsize;
  workstate->stage_elapsed_time = cstate.stage_elapsed_time + stepsize;
  /// If start_time is later changed to index off of stage start
  /// instead of run start, then the RHS above is just end_time.

  // Write new spin configuration into next_state
  ThreeVector tempspin;
  unsigned int i,j;
  real* yout_data = N_VDATA(yout);
  OC_REAL8m maxdrift = 0;
  for(i=j=0;i<size;++i,j+=3) {
    tempspin.x = yout_data[j];
    tempspin.y = yout_data[j+1];
    tempspin.z = yout_data[j+2];

    // Check size
    OC_REAL8m mag = sqrt(tempspin.MagSq());
    if(mag>0.03125) tempspin *= 1.0/mag;
    else            tempspin.SetMag(1.0);
    OC_REAL8m drift = fabs(1.0-mag);
    if(drift>maxdrift) maxdrift=drift;
    workstate->spin[i] = tempspin;
  }
  // If solver spins have drifted farther than renormtol from
  // unity, then on next pass reinitialize solver with unit spins.
  if(maxdrift>renormtol) cvode_reset_request=1;

  // Derived data
  if(state_cheat) {
    fixedstate = workstate;
    workstate = NULL;
    fixedstate->AddDerivedData("Max dm/dt",
			      datablock.max_dm_dt);
    fixedstate->AddDerivedData("Delta E",
			      datablock.total_energy - old_energy);
    fixedstate->AddDerivedData("CVODE reset count",cvode_reset_count);
    // Check for stage done event
    if(driver->IsStageDone(*fixedstate)) {
      // Stage done.  Uncheat.
      state_cheat = 0;
      workstate = &(next_state.GetWriteReference()); // Toggle state id
      workstate->ClearDerivedData();
      fixedstate=NULL; // Safety
    }
  }
  workstate = NULL;
  if(!state_cheat) {
    fixedstate = &(next_state.GetReadReference());
    OC_REAL8m fixed_max_dm_dt,fixed_energy;
    ComputeEnergyAndMaxDmDt(*fixedstate,fixed_max_dm_dt,fixed_energy);
    fixedstate->AddDerivedData("Max dm/dt",fixed_max_dm_dt);
    fixedstate->AddDerivedData("Delta E",fixed_energy - old_energy);
    fixedstate->AddDerivedData("CVODE reset count",cvode_reset_count);
    /// If state_cheat was previously true, then the preceding
    /// AddDerivedData calls will overwrite the old approximate
    /// values.
  }

  state_id = fixedstate->Id();
  stage_number = fixedstate->stage_number;

  return 1;  // Good step
}

void Tuw_CvodeEvolve::UpdateDerivedOutputs(const Oxs_SimState& state)
{ // This routine fills all the Tuw_CvodeEvolve Oxs_ScalarOutput's to
  // the appropriate value based on the import "state".  Also any of
  // Oxs_VectorOutput's that have CacheRequest enabled are filled. 
  max_dm_dt_output.cache.state_id
    = cvode_reset_count_output.cache.state_id
    = Delta_E_output.cache.state_id
    = 0;

  if(!state.GetDerivedData("Max dm/dt",
			   max_dm_dt_output.cache.value) ||
     !state.GetDerivedData("Delta E",
			   Delta_E_output.cache.value)) {
    // Missing at least one value.  Recompute
    OC_REAL8m max_dm_dt,energy;
    ComputeEnergyAndMaxDmDt(state,max_dm_dt,energy);
    if(!state.GetDerivedData("Max dm/dt",
			     max_dm_dt_output.cache.value)) {
      max_dm_dt_output.cache.value = max_dm_dt;
      state.AddDerivedData("Max dm/dt",max_dm_dt);
    }
    if(!state.GetDerivedData("Delta E",
			   Delta_E_output.cache.value)) {
      // We don't have last energy data, so punt---assume 0.
      // The proper way to fix this is probably to store
      // energy and old_energy in each state, the former
      // (and possibly the latter) as DerivedData.  OTOH,
      // maybe eventually *all* outputs get their caches
      // moved into the state.
      Delta_E_output.cache.value = energy;
      state.AddDerivedData("Delta E",Delta_E_output.cache.value);
    }
  }
  max_dm_dt_output.cache.value *= (180e-9/PI);
  /// Convert from radians/second to deg/ns

  if(!state.GetDerivedData("CVODE reset count",
			   cvode_reset_count_output.cache.value)) {
    // cvode_reset_count not stored for this state.  Punt: use
    // current value.
    cvode_reset_count_output.cache.value = 
      static_cast<OC_REAL8m>(cvode_reset_count);
  }

  max_dm_dt_output.cache.state_id
    = cvode_reset_count_output.cache.state_id
    = Delta_E_output.cache.state_id
    = state.Id();
}
