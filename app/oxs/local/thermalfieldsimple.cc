/* FILE: thermalfieldsimple.cc
*
* Example anisotropy class implementation.
* This class is derived from the Oxs_Energy class.
*
*/
#include "thermalfieldsimple.h"

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(Oxs_ThermalFieldSimple);
/* End includes */

//#define MU0 12.56637061435917295385e-7
/* 4 PI 10^7 */
// Constructor
Oxs_ThermalFieldSimple::Oxs_ThermalFieldSimple(
const char* name,
// Child instance id
Oxs_Director* newdtr, // App director
//Tcl_Interp* safe_interp, // Safe interpreter
const char* argstr)
// MIF input block parameters
: Oxs_Energy(name,newdtr,argstr),mesh_id(0)
{
  // Process arguments
  if(HasInitValue("T")){
    OXS_GET_INIT_EXT_OBJECT("T",Oxs_ScalarField,T_init);
  } else {
      // Set default temperature scalar field to 300 K
      //Oxs_Ext* foo=MakeNew("Oxs_UniformScalarField",director,"value 300.0");
      //T_init.SetAsOwner(dynamic_cast<Oxs_ScalarField*>(foo));
      throw Oxs_ExtError(this,"Error Oxs_ThermalFieldSimple: No temperature T defined.\n");
  }
 if(HasInitValue("alpha")){
  OXS_GET_INIT_EXT_OBJECT("alpha",Oxs_ScalarField,alpha_init);
  } else {
      // Set default alpha to 0.01
       //Oxs_Ext* foo=MakeNew("Oxs_UniformScalarField",director,"value 0.01");
       //alpha_init.SetAsOwner(dynamic_cast<Oxs_ScalarField*>(foo));
      throw Oxs_ExtError(this,"Error Oxs_ThermalFieldSimple: No damping coefficient alpha defined.\n");
  } 
//allows boxwise definition
 //if(HasInitValue("D")){
 // OXS_GET_INIT_EXT_OBJECT("D",Oxs_ScalarField,D_init);
 // } else {
 //     // Set default D to 1.0
 //      Oxs_Ext* foo=MakeNew("Oxs_UniformScalarField",director,"value 1.");
 //      D_init.SetAsOwner(dynamic_cast<Oxs_ScalarField*>(foo));
 //     //throw Oxs_ExtError(this,"Error Oxs_ThermalField: No thermal fluctuations strength D defined.\n");
 // } 

 // D = GetRealInitValue("D",1.0);
 // min_time_step = GetRealInitValue("min_time_step",1.0e-14);
  time_step = GetRealInitValue("time_step",1.0e-14);
  VerifyAllInitArgsUsed();
 // start the clock used for seeding
 beginning = myclock::now();
}

OC_BOOL Oxs_ThermalFieldSimple::Init()
{
mesh_id = 0;
T.Release();
alpha.Release();
beginning = myclock::now();
return Oxs_Energy::Init(); 
}

void Oxs_ThermalFieldSimple::GetEnergy
(const Oxs_SimState& state,
Oxs_EnergyData& oed
) const
{


  //beginning = myclock::now();
  const Oxs_Mesh* mesh = state.mesh;
  const Oxs_MeshValue<OC_REAL8m>& Ms_inverse = *(state.Ms_inverse);
  const Oxs_MeshValue<OC_REAL8m>& Ms = *(state.Ms);
  //	OC_REAL8m          dt =  state.last_timestep;
  const Oxs_MeshValue<ThreeVector>& spin = state.spin;
  const OC_REAL8m KB     = 1.3806488e-23;  // J/K
  const OC_REAL8m GAMMAE = 1.760859708e11; // 1/(s*T)
        OC_UINT4m size   =  mesh->Size();
  const OC_REAL8m sign   = -1.;
  // Use supplied buffer space, and reflect that use in oed.
  oed.energy = oed.energy_buffer;
  oed.field = oed.field_buffer;
  Oxs_MeshValue<OC_REAL8m>& energy = *oed.energy_buffer;
  Oxs_MeshValue<ThreeVector>& field = *oed.field_buffer;

 // random number generators and normal distributions for x, y, and z direction
  std::default_random_engine generatorx;
  std::default_random_engine generatory;
  std::default_random_engine generatorz;
  std::normal_distribution<OC_REAL8m> distributionx (0.,1.);
  std::normal_distribution<OC_REAL8m> distributiony (0.,1.);
  std::normal_distribution<OC_REAL8m> distributionz (0.,1.);

  //if(dt<min_time_step) { dt = min_time_step;}
  OC_REAL8m dt = time_step;
  alpha_init->FillMeshValue(mesh,alpha);
  T_init->FillMeshValue(mesh,T);
  //D_init->FillMeshValue(mesh,D);

    //OC_CHAR msg[100];
    //Oc_Snprintf(msg,100,"KB: %g   GAMMAE: %g   MU0 : %g   dt: %g ",KB, GAMMAE, MU0, dt);
    //Oc_Snprintf(msg,100,"last_time_step: %g", state.last_timestep);
    //throw Oxs_ExtError(this,msg);
    //Oc_Snprintf(stdout,"last_time_step: %g", state.last_timestep);

if(mesh_id != state.mesh->Id()){
 // obtain a seed from the timer and initilaize random number setup
  myclock::duration dx = myclock::now() - beginning;
  unsigned seedx = dx.count();
  myclock::duration dy = myclock::now() - beginning;
  unsigned seedy = dy.count();
  myclock::duration dz = myclock::now() - beginning;
  unsigned seedz = dz.count();
  generatorx.seed(seedx);
  generatory.seed(seedy);
  generatorz.seed(seedz);
  distributionx.reset();
  distributiony.reset();
  distributionz.reset();
  mesh_id=0;
// These values do not change wihtin the subsequent loop so it is more efficient to calculate only once
  OC_REAL8m pre_factor  = 2.0;
  pre_factor           *= KB;
  pre_factor           /= MU0;
  pre_factor           /= GAMMAE;
  pre_factor           /= dt;
  pre_factor            = sqrt(pre_factor);

 for(OC_INDEX i=0 ;i<size;++i) {


    	OC_REAL8m cell_alpha  = alpha[i];
    	OC_REAL8m cell_T      = T[i]; 
    	OC_REAL8m cell_Volume = mesh->Volume(i); 
    	ThreeVector         m = spin[i]; 
    	ThreeVector       mxH = m;

    	OC_REAL8m cell_factor    = cell_alpha;
        cell_factor   *= cell_T;
       	cell_factor   /= cell_Volume;
       	cell_factor   *= Ms_inverse[i];
       	cell_factor   /= (1.+cell_alpha*cell_alpha);
       	cell_factor    = sqrt(cell_factor);
       	cell_factor   *= pre_factor;

    	//if(cell_factor==0.0) {
	energy[i] = 0.0;
	field[i].Set(0.0,0.0,0.0);
      	//continue;
    	//}else{
      
    	ThreeVector xi = ThreeVector(distributionx(generatorx),distributiony(generatory),distributionz(generatorz)); 

    	ThreeVector     H  =  xi*cell_factor;
        //              mxH ^=  H;
    	OC_REAL8m   mdotH  = (m*H); 
	//OC_REAL8m     ei   = H.MagSq();
    	OC_REAL8m   ei     = mdotH;
                    ei    *= Ms[i];
                    ei    *= MU0;
                    ei    *= sign;
	energy[i] = ei;
        field[i]  = H;
   //}
  }
   mesh_id   = state.mesh->Id();
 }
}

