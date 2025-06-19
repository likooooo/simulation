#include <iostream>
#include "RS.h"

int main()
{
	// ------------------init new simulation----------------------------
	std::cout << "Here begin"  << std::endl;
	
	// set the lattice vectors (1,0) and (0,0) -- 1D lattice
	RS_real Lr[4] = { 1, 0, 0, 1 };
	unsigned int nG = 1;
	unsigned int argflag = 0;
	RS_Simulation PS0;
	RS_Simulation *ps1 = &PS0;
	RS_Simulation** pS ;
	pS= &ps1;
	*pS = RS_Simulation_New(Lr, nG, NULL);
	RS_Simulation *S = *pS;
	S->options.use_less_memory = true;
	RS_Simulation_SetLattice(S, Lr);
	std::cout << "Here SetLattice: Lr is    " << Lr[0] << "," << Lr[1] << "," << Lr[2]<< "," << Lr[3]<<  std::endl;

    // -------------------Set NumG ------------------------------------------------
	Simulation_SetNumG(S, 27);
	std::cout << "Here SetNumG"  << std::endl;

    //-------------------- Material definition--------------------------------------
	//AddMaterial("Silicon", {12,0}) -- real and imag parts
    const char *name1 = "Silicon";
	int type1 = RS_MATERIAL_TYPE_SCALAR_COMPLEX; // epsilon type
	RS_real eps1[2] =  {12, 0};
	RS_MaterialID M1 = RS_Simulation_SetMaterial(S, -1, name1, type1, eps1);
	std::cout << "Here AddMaterial: Silicon"  << std::endl;

    //AddMaterial("Vacuum", {1,0})
	const char *name2 = "Vacuum";
	int type2 = RS_MATERIAL_TYPE_SCALAR_COMPLEX; // epsilon type
	RS_real eps2[2] =  {1,0};
	RS_MaterialID M2 = RS_Simulation_SetMaterial(S, -1, name2, type2, eps2);
	std::cout << "Here AddMaterial: Vacuum"  << std::endl;

    //------------------------- Layer Add----------------------------------------
	// AddLayer(
	//  'AirAbove', --name
	//  0,          --thickness
	//  'Vacuum')   --background material 
	RS_LayerID layer1;
	const char *layername1 = "AirAbove";
	RS_real thickness1 = 0;
	const char *matname1 = "Vacuum";
	RS_MaterialID LM1 = RS_Simulation_GetMaterialByName(S, matname1);
	layer1 = RS_Simulation_SetLayer(S, -1, layername1, &thickness1, -1, LM1);
	std::cout << "Here AddLayer 1"  << std::endl;

	// AddLayer('Slab', 0.5, 'Vacuum')
	RS_LayerID layer2;
	const char *layername2 = "Slab";
	RS_real thickness2 = 0.5;
	const char *matname2 = "Vacuum";
	RS_MaterialID LM2 = RS_Simulation_GetMaterialByName(S, matname2);
	layer2 = RS_Simulation_SetLayer(S, -1, layername2, &thickness2, -1, LM2);
	std::cout << "Here AddLayer 2"  << std::endl;

    // Set Layer Pattern as Rectangle('Slab',        -- which layer to alter
    //                        'Silicon',     -- material in rectangle
	//                        {0,0},         -- center
	//                        0,             -- tilt angle (degrees)
	//                        {0.25, 0.5}) -- half-widths
    const char *layer_name2 = "Slab";
	const char *material_name2 = "Silicon";
	RS_LayerID patternlayer2 = RS_Simulation_GetLayerByName(S, layer_name2);
	double center2[2] = {0,0}, halfwidths2[2] = {0.25, 0.5};
    RS_MaterialID PM2 = RS_Simulation_GetMaterialByName(S, material_name2);
	RS_real angle2 = 0 / 360.;
	int ret2;
    if(0 == Lr[1] && 0 == Lr[2] && 0 == Lr[3]){
		ret2 = RS_Layer_SetRegionHalfwidths(S, patternlayer2, PM2, RS_REGION_TYPE_INTERVAL, halfwidths2, center2, &angle2);
		std::cout << "Here SetLayerPatternRectangle 2 : enter if"  << std::endl;
	}else{
		ret2 = RS_Layer_SetRegionHalfwidths(S, patternlayer2, PM2, RS_REGION_TYPE_RECTANGLE, halfwidths2, center2, &angle2);
		std::cout << "Here SetLayerPatternRectangle 2 : enter else"  << std::endl;
	}

	// AddLayerCopy('AirBelow',         --new layer name
	//              0,                  --thickness
	//             'AirAbove')          --layer to copy
	RS_LayerID layer5;
	const char *layername5 = "AirBelow";
	RS_real thickness5 = 0;
	const char *copyname5 = "AirAbove";
	RS_LayerID Lcopy5 = RS_Simulation_GetLayerByName(S, copyname5);
	layer5 = RS_Simulation_SetLayer(S, -1, layername5, &thickness5, Lcopy5, -1);
	std::cout << "Here AddLayer 5"  << std::endl;

	//-------------------------------------Set excitation------------------------------------------- 
	// -- E polarized along the grating "rods"
	// SetExcitationPlanewave(
	// 	{0,0},  -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
	// 	{0,0},  -- s-polarization amplitude and phase (in degrees);TE
	// 	{1,0})  -- p-polarization amplitude and phase;TM
	double angleE[2] = {M_PI/180.* 0, M_PI/180.* 0};
	double pol_s[2] = {0, M_PI/180.* 0}; /* s polarization; E out of plane */
	double pol_p[2] = {1, M_PI/180.* 0}; /* p polarization; E in plane */
	int order = 0;
	int retE = Simulation_MakeExcitationPlanewave(S, angleE, pol_s, pol_p, order);
	std::cout << "Here MakeExcitationPlanewave"  << std::endl;

	//-------------------------------------Set Frequency--------------------------------------------- 
	// SetFrequency(0.64)
	RS_real freq[2] = {0.01, 0}; // real part and img part 
	RS_Simulation_SetFrequency(S, freq);
	std::cout << "Here SetFrequency"  << std::endl;

	//-------------------------------------Calculate--------------------------------------------------
	//for x=-0.5,3.5,0.02 do
	//  for z=-1,1.5,0.02 do
	//	   Ex,Ey,Ez = S:GetEField({x,0,z})
	//	   print(x..'\t'..z..'\t'..Ex..'\t'..Ey..'\t'..Ez)
	//  end
	//  print('')
    //end
	int ret6;
	double r[3], fE[6];
	for (double x=-0.5; x < 3.5; x += 0.02){
		double z = -1;
		// for (double z=-1; z < 1.5; z += 0.02)
		{
			r[0] = x;
			r[1] = 0;
			r[2] = z;
			ret6 = Simulation_GetField(S, r, fE, NULL);
			std::cout << "x:  " << x << ",   z:  " << z << ",    Ex: " << fE[0]<< ",    Ey: " << fE[1]<< ",    Ez: " << fE[2] << std::endl;

		}
	}
	
}