#ifndef CONSTANTS_H_GUARD
#define CONSTANTS_H_GUARD

#include <unordered_map>
#include "TMath.h"
#include "TLorentzVector.h"

//static int num_mixed_p_pip = 0; 

static const int MAX_PARTS = 100; 
static const int NUM_THREADS = 1;

static const float c_special = 29.9792458; //speed of light in cm/ns
static const float c_convert = 10000000; //Convert c_special to m/s

//Beam energies in GeV
static const float energy_e16 = 5.7696;
//3_14 was 5.754, but led to bad results: updated to 5.759 after reading Paremuzyan's thesis 
//Check 4.794 GeV for e16
static const float energy_e1f = 5.499;

static const float fine_structure = 1.0/137.0; 

//Masses of relevant particles
static const	float me = 0.0005109989; //mass of electron in GeV
static const	float mp = 0.93828;	//Mass of proton in GeV
static const	float mpi = 0.1395;	//Mass of pion in GeV

//Range of W/Q2
static const float WminAna = 1.4;
static const float WmaxAna = 2.125;
static const float Q2minAna = 3.0;//2.0; Moved up to 3.0 in order to compare with simulation in that range
static const float Q2maxAna = 5.0; 

//Initial Four Vectors

const TLorentzVector k_mu_e16(0.0,0.0,sqrt(energy_e16*energy_e16-me*me),energy_e16);
const TLorentzVector k_mu_e1f(0.0,0.0,sqrt(energy_e16*energy_e1f-me*me),energy_e1f);
const TLorentzVector p_mu(0.0,0.0,0.0,mp);

//e1-6 Luminosity Values
static const float lt_e16 = 5.0; //Target length in cm
static const float Dt_e16 = 0.073; //Density of target in g/cm^3
static const float NA = 6.022 * pow(10.0,23); //Avogadro's number
static const float qe = 1.602 * pow(10.0,-19); // fundamental Coulomb charge 
static const float Mt_e16 = 1.007; //Molar mass of target in g/mole
static const float Qt_e16 = 21.32*pow(10.0,-3); //Total charge incident on target from Arjun in Coulombs //Will need to verify myself
static const float L_e16 = Qt_e16*lt_e16*Dt_e16*NA/(qe*Mt_e16);
static const float cm2_to_mbarn = pow(10.0,-27);

//Particle ID Numbers
//anti particles are negative
static const 	int PROTON = 2212;
static const	int ELECTRON = 11;
static const	int PION = 211;
static const	int PION_0 = 111;

// Beam spot from Arjun's Fortran vertex_e1f.f
static const float x_beam = 0.090;
static const float y_beam = -0.345;

//Fiducial Cut Parameters (electrons have e, hadrons have h)
//Arjun's cut_fid_e1f.f
static const float c1e = 12.0;
static const float c2e = 18.5;
static const float c3e = 0.25;
static const float c4e = 15.0;
static const float factor_e = 0.416667;
static const float p_shift_e = 0.14;
static const float a0xh[6] = {24.0,24.0,23.0,23.5,24.5,24.5};
static const float a0mh[6] = {25.0,26.0,26.0,25.5,27.0,26.0};
static const float a1xh[6] = {0.22,0.23,0.20,0.20,0.22,0.22};
static const float a1mh[6] = {0.22,0.22,0.22,0.22,0.16,0.16};
static const float a2xh[6] = {8.0,8.0,8.0,8.0,8.0,8.0};
static const float a2mh[6] = {8.0,8.0,8.0,8.0,8.0,8.0};
static const float a3xh[6] = {1.0,1.0,1.0,1.0,1.0,1.0};
static const float a3mh[6] = {1.0,1.0,1.0,1.0,1.0,1.0};

//Missing Mass Cuts
static const float pim_center = 0.164369;//These are all determined through fitting of b_wig fitting.h
static const float pim_sig = 0.0862474;
static const float pip_center = 0.157301;
static const float pip_sig = 0.0752388;
static const float p_center = 0.946847;
static const float p_sig = 0.0420984;
static const float pim_center2 = 0.022;
static const float pim_sig2 = 0.022;
static const float pip_center2 = 0.022;
static const float pip_sig2 = 0.022;
static const float p_center2 = 0.89;
static const float p_sig2 = 0.05;

static const float pim_bot_MM[11] = {0.0781216,0.0781216,0.0781216,0.0781216,0.0781216,0.0781216,0.0781216,0.0781216,0.083,0.109,0.109};
static const float pim_top_MM[11] = {0.2506164,0.2506164,0.2506164,0.2506164,0.2506164,0.2506164,0.2506164,0.2506164,0.34,0.415,0.415};


//For Project
static const float MM_n_center = 0.944;
static const float MM_n_sigma = 0.07;
static const float MM_piz_center = 0.140;
static const float MM_piz_sigma = 0.1;
static const float MM_D_center = 1.232;
static const float MM_D_sigma = 0.08;

//My Own MM cut parameters
static const float MM_zero_center = 0.0;
static const float MM_zero_sigma = 0.02;
static const float MM_zero_center2 = 0.0;
static const float MM_zero_sigma2 = 0.004;

//Arjun's cut parameters: 3rd order polynomial
static const float DTL[4] = {-0.778903, 0.027350, 0.047947, -0.009641};
static const float DTH[4] = {0.758057, -0.147383, 0.034343, -0.002367};

//attempt to remove slice of electrons from Pi- delta t
static const float dt_e_sig = 0.05;
static const float dt_e_A = 20.0;
static const float dt_e_a = 8.55;
static const float dt_e_b = 0.31;

static const int MinCC_Cut[6][18] = {{35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35},{35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35},{35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35},{35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35},{35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35},{35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35}};



//Kinematic Cuts
// Margin of error allowed for kinematic fitting 
const float px_dev = 0.02;
const float py_dev = 0.02;
const float pz_dev = 0.02;

//Electron EC cut parameters
//currently from Arjun and not verified
//High and Low functions are third order polynomials separtated by sector
const Float_t p_min_e16 = 0.70; //in GeV
const Float_t ec_min_e16 = 0.06; //in GeV 
const float sf_high_e16[6][4] = {
		{0.380401, -0.019463, 0.004609, -0.000359},
		{0.428533, -0.047554, 0.016350, -0.001938},
		{0.420563, -0.064622, 0.025420, -0.003105},
		{0.411866, -0.062074, 0.022444, -0.002734},
		{0.383041, -0.020678, 0.007576, -0.000897},
		{0.394516, -0.023219, 0.010402, -0.001431}
		};
const float sf_low_e16[6][4] = {
		{0.141603, 0.063643, -0.012677, 0.000918},
		{0.149598, 0.053507, -0.007504, 0.000323},
		{0.172657, 0.038093, -0.002066, -0.000355},
		{0.149623, 0.060981, -0.012957, 0.001054},
		{0.099589, 0.102036, -0.028452, 0.002751},
		{0.117104, 0.097765, -0.023769, 0.002149}
		};
const Float_t p_min_e1f = 0.64;
const Float_t ec_min_e1f_exp[6] = {0.058,0.064,0.060,0.056,0.058,0.056}; 
const Float_t ec_min_e1f_sim[6] = {0.063,0.063,0.063,0.063,0.063,0.063}; 

/*
//Binning
//Plot Formation Constants
//W Q2
const int WQxres = 300;
const int WQyres = 300;
const float WQxmin = -0.01;
const float WQymin = -0.01;
const float WQxmax = 3.99;
const float WQymax = 8.99;
//For Topology Selection
const int WQ2xres = 20;
const int WQ2yres = 20;

//Electron Sampling Fraction
const int SFxres = 100;
const int SFyres = 100;
const float SFxmin = 0.0;
const float SFymin = 0.0;
const float SFxmax = 6.0;
const float SFymax = 1.0;
//Minimum Electron Energy
const int MEExres = 100;
const int MEEyres = 100;
const float MEExmin = 0.0;
const float MEEymin = 0.0;
const float MEExmax = 6.0;
const float MEEymax = 1.0;
//Fiducial Cuts
const int FIDxres = 200;
const int FIDyres = 200;
const float FIDxmin = -30.0;
const float FIDymin = 0.0;
const float FIDxmax = 30.0;
const float FIDymax = 180.0;
//Delta_t
const int DTxres = 600;//300;
const int DTyres = 400;
const float DTxmin = 0.0;
const float DTymin = -5.0;//-4.0;
const float DTxmax = 4.5;//7.0;
const float DTymax = 5.0;//4.0;*/

//Missing Mass
static const int MMxres = 1500;//[2][4] = {
	//{1500,500,500,500},
	//{1500,500,500,500}
//};
static const float MMxmin2[4] = {-0.2,-0.2,-0.2,-0.2};
/*[2][4] = {
	{-0.2,-0.2,-0.2,-0.2},
	{-0.2,-0.2,-0.2,-0.2}
};*/
static const float MMxmax2[4] = {1.6,0.6,0.6,0.3};
/*[2][4] = {
	{3.0,1.0,1.0,0.4},
	{2.0,0.4,0.4,0.2}
};*/
/*
//Alpha
const int alphaxres = 100;
const float alphaxmin = 0.0;
const float alphaxmax = 3.2;
//Binning
const float Wmin = 1.4;
const float Wres = 0.5;
const float Wmax = 2.125;
const float Q2min = 2.0;
const float Q2max = 5.0;
const float Q2res = 0.5;
//CC Min
const float MinCCmin = -0.5;
const float MinCCmax = 501.5;
const int MinCCres = 502;
*/
//Project
const int p_MMxres = 400;
const float p_MMxmin = -0.2;
const float p_MMxmax = 8.0;

//binning
const float Wbin_res = 0.025;//The width of a W bin //30 steps
const float Wbin_start = 1.4;//The starting of W bins

const float Q2bin_res = 0.5;//6 steps
const float Q2bin_start = 2.0; 

const float pbin_res = 0.5;//Range: 0-6.0
const float pbin_start = 0.5;


//Fun names for file lists
static const std::string list1 = "one";
static const std::string list2 = "two";
static const std::string list3 = "three";
static const std::string list4 = "four";
static const std::string list3p = "three+";
static const std::string list3n = "three-";
static const std::string list3h = "threeh";
static const std::string lists1 = "sim_e16_1";
static const std::string lists2 = "sim_pre";
static const std::string lists3n = "sim_gsim_ngpp_nmcdata";
static const std::string lists3y = "sim_gsim_ngpp_mcdata";
static const std::string lists3 = "sim_pair";

//Paths for file names
static const std::string path1 = "/home/mclauchlinc/Desktop/analysis/nick_convert_e16.txt";
static const std::string path2 = "/Users/cmc/Desktop/analysis/e16_10_18_17_ntple.txt";
static const std::string path3 = "/Users/cmc/Desktop/analysis/Nick_skim_e16.txt";
static const std::string path4 = "/Users/cmc/Desktop/analysis/arjun_sim.txt";
static const std::string path3p = "/Users/cmc/Desktop/analysis/NickSkim_e16_PlateIN.txt";
static const std::string path3n = "/Users/cmc/Desktop/analysis/NickSkim_e16_PlateOUT.txt";
static const std::string paths1 = "/Users/cmc/Desktop/analysis/simulation/sim_e16_group1.txt";
static const std::string paths2 = "/Users/cmc/Desktop/analysis/simulation/sim_e16_pre_gpp.txt";
static const std::string paths3n = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/gsim_no_gpp_nmcdata.txt";
static const std::string paths3y = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/gsim_no_gpp_mcdata.txt";
static const std::string paths3 = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/full_sim_pair.txt";

//Run type for file lists {e16,e1f,e16sim,e1fsim} -> {1,2,3,4}
static const int type1 = 1;
static const int type2 = 1;
static const int type3 = 1;
static const int type4 = 3;
static const int type3p = 1;
static const int type3n = 1; 
static const int types1 = 3; 
static const int types2 = 5; 
static const int types3n = 5; 
static const int types3y = 5; 
static const int types3 = 1; 

//Map of file lists to run type
static std::unordered_map<std::string, std::string> filepath_map = 	{	{list1,path1},
																		{list2,path2},
																		{list3,path3},
																		{list4,path4},
																		{list3p,path3p},
																		{list3n,path3n},
																		{lists1,paths1},
																		{lists2,paths2},
																		{lists3n,paths3n},
																		{lists3y,paths3y},
																		{lists3,paths3}};

static std::unordered_map<std::string, int> filetype_map = 	{	{list1,type1},
																{list2,type2},
																{list3,type3},
																{list4,type4},
																{list3p,type3p},
																{list3n,type3n},
																{lists1,types1},
																{lists2,types2},
																{lists3n,types3n},
																{lists3y,types3y},
																{lists3,types3}};

static const char * species[] = {"ele","pro","pip","pim"};//4
static const char * eid_cut[] = {"pre","sanity","fid","sf","min_cc","fid+sf","fid+cc","sf+cc","eid","bank","event"}; //11
static const char * cut_ver[] = {"cut","anticut"};
static const char * hid_cut[] = {"pre","sanity","fid","dt","hid","bank","event"};//,"pWQ2"}; //7
static const char * topologies[] = {"None","Pmiss","PIPmiss","PIMmiss","Zeromiss","ALLmiss"}; //6
static const char * sec_list[] = {"all_sectors","sec1","sec2","sec3","sec4","sec5","sec6"}; //7`
static const char * W_dep_list[] = {"No_W_Dep","W_Dep"};
static const char * CC_det_side[] = {"Left","Coince","Right","All"};
static const char * basic_cut[] = {"pre","cut","anti"};
static const char * MM_sq[] = {"linear","squared"};
static const char * par_cut[4][11] = {{"pre","sanity","fid","sf","min_cc","fid+sf","fid+cc","sf+cc","eid","bank","event"},{"pre","sanity","fid","dt","hid","bank","event"},{"pre","sanity","fid","dt","hid","bank","event"},{"pre","sanity","fid","dt","hid","bank","event"}};
static const char * fit_q[] = {"4fit","4show"};
static const char * thrown[] = {"recon","thrown"};//For reconstructed vs. thrown events

#endif
