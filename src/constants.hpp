#ifndef CONSTANTS_H_GUARD
#define CONSTANTS_H_GUARD

#include <unordered_map>
#include "TMath.h"
#include "TLorentzVector.h"


static const int MAX_PARTS = 100; 
static const int NUM_THREADS = 4;

static const double c_special = 29.9792458; //speed of light in cm/ns
static const double c_convert = 10000000; //Convert c_special to m/s

//Beam energies in GeV
static const double energy_e16 = 5.7696;
//3_14 was 5.754, but led to bad results: updated to 5.759 after reading Paremuzyan's thesis 
//Check 4.794 GeV for e16
static const double energy_e1f = 5.499;

static const double fine_structure = 1.0/137.0; 

//Masses of relevant particles
static const	double me = 0.0005109989; //mass of electron in GeV
static const	double mp = 0.93828;	//Mass of proton in GeV
static const	double mpi = 0.1395;	//Mass of pion in GeV

//Range of W/Q2
static const double WminAna = 1.4;
static const double WmaxAna = 2.125;
static const double Q2minAna = 3.0;//2.0; Moved up to 3.0 in order to compare with simulation in that range
static const double Q2maxAna = 5.0; 

//Initial Four Vectors

//TLorentzVector k_mu_e16 = Make_4Vector(energy_e16,0.0,0.0,1.0,me);
//TLorentzVector k_mu_e1f = Make_4Vector(energy_e1f,0.0,0.0,1.0,me);
//TLorentzVector p_mu = Make_4Vector(0.0,0.0,0.0,0.0,mp);

//e1-6 Luminosity Values
static const double lt_e16 = 5.0; //Target length in cm
static const double Dt_e16 = 0.073; //Density of target in g/cm^3
static const double NA = 6.022 * pow(10.0,23); //Avogadro's number
static const double qe = 1.602 * pow(10.0,-19); // fundamental Coulomb charge 
static const double Mt_e16 = 1.007; //Molar mass of target in g/mole
static const double Qt_e16 = 21.32*pow(10.0,-3); //Total charge incident on target from Arjun in Coulombs //Will need to verify myself
static const double L_e16 = Qt_e16*lt_e16*Dt_e16*NA/(qe*Mt_e16);
static const double cm2_to_mbarn = pow(10.0,-27);

//Particle ID Numbers
//anti particles are negative
static const 	int PROTON = 2212;
static const	int ELECTRON = 11;
static const	int PION = 211;
static const	int PION_0 = 111;

// Beam spot from Arjun's Fortran vertex_e1f.f
static const double x_beam = 0.090;
static const double y_beam = -0.345;


//Fun names for file lists
static const std::string list1 = "one";
static const std::string list2 = "two";
static const std::string list3 = "three";
static const std::string list4 = "four";
static const std::string list3p = "three+";
static const std::string list3n = "three-";

//Paths for file names
static const std::string path1 = "/home/mclauchlinc/Desktop/analysis/nick_convert_e16.txt";
static const std::string path2 = "/Users/cmc/Desktop/analysis/e16_10_18_17_ntple.txt";
static const std::string path3 = "/Users/cmc/Desktop/analysis/Nick_skim_e16.txt";
static const std::string path4 = "/Users/cmc/Desktop/analysis/arjun_sim.txt";
static const std::string path3p = "/Users/cmc/Desktop/analysis/NickSkim_e16_PlateIN.txt";
static const std::string path3n = "/Users/cmc/Desktop/analysis/NickSkim_e16_PlateOUT.txt";

//Run type for file lists {e16,e1f,e16sim,e1fsim} -> {1,2,3,4}
static const int type1 = 1;
static const int type2 = 1;
static const int type3 = 1;
static const int type4 = 3;
static const int type3p = 1;
static const int type3n = 1; 

//Map of file lists to run type
static std::unordered_map<std::string, std::string> filepath_map = 	{{list1,path1},
																					{list2,path2},
																					{list3,path3},
																					{list4,path4},
																					{list3p,path3p},
																					{list3n,path3n}};

static std::unordered_map<std::string, int> filetype_map = 	{{list1,type1},
																			{list2,type2},
																			{list3,type3},
																			{list4,type4},
																			{list3p,type3p},
																			{list3n,type3n}};

#endif
