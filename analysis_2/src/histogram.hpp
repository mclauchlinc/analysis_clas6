#ifndef HISTOGRAM_H_GUARD
#define HISTOGRAM_H_GUARD

#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "constants.hpp"
#include "functions.hpp"
#include "CartesianGenerator.hpp"
#include "physics.hpp"
#include "detectors.hpp"
//#include "particle.hpp"
//#include "variables.h"
//#include "CartesianGenerator.hh"


using TH2F_ptr = std::shared_ptr<TH2F>;
using TH1F_ptr = std::shared_ptr<TH1F>;
using THn_ptr = std::shared_ptr<ThnSparseF>;


class Histogram {
protected:
	std::shared_ptr<TFile> RootOutputFile;
	TCanvas* def; 

	//Plot Formation Constants
	//W Q2
	 int WQxres = 300;
	 int WQyres = 300;
	 double WQxmin = -0.01;
	 double WQymin = -0.01;
	 double WQxmax = 3.99;
	 double WQymax = 8.99;
	//For Topology Selection
	 int WQ2xres = 20;
	 int WQ2yres = 20;

	//Electron Sampling Fraction
	 int SFxres = 300;
	 int SFyres = 200;
	 double SFxmin = 0.0;
	 double SFymin = 0.0;
	 double SFxmax = 6.0;
	 double SFymax = 1.0;
	//Minimum Electron Energy
	 int MEExres = 100;
	 int MEEyres = 100;
	 double MEExmin = 0.0;
	 double MEEymin = 0.0;
	 double MEExmax = 6.0;
	 double MEEymax = 1.0;
	//Fiducial Cuts
	 int FIDxres = 200;
	 int FIDyres = 200;
	 double FIDxmin = -30.0;
	 double FIDymin = 0.0;
	 double FIDxmax = 30.0;
	 double FIDymax = 180.0;
	//Delta_t
	 int DTxres = 100;//300;
	 int DTyres = 200;
	 double DTxmin = 0.0;
	 double DTymin = -5.0;//-4.0;
	 double DTxmax = 4.5;//7.0;
	 double DTymax = 5.0;//4.0;
	//Missing Mass
	 int MMxres = 1500;
	 double MMxmin = -0.2;
	 double MMxmax = 3.0;
	//Alpha
	 int alphaxres = 100;
	 double alphaxmin = 0.0;
	 double alphaxmax = 3.2;
	//Binning for Analysis
	 double Wmin = 1;
	 double Wres = 0.5;
	 double Wmax = 3;
	 double Q2min = 1.5;
	 double Q2max = 5.0;
	 double Q2res = 0.5;
	 //CC Min
	 double MinCCmin = -0.5;
	 double MinCCmax = 501.5;
	 int MinCCres = 502;

	//binning
	int bins[3][6] = {{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}}; //{W, Q2, MM, theta, alpha, helicity}
	float xmin[3][6] = {{NAN,NAN,NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN,NAN,NAN}};
	float xmax[3][6] = {{NAN,NAN,NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN,NAN,NAN}Zaq	ww5678\}; 


	 float Wbin_res = 0.025;//The width of a W bin //30 steps
	 float Wbin_start = 1.4;//The starting of W bins

	 float Q2bin_res = 0.5;//6 steps
	 float Q2bin_start = 2.0; 

	double Yth_start = 9.0;
	double Yth_res = 18.0;
	double Yal_start = 18.0;//10.0; 0-36, 36-72, 72-108,  
	double Yal_res = 36.0;//40.0;
	double YM_start[3] = {1.1,1.1,0.3};
	double YM_res[3] = {0.06,0.06,0.06};

	 //Making the Histograms
	THn_ptr egg[4][3]; //{e16,e1f,e16sim,e1fsim},{pim,p,pip}
	//projections
	TH1F_ptr hist_mm[3][]
	


public:
	Histogram(const std::string& output_file);
	~Histogram();
	void Write();
	//W Qsquared plots
	void Set_THn(); 
	//int Set_Bins(int bin, std::vector<float> bin);
};





#endif