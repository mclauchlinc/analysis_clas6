#ifndef HISTOGRAM_H_GUARD
#define HISTOGRAM_H_GUARD

#include "TH1.h"
#include "TH2.h"
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
	TH2F_ptr WQ2_hist[11][6];//electron cuts, topologies (including pre)
	TH2F_ptr Fid_hist[7][4][11][30][12][6][2];//sector, species, cut, W binning, p binning, topology
	TH2F_ptr SF_hist[10][30][7][6][2];//cuts, W Binning, Sector, topology
	TH2F_ptr DT_hist[3][7][30][7][6][2]; //particle, cuts, W binning, sector, topology
	TH1F_ptr CC_hist[6][18][11][4][6][2]; //Sector, segment, cut, side of detector, topology
	TH1F_ptr MM_hist[5][3][2];//topology, cut, squared v linear

	bool Fid_made_hist[7][4][11][30][12][6][2];
	bool Fid_fill_hist[7][4][11][30][12][6][2];
	bool Fid_write_hist[7][4][11][30][12][6][2];

	bool CC_made_hist[6][18][11][4][6][2];
	bool CC_fill_hist[6][18][11][4][6][2];
	bool CC_write_hist[6][18][11][4][6][2];

	bool DT_made_hist[3][7][30][7][6][2];
	bool DT_fill_hist[3][7][30][7][6][2];
	bool DT_dir_hist[3][7][30][7][6][2];
	bool DT_dir_made[3][8][2][8][6];

	bool WQ2_made_hist[11][6];
	bool WQ2_dir_made[11][6];


public:
	Histogram(const std::string& output_file);
	~Histogram();
	void Write();
	//W Qsquared plots
	int W_binning(float W_);
	int p_binning(float p_);
	char Part_cut(int species, int cut);
	void WQ2_Make();
	void WQ2_Fill(int top, int cut, float W_, float Q2_);
	void WQ2_Write();
	//Fiducial Cuts
	void Fid_Make();
	void Fid_Fill(int top, float theta, float phi, int part, int cut, int cutvanti, float W_, float p);
	void Fid_Write();
	//Sampling Fraction Cuts
	void SF_Make();
	void SF_Fill(int top, float p, float en, int cut, int cva, float W_, int sec);
	void SF_Write();
	//Delta T Cuts
	void DT_Make();
	void DT_Fill(int top, int part, float p, float d, float t, float d0, float t0, int cut, int anti, float W_, int sec);
	void DT_Fill(int top, int part, float p, float dt, int cut, int anti, float W_, int sec);
	void DT_Write();
	//Min CC Cuts
	void CC_Make();
	void CC_Fill(int top, int sec, int segm, int nphe, int cut, int anti);
	void CC_Write();
	//Missing Mass Cuts
	void MM_Make();
	void MM_Fill(int top, float mm, int cut, int square);
	void MM_Write();
	/*//Signature Plots //We'll get there
	void MM2_Make();
	void MM2_Fill(int top, float mm, float W_, float Q2_);
	void MM2_Write();
	void Th2_Make();
	void Th2_Fill(int top, float theta, float W_, float Q2_);
	void Th2_Write();
	void Alpha_Make();
	void Alpha_Fill(int top, float alpha, float W_, float Q2_);
	void Alpha_Write();
	*/
	//void Fill_EID(std::shared_ptr<Particle> par);
	//void Fill_HID(std::shared_ptr<Particle>  par);
};





#endif