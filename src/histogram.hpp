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
#include "environment.hpp"
#include "THnSparse.h"
//#include "particle.hpp"
//#include "variables.h"
//#include "CartesianGenerator.hh"


using TH2F_ptr = std::shared_ptr<TH2F>;
using TH1F_ptr = std::shared_ptr<TH1F>;
using THn_ptr = std::shared_ptr<THnSparseD>;
//using TH1I_ptr = std::shared_ptr<TH1I>;


class Histogram {
protected:
	std::shared_ptr<TFile> RootOutputFile;
	TCanvas* def; 

	//Plot Formation Constants
	//W Q2
	 int WQxres = 200;
	 int WQyres = 200;
	 double WQxmin = -0.01;
	 double WQymin = -0.01;
	 double WQxmax = 3.99;
	 double WQymax = 8.99;
	//For Topology Selection
	 int WQ2xres = 20;
	 int WQ2yres = 20;

	//Electron Sampling Fraction
	 int SFxres = 500;
	 int SFyres = 300;
	 double SFxmin = 0.0;
	 double SFymin = 0.0;
	 double SFxmax = 6.0;
	 double SFymax = 1.0;
	//Minimum Electron Energy
	 int MEExres = 300;
	 int MEEyres = 300;
	 double MEExmin = 0.0;
	 double MEEymin = 0.0;
	 double MEExmax = 6.0;
	 double MEEymax = 1.0;
	//Fiducial Cuts
	 int FIDxres = 400;
	 int FIDyres = 300;
	 double FIDxmin = -30.0;
	 double FIDymin = 0.0;
	 double FIDxmax = 30.0;
	 double FIDymax = 180.0;
	//Delta_t
	 int DTxres = 1000;//300;
	 int DTyres = 400;
	 double DTxmin = 0.0;
	 double DTymin = -10.0;//-4.0;
	 double DTxmax = 4.5;//7.0;
	 double DTymax = 10.0;//4.0;
	//Missing Mass
	 int MMxres = 400;
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
	TH2F_ptr WQ2_hist[11][6][2];//electron cuts, topologies (including pre), Recon vs. Raw (for data this should always be "Recon")
	TH2F_ptr Fid_hist[7][4][11][30][12][6][2];//sector, species, cut, W binning, p binning, topology, anti
	TH2F_ptr SF_hist[10][30][7][6][2];//cuts, W Binning, Sector, topology
	TH2F_ptr DT_hist[4][7][30][7][6][2]; //particle, cuts, W binning, sector, topology
	TH1F_ptr CC_hist[6][18][11][4][6][2]; //Sector, segment, cut, side of detector, topology, anti
	TH1F_ptr MM_hist[4][3][2][2];//topology, cut, squared v linear, fitting vs. not fitting plots
	TH1F_ptr Cross_hist[2]; //Showing how many events counted in mulitiple topologies

	bool Fid_made_hist[7][4][11][30][12][6][2];
	bool Fid_fill_hist[7][4][11][30][12][6][2];
	bool Fid_write_hist[7][4][11][30][12][6][2];

	bool CC_made_hist[6][18][11][4][6][2];
	bool CC_fill_hist[6][18][11][4][6][2];
	bool CC_write_hist[6][18][11][4][6][2];

	bool DT_made_hist[4][7][30][7][6][2];//Added one to the second bin for cuts to all for the electron WQ2
	bool DT_fill_hist[4][7][30][7][6][2];
	bool DT_dir_hist[4][7][30][7][6][2];
	bool DT_dir_made[4][8][2][8][6];

	bool WQ2_made_hist[11][6][2];
	bool WQ2_dir_made[11][6][2];

	THn_ptr Friend[3];//This will be the 7 dimensional histogram from which I can project out different pieces
	Int_t _Friend_bins[7] = {5,29,5,14,10,10,10}; //topology, W, Q2, MM, Theta, Alpha, Phi
	float _W_min = 1.4;
	float _W_max = 2.125;
	float _Q2_min = 2.0;
	float _Q2_max = 5.0;
	float _MM_min[3] = {1.1,0.3,1.1};
	float _MM_max[3] = {2.0,1.2,2.0};
	float _theta_min = 0.0;
	float _theta_max = 180.0;
	float _alpha_min = 0.0;
	float _alpha_max = 360.;
	float _phi_min = 0.0; 
	float _phi_max = 360.0;



public:
	Histogram(std::shared_ptr<Environment> _envi, const std::string& output_file);
	//~Histogram();
	void Write(std::shared_ptr<Environment> _envi);
	//W Qsquared plots
	int W_binning(float W_);
	int p_binning(float p_);
	char Part_cut(int species, int cut);
	void WQ2_Make(std::shared_ptr<Environment> _envi);
	void WQ2_Fill(std::shared_ptr<Environment> _envi, int top, int cut, float W_, float Q2_, int thr = 0);
	void WQ2_Write(std::shared_ptr<Environment> _envi);
	//Fiducial Cuts
	void Fid_Make(std::shared_ptr<Environment> _envi );
	void Fid_Fill(std::shared_ptr<Environment> _envi, int top, float theta, float phi, int part, int cut, int cutvanti, float W_, float p);
	//void Fid_Fill(std::shared_ptr<Environment> _envi, Particle par_, float W_);
	void Fid_Write(std::shared_ptr<Environment> _envi);
	//Sampling Fraction Cuts
	void SF_Make(std::shared_ptr<Environment> _envi );
	void SF_Fill(std::shared_ptr<Environment> _envi, int top, float p, float en, int cut, int cva, float W_, int sec);
	//void SF_Fill(std::shared_ptr<Environment> _envi, Particle par_, float W_);
	void SF_Write(std::shared_ptr<Environment> _envi);
	//Delta T Cuts
	void DT_Make(std::shared_ptr<Environment> _envi );
	void DT_Fill(std::shared_ptr<Environment> _envi, int top, int part, float p, float d, float t, float d0, float t0, int cut, int anti, float W_, int sec);
	void DT_Fill(std::shared_ptr<Environment> _envi, int top, int part, float p, float dt, int cut, int anti, float W_, int sec);
	//void DT_Fill(std::shared_ptr<Environment> _envi, Particle par_, float W_);
	void DT_Write(std::shared_ptr<Environment> _envi);
	//Min CC Cuts
	void CC_Make(std::shared_ptr<Environment> _envi );
	void CC_Fill(std::shared_ptr<Environment> _envi, int top, int sec, int segm, int nphe, int cut, int anti);
	//void CC_Fill(std::shared_ptr<Environment> _envi, Particle par_);
	void CC_Write(std::shared_ptr<Environment> _envi);
	//Missing Mass Cuts
	void MM_Make(std::shared_ptr<Environment> _envi );
	void MM_Fill(std::shared_ptr<Environment> _envi, int top, float mm, int cut, int square, bool fit);
	//void MM_Fill(std::shared_ptr<Environment> _envi, Particle _par);
	void MM_Write(std::shared_ptr<Environment> _envi);

	void Friend_Make(std::shared_ptr<Environment> _envi);
	int Friend_W_binning(float W_);
	int Friend_Q2_binning(float Q2_);
	int Friend_MM_binning(float MM_, int chan);
	int Friend_theta_binning(float theta_);
	int Friend_alpha_binning(float alpha_);
	int Friend_phi_binning(float phi_);
	int * Friend_binning(int top, float W_, float Q2_, float MM_, float theta_, float alpha_, float phi_ , int channel);
	void Friend_Fill(std::shared_ptr<Environment> _envi, int top_, float W_, float Q2_, float MM_, float theta_, float alpha_, float phi_ , int chan_, float weight_);
	void Friend_Write(std::shared_ptr<Environment> _envi);

	void Cross_Make(std::shared_ptr<Environment> envi_);
	void Cross_Fill(std::shared_ptr<Environment> envi_, int gevnt_[4], float weight_);
	void Cross_Write(std::shared_ptr<Environment> envi_);

	//void Event_Particle_Hist(std::shared_ptr<Environment> envi_, const Particle p1, float W_, int top_, int par_, bool pass_);


	//**Include functions to fill histograms from events and particle




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