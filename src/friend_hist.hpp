#ifndef FRIEND_HIST_H_GUARD
#define FRIEND_HIST_H_GUARD


#include "histogram.hpp"
#include "constants.hpp"
#include "environment.hpp"

using TH2F_ptr = std::shared_ptr<TH2F>;
using TH1F_ptr = std::shared_ptr<TH1F>;
using THn_ptr = std::shared_ptr<THnSparseD>;

class F_Hist{
private:
	//static Int_t _bins[7] = {5,29,5,14,10,10,10};//Be sure to check that these are identical to the ones in _Friend_bins in histograms.hpp
	///topology, W, Q2, MM, Theta, Alpha, Phi
	/*
	Topology: [_bin[0]]
	W: [_bin[1]]
	Q2: [_bin[2]]
	MM: [_bin[3]]
	theta: [_bin[4]]
	alpha: [_bin[5]]
	phi: [_bin[6]]
	[_bins[0]][_bin[1]][_bin[2]][_bins[3]][_bins[4]][_bins[5]][_bins[6]];
	*/
	//The Variable Plots should only depend on W and Qsquared, and be integrated over the others 
	TH2F_ptr _WQ2_plots[5][14][10][10][3];//the 3 is for the different channels
	TH1F_ptr _MM_plots[5][29][5][3];
	TH1F_ptr _Theta_plots[5][29][5][3];
	TH1F_ptr _Alpha_plots[5][29][5][3];
	//TH1F _W_plots[_bins[0]][_bin[2]][_bins[3]][_bins[4]][_bins[5]][_bins[6]];
	//TH1F _WQ2_plots[_bins[0]][_bin[1]][_bins[3]][_bins[4]][_bins[5]][_bins[6]];
public: 
	F_Hist(std::shared_ptr<Environment> envi_, std::shared_ptr<THnSparseD> hist_[3]);//Channel being which is our "special" particle
};

#endif