#include "friend_hist.hpp"


F_Hist::F_Hist(std::shared_ptr<Environment> envi_, std::shared_ptr<THnSparseD> hist_[3]){
	//TH2F _WQ2_plots[_bins[0]][_bins[3]][_bins[4]][_bins[5]][3];//the 3 is for the different channels
	//TH1F _MM_plots[_bins[0]][_bins[1]][_bins[2]][3];
	//TH1F _Theta_plots[_bins[0]][_bin[1]][_bins[2]][3];
	//TH1F _Alpha_plots[_bins[0]][_bin[1]][_bins[2]][3];
	//Histogram::Friend_bin_reverse(var,bin,channel)
	//Wq2 Plots
	Int_t *x; //Getting the bin really only wants an Int_t * as an input and I cannot for the life of me figure out another way of doing this
	for(int f = 0; f < 3; f++){//Channel
		for(int i = 0; i < xc_bins[0]; i++){//Top
			x[0] = i; 
			for(int j = 0; j < xc_bins[1]; j++){//W
				x[1] = j;
				for(int k = 0; k < xc_bins[2]; k++){//Q2
					x[2] = k;
					for(int l = 0; l < xc_bins[3]; l++){//MM
						x[3] = l;
						for(int m = 0; m < xc_bins[4]; m++){//theta
							x[4] = m;
							for(int n = 0; n < xc_bins[5]; n++){//alpha
								x[5] = n;
								for(int o = 0; o < xc_bins[6]; o++){//phi
									x[6] = o; 
									//_WQ2_plots[i][l][m][n][f]->Fill(hist_[f]->GetBinContent({i,j,k,l,m,n,o}),hist_[f]->Histogram::Friend_bin_reverse(2,j,f));
									_MM_plots[i][j][k][f]->Fill(hist_[f]->GetBinContent(hist_[f]->GetBin(x)));
									_Theta_plots[i][j][k][f]->Fill(hist_[f]->GetBinContent(hist_[f]->GetBin(x)));
									_Alpha_plots[i][j][k][f]->Fill(hist_[f]->GetBinContent(hist_[f]->GetBin(x)));
								}
							}
						}
					}
				}
			}
		}
	}
	
	//MM Plots

	//Theta Plots

	//Alpha Plots

}
