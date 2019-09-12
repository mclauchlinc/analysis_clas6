#include "histogram.hpp"


Histogram::Histogram(const std::string& output_file){
	RootOutputFile = fun::Name_File(output_file);
	def = new TCanvas("def");
	Histogram::WQ2_Make();
	Histogram::Fid_Make();
	Histogram::SF_Make();
	//Histogram::DT_Make();
	//Histogram::CC_Make();
	//Histogram::MM_Make();
}

Histogram::~Histogram() { this->Write(); }

void Histogram::Write(){
	std::cout<< "Writing" <<std::endl;
	RootOutputFile->cd();
	WQ2_Write();
	Fid_Write();
	SF_Write();
	//DT_Write();
	//CC_Write();
	//MM_Write();
	RootOutputFile->Close();
}


//W Qsquared plots
int Histogram::W_binning(float W_){
  int j = 0;
  float top, bot; 
  for(int i = 1; i < 30; i++){
    top = Wbin_start + i*Wbin_res;//constants.hpp
    bot = top - Wbin_res; 
    if(W_ < top && W_ > bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::p_binning(float p_){
  int j = 0;
  float top, bot; 
  for(int i = 1; i < 12; i++){
    top = pbin_start + i*pbin_res;//constants.hpp
    bot = top - pbin_res; 
    if(p_ < top && p_ > bot){
      j = i; 
    }
  }
  return j; 
}

char Histogram::Part_cut(int species, int cut){
	char the_cut[100]; 
	switch(species){
		case 0:
		sprintf(the_cut,"%s",eid_cut[cut]);
		break;
		case 1:
		if(cut >= 6){
			sprintf(the_cut,"INVALD_CUT");
		}else{
			sprintf(the_cut,"%s",hid_cut[cut]);
		}
		break;
	}
	//return the_cut; 
}

void Histogram::WQ2_Make(){
	std::vector<long> space_dims(2);
	space_dims[0] = 11; //Electron Cuts
	space_dims[1] = 6; //Topologies

	CartesianGenerator cart(space_dims);//CartesianGenerator.hpp
	char hname[100]; 

	while(cart.GetNextCombination()){
		sprintf(hname,"W_Q2_%s_%s",eid_cut[cart[0]],topologies[cart[1]]); //constants.h and otherwise writing the specific cut to the right plot
    	WQ2_hist[cart[0]][cart[1]] = std::make_shared<TH2D>( hname, hname, WQxres, WQxmin, WQxmax, WQyres, WQymin, WQymax); // constants.h
	}
}

void Histogram::WQ2_Fill(int top, int cut, float W_, float Q2_){
	WQ2_hist[cut][top]->Fill(W_,Q2_);
}

void Histogram::WQ2_Write(){
	TDirectory* dir_WQ2 = RootOutputFile->mkdir("W vs. Q2");
	dir_WQ2->cd();
	for(int top = 0; top < 6; top ++){
		for(int cut = 0; cut < 10; cut++){
			WQ2_hist[top][cut]->SetXTitle("W (GeV)");
			WQ2_hist[top][cut]->SetYTitle("Q^{2} (GeV^{2}");
			WQ2_hist[top][cut]->SetOption("Colz");
			WQ2_hist[top][cut]->Write();
		}
	}
}

//Fiducial Cuts
void Histogram::Fid_Make(){
	char hname[100];
	char par_cut[100]; 
	float wtop, wbot,ptop,pbot; 
	std::vector<long> space_dims(6);
	space_dims[0] = 7; //Sector
	space_dims[1] = 4; //species
	space_dims[2] = 11;//cuts
	space_dims[3] = 30;//W Binning
	space_dims[4] = 12;//Momentum binning
	space_dims[5] = 6; //Topology

	CartesianGenerator cart(space_dims);

	//char * fid_cuts[];
	int length = 0;
	char p_cut[100]; 
	int num_hist = 0; 
	//std::cout<<std::endl <<"Here we are. Making Fiducial histograms" <<std::endl;

	while(cart.GetNextCombination()){
		//std::cout<<"In the loop: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<" Histogram Number: " <<num_hist <<std::endl;
		if(cart[1] == 0){ //Dealing with the different numbers of cuts for electrons vs. hadrons
			//fid_cuts = eid_cut;
			length = 10; 
		}else{
			//fid_cuts = hid_cut; 
			length = 6; 
		}
		if(cart[4] ==0){
			sprintf(p_cut,"p_range:All");
		}else{
			ptop = pbin_start + (cart[4]*pbin_res);//constants.hpp
			pbot = ptop - pbin_res;//constants.hpp
			sprintf(p_cut,"p_range:%f-%f",pbot,ptop);
		}
		if(cart[2]<length){//Here we go
			if(cart[1]==0){
				sprintf(par_cut,"%s",eid_cut[cart[2]]);
			}else{
				sprintf(par_cut,"%s",hid_cut[cart[2]]);
			}
			if(cart[3] == 0){//All W bins
				if(cart[5] == 0 && cart[1] != (length-1) && cart[4]!=0){//Specific Momentum Bins
					if((cart[1]!=0 && (cart[1] == 0 || cart[1]==3)) || (cart[1]==0 && (cart[1]==0 || cart[1] == 7))){//Isolated cuts for particles
						sprintf(hname,"%s_Fid_%s_%s_W:ALL_%s_%s",species[cart[1]],sec_list[cart[0]],par_cut,p_cut,topologies[cart[5]]);//constants.hpp
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = std::make_shared<TH2D>(hname,hname, FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
						num_hist++;
					}
				}
				if(cart[4]==0 && (((cart[1] != (length-1) || cart[5]==0)) || (cart[1]== (length-1) && cart[5]!=0 ))) {//No momentum dependence and correct event stuff
					sprintf(hname,"%s_Fid_%s_%s_W:ALL_%s_%s",species[cart[1]],sec_list[cart[0]],par_cut,p_cut,topologies[cart[5]]);//constants.hpp
					Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = std::make_shared<TH2D>(hname,hname, FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
					num_hist++;
				}
			}else{//Specific W Bins
				if(cart[5]==0 && cart[4]==0 && cart[1] != (length-1) ){ //No event selection or momentum dependence
					if((cart[1]!=0 && (cart[1] == 0 || cart[1]==3)) || (cart[1]==0 && (cart[1]==0 || cart[1] == 7))){//Isolated cuts for particles
						wtop = Wbin_start + (cart[3]*Wbin_res);//constants.hpp
						wbot = wtop - Wbin_res;//constants.hpp
						sprintf(hname,"%s_Fid_%s_%s_W:%f-%f_%s_%s",species[cart[1]],sec_list[cart[0]],par_cut,wbot,wtop,p_cut,topologies[cart[5]]);//constants.hpp
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = std::make_shared<TH2D>(hname,hname, FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
						num_hist++;
					}
				}
			}
		}
	}
}

void Histogram::Fid_Fill(int top, float theta, float phi, int part, int cut, float W_, float p_){
	float phic = physics::phi_center(phi);
	if(p_binning(p_) + W_binning(W_) != 0){//Both p and W dependence
		Fid_hist[0][part][cut][Histogram::W_binning(W_)][Histogram::p_binning(p_)][top]->Fill(phic,theta);//All Sectors
		Fid_hist[physics::get_sector(phi)][part][cut][Histogram::W_binning(W_)][Histogram::p_binning(p_)][top]->Fill(phic,theta);//Individual Sectors
	}
	if(p_binning(p_) != 0){ //All W with p binning
		Fid_hist[0][part][cut][0][Histogram::p_binning(p_)][top]->Fill(phic,theta);//All Sectors
		Fid_hist[physics::get_sector(phi)][part][cut][0][Histogram::p_binning(p_)][top]->Fill(phic,theta);//Individual Sectors
	}
	if(W_binning(W_) != 0){ //All p with W binning
		Fid_hist[0][part][cut][Histogram::W_binning(W_)][0][top]->Fill(phic,theta);//All sectors
		Fid_hist[physics::get_sector(phi)][part][cut][Histogram::W_binning(W_)][0][top]->Fill(phic,theta);//Individual Sectors
	}
	//All W and All Momentum
	Fid_hist[0][part][cut][0][0][top]->Fill(phic,theta);//All sectors, all W, all P
	Fid_hist[physics::get_sector(phi)][part][cut][0][0][top]->Fill(phic,theta);//Individual sectors all W, all p
}

void Histogram::Fid_Write(){
	TDirectory* dir_Fid = RootOutputFile->mkdir("Fiducial Plots");
	dir_Fid->cd();
	TDirectory* ele_fid = dir_Fid->mkdir("Electron Fiducial Plots");
	TDirectory* pro_fid = dir_Fid->mkdir("Proton Fiducial Plots");
	TDirectory* pip_fid = dir_Fid->mkdir("Pi+ Fiducial Plots");
	TDirectory* pim_fid = dir_Fid->mkdir("Pi- Fiducial Plots");
	TDirectory* ele_fid_w = ele_fid->mkdir("Electron Fid W-Dependence");
	TDirectory* ele_fid_p = ele_fid->mkdir("Electron Fid P-Dependence");
	TDirectory* ele_fid_wp = ele_fid->mkdir("Electron Fid W-P-Dependence");
	TDirectory* pro_fid_w = pro_fid->mkdir("Proton Fid W-Dependence");
	TDirectory* pro_fid_p = pro_fid->mkdir("Proton Fid P-Dependence");
	TDirectory* pro_fid_wp = pro_fid->mkdir("Proton Fid W-P-Dependence");
	TDirectory* pip_fid_w = pip_fid->mkdir("Pi+ Fid W-Dependence");
	TDirectory* pip_fid_p = pip_fid->mkdir("Pi+ Fid P-Dependence");
	TDirectory* pip_fid_wp = pip_fid->mkdir("Pi+ Fid W-P-Dependence");
	TDirectory* pim_fid_w = pim_fid->mkdir("Pi- Fid W-Dependence");
	TDirectory* pim_fid_p = pim_fid->mkdir("Pi- Fid P-Dependence");
	TDirectory* pim_fid_wp = pim_fid->mkdir("Pi- Fid W-P-Dependence");

	
	std::vector<long> space_dims(6);
	space_dims[0] = 7; //Sector
	space_dims[1] = 4; //species
	space_dims[2] = 10;//cuts
	space_dims[3] = 30;//W Binning
	space_dims[4] = 12;//Momentum binning
	space_dims[5] = 6; //Topology

	CartesianGenerator cart(space_dims);

	//char * fid_cuts[];
	int length = 0;
	char p_cut[100]; 

	bool p_dep = false;
	bool w_dep = false; 

	while(cart.GetNextCombination()){
		p_dep = false;
		w_dep = false;
		if(cart[3] == 0){
			w_dep = false;
		}else {
			w_dep = true; 
		}
		if(cart[4] == 0 || w_dep){
			p_dep = false;
		}else{
			p_dep = true;
		}
		switch(cart[1]){
			case 0:
				ele_fid->cd();
				if(p_dep){
					ele_fid_p->cd();
				}else if(w_dep){
					ele_fid_w->cd();
				}
			break;
			case 1:
				pro_fid->cd();
				if(p_dep){
					pro_fid_p->cd();
				}else if(w_dep){
					pro_fid_w->cd();
				}
			break;
			case 2:
				pip_fid->cd();
				if(p_dep){
					pip_fid_p->cd();
				}else if(w_dep){
					pip_fid_w->cd();
				}
			break;
			case 3:
				pim_fid->cd();
				if(p_dep){
					pim_fid_p->cd();
				}else if(w_dep){
					pim_fid_w->cd();
				}
			break;
		}

		//Just making sure we aren't looking in places that don't exist 
		if(cart[1] == 0){ //Dealing with the different numbers of cuts for electrons vs. hadrons
			//fid_cuts = eid_cut;
			length = 10; 
		}else{
			//fid_cuts = hid_cut; 
			length = 6; 
		}
		if(cart[2]<length){
			if(cart[3] == 0){//All W bins
				if(cart[5] == 0 && cart[1] != (length-1) && cart[4]!=0){//Specific Momentum Bins
					if((cart[1]!=0 && (cart[1] == 0 || cart[1]==3)) || (cart[1]==0 && (cart[1]==0 || cart[1] == 7))){//Isolated cuts for particles
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetXTitle("{phi} (degrees)");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetYTitle("{theta} (degrees)");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetOption("COLZ");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->Write();
					}
				}
				if(cart[4]==0 && (((cart[1] != (length-1) || cart[5]==0)) || (cart[1]== (length-1) && cart[5]!=0 ))) {//No momentum dependence and correct event stuff
					Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetXTitle("{phi} (degrees)");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetYTitle("{theta} (degrees)");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetOption("COLZ");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->Write();
				}
			}else{//Specific W Bins
				if(cart[5]==0 && cart[4]==0 && cart[1] != (length-1) ){ //No event selection or momentum dependence
					if((cart[1]!=0 && (cart[1] == 0 || cart[1]==3)) || (cart[1]==0 && (cart[1]==0 || cart[1] == 7))){//Isolated cuts for particles
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetXTitle("{phi} (degrees)");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetYTitle("{theta} (degrees)");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetOption("COLZ");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->Write();					
					}
				}
			}
		}
	}
}

//Sampling Fraction Cuts
void Histogram::SF_Make(){
	char hname[100];
	std::vector<long> space_dims(4);
	space_dims[0] = 11; //Electron Cuts
	space_dims[1] = 30; //W binning
	space_dims[2] = 7;  //Sector
	space_dims[3] = 6;  //Topology

	float top,bot; 

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		if((cart[0] == 10 && cart[3] != 0) || (cart[0] != 0 && cart[3] == 0) ){//Topology only matters for event selection cut
			if(cart[1] == 0 ){ //All W 
				sprintf(hname,"SF_%s_%s_W:ALL_%s",eid_cut[cart[0]],sec_list[cart[2]],topologies[cart[3]]);
			}else{	//Specific W Bins
				top = Wbin_start + cart[1]*Wbin_res;
				bot = top - Wbin_res;
				sprintf(hname,"SF_%s_%s_W:%f-%f_%s",eid_cut[cart[0]],sec_list[cart[2]],bot,top,topologies[cart[3]]);
			}
		}
		SF_hist[cart[0]][cart[1]][cart[2]][cart[3]] = std::make_shared<TH2D>(hname,hname, SFxres, SFxmin, SFxmax, SFyres, SFymin, SFymax);
	}
}

void Histogram::SF_Fill(int top, float p, float en, int cut, float W_, int sec){
	SF_hist[cut][Histogram::W_binning(W_)][sec][top]->Fill(p,en/p);
	SF_hist[cut][Histogram::W_binning(W_)][0][top]->Fill(p,en/p);
}

void Histogram::SF_Write(){
	TDirectory* dir_SF = RootOutputFile->mkdir("SF Plots");
	dir_SF->cd();
	TDirectory* sf_cut_all = dir_SF->mkdir("SF: All W,Sector, Topologies");
	TDirectory* sf_cut1 = dir_SF->mkdir("SF Pre Cut");
	TDirectory* sf_cut2 = dir_SF->mkdir("SF Sanity Cut");
	TDirectory* sf_cut3 = dir_SF->mkdir("SF Fid Cut");
	TDirectory* sf_cut4 = dir_SF->mkdir("SF SF Cut");
	TDirectory* sf_cut5 = dir_SF->mkdir("SF CC Cut");
	TDirectory* sf_cut6 = dir_SF->mkdir("SF Fid+SF Cut");
	TDirectory* sf_cut7 = dir_SF->mkdir("SF Fid+CC Cut");
	TDirectory* sf_cut8 = dir_SF->mkdir("SF SF+CC Cut");
	TDirectory* sf_cut9 = dir_SF->mkdir("SF EID Cut");
	TDirectory* sf_cut10 = dir_SF->mkdir("SF Bank Cut");
	TDirectory* sf_cut11 = dir_SF->mkdir("SF Event Cut");
	TDirectory* sf_cut1_w = sf_cut1->mkdir("SF Pre Cut W Dep");
	TDirectory* sf_cut2_w = sf_cut2->mkdir("SF Sanity Cut W Dep");
	TDirectory* sf_cut3_w = sf_cut3->mkdir("SF Fid Cut W Dep");
	TDirectory* sf_cut4_w = sf_cut4->mkdir("SF SF Cut W Dep");
	TDirectory* sf_cut5_w = sf_cut5->mkdir("SF CC Cut W Dep");
	TDirectory* sf_cut6_w = sf_cut6->mkdir("SF Fid+SF Cut W Dep");
	TDirectory* sf_cut7_w = sf_cut7->mkdir("SF Fid+CC Cut W Dep");
	TDirectory* sf_cut8_w = sf_cut8->mkdir("SF SF+CC Cut W Dep");
	TDirectory* sf_cut9_w = sf_cut9->mkdir("SF EID Cut W Dep");
	TDirectory* sf_cut10_w = sf_cut10->mkdir("SF Bank Cut W Dep");
	TDirectory* sf_cut11_w = sf_cut11->mkdir("SF Event Cut W Dep");


	std::vector<long> space_dims(4);
	space_dims[0] = 11; //Electron Cuts
	space_dims[1] = 30; //W binning
	space_dims[2] = 7;  //Sector
	space_dims[3] = 6;  //Topology

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		if(cart[1] == 0 && cart[3]==0 && cart[2] == 0){
			sf_cut_all->cd();
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetXTitle("Momentum (GeV)");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetYTitle("Sampling Fraction");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetOption("COLZ");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->Write();
		}else{
			switch(cart[0]){
				case 0:
					sf_cut1->cd();
					if(cart[2] !=0){
						sf_cut1_w->cd();
					}
				break;
				case 1:
					sf_cut2->cd();
					if(cart[2] !=0){
						sf_cut2_w->cd();
					}
				break;
				case 2:
					sf_cut3->cd();
					if(cart[2] !=0){
						sf_cut3_w->cd();
					}
				break;
				case 3:
					sf_cut4->cd();
					if(cart[2] !=0){
						sf_cut4_w->cd();
					}
				break;
				case 4:
					sf_cut5->cd();
					if(cart[2] !=0){
						sf_cut5_w->cd();
					}
				break;
				case 5:
					sf_cut6->cd();
					if(cart[2] !=0){
						sf_cut6_w->cd();
					}
				break;
				case 6:
					sf_cut7->cd();
					if(cart[2] !=0){
						sf_cut7_w->cd();
					}
				break;
				case 7:
					sf_cut8->cd();
					if(cart[2] !=0){
						sf_cut8_w->cd();
					}
				break;
				case 8:
					sf_cut9->cd();
					if(cart[2] !=0){
						sf_cut9_w->cd();
					}
				break;
				case 9:
					sf_cut10->cd();
					if(cart[2] !=0){
						sf_cut10_w->cd();
					}
				break;
				case 10:
					sf_cut11->cd();
					if(cart[2] !=0){
						sf_cut11_w->cd();
					}
				break;
			}
		}
		SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetXTitle("Momentum (GeV)");
		SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetYTitle("Sampling Fraction");
		SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetOption("COLZ");
		SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->Write();
	}

}


//Delta T Cuts
void Histogram::DT_Make(){
	char hname[100];
	std::vector<long> space_dims(5);
	space_dims[0] = 3;  //species
	space_dims[1] = 7;  //Cuts
	space_dims[2] = 30; //W Binning
	space_dims[3] = 7; //Sector
	space_dims[4] = 6; //topology

	float bot,top;

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		if((cart[1] == 6 && cart[4]!=0) || (cart[1]!=6 && cart[4] ==0)){
			if(cart[2] == 0){
				sprintf(hname,"%s_DeltaT_%s_%s_W:ALL_%s",species[cart[0]+1],hid_cut[cart[1]],sec_list[cart[3]],topologies[cart[4]]);
			}else if(cart[3]==0 & cart[4] ==0 && cart[1]!=6){
				top = Wbin_start + cart[2]*Wbin_res;
				bot = top - Wbin_res;
				sprintf(hname,"%s_DeltaT_%s_%s_W:%f-%f_%s",species[cart[0]+1],hid_cut[cart[1]],sec_list[cart[3]],bot,top,topologies[cart[4]]);
			}
		}
		DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]] = std::make_shared<TH2D>(hname,hname, DTxres, DTxmin, DTxmax, DTyres, DTymin, DTymax);
	}
}

void Histogram::DT_Fill(int top, int part, float p, float d, float t, float d0, float t0, int cut, float W_, int sec){
	float dt = physics::delta_t(part, p, d, t, d0, t0);
	DT_hist[part][cut][Histogram::W_binning(W_)][sec][top]->Fill(p,dt);
	DT_hist[part][cut][0][sec][top]->Fill(p,dt);
	DT_hist[part][cut][Histogram::W_binning(W_)][0][top]->Fill(p,dt);
	DT_hist[part][cut][0][0][top]->Fill(p,dt);
}
void Histogram::DT_Write(){
	char dir_name[100]; 
	TDirectory* DT_plot = RootOutputFile->mkdir("DT_plots");
	TDirectory* par_dt[3][8][3][8][7][3];
	DT_plot->cd(); 
	for(int i = 1; i<4; i++){
		sprintf(dir_name,"%s_DT_plots",species[i]);
		par_dt[i][0][0][0][0][0]= DT_plot->mkdir(dir_name);
		for(int j = 1; j < 8; j++){//cut
			sprintf(dir_name,"%s_DT_%s_%s_%s_%s",species[i],hid_cut[j-1]);
			par_dt[i][j][0][0][0][0] = par_dt[i][0][0][0][0][0]->mkdir(dir_name);
			for( int k = 1; k < 3; k++){//W Dep 
				sprintf(dir_name,"%s_DT_%s_%s",species[i],hid_cut[j-1],W_dep_list[k-1]);
				par_dt[i][j][k][0][0][0] = par_dt[i][j][0][0][0][0]->mkdir(dir_name);
				for( int l = 1; l < 8; l++){//Sector
					sprintf(dir_name,"%s_DT_%s_%s_%s",species[i],hid_cut[j-1],W_dep_list[k-1],sec_list[l-1]);
					par_dt[i][j][k][l][0][0] = par_dt[i][j][k][0][0][0]->mkdir(dir_name);
					for( int m = 1; m <7; m++){//Topology
						sprintf(dir_name,"%s_DT_%s_%s_%s_%s",species[i],hid_cut[j-1],W_dep_list[k-1],sec_list[l-1],topologies[m-1]);
						par_dt[i][j][k][l][m][0] = par_dt[i][j][k][l][0][0]->mkdir(dir_name);
					}
				}
			}
			for(int k = 1; k < 8; k++){ //Sector 
				sprintf(dir_name,"%s_DT_%s_%s",species[i],hid_cut[j-1],sec_list[k-1]);
				par_dt[i][j][0][k][0][1] = par_dt[i][j][0][0][0][0]->mkdir(dir_name);
				for(int l = 1; l <7; l++){ //topology
					sprintf(dir_name,"%s_DT_%s_%s_%s",species[i],hid_cut[j-1],sec_list[k-1],topologies[l-1]);
					par_dt[i][j][0][k][l][1] = par_dt[i][j][0][k][0][1]->mkdir(dir_name);
					for(int m = 1; m < 3; m++){//W dep
						sprintf(dir_name,"%s_DT_%s_%s",species[i],hid_cut[j-1],sec_list[k-1],topologies[l-1],W_dep_list[m-1]);
						par_dt[i][j][m][k][l][1] = par_dt[i][j][0][k][l][1]->mkdir(dir_name);
					}
				}
			}
			for(int k = 1; k < 7; k++){ //topology 
				sprintf(dir_name,"%s_DT_%s_%s",species[i],hid_cut[j-1],topologies[k-1]);
				par_dt[i][j][0][0][k][2] = par_dt[i][j][0][0][0][0]->mkdir(dir_name);
				for(int l = 1; l <8; l++){ //sector
					sprintf(dir_name,"%s_DT_%s_%s_%s",species[i],hid_cut[j-1],topologies[k-1],sec_list[l-1]);
					par_dt[i][j][0][l][k][2] = par_dt[i][j][0][k][0][2]->mkdir(dir_name);
					for(int m = 1; m < 3; m++){//W dep
						sprintf(dir_name,"%s_DT_%s_%s_%s",species[i],hid_cut[j-1],topologies[k-1],sec_list[l-1],W_dep_list[m-1]);
						par_dt[i][j][m][l][k][2] = par_dt[i][j][0][k][l][2]->mkdir(dir_name);
					}
				}
			}
		}
	}
	TDirectory* pro_dt = DT_plot->mkdir("Proton DT Plots");
	TDirectory* pip_dt = DT_plot->mkdir("Pi+ DT Plots");
	TDirectory* pim_dt = DT_plot->mkdir("Pi- DT Plots");
	TDirectory* pro_dt_w = pro_dt->mkdir("Proton DT Plots W Dep");
	TDirectory* pip_dt_w = pip_dt->mkdir("Pi+ DT Plots W Dep");
	TDirectory* pim_dt_w = pim_dt->mkdir("Pi- DT Plots W Dep");

	std::vector<long> space_dims(5);
	space_dims[0] = 3;  //species
	space_dims[1] = 7;  //Cuts
	space_dims[2] = 30; //W Binning
	space_dims[3] = 7; //Sector
	space_dims[4] = 6; //topology

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		par_dt[cart[0]][0][0][0][0][0]->cd();//Main folder for particles
		if(cart[2]==0 && cart[3] == 0 && ((cart[1]!=6 && cart[4] == 0) || (cart[1] ==6 && cart[4]!=0))){ //These are all W, All Sectors, Combined Topology
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Delta T (ns)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
		}
		par_dt[cart[0]][cart[1]+1][0][0][0][0]->cd();//Go into the particle cuts
		if(cart[2] == 0 && cart[3] ==0 && ((cart[1]==6 && cart[4] !=0) || (cart[1]!=6 && cart[4]==0))){
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Delta T (ns)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
		}

		// W Dependence
		par_dt[cart[0]][cart[1]][cart[2]+1][0][0][0]->cd();
		if(cart[3] == 0 && ((cart[1]!=6 &&cart[4] == 0) || (cart[1]==6 && cart[4] == 5))){//All sec & all/No tops
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Delta T (ns)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
		}
		par_dt[cart[0]][cart[1]][cart[2]+1][cart[3]+1][0][0]->cd(); 
		if((cart[1]!=6 &&cart[4] == 0) || (cart[1]==6 && cart[4]==5) ){//all or no tops
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Delta T (ns)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
		}
		if(cart[1]==6){//Make sure we're looking at event selected things
			par_dt[cart[0]][cart[1]][cart[2]+1][cart[3]+1][cart[4]+1][0]->cd();
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Delta T (ns)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
		}
		
		// Sector
		par_dt[cart[0]][cart[1]][0][cart[3]+1][0][1]->cd();
		if(cart[3] == 0 && ((cart[1]!=6 &&cart[4] == 0) || (cart[1]==6 && cart[4] == 5))){//All W & all/No tops
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Delta T (ns)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
		}
		if(cart[1] == 6){
			par_dt[cart[0]][cart[1]][0][cart[3]+1][cart[4]+1][1]->cd(); 
			if(cart[2]==0){//all W
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Delta T (ns)");
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
			}
			par_dt[cart[0]][cart[1]][cart[2]+1][cart[3]+1][cart[4]+1][1]->cd();
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Delta T (ns)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
		}

		// Topology
		par_dt[cart[0]][cart[1]][0][0][cart[4]+1][2]->cd();
		if(cart[1] == 6){
			if(cart[3] == 0 && cart[2]==0){//All W & all sec
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Delta T (ns)");
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
			}
			par_dt[cart[0]][cart[1]][0][cart[3]+1][cart[4]+1][2]->cd(); 
			if(cart[2]==0){//all W
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Delta T (ns)");
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
			}
			par_dt[cart[0]][cart[1]][cart[2]+1][cart[3]+1][cart[4]+1][2]->cd();
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Delta T (ns)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
		}
	}
}

//Min CC Cuts
void Histogram::CC_Make(){
	char hname[100];
	std::vector<long> space_dims(5);
	space_dims[0] = 6;  //Sector
	space_dims[1] = 18; //Segment
	space_dims[2] = 5;  //Cut
	space_dims[3] = 4;  //Side of detector
	space_dims[4] = 6;  //Topology
}
void Histogram::CC_Fill(int top, int sec, int segm, int nphe, int cut){

}
void Histogram::CC_Write(){

}
//Missing Mass Cuts
void Histogram::MM_Make(){
	char hname[100];
	std::vector<long> space_dims(4);
	space_dims[0] = 5;  //Topology
	space_dims[1] = 3;  //Cuts
	space_dims[2] = 30; //W Binning
	space_dims[3] = 6;  //Topology
}
void Histogram::MM_Fill(int top, float mm, int cut){

}
void Histogram::MM_Write(){

}

