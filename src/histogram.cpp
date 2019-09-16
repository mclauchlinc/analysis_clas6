#include "histogram.hpp"


Histogram::Histogram(const std::string& output_file){
	RootOutputFile = fun::Name_File(output_file);
	def = new TCanvas("def");
	Histogram::WQ2_Make();
	Histogram::Fid_Make();
	Histogram::SF_Make();
	Histogram::DT_Make();
	Histogram::CC_Make();
	//Histogram::MM_Make();
}

Histogram::~Histogram() { this->Write(); }

void Histogram::Write(){
	std::cout<< "Writing" <<std::endl;
	RootOutputFile->cd();
	std::cout<<"Writing WQ2 Plots"; 
	Histogram::WQ2_Write();
	std::cout <<"Done!" <<std::endl <<"Writing Fid Plots"; 
	Histogram::Fid_Write();
	std::cout <<"Done!" <<std::endl <<"Writing SF Plots"; 
	Histogram::SF_Write();
	std::cout <<"Done!" <<std::endl <<"Writing DT Plots"; 
	Histogram::DT_Write();
	std::cout<<"Done!";
	CC_Write();
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
    	WQ2_hist[cart[0]][cart[1]] = std::make_shared<TH2F>( hname, hname, WQxres, WQxmin, WQxmax, WQyres, WQymin, WQymax); // constants.h
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
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = std::make_shared<TH2F>(hname,hname, FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
						num_hist++;
					}
				}
				if(cart[4]==0 && (((cart[1] != (length-1) || cart[5]==0)) || (cart[1]== (length-1) && cart[5]!=0 ))) {//No momentum dependence and correct event stuff
					sprintf(hname,"%s_Fid_%s_%s_W:ALL_%s_%s",species[cart[1]],sec_list[cart[0]],par_cut,p_cut,topologies[cart[5]]);//constants.hpp
					Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = std::make_shared<TH2F>(hname,hname, FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
					num_hist++;
				}
			}else{//Specific W Bins
				if(cart[5]==0 && cart[4]==0 && cart[1] != (length-1) ){ //No event selection or momentum dependence
					if((cart[1]!=0 && (cart[1] == 0 || cart[1]==3)) || (cart[1]==0 && (cart[1]==0 || cart[1] == 7))){//Isolated cuts for particles
						wtop = Wbin_start + (cart[3]*Wbin_res);//constants.hpp
						wbot = wtop - Wbin_res;//constants.hpp
						sprintf(hname,"%s_Fid_%s_%s_W:%f-%f_%s_%s",species[cart[1]],sec_list[cart[0]],par_cut,wbot,wtop,p_cut,topologies[cart[5]]);//constants.hpp
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = std::make_shared<TH2F>(hname,hname, FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
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
	TDirectory* pro_fid_w = pro_fid->mkdir("Proton Fid W-Dependence");
	TDirectory* pro_fid_p = pro_fid->mkdir("Proton Fid P-Dependence");
	TDirectory* pip_fid_w = pip_fid->mkdir("Pi+ Fid W-Dependence");
	TDirectory* pip_fid_p = pip_fid->mkdir("Pi+ Fid P-Dependence");
	TDirectory* pim_fid_w = pim_fid->mkdir("Pi- Fid W-Dependence");
	TDirectory* pim_fid_p = pim_fid->mkdir("Pi- Fid P-Dependence");
	TDirectory* ele_fid_s = ele_fid->mkdir("Electron Fid Sec-Dependence");
	TDirectory* pro_fid_s = pro_fid->mkdir("Proton Fid Sec-Dependence");
	TDirectory* pip_fid_s = pip_fid->mkdir("Pi+ Fid Sec-Dependence");
	TDirectory* pim_fid_s = pim_fid->mkdir("Pi- Fid Sec-Dependence");

	
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
	bool s_dep = false; 

	while(cart.GetNextCombination()){
		p_dep = false;
		w_dep = false;
		s_dep = false;
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
		if( cart[0] == 0 || w_dep || p_dep){
			s_dep = false;
		} else{
			s_dep = true; 
		}
		switch(cart[1]){
			case 0:
				ele_fid->cd();
				if(p_dep){
					ele_fid_p->cd();
				}else if(w_dep){
					ele_fid_w->cd();
				}else if(s_dep){
					ele_fid_s->cd();
				}
			break;
			case 1:
				pro_fid->cd();
				if(p_dep){
					pro_fid_p->cd();
				}else if(w_dep){
					pro_fid_w->cd();
				}
				else if(s_dep){
					pro_fid_s->cd();
				}
			break;
			case 2:
				pip_fid->cd();
				if(p_dep){
					pip_fid_p->cd();
				}else if(w_dep){
					pip_fid_w->cd();
				}
				else if(s_dep){
					pip_fid_s->cd();
				}
			break;
			case 3:
				pim_fid->cd();
				if(p_dep){
					pim_fid_p->cd();
				}else if(w_dep){
					pim_fid_w->cd();
				}else if(s_dep){
					pim_fid_s->cd();
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
		if((cart[0] == 10 && cart[3] != 0) || (cart[0] != 10 && cart[3] == 0) ){//Topology only matters for event selection cut
			if(cart[1] == 0 ){ //All W 
				sprintf(hname,"SF_%s_%s_W:ALL_%s",eid_cut[cart[0]],sec_list[cart[2]],topologies[cart[3]]);
				SF_hist[cart[0]][cart[1]][cart[2]][cart[3]] = std::make_shared<TH2F>(hname,hname, SFxres, SFxmin, SFxmax, SFyres, SFymin, SFymax);
				//std::cout<<std::endl <<"Created plot: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3];

			}else{	//Specific W Bins
				top = Wbin_start + cart[1]*Wbin_res;
				bot = top - Wbin_res;
				sprintf(hname,"SF_%s_%s_W:%f-%f_%s",eid_cut[cart[0]],sec_list[cart[2]],bot,top,topologies[cart[3]]);
				SF_hist[cart[0]][cart[1]][cart[2]][cart[3]] = std::make_shared<TH2F>(hname,hname, SFxres, SFxmin, SFxmax, SFyres, SFymin, SFymax);
				//std::cout<<std::endl <<"Created plot: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3];
			}
		}
	}
}

void Histogram::SF_Fill(int top, float p, float en, int cut, float W_, int sec){
	SF_hist[cut][Histogram::W_binning(W_)][sec][top]->Fill(p,en/p);
	SF_hist[cut][Histogram::W_binning(W_)][0][top]->Fill(p,en/p);
}

void Histogram::SF_Write(){
	TDirectory* dir_SF = RootOutputFile->mkdir("SF Plots");
	dir_SF->cd();
	TDirectory* sf_dir[11];//Cut
	TDirectory* sf_dir_w[11];
	TDirectory* sf_dir_sec[11][7];
	TDirectory* sf_dir_top[11][5];
	//[8][6];// W binning, Sector, Topology
	char dir_name[100];
	std::cout<<std::endl;
	for(int cut = 0 ; cut < 11; cut++){
		sprintf(dir_name,"SF_%s",eid_cut[cut]);
		sf_dir[cut] = dir_SF->mkdir(dir_name);
		//std::cout<<"Made Directory: " <<cut <<" 0 0 0" <<std::endl; 
		sprintf(dir_name,"SF_%s_%s",eid_cut[cut],W_dep_list[1]);
		sf_dir_w[cut] = sf_dir[cut]->mkdir(dir_name);
		//std::cout<<"Made Directory: " <<cut <<" " <<1 <<" 0 0" <<std::endl;
		for(int sec = 1; sec < 8 ; sec++){
			sprintf(dir_name,"SF_%s_%s",eid_cut[cut],sec_list[sec-1]);
			sf_dir_sec[cut][sec] = sf_dir[cut]->mkdir(dir_name);
			//std::cout<<"Made Directory: " <<cut <<" " <<0 <<" " <<sec <<" 0" <<std::endl;
		}
		for(int top = 1; top < 6; top++){
			sprintf(dir_name,"SF_%s_%s",eid_cut[cut],topologies[top]);
			sf_dir_top[cut][top] = sf_dir[cut]->mkdir(dir_name);
			//std::cout<<"Made Directory: " <<cut <<" " <<0 <<" " <<0 <<" "<<top <<std::endl;
		}
	}


	std::vector<long> space_dims(4);
	space_dims[0] = 11; //Electron Cuts
	space_dims[1] = 30; //W binning
	space_dims[2] = 7;  //Sector
	space_dims[3] = 6;  //Topology

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		dir_SF->cd();
		//General Entry
		//std::cout<<"Curr Vals: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3]<<std::endl ;
		sf_dir[cart[0]]->cd();
		//std::cout<<"    Now In: " <<cart[0] /*<<" " <<"0" <<" " <<"0" <<" " <<"0"*/<<std::endl ; 
		//Want just cuts, so all Sectors, W, and use Combined Topology
		if(cart[2]==0 && cart[1]==0 && ((cart[0]!=10 && cart[3]==0)||(cart[0]==10 && cart[3]==5))){
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetXTitle("Momentum (GeV)");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetYTitle("Sampling Fraction");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetOption("COLZ");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->Write();
		}
		//W binning
		if(cart[2]==0 && cart[1]!=0 && ((cart[0]!=10 && cart[3]==0)||(cart[0]==10 && cart[3]==5))){
			//std::cout<<"    Try In: " <<cart[0] <<" " <<"1" <<" " <<"0" <<" " <<"0"<<std::endl ;
			sf_dir_w[cart[0]]->cd();
			//std::cout<<"    Now In: " <<cart[0] <<" " <<"1" <<" " <<"0" <<" " <<"0"<<std::endl ; 
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetXTitle("Momentum (GeV)");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetYTitle("Sampling Fraction");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetOption("COLZ");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->Write();
		}
		//Sector Binning
		if(cart[2]!=0 && cart[1]==0 && ((cart[0]!=10 && cart[3]==0)||(cart[0]==10 && cart[3]==5))){
			//std::cout<<"    Try In: " <<cart[0] <<" " <<"0" <<" " <<cart[2]+1 <<" " <<"0"<<std::endl ; 
			sf_dir_sec[cart[0]][cart[2]]->cd();
			//std::cout<<"    Now In: " <<cart[0] <<" " <<"0" <<" " <<cart[2]+1 <<" " <<"0"<<std::endl ; 
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetXTitle("Momentum (GeV)");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetYTitle("Sampling Fraction");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetOption("COLZ");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->Write();
		}
		//Topology Binning
		if(cart[2]==0 && cart[1]==0 && (cart[0]==10 && cart[3]!=0)){
			//std::cout<<"    Now In: " <<cart[0] <<" " <<"0" <<" " <<"0" <<" " <<cart[3]<<std::endl ;
			sf_dir_top[cart[0]][cart[3]]->cd();
			//std::cout<<"    Now In: " <<cart[0] <<" " <<"0" <<" " <<"0" <<" " <<cart[3]<<std::endl ; 
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetXTitle("Momentum (GeV)");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetYTitle("Sampling Fraction");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->SetOption("COLZ");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]]->Write();
		}
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
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]] = std::make_shared<TH2F>(hname,hname, DTxres, DTxmin, DTxmax, DTyres, DTymin, DTymax);
			}else if(cart[3]==0 & cart[4] == 0 && cart[1]!=6){//Looking at specific Cuts on fiducial and pre cut regimes 
				top = Wbin_start + cart[2]*Wbin_res;
				bot = top - Wbin_res;
				sprintf(hname,"%s_DeltaT_%s_%s_W:%f-%f_%s",species[cart[0]+1],hid_cut[cart[1]],sec_list[cart[3]],bot,top,topologies[cart[4]]);
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]] = std::make_shared<TH2F>(hname,hname, DTxres, DTxmin, DTxmax, DTyres, DTymin, DTymax);
			}
		}
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
	TDirectory* par_dt[3][8][2][8][6];
	//std::cout<<"Did I get here?" <<std::endl; 
	DT_plot->cd(); 
	for(int i = 0; i<3; i++){
		sprintf(dir_name,"%s_DT_plots",species[i+1]);
		par_dt[i][0][0][0][0]= DT_plot->mkdir(dir_name);
		//std::cout<<"Made pointer:" <<i <<" 0 0 0 0" <<std::endl;
		for(int j = 1; j < 8; j++){//cut
			sprintf(dir_name,"%s_DT_%s",species[i+1],hid_cut[j-1]);
			par_dt[i][j][0][0][0] = par_dt[i][0][0][0][0]->mkdir(dir_name);
			//std::cout<<"    Made Pointers" <<i <<j <<" 0 0 0"<<std::endl;
			//W Dependence
			sprintf(dir_name,"%s_DT_%s_%s",species[i+1],hid_cut[j-1],W_dep_list[1]);
			par_dt[i][j][1][0][0] = par_dt[i][j][0][0][0]->mkdir(dir_name);
			//std::cout<<"    Made Pointers" <<i <<j <<" 1 0 0"<<std::endl;
			for(int k = 1; k < 8; k++){ //Sector 
				sprintf(dir_name,"%s_DT_%s_%s",species[i+1],hid_cut[j-1],sec_list[k-1]);
				par_dt[i][j][0][k][0] = par_dt[i][j][0][0][0]->mkdir(dir_name);
				//std::cout<<"    Made Pointers" <<i <<j <<" 0 " <<k <<" 0"<<std::endl;
				//std::cout<<"Sector Pointer" <<std::endl;
				}
			if(j == 7){//Event cut
				for(int k = 1; k < 6; k++){ //topology 
					sprintf(dir_name,"%s_DT_%s_%s",species[i+1],hid_cut[j-1],topologies[k]);
					par_dt[i][j][0][0][k] = par_dt[i][j][0][0][0]->mkdir(dir_name);
					//std::cout<<"    Made Pointers" <<i <<j <<" 0 0 " <<k<<std::endl;
					//std::cout<<"Topology Pointer" <<std::endl;

				}
			}
		}
	}

	std::vector<long> space_dims(5);
	space_dims[0] = 3;  //species
	space_dims[1] = 7;  //Cuts
	space_dims[2] = 30; //W Binning
	space_dims[3] = 7; //Sector
	space_dims[4] = 6; //topology

	CartesianGenerator cart(space_dims);
	std::cout<<std::endl;
	while(cart.GetNextCombination()){
		par_dt[cart[0]][0][0][0][0]->cd();//Main folder for particles
		//std::cout<<std::endl <<"Current Vals: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4]; 
		//std::cout<<" Now Writing in " <<cart[0] <<" 0 0 0 0"<<std::endl;
		par_dt[cart[0]][cart[1]+1][0][0][0]->cd();//Get into those CUTS
		//std::cout <<"      Now Writing in " <<cart[0] <<" " <<cart[1]+1 <<" 0 0 0"<<std::endl;
		//All W, Sectors, and Combined Topology for Event selection, but still by Cut
		if(cart[2] ==0 && cart[3] == 0 && ((cart[1]!=6 && cart[4]==0)||(cart[1]==6 && cart[4]==5))){
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Delta T (ns)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();			
		}
		
		//For W Range, but all sectors, combine topology for event selection, all sectors
		if(cart[2]!=0 && cart[3] == 0 && ((cart[1]!=6 && cart[4]==0))){//||(cart[1]==6 && cart[4]==5))){ //Issue when trying to look at event selection for this. Not sure why, but getting rid of it solved it *shrug* 9/12/19
			//std::cout <<"          trying to Write in " <<cart[0] <<" " <<cart[1]+1 <<" 1 0 0"<<std::endl;
			par_dt[cart[0]][cart[1]+1][0][0][0]->cd();//Get into those CUTS
			par_dt[cart[0]][cart[1]+1][1][0][0]->cd();
			//std::cout <<"   We are writing" <<std::endl;
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Delta T (ns)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();			
		}
		//For Sector Range, but all W, combine topology for event selection
		if(cart[2] ==0 && ((cart[1]!=6 && cart[4]==0)||(cart[1]==6 && cart[4]==5))){
			//std::cout <<"          Trying to Write in " <<cart[0] <<" " <<cart[1]+1 <<" 0 " <<cart[3]+1 <<" 0"<<std::endl;
			par_dt[cart[0]][cart[1]+1][0][0][0]->cd();//Get into those CUTS
			par_dt[cart[0]][cart[1]+1][0][cart[3]+1][0]->cd();
			//std::cout <<"  We are writing "<<std::endl;
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Delta T (ns)");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();			
		}
		//For topology Range, but all W, all sector
		if(cart[2] ==0 && cart[3] == 0 && cart[1] == 6 && cart[4]!=0){
			//std::cout <<"          Trying to Write in " <<cart[0] <<" " <<cart[1]+1 <<" 0 0 " <<cart[4]+1<<std::endl;
			par_dt[cart[0]][cart[1]+1][0][0][0]->cd();//Get into those CUTS
			par_dt[cart[0]][cart[1]+1][0][0][cart[4]]->cd();
			//std::cout <<"  we are writing"<<std::endl;
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
	space_dims[2] = 11;  //Cut
	space_dims[3] = 4;  //Side of detector
	space_dims[4] = 6;  //Topology

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		if((cart[2] == 10 && cart[4]!=0) || (cart[2]!=10 && cart[4] ==0)){
				sprintf(hname,"MinCC_%s_%s_Seg%d_%s_%s",eid_cut[cart[2]],sec_list[cart[0]+1],cart[1]+1,CC_det_side[cart[3]],topologies[cart[4]]);
				CC_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]] = std::make_shared<TH1F>(hname,hname, MinCCres, MinCCmin, MinCCmax);
		}
	}
}

void Histogram::CC_Fill(int top, int sec, int segm, int nphe, int cut){
	CC_hist[sec][segm][cut][detect::cc_lrc(segm)][top]->Fill(nphe);

}
void Histogram::CC_Write(){
	char dir_name[100]; 
	TDirectory* CC_plot = RootOutputFile->mkdir("CC_plots");
	TDirectory* par_cc[6][19][12][5][6];
	for(int i = 0; i< 6; i++){//Sectors
		sprintf(dir_name,"CC_%s",sec_list[i+1]);
		par_cc[i][0][0][0][0] = CC_plot->mkdir(dir_name);
		for(int j = 1; j < 19 ; j++){//Segments
			sprintf(dir_name,"CC_%s_Segm%d",sec_list[i+1],j);
			par_cc[i][j][0][0][0] = par_cc[i][0][0][0][0]->mkdir(dir_name);
			for(int k = 1; k < 5; k++){//Part of CC hit
				sprintf(dir_name,"CC_%s_Segm%d_%s",sec_list[i+1],j,CC_det_side[k-1]);
				par_cc[i][j][0][k][0] = par_cc[i][j][0][0][0]->mkdir(dir_name);
				for(int l = 1; l<12; l++){//EID Cuts
					sprintf(dir_name,"CC_%s_Segm%d_%s_%s",sec_list[i+1],j,CC_det_side[k-1],eid_cut[l-1]);
					par_cc[i][j][l][k][0] = par_cc[i][j][0][k][0]->mkdir(dir_name);
					if(l == 11){
						for(int m = 1; m<6; m++){//Topologies
							sprintf(dir_name,"CC_%s_Segm%d_%s_%s_%s",sec_list[i+1],j,CC_det_side[k-1],eid_cut[l-1],topologies[m]);
							par_cc[i][j][l][k][m] = par_cc[i][j][l][k][0]->mkdir(dir_name);
						}
					}
				}
			}
		}
	}

	CC_plot->cd();
	for(int i = 0; i< 6; i++){//Sectors
		//std::cout<<" Trying to Enter: " <<i <<std::endl;
		par_cc[i][0][0][0][0]->cd();
		//std::cout<<" Did do it Enter: " <<i <<std::endl;
		for(int j = 1; j < 19 ; j++){//Segments
			//std::cout<<" Trying to Enter: " <<i <<" " <<j <<std::endl;
			par_cc[i][j][0][0][0]->cd();
			//std::cout<<" Did do it Enter: " <<i <<" " <<j <<std::endl;
			for(int k = 1; k < 5; k++){//Part of CC hit
				//std::cout<<" Trying to Enter: " <<i <<" " <<j <<" 0 " <<k <<std::endl;
				par_cc[i][j][0][k][0]->cd();
				//std::cout<<" Did do it Enter: " <<i <<" " <<j <<" 0 " <<k <<std::endl;
				for(int l = 1; l<12; l++){//EID Cuts
					//std::cout<<" Trying to Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<std::endl;
					par_cc[i][j][l][k][0]->cd();
					//std::cout<<" Did do it Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<std::endl;
					for(int m = 0; m<6; m++){
						if(l == 11 && m!=0){//Event selected
							//std::cout<<" Trying to Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<" " <<m <<std::endl;
							par_cc[i][j][l][k][m]->cd();
							//std::cout<<" Did do it Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<" " <<m <<std::endl;
							CC_hist[i][j-1][l-1][k-1][m]->SetXTitle("num photoelectrons");
							CC_hist[i][j-1][l-1][k-1][m]->SetYTitle("Counts");
							CC_hist[i][j-1][l-1][k-1][m]->Write();
						}else if( l != 11 && m==0){//Not event selected
							//std::cout<<" Trying to Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<" " <<m <<std::endl;
							par_cc[i][j][l][k][0]->cd();
							//std::cout<<" Did do it Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<" " <<m <<std::endl;
							CC_hist[i][j-1][l-1][k-1][m]->SetXTitle("num photoelectrons");
							CC_hist[i][j-1][l-1][k-1][m]->SetYTitle("Counts");
							CC_hist[i][j-1][l-1][k-1][m]->Write();
						}
					}
				}
			}
		}
	}

	

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

