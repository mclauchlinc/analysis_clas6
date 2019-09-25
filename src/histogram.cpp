#include "histogram.hpp"


Histogram::Histogram(const std::string& output_file){
	RootOutputFile = fun::Name_File(output_file);
	def = new TCanvas("def");
	Histogram::WQ2_Make();
	Histogram::Fid_Make();
	Histogram::SF_Make();
	Histogram::DT_Make();
	Histogram::CC_Make();
	Histogram::MM_Make();
}

Histogram::~Histogram() { this->Write(); }

void Histogram::Write(){
	std::cout<< "Writing Plots" <<std::endl;
	RootOutputFile->cd();
	std::cout<<"Writing WQ2 Plots"; 
	Histogram::WQ2_Write();
	std::cout <<" Done!" <<std::endl <<"Writing Fid Plots"; 
	Histogram::Fid_Write();
	std::cout <<" Done!" <<std::endl <<"Writing SF Plots"; 
	Histogram::SF_Write();
	std::cout <<" Done!" <<std::endl <<"Writing DT Plots"; 
	Histogram::DT_Write();
	std::cout<<" Done!" <<std::endl <<"Writing CC Plots";
	CC_Write();
	std::cout<<" Done!" <<std::endl <<"Writing MM Plots";
	MM_Write();
	std::cout<<" Done!" <<std::endl <<"Closing Event Selection RootFile";
	RootOutputFile->Close();
	std::cout<<" Done!" <<std::endl;
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
		if(cart[0]==10 && cart[1]!=0){
			sprintf(hname,"W_Q2_%s_%s",eid_cut[cart[0]],topologies[cart[1]]); //constants.h and otherwise writing the specific cut to the right plot
	    	WQ2_hist[cart[0]][cart[1]] = std::make_shared<TH2F>( hname, hname, WQxres, WQxmin, WQxmax, WQyres, WQymin, WQymax); // constants.h
	    	WQ2_made_hist[cart[0]][cart[1]]=true;
		}else if( cart[0]!= 10 && cart[1]==0){
			sprintf(hname,"W_Q2_%s_%s",eid_cut[cart[0]],topologies[cart[1]]); //constants.h and otherwise writing the specific cut to the right plot
	    	WQ2_hist[cart[0]][cart[1]] = std::make_shared<TH2F>( hname, hname, WQxres, WQxmin, WQxmax, WQyres, WQymin, WQymax); // constants.h	
	    	WQ2_made_hist[cart[0]][cart[1]]=true;		
	    }
	}
}

void Histogram::WQ2_Fill(int top, int cut, float W_, float Q2_){
	WQ2_hist[cut][top]->Fill(W_,Q2_);
}

void Histogram::WQ2_Write(){
	char dir_name[100];
	TDirectory* dir_WQ2 = RootOutputFile->mkdir("W vs. Q2");
	dir_WQ2->cd();
	TDirectory* dir_WQ2_sub[6];
	for(int top = 0; top < 6; top ++){
		sprintf(dir_name,"W Q2 %s",topologies[top]);
		dir_WQ2_sub[top] = dir_WQ2->mkdir(dir_name);
		dir_WQ2_sub[top]->cd();
		for(int cut = 0; cut < 11; cut++){
			if((top==0 && cut !=10) || (top!=0 && cut==10)){
				if(WQ2_made_hist[cut][top]){
					WQ2_hist[cut][top]->SetXTitle("W (GeV)");
					WQ2_hist[cut][top]->SetYTitle("Q^{2} (GeV^{2}");
					WQ2_hist[cut][top]->SetOption("Colz");
					WQ2_hist[cut][top]->Write();
				}else{
					std::cout<<std::endl <<"WQ2 would have segfaulted: " <<cut <<" " <<top;
				}
				
			}
		}
		dir_WQ2_sub[top]->Close();
	}
	dir_WQ2->Close();
}

//Fiducial Cuts
void Histogram::Fid_Make(){
	char hname[100];
	char par_cut[100]; 
	float wtop, wbot,ptop,pbot; 
	std::vector<long> space_dims(7);
	space_dims[0] = 7; //Sector
	space_dims[1] = 4; //species
	space_dims[2] = 11;//cuts
	space_dims[3] = 30;//W Binning
	space_dims[4] = 12;//Momentum binning
	space_dims[5] = 6; //Topology
	space_dims[6] = 2; //Cut vs. anti-cut

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
			length = 11; 
		}else{
			//fid_cuts = hid_cut; 
			length = 7; 
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
				sprintf(par_cut,"%s_%s",eid_cut[cart[2]],cut_ver[cart[6]]);
			}else{
				sprintf(par_cut,"%s_%s",hid_cut[cart[2]],cut_ver[cart[6]]);
			}
			if(cart[3] == 0){//No W dependence
				if((cart[5] == 0) && (cart[2] != (length-1)) && (cart[4]!=0)){//Specific Momentum Bins and no topolgy dependence 
					//if((cart[1]!=0 && (cart[2] == 0 || cart[2]==3)) || (cart[1]==0 && (cart[2]==0 || cart[2] == 7))){//Isolated cuts for particles
						sprintf(hname,"%s_Fid_%s_%s_W:ALL_%s_%s",species[cart[1]],sec_list[cart[0]],par_cut,p_cut,topologies[cart[5]]);//constants.hpp
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]] = std::make_shared<TH2F>(hname,hname, FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
						num_hist++;
						Fid_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]=true;
					//}
				}
				if(cart[4]==0) {//No momentum dependence
					if(cart[2]!= (length-1)){//No event selection 
						if(cart[5]==0){
							sprintf(hname,"%s_Fid_%s_%s_W:ALL_%s_%s",species[cart[1]],sec_list[cart[0]],par_cut,p_cut,topologies[cart[5]]);//constants.hpp
							Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]] = std::make_shared<TH2F>(hname,hname, FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
							Fid_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]=true;
							num_hist++;
						}
					}else{//Event selection 
						if(cart[5]!=0){
							sprintf(hname,"%s_Fid_%s_%s_W:ALL_%s_%s",species[cart[1]],sec_list[cart[0]],par_cut,p_cut,topologies[cart[5]]);//constants.hpp
							Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]] = std::make_shared<TH2F>(hname,hname, FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
							Fid_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]=true;
							num_hist++;
						}
					}
				}
			}else{//Specific W Bins
				if(cart[5]==0 && cart[4]==0 && cart[1] != (length-1) ){ //No event selection or momentum dependence
					//if((cart[1]!=0 && (cart[2] == 0 || cart[2]==3)) || (cart[1]==0 && (cart[2]==0 || cart[2] == 7))){//Isolated cuts for particles
						wtop = Wbin_start + (cart[3]*Wbin_res);//constants.hpp
						wbot = wtop - Wbin_res;//constants.hpp
						sprintf(hname,"%s_Fid_%s_%s_W:%f-%f_%s_%s",species[cart[1]],sec_list[cart[0]],par_cut,wbot,wtop,p_cut,topologies[cart[5]]);//constants.hpp
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]] = std::make_shared<TH2F>(hname,hname, FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
						Fid_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]=true;
						num_hist++;
					//}
				}
			}
		}
	}
}

void Histogram::Fid_Fill(int top, float theta, float phi, int part, int cut, int cutvanti, float W_, float p_){
	//std::cout<< std::endl <<"We are inside the Fid Fill function ";
	//std::cout<< std::endl <<"top: " <<top <<" theta: " <<theta <<" phi: " <<phi <<" part: " <<part <<" cut: " <<cut <<" cva: " <<cutvanti <<" W: " <<W_ <<" p: " <<p_ <<std::endl; 
	//std::cout<< std::endl <<"filling: " <<physics::get_sector(phi) <<" " <<part <<" " <<cut <<" " <<0 <<" " <<Histogram::p_binning(p_) <<" " <<top <<" " <<cutvanti <<std::endl;

	float phic = physics::phi_center(phi);
	//std::cout<<"Filling p binning" <<std::endl;
	if(p_binning(p_) != 0 && top==0){ //All W with p binning
		Fid_fill_hist[physics::get_sector(phi)][part][cut][0][Histogram::p_binning(p_)][0][cutvanti]=true;
		if(Fid_made_hist[physics::get_sector(phi)][part][cut][0][Histogram::p_binning(p_)][0][cutvanti]){
			Fid_hist[0][part][cut][0][Histogram::p_binning(p_)][0][cutvanti]->Fill(phic,theta);//All Sectors
			Fid_hist[physics::get_sector(phi)][part][cut][0][Histogram::p_binning(p_)][0][cutvanti]->Fill(phic,theta);//Individual Sectors
			//std::cout<<"We did it! ";
		}else{
			std::cout<<"would have segfaulted ";
			std::cout<<std::endl <<"For Fid Plot: " <<physics::get_sector(phi) <<" " <<part <<" " <<cut <<" " <<0 <<" " <<Histogram::p_binning(p_) <<" " <<top <<" " <<cutvanti <<std::endl;
		}
	}
	//std::cout<<"Filling W binning" <<std::endl;
	if(W_binning(W_) != 0 && top==0){ //All p with W binning
		Fid_fill_hist[physics::get_sector(phi)][part][cut][Histogram::W_binning(W_)][0][0][cutvanti]=true;
		if(Fid_made_hist[physics::get_sector(phi)][part][cut][Histogram::W_binning(W_)][0][0][cutvanti]){
			Fid_hist[0][part][cut][Histogram::W_binning(W_)][0][0][cutvanti]->Fill(phic,theta);//All sectors
			Fid_hist[physics::get_sector(phi)][part][cut][Histogram::W_binning(W_)][0][0][cutvanti]->Fill(phic,theta);//Individual Sectors
			//std::cout<<"We did it! ";
		}else{
			std::cout<<"would have segfaulted ";
			std::cout<<std::endl <<"For Fid Plot: " <<physics::get_sector(phi) <<" " <<part <<" " <<cut <<" " <<Histogram::W_binning(W_) <<" " <<0 <<" " <<top <<" " <<cutvanti <<std::endl;
		}
		
	}
	//std::cout<<"Filling all p and w" <<std::endl;
	//All W and All Momentum
	Fid_fill_hist[physics::get_sector(phi)][part][cut][0][0][top][cutvanti]=true;
	if(Fid_made_hist[physics::get_sector(phi)][part][cut][0][0][top][cutvanti]){
		Fid_hist[0][part][cut][0][0][top][cutvanti]->Fill(phic,theta);//All sectors, all W, all P
		Fid_hist[physics::get_sector(phi)][part][cut][0][0][top][cutvanti]->Fill(phic,theta);//Individual sectors all W, all p
		//std::cout<<"We did it! ";
	}else{
		std::cout<<"would have segfaulted ";
		std::cout<<std::endl <<"For Fid Plot: " <<physics::get_sector(phi) <<" " <<part <<" " <<cut <<" " <<0 <<" " <<0 <<" " <<top <<" " <<cutvanti <<std::endl;
	}
	if(cut == 6 && (top == 3 || top == 4) && cutvanti ==1 ){
		//std::cout<<std::endl <<"For Fid plots sector is: " <<physics::get_sector(phi) <<" from phi: " <<phi <<std::endl;
	}
}

void Histogram::Fid_Write(){
	char dir_name[100];
	TDirectory* dir_Fid = RootOutputFile->mkdir("Fiducial Plots");
	dir_Fid->cd();
	TDirectory* par_fid[4];
	TDirectory* par_fid_p[4][10];
	TDirectory* par_fid_W[4][10];
	TDirectory* par_fid_cut[4][11];
	TDirectory* par_fid_top[4][5];
	TDirectory* par_fid_sec[4][11][8][7];
	int dir_len = 0; 

	for(int i = 0; i<4; i++){
		sprintf(dir_name,"%s Fid Plots",species[i]);
		par_fid[i]= dir_Fid->mkdir(dir_name);
		if(i ==0){
			dir_len = 11;
		}else{
			dir_len = 7; 
		}
		for(int j = 0; j<dir_len; j++){
			sprintf(dir_name,"%s Fid %s",species[i],par_cut[i][j]);
			par_fid_cut[i][j]=par_fid[i]->mkdir(dir_name);
			for(int l = 0; l<7; l++){
				sprintf(dir_name,"%s Fid %s %s",species[i],par_cut[i][j],sec_list[l]);
				par_fid_sec[i][j][0][l] = par_fid_cut[i][j]->mkdir(dir_name);
			}
			if(j!=(dir_len-1)){
				sprintf(dir_name,"%s Fid %s W Dep",species[i],par_cut[i][j]);
				par_fid_W[i][j] = par_fid_cut[i][j]->mkdir(dir_name);
				
				sprintf(dir_name,"%s Fid %s P Dep",species[i],par_cut[i][j]);
				par_fid_p[i][j] = par_fid_cut[i][j]->mkdir(dir_name);
				
				for(int l = 0; l<7; l++){
					sprintf(dir_name,"%s Fid %s P Dep %s",species[i],par_cut[i][j],sec_list[l]);
					par_fid_sec[i][j][2][l] = par_fid_p[i][j]->mkdir(dir_name);
					sprintf(dir_name,"%s Fid %s W Dep %s",species[i],par_cut[i][j],sec_list[l]);
					par_fid_sec[i][j][1][l] = par_fid_W[i][j]->mkdir(dir_name);
				}
			}else{
				for(int k = 0; k < 5; k++){
					sprintf(dir_name,"%s Fid Event: %s",species[i],topologies[k+1]);
					par_fid_top[i][k] = par_fid_cut[i][k]->mkdir(dir_name);
					for(int l = 0; l<7; l++){
						sprintf(dir_name,"%s Fid %s %s %s",species[i],par_cut[i][j],sec_list[l],topologies[k+1]);
						par_fid_sec[i][j][3+k][l] = par_fid_cut[i][j]->mkdir(dir_name);
					}
				}
			}	
		}
	}

	
	std::vector<long> space_dims(7);
	space_dims[0] = 7; //Sector
	space_dims[1] = 4; //species
	space_dims[2] = 11;//cuts
	space_dims[3] = 30;//W Binning
	space_dims[4] = 12;//Momentum binning
	space_dims[5] = 6; //Topology
	space_dims[6] = 2;//Cut vs anti cut

	CartesianGenerator cart(space_dims);

	//char * fid_cuts[];
	int length = 0;
	char p_cut[100]; 

	bool p_dep = false;
	bool w_dep = false; 
	bool s_dep = false; 

	while(cart.GetNextCombination()){
		par_fid[cart[1]]->cd();

		//Just making sure we aren't looking in places that don't exist 
		if(cart[1] == 0){ //Dealing with the different numbers of cuts for electrons vs. hadrons
			//fid_cuts = eid_cut;
			length = 10; //number of cuts outside of event selection 
		}else{
			//fid_cuts = hid_cut; 
			length = 6;  //number of cuts outside of event selection 
		}

		
		if(cart[2]<(length)){
			par_fid_cut[cart[1]][cart[2]]->cd();//Just the cuts
			if(cart[3] == 0){//All W bins
				if(cart[4]!=0){//Specific Momentum Bins
					par_fid_p[cart[1]][cart[2]]->cd();
					if(cart[5]==0){//Isolated cuts for particles
						if(Fid_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]){
							par_fid_sec[cart[1]][cart[2]][2][cart[0]]->cd();
							Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetXTitle("{phi} (degrees)");
							Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetYTitle("{theta} (degrees)");
							Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetOption("COLZ");
							Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->Write();
						}else{
							std::cout<<std::endl <<"Would have segfaulted Fid plot " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<" " <<cart[6] <<std::endl;
						}
					}
				}else if(cart[5]==0){//All P bins
					if(Fid_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]){
						par_fid_sec[cart[1]][cart[2]][0][cart[0]]->cd();
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetXTitle("{phi} (degrees)");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetYTitle("{theta} (degrees)");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetOption("COLZ");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->Write();
					}else{
						std::cout<<std::endl <<"Would have segfaulted Fid plot " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<" " <<cart[6] <<std::endl;
					}
				}
			}else{//Specific W Bins
				par_fid_W[cart[1]][cart[2]]->cd();
				if(cart[4]==0 ){ //No momentum dependence
					if(cart[5]==0){//No event selection 
						if(Fid_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]){
							par_fid_sec[cart[1]][cart[2]][1][cart[0]]->cd();
							Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetXTitle("{phi} (degrees)");
							Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetYTitle("{theta} (degrees)");
							Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetOption("COLZ");
							Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->Write();
						}else{
							std::cout<<std::endl <<"Would have segfaulted Fid plot " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<" " <<cart[6] <<std::endl;
						}					
					}
				}
			}
		}else if(cart[2]==(length)){
			par_fid_cut[cart[1]][cart[2]]->cd();//Just the cuts
			if(cart[5]!=0){//Event Selection topologies
				if(cart[4]==0 && cart[3]==0){
					if(Fid_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]){
						par_fid_sec[cart[1]][cart[2]][2+cart[5]][cart[0]]->cd();
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetXTitle("{phi} (degrees)");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetYTitle("{theta} (degrees)");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetOption("COLZ");
						Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->Write();
					}else{
						std::cout<<std::endl <<"Would have segfaulted Fid plot " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<" " <<cart[6] <<std::endl;
					}
				}
			}
		}
	}
}

//Sampling Fraction Cuts
void Histogram::SF_Make(){
	char hname[100];
	std::vector<long> space_dims(5);
	space_dims[0] = 11; //Electron Cuts
	space_dims[1] = 30; //W binning
	space_dims[2] = 7;  //Sector
	space_dims[3] = 6;  //Topology
	space_dims[4] = 2; //cut v anti

	float top,bot; 

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		if((cart[0] == 10 && cart[3] != 0) || (cart[0] != 10 && cart[3] == 0) ){//Topology only matters for event selection cut
			if(cart[1] == 0 ){ //All W 
				sprintf(hname,"SF_%s_%s_%s_W:ALL_%s",eid_cut[cart[0]],cut_ver[cart[4]],sec_list[cart[2]],topologies[cart[3]]);
				SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]] = std::make_shared<TH2F>(hname,hname, SFxres, SFxmin, SFxmax, SFyres, SFymin, SFymax);
				//std::cout<<std::endl <<"Created plot: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3];

			}else{	//Specific W Bins
				top = Wbin_start + cart[1]*Wbin_res;
				bot = top - Wbin_res;
				sprintf(hname,"SF_%s_%s_%s_W:%f-%f_%s",eid_cut[cart[0]],cut_ver[cart[4]],sec_list[cart[2]],bot,top,topologies[cart[3]]);
				SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]] = std::make_shared<TH2F>(hname,hname, SFxres, SFxmin, SFxmax, SFyres, SFymin, SFymax);
				//std::cout<<std::endl <<"Created plot: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3];
			}
		}
	}
}

void Histogram::SF_Fill(int top, float p, float en, int cut, int cva, float W_, int sec){
	SF_hist[cut][Histogram::W_binning(W_)][sec][top][cva]->Fill(p,en/p);
	SF_hist[cut][Histogram::W_binning(W_)][0][top][cva]->Fill(p,en/p);
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


	std::vector<long> space_dims(5);
	space_dims[0] = 11; //Electron Cuts
	space_dims[1] = 30; //W binning
	space_dims[2] = 7;  //Sector
	space_dims[3] = 6;  //Topology
	space_dims[4] = 2; //Cut v anti

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		dir_SF->cd();
		//General Entry
		//std::cout<<"Curr Vals: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3]<<std::endl ;
		sf_dir[cart[0]]->cd();
		//std::cout<<"    Now In: " <<cart[0] /*<<" " <<"0" <<" " <<"0" <<" " <<"0"*/<<std::endl ; 
		//Want just cuts, so all Sectors, W, and use Combined Topology
		if(cart[2]==0 && cart[1]==0 && ((cart[0]!=10 && cart[3]==0)||(cart[0]==10 && cart[3]==5))){
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Sampling Fraction");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
		}
		//W binning
		if(cart[2]==0 && cart[1]!=0 && ((cart[0]!=10 && cart[3]==0)||(cart[0]==10 && cart[3]==5))){
			//std::cout<<"    Try In: " <<cart[0] <<" " <<"1" <<" " <<"0" <<" " <<"0"<<std::endl ;
			sf_dir_w[cart[0]]->cd();
			//std::cout<<"    Now In: " <<cart[0] <<" " <<"1" <<" " <<"0" <<" " <<"0"<<std::endl ; 
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Sampling Fraction");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
		}
		//Sector Binning
		if(cart[2]!=0 && cart[1]==0 && ((cart[0]!=10 && cart[3]==0)||(cart[0]==10 && cart[3]==5))){
			//std::cout<<"    Try In: " <<cart[0] <<" " <<"0" <<" " <<cart[2]+1 <<" " <<"0"<<std::endl ; 
			sf_dir_sec[cart[0]][cart[2]]->cd();
			//std::cout<<"    Now In: " <<cart[0] <<" " <<"0" <<" " <<cart[2]+1 <<" " <<"0"<<std::endl ; 
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Sampling Fraction");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
		}
		//Topology Binning
		if(cart[2]==0 && cart[1]==0 && (cart[0]==10 && cart[3]!=0)){
			//std::cout<<"    Now In: " <<cart[0] <<" " <<"0" <<" " <<"0" <<" " <<cart[3]<<std::endl ;
			sf_dir_top[cart[0]][cart[3]]->cd();
			//std::cout<<"    Now In: " <<cart[0] <<" " <<"0" <<" " <<"0" <<" " <<cart[3]<<std::endl ; 
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Sampling Fraction");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
		}
	}
}


//Delta T Cuts
void Histogram::DT_Make(){
	char hname[100];
	std::vector<long> space_dims(6);
	space_dims[0] = 3;  //species
	space_dims[1] = 7;  //Cuts
	space_dims[2] = 30; //W Binning
	space_dims[3] = 7; //Sector
	space_dims[4] = 6; //topology
	space_dims[5] = 2; //cut v anti

	float bot,top;

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		if((cart[1] == 6 && cart[4]!=0) || (cart[1]!=6 && cart[4] ==0)){
			if(cart[2] == 0){
				sprintf(hname,"%s_DeltaT_%s_%s_%s_W:ALL_%s",species[cart[0]+1],hid_cut[cart[1]],cut_ver[cart[5]],sec_list[cart[3]],topologies[cart[4]]);
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = std::make_shared<TH2F>(hname,hname, DTxres, DTxmin, DTxmax, DTyres, DTymin, DTymax);
				DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = true; 	
			}else if(cart[3]==0 & cart[4] == 0 && cart[1]!=6){//Looking at specific Cuts on fiducial and pre cut regimes 
				top = Wbin_start + cart[2]*Wbin_res;
				bot = top - Wbin_res;
				sprintf(hname,"%s_DeltaT_%s_%s_%s_W:%f-%f_%s",species[cart[0]+1],hid_cut[cart[1]],cut_ver[cart[5]],sec_list[cart[3]],bot,top,topologies[cart[4]]);
				DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = std::make_shared<TH2F>(hname,hname, DTxres, DTxmin, DTxmax, DTyres, DTymin, DTymax);
				DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = true;
			}
		}
	}
}

void Histogram::DT_Fill(int top, int part, float p, float d, float t, float d0, float t0, int cut, int anti, float W_, int sec){
	float dt = physics::delta_t(part, p, d, t, d0, t0);
	if(sec < 7){
		if(cut!=6){//No event selection in the W variance. 
			if(DT_made_hist[part][cut][Histogram::W_binning(W_)][0][top][anti]){
				DT_hist[part][cut][Histogram::W_binning(W_)][0][top][anti]->Fill(p,dt);
			}else{
				std::cout<<"DT Would have segfaulted filling: " <<part <<" " <<cut <<" " <<Histogram::W_binning(W_) <<" " <<sec <<" " <<top <<" " <<anti <<std::endl;
			}
		}
		if(DT_made_hist[part][cut][0][sec][top][anti]){
			DT_hist[part][cut][0][sec][top][anti]->Fill(p,dt);
		}else{
			std::cout<<"DT Would have segfaulted filling: " <<part <<" " <<cut <<" " <<0 <<" " <<sec <<" " <<top <<" " <<anti <<std::endl;
		}
		if(DT_made_hist[part][cut][0][0][top][anti]){
			DT_hist[part][cut][0][0][top][anti]->Fill(p,dt);
		}else{
			std::cout<<"DT Would have segfaulted filling: " <<part <<" " <<cut <<" " <<0 <<" " <<0 <<" " <<top <<" " <<anti <<std::endl;
		}
	}else{
		std::cout<<"DT would have had weird nonexistant sectors " <<std::endl;
	}
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
		DT_dir_made[i][0][0][0][0]=true;
		//std::cout<<"Made pointer:" <<i <<" 0 0 0 0" <<std::endl;
		for(int j = 1; j < 8; j++){//cut
			sprintf(dir_name,"%s_DT_%s",species[i+1],hid_cut[j-1]);
			par_dt[i][j][0][0][0] = par_dt[i][0][0][0][0]->mkdir(dir_name);
			DT_dir_made[i][j][0][0][0]=true;
			//std::cout<<"    Made Pointers" <<i <<j <<" 0 0 0"<<std::endl;
			//W Dependence
			sprintf(dir_name,"%s_DT_%s_%s",species[i+1],hid_cut[j-1],W_dep_list[1]);
			par_dt[i][j][1][0][0] = par_dt[i][j][0][0][0]->mkdir(dir_name);
			DT_dir_made[i][j][1][0][0]=true;
			//std::cout<<"    Made Pointers" <<i <<j <<" 1 0 0"<<std::endl;
			for(int k = 1; k < 8; k++){ //Sector 
				sprintf(dir_name,"%s_DT_%s_%s",species[i+1],hid_cut[j-1],sec_list[k-1]);
				par_dt[i][j][0][k][0] = par_dt[i][j][0][0][0]->mkdir(dir_name);
				DT_dir_made[i][j][0][k][0]=true;
				//std::cout<<"    Made Pointers" <<i <<j <<" 0 " <<k <<" 0"<<std::endl;
				//std::cout<<"Sector Pointer" <<std::endl;
				}
			if(j == 7){//Event cut
				for(int l = 1; l < 6; l++){ //topology 
					sprintf(dir_name,"%s_DT_%s_%s",species[i+1],hid_cut[j-1],topologies[l]);
					par_dt[i][j][0][0][l] = par_dt[i][j][0][0][0]->mkdir(dir_name);
					DT_dir_made[i][j][0][0][l]=true;
					//std::cout<<"    Made Pointers" <<i <<j <<" 0 0 " <<k<<std::endl;
					//std::cout<<"Topology Pointer" <<std::endl;

				}
			}
		}
	}

	std::vector<long> space_dims(6);
	space_dims[0] = 3;  //species
	space_dims[1] = 7;  //Cuts
	space_dims[2] = 30; //W Binning
	space_dims[3] = 7; //Sector
	space_dims[4] = 6; //topology
	space_dims[5] = 2; //Cut v anti

	CartesianGenerator cart(space_dims);
	while(cart.GetNextCombination()){
		if(DT_dir_made[0][0][0][0][0]){
			par_dt[cart[0]][0][0][0][0]->cd();//Main folder for particles

			//std::cout<<std::endl <<"Current Vals: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4]; 
			//std::cout<<" Now Writing in " <<cart[0] <<" 0 0 0 0"<<std::endl;
			if(DT_dir_made[cart[0]][cart[1]+1][0][0][0]){
				par_dt[cart[0]][cart[1]+1][0][0][0]->cd();//Get into those CUTS
				//std::cout <<"      Now Writing in " <<cart[0] <<" " <<cart[1]+1 <<" 0 0 0"<<std::endl;
				//All W, Sectors, and Combined Topology for Event selection, but still by Cut
				if(cart[2] ==0 && cart[3] == 0 && ((cart[1]!=6 && cart[4]==0)||(cart[1]==6 && cart[4]==5))){
					if(DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]){
						DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetXTitle("Momentum (GeV)");
						DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetYTitle("Delta T (ns)");
						DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetOption("COLZ");
						DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->Write();	
					}else{
						std::cout<<"DT Would have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
					}		
				}
				
				//For W Range, but all sectors, combine topology for event selection, all sectors
				if(cart[2]!=0 && cart[3] == 0 && ((cart[1]!=6 && cart[4]==0))){//||(cart[1]==6 && cart[4]==5))){ //Issue when trying to look at event selection for this. Not sure why, but getting rid of it solved it *shrug* 9/12/19
					//std::cout <<"          trying to Write in " <<cart[0] <<" " <<cart[1]+1 <<" 1 0 0"<<std::endl;
					if(DT_dir_made[cart[0]][cart[1]+1][1][0][0]){
						par_dt[cart[0]][cart[1]+1][0][0][0]->cd();//Get into those CUTS
						par_dt[cart[0]][cart[1]+1][1][0][0]->cd();
						//std::cout <<"   We are writing" <<std::endl;
						if(DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]){
							DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetXTitle("Momentum (GeV)");
							DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetYTitle("Delta T (ns)");
							DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetOption("COLZ");
							DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->Write();	
						}else{
							std::cout<<"DT Would have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
						}	
					}else{
						std::cout<<"cart values: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
						std::cout<<"DTWould have filled nonexistant dir on: " <<cart[0] <<" " <<cart[1]+1 <<" " <<1 <<" " <<0 <<" " <<0 <<" " <<0 <<std::endl;
					}
							
				}
				//For Sector Range, but all W, combine topology for event selection
				if(cart[2] ==0 && ((cart[1]!=6 && cart[4]==0)||(cart[1]==6 && cart[4]==5))){
					//std::cout <<"          Trying to Write in " <<cart[0] <<" " <<cart[1]+1 <<" 0 " <<cart[3]+1 <<" 0"<<std::endl;
					if(DT_dir_made[cart[0]][cart[1]+1][0][cart[3]+1][0]){
						par_dt[cart[0]][cart[1]+1][0][0][0]->cd();//Get into those CUTS
						par_dt[cart[0]][cart[1]+1][0][cart[3]+1][0]->cd();
						//std::cout <<"  We are writing "<<std::endl;
						if(DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]){
							DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetXTitle("Momentum (GeV)");
							DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetYTitle("Delta T (ns)");
							DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetOption("COLZ");
							DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->Write();	
						}else{
							std::cout<<"DT Would have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
						}
					}else{
				std::cout<<"cart values: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
				std::cout<<"DTWould have filled nonexistant dir on: " <<cart[0] <<" " <<cart[1]+1 <<" " <<0 <<" " <<cart[3]+1 <<" " <<0 <<" " <<0 <<std::endl;
			}		
				}
				//For topology Range, but all W, all sector
				if(cart[2] ==0 && cart[3] == 0 && cart[1] == 6 && cart[4]!=0){
					//std::cout <<"          Trying to Write in " <<cart[0] <<" " <<cart[1]+1 <<" 0 0 " <<cart[4]+1<<std::endl;
					if(DT_dir_made[cart[0]][cart[1]+1][0][0][cart[4]]+1){
						par_dt[cart[0]][cart[1]+1][0][0][0]->cd();//Get into those CUTS
						par_dt[cart[0]][cart[1]+1][0][0][cart[4]+1]->cd();
						//std::cout <<"  we are writing"<<std::endl;
						if(DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]){
							DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetXTitle("Momentum (GeV)");
							DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetYTitle("Delta T (ns)");
							DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetOption("COLZ");
							DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->Write();	
						}else{
							std::cout<<"DT Would have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
						}
					}else{
						std::cout<<"cart values: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
						std::cout<<"DTWould have filled nonexistant dir on: " <<cart[0] <<" " <<cart[1]+1 <<" " <<1 <<" " <<0 <<" " <<0 <<" " <<0 <<std::endl;
					}		
				}
			}else{
				std::cout<<"cart values: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
				std::cout<<"DTWould have filled nonexistant dir on: " <<cart[0] <<" " <<cart[1]+1 <<" " <<0 <<" " <<0 <<" " <<0 <<" " <<cart[4]+1 <<std::endl;
			}
						
		}else{
				std::cout<<"cart values: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
				std::cout<<"DTWould have filled nonexistant dir on: " <<0 <<" " <<0 <<" " <<0 <<" " <<0 <<" " <<0 <<" " <<cart[4]+1 <<std::endl;
			}

	}
}

//Min CC Cuts
void Histogram::CC_Make(){
	char hname[100];
	std::vector<long> space_dims(6);
	space_dims[0] = 6;  //Sector
	space_dims[1] = 18; //Segment
	space_dims[2] = 11;  //Cut
	space_dims[3] = 4;  //Side of detector
	space_dims[4] = 6;  //Topology
	space_dims[5] = 2; //cut v anti

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		if((cart[2] == 10 && cart[4]!=0) || (cart[2]!=10 && cart[4] ==0)){//Making sure event selection lines up with topologies
				sprintf(hname,"MinCC_%s_%s_%s_Seg%d_%s_%s",eid_cut[cart[2]],cut_ver[cart[5]],sec_list[cart[0]+1],cart[1]+1,CC_det_side[cart[3]],topologies[cart[4]]);
				CC_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = std::make_shared<TH1F>(hname,hname, MinCCres, MinCCmin, MinCCmax);
				CC_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]=true;//std::cout<<std::endl<<"Making CC plot: ";
				//std::cout<<" Made CC plot: "<<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
		}
		if(!CC_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]){
			//std::cout<<"Didn't Make CC plot: "<<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
		}
	}
}

void Histogram::CC_Fill(int top, int sec, int segm, int nphe, int cut, int anti){
	//std::cout<<std::endl <<"top: " <<top <<" sec: " <<sec <<" segm: " <<segm <<" nphe: " <<nphe <<" cut: " <<cut <<" anti: " <<anti <<std::endl;
	//std::cout<<"filling plot";
	//std::cout<<std::endl <<sec <<" " <<detect::cc_segment(segm) <<" " <<cut <<" " <<detect::cc_lrc(segm) <<" " <<top <<" " <<anti <<std::endl;
	
	if((cut == 10 && top!=0) || (cut!=10 && top ==0)){
		//if(((sec-1) < 6) && ((detect::cc_segment(segm)) <18)&& (cut < 11) && (detect::cc_lrc(segm)<4) && (top < 6) && (anti < 2)){
		if(CC_made_hist[sec-1][detect::cc_segment(segm)][cut][detect::cc_lrc(segm)][top][anti]){
			//CC_fill_hist[sec-1][detect::cc_segment(segm)][cut][detect::cc_lrc(segm)][top][anti]=true;
			//if(CC_fill_hist[sec-1][detect::cc_segment(segm)][cut][detect::cc_lrc(segm)][top][anti] && ){
				CC_hist[sec-1][detect::cc_segment(segm)][cut][detect::cc_lrc(segm)][top][anti]->Fill(nphe);
				CC_hist[sec-1][detect::cc_segment(segm)][cut][3][top][anti]->Fill(nphe);
			
			//}
			
		}else{
				std::cout<<std::endl <<"Would have segfaulted in CC filling";
				std::cout<<std::endl <<"for CC Plot: " <<sec-1 <<" " <<detect::cc_segment(segm) <<" " <<cut <<" " <<detect::cc_lrc(segm) <<" " <<top <<" " <<anti <<std::endl;
		}
	}else{
		//std::cout<<std::endl <<"You're trying to fill an event thing that doesn't correspond to the correct cut";
	}

}
void Histogram::CC_Write(){
	char dir_name[100]; 
	TDirectory* CC_plot = RootOutputFile->mkdir("CC_plots");
	TDirectory* par_cc[6][12][5][6];
	for(int i = 0; i< 6; i++){//Sectors
		sprintf(dir_name,"CC_%s",sec_list[i+1]);
		par_cc[i][0][0][0] = CC_plot->mkdir(dir_name);
		//for(int j = 1; j < 19 ; j++){//Segments
			//sprintf(dir_name,"CC_%s_Segm%d",sec_list[i+1],j);
			//par_cc[i][j][0][0][0] = par_cc[i][0][0][0][0]->mkdir(dir_name);
			for(int k = 1; k < 5; k++){//Part of CC hit
				sprintf(dir_name,"CC_%s_%s",sec_list[i+1],CC_det_side[k-1]);
				par_cc[i][0][k][0] = par_cc[i][0][0][0]->mkdir(dir_name);
				for(int l = 1; l<12; l++){//EID Cuts
					sprintf(dir_name,"CC_%s_%s_%s",sec_list[i+1],CC_det_side[k-1],eid_cut[l-1]);
					par_cc[i][l][k][0] = par_cc[i][0][k][0]->mkdir(dir_name);
					if(l == 11){
						for(int m = 1; m<6; m++){//Topologies
							sprintf(dir_name,"CC_%s_%s_%s_%s",sec_list[i+1],CC_det_side[k-1],eid_cut[l-1],topologies[m]);
							par_cc[i][l][k][m] = par_cc[i][l][k][0]->mkdir(dir_name);
						}
					}
				}
			}
		//}
	}

	CC_plot->cd();
	for(int i = 0; i< 6; i++){//Sectors
		//std::cout<<" Trying to Enter: " <<i <<std::endl;
		par_cc[i][0][0][0]->cd();
		//std::cout<<" Did do it Enter: " <<i <<std::endl;
		for(int j = 0; j < 18 ; j++){//Segments
			//std::cout<<" Trying to Enter: " <<i <<" " <<j <<std::endl;
			//par_cc[i][0][0][0]->cd();
			//std::cout<<" Did do it Enter: " <<i <<" " <<j <<std::endl;
			for(int k = 1; k < 5; k++){//Part of CC hit
				//std::cout<<" Trying to Enter: " <<i <<" " <<j <<" 0 " <<k <<std::endl;
				par_cc[i][0][k][0]->cd();
				//std::cout<<" Did do it Enter: " <<i <<" " <<j <<" 0 " <<k <<std::endl;
				for(int l = 1; l<12; l++){//EID Cuts
					//std::cout<<" Trying to Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<std::endl;
					par_cc[i][l][k][0]->cd();
					//std::cout<<" Did do it Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<std::endl;
					for(int m = 0; m<6; m++){
						for(int n = 0; n<2; n++){//cut v anti-cut
							if(l == 11 && m!=0){//Event selected
								//std::cout<<" Trying to Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<" " <<m <<std::endl;
								par_cc[i][l][k][m]->cd();
								//std::cout<<" Did do it Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<" " <<m <<std::endl;
								CC_hist[i][j][l-1][k-1][m][n]->SetXTitle("num photoelectrons");
								CC_hist[i][j][l-1][k-1][m][n]->SetYTitle("Counts");
								CC_hist[i][j][l-1][k-1][m][n]->Write();
							}else if( l != 11 && m==0){//Not event selected
								//std::cout<<" Trying to Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<" " <<m <<std::endl;
								par_cc[i][l][k][m]->cd();
								//std::cout<<" Did do it Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<" " <<m <<std::endl;
								CC_hist[i][j][l-1][k-1][m][n]->SetXTitle("num photoelectrons");
								CC_hist[i][j][l-1][k-1][m][n]->SetYTitle("Counts");
								CC_hist[i][j][l-1][k-1][m][n]->Write();
							}
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
	std::vector<long> space_dims(3);
	space_dims[0] = 4;  //Topology
	space_dims[1] = 3;  //Cuts {pre, cut, anti}
	space_dims[2] = 2; //squared vs linear

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		sprintf(hname,"%s_MM_%s_%s",topologies[cart[0]+1],basic_cut[cart[1]],MM_sq[cart[2]]);
		MM_hist[cart[0]][cart[1]][cart[2]] = std::make_shared<TH1F>(hname,hname,MMxres,MMxmin,MMxmax);
	}
}
void Histogram::MM_Fill(int top, float mm, int cut, int square){
	MM_hist[top][cut][square]->Fill(mm);
}
void Histogram::MM_Write(){
	char dir_name[100];
	TDirectory * MM_plot = RootOutputFile->mkdir("MM plots");
	MM_plot->cd();
	TDirectory * MM_dir[4][3];//top, cut
	for(int i = 0; i < 4; i++ ){
		sprintf(dir_name,"MM %s",topologies[i+1]);
		MM_dir[i][0] = MM_plot->mkdir(dir_name);
		for(int j = 1; j <3; j++){
			sprintf(dir_name,"MM %s %s",topologies[i+1],MM_sq[j-1]);
			MM_dir[i][j] = MM_dir[i][0]->mkdir(dir_name);
		}
	}
	for(int i = 0; i < 4; i++ ){//Topology without "All"
		MM_dir[i][0]->cd();
		for(int j = 1; j <3; j++){//linear vs squared
			MM_dir[i][j]->cd();
			for( int k = 0 ; k<3; k++){//Cuts
				sprintf(dir_name,"MM %s %s %s",topologies[i+1],MM_sq[j-1],basic_cut[k]);
				MM_hist[i][k][j-1]->SetXTitle(dir_name);
				MM_hist[i][k][j-1]->SetYTitle("Events");
				MM_hist[i][k][j-1]->Write();
			}
		}
	}

}

