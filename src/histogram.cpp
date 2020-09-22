#include "histogram.hpp"


Histogram::Histogram(std::shared_ptr<Environment> _envi, const std::string& output_file){
	//RootOutputFile = fun::Name_File(output_file);
	//def = new TCanvas("def");
	Histogram::WQ2_Make(_envi);
	Histogram::Fid_Make(_envi);
	Histogram::SF_Make(_envi);
	Histogram::DT_Make(_envi);
	Histogram::CC_Make(_envi);
	Histogram::MM_Make(_envi);
	Histogram::Friend_Make(_envi);
	Histogram::Cross_Make(_envi);
	Histogram::XY_Make(_envi);
	Histogram::Fid_Det_Make(_envi);
}

//Histogram::~Histogram() { this->Write(); }

void Histogram::Write(const std::string& output_file, std::shared_ptr<Environment> _envi){
	RootOutputFile = fun::Name_File(output_file);
	def = new TCanvas("def");
	std::cout<< "Writing Plots" <<std::endl;
	RootOutputFile->cd();
	Histogram::WQ2_Write( _envi);
	Histogram::Fid_Write( _envi);
	Histogram::SF_Write( _envi);
	Histogram::DT_Write( _envi);
	Histogram::CC_Write( _envi);
	Histogram::MM_Write( _envi);
	Histogram::XY_Write(_envi);
	Histogram::Fid_Det_Write(_envi);
	Friend_Write(_envi);
	Histogram::Cross_Write(_envi);
	RootOutputFile->Close();
	std::cout<<"Histograms Done!" <<std::endl;
}


void Histogram::Print(const std::string& output_dir, std::shared_ptr<Environment> envi_){
	int _check_ = -1;
	
	//if(envi_->was_print()>0){
		std::string _curr_dir_ = fun::get_current_dir();
		std::string _out_dir_ = "$curr/$name/Plots";
		std::cout<<"Printing Plots \n 	Making Image File: ";
		//HistImageFile = fun::Name_Image_File(output_dir);
		std::cout<<"Making Plot Directory:";
		fun::replace(_out_dir_, "$curr", _curr_dir_);
		fun::replace(_out_dir_, "$name", output_dir);
		_check_ = fun::Make_Dir(_out_dir_);
		//chdir(_out_dir_.c_str());
		//std::cout<<" Done \n 	Printing Histograms";
	//}
	//if(_check_ == 0){
		//std::cout<<"Pass the check?";
		Histogram::WQ2_Print(output_dir, envi_);
		//Histogram::Fid_Print(const std::string& output_dir, envi_);
		//Histogram::SF_Print(const std::string& output_dir, envi_);
		//Histogram::DT_Print(const std::string& output_dir, envi_);
		//Histogram::CC_Print(const std::string& output_dir, envi_);
		//Histogram::MM_Print(const std::string& output_dir, envi_);
		//Histogram::XY_Print(const std::string& output_dir,envi_);
		//Histogram::Fid_Det_Print(const std::string& output_dir,envi_);
		//Friend_Print(const std::string& output_dir,envi_);
		//Histogram::Cross_Print(const std::string& output_dir,envi_);
	//}	
}

//W Qsquared plots
int Histogram::W_binning(float W_){
  int j = 0;
  float top, bot; 
  for(int i = 1; i < 30; i++){
    top = Wbin_start + i*Wbin_res;//constants.hpp
    bot = top - Wbin_res; 
    if(W_ < top && W_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::p_binning(float p_){
  int j = 0;
  float top, bot; 
  for(int i = 1; i < 26; i++){//Changed from 12 to 26 7/7/20
    top = pbin_start + i*pbin_res;//constants.hpp
    bot = top - pbin_res; 
    if(p_ < top && p_ >= bot){
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

void Histogram::WQ2_Make(std::shared_ptr<Environment> _envi ){
	if(_envi->was_WQ2_plot() && _envi->was_fit_type()!=1){
		std::vector<long> space_dims(4);
		space_dims[0] = 11; //Electron Cuts
		space_dims[1] = 6; //Topologies
		space_dims[2] = 2; //Recon vs. Thrown
		space_dims[3] = 2; //nweight vs. weight

		CartesianGenerator cart(space_dims);//CartesianGenerator.hpp
		char hname[100]; 

		while(cart.GetNextCombination()){
			if(_envi->Environment::was_sim()){
				if(cart[0]==10 && cart[1]!=0){
					sprintf(hname,"W_Q2_%s_%s_%s_%s",eid_cut[cart[0]],topologies[cart[1]],throw_stat[cart[2]],w_stat[cart[3]]); //constants.h and otherwise writing the specific cut to the right plot
			    	WQ2_hist[cart[0]][cart[1]][cart[2]][cart[3]] = std::make_shared<TH2F>( hname, hname, WQxres, WQxmin, WQxmax, WQyres, WQymin, WQymax); // constants.h
			    	WQ2_made_hist[cart[0]][cart[1]][cart[2]][cart[3]]=true;
				}else if( cart[0]!= 10 && cart[1]==0){
					sprintf(hname,"W_Q2_%s_%s_%s_%s",eid_cut[cart[0]],topologies[cart[1]],throw_stat[cart[2]],w_stat[cart[3]]); //constants.h and otherwise writing the specific cut to the right plot
			    	WQ2_hist[cart[0]][cart[1]][cart[2]][cart[3]] = std::make_shared<TH2F>( hname, hname, WQxres, WQxmin, WQxmax, WQyres, WQymin, WQymax); // constants.h	
			    	WQ2_made_hist[cart[0]][cart[1]][cart[2]][cart[3]]=true;		
		    	}
			}else if(cart[2] == 0){
				if(cart[0]==10 && cart[1]!=0){
					sprintf(hname,"W_Q2_%s_%s_%s",eid_cut[cart[0]],topologies[cart[1]],w_stat[cart[3]]); //constants.h and otherwise writing the specific cut to the right plot
			    	WQ2_hist[cart[0]][cart[1]][cart[2]][cart[3]] = std::make_shared<TH2F>( hname, hname, WQxres, WQxmin, WQxmax, WQyres, WQymin, WQymax); // constants.h
			    	WQ2_made_hist[cart[0]][cart[1]][cart[2]][cart[3]]=true;
				}else if( cart[0]!= 10 && cart[1]==0){
					sprintf(hname,"W_Q2_%s_%s_%s",eid_cut[cart[0]],topologies[cart[1]],w_stat[cart[3]]); //constants.h and otherwise writing the specific cut to the right plot
			    	WQ2_hist[cart[0]][cart[1]][cart[2]][cart[3]] = std::make_shared<TH2F>( hname, hname, WQxres, WQxmin, WQxmax, WQyres, WQymin, WQymax); // constants.h	
			    	WQ2_made_hist[cart[0]][cart[1]][cart[2]][cart[3]]=true;		
		    	}
			}
		}
	}
}
//hist_->Histogram::WQ2_Fill(envi_, 0, 0, W_, Q2, 1.0,thr);
void Histogram::WQ2_Fill(std::shared_ptr<Environment> _envi, int top, int cut, float W_, float Q2_, float weight_, int thr){
	if(_envi->was_WQ2_plot() && _envi->was_fit_type()!=1){
		if(WQ2_made_hist[cut][top][thr][0]){
			WQ2_hist[cut][top][thr][0]->Fill(W_,Q2_,1.0);
		}else{
			std::cout<<"Incorrectly filled WQ2 at: " <<cut <<" " <<top <<" " <<thr <<" " <<0 <<std::endl;
		}
		if(WQ2_made_hist[cut][top][thr][1]){
			WQ2_hist[cut][top][thr][1]->Fill(W_,Q2_,weight_);
		}else{
			std::cout<<"Incorrectly filled WQ2 at: " <<cut <<" " <<top <<" " <<thr <<" " <<1 <<std::endl;
		}
	}
}

void Histogram::WQ2_Write(std::shared_ptr<Environment> _envi){
	if(_envi->was_WQ2_plot() && _envi->was_fit_type()!=1){
		std::cout<<"WQ2 Plots: "; 
		char dir_name[100];
		TDirectory* dir_WQ2 = RootOutputFile->mkdir("W vs. Q2");
		dir_WQ2->cd();
		TDirectory* dir_WQ2_sub[6];
		for(int top = 0; top < 6; top ++){
			sprintf(dir_name,"W Q2 %s",topologies[top]);
			dir_WQ2_sub[top] = dir_WQ2->mkdir(dir_name);
			dir_WQ2_sub[top]->cd();
			for(int cut = 0; cut < 11; cut++){
				if(_envi->was_sim()){
					for(int thr = 0; thr < 2; thr++){
						if((top==0 && cut !=10) || (top!=0 && cut==10)){
							if(WQ2_made_hist[cut][top]){
								WQ2_hist[cut][top][thr][0]->SetXTitle("W (GeV)");
								WQ2_hist[cut][top][thr][0]->SetYTitle("Q^{2} (GeV^{2}");
								WQ2_hist[cut][top][thr][0]->SetOption("Colz");
								WQ2_hist[cut][top][thr][0]->Write();
								WQ2_hist[cut][top][thr][1]->SetXTitle("W (GeV)");
								WQ2_hist[cut][top][thr][1]->SetYTitle("Q^{2} (GeV^{2}");
								WQ2_hist[cut][top][thr][1]->SetOption("Colz");
								WQ2_hist[cut][top][thr][1]->Write();
							}else{
								std::cout<<std::endl <<"WQ2 would have segfaulted: " <<cut <<" " <<top;
							}
						}
					}
				}else{
					if((top==0 && cut !=10) || (top!=0 && cut==10)){
						if(WQ2_made_hist[cut][top]){
							WQ2_hist[cut][top][0][0]->SetXTitle("W (GeV)");
							WQ2_hist[cut][top][0][0]->SetYTitle("Q^{2} (GeV^{2}");
							WQ2_hist[cut][top][0][0]->SetOption("Colz");
							WQ2_hist[cut][top][0][0]->Write();
							WQ2_hist[cut][top][0][1]->SetXTitle("W (GeV)");
							WQ2_hist[cut][top][0][1]->SetYTitle("Q^{2} (GeV^{2}");
							WQ2_hist[cut][top][0][1]->SetOption("Colz");
							WQ2_hist[cut][top][0][1]->Write();
						}else{
							std::cout<<std::endl <<"WQ2 would have segfaulted: " <<cut <<" " <<top;
						}
					}
				}	
			}
			dir_WQ2_sub[top]->Close();
		}
		dir_WQ2->Close();
		std::cout<<"Done" <<std::endl;
	}
}

void Histogram::WQ2_Print(const std::string& output_dir, std::shared_ptr<Environment> envi_){
	Double_t w = 600;
	Double_t h = 1200;
	std::string _curr_dir_ = fun::get_current_dir();
	std::string print_name = "$curr/$name/Plots/WQ2.png";
	//std::string print_name = "/Users/cmc/Desktop/analysis/analysis_clas6/bin/$name/Plots/WQ2.png";// = "WQ2.png"; 
	//if(envi_->was_print() > 0){
		fun::replace(print_name, "$curr", _curr_dir_);
		fun::replace(print_name,"$name",output_dir);
		//sprintf(print_name,"/Users/cmc/Desktop/analysis/analysis_clas6/bin/%s/WQ2.png",output_dir.c_str());//Users/cmc/Desktop/analysis/analysis_clas6/bin/turst/
		//TDirectory * iWQ2 = HistImageFile->mkdir("WQ2 images");
		//iWQ2->cd();

		def = new TCanvas("cwq","cwq",w,h);
		def->Divide(1,2);
		def->cd(1);
		WQ2_hist[0][0][0][0]->Draw("colz");
		def->cd(2);
		WQ2_hist[1][0][0][0]->Draw("colz");
		
		def->SaveAs(print_name.c_str());
		//def->SaveAs("WQ2.png");
		delete def;

	//}
}

//Fiducial Cuts
void Histogram::Fid_Make(std::shared_ptr<Environment> _envi){
	char hname[100];
	char par_cut[100]; 
	float wtop, wbot,ptop,pbot; 
	std::vector<long> space_dims(7);
	space_dims[0] = 7; //Sector
	space_dims[1] = 4; //species
	space_dims[2] = 11;//cuts
	space_dims[3] = 30;//W Binning
	space_dims[4] = 26;//Momentum binning
	space_dims[5] = 6; //Topology
	space_dims[6] = 2; //Cut vs. anti-cut

	CartesianGenerator cart(space_dims);

	//char * fid_cuts[];
	int length = 0;
	char p_cut[100]; 
	int num_hist = 0; 
	//std::cout<<std::endl <<"Here we are. Making Fiducial histograms" <<std::endl;

	while(cart.GetNextCombination()){
		if(_envi->was_fid_plot(cart[1]) && fun::hist_fitting(cart[1],cart[2],cart[3],cart[4],_envi->was_fit_type())){
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
}

void Histogram::Fid_Fill(std::shared_ptr<Environment> _envi,int top, float theta, float phi, int part, int cut, int cutvanti, float W_, float p_){
	//std::cout<< std::endl <<"We are inside the Fid Fill function ";
	//std::cout<< std::endl <<"top: " <<top <<" theta: " <<theta <<" phi: " <<phi <<" part: " <<part <<" cut: " <<cut <<" cva: " <<cutvanti <<" W: " <<W_ <<" p: " <<p_ <<std::endl; 
	//std::cout<< std::endl <<"filling: " <<physics::get_sector(phi) <<" " <<part <<" " <<cut <<" " <<0 <<" " <<Histogram::p_binning(p_) <<" " <<top <<" " <<cutvanti <<std::endl;
	if(_envi->was_fid_plot(part) && fun::hist_fitting(part,cut,0,0,_envi->was_fit_type())){
		float phic = physics::phi_center(phi);
		//std::cout<<"Filling p binning" <<std::endl;
		if(p_binning(p_) != 0 && top==0 && _envi->was_fitting()){ //All W with p binning
			Fid_fill_hist[physics::get_sector(phi)][part][cut][0][Histogram::p_binning(p_)][0][cutvanti]=true;
			if(Fid_made_hist[physics::get_sector(phi)][part][cut][0][Histogram::p_binning(p_)][0][cutvanti]){
				Fid_hist[0][part][cut][0][Histogram::p_binning(p_)][0][cutvanti]->Fill(phic,theta);//All Sectors
				Fid_hist[physics::get_sector(phi)][part][cut][0][Histogram::p_binning(p_)][0][cutvanti]->Fill(phic,theta);//Individual Sectors			
			}else{
				std::cout<<"would have segfaulted ";
				std::cout<<std::endl <<"For Fid Plot: " <<physics::get_sector(phi) <<" " <<part <<" " <<cut <<" " <<0 <<" " <<Histogram::p_binning(p_) <<" " <<top <<" " <<cutvanti <<std::endl;
			}
		}
		//std::cout<<"Filling W binning" <<std::endl;
		if(W_binning(W_) != 0 && top==0 && _envi->was_fitting()){ //All p with W binning
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
}

void Histogram::Fid_Write(std::shared_ptr<Environment> _envi){
	if(_envi->was_fid_plot(0) || _envi->was_fid_plot(1) || _envi->was_fid_plot(2)|| _envi->was_fid_plot(3)){//If nothing has being plotted then don't make the directories
		std::cout<<"Fid Plots: "; 
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
			if(_envi->was_fid_plot(i)){
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
		}

		
		std::vector<long> space_dims(7);
		space_dims[0] = 7; //Sector
		space_dims[1] = 4; //species
		space_dims[2] = 11;//cuts
		space_dims[3] = 30;//W Binning
		space_dims[4] = 26;//Momentum binning
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
			if(_envi->was_fid_plot(cart[1]) && fun::hist_fitting(cart[1],cart[2],cart[3],cart[4],_envi->was_fit_type())){
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
									Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetXTitle("#phi (degrees)");
									Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetYTitle("#theta (degrees)");
									Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetOption("COLZ");
									Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->Write();
								}else{
									std::cout<<std::endl <<"Would have segfaulted Fid plot " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<" " <<cart[6] <<std::endl;
								}
							}
						}else if(cart[5]==0){//All P bins
							if(Fid_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]){
								par_fid_sec[cart[1]][cart[2]][0][cart[0]]->cd();
								Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetXTitle("#phi (degrees)");
								Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetYTitle("#theta (degrees)");
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
									Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetXTitle("#phi (degrees)");
									Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetYTitle("#theta (degrees)");
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
								Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetXTitle("#phi (degrees)");
								Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]][cart[6]]->SetYTitle("#theta (degrees)");
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
		std::cout<<"Done" <<std::endl;
	}
}

//Sampling Fraction Cuts
void Histogram::SF_Make(std::shared_ptr<Environment> _envi){
	if(_envi->was_sf_plot()){
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
			if((cart[0] == 10 && cart[3] != 0) || (cart[0] != 10 && cart[3] == 0) && fun::hist_fitting(0,cart[0],cart[1],0,_envi->was_fit_type())){//Topology only matters for event selection cut
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
}

void Histogram::SF_Fill(std::shared_ptr<Environment> _envi,int top, float p, float en, int cut, int cva, float W_, int sec){
	if(_envi->was_sf_plot() && fun::hist_fitting(0,cut,Histogram::W_binning(W_),0,_envi->was_fit_type())){
		SF_hist[cut][Histogram::W_binning(W_)][sec][top][cva]->Fill(p,en/p);
		SF_hist[cut][Histogram::W_binning(W_)][0][top][cva]->Fill(p,en/p);
	}
}

void Histogram::SF_Write(std::shared_ptr<Environment> _envi){
	if(_envi->was_sf_plot()){
		std::cout<<"SF Plots: "; 
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
			if(fun::hist_fitting(0,cart[0],cart[1],0,_envi->was_fit_type())){
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
		std::cout<<"Done" <<std::endl;
	}
}


//Delta T Cuts
void Histogram::DT_Make(std::shared_ptr<Environment> _envi){
	char hname[100];
	std::vector<long> space_dims(6);
	space_dims[0] = 4;  //species
	space_dims[1] = 7;  //Cuts
	space_dims[2] = 30; //W Binning
	space_dims[3] = 7; //Sector
	space_dims[4] = 6; //topology
	space_dims[5] = 2; //cut v anti


	float bot,top;

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		if(_envi->was_dt_plot(cart[0]) && fun::hist_fitting(cart[0],cart[1],cart[2],0,_envi->was_fit_type())){
			if((cart[1] == 6 && cart[4]!=0) || (cart[1]!=6 && cart[4] ==0)){
				if(cart[2] == 0){
					sprintf(hname,"%s_DeltaT_%s_%s_%s_W:ALL_%s",species[cart[0]],hid_cut[cart[1]],cut_ver[cart[5]],sec_list[cart[3]],topologies[cart[4]]);
					DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = std::make_shared<TH2F>(hname,hname, DTxres, DTxmin, DTxmax, DTyres, DTymin, DTymax);
					DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = true; 	
				}else if(cart[3]==0 & cart[4] == 0 && cart[1]!=6){//Looking at specific Cuts on fiducial and pre cut regimes 
					top = Wbin_start + cart[2]*Wbin_res;
					bot = top - Wbin_res;
					sprintf(hname,"%s_DeltaT_%s_%s_%s_W:%f-%f_%s",species[cart[0]],hid_cut[cart[1]],cut_ver[cart[5]],sec_list[cart[3]],bot,top,topologies[cart[4]]);
					DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = std::make_shared<TH2F>(hname,hname, DTxres, DTxmin, DTxmax, DTyres, DTymin, DTymax);
					DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = true;
				}else{
					DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = false;
				}
			}
		}
	}
}

void Histogram::DT_Fill(std::shared_ptr<Environment> _envi,int top, int part, float p, float d, float t, float d0, float t0, int cut, int anti, float W_, int sec){
	if(_envi->was_dt_plot(part) && fun::hist_fitting(part,cut,Histogram::W_binning(W_),0,_envi->was_fit_type())){
		float dt = physics::delta_t(part, p, d, t, d0, t0);
		if(sec < 7){//If sector dependent
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
}
			
void Histogram::DT_Fill(std::shared_ptr<Environment> _envi,int top, int part, float p, float dt, int cut, int anti, float W_, int sec){
	//std::cout<<"Filling DT Plot" <<std::endl;
	if(_envi->was_dt_plot(part) && fun::hist_fitting(part,cut,Histogram::W_binning(W_),0,_envi->was_fit_type())){
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
}

void Histogram::DT_Write(std::shared_ptr<Environment> _envi){
	bool do_dt_plots = false;
	for(int i = 0; i< 4; i++){
		if(_envi->was_dt_plot(i)){
			do_dt_plots = true;
		}
	}
	if(do_dt_plots){
		std::cout <<"Writing DT Plots: ";
		char dir_name[100]; 
		TDirectory* DT_plot = RootOutputFile->mkdir("DT_plots");
		TDirectory* par_dt[4][9][2][8][6];//Particle, cut, anti
		//std::cout<<"Did I get here?" <<std::endl; 
		DT_plot->cd(); 
		for(int i = 0; i<4; i++){//Species
			if(_envi->was_dt_plot(i)){
				sprintf(dir_name,"%s_DT_plots",species[i]);
				par_dt[i][0][0][0][0]= DT_plot->mkdir(dir_name);
				DT_dir_made[i][0][0][0][0]=true;
				//std::cout<<"Made pointer:" <<i <<" 0 0 0 0" <<std::endl;
				for(int j = 1; j < 8; j++){//cut
					sprintf(dir_name,"%s_DT_%s",species[i],hid_cut[j-1]);
					par_dt[i][j][0][0][0] = par_dt[i][0][0][0][0]->mkdir(dir_name);
					DT_dir_made[i][j][0][0][0]=true;
					//std::cout<<"    Made Pointers" <<i <<j <<" 0 0 0"<<std::endl;
					//W Dependence
					sprintf(dir_name,"%s_DT_%s_%s",species[i],hid_cut[j-1],W_dep_list[1]);
					par_dt[i][j][1][0][0] = par_dt[i][j][0][0][0]->mkdir(dir_name);
					DT_dir_made[i][j][1][0][0]=true;
					//std::cout<<"    Made Pointers" <<i <<j <<" 1 0 0"<<std::endl;
					for(int k = 1; k < 8; k++){ //Sector 
						sprintf(dir_name,"%s_DT_%s_%s",species[i],hid_cut[j-1],sec_list[k-1]);
						par_dt[i][j][0][k][0] = par_dt[i][j][0][0][0]->mkdir(dir_name);
						DT_dir_made[i][j][0][k][0]=true;
						//std::cout<<"    Made Pointers" <<i <<j <<" 0 " <<k <<" 0"<<std::endl;
						//std::cout<<"Sector Pointer" <<std::endl;
						}
					if(j == 7){//Event cut
						for(int l = 1; l < 6; l++){ //topology 
							sprintf(dir_name,"%s_DT_%s_%s",species[i],hid_cut[j-1],topologies[l]);
							par_dt[i][j][0][0][l] = par_dt[i][j][0][0][0]->mkdir(dir_name);
							DT_dir_made[i][j][0][0][l]=true;
							//std::cout<<"    Made Pointers" <<i <<j <<" 0 0 " <<k<<std::endl;
							//std::cout<<"Topology Pointer" <<std::endl;

						}
					}
				}
			}
		}
		//std::cout<<"Made it through making directories" <<std::endl;

		std::vector<long> space_dims(6);
		space_dims[0] = 4;  //species
		space_dims[1] = 7;  //Cuts
		space_dims[2] = 30; //W Binning
		space_dims[3] = 7; //Sector
		space_dims[4] = 6; //topology
		space_dims[5] = 2; //Cut v anti

		CartesianGenerator cart(space_dims);
		while(cart.GetNextCombination()){
			if(_envi->was_dt_plot(cart[0]) && fun::hist_fitting(cart[0],cart[1],cart[2],0,_envi->was_fit_type())){
				if(DT_dir_made[0][0][0][0][0]){
					par_dt[cart[0]][0][0][0][0]->cd();//Main folder for particles
					if(DT_dir_made[cart[0]][cart[1]+1][0][0][0]){
						par_dt[cart[0]][cart[1]+1][0][0][0]->cd();//Get into those CUTS
						//std::cout <<"      Now Writing in " <<cart[0] <<" " <<cart[1]+1 <<" 0 0 0"<<std::endl;
						//All W, Sectors, and Combined Topology for Event selection, but still by Cut
						if(cart[2] ==0 && cart[3] == 0 && ((cart[1]!=6 && cart[4]==0)||(cart[1]==6 && cart[4]==5))){
							if(DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]){
								//std::cout<<"DT shouldn't have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
								DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetXTitle("Momentum (GeV)");
								DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetYTitle("Delta T (ns)");
								DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetOption("COLZ");
								DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->Write();	
							}else{
								std::cout<<"DT Would have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
							}		
						}
						
						//For W Range, but all sectors, combine topology for event selection, all sectors
						//if(_envi->was_W_dep_plot()){	
							if(cart[2]!=0 && cart[3] == 0 && ((cart[1]!=6 && cart[4]==0))){//||(cart[1]==6 && cart[4]==5))){ //Issue when trying to look at event selection for this. Not sure why, but getting rid of it solved it *shrug* 9/12/19
								//std::cout <<"          trying to Write in " <<cart[0] <<" " <<cart[1]+1 <<" 1 0 0"<<std::endl;
								if(DT_dir_made[cart[0]][cart[1]+1][1][0][0]){
									par_dt[cart[0]][cart[1]+1][0][0][0]->cd();//Get into those CUTS
									par_dt[cart[0]][cart[1]+1][1][0][0]->cd();
									//std::cout <<"   We are writing" <<std::endl;
									if(DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]){
										//std::cout<<"DT shouldn't have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
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
						//}
						//For Sector Range, but all W, combine topology for event selection
						if(cart[2] ==0 && ((cart[1]!=6 && cart[4]==0)||(cart[1]==6 && cart[4]==5))){
							//std::cout <<"          Trying to Write in " <<cart[0] <<" " <<cart[1]+1 <<" 0 " <<cart[3]+1 <<" 0"<<std::endl;
							if(DT_dir_made[cart[0]][cart[1]+1][0][cart[3]+1][0]){
								par_dt[cart[0]][cart[1]+1][0][0][0]->cd();//Get into those CUTS
								par_dt[cart[0]][cart[1]+1][0][cart[3]+1][0]->cd();
								//std::cout <<"  We are writing "<<std::endl;
								if(DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]){
									//std::cout<<"DT shouldn't have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
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
							if(DT_dir_made[cart[0]][cart[1]+1][0][0][cart[4]]){//Changed cart[4] + 1 to just cart[4] 5/5/20
								par_dt[cart[0]][cart[1]+1][0][0][0]->cd();//Get into those CUTS
								par_dt[cart[0]][cart[1]+1][0][0][cart[4]]->cd();
								//std::cout <<"  we are writing"<<std::endl;
								if(DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]){
									//std::cout<<"DT shouldn't have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
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
		std::cout<<" Done" <<std::endl;
	}
}

//Min CC Cuts
void Histogram::CC_Make(std::shared_ptr<Environment> _envi){
	
	if(_envi->was_cc_plot() && !(_envi->was_sim())){
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
}

void Histogram::CC_Fill(std::shared_ptr<Environment> _envi,int top, int sec, int segm, int nphe, int cut, int anti){
	//std::cout<<std::endl <<"top: " <<top <<" sec: " <<sec <<" segm: " <<segm <<" nphe: " <<nphe <<" cut: " <<cut <<" anti: " <<anti <<std::endl;
	//std::cout<<"filling plot";
	//std::cout<<std::endl <<sec <<" " <<detect::cc_segment(segm) <<" " <<cut <<" " <<detect::cc_lrc(segm) <<" " <<top <<" " <<anti <<std::endl;
	if(_envi->was_cc_plot() && !(_envi->was_sim())){
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

}
void Histogram::CC_Write(std::shared_ptr<Environment> _envi){
	if(_envi->was_cc_plot() && !(_envi->was_sim())){
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

	

}
//Missing Mass Cuts
void Histogram::MM_Make(std::shared_ptr<Environment> _envi){
	bool do_mm_plots = false;
	for(int i = 0; i< 4; i++){
		if(_envi->was_MM_plot(i)){
			do_mm_plots = true;
		}
	}
	if(do_mm_plots){
		char hname[100];
		std::vector<long> space_dims(4);
		space_dims[0] = 4;  //Topology
		space_dims[1] = 3;  //Cuts {pre, cut, anti}
		space_dims[2] = 2; //squared vs linear
		space_dims[3] = 2; //Fitting vs. not fitting plots

		float MMmin = NAN; 
		float MMmax = NAN; 

		CartesianGenerator cart(space_dims);

		while(cart.GetNextCombination()){
			if(_envi->was_MM_plot(cart[0])){
				sprintf(hname,"%s_MM_%s_%s_%s",topologies[cart[0]+1],basic_cut[cart[1]],MM_sq[cart[2]],fit_q[cart[3]]);
				MMmin = MMxmin2[cart[0]];
				MMmax = MMxmax2[cart[0]];
				MM_hist[cart[0]][cart[1]][cart[2]][cart[3]] = std::make_shared<TH1F>(hname,hname,MMxres,MMmin,MMmax);
			}
		}
	}
}
void Histogram::MM_Fill(std::shared_ptr<Environment> _envi,int top, float mm, int cut, int square,bool fit){
	if(_envi->was_MM_plot(top)){
		MM_hist[top][cut][square][1]->Fill(mm);
		if(fit){
			MM_hist[top][cut][square][0]->Fill(mm);
		}
	}
}
void Histogram::MM_Write(std::shared_ptr<Environment> _envi){
	bool do_mm_plots = false;
	for(int i = 0; i< 4; i++){
		if(_envi->was_MM_plot(i)){
			do_mm_plots = true;
		}
	}
	if(do_mm_plots){
		char dir_name[100];
		TDirectory * MM_plot = RootOutputFile->mkdir("MM plots");
		MM_plot->cd();
		TDirectory * MM_dir[4][3][2];//top, cut, fit
		for(int k = 0; k<2; k++){
			for(int i = 0; i < 4; i++ ){
				if(_envi->was_MM_plot(i)){
					sprintf(dir_name,"MM %s %s",topologies[i+1],fit_q[k]);
					MM_dir[i][0][k] = MM_plot->mkdir(dir_name);
					for(int j = 1; j <3; j++){
						sprintf(dir_name,"MM %s %s %s",topologies[i+1],MM_sq[j-1],fit_q[k]);
						MM_dir[i][j][k] = MM_dir[i][0][k]->mkdir(dir_name);
					}
				}
			}
		}
		for(int p = 0; p<2 ; p++){
			for(int i = 0; i < 4; i++ ){//Topology without "All"
				if(_envi->was_MM_plot(i)){
					MM_dir[i][0][p]->cd();
					for(int j = 1; j <3; j++){//linear vs squared
						MM_dir[i][j][p]->cd();
						for( int k = 0 ; k<3; k++){//Cuts
							sprintf(dir_name,"MM %s (GeV %s)",MM_sq[j-1],MM_sq[j-1]);
							MM_hist[i][k][j-1][p]->SetXTitle(dir_name);
							MM_hist[i][k][j-1][p]->SetYTitle("Events");
							MM_hist[i][k][j-1][p]->Write();
						}
					}
				}
			}
		}
	}
	

}

void Histogram::XY_Make(std::shared_ptr<Environment> envi_){
	char hname[100];
	std::vector<long> space_dims(3);
	space_dims[0] = 3;  //Detectors {CC,SC,EC}
	space_dims[1] = 2;  //Simulation {Recon,Thrown}
	space_dims[2] = 4; //species

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		if(cart[1]==0){
			sprintf(hname,"XY_%s_%s",detectors[cart[0]+1],species[cart[2]]);
			XY_hist[cart[0]][cart[1]][cart[2]]= std::make_shared<TH2F>(hname,hname,XYres,-XYmax[cart[0]],XYmax[cart[0]],XYres,-XYmax[cart[0]],XYmax[cart[0]]);
		}else if(envi_->was_sim() && cart[0]==0){
			sprintf(hname,"XY_%s_Thrown",species[cart[2]]);
			XY_hist[cart[0]][cart[1]][cart[2]]= std::make_shared<TH2F>(hname,hname,XYres,-XYmax[cart[0]],XYmax[cart[0]],XYres,-XYmax[cart[0]],XYmax[cart[0]]);
		}
	}
}

void Histogram::XY_Fill(std::shared_ptr<Environment> envi_, int species_, float x_, float y_, int detector_ , bool thrown_){
	//{1,2,3} -> {cc,sc,ed}
	if(thrown_ && envi_->was_sim()){
		XY_hist[0][1][species_]->Fill(x_,y_);
	}else if(!thrown_){
		XY_hist[detector_][0][species_]->Fill(x_,y_);
	}
}

void Histogram::XY_Write(std::shared_ptr<Environment> envi_){
	char dir_name[100];
	TDirectory * XY_plot = RootOutputFile->mkdir("XY plots");
	XY_plot->cd();
	for(int l = 0; l < 4; l++){
		for(int k = 0; k<3; k++){
			XY_hist[k][0][l]->SetXTitle("X Position (mm)");//Check what these actual units are...
			XY_hist[k][0][l]->SetYTitle("Y Position (mm)");
			XY_hist[k][0][l]->Write();
		}
		if(envi_->was_sim()){
			XY_hist[0][1][l]->SetXTitle("X Position (mm)");//Check what these actual units are...
			XY_hist[0][1][l]->SetYTitle("Y Position (mm)");
			XY_hist[0][1][l]->Write();
		}
	}
}

void Histogram::Fid_Det_Make(std::shared_ptr<Environment> envi_){
	char hname[100];
	std::vector<long> space_dims(3);//[3][4][7][2]
	space_dims[0] = 3;  //Detectors {CC,SC,EC}
	space_dims[1] = 4;  //Species
	space_dims[2] = 7; //All, Sector

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		if(cart[2]==0){
			if(envi_->was_sim()){
				sprintf(hname,"%sFid_%s_Sector:All_Sim",detectors[cart[0]+1],species[cart[1]]);
				Fid_Det_hist[cart[0]][cart[1]][cart[2]]= std::make_shared<TH2F>(hname,hname,FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
			}else{
				sprintf(hname,"%sFid_%s_Sector:All",detectors[cart[0]+1],species[cart[1]]);
				Fid_Det_hist[cart[0]][cart[1]][cart[2]]= std::make_shared<TH2F>(hname,hname,FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
			}
		}else{
			if(envi_->was_sim()){
				sprintf(hname,"%sFid_%s_Sector:%d_Sim",detectors[cart[0]+1],species[cart[1]],cart[2]);
				Fid_Det_hist[cart[0]][cart[1]][cart[2]]= std::make_shared<TH2F>(hname,hname,FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
			}else{
				sprintf(hname,"%sFid_%s_Sector:%d",detectors[cart[0]+1],species[cart[1]],cart[2]);
				Fid_Det_hist[cart[0]][cart[1]][cart[2]]= std::make_shared<TH2F>(hname,hname,FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
			}
		}
		
	}
}

void Histogram::Fid_Det_Fill(std::shared_ptr<Environment> envi_, int species_, float theta_, float phi_, int sector_, int detector_){
	Fid_Det_hist[detector_][species_][sector_]->Fill(phi_,theta_);
	Fid_Det_hist[detector_][species_][0]->Fill(phi_,theta_);
}

void Histogram::Fid_Det_Write(std::shared_ptr<Environment> envi_){
	char dir_name[100];
	TDirectory * FidD_plot = RootOutputFile->mkdir("Fid Detector plots");
	FidD_plot->cd();
	TDirectory * FidD_dir[4][4];
	for(int spe = 0; spe < 4; spe++){
		FidD_dir[spe][0] = FidD_plot->mkdir(species[spe]);
		FidD_dir[spe][0]->cd();
		for(int det = 0; det < 3; det++){
			sprintf(dir_name,"Det_%s_%s",detectors[det+1],species[spe]);
			FidD_dir[spe][det+1] = FidD_dir[spe][0]->mkdir(dir_name);
			FidD_dir[spe][det+1]->cd();
			for(int sec = 0; sec < 7; sec++){
			Fid_Det_hist[det][spe][sec]->SetXTitle("Phi (deg)");
			Fid_Det_hist[det][spe][sec]->SetYTitle("Theta (deg)");
			Fid_Det_hist[det][spe][sec]->Write();
			}
		}
	}
}




void Histogram::Friend_Make(std::shared_ptr<Environment> _envi){
	//If we decide we want to fill this stuff
	if(_envi->was_Friend_plot()){
		char hname[100];
		sprintf(hname,"2#pi_off_proton_#Delta^{++}");
		Double_t xmin1[7] = {0.5,_W_min,_Q2_min,_MM_min[0],_theta_min,_alpha_min,_phi_min};
		Double_t xmax1[7] = {5.5,_W_max,_Q2_max,_MM_max[0],_theta_max,_alpha_max,_phi_max};
		Double_t xmin2[7] = {-0.5,_W_min,_Q2_min,_MM_min[1],_theta_min,_alpha_min,_phi_min};
		Double_t xmax2[7] = {5.5,_W_max,_Q2_max,_MM_max[1],_theta_max,_alpha_max,_phi_max};
		Double_t xmin3[7] = {-0.5,_W_min,_Q2_min,_MM_min[2],_theta_min,_alpha_min,_phi_min};
		Double_t xmax3[7] = {5.5,_W_max,_Q2_max,_MM_max[2],_theta_max,_alpha_max,_phi_max};
		Friend[0] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin1,xmax1);
		//Friend[0]->Sumw2();//Must be called before error bars can be placed. 
		sprintf(hname,"2#pi_off_proton_#rho");
		Friend[1] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin2,xmax2);
		//Friend[1]->Sumw2();//Must be called before error bars can be placed. 
		sprintf(hname,"2#pi_off_proton_#Delta^{0}");
		Friend[2] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin3,xmax3);
		//Friend[2]->Sumw2();//Must be called before error bars can be placed. 
	}
}

int Histogram::Friend_W_binning(float W_){
	int Wbin = W_binning(W_)-1;
	return Wbin;
}

int Histogram::Friend_Q2_binning(float Q2_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _Friend_bins[2]; i++){
    top = _Q2_min + (i+1)*((_Q2_max-_Q2_min)/_Friend_bins[2]);//constants.hpp
    bot = top - ((_Q2_max-_Q2_min)/_Friend_bins[2]); 
    if(Q2_ < top && Q2_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::Friend_MM_binning(float MM_, int chan){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _Friend_bins[3]; i++){
    top = _MM_min[chan] + (i+1)*((_MM_max[chan]-_MM_min[chan])/_Friend_bins[3]);//constants.hpp
    bot = top - ((_MM_max[chan]-_MM_min[chan])/_Friend_bins[3]); 
    if(MM_ < top && MM_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::Friend_theta_binning(float theta_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _Friend_bins[4]; i++){
    top = _theta_min + (i+1)*((_theta_max-_theta_min)/_Friend_bins[4]);//constants.hpp
    bot = top - ((_theta_max-_theta_min)/_Friend_bins[4]); 
    if(theta_ < top && theta_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::Friend_alpha_binning(float alpha_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _Friend_bins[5]; i++){
    top = _alpha_min + (i+1)*((_alpha_max-_alpha_min)/_Friend_bins[5]);//constants.hpp
    bot = top - ((_alpha_max-_alpha_min)/_Friend_bins[5]); 
    if(alpha_ < top && alpha_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::Friend_phi_binning(float phi_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _Friend_bins[6]; i++){
    top = _phi_min + (i+1)*((_phi_max-_phi_min)/_Friend_bins[6]);//constants.hpp
    bot = top - ((_phi_max-_phi_min)/_Friend_bins[6]); 
    if(phi_ < top && phi_ >= bot){
      j = i; 
    }
  }
  return j; 
}


int * Histogram::Friend_binning(int top, float W_, float Q2_, float MM_, float theta_, float alpha_, float phi_ , int channel){
	int x[7]; 
	int test = 0; 
	bool in_region = false; 
	x[0] = top; //{pmiss,pipmiss,pimmiss,zeromiss,all}
	x[1] = Friend_W_binning(W_);
	x[2] = Friend_Q2_binning(Q2_);
	x[3] = Friend_MM_binning(MM_,channel);
	x[4] = Friend_theta_binning(theta_);
	x[5] = Friend_alpha_binning(alpha_);
	x[6] = Friend_phi_binning(phi_);
	/*std::cout<<std::endl <<"Filling Friend in channel " <<channel <<std::endl; 
	std::cout<<"top = " <<top <<" bin: " <<x[0] <<std::endl;
	std::cout<<"W = " <<W_ <<" bin: " <<x[1] <<std::endl;
	std::cout<<"Q2 = " <<Q2_ <<" bin: " <<x[2] <<std::endl;
	std::cout<<"MM = " <<MM_ <<" bin: " <<x[3] <<std::endl;
	std::cout<<"theta = " <<theta_ <<" bin: " <<x[4] <<std::endl;
	std::cout<<"alpha = " <<alpha_ <<" bin: " <<x[5] <<std::endl;
	std::cout<<"phi = " <<phi_ <<" bin: " <<x[6] <<std::endl;
	*/
	for(int i = 0; i < 7; i++){
		if(x[i] >= 0){
			test++;
		}
	}
	if(test < 7){
		x[0] = -1;
	}
	return x; 
}

void Histogram::Friend_Fill(std::shared_ptr<Environment> _envi, int top_, float W_, float Q2_, float MM_, float theta_, float alpha_, float phi_ , int chan_, float weight_){
	int *y = Friend_binning(top_,W_,Q2_,MM_,theta_,alpha_,phi_,chan_);
	Double_t x[7] = {(float)top_, W_, Q2_, MM_, theta_, alpha_, phi_};
	if(y[0]>=0){
		Friend[chan_]->Fill(x,weight_);
	}
}

void Histogram::Friend_Write(std::shared_ptr<Environment> _envi){
	if(_envi->Environment::was_Friend_plot()){
		char dir_name[100];
		TDirectory * Friend_plot = RootOutputFile->mkdir("Friend plots");
		Friend_plot->cd();
		TH1D* MM_proj[3]; 
		TH1D* alpha_proj[3];
		TH1D* theta_proj[3];

		for(int i = 0; i< 3; i++){
			MM_proj[i] = Friend[i]->Projection(3);
			alpha_proj[i] = Friend[i]->Projection(5);
			theta_proj[i] = Friend[i]->Projection(4);
			Friend[i]->Write(); 
			MM_proj[i]->Draw();
			MM_proj[i]->Write();
			alpha_proj[i]->Draw();
			alpha_proj[i]->Write();
			theta_proj[i]->Draw();
			theta_proj[i]->Write();
		}
	}
	
}

float Histogram::Friend_bin_reverse(int var_, int bin_, int channel_){//This only works for equally spaced bins. Will need rework when bins have varied sizes
	float val = NAN;
	float max = NAN;
	float min = NAN;
	int bins = -1; 
	switch(var_){
		case 0: 
			max = 4.5;
			min = -0.5;
			bins = _Friend_bins[var_]; 
		break; //Top
		case 1:
			max = _W_max;
			min = _W_min;
			bins = _Friend_bins[var_] ; 
		break;//W
		case 2: 
			max = _Q2_max;
			min = _Q2_min;
			bins = _Friend_bins[var_]; 
		break;//Q2
		case 3: 
			max = _MM_max[channel_];
			min = _MM_min[channel_];
			bins = _Friend_bins[var_]; 
		break;//MM
		case 4: 
			max = _theta_max;
			min = _theta_min;
			bins = _Friend_bins[var_]; 
		break;//Theta
		case 5: 
			max = _alpha_max;
			min = _alpha_min;
			bins = _Friend_bins[var_]; 
		break;//Alpha
		case 6: 
			max = _phi_max;
			min = _phi_min;
			bins = _Friend_bins[var_]; 
		break;//phi
	}
	for(int i = 0; i < bins; i++){
		if(bin_ == i){
			val = ((max - min)/bins)*(i + 0.5) + min; 
		}
	}
	return val;
}

void Histogram::Cross_Make(std::shared_ptr<Environment> envi_){
	char hname[100];
	
	/*
	0 - Pmiss Only
	1 - Pipmiss Only
	2 - Pimmiss Only
	3 - Zeromiss Only
	4 - Zeromiss + 3
	5 - Zeromiss + 2
	6 - Zeromiss + 1
	7 - Pmiss + Pipmiss
	8 - Pmiss + Pimmiss
	9 - Pipmiss + Pimmiss
	10 - No Zeromiss + 3
	11 - Multiples Pmiss
	12 - Multiples Pipmiss
	13 - Multiples Pimmiss
	14 - Multiples Zeromiss
	*/
	sprintf(hname,"Topology_Crossing_no_weights");
	Cross_hist[0] = std::make_shared<TH1F>(hname,hname, 15, -0.5, 14.5);
	sprintf(hname,"Topology_Crossing__weights");
	Cross_hist[1] = std::make_shared<TH1F>(hname,hname, 15, -0.5, 14.5);
}

void Histogram::Cross_Fill(std::shared_ptr<Environment> envi_, int gevnt_[4], float weight_){
	bool top_pass[4] = {false,false,false,false};
	bool top_mult[4] = {false,false,false,false};
	int num_hadmiss = 0;
	int tot_evnt = 0; 
	for(int i = 0; i< 4; i++){
		if(gevnt_[i] > 0){
			top_pass[i] = true;
			if(i != 3){
				num_hadmiss += 1; 
			}
			if(gevnt_[i] > 1){
				top_mult[i] = true;
			}
		}
		tot_evnt += gevnt_[i];
	}
	
	if(top_pass[0] && !top_pass[1] && !top_pass[2] && !top_pass[3] && !top_mult[0]){//Pmiss Only
		Cross_hist[0]->Fill(0.0);
		Cross_hist[1]->Fill(0.0,weight_);
	}else if(!top_pass[0] && top_pass[1] && !top_pass[2] && !top_pass[3] && !top_mult[1]){//Pipmiss only
		Cross_hist[0]->Fill(1.0);
		Cross_hist[1]->Fill(1.0,weight_);
	}else if(!top_pass[0] && !top_pass[1] && top_pass[2] && !top_pass[3] && !top_mult[2]){//Pimmiss only
		Cross_hist[0]->Fill(2.0);
		Cross_hist[1]->Fill(2.0,weight_);
	}else if(!top_pass[0] && !top_pass[1] && !top_pass[2] && top_pass[3] && !top_mult[3]){//Zeromiss only
		Cross_hist[0]->Fill(3.0);
		Cross_hist[1]->Fill(3.0,weight_);
	}else if(top_pass[0] && top_pass[1] && top_pass[2] && top_pass[3] && !top_mult[3]){//Zeromiss +3
		Cross_hist[0]->Fill(4.0);
		Cross_hist[1]->Fill(4.0,weight_);
	}else if((num_hadmiss == 2) && top_pass[3] && !top_mult[3]){//Zeromiss +2
		Cross_hist[0]->Fill(5.0);
		Cross_hist[1]->Fill(5.0,weight_);
	}else if((num_hadmiss == 1) && top_pass[3] && !top_mult[3]){//Zeromiss +1
		Cross_hist[0]->Fill(6.0);
		Cross_hist[1]->Fill(6.0,weight_);
	}else if(top_pass[0] && top_pass[1]){//P + Pip
		Cross_hist[0]->Fill(7.0);
		Cross_hist[1]->Fill(7.0,weight_);
	}else if(top_pass[0] && top_pass[2]){//P + Pip
		Cross_hist[0]->Fill(8.0);
		Cross_hist[1]->Fill(8.0,weight_);
	}else if(top_pass[1] && top_pass[2]){//P + Pip
		Cross_hist[0]->Fill(9.0);
		Cross_hist[1]->Fill(9.0,weight_);
	}else if(!top_pass[3] && (num_hadmiss == 3)){//P + Pip
		Cross_hist[0]->Fill(10.0);
		Cross_hist[1]->Fill(10.0,weight_);
	}else if(top_mult[0]){//P + Pip
		Cross_hist[0]->Fill(11.0);
		Cross_hist[1]->Fill(11.0,weight_);
	}else if(top_mult[1]){//P + Pip
		Cross_hist[0]->Fill(12.0);
		Cross_hist[1]->Fill(12.0,weight_);
	}else if(top_mult[2]){//P + Pip
		Cross_hist[0]->Fill(13.0);
		Cross_hist[1]->Fill(13.0,weight_);
	}else if(top_mult[3]){//P + Pip
		Cross_hist[0]->Fill(14.0);
		Cross_hist[1]->Fill(14.0,weight_);
	}

}

void Histogram::Cross_Write(std::shared_ptr<Environment> envi_){
	//char dir_name[100];
	std::cout<<"Cross Plots: ";
	TDirectory * Cross_plot = RootOutputFile->mkdir("Cross Plots");
	Cross_plot->cd();
	Cross_hist[0]->SetXTitle("Event Topology Mixing");
	Cross_hist[0]->SetYTitle("Events");
	Cross_hist[0]->Write();
	Cross_hist[1]->SetXTitle("Event Topology Mixing");
	Cross_hist[1]->SetYTitle("Events");
	Cross_hist[1]->Write();
	std::cout<<"Done" <<std::endl;
}


void Histogram::Acceptance_Make(std::shared_ptr<Environment> _envi){
	//If we decide we want to fill this stuff
	if(_envi->was_Friend_plot() && _envi->was_sim()){
		char hname[100];
		sprintf(hname,"2#pi_off_proton_#Delta^{++}");
		Double_t xmin1[7] = {-0.5,_W_min,_Q2_min,_MM_min[0],_theta_min,_alpha_min,_phi_min};
		Double_t xmax1[7] = {5.5,_W_max,_Q2_max,_MM_max[0],_theta_max,_alpha_max,_phi_max};
		Double_t xmin2[7] = {-0.5,_W_min,_Q2_min,_MM_min[1],_theta_min,_alpha_min,_phi_min};
		Double_t xmax2[7] = {5.5,_W_max,_Q2_max,_MM_max[1],_theta_max,_alpha_max,_phi_max};
		Double_t xmin3[7] = {-0.5,_W_min,_Q2_min,_MM_min[2],_theta_min,_alpha_min,_phi_min};
		Double_t xmax3[7] = {5.5,_W_max,_Q2_max,_MM_max[2],_theta_max,_alpha_max,_phi_max};
		Friend[0] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin1,xmax1);
		//Friend[0]->Sumw2();//Must be called before error bars can be placed. 
		sprintf(hname,"2#pi_off_proton_#rho");
		Friend[1] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin2,xmax2);
		//Friend[1]->Sumw2();//Must be called before error bars can be placed. 
		sprintf(hname,"2#pi_off_proton_#Delta^{0}");
		Friend[2] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin3,xmax3);
		//Friend[2]->Sumw2();//Must be called before error bars can be placed. 
	}
}

//Fill this with Thrown data after having already made the Friend
void Histogram::Acceptance_Fill(std::shared_ptr<Environment> _envi, int top_, float W_, float Q2_, float MM_, float theta_, float alpha_, float phi_ , int chan_, float weight_){
	int *y = Friend_binning(top_,W_,Q2_,MM_,theta_,alpha_,phi_,chan_);
	Double_t x[7] = {(float)top_, W_, Q2_, MM_, theta_, alpha_, phi_};
	if(y[0]>=0){
		Friend[chan_]->Fill(x,weight_);
	}
}
/*
void Histogram::Event_Particle_Hist(std::shared_ptr<Environment> envi_, const Particle p1, float W_, int top_, int par_, bool pass_){
	std::cout<<"		Filling Particle Event" <<std::endl;
	int pass = -1; 
	if(pass_){
		pass = 0; 
	}else{
		pass = 1; 
	}
	//Electron
	if(pass != -1){
		if(par_ == 0 && _pid[0]){
			std::cout<<"			Electron cc_segm: " <<_cc_seg <<std::endl;
			Fid_Fill(envi_,top_+1,_theta,_phi,0,10,pass,W_,_p);
			SF_Fill(envi_,top_+1,_p,_etot,10,pass,W_,physics::get_sector(_phi));
			CC_Fill(envi_,top_+1,physics::get_sector(_phi),_cc_seg,_nphe,10,pass);
		}else if(par_ ==0){
			std::cout<<"			Electron issue, friend" <<std::endl;
		}
		if(par_ != 0 && _pid[par_]){
			Fid_Fill(envi_,top_+1,_theta,_phi,par_,6,pass,W_,_p);
			DT_Fill(envi_,top_+1,par_,_p,_dt[par_],6,pass,W_,physics::get_sector(_phi));
		}else if(par_ !=0){
			std::cout<<"			Hadron issue, friend" <<std::endl;
		}
	}
}
/*
void Histogram::Fill_EID(std::shared_ptr<Particle> par, float W_, float Q2_){
	bool san, fid, cc, sf;
	if(par->Particle::Is_sanity_pass(0)){
		Histogram::WQ2_Fill(0,1,W_,Q2_);
		Histogram::Fid_Fill(0,par->Particle::par_theta(),par->Particle::par_phi(),0,1,0,W_,par->Particle::par_p());
		Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),1,0,W_,sector[0]);
		Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),1,0);
		if(par->Particle::Is_fid_pass(0)){
			fid = true;
			Histogram::WQ2_Fill(0,2,W_,Q2_);
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,2,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),2,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),2,0);
		}else{
			fid = false;
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,2,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),2,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),2,1);
		}
		if(par->Particle::Is_cc_pass()){
			cc = true;
			Histogram::WQ2_Fill(0,4,W_,Q2_);
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,4,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),4,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),4,0);
		}else{
			cc = false;
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,4,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),4,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),4,1);
		}
		if(par->Is_sf_pass() && par->Is_min_cc_pass()){
			sf = true;
			Histogram::WQ2_Fill(0,3,W_,Q2_);
			Histogram::Fid_Fill(0,par->Particle::par_theta(),par->Particle::par_phi(),0,3,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),3,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),3,0);
		}else{
			sf = false;
			Histogram::Fid_Fill(0,par->Particle::par_theta(),par->Particle::par_phi(),0,3,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),3,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),3,1);
		}
		if(fid && sf){
			Histogram::WQ2_Fill(0,5,W_,Q2_);
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,5,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),5,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),5,0);
		}else{
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,5,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),5,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),5,1);
		}
		if(fid && cc){
			Histogram::WQ2_Fill(0,6,W_,Q2_);
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,6,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),6,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),6,0);
		}else{
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,6,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),6,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),6,1);
		}
		if(sf && cc){
			Histogram::WQ2_Fill(0,7,W_,Q2_);
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,7,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),7,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),7,0);
		}else{
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,7,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),7,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),7,1);
		}
		if(sf && fid && cc){
			Histogram::WQ2_Fill(0,8,W_,Q2_);
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,8,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),8,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),8,0);
			_elec = physics::Make_4Vector(par->Particle::par_p(),data->Branches::cx(0),data->Branches::cy(0),data->Branches::cz(0),me);
			good_electron++;
		}else{
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,8,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),8,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),8,1);
		}
		if(par->Particle::Bank()==ELECTRON){
			Histogram::WQ2_Fill(0,9,W_,Q2_);
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,9,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),9,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),9,0);
		}
		else{
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,9,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),9,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),9,1);
		}
	}else{
		san = false;
		Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,1,1,W_,par->Particle::par_p());
		Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),1,1,W_,sector[0]);
		Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),1,1);
	}
}

void Histogram::Fill_HID(std::shared_ptr<Particle> par){

}*/

