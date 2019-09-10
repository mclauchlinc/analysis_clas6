#include "histogram.hpp"


Histogram::Histogram(const std::string& output_file){
	RootOutputFile = fun::Name_File(output_file);
	def = new TCanvas("def");

}

Histogram::~Histogram() { this->Write(); }

void Histogram::Write(){
	std::cout<< "Writing" <<std::endl;
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

}

void Histogram::WQ2_Make(){
	std::vector<long> space_dims(2);
	space_dims[0] = 10; //Electron Cuts
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
	float wtop, wbot,ptop,pbot; 
	std::vector<long> space_dims(5);
	space_dims[0] = 7; //Sector
	space_dims[1] = 4; //species
	space_dims[2] = 10;//cuts
	space_dims[3] = 30;//W Binning
	space_dims[4] = 12;//Momentum binning

	CartesianGenerator cart(space_dims);

	//char * fid_cuts[];
	int length = 0;
	char p_cut[100]; 

	while(cart.GetNextCombination()){
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
		if(cart[2]<length){
			if(cart[3] == 0){//All W bins
				sprintf(hname,"%s_Fid_%s_%s_W:ALL_%s",species[cart[1]],sec_list[cart[0]],Histogram::Part_cut(cart[1],cart[2]),p_cut);//constants.hpp
			}else{//Specific W Bins
				wtop = Wbin_start + (cart[3]*Wbin_res);//constants.hpp
				wbot = wtop - Wbin_res;//constants.hpp
				sprintf(hname,"%s_Fid_%s_%s_W:%f-%f_%s",species[cart[1]],sec_list[cart[0]],Histogram::Part_cut(cart[1],cart[2]),wbot,wtop,p_cut);//constants.hpp
			}
			Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]] = std::make_shared<TH2D>(hname,hname, FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
		}
	}
}

void Histogram::Fid_Fill(float theta, float phi, int part, int cut, float W_, float p_){
	float phic = physics::phi_center(phi);
	if(p_binning(p_) + W_binning(W_) != 0){//Both p and W dependence
		Fid_hist[0][part][cut][Histogram::W_binning(W_)][Histogram::p_binning(p_)]->Fill(phic,theta);//All Sectors
		Fid_hist[physics::get_sector(phi)][part][cut][Histogram::W_binning(W_)][Histogram::p_binning(p_)]->Fill(phic,theta);//Individual Sectors
	}
	if(p_binning(p_) != 0){ //All W with p binning
		Fid_hist[0][part][cut][0][Histogram::p_binning(p_)]->Fill(phic,theta);//All Sectors
		Fid_hist[physics::get_sector(phi)][part][cut][0][Histogram::p_binning(p_)]->Fill(phic,theta);//Individual Sectors
	}
	if(W_binning(W_) != 0){ //All p with W binning
		Fid_hist[0][part][cut][Histogram::W_binning(W_)][0]->Fill(phic,theta);//All sectors
		Fid_hist[physics::get_sector(phi)][part][cut][Histogram::W_binning(W_)][0]->Fill(phic,theta);//Individual Sectors
	}
	//All W and All Momentum
	Fid_hist[0][part][cut][0][0]->Fill(phic,theta);//All sectors, all W, all P
	Fid_hist[physics::get_sector(phi)][part][cut][0][0]->Fill(phic,theta);//Individual sectors all W, all p
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

	
	std::vector<long> space_dims(5);
	space_dims[0] = 7; //Sector
	space_dims[1] = 4; //species
	space_dims[2] = 10;//cuts
	space_dims[3] = 30;//W Binning
	space_dims[4] = 12;//Momentum binning

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
		if(cart[4] ==0 ){
			p_dep = false;
		}else{
			p_dep = true;
		}
		/*switch(cart[1]){
			case 0:
			ele_fid->cd();
			if(p_dep){
				if(w_dep){

				}else{}
			}else{
				if(w_dep){

				}
			}
			break;
			case 0:
			pro_fid->cd();
			break;
			case 0:
			pip_fid->cd();
			break;
			case 0:
			pim_fid->cd();
			break;
		}*/


		if(cart[1] == 0){ //Dealing with the different numbers of cuts for electrons vs. hadrons
			//fid_cuts = eid_cut;
			length = 10; 
		}else{
			//fid_cuts = hid_cut; 
			length = 6; 
		}
		if(cart[2]<length){
			if(cart[3] == 0){//All W bins
			}else{//Specific W Bins
			}
			Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("{phi} (degrees)");
			Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("{theta} (degrees)");
			Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
			Fid_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
		}
	}
}
//Sampling Fraction Cuts
void Histogram::SF_Make(){
	char hname[100];
	std::vector<long> space_dims(3);
	space_dims[0] = 10; //Electron Cuts
	space_dims[1] = 30; //W binning
	space_dims[2] = 7;  //Sector
}
void Histogram::SF_Fill(float p, float en, int cut){

}
void Histogram::SF_Write(){

}
//Delta T Cuts
void Histogram::DT_Make(){
	char hname[100];
	std::vector<long> space_dims(4);
	space_dims[0] = 3;  //species
	space_dims[1] = 6;  //Cuts
	space_dims[2] = 30; //W Binning
	space_dims[3] = 7; //Sector
}
void Histogram::DT_Fill(int part, float p, float d, float t, float d0, float t0, int cut, float W_){

}
void Histogram::DT_Write(){

}
//Min CC Cuts
void Histogram::CC_Make(){
	char hname[100];
	std::vector<long> space_dims(4);
	space_dims[0] = 6;  //Sector
	space_dims[1] = 18; //Segment
	space_dims[2] = 5;  //Cut
	space_dims[3] = 4;  //Side of detector
}
void Histogram::CC_Fill(int sec, int segm, int nphe, int cut){

}
void Histogram::CC_Write(){

}
//Missing Mass Cuts
void Histogram::MM_Make(){
	char hname[100];
	std::vector<long> space_dims(3);
	space_dims[0] = 5;  //Topology
	space_dims[1] = 3;  //Cuts
	space_dims[2] = 30; //W Binning
}
void Histogram::MM_Fill(int top, float mm, int cut){

}
void Histogram::MM_Write(){

}

