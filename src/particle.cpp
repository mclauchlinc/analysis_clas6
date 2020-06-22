#include "particle.hpp"

//Initialize the particle by filling it with everything we want
Particle::Particle(){

}

void Particle::Fill_Particle(std::shared_ptr<Branches> data_, int par_idx_, int set_, bool sim_, bool thrown_){
	float m = -99; 
	_idx = par_idx_;
	_set = set_;
	_sim = sim_;
	_thrown = thrown_;
	if(data_->Branches::dc(_idx) > 0 && data_->Branches::sc(_idx) > 0){//These are things necessary for all events
		for(int i = 0; i< 4; i++){
			switch(i){
				case 0:
					m = me; 
				break;
				case 1: 
					m = mp; 
				break; 
				case 2: 
					m = mpi; 
				break;
				case 3: 
					m = mpi; 
				break;
			}
			_dt[i] = physics::delta_t(i,data_,_idx);
		}
		_p = data_->Branches::p(_idx);
		_theta = physics::get_theta(_idx,data_);
		_phi = physics::get_phi(_idx,data_);
		_sf = data_->Branches::etot(_idx)/data_->Branches::p(_idx);
		_q = data_->Branches::q(_idx);
		//std::cout<<"	Particle " <<_idx <<" has charge of " <<_q <<std::endl;
	}else if(_thrown){
		_p = data_->Branches::mcp(_idx);
		_theta = data_->Branches::mctheta(_idx);
		_phi = data_->Branches::mcphi(_idx);
	}
	if(!_sim && data_->Branches::cc(_idx) > 0){
		_cc_seg = data_->Branches::cc_segm(_idx);
		_nphe = data_->Branches::nphe(_idx);
		//std::cout<<"		Has a CC value " <<std::endl;
		if(_idx == 0){
			//std::cout<<"		Here is the cc_segm for this particle: " <<_cc_seg <<std::endl;
		}
	}
	if(data_->Branches::ec(_idx)>0){
		_etot = data_->Branches::etot(_idx);
	}
}


void Particle::EID(std::shared_ptr<Branches> data_, std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, float W_){
	bool _trigger_ = false;
	int thr = 0; 
	float Q2 = physics::Qsquared(_set, data_, _thrown);
	//Sanity Cuts
	if(_idx == 0){
		hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,0,0,W_,_p);
		hist_->Histogram::SF_Fill(envi_,0,_p,_etot,0,0,W_,physics::get_sector(_phi));
		hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,0,0);
		hist_->Histogram::WQ2_Fill(envi_, 0, 0, W_, Q2, thr);
		_trigger_ = true;
	}
	if(cuts::e_sanity(data_,envi_,_idx)){
		_sanity_pass[0] = true;
		if(_trigger_){
			hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,1,0,W_,_p);
			hist_->Histogram::SF_Fill(envi_,0,_p,_etot,1,0,W_,physics::get_sector(_phi));
			hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,1,0);
			hist_->Histogram::WQ2_Fill(envi_, 0, 1, W_, Q2, thr);
		}
		//Fiducial Cuts
		if(cuts::e_fid(data_,envi_,_idx)){
			_fid_pass[0] = true; 
			if(_trigger_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,2,0,W_,_p);
				hist_->Histogram::SF_Fill(envi_,0,_p,_etot,2,0,W_,physics::get_sector(_phi));
				hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,2,0);
				hist_->Histogram::WQ2_Fill(envi_, 0, 2, W_, Q2, thr);
			}
		}else{
			_fid_pass[0] = false;
			if(_trigger_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,2,1,W_,_p);
				hist_->Histogram::SF_Fill(envi_,0,_p,_etot,2,1,W_,physics::get_sector(_phi));
				hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,2,1);

			}
		}
		//Sampling Fraction Cuts
		if(cuts::e_sf(data_,envi_,_idx)){
			_sf_pass = true;
			if(_trigger_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,3,0,W_,_p);
				hist_->Histogram::SF_Fill(envi_,0,_p,_etot,3,0,W_,physics::get_sector(_phi));
				hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,3,0);
				hist_->Histogram::WQ2_Fill(envi_, 0, 3, W_, Q2, thr);
			}
		}
		else{
			_sf_pass = false;
			if(_trigger_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,3,1,W_,_p);
				hist_->Histogram::SF_Fill(envi_,0,_p,_etot,3,1,W_,physics::get_sector(_phi));
				hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,3,1);
			}
		}
		//Min CC cuts
		if(cuts::e_cc(data_,envi_,_idx)){
			_cc_pass = true; 
			if(_trigger_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,4,0,W_,_p);
				hist_->Histogram::SF_Fill(envi_,0,_p,_etot,4,0,W_,physics::get_sector(_phi));
				hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,4,0);
				hist_->Histogram::WQ2_Fill(envi_, 0, 4, W_, Q2, thr);
			}
		}else{
			_cc_pass = false; 
			if(_trigger_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,4,1,W_,_p);
				hist_->Histogram::SF_Fill(envi_,0,_p,_etot,4,1,W_,physics::get_sector(_phi));
				hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,4,1);
			}
		}
		//Min EC cuts
		if(cuts::min_ec(data_->Branches::etot(_idx))){
			_min_ec_pass = true; 
		}else{
			_min_ec_pass = false; 
		}
		//Combo cuts
		if(_fid_pass[0] && _sf_pass){
			if(_trigger_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,5,0,W_,_p);
				hist_->Histogram::SF_Fill(envi_,0,_p,_etot,5,0,W_,physics::get_sector(_phi));
				hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,5,0);
				hist_->Histogram::WQ2_Fill(envi_, 0, 5, W_, Q2, thr);
			}
		}else{
			if(_trigger_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,5,1,W_,_p);
				hist_->Histogram::SF_Fill(envi_,0,_p,_etot,5,1,W_,physics::get_sector(_phi));
				hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,5,1);
			}
		}
		if(_fid_pass[0] && _cc_pass){
			if(_trigger_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,6,0,W_,_p);
				hist_->Histogram::SF_Fill(envi_,0,_p,_etot,6,0,W_,physics::get_sector(_phi));
				hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,6,0);
				hist_->Histogram::WQ2_Fill(envi_, 0, 6, W_, Q2, thr);
			}
		}else{
			if(_trigger_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,6,1,W_,_p);
				hist_->Histogram::SF_Fill(envi_,0,_p,_etot,6,1,W_,physics::get_sector(_phi));
				hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,6,1);
			}
		}
		if(_sf_pass && _cc_pass){
			if(_trigger_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,7,0,W_,_p);
				hist_->Histogram::SF_Fill(envi_,0,_p,_etot,7,0,W_,physics::get_sector(_phi));
				hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,7,0);
				hist_->Histogram::WQ2_Fill(envi_, 0, 7, W_, Q2, thr);
			}
		}else{
			if(_trigger_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,7,1,W_,_p);
				hist_->Histogram::SF_Fill(envi_,0,_p,_etot,7,1,W_,physics::get_sector(_phi));
				hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,7,1);
			}
		}
	} else{
		if(_trigger_){
			hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,1,1,W_,_p);
			hist_->Histogram::SF_Fill(envi_,0,_p,_etot,1,1,W_,physics::get_sector(_phi));
			hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,1,1);
		}
	}
	//Full EID pass
	if(cuts::eid(data_,envi_,_idx)){
		_pid[0] = true;
		_ided = true;
		if(_trigger_){
			hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,8,0,W_,_p);
			hist_->Histogram::SF_Fill(envi_,0,_p,_etot,8,0,W_,physics::get_sector(_phi));
			hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,8,0); 
			hist_->Histogram::WQ2_Fill(envi_, 0, 8, W_, Q2, thr);
		}
	}else{
		_pid[0] = false;
		if(_trigger_){
			hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,8,1,W_,_p);
			hist_->Histogram::SF_Fill(envi_,0,_p,_etot,8,1,W_,physics::get_sector(_phi));
			hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,8,1);
		}
	}
	//Bank ID
	if(data_->Branches::id(_idx) == ELECTRON && _idx == 0){
		if(_trigger_){
			hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,9,0,W_,_p);
			hist_->Histogram::SF_Fill(envi_,0,_p,_etot,9,0,W_,physics::get_sector(_phi));
			hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,9,0);
			hist_->Histogram::WQ2_Fill(envi_, 0, 9, W_, Q2, thr);
		}
	}else{
		if(_trigger_){
			hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,0,9,1,W_,_p);
			hist_->Histogram::SF_Fill(envi_,0,_p,_etot,9,1,W_,physics::get_sector(_phi));
			hist_->Histogram::CC_Fill(envi_,0,physics::get_sector(_phi),_cc_seg,_nphe,9,1);
		}
	}
}


void Particle::HID(std::shared_ptr<Branches> data_, std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, float W_){
	bool _ntrig_ = false;//If not a trigger particle
	if(_idx != 0){
		_ntrig_ = true;
	}
	//Cycle between all relevant hadrons
	int _pid_[3] = {PROTON,PION,-PION}; 
	//std::cout<<"	Particle " <<_idx <<" charge: " <<data_->Branches::q(_idx) <<std::endl;
	for(int i = 0; i<3; i++){
		if(_ntrig_){
			hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,i+1,0,0,W_,_p);
			hist_->Histogram::DT_Fill(envi_,0,i+1,_p,_dt[i+1],0,0,W_,physics::get_sector(_phi));
		}
		//Sanity Cuts
		
		if(cuts::h_sanity(data_,envi_,_idx,i)){
				//std::cout<<"	Particle " <<_idx <<" Pass H sanity for " <<species[i+1]  <<" with charge: " <<data_->Branches::q(_idx) <<std::endl;
			_sanity_pass[i+1] = true;
			if(_ntrig_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,i+1,1,0,W_,_p);
				hist_->Histogram::DT_Fill(envi_,0,i+1,_p,_dt[i+1],1,0,W_,physics::get_sector(_phi));
			}
			//Fiducial Cut
			if(cuts::h_fid(data_,envi_,_idx,i)){
				_fid_pass[i+1] = true; 
				if(_ntrig_){
					hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,i+1,2,0,W_,_p);
					hist_->Histogram::DT_Fill(envi_,0,i+1,_p,_dt[i+1],2,0,W_,physics::get_sector(_phi));
				}
			}else{
				_fid_pass[i+1] = false; 
				if(_ntrig_){
					hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,i+1,2,1,W_,_p);
					hist_->Histogram::DT_Fill(envi_,0,i+1,_p,_dt[i+1],2,1,W_,physics::get_sector(_phi));
				}
			}
			//Delta T cut
			//std::cout<<"Delta t cut had " <<species[i+1] <<": " <<cuts::delta_t_cut(i+1,_p,_dt[i+1]) <<std::endl;
			if(cuts::delta_t_cut(i+1,_p,_dt[i+1])){//cuts::h_dt(data_,envi_,_idx,i)){
				_dt_pass[i+1] = true;
				if(_ntrig_){
					hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,i+1,3,0,W_,_p);
					hist_->Histogram::DT_Fill(envi_,0,i+1,_p,_dt[i+1],3,0,W_,physics::get_sector(_phi));
				}
			}else{
				_dt_pass[i+1] = false;
				if(_ntrig_){
					hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,i+1,3,1,W_,_p);
					hist_->Histogram::DT_Fill(envi_,0,i+1,_p,_dt[i+1],3,1,W_,physics::get_sector(_phi));
				}
			}
		}
		else{
			_sanity_pass[i+1] = false;
			if(_ntrig_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,i+1,1,1,W_,_p);
				hist_->Histogram::DT_Fill(envi_,0,i+1,_p,_dt[i+1],1,1,W_,physics::get_sector(_phi));
			}
		}
		if(_dt_pass[i+1] && _fid_pass[i+1] && _sanity_pass[i+1]){//cuts::hid(_set, data_, envi_, _idx, i)){
			_pid[i+1] = true;
			_ided = true; 
			if(_ntrig_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,i+1,4,0,W_,_p);
				hist_->Histogram::DT_Fill(envi_,0,i+1,_p,_dt[i+1],4,0,W_,physics::get_sector(_phi));
			}
		}else{
			_pid[i+1] = false;
			if(_ntrig_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,i+1,4,1,W_,_p);
				hist_->Histogram::DT_Fill(envi_,0,i+1,_p,_dt[i+1],4,1,W_,physics::get_sector(_phi));
			}
		}
		if(data_->Branches::id(_idx) == _pid_[i]){
			if(_ntrig_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,i+1,5,0,W_,_p);
				hist_->Histogram::DT_Fill(envi_,0,i+1,_p,_dt[i+1],5,0,W_,physics::get_sector(_phi));
			}

		}else{
			if(_ntrig_){
				hist_->Histogram::Fid_Fill(envi_,0,_theta,_phi,i+1,5,1,W_,_p);
				hist_->Histogram::DT_Fill(envi_,0,i+1,_p,_dt[i+1],5,1,W_,physics::get_sector(_phi));
			}
		}
	}
}

void Particle::PID(std::shared_ptr<Branches> data_, std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, float W_){
	//std::cout<<"Doing Particle ID" <<std::endl;
	if(_thrown){
		switch(data_->Branches::mcid(_idx)){
			case ELECTRON:
				_pid[0] = true;
				_ided = true;
				hist_->Histogram::WQ2_Fill(envi_, 0, 0, W_, physics::Qsquared(envi_->Environment::was_data_set(), data_, _thrown), 1);
			break;
			case PROTON:
				_pid[1] = true;
				_ided = true;
			break;
			case PION:
				_pid[2] = true;
				_ided = true;
			break;
			case -PION:
				_pid[3] = true;
				_ided = true;
			break;
		}
	}else{
		Particle::EID(data_,envi_,hist_,W_);
		Particle::HID(data_,envi_,hist_,W_);
	}
	//Negative Particles
	if(_pid[0] && _pid[3] && !_pid[1] && !_pid[2]){
		_id_crisis = true;
	}else{
		_id_crisis = false;
	}
	//Positive Particles
	if(!_pid[0] && !_pid[3] && _pid[1] && _pid[2]){
		_id_crisis = true;
	}else{
		_id_crisis = false;
	}
}



bool Particle::Pass_Sanity(int i){
	return _sanity_pass[i];
}
bool Particle::Pass_ec(){
	return _min_ec_pass;
}
bool Particle::Pass_fid(int i){
	return _fid_pass[i];
}
bool Particle::Pass_sf(){
	return _sf_pass;
}
bool Particle::Pass_cc(){
	return _cc_pass;
}
bool Particle::Pass_dt(int i){
	return _dt_pass[i];
}
bool Particle::Corr_p(){
	return _p_corr; 
}
int Particle::ID_crisis(){
	return _id_crisis; 
}
bool Particle::IDed(){
	return _ided;
}
bool Particle::Is_Sim(){
	return _sim;
}

bool Particle::Is_Thrown(){
	return _thrown;
}


bool Particle::Is_Elec(){
	return _pid[0];
}
bool Particle::Is_Pro(){
	return _pid[1];
}
bool Particle::Is_Pip(){
	return _pid[2];
}
bool Particle::Is_Pim(){
	return _pid[3];
}

float Particle::Get_p(){
	return _p; 
}
float Particle::Get_theta(){
	return _theta; 
}
float Particle::Get_phi(){
	return _phi; 
}

int Particle::Get_set(){
	return _set;
}

int Particle::Get_idx(){
	return _idx;
}

void Particle::Fill_Par_Event(std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, float W_, int top_, int par_, bool pass_){
	//std::cout<<"		Filling Particle Event" <<std::endl;
	int _pass_ = -1; 
	if(pass_){
		_pass_ = 0; 
	}else{
		_pass_ = 1; 
	}
	//Electron
	if(_pass_ != -1){
		if(par_ == 0 && _pid[0]){
			//std::cout<<"			Electron cc_segm: " <<_cc_seg <<std::endl;
			hist_->Histogram::Fid_Fill(envi_,top_+1,_theta,_phi,0,10,_pass_,W_,_p);
			hist_->Histogram::SF_Fill(envi_,top_+1,_p,_etot,10,_pass_,W_,physics::get_sector(_phi));
			hist_->Histogram::CC_Fill(envi_,top_+1,physics::get_sector(_phi),_cc_seg,_nphe,10,_pass_);
		}else if(par_ ==0){
			std::cout<<"			Electron issue, friend" <<std::endl;
		}
		if(par_ != 0 && _pid[par_]){
			hist_->Histogram::Fid_Fill(envi_,top_+1,_theta,_phi,par_,6,_pass_,W_,_p);
			hist_->Histogram::DT_Fill(envi_,top_+1,par_,_p,_dt[par_],6,_pass_,W_,physics::get_sector(_phi));
		}else if(par_ !=0){
			std::cout<<"			Hadron issue, friend" <<std::endl;
		}
	}
}

void Particle::Check_Particle(){
	std::cout<<"			Looking inside Particle " <<_idx <<std::endl;
	std::cout<<"				From set: " <<_set <<std::endl;
	std::cout<<"				Sim Status: " <<_sim <<std::endl;
	std::cout<<"				Thrown Status: " <<_thrown <<std::endl;
	std::cout<<"				p: " <<_p <<std::endl;
	std::cout<<"				q: " <<_q <<std::endl;
	std::cout<<"				dt: {" <<_dt[0] <<", " <<_dt[1] <<", " <<_dt[2] <<", " <<_dt[3] <<"}" <<std::endl; 
	std::cout<<"				theta: " <<_theta <<std::endl;
	std::cout<<"				phi: " <<_phi <<std::endl;
	std::cout<<"				sf: " <<_sf <<std::endl;
	std::cout<<"				etot: " <<_etot <<std::endl;
	std::cout<<"				cc_segm: " <<_cc_seg <<std::endl;
	std::cout<<"				nphe: " <<_nphe <<std::endl;
	std::cout<<"				sanity: {" <<_sanity_pass[0] <<", " <<_sanity_pass[1] <<", " <<_sanity_pass[2] <<", " <<_sanity_pass[3] <<"}" <<std::endl;
	std::cout<<"				ec min: " <<_min_ec_pass <<std::endl;
	std::cout<<"				fid pass: {" <<_fid_pass[0] <<", " <<_fid_pass[1] <<", " <<_fid_pass[2] <<", " <<_fid_pass[3] <<"}" <<std::endl;
	std::cout<<"				sf pass: " <<_sf_pass <<std::endl;
	std::cout<<"				cc pass: " <<_cc_pass <<std::endl;
	std::cout<<"				dt pass: {" <<_dt_pass[0] <<", " <<_dt_pass[1] <<", " <<_dt_pass[2] <<", " <<_dt_pass[3] <<"}" <<std::endl;
	std::cout<<"				p corr: " <<_p_corr <<std::endl;
	std::cout<<"				ID Crisis: " <<_id_crisis <<std::endl;
	std::cout<<"				PID: {" <<_pid[0] <<", " <<_pid[1] <<", " <<_pid[2] <<", " <<_pid[3] <<"}" <<std::endl;
	std::cout<<"				IDed: " <<_ided <<std::endl;
}
