#include "particle.hpp"

Particle::Particle(){
}


void Particle::Fill_Particle(std::shared_ptr<Branches> data, int par_idx){
	idx = par_idx; 
	if(idx > 0){
		hadron = true;
	}else{
		hadron = false;
	}
	//Hits in Detector Systems
	dc = data->Branches::dc(idx);
	cc = data->Branches::dc(idx);
	sc = data->Branches::dc(idx);
	ec = data->Branches::dc(idx);
	//Particle Information 
	cx = data->Branches::cx(idx);
	cy = data->Branches::cy(idx);
	cz = data->Branches::cz(idx);
	p = data->Branches::p(idx);
	etot = data->Branches::etot(idx);
	bank = data->Branches::id(idx);
	
	cc_segm = data->Branches::cc_segm(idx);
	cc_sect = data->Branches::cc_sect(idx);
	//Conversions to more useful values
	phi = physics::get_phi(cx,cy);
	theta = physics::get_theta(cz);
	sector = physics::get_sector(phi);
	sf = etot/p; 
	if(hadron){
		sc_r = data->Branches::sc_r(idx);
		sc_t = data->Branches::sc_t(idx);
		sc_r0 = data->Branches::sc_r(0);
		sc_t0 = data->Branches::sc_t(0);
		for(int i = 0; i<3; i++){
			dt_t[i] = physics::delta_t(i,p,sc_r,sc_t,sc_r0,sc_t0);
		}
	}else{
		pim_ele = false;
		pro_pip = false; 
	}
}

void Particle::EID_Cuts(){
	//Sanity Cuts
	if(dc > 0 && cc > 0 && sc > 0 && ec > 0 && q == -1){
		sanity_pass[0] = true;
		//Fiducial Cuts
		if(cuts::fid_cut(0,p,cx,cy,cz)){
			fid_pass[0] = true; 
		}else{
			fid_pass[0] = false; 
		}
		//Sampling Fraction Cuts
		if(cuts::sf_cut(p,etot,cx,cy)){
			sf_pass = true;
		}
		else{
			sf_pass = false;
		}
		//Min CC cuts
		if(cuts::min_cc(cc_segm,cc_sect,nphe)){
			cc_pass = true; 
		}else{
			cc_pass = false; 
		}
		//Min EC cuts
		if(cuts::min_ec(etot)){
			min_ec_pass = true; 
		}else{
			min_ec_pass = false; 
		}
	}
	//Full EID pass
	if(sanity_pass[0] && fid_pass[0] && sf_pass && cc_pass && min_ec_pass){
		pid[0] = true; 
	}else{
		pid[0] = false;
	}
}

void Particle::Ele_Pim_Resolve(){//Requires EID_Cuts to have been performed 
	if(cc_pass && min_ec_pass && sf_pass){//Don't need the fiducial one
		double_id_resolved = false;
	}else{
		double_id_resolved = true;
		pim_ele = false;
	}
}

void Particle::HID_Cuts(){
	int q_goal = 0; 
	//Cycle between all relevant hadrons
	for(int i = 0; i<3; i++){
		if(i == 2){
			q_goal = -1;
			if(i == 2){
				dt = dt_t[2];
			}
		}else{
			q_goal = 1;
		}
		//Sanity Cuts
		if(dc > 0 && sc > 0 && q==q_goal){
			sanity_pass[i+1] = true;
			//Fiducial Cut
			if(cuts::fid_cut(i+1,p,cx,cy,cz)){
				fid_pass[i+1] = true; 
			}else{
				fid_pass[i+1] = false; 
			}
			//Delta T cut
			if(cuts::delta_t_cut(i,p,sc_r,sc_t,sc_r0,sc_t0)){
				dt_pass[i] = true;
			}else{
				dt_pass[i] = false;
			}
		}
		else{
			sanity_pass[i+1] = false;
		}
		if(fid_pass[i+1] && sanity_pass[i+1] && dt_pass[i]){
			pid[i+1] = true;
		}else{
			pid[i+1] = false;
		}
	}
	if(pid[1] && pid[2]){
		double_id_resolved = false; 
		pro_pip = true;
	}else if(!pid[3]){
		double_id_resolved = true;
	}else{
		pim_ele = true;
		double_id_resolved = false;
	}
}

bool Particle::Is_Elec(){
	return pid[0];
}
bool Particle::Is_Pro(){
	return pid[1];
}
bool Particle::Is_Pip(){
	return pid[2];
}
bool Particle::Is_Pim(){
	return pid[3];
}

void Particle::Pro_Pip_Resolved(int par){//Requires HID_Cuts to have been performed
	switch(par){
		case 0:
			pid[1] = true;
			pid[2] = false;
			double_id_resolved = true;
			pro_pip = false;
			dt = dt_t[0];
		break;
		case 1:
			pid[1] = false;
			pid[2] = true;
			double_id_resolved = true;
			pro_pip = false;
			dt = dt_t[1];
		break;
		default:
		std::cout<<"Improperly Entered a parameter";
		break; 
	}
}

bool Particle::Is_fid_pass(int par){
	return fid_pass[par];
}

bool Particle::Is_cc_pass(){
	return cc_pass;
}

bool Particle::Is_min_ec_pass(){
	return min_ec_pass;
}

bool Particle::Is_dt_pass(int par){
	return dt_pass[par];
}

bool Particle::Is_sanity_pass(int par){
	return sanity_pass[par];
}

bool Particle::Is_sf_pass(){
	return sf_pass;
}

bool Particle::Pro_Pip_Prob(){
	return pro_pip;
}

int Particle::par_idx(){
	return idx;
}
int Particle::par_nphe(){
	return nphe;
}
int Particle::par_id(){
	return id;
}
int Particle::par_sector(){
	return sector;
}
int Particle::par_cc_sector(){
	return cc_sect;
}
int Particle::par_cc_segm(){
	return cc_segm;
}

float Particle::par_dt(int par){
	return dt_t[par];
}
float Particle::par_theta(){
	return theta;
}
float Particle::par_phi(){
	return phi;
}
float Particle::par_etot(){
	return etot;
}
float Particle::par_p(){
	return p;
}

TLorentzVector Particle::Par_4Vec(){
	return k1;
}
