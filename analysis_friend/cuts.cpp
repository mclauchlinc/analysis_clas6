#include "cuts.hpp"



bool cuts::fid_e ( float p, float cx, float cy, float cz)
{
	bool in_fid = kFALSE;

	//degree Conversion
	float degree = 180.0/TMath::Pi();

	//Calculate angles and sector
	float theta = physics::get_theta(cz);
	float phi_c = physics::phi_center(physics::get_phi(cx, cy));
	int sector = physics::get_sector(physics::get_phi(cx, cy));

	//Creating the cut function
	float theta_cut = c1e + c2e / ((float)p+p_shift_e);
	float expon = c3e * TMath::Power((float)p,factor_e);
	float del_phi = c4e * TMath::Power((TMath::Sin((theta-theta_cut)/degree)),expon);
	

	//Actual application of the cut
	if(TMath::Abs(phi_c)<=del_phi && theta>=theta_cut)
	{
		in_fid = kTRUE;
	}

	return in_fid;
}

float cuts::phi_min(float theta, int idx){
	return -(a0mh[idx]*(1.0-TMath::Exp(-a1mh[idx]*(theta-a2mh[idx])))-a3mh[idx]);
}

float cuts::phi_max(float theta, int idx){
	return (a0xh[idx]*(1.0-TMath::Exp(-a1xh[idx]*(theta-a2xh[idx])))+a3xh[idx]);
}

bool cuts::fid_h ( float p, float cx, float cy, float cz)
{
	bool in_fid = kFALSE;

	//degree Conversion
	float degree = 180.0/TMath::Pi();

	//Calculate angles and sector
	//Phi is centered for the sector
	float theta = physics::get_theta(cz);
	float phi_c = physics::phi_center(physics::get_phi(cx, cy));
	int sector = physics::get_sector(physics::get_phi(cx, cy));

	int sec_indx = sector -1;

	//Actual application of the cut
	if(phi_c>=cuts::phi_min(theta , sec_indx) && phi_c<=cuts::phi_max(theta, sec_indx))
	{
		in_fid = kTRUE;
	}

	return in_fid;
}

bool cuts::fid_cut(int part, float p, float cx, float cy, float cz){
	bool pass = false; 
	if(part ==0){//Electron
		pass = fid_e(p,cx,cy,cz);
	}else{
		pass = fid_h(p,cx,cy,cz);
	}
	return pass;
}

//Proton Formulae from Arjun
float cuts::dt_p_low(float p){
	return DTL[0]+DTL[1]*p+DTL[2]*p*p+DTL[3]*p*p*p;
}

float cuts::dt_p_high(float p){
	return DTH[0]+DTH[1]*p+DTH[2]*p*p+DTH[3]*p*p*p;
}

bool cuts::delta_t_cut_iso(int part, int part_iso, float p, float d0, float d, float t0, float t) //Note: d and t need the sc_index 
{
	bool pass = false;

	float dt = physics::delta_t(part, p, d, t, d0, t0);
	//std::cout<<std::endl <<"delta t cut: p: " <<p <<" dt: " <<dt <<" cut_low: " <<cuts::dt_p_low(p) <<" cut_high:" <<cuts::dt_p_high(p) <<std::endl;  
	if(dt>cuts::dt_p_low(p) && dt<cuts::dt_p_high(p) )
	{
		pass = true;
	}else if(part == part_iso && TMath::Abs(dt) < 10){//part allowing a wiiiide delta t cut
		pass = true; 
	}
	return pass;
}

bool cuts::delta_t_cut(int part, float p, float d0, float d, float t0, float t) //Note: d and t need the sc_index 
{
	bool pass = false;

	float dt = physics::delta_t(part, p, d, t, d0, t0);
	//std::cout<<std::endl <<"delta t cut: p: " <<p <<" dt: " <<dt <<" cut_low: " <<cuts::dt_p_low(p) <<" cut_high:" <<cuts::dt_p_high(p) <<std::endl;  
	if(dt>cuts::dt_p_low(p) && dt<cuts::dt_p_high(p) )
	{
		pass = true;
	}else //if(part == 1 && TMath::Abs(dt) < 10){//part allowing a wiiiide delta t cut
	//	pass = true; 
	//}
	return pass;
}

bool cuts::delta_t_cut(int species_, float p_, float dt_){
	bool pass = false;
	//std::cout<<std::endl <<"delta t cut: p: " <<p <<" dt: " <<dt <<" cut_low: " <<cuts::dt_p_low(p) <<" cut_high:" <<cuts::dt_p_high(p) <<std::endl;  
	if(dt_>cuts::dt_p_low(p_) && dt_<cuts::dt_p_high(p_)){
		//std::cout<<"Low: " <<cuts::dt_p_low(p_) <<"	| Value: " <<dt_ <<"	|High: " <<cuts::dt_p_high(p_);
		pass = true;
		//std::cout<<"	| " <<species_ <<" Dt Pass: " <<pass <<std::endl;; 
	} //if(part == 1 && TMath::Abs(dt) < 10){//part allowing a wiiiide delta t cut
	//	pass = true; 
	//}
	return pass;
}

bool cuts::min_cc(int cc_segm, int cc_sect, int nphe){
	bool pass = false;
	int seg = detect::cc_segment(cc_segm);
	if(nphe > 35){//Just a flat cut for now
		pass = true;
	}
	return pass; 
}


bool cuts::min_ec(Float_t etot)
{
	bool pass_min_ec = kFALSE;
	if(etot >= ec_min_e16)
	{
		pass_min_ec = kTRUE;
	}
	return pass_min_ec;
}

float cuts::sf(Float_t etot, Float_t p){
	return (float)etot/(float)p;
}

float cuts::sf_low(Float_t p, int sidx, int r ){
	//If no run is given assume e16
	//1 = e16
	//2 = e1f
	/* haven't defined e1f parameters yet
	if(r = 2){
		return sf_low_e1f[sidx][0] + sf_low_e1f[sidx][1]*p + sf_low_e1f[sidx][2]*p*p + sf_low_e1f[sidx][3]*p*p*p;
	}
	else
	*/
	return sf_low_e16[sidx][0] + sf_low_e16[sidx][1]*p + sf_low_e16[sidx][2]*p*p + sf_low_e16[sidx][3]*p*p*p;
}

float cuts::sf_high(Float_t p, int sidx, int r){
	//If no run is given assume e16
	//1 = e16
	//2 = e1f
	/*haven't defined e1f parameters yet
	if(r = 2){
		return sf_high_e1f[sidx][0] + sf_high_e1f[sidx][1]*p + sf_high_e1f[sidx][2]*p*p + sf_high_e1f[sidx][3]*p*p*p;
	}
	else
		*/
	return sf_high_e16[sidx][0] + sf_high_e16[sidx][1]*p + sf_high_e16[sidx][2]*p*p + sf_high_e16[sidx][3]*p*p*p;
}
//Electron Sampling Fraction cut using the EC
bool cuts::sf_cut(Float_t p, Float_t etot, Float_t cx, Float_t cy, int run, bool sim)
{
	bool pass_sf = kFALSE;
	float s_f = sf(etot,p);
	int sector = physics::get_sector(physics::get_phi(cx,cy));
	int sidx = sector -1;
	float low, high;
	if(!sim){
		low = sf_low( p, sidx, run);
		high = sf_high(p, sidx, run);
	}else{
		low = 0.0;
		high = 0.18;
	}
	if(s_f >= low && s_f <= high)
	{
		pass_sf = kTRUE;
	}
	return pass_sf;
}

bool cuts::MM_cut(int top_, float MM){
	bool pass = false;
	float top = NAN;
	float bot = NAN;
	switch(top_){
		case 0:
			top = p_center + p_sig;
			bot = p_center - p_sig; 
		break;
		case 1:
			top = pip_center + pip_sig;
			bot = pip_center - pip_sig; 
		break;
		case 2:
			top = pim_center + pim_sig;
			bot = pim_center - pim_sig; 
		break;
		case 3:
			top = MM_zero_center + MM_zero_sigma;
			bot = MM_zero_center - MM_zero_sigma; 
		break;
	}
	if((MM > bot) && (MM < top)){
		pass = true;
	}
	return pass; 
}

//Putting together the e_sanity cuts based on environment we set
bool cuts::in_range(float W_, float Q2_, std::shared_ptr<Environment> envi){
	bool pass = false;
	if(W_ >= envi->was_Wmin() && W_ <= envi->was_Wmax() ){//Checking to see if the particle is in the relevant W Q2 region 
		if(Q2_ >= envi->was_Qmin() && Q2_ <= envi->was_Qmax()){
			pass = true;
		}
	}
	return pass; 
}

bool cuts::e_sanity(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int idx_){
  bool pass = false; 
  bool qcheck = false;
  bool dc = false;
  bool sc = false;
  bool ec = false;
  bool cc = false;
  if(data->Branches::q(idx_) == -1){
  	qcheck = true;
  }
  //DC Cut
  if(envi->was_dc_hit() && data->Branches::dc(idx_)){
    dc = true;
  }else if(!(envi->was_dc_hit())){
    dc = true;
  }
  if(envi->was_sc_hit() && data->Branches::sc(idx_)){
    sc = true;
  }else if(!(envi->was_sc_hit())){
    sc = true;
  }
  if(envi->was_ec_hit() && data->Branches::ec(idx_)){
    ec = true;
  }else if(!(envi->was_ec_hit())){
    ec = true;
  }
  if(envi->was_cc_hit() && data->Branches::cc(idx_) && !(envi->was_sim())){
    cc = true;
  }else if(!(envi->was_cc_hit()) || envi->was_sim()){
    cc = true;
  }
  pass = qcheck && dc && sc && ec && cc; 
  return pass; 
}
bool cuts::h_sanity(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int idx_, int par){
  bool pass = false;
  bool dc = false;
  bool charge = false;
  bool sc = false;
  //DC Cut
  if(envi->was_dc_hit() && data->Branches::dc(idx_)){
    dc = true;
  }else{
  	dc = false;
  }
  //SC
  if(envi->was_sc_hit() && data->Branches::sc(idx_)){
   	sc = true;
  } else{
  	sc = false;
  }
  //Charge Match
  switch(par){
  	case 0: 
  		if(data->Branches::q(idx_) == 1){
  			charge = true; 
  		}else{
  			charge = false;
  		}
  	break;
  	case 1: 
  		if(data->Branches::q(idx_) == 1){
  			charge = true; 
  		}else{
  			charge = false;
  		}
  	break;
  	case 2: 
  		if(data->Branches::q(idx_) == -1){
  			charge = true; 
  		}else{
  			charge = false;
  		}
  	break;
  }
  if(charge && dc && sc){
  	pass = true;
  }
  return pass; 
}
bool cuts::e_cc(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int idx_){
	bool cc = false;
	if(!(envi->was_sim())){
		if(envi->was_eid_cc() && cuts::min_cc(data->Branches::cc_segm(idx_),data->Branches::cc_sect(idx_),data->Branches::nphe(idx_))){
	    	cc = true;
	  	}
	}
  	return cc; 
}
bool cuts::e_ec(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int idx_){
	bool ec = false;
	if(envi->was_eid_ec() && cuts::min_ec(data->Branches::etot(idx_))){
    	ec = true;
  	}
  	return ec; 
}
bool cuts::e_sf(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int idx_){
	bool sf = false;
	if(envi->was_eid_sf() && cuts::sf_cut(data->Branches::p(idx_),data->Branches::etot(idx_),data->Branches::cx(idx_),data->Branches::cy(idx_),envi->Environment::was_data_set(),envi->Environment::was_sim())){
    	sf = true;
  	}
  	return sf; 
}
bool cuts::e_fid(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int idx_){
	bool fid = false;
	if(envi->was_eid_fid() && cuts::fid_cut(idx_,data->Branches::p(idx_),data->Branches::cx(idx_),data->Branches::cy(idx_),data->Branches::cz(idx_))){
    	fid = true;
  	}
  	return fid; 
}

bool cuts::h_fid(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int par, int had){
	bool pass = false;
	if(envi->was_hid_fid(had) && cuts::fid_cut(had+1,data->Branches::p(par),data->Branches::cx(par),data->Branches::cy(par),data->Branches::cz(par))){
    	pass = true;
  	}
  	return pass; 
}


bool cuts::h_dt(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int par, int had){
	bool pass = false;
	if(envi->was_hid_dt(had) && cuts::delta_t_cut(had+1, data->Branches::p(par), data->Branches::sc_r(0), data->Branches::sc_r(par), data->Branches::sc_t(0), data->Branches::sc_t(par))){
    	pass = true;
  	}
  	return pass; 
}

bool cuts::pim_e_sep(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int par, int had){
	bool pass = false; 
	if(envi->was_hid_e()){
		if(had == 2){
			if(data->cc(par)>0){//If it has a registerd hit in the CC 
				if(cuts::min_cc(data->Branches::cc_segm(par),data->Branches::cc_sect(par),data->Branches::nphe(par))){//Min CC Cut
					//if(cuts::sf_cut(data->Branches::p(par),data->Branches::etot(par),data->Branches::cx(par),data->Branches::cy(par))){//SF Cut
						//if(cuts::fid_cut(h,data->Branches::p(h),data->Branches::cx(h),data->Branches::cy(h),data->Branches::cz(h))){//Electron Fiducial //I don't know that this should be necessary because it doesn't work towards e-pim separation 
							pass = true;//Defined quantity for Event object  
						//}
					//}
				}
			}else{
				pass = true;
			}
		}else{
			pass = true;
		}
	}
	return pass; 
}

/*
bool cuts::elec_p_cut(int set, float W_, float Q2_, float p, int had){
	bool pass = false;
	float W_width = 0.05;
	float Q2_width = 0.05;
	float p_beam = NAN; 
	switch(set){
		case 1:
			p_beam = energy_e16;
		break;
		case 0:
			p_beam = energy_e1f;
		break;
	}
	if(had!=2){
		pass = true;
	}else{
		//W width
		if(p >= (-(W_-W_width)*(W_-W_width)+Q2_-mp*mp-2*mp*p_beam)/(2*mp) || p <= (-(W_+W_width)*(W_+W_width)+Q2_-mp*mp-2*mp*p_beam)/(2*mp)){
			//Q2 Width
			if(p >= (-(W_)*(W_-W_width)+(Q2_+Q2_width)-mp*mp-2*mp*p_beam)/(2*mp) || p <= (-(W_)*(W_)+(Q2_-Q2_width)-mp*mp-2*mp*p_beam)/(2*mp)){
				pass = true; 
			}
		}
	}
	return pass;
}*/


bool cuts::elec_p_cut(int set, std::shared_ptr<Branches> data, int part, int had){
	bool pass = false;
	float W_width = 0.05;
	float Q2_width = 0.05;
	float p_beam = NAN; 
	float W_ = physics::WP(set,data);
	float Q2_ = physics::Qsquared(set,data);
	float p = data->Branches::p(part);
	switch(set){
		case 1:
			p_beam = energy_e16;
		break;
		case 0:
			p_beam = energy_e1f;
		break;
	}
	if(had!=2){
		pass = true;
	}else{
		//W width
		if(p >= (-(W_-W_width)*(W_-W_width)-(Q2_-Q2_width)+mp*mp+2*mp*p_beam)/(2*mp) || p <= (-(W_+W_width)*(W_+W_width)-(Q2_+Q2_width)+mp*mp+2*mp*p_beam)/(2*mp)){
			pass == true;
		}
	}
	return pass;
}

bool cuts::e_dt(int set, std::shared_ptr<Branches> data, int idx){
	bool pass = false;
	float W_width = 0.05;
	float Q2_width = 0.05;
	float p_beam = NAN; 
	float W_ = physics::WP(set,data);
	float Q2_ = physics::Qsquared(set,data);
	float p = data->Branches::p(idx);
	switch(set){
		case 1:
			p_beam = energy_e16;
		break;
		case 0:
			p_beam = energy_e1f;
		break;
	}		
	//W width
	if(p <= (-(W_-W_width)*(W_-W_width)-(Q2_-Q2_width)+mp*mp+2*mp*p_beam)/(2*mp) && p >= (-(W_+W_width)*(W_+W_width)-(Q2_+Q2_width)+mp*mp+2*mp*p_beam)/(2*mp)){
		pass == true;
	}
	return pass;

}

bool cuts::eid(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int idx_){
	bool pass = true;
	if(envi->was_dc_hit()){
		if(data->Branches::dc(idx_)>0){
			pass &= true;
		}
	}
	if(envi->was_ec_hit()){
		if(data->Branches::ec(idx_)>0){
			pass &= true;
		}
	}
	if(envi->was_sc_hit()){
		if(data->Branches::sc(idx_)>0){
			pass &= true;
		}
	}
	if(envi->was_cc_hit()){
		if(data->Branches::cc(idx_)>0){
			pass &= true;
		}
	}
	if(envi->was_eid_fid()){
		pass &= cuts::e_fid(data,envi,idx_);
	}
	if(envi->was_eid_sf()){
		pass &= cuts::e_sf(data,envi,idx_);
	}
	if(envi->was_eid_ec()){
		pass &= cuts::e_ec(data,envi,idx_);
	}
	if(envi->was_eid_cc()){
		pass &= cuts::e_cc(data,envi,idx_);
	}
	return pass; 
}



bool cuts::hid(int set, std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int par, int had){
	bool pass = false;
	bool fid = false;
	bool san = false;
	bool dt = false;
	san = cuts::h_sanity(data,envi,par,had);
	if(envi->was_hid_fid(had)){
		fid = cuts::h_fid(data,envi,par,had);
	}
	if(envi->was_hid_dt(had)){
		dt = cuts::h_dt(data,envi,par,had);
	}
	if(envi->was_hid_e()){
		pass &= cuts::pim_e_sep(data,envi,par,had);
	}
	//pass &= cuts::elec_p_cut(set,data,par,had);
	if(san && fid && dt){
		pass = true;
	}
	return pass; 
}

bool cuts::p_miss(std::shared_ptr<Environment> envi){
	bool pass = false; 
	if(envi->was_top(0)){
		pass = true; 
	}
	return pass; 
}
bool cuts::pip_miss(std::shared_ptr<Environment> envi){
	bool pass = false; 
	if(envi->was_top(1)){
		pass = true; 
	}
	return pass;
}
bool cuts::pim_miss(std::shared_ptr<Environment> envi){
	bool pass = false; 
	if(envi->was_top(2)){
		pass = true; 
	}
	return pass;
}
bool cuts::z_miss(std::shared_ptr<Environment> envi){
	bool pass = false; 
	if(envi->was_top(3)){
		pass = true; 
	}
	return pass;
}

bool cuts::p_corr(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi){
	bool pass = false; 
	return pass; 
}
bool cuts::eff_cut(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi){
	bool pass = false; 
	return pass; 
}



