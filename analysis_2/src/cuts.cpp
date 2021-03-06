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

bool cuts::delta_t_cut(int part, float p, float d0, float d, float t0, float t) //Note: d and t need the sc_index 
{
	bool pass = false;

	float dt = physics::delta_t(part, p, d, t, d0, t0);
	//std::cout<<std::endl <<"delta t cut: p: " <<p <<" dt: " <<dt <<" cut_low: " <<cuts::dt_p_low(p) <<" cut_high:" <<cuts::dt_p_high(p) <<std::endl;  
	if(dt>cuts::dt_p_low(p) && dt<cuts::dt_p_high(p) )
	{
		pass = true;
	}
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
bool cuts::sf_cut(Float_t p, Float_t etot, Float_t cx, Float_t cy, int r)
{
	bool pass_sf = kFALSE;
	float s_f = sf(etot,p);
	int sector = physics::get_sector(physics::get_phi(cx,cy));
	int sidx = sector -1;
	float low = sf_low( p, sidx, r);
	float high = sf_high(p, sidx, r);
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
			top = MM_zero_center2 + MM_zero_sigma2;
			bot = MM_zero_center2 - MM_zero_sigma2; 
		break;
	}
	if((MM > bot) && (MM < top)){
		pass = true;
	}
	return pass; 
}

