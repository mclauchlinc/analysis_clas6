#include "cuts.hpp"


double cuts::phi_min(double theta, int idx){
	return -(a0mh[idx]*(1.0-TMath::Exp(-a1mh[idx]*(theta-a2mh[idx])))-a3mh[idx]);
}

double cuts::phi_max(double theta, int idx){
	return (a0xh[idx]*(1.0-TMath::Exp(-a1xh[idx]*(theta-a2xh[idx])))+a3xh[idx]);
}

bool cuts::fid_h ( double p, double cx, double cy, double cz)
{
	bool in_fid = kFALSE;

	//degree Conversion
	double degree = 180.0/TMath::Pi();

	//Calculate angles and sector
	//Phi is centered for the sector
	double theta = physics::get_theta(cz);
	double phi_c = physics::phi_center(physics::get_phi(cx, cy));
	int sector = physics::get_sector(physics::get_phi(cx, cy));

	int sec_indx = sector -1;

	//Actual application of the cut
	if(phi_c>=cuts::phi_min(theta , sec_indx) && phi_c<=cuts::phi_max(theta, sec_indx))
	{
		in_fid = kTRUE;
	}

	return in_fid;
}

//Proton Formulae from Arjun
double cuts::dt_p_low(double p){
	return DTL[0]+DTL[1]*p+DTL[2]*p*p+DTL[3]*p*p*p;
}

double cuts::dt_p_high(double p){
	return DTH[0]+DTH[1]*p+DTH[2]*p*p+DTH[3]*p*p*p;
}

bool cuts::delta_t(int part, double p, double d0, double d, double t0, double t) //Note: d and t need the sc_index 
{
	bool pass = false;
	
	double dt = physics::delta_t(part, p, d, d0, t, t0);
	if(dt>cuts::dt_p_low(p) && dt<cuts::dt_p_high(p) )
	{
		pass = true;
	}
	return false;
}

bool cuts::min_cc(int cc_segm, int cc_sect, int nphe){
	bool pass = false;
	int seg = detect::cc_segment(cc_segm);
	if(nphe > 35){//Just a flat cut for now
		pass = true;
	}
	return pass; 
}