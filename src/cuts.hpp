#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD

#include "functions.hpp"
#include "physics.hpp"
#include "detectors.hpp"

namespace cuts{
	//Fiducial Cuts
	bool fid_e (float p, float cx, float cy, float cz);
	float phi_min(float theta, int idx);
	float phi_max(float theta, int idx);
	bool fid_h ( float p, float cx, float cy, float cz);
	bool fid_cut( int part, float p, float cx, float cy, float cz);

	//DT Cuts
	//Proton Formulae from Arjun
	float dt_p_low(float p);
	float dt_p_high(float p);
	bool delta_t_cut(int part, float p, float d0, float d, float t0, float t);
	
	//CC Cuts
	bool min_cc(int cc_segm, int cc_sect, int nphe);
	
	//Sampling Fraction Cuts
	bool min_ec(Float_t etot);
	float sf(Float_t etot, Float_t p);
	float sf_low(Float_t p, int sidx, int r =1);
	float sf_high(Float_t p, int sidx, int r=1);
	bool sf_cut(Float_t p, Float_t etot, Float_t cx, Float_t cy, int r = 1);


}


#endif