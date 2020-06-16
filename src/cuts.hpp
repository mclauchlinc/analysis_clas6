#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD

#include "functions.hpp"
#include "physics.hpp"
#include "detectors.hpp"
#include "branches.hpp"
#include "environment.hpp"
#include "TMath.h"

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
	bool delta_t_cut_iso(int part, int part_iso, float p, float d0, float d, float t0, float t);
	bool delta_t_cut(int part, float p, float d0, float d, float t0, float t);
	bool delta_t_cut(int species_, float p_, float dt_);

	//CC Cuts
	bool min_cc(int cc_segm, int cc_sect, int nphe);
	
	//Sampling Fraction Cuts
	bool min_ec(Float_t etot);
	float sf(Float_t etot, Float_t p);
	float sf_low(Float_t p, int sidx, int r =1);
	float sf_high(Float_t p, int sidx, int r=1);
	bool sf_cut(Float_t p, Float_t etot, Float_t cx, Float_t cy, int r = 1);

	//Missing Mass Cuts
	bool MM_cut(int top, float MM);

	//Environment cuts
	//Environment things
	bool in_range(float W_, float Q2_, std::shared_ptr<Environment> envi);
	bool e_sanity(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int idx_);
	bool h_sanity(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int idx_, int par);
	bool e_cc(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int idx_);
	bool e_ec(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int idx_);
	bool e_sf(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int idx_);
	bool e_fid(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int idx_);
	//bool e_dt(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi);
	bool h_fid(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int par, int had);
	bool h_dt(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int par, int had);
	bool pim_e_sep(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int par, int had);
	bool elec_p_cut(int set, std::shared_ptr<Branches>, int part, int had);
	bool e_dt(int set, std::shared_ptr<Branches>, int part);
	bool eid(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int idx_);
	bool hid(int set, std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi, int par, int had);

	bool p_miss(std::shared_ptr<Environment> envi);
	bool pip_miss(std::shared_ptr<Environment> envi);
	bool pim_miss(std::shared_ptr<Environment> envi);
	bool z_miss(std::shared_ptr<Environment> envi);

	bool p_corr(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi);
	bool eff_cut(std::shared_ptr<Branches> data, std::shared_ptr<Environment> envi);



}


#endif