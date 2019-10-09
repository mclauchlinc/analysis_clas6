#ifndef PARTICLE_H_GUARD
#define PARTICLE_H_GUARD

#include "physics.hpp"
#include "Branches.hpp"
#include "cuts.hpp"
#include "event_class.hpp"
#include "histogram.hpp"

class Particle{
private:
	bool hadron = false; 

	int q = 0; 
	int idx = -99;
	int nphe = -99;
	int id = 0; 
	int sector = 0; 
	int cc_segm = NAN; 
	int cc_sect = NAN;
	int bank = -99;//bank id

	float sc_r = NAN; //Distance through SC 
	float sc_r0 = NAN; //Distance through SC for tagged electron
	float sc_t = NAN;	//Time through SC 
	float sc_t0 = NAN; //TIme through SC for tagged electron
	float dt_t[3] = {NAN,NAN,NAN}; //Delta T for all particle types
	float dt = NAN; //Assigned post particle determination 
	float cx = NAN; 
	float cy = NAN; 
	float cz = NAN; 
	float p = NAN;
	float theta = NAN; 
	float phi = NAN;
	float etot = NAN; 
	float sf = NAN; 

	int dc = 0;
	int sc = 0; 
	int ec = 0;
	int cc = 0; 
	

	TLorentzVector k1 = {NAN,NAN,NAN,NAN};  

	bool pim_ele = false; //Double ID
	bool pro_pip = false; //Double ID

	bool sanity_pass[4] = {false,false,false,false}; 
	bool min_ec_pass = false;
	bool fid_pass[4] = {false,false,false,false}; 
	bool sf_pass = false;
	bool cc_pass = false; 
	bool dt_pass[3] = {false,false,false}; 
	bool double_id_resolved = false;
	bool pid[4] = {false,false,false,false};


public:
	Particle();

	void Fill_Particle(std::shared_ptr<Branches> data, int par_idx);

	void EID_Cuts();//std::shared_ptr<Histogram> hist_);
	void Ele_Pim_Resolve();
	void HID_Cuts();//std::shared_ptr<Histogram> hist_);

	bool Is_Elec();
	bool Is_Pro();
	bool Is_Pip();
	bool Is_Pim();


	void Pro_Pip_Resolved(int par);

	bool Is_fid_pass(int par);
	bool Is_cc_pass();
	bool Is_min_ec_pass();
	bool Is_dt_pass(int par);
	bool Is_sanity_pass(int par);
	bool Is_sf_pass();

	bool Pro_Pip_Prob();

	void Fill_Fid();

	void Bank();

	int par_idx();
	int par_nphe();
	int par_id();
	int par_sector();
	int par_cc_sector();
	int par_cc_segm();

	float par_dt(int par);
	float par_theta();
	float par_phi();
	float par_etot();
	float par_p();

	TLorentzVector Par_4Vec();
};


#endif