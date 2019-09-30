#ifndef PARTICLE_H_GUARD
#define PARTICLE_H_GUARD

#include "physics.hpp"
#include "Branches.hpp"
#include "cuts.hpp"
#include "event_class.hpp"

class Particle: public Event_Class{
private:
	bool hadron = false; 

	int q = 0; 
	int idx = -99;
	int nphe = -99;
	int id = 0; 
	int sector = 0; 
	int cc_segm = NAN; 
	int cc_sect = NAN

	float sc_r = NAN; //Distance through SC 
	float sc_r0 = NAN; //Distance through SC for tagged electron
	float sc_t = NAN;	//Time through SC 
	float sc_t0 = NAN; //TIme through SC for tagged electron
	float dt_t[3] = NAN; //Delta T for all particle types
	float dt = NAN; //Assigned post particle determination 
	float cx = NAN; 
	float cy = NAN; 
	float cz = NAN; 
	float p = NAN;
	float theta = NAN; 
	float phi = NAN;
	float ec_energy = NAN; 
	float sf = NAN; 
	

	TLorentzVector k1;  

	bool pim_ele = false; //Double ID
	bool pro_pip = false; //Double ID

	bool sanity_pass[4] = false; 
	bool min_ec_pass = false;
	bool fid_pass[4] = false; 
	bool sf_pass = false;
	bool cc_pass = false; 
	bool dt_pass[3] = false; 
	bool double_id_resolved = false;
	bool pid_pass[4] = {false,false,false,false};


public:
	Particle(bool is_hadron);

	void EID_Cuts();
	void Fill_Particle(std::shared_ptr<Branches> data);
}


#endif