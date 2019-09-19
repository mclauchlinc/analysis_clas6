#ifndef EVENT_CLASS_H_GUARD
#define EVENT_CLASS_H_GUARD

#include "TLorentzVector.h"
#include <iostream>
#include "histogram.hpp"
#include "constants.hpp"
#include "branches.hpp"
#include "eid.hpp"//Electron ID Cuts
#include "cuts.hpp"
//#include "hid.hpp"//Hadron ID Cuts

//#include "ntuple.hpp"

/*
class Particle{
private:
	int _pid = 0; //final particle ID
	int _pidn = 0; //number of particles IDed as. Will need to be resolved if >1, and considered empty when 0

	TLorentzVector _4vec; //The four vector for this particle in lab frame

	bool _epass = false;
	bool _ppass = false;
	bool _pippass = false;
	bool _pimpass = false;

	float _px = 0.0;
	float _py = 0.0;
	float _pz = 0.0;

public:
	void SetPID(int i);

	int Get_PIDn(){
		return _pidn; 
	}

};*/

class Event_Class{
private:
	//int _gpart = 0; //Number of tagged particles in the Event_Class
	bool _valid = false; //Valid trigger electron from eid 
	int _top = 0; //Topology this fits under {bad, pmiss, pipmiss,pimmiss,zeromiss }->{0,1,2,3,4}

	//All four vectors in the center of mass frame 
	TLorentzVector _beam;
	TLorentzVector _elec;
	TLorentzVector _gamma;
	TLorentzVector _target;
	TLorentzVector _prot;
	TLorentzVector _pip;
	TLorentzVector _pim;

	int good_electron = 0; 
	int good_proton = 0;
	int good_pip = 0; 
	int good_pim = 0; 

	int check_idx[3]={-99,-99,-99};

	int pro_idx[12];
	int pip_idx[12];
	int pim_idx[12];

	float _alpha1 = -99;
	float _alpha2 = -99;
	float _alpha3 = -99;
	float _theta1 = -99;
	float _theta2 = -99;
	float _theta3 = -99;
	float _MMt1 = -99; //MM proton/pip
	float _MMt2 = -99;	//MM proton/pim
	float _MMt3 = -99;	//MM Pip/pim

	float _MM = -99; 
	float _MM2 = -99; 

	float _W = -99;
	float _Q2 = -99; 

	int _helicity = 0; 

	float MM_p = -99;
	float MM_p2 = -99;
	float MM_pip = -99;
	float MM_pip2 = -99;
	float MM_pim = -99;
	float MM_pim2 = -99;
	float MM_z = -99;
	float MM_z2 = -99;

	int ppip = 0; 

	float d[3] = {-99.0,-99.0,-99.0};
	float t[3] = {-99.0,-99.0,-99.0} ;
	float cx[3] = {-99.0,-99.0,-99.0} ;
	float cy[3] = {-99.0,-99.0,-99.0};
	float cz[3] = {-99.0,-99.0,-99.0};
	float _p[3] = {-99.0,-99.0,-99.0};
	int h_sec[3] = {-99,-99,-99};
	
public:
	Event_Class(std::shared_ptr<Branches> data, std::shared_ptr<Histogram> _hists, int run_type);
	//~Event_Class();

	float Get_px(int i);
	float Get_py(int i);
	float Get_pz(int i);
	float Get_p0(int i);
	float Get_hel();
	float Get_top();
	float Get_pid(int i);
	bool is_valid();

	int Get_ppip();
	//void Fill_Tree(forest tree, int event_n);
	/*
	void Assign_electron(float p, float cx, float cy, float cz);
	void Assign_proton(float p, float cx, float cy, float cz);
	void Assign_pip(float p, float cx, float cy, float cz);
	void Assign_pim(float p, float cx, float cy, float cz);
	*/


};

/*
//Overall reaction with information from all the given events
class Reaction{
private:
	size_t _total_yield = 0;//total number of measured events 
	size_t _top_yield[4] = {0,0,0,0};//yields for the various topologies

	/*
	What I want bins for
	alpha[3]
	theta[3]


	//

public:
	//Constructor that adds yields to the overall reaction 
	Reaction(std::unique_ptr<Event> event);
	size_t Get_total_yield();
	size_t Get_top_yield(int i);

};*/


#endif

