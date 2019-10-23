#ifndef EVENT_H_GUARD
#define EVENT_H_GUARD

#include "TLorentzVector.h"
#include <iostream>
#include "histogram.hpp"
#include "constants.hpp"
#include "branches.hpp"
#include "eid.hpp"//Electron ID Cuts
#include "cuts.hpp"
//#include "particle.hpp"
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

class Event{
private:
	int _top = 0; //Topology this fits under {bad, pmiss, pipmiss,pimmiss,zeromiss }->{0,1,2,3,4}

	//All four vectors in the center of mass frame 
	TLorentzVector _beam = {NAN,NAN,NAN,NAN};
	TLorentzVector _elec= {NAN,NAN,NAN,NAN};
	TLorentzVector _gamma= {NAN,NAN,NAN,NAN};
	TLorentzVector _target= {NAN,NAN,NAN,NAN};
	TLorentzVector _pro= {NAN,NAN,NAN,NAN};
	TLorentzVector _pip= {NAN,NAN,NAN,NAN};
	TLorentzVector _pim= {NAN,NAN,NAN,NAN};

	int _fc_curr = 0; 

	float _alpha1 = NAN;	//{}
	float _alpha2 = NAN;	//{}
	float _alpha3 = NAN;	//{}
	float _theta1 = NAN; 	//proton
	float _theta2 = NAN;	//pip
	float _theta3 = NAN;	//pim
	float _MMt1 = NAN; //MM proton/pip
	float _MMt2 = NAN;	//MM proton/pim
	float _MMt3 = NAN;	//MM Pip/pim

	float _W = NAN;
	float _Q2 = NAN; 

	int _helicity = 0; 

	
public:
	Event(std::shared_ptr<Branches> data, std::shared_ptr<Histogram> _hists, int run_type, int plate_info, int data_set = 0);//default to e16
	//Run type reffers to simulation vs data
	//~Event();

	float Get_px(int i);
	float Get_py(int i);
	float Get_pz(int i);
	float Get_p0(int i);
	float Get_hel();
	float Get_top();
	float Get_pid(int i);
	bool is_valid();

	int Get_ppip(int idx);
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

