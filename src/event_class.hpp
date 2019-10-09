#ifndef EVENT_CLASS_H_GUARD
#define EVENT_CLASS_H_GUARD

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

class Event_Class{
private:
	//int _gpart = 0; //Number of tagged particles in the Event_Class
	bool _valid = false; //Valid trigger electron from eid 
	int _top = 0; //Topology this fits under {bad, pmiss, pipmiss,pimmiss,zeromiss }->{0,1,2,3,4}
	bool _assigned_4vecs = false; 

	//All four vectors in the center of mass frame 
	TLorentzVector _beam;
	TLorentzVector _elec;
	TLorentzVector _gamma;
	TLorentzVector _target;
	TLorentzVector _pro;
	TLorentzVector _pip;
	TLorentzVector _pim;

	int good_electron = 0; 
	int good_pro = 0;
	int good_pip = 0; 
	int good_pim = 0; 

	bool pim_is_ele = false; 

	int check_idx[3] = {-99,-99,-99};

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


	//Variables for multiple particles being Identified
	float dpro[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	float tpro[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN} ;
	float cxpro[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN} ;
	float cypro[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	float czpro[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	float ppro[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	int h_secpro[20] = {-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99};
	int pro_idx[20] = {-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99};

	float dpip[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	float tpip[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	float cxpip[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	float cypip[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	float czpip[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	float ppip[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	int h_secpip[20] = {-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99};
	int pip_idx[20] = {-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99};

	float dpim[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	float tpim[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	float cxpim[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	float cypim[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	float czpim[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	float ppim[20] = {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	int h_secpim[20] = {-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99};
	int pim_idx[20] = {-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99};

	float d[3]= {NAN,NAN,NAN}; 
	float t[3]= {NAN,NAN,NAN};
	float cx[3]= {NAN,NAN,NAN}; 
	float cy[3]= {NAN,NAN,NAN};
	float cz[3]= {NAN,NAN,NAN};
	float _p[3]= {NAN,NAN,NAN};
	int h_sec[3]= {-99,-99,-99};

	bool top_possible[4]= {false,false,false,false};
	
public:
	Event_Class(std::shared_ptr<Branches> data, std::shared_ptr<Histogram> _hists, int run_type, int data_set = 0);
	//~Event_Class();

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

