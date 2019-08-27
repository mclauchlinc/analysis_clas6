#ifndef EVENT_CLASS_H_GUARD
#define EVENT_CLASS_H_GUARD

#include "TLorentzVector.h"
#include <iostream>
#include "Branches.hpp"
#include "histogram.hpp"



class Particle{
private:
	int _pid = 0; //final particle ID
	int _pidn = 0; //number of particles IDed as. Will need to be resolved if >1, and considered empty when 0

	TLorentzVector _4vec; //The four vector for this particle in lab frame

	bool _epass = false;
	bool _ppass = false;
	bool _pippass = false;
	bool _pimpass = false;

public:
	void SetPID(int i);

	int Get_PIDn(){
		return _pidn; 
	}

};

class Event{
private:
	//int _gpart = 0; //Number of tagged particles in the event
	bool _valid = 0; //Valid trigger electron from eid 
	int _top = 0; //Topology this fits under {bad, pmiss, pipmiss,pimmiss,zeromiss }->{0,1,2,3,4}

	//All four vectors in the center of mass frame 
	TLorentzVector _beam;
	TLorentzVector _elec;
	TLorentzVector _gamma;
	TLorentzVector _target;
	TLorentzVector _prot;
	TLorentzVector _pip;
	TLorentzVector _pim;

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
	
public:
	Event(std::shared_ptr<Branches> data, std::shared_ptr<Histogram> hists, int run_type);

};


#endif

