#ifndef EVENT_H_GUARD
#define EVENT_H_GUARD

#include "particle.hpp"
#include "cuts.hpp"
#include "histogram.hpp"
#include "physics.hpp"
#include "constants.hpp"
#include "branches.hpp"
#include "event_selection.hpp"

class Event{
private:
	//Event Information
	bool _valid = false; 
	bool _sim = false;
	int topology = 0; //{0,1,2,3,4}->{none,pMiss,pipMiss,pimMiss,All Detected}
	float _W = NAN;
	float _Q2 = NAN;
	float _MM[4] = {NAN,NAN,NAN,NAN};//For the four topologies
	float _MM2[4] = {NAN,NAN,NAN,NAN}; //For the four topologies

	Particle _elec;
	Particle _pro;
	Particle _pip; 
	Particle _pim; 

	TLorentzVector _k_mu = {NAN,NAN,NAN,NAN};
	TLorentzVector _p_mu = p_mu;

	int _meas_had[4] = {0,0,0,0}; //{pro,pip,pim,pro/pip}


public:
	Event(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_, int data_set_, bool sim_ = false);
	~Event();


};



#endif