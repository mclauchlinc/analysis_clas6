#ifndef EVENT_ANALYSIS_H
#define EVENT_ANALYSIS_H

#include "TLorentzVector.h"
#include <iostream>
#include "histogram.hpp"
#include "constants.hpp"
#include "branches.hpp"
#include "cuts.hpp"
#include "physics.hpp"
#include "ntuple.hpp"
#include "environment.hpp"
#include "functions.hpp"
#include "particle.hpp"
#include "event.hpp"

class Analysis{
private:
	float _gevts = 0.0; //Number of good events
	int _ntop[4] = {0,0,0,0}; //Number of good events for each topology
	int _mtop[4] = {0,0,0,0}; //Number of candidate events for each topology
	int _npart = 0; //Number of particles measured
	int _gpart[4] = {0,0,0,0}; //Number of identified particles of each species

	int _gElec[20] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	int _gPro[20] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	int _gPip[20] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	int _gPim[20] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

	float _W[2] = {NAN,NAN};
	float _Q2[2] = {NAN,NAN};

	int _set = -1; //{0,1} -> {e16,e1f}
	bool _sim = false;
	bool _thrown = false;
	bool _recon = false;




public:
	Analysis(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Environment> envi_, int run_type_);


	int Gevts();
	int Ntop(int i);

	float Get_W(int i);
	float Get_Q2(int i);
	int Get_set();
	bool Is_Sim();
	bool Is_Thrown();



};


#endif