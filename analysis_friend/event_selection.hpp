#ifndef EVENT_SELECTION_H_GUARD
#define EVENT_SELECTION_H_GUARD

//#include "histogram.hpp"
//#include "event_class.hpp"
#include "physics.hpp"
//#include "functions.hpp"
#include "particle.hpp"
#include "cuts.hpp"

namespace Selection{
	//bool Event_Selection(TLorentzVector k_mu_, Particle elec_, Particle h1_, Particle h2_);//Missing Particle
	//bool Event_Selection(TLorentzVector k_mu_, Particle elec_, Particle h1_, Particle h2_, Particle h3_);//All Detected

	bool Event_Selection( int top_, TLorentzVector k_mu_, TLorentzVector elec_, TLorentzVector h1_, TLorentzVector h2_);
	bool Event_Selection( int top_, TLorentzVector k_mu_, TLorentzVector elec_, TLorentzVector h1_, TLorentzVector h2_, TLorentzVector h3_);
	bool Event_Selection( int top_, float MM_);
	//bool Topology(bool top[4], std::shared_ptr<Histogram> hist_, std::shared_ptr<Branches> data, TLorentzVector k0, TLorentzVector k1, TLorentzVector k2, TLorentzVector k3, int idx1, int idx2, int idx3);
}




#endif