#ifndef EVENT_SELECTION_H_GUARD
#define EVENT_SELECTION_H_GUARD

#include "histogram.hpp"
#include "event_class.hpp"
#include "physics.hpp"
#include "functions.hpp"

namespace Selection{
	bool Topology(bool top[4], std::shared_ptr<Histogram> hist_, std::shared_ptr<Branches> data, TLorentzVector k0, TLorentzVector k1, TLorentzVector k2, TLorentzVector k3, int idx1, int idx2, int idx3);
}




#endif