#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD

#include "functions.hpp"
#include "physics.hpp"
#include "detectors.hpp"

namespace cuts{
	//Fiducial Cuts
	double phi_min(double theta, int idx);
	double phi_max(double theta, int idx);
	bool fid_h ( double p, double cx, double cy, double cz);

	//DT Cuts
	//Proton Formulae from Arjun
	double dt_p_low(double p);
	double dt_p_high(double p);
	bool delta_t(int part, double p, double d0, double d, double t0, double t);

	//CC Cuts
	bool min_cc(int cc_segm, int cc_sect, int nphe);
	


}


#endif