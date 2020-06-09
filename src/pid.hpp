#ifndef PID_H
#define PID_H

#include "cuts.hpp"
#include "constants.hpp"
#include "branches.hpp"
#include "histograms.hpp"


class Pid{
private:
	int _pid_n = 5 //Number of identified particles allowed per event
	//Particle Characteristics
	//Fermion
	float _p_e = NAN; 
	float _theta_e = NAN;
	int 

	//Hadron
	float _p_pro[_pid_n];
	float _p_pip[_pid_n];
	float _p_pim[_pid_n];
	float _theta_pro[_pid_n];
	float _theta_pip[_pid_n];
	float _theta_pim[_pid_n];
	float _phi_pro[_pid_n];
	float _phi_pip[_pid_n];
	float _phi_pim[_pid_n];
	float _dt_pro[_pid_n];
	float _dt_pip[_pid_n];
	float _dt_pim[_pid_n];
	int _idx_pro[_pid_n];
	int _idx_pip[_pid_n]; 
	int _idx_pim[_pid_n];

	int _n_pro = 0;
	int _n_pip = 0;
	int _n_pim = 0; 


public:
	PID();
}



#endif