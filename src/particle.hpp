#ifndef PARTICLE_H_GUARD
#define PARTICLE_H_GUARD

#include "physics.hpp"
#include "Branches.hpp"
#include "cuts.hpp"
//#include "event_class.hpp"
#include "histogram.hpp"
#include "physics.hpp"

class Particle{
private:
	int _idx = -1; 
	int _set = -1; //{1,0} -> {e16,e1f}
	bool _sim = false;
	bool _thrown = false;
	float _weight = NAN;

	float _p = NAN;//In lab frame
	int _q = 0; 
	float _dt[4] = {NAN,NAN,NAN,NAN};
	float _theta = NAN;//In lab frame
	float _phi = NAN;//In lab frame
	float _sf = NAN;//In lab frame
	float _etot = NAN;//Energy deposited in EC
	int _cc_seg = -1; //Segment of CC hit
	int _nphe = -1;//number photo electrons generated in CC

	bool _sanity_pass[4] = {false,false,false,false}; 
	bool _min_ec_pass = false;
	bool _fid_pass[4] = {false,false,false,false}; 
	bool _sf_pass = false;
	bool _cc_pass = false; 
	bool _dt_pass[4] = {false,false,false,false}; 
	bool _sp_dt_pass[3] = {false,false, false};
	bool _p_corr = false; 
	int _id_crisis = 0;//{0,1,2} = {none, pro/pip, pim/e}
	bool _pid[4] = {false,false,false,false};
	bool _ided = false; 

	float _x[3] = {NAN,NAN,NAN};//cc,sc,ec 
	float _y[3] = {NAN,NAN,NAN};//cc,sc,ec
	float _dtheta[3] = {NAN,NAN,NAN};
	float _dphi[3] = {NAN,NAN,NAN};


public:
	Particle();
	void Fill_Particle(std::shared_ptr<Branches> data_, int par_idx_, int set_, bool sim_, float _weight, bool thrown_ = false);


	void EID(std::shared_ptr<Branches> data_, std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, float W_);//Electron ID
	void HID(std::shared_ptr<Branches> data_, std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, float W_);//Hadron ID
	void PID(std::shared_ptr<Branches> data_, std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, float W_);//Acutal particle ID cut

	//Pass checks for all relevant pieces
	bool Pass_Sanity(int i);//i = {0,1,2,3} -> {e,pro,pip,pim}
	bool Pass_ec();
	bool Pass_fid(int i);//i = {0,1,2,3} -> {e,pro,pip,pim}
	bool Pass_sf();
	bool Pass_cc();
	bool Pass_dt(int i);//i = {0,1,2,3} -> {e,pro,pip,pim}
	bool Corr_p();
	int ID_crisis();
	bool IDed();//Has the particle been identified as a particle of interest?
	bool Is_Sim();
	bool Is_Thrown();

	bool Is_Elec();
	bool Is_Pro();
	bool Is_Pip();
	bool Is_Pim();

	float Get_p();
	float Get_theta();
	float Get_phi();

	int Get_set();
	int Get_idx();

	void Fill_Par_Event(std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, float W_, int top_, int par_, bool pass_);

	void Check_Particle();
	
};


#endif