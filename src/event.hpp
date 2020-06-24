#ifndef EVENT_H_GUARD
#define EVENT_H_GUARD

#include "histogram.hpp"
#include "particle.hpp"
#include "cuts.hpp"
#include "physics.hpp"
#include "constants.hpp"
#include "branches.hpp"
#include "event_selection.hpp"


class Event{
private:
	//Topology
	bool _top[4] = {false,false,false,false}; //Which topology was looked at {pmiss,pipmiss,pimmiss,zeromiss}
	bool _pass = false; //Did it pass this topology's MM cut? 
	float _weight = NAN; 
	bool _sim = false; 
	bool _thrown;
	bool _filled_correctly = false; //Was this event filled correctly?
	int _hel = 0;
	bool _COM = false;
	int _set = -1; 
	float _MM = NAN;
	float _MM2 = NAN;

	bool _part[4] = {false,false,false,false}; //Which particles are present?

	//Particle Attributes {ele',pro',pip',pim'}
	float _p_lab[4] = {NAN,NAN,NAN,NAN};//Particle momentum in lab frame
	float _theta_lab[4] = {NAN,NAN,NAN,NAN};//Particle theta  in lab frame
	float _phi_lab[4] = {NAN,NAN,NAN,NAN}; //Particle phi in lab frame
	float _p[4] = {NAN,NAN,NAN,NAN};//Particle momentum in COM frame
	float _theta[4] = {NAN,NAN,NAN,NAN};//Particle theta in COM frame
	float _phi[4] = {NAN,NAN,NAN,NAN}; //Particle phi in COM frame

	TLorentzVector _k1_lab = physics::Make_4Vector(energy_e16,0.0,0.0,1.0,me);//{NAN,NAN,NAN,NAN};
	TLorentzVector _k1 = {NAN,NAN,NAN,NAN};
	TLorentzVector _p1 = {NAN,NAN,NAN,NAN};
	TLorentzVector _vec_lab[4] = {{NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN}};//4 vectors in lab frame
	TLorentzVector _vec[4] = {{NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN}};//4 vectors in COM frame

	//Event Binning
	float _W = NAN; 
	float _Q2 = NAN; 
	float _MMb[3] = {NAN,NAN,NAN};//Combined missing mass for given expression of data 
	float _thetab[3] = {NAN,NAN,NAN};
	float _alphab[3] = {NAN,NAN,NAN};


public:
	Event();
	void COM_4Vec();
	void Vars();
	void Fill_Event(std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, int top_, float W_, float Q2_,  Particle p1,  Particle p2,  Particle p3, int hel_ = 1);
	void Fill_Event(std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, int top_, float W_, float Q2_,  Particle p1,  Particle p2,  Particle p3,  Particle p4, int hel_ = 1);
	void Fill_Event_Hists(std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, bool fit_);
	
	void Assign_Weight(float weight_);

	bool Top(int i);
	int Get_Top();
	bool Gevnt();
	float Get_Weight();
	bool Part_Present(int i);
	float Get_P(int i, bool COM_ = false);
	float Get_Theta(int i, bool COM_ = false);
	float Get_Phi(int i, bool COM_ = false); 
	float Get_Px(int i, bool COM_ = false);
	float Get_Py(int i, bool COM_ = false);
	float Get_Pz(int i, bool COM_ = false);
	float Get_P0(int i, bool COM_ = false);
	TLorentzVector Get_Beam(bool COM_ = false);
	TLorentzVector Get_Target(bool COM_ = false);
	float Get_Beam_Comp(int i, bool COM_ = false);
	float Get_Target_Comp(int i, bool COM_ = false);
	TLorentzVector Get_4Vec(int i, bool COM_ = false);
	float Get_W();
	float Get_Q2();
	float Get_MMb(int i);
	float Get_Thetab(int i);
	float Get_Alphab(int i);
	bool Get_COM();
	int Get_PID(int i);
	int Get_Set();
	int Get_Hel();

	//~Event();
};



#endif