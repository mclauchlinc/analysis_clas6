#ifndef PHYSICS_H_GUARD
#define PHYSICS_H_GUARD

#include "TMath.h"
#include <TLorentzVector.h>
#include "TVector3.h"
#include <cmath>
#include "constants.hpp"
#include "branches.hpp"

namespace physics{
	void Print_4Vec(TLorentzVector k1);
	void Print_3Vec(TVector3 p1);
	bool Check_4Vec(TLorentzVector k1);
	TLorentzVector Make_4Vector(float p, float cx, float cy, float cz, float m);
	TLorentzVector Make_4Vector(float px, float py, float pz, float m);
	TLorentzVector Make_4Vector(bool here, float p, float theta, float phi, float m);
	TLorentzVector Set_k_mu(int set);
	int event_helicity(std::shared_ptr<Branches> data, int plate_stat);
	float Qsquared(int set, std::shared_ptr<Branches> data, int thr = 0);
	float Qsquared(int set, std::shared_ptr<Branches> data, bool thrown_);
	float WP(int set, std::shared_ptr<Branches> data, int thr = 0);

	float beta_calc(float m, std::shared_ptr<Branches> data, int i);
	float MM_event(int set, int squared, TLorentzVector k1_mu, TLorentzVector k2_mu, TLorentzVector k3_mu, TLorentzVector k4_mu={0.0,0.0,0.0,0.0});
	float MM_event(int squared, TLorentzVector k0_mu, TLorentzVector k1_mu, TLorentzVector k2_mu, TLorentzVector k3_mu, TLorentzVector k4_mu);
	float MM_event(int squared, TLorentzVector k0_mu, TLorentzVector k1_mu, TLorentzVector k2_mu, TLorentzVector k3_mu);
	float get_theta(float cz_);//in lab frame
	float get_theta(int part, std::shared_ptr<Branches> data);
	float get_phi(float cx_, float cy_);//in lab frame
	float get_phi(int part, std::shared_ptr<Branches> data);//in lab frame
	float get_phi_pos(float cx_, float cy_);//in lab frame, but going to 360 rather than -180->180. Keeping sector 1 centered on 0 degrees. 
	int get_sector(float phi_);//Phi must be from lab frame
	float phi_center(float phi_);//Center phi within the sector

	//Delta T
	float vert_e(float d, float t);
	float vert_h(float p, float d, float t, float m);
	float delta_t(int part, float p, float d, float t, float d0, float t0);
	float delta_t(int part, std::shared_ptr<Branches> data, int idx);


	//Math
	TVector3 V4_to_V3(TLorentzVector p1);//Get just the three vector part out of a four vector
	float Vec3_Mag(TVector3 v1);
	float Vec3_Mag(TLorentzVector p1);
	TVector3 Cross_Product(TLorentzVector p1, TLorentzVector p2);//Get the cross product between two vectors
	TVector3 Cross_Product(TVector3 v1, TVector3 v2); 
	float Dot_Product(TVector3 v1, TVector3 v2);//Get the dot product between two vectors
	float Dot_Product(TLorentzVector p1, TLorentzVector p2);
	float Cos_Vecs(TVector3 v1, TVector3 v2); //Get the Cosine between two vectors 
	float Cos_Vecs(TLorentzVector p1, TLorentzVector p2);
	float Sin_Vecs(TVector3 v1, TVector3 v2);
	float Sin_Vecs(TLorentzVector p1, TLorentzVector p2); //Get the Sin between two vectors 
	float Get_phie(int set, TLorentzVector p0);
	float Get_phie(TLorentzVector k0, TLorentzVector p0);
	void Rotate_4Vec(int set, float theta, float phi, float phie, TLorentzVector &p1); //Rotate Four vectors along the theta and phi angles
	void Rotate_4Vec(float theta, float phi, float phie, TLorentzVector &p1);
	void Boost_4Vec(float beta, TLorentzVector &p1 );// Boost a four vector in the z direction 
	void COM_gp(int set, TLorentzVector &p0, TLorentzVector &p1, TLorentzVector &p2, TLorentzVector &p3); //Bring four vectors into the COM reference frame for excited nucleon 
	void COM_gp(TLorentzVector &k0, TLorentzVector &p0, TLorentzVector &p1, TLorentzVector &p2, TLorentzVector &p3);
	TLorentzVector COM_gp(int par, TLorentzVector k0, TLorentzVector p0, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3);
	float alpha(int top, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, int set);
	float alpha(int top, TLorentzVector k0, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, bool COM); //Alpha angle between scattering planes
	float epsilon(int set, float Energy, float Q_2); //Virtual photon transverse polarization
	float epsilon(TLorentzVector k0, float Energy, float Q_2); 
	float MM_2(TLorentzVector p1, TLorentzVector p2);//Get the MM of a two particle system
	float gamma_nu(int set, float Ep, float Q_2, float W_);//Virtual Photon Flux
	int Qfaraday(int q_last, int q_next, int q_tot);//Faraday Cup counting 
	float Ev_Theta(int top, TLorentzVector k0, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, bool COM);
	float Ev_MM(int top, TLorentzVector k0, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, bool COM);


	


}

#endif