#ifndef PHYSICS_H_GUARD
#define PHYSICS_H_GUARD

#include "TMath.h"
#include <TLorentzVector.h>
#include "TVector3.h"
#include <cmath>
#include "constants.hpp"
#include "branches.hpp"

namespace physics{
	TLorentzVector Make_4Vector(float p, float cx, float cy, float cz, float m);
	TLorentzVector Set_k_mu(int set);
	int event_helicity(std::shared_ptr<Branches> data, int plate_stat);
	float Qsquared(int set, std::shared_ptr<Branches> data);
	float WP(int set, std::shared_ptr<Branches> data);

	float beta_calc(float m, std::shared_ptr<Branches> data, int i);
	float MM_event(int set, TLorentzVector k1_mu, TLorentzVector k2_mu, TLorentzVector k3_mu, TLorentzVector k4_mu , int squared = 0);
	float get_theta(float cz_);//in lab frame
	float get_phi(float cx_, float cy_);//in lab frame
	int get_sector(float phi_);//Phi must be from lab frame
	float phi_center(float phi_);//Center phi within the sector


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
	void Rotate_4Vec(int set, float theta, float phi, float phie, TLorentzVector &p1); //Rotate Four vectors along the theta and phi angles
	void Boost_4Vec(float beta, TLorentzVector &p1 );// Boost a four vector in the z direction 
	void COM_gp(int set, TLorentzVector &p0, TLorentzVector &p1, TLorentzVector &p2, TLorentzVector &p3); //Bring four vectors into the COM reference frame for excited nucleon 
	float alpha(int top, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4); //Alpha angle between scattering planes
	float epsilon(int set, float Energy, float Q_2); //Virtual photon transverse polarization 
	float MM_2(TLorentzVector p1, TLorentzVector p2);//Get the MM of a two particle system
	float gamma_nu(int set, float Ep, float Q_2, float W_);//Virtual Photon Flux
	int Qfaraday(int q_last, int q_next, int q_tot);//Faraday Cup counting 


	


}

#endif