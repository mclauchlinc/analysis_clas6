#include "physics.hpp"


TLorentzVector physics::Make_4Vector(float p, float cx, float cy, float cz, float m){
	TVector3 k_mu_3(p*cx, p*cy, p*cz);
	TLorentzVector k_mu;
	k_mu.SetVectM(k_mu_3,m);
	return k_mu;
}

TLorentzVector physics::Set_k_mu(int set){
	TLorentzVector k_mu;
	switch(set){
		case 0:
		k_mu = k_mu_e16; //Constants.hpp
		break;
		case 1:
		k_mu = k_mu_e1f; //Constants.hpp
		break;
	}
	return k_mu; 
}

int physics::event_helicity(std::shared_ptr<Branches> data, int plate_stat){
	int eh = 0; 
	if(data->evntclas2() >= 1000) eh = 1; 
	if(data->evntclas2() <= -1000) eh = -1; 
	if(data->evntclas2() < 1000 && data->evntclas2() > -1000) eh = 0; 
	if(plate_stat == 0 ) eh = 1; 
	return plate_stat*eh; 
}

float physics::Qsquared(int set, std::shared_ptr<Branches> data){
	TVector3 k_mu_3(data->p(0)*data->cx(0),data->p(0)*data->cy(0),data->p(0)*data->cz(0));
	TLorentzVector k_mu = physics::Set_k_mu(set);
	TLorentzVector k_mu_prime; 
	k_mu_prime.SetVectM(k_mu_3,me); //Constants.hpp
	return -(k_mu - k_mu_prime).Mag2();
}

float physics::WP(int set, std::shared_ptr<Branches> data){
	TLorentzVector k_mu = physics::Set_k_mu(set);
	TVector3 k_mu_3(data->p(0)*data->cx(0),data->p(0)*data->cy(0),data->p(0)*data->cz(0));
	TLorentzVector k_mu_prime.SetVectM(k_mu_3,me);
	TLorentzVector q_mu; 
	q_mu = k_mu - k_mu_prime;
	return (p_mu + q_mu).Mag();
}

float physics::beta_calc(float m, std::shared_ptr<Branches> data, int i){
	return data->p(i)/TMath::Sqrt(m*m+data->p(i)*data->p(i));
}

float physics::MM_event(int set, TLorentzVector k1_mu, TLorentzVector k2_mu, TLorentzVector k3_mu, TLorentzVector k4_mu = {0.0,0.0,0.0,0.0}, int squared){
	TLorentzVector k_mu = physics::Set_k_mu(set);
	float MM = -99; 
	if(squared == 0 ){
		MM = (k_mu + p_mu - k1_mu - k2_mu - k3_mu - k4_mu).Mag();
	} else{
		MM = (k_mu + p_mu - k1_mu - k2_mu - k3_mu - k4_mu).Mag2();
	}
	return  MM;
}


//Math
TVector3 physics::V4_to_V3(TLorentzVector p1){
	TVector3 v1;
	for(int i = 0; i < 3; i++){
		v1[i]=p1[i];
	}
	return v1;
}

float physics::Vec3_Mag(TVector3 v1){
	return sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
}

float physics::Vec3_Mag(TLorentzVector p1){
	return sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]);
}

TVector3 physics::Cross_Product(TLorentzVector p1, TLorentzVector p2){
	TVector3 product(p1[1]*p2[2]-p1[2]*p2[1],p1[2]*p2[0]-p1[0]*p2[2],p1[0]*p2[1]-p1[1]*p2[0]); 
	return product;
}//Get the cross product between two vectors

TVector3 physics::Cross_Product(TVector3 v1, TVector3 v2){
	TVector3 product(v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0]); 
	return product;
} 
float physics::Dot_Product(TVector3 v1, TVector3 v2){
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}//Get the dot product between two vectors

float physics::Dot_Product(TLorentzVector p1, TLorentzVector p2){
	return p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2];
}

float physics::Cos_Vecs(TVector3 v1, TVector3 v2){
	double other1mag;
	double other2mag;
	double dotpro;
	dotpro = v1 * v2;
	other1mag = Vec3_Mag(v1);
	other2mag = Vec3_Mag(v2);
	return dotpro/(other2mag*other1mag);
} //Get the Cosine between two vectors 

float physics::Cos_Vecs(TLorentzVector p1, TLorentzVector p2){
	float other1mag;
	float other2mag;
	float dotpro;
	dotpro = physics::Dot_Product(p1,p2);
	other1mag = Vec3_Mag(p1);
	other2mag = Vec3_Mag(p2);
	return dotpro/(other2mag*other1mag);
}

float physics::Sin_Vecs(TVector3 v1, TVector3 v2){
	float other1mag;
	float other2mag;
	float Cross_Mag;
	TVector3 product = Cross_Product(v1,v2);
	other1mag = Vec3_Mag(v1);
	other2mag = Vec3_Mag(v2);
	Cross_Mag = Vec3_Mag(product);
	return Cross_Mag/(other1mag*other2mag);
}

float physics::Sin_Vecs(TLorentzVector p1, TLorentzVector p2){
	float other1mag;
	float other2mag;
	float Cross_Mag;
	TVector3 product = Cross_Product(p1,p2);
	other1mag = Vec3_Mag(p1);
	other2mag = Vec3_Mag(p2);
	Cross_Mag = Vec3_Mag(product);
	return Cross_Mag/(other1mag*other2mag);
} //Get the Sin between two vectors 

float physics::Get_phie(int set, TLorentzVector p0){
	TLorentzVector k_mu = physics::Get_k_mu(set);
	TVector3 nE = (1/(Vec3_Mag(k_mu)*Vec3_Mag(p0)))*Cross_Product(k_mu,p0);
	float phie = TMath::ATan2(nE[0],nE[1]);
	return phie; 
}

void physics::Rotate_4Vec(int set, float theta, float phi, float phie, TLorentzVector &p1){
	TLorentzVector k_mu = physics::Set_k_mu(set);
	p1.RotateZ(-phi);
	p1.RotateY(-theta);
	p1.RotateZ(-phie);
} //Rotate Four vectors along the theta and phi angles

void physics::Boost_4Vec(float beta, TLorentzVector &p1 ){
	p1.Boost(0.0,0.0,beta);
}// Boost a four vector in the z direction 

void physics::COM_gp(int set, TLorentzVector &p0, TLorentzVector &p1, TLorentzVector &p2, TLorentzVector &p3){
	TLorentzVector k_mu = physics::Set_k_mu(set);//Establish set
	TLorentzVector q_mu = k_mu - p0;//Four vector for virtual particle
	TLorentzVector nstar_mu = p_mu + q_mu; //Combined photon-target system
	float phigp = TMath::ATan2(nstar_mu[1],nstar_mu[0]);//Phi angle out of the x-plane
	nstar_mu.RotateZ(-phigp);//Get all horizontal momentum on x axis by roating around z axis
	float thgp - TMath::ATan2(nstar_mu[0],nstar_mu[2]); //Theta angle away from z-axis
	nstar_mu.RotateY(-thgp);//Rotate towards z-axis so all momentum is in the z direction
	float b = nstar_mu.Beta();//Get the beta to boost to the rest frame for the center of mass
	nstar_mu.Boost(0.0,0.0,-b);
	float phie = physis::Get_phie(set,p0);//Get the angle for the scattering plane of the electrons just to have a consistent definition of phi 
	physics::Rotate_4Vec(set, thgp,phigp,phie,p0);
	physics::Rotate_4Vec(set, thgp,phigp,phie,p1);
	physics::Rotate_4Vec(set, thgp,phigp,phie,p2);
	physics::Rotate_4Vec(set, thgp,phigp,phie,p3);
	physics::Boost_4Vec(-b,p0);
	physics::Boost_4Vec(-b,p1);
	physics::Boost_4Vec(-b,p2);
	physics::Boost_4Vec(-b,p3);
} //Bring four vectors into the COM reference frame for excited nucleon 

float physics::alpha(int top, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4){
	float dotpro;
	float alph; 
	float sin, cos; 
	float theta_b, phi_b, phi_c; 
	TVector3 norm1, norm2, delta_v, gamma_v, beta_v, v1, v2, v3, v4, norm3;
	switch(top){
		//v1 = Target particle (the particle whose theta and phi angles are being measured) {pi-, p, pi+}
		//v2 = Paired particle for  scattering plane {p, pp, p}
		// v3= First particle in other scattering plane {pp, pi+, pp}
		// v4 = second particle in other scattering plane{pi+, pi-, pi-}
		case 0://{pi-,p},{pp,pi+}
		v1 = physics::V4_to_V3(p3);
		v2 = physics::V4_to_V3(p4);
		v3 = physics::V4_to_V3(p1);
		v4 = physics::V4_to_V3(p2);
		break;
		case 1://{p,pp},{pi+,pi-}
		v1 = physics::V4_to_V3(p4);
		v2 = physics::V4_to_V3(p1);
		v3 = physics::V4_to_V3(p2);
		v4 = physics::V4_to_V3(p3);
		break;
		case 2://{pi+,p},{pp,pi-}
		v1 = physics::V4_to_V3(p2);
		v2 = physics::V4_to_V3(p1);
		v3 = physics::V4_to_V3(p4);
		v4 = physics::V4_to_V3(p3);
		break;
	}
	delta_v = -(1/physics::Vec3_Mag(v1))*v1;
	beta_v = (1/(physics::Vec3_Mag(physics::Cross_Product(physics::Cross_Product(v2,v4),delta_v))))*physics::Cross_Product(physics::Cross_Product(v2,v4),delta_v);
	gamma_v = -(1/(physics::Vec3_Mag(physics::Cross_Product(physics::Cross_Product(v3,v1),delta_v))))*physics::Cross_Product(physics::Cross_Product(v3,v1),delta_v);
	
	phi_b = TMath::ATan2(delta_v[1],delta_v[0]);//Angle from x axis to y-axis
	
	delta_v.RotateZ(-phi_b);
	beta_v.RotateZ(-phi_b);
	gamma_v.RotateZ(-phi_b);
	
	theta_b = TMath::ATan2(delta_v[0],delta_v[2]);
	delta_v.RotateY(-theta_b+TMath::Pi());
	beta_v.RotateY(-theta_b+TMath::Pi());
	gamma_v.RotateY(-theta_b+TMath::Pi());
	
	phi_c = TMath::ATan2(beta_v[1],beta_v[0]);
	delta_v.RotateZ(-phi_c);
	beta_v.RotateZ(-phi_c);
	gamma_v.RotateZ(-phi_c);
	alph = TMath::ATan2(gamma_v[1],gamma_v[0])*180/TMath::Pi();
	if(alph < 0){
		alph = 360 + alph; 
	} 
	
	//std::cout<<"alpha post = " <<alph <<std::endl;
	//std::cout<<std::endl <<"aCos(-.5)" <<TMath::ACos(-0.5) <<std::endl;

	
	return alph; 
} //Alpha angle between scattering planes

float physics::epsilon(int set, float Energy, float Q_2){
	float event_energy, omega; 
	switch(set){
		case 0:
		event_energy = energy_e16; //Constants.hpp
		break;
		case 1:
		event_energy = energy_e1f; //Constants.hpp
		break;
	}
	omega = event_energy - Energy;
	return 1.0/(1+2*(Q_2+omega)/(4*event_energy*Energy-Q_2));
} //Virtual photon transverse polarization 

float physics::MM_2(TLorentzVector p1, TLorentzVector p2){
	return (p1+p2).Mag();
}
//Get the MM of a two particle system
float physics::gamma_nu(int set, float Ep, float Q_2, float W_){
	float event_energy;
	switch(set){
		case 0:
		event_energy = energy_e16; //Constants.hpp
		break;
		case 1:
		event_energy = energy_e1f; //Constants.hpp
		break;
	}
	return fine_structure*W_*(W_*W_-mp*mp)/(4*TMath::Pi()*event_energy*event_energy*mp*mp*(1-epsilon(set,Ep,Q_2))*Q_2);
}
//Virtual Photon Flux
int physics::Qfaraday(int q_last, int q_next, int q_tot){
	if(q_next < q_last){
		q_tot = q_tot + q_next; 
	}
	else{
		q_tot = q_tot + (q_next - q_last);
	}
	return q_tot;
}
//Faraday Cup counting 

