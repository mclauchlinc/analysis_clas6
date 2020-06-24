#include "physics.hpp"


void physics::Print_4Vec(TLorentzVector k1){
	std::cout<<std::endl <<"Four Vector= px: " <<k1[0] <<" py: " <<k1[1] <<" pz: " <<k1[2]; 
}
void physics::Print_3Vec(TVector3 p1){
	std::cout<<std::endl <<"Three Vector= px: " <<p1[0] <<" py: " <<p1[1] <<" pz: " <<p1[2] ;
}

bool physics::Check_4Vec(TLorentzVector k1){
	bool pass = false; 
	if(k1[3]>0){//A non zero energy which all particles should have
		pass = true; 
	}
	return pass; 
}

TLorentzVector physics::Make_4Vector(float p, float cx, float cy, float cz, float m){
	TVector3 k_mu_3(p*cx, p*cy, p*cz);
	TLorentzVector k_mu;
	k_mu.SetVectM(k_mu_3,m);
	return k_mu;
}

TLorentzVector physics::Make_4Vector(float px, float py, float pz, float m){
	TVector3 k_mu_3(px, py, pz);
	TLorentzVector k_mu;
	k_mu.SetVectM(k_mu_3,m);
	return k_mu;
}

TLorentzVector physics::Make_4Vector(bool here, float p, float theta, float phi, float m){
	TLorentzVector k_mu;
	if(here){
		TVector3 k_mu_3(p*TMath::Sin(theta*TMath::Pi()/180.0)*TMath::Cos(phi*TMath::Pi()/180.0), p*TMath::Sin(theta*TMath::Pi()/180.0)*TMath::Sin(phi*TMath::Pi()/180.0), p*TMath::Cos(theta*TMath::Pi()/180.0));
		k_mu.SetVectM(k_mu_3,m);
	}
	return k_mu;
}

TLorentzVector physics::Set_k_mu(int set){
	TLorentzVector k_mu;
	switch(set%2){
		case 1:
		k_mu = k_mu_e16; //Constants.hpp
		break;
		case 0:
		k_mu = k_mu_e1f; //Constants.hpp
		break;
		default:
		k_mu = k_mu_e16;
		break;
	}
	return k_mu; 
}

int physics::event_helicity(std::shared_ptr<Branches> data, int plate_stat){
	int eh = 0; 
	if(data->evntclas2() >= 1000) eh = 1; 
	if(data->evntclas2() <= -1000) eh = -1; 
	if(data->evntclas2() < 1000 && data->evntclas2() > -1000) eh = 0; 
	//if(plate_stat == 0 ) eh = 1; 
	return plate_stat*eh; 
}

float physics::Qsquared(int set, std::shared_ptr<Branches> data, int thr){
	TLorentzVector k_mu_prime; 
	//TLorentzVector k_mu_prime; 
	if(thr == 1){
		k_mu_prime = Make_4Vector(true,data->mcp(0),data->mctheta(0),data->mcphi(0),me);
	}else{
		k_mu_prime = Make_4Vector(data->p(0)*data->cx(0),data->p(0)*data->cy(0),data->p(0)*data->cz(0),me);
	}
	TLorentzVector k_mu = physics::Set_k_mu(set);
	return -(k_mu - k_mu_prime).Mag2();
}

float physics::Qsquared(int set, std::shared_ptr<Branches> data, bool thrown_){
	TLorentzVector k_mu_prime; 
	if(thrown_){
		k_mu_prime = Make_4Vector(true,data->mcp(0),data->mctheta(0),data->mcphi(0),me);
	}else{
		k_mu_prime = Make_4Vector(data->p(0)*data->cx(0),data->p(0)*data->cy(0),data->p(0)*data->cz(0),me);
	}
	TLorentzVector k_mu = physics::Set_k_mu(set);
	return -(k_mu - k_mu_prime).Mag2();
}


float physics::WP(int set, std::shared_ptr<Branches> data, int thr){
	TLorentzVector k_mu = physics::Set_k_mu(set);
	TLorentzVector k_mu_prime; 
	if(thr==1){
		k_mu_prime = Make_4Vector(true,data->mcp(0),data->mctheta(0),data->mcphi(0),me);
			//Print_4Vec(k_mu_prime);
			//std::cout<<std::endl <<"again particle ID: " <<data->pidpart(0) <<" px: " <<data->pxpart(0) <<" py: " <<data->pypart(0) <<" pz: " <<data->pzpart(0);  
	}else{
			k_mu_prime = Make_4Vector(data->p(0)*data->cx(0),data->p(0)*data->cy(0),data->p(0)*data->cz(0),me);
	}
	TLorentzVector q_mu; 
	q_mu = k_mu - k_mu_prime;
	//std::cout<<" W: " <<(p_mu + q_mu).Mag();
	return (p_mu + q_mu).Mag();
}


float physics::beta_calc(float m, std::shared_ptr<Branches> data, int i){
	return data->p(i)/TMath::Sqrt(m*m+data->p(i)*data->p(i));
}

float physics::MM_event(int set, int squared, TLorentzVector k1_mu, TLorentzVector k2_mu, TLorentzVector k3_mu, TLorentzVector k4_mu){
	TLorentzVector k_mu = physics::Set_k_mu(set);
	float MM = -99; 
	if(squared == 0 ){
		MM = (k_mu + p_mu - k1_mu - k2_mu - k3_mu - k4_mu).Mag();
	} else{
		MM = (k_mu + p_mu - k1_mu - k2_mu - k3_mu - k4_mu).Mag2();
	}
	return  MM;
}

float physics::MM_event(int squared, TLorentzVector k0_mu, TLorentzVector k1_mu, TLorentzVector k2_mu, TLorentzVector k3_mu, TLorentzVector k4_mu){
	float MM = NAN; 
	if(squared == 0 ){
		MM = (k0_mu + p_mu - k1_mu - k2_mu - k3_mu - k4_mu).Mag();
	} else{
		MM = (k0_mu + p_mu - k1_mu - k2_mu - k3_mu - k4_mu).Mag2();
	}
	return  MM;
}

float physics::MM_event(int squared, TLorentzVector k0_mu, TLorentzVector k1_mu, TLorentzVector k2_mu, TLorentzVector k3_mu){
	float MM = NAN; 
	if(squared == 0 ){
		MM = (k0_mu + p_mu - k1_mu - k2_mu - k3_mu).Mag();
	} else{
		MM = (k0_mu + p_mu - k1_mu - k2_mu - k3_mu).Mag2();
	}
	return  MM;
}

float physics::get_theta(float cz_){
	float degree = 180.0/TMath::Pi();
	return TMath::ACos(cz_)*degree;
}

float physics::get_theta(int part, std::shared_ptr<Branches> data){
	return physics::get_theta(data->Branches::cz(part));
}

float physics::get_phi(float cx_, float cy_){
	float degree = 180.0/TMath::Pi();
	return TMath::ATan2(cy_,cx_)*degree;
}

float physics::get_phi(int part, std::shared_ptr<Branches> data){
	return physics::get_phi(data->Branches::cx(part), data->Branches::cy(part));
}

float physics::get_phi_pos(float cx_, float cy_){
	float phi = physics::get_phi(cx_,cy_);
	if(phi<0){
		phi = 180.0 - phi; 
	}
	return phi;
}

int physics::get_sector(float phi_){
	int sector;
	if(phi_>=-30 && phi_ <=30)
	{
		sector = 1;
	}else
	if(phi_>=30 && phi_<=90)
	{
		sector = 2;
	}else
	if(phi_>=90 && phi_ <=150)
	{
		sector = 3;
	}else
	if(phi_>=150 || phi_<=-150)
	{
		sector = 4;
	}else
	if(phi_>=-150 && phi_<=-90)
	{
		sector = 5;
	}else
	if(phi_>=-90 && phi_<=-30)
	{
		sector = 6;
	}//got rid of pointless "else" statement 3/10/2017
	return sector;
}

float physics::phi_center( float phi_)
{
	double phi_corr;
	int sector = get_sector(phi_);

	if(sector ==1)
	{
		phi_corr = phi_;
	}else
	if(sector==2)
	{
		phi_corr = phi_-60;
	}else
	if(sector==3)
	{
		phi_corr = phi_-120;
	}else
	if(sector == 4)
	{
		if(phi_<=-150)
		{
			phi_corr = phi_+180;
		}
		if(phi_>=150)
		{
			phi_corr = phi_-180;
		}
	}else
	if(sector==5)
	{
		phi_corr = phi_+120;
	}else
	if(sector==6)
	{
		phi_corr = phi_+60;
	}

	//Not working, but elegant. Adjust for elegance later
	/*phi0 = phi0 +180;
	int mod6 = ((int)phi0+210)/60;
	phi_corr = phi0 - (double)mod6*60.0;
	*/
	return phi_corr;
}

//Delta t
float physics::vert_e(float d, float t){
	return t-(d/c_special);
}

float physics::vert_h(float p, float d, float t, float m){
	return t-((d/c_special)*sqrt(1.0 + m*m/(p*p)));
}

float physics::delta_t(int part, float p, float d, float t, float d0, float t0){
	float mass = -99; 
	switch(part){
		case 0:
			mass = mp; 
		break;
		case 1:
			mass = mpi;
		break;
		case 2:
			mass = mpi;
		break;
	}
	float vertex_e = physics::vert_e(d0,t0);
	float vertex_h = physics::vert_h(p,d,t,mass);
	return vertex_e - vertex_h;
}

float physics::delta_t(int part, std::shared_ptr<Branches> data, int idx){
	float mass = -99;
	switch(part){
		case 0:
			mass = me; 
		break;
		case 1:
			mass = mp; 
		break;
		case 2:
			mass = mpi;
		break;
		case 3:
			mass = mpi;
		break;
	}
	float vertex_e = physics::vert_e(data->Branches::sc_r(0),data->Branches::sc_t(0));
	float vertex_h = physics::vert_h(data->Branches::p(idx),data->Branches::sc_r(idx),data->Branches::sc_t(idx),mass);
	return vertex_e - vertex_h;
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
	TLorentzVector k_mu = physics::Set_k_mu(set);
	TVector3 nE = (1/(Vec3_Mag(k_mu)*Vec3_Mag(p0)))*Cross_Product(k_mu,p0);
	float phie = TMath::ATan2(nE[0],nE[1]);
	return phie; 
}

float physics::Get_phie(TLorentzVector k0, TLorentzVector p0){
	TVector3 nE = (1/(Vec3_Mag(k0)*Vec3_Mag(p0)))*Cross_Product(k0,p0);
	float phie = TMath::ATan2(nE[0],nE[1]);
	return phie; 
}

void physics::Rotate_4Vec(int set, float theta, float phi, float phie, TLorentzVector &p1){
	TLorentzVector k_mu = physics::Set_k_mu(set);
	p1.RotateZ(-phi);
	p1.RotateY(-theta);
	p1.RotateZ(-phie);
} //Rotate Four vectors along the theta and phi angles

void physics::Rotate_4Vec(float theta, float phi, float phie, TLorentzVector &p1){
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
	float thgp = TMath::ATan2(nstar_mu[0],nstar_mu[2]); //Theta angle away from z-axis
	nstar_mu.RotateY(-thgp);//Rotate towards z-axis so all momentum is in the z direction
	float b = nstar_mu.Beta();//Get the beta to boost to the rest frame for the center of mass
	nstar_mu.Boost(0.0,0.0,-b);
	float phie = physics::Get_phie(set,p0);//Get the angle for the scattering plane of the electrons just to have a consistent definition of phi 
	physics::Rotate_4Vec(set, thgp,phigp,phie,p0);
	physics::Rotate_4Vec(set, thgp,phigp,phie,p1);
	physics::Rotate_4Vec(set, thgp,phigp,phie,p2);
	physics::Rotate_4Vec(set, thgp,phigp,phie,p3);
	physics::Boost_4Vec(-b,p0);
	physics::Boost_4Vec(-b,p1);
	physics::Boost_4Vec(-b,p2);
	physics::Boost_4Vec(-b,p3);
} //Bring four vectors into the COM reference frame for excited nucleon 

void COM_gp(TLorentzVector &k0, TLorentzVector &p0, TLorentzVector &p1, TLorentzVector &p2, TLorentzVector &p3){
	TLorentzVector q_mu = k0 - p0;//Four vector for virtual particle
	TLorentzVector nstar_mu = p_mu + q_mu; //Combined photon-target system
	float phigp = TMath::ATan2(nstar_mu[1],nstar_mu[0]);//Phi angle out of the x-plane
	nstar_mu.RotateZ(-phigp);//Get all horizontal momentum on x axis by roating around z axis
	float thgp = TMath::ATan2(nstar_mu[0],nstar_mu[2]); //Theta angle away from z-axis
	nstar_mu.RotateY(-thgp);//Rotate towards z-axis so all momentum is in the z direction
	float b = nstar_mu.Beta();//Get the beta to boost to the rest frame for the center of mass
	nstar_mu.Boost(0.0,0.0,-b);
	float phie = physics::Get_phie(k0,p0);//Get the angle for the scattering plane of the electrons just to have a consistent definition of phi 
	physics::Rotate_4Vec(thgp,phigp,phie,p0);
	physics::Rotate_4Vec(thgp,phigp,phie,p1);
	physics::Rotate_4Vec(thgp,phigp,phie,p2);
	physics::Rotate_4Vec(thgp,phigp,phie,p3);
	physics::Boost_4Vec(-b,p0);
	physics::Boost_4Vec(-b,p1);
	physics::Boost_4Vec(-b,p2);
	physics::Boost_4Vec(-b,p3);
}

TLorentzVector physics::COM_gp(int par, TLorentzVector k0, TLorentzVector p0, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3){
	TLorentzVector output;
	TLorentzVector q_mu = k0 - p0;//Four vector for virtual particle
	TLorentzVector nstar_mu = p_mu + q_mu; //Combined photon-target system
	switch(par){
		case 0: output = p0; break;
		case 1: output = p1; break;
		case 2: output = p2; break;
		case 3: output = p3; break;
		case 4: output = k0; break;
		case 5: output = p_mu; break;
	}
	float phigp = TMath::ATan2(nstar_mu[1],nstar_mu[0]);//Phi angle out of the x-plane
	nstar_mu.RotateZ(-phigp);//Get all horizontal momentum on x axis by roating around z axis
	float thgp = TMath::ATan2(nstar_mu[0],nstar_mu[2]); //Theta angle away from z-axis
	nstar_mu.RotateY(-thgp);//Rotate towards z-axis so all momentum is in the z direction
	float b = nstar_mu.Beta();//Get the beta to boost to the rest frame for the center of mass
	nstar_mu.Boost(0.0,0.0,-b);
	float phie = physics::Get_phie(k0,output);//Get the angle for the scattering plane of the electrons just to have a consistent definition of phi 
	physics::Rotate_4Vec(thgp,phigp,phie,output);
	physics::Boost_4Vec(-b,output);
	return output;
}

float physics::alpha(int top, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, int set){
	float dotpro;
	float alph; 
	float sin, cos; 
	float theta_b, phi_b, phi_c; 
	TVector3 norm1, norm2, delta_v, gamma_v, beta_v, v1, v2, v3, v4, norm3;
	physics::COM_gp(set,p1,p2,p3,p4);
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
	
	//Make them all unit vectors. 
	v1 = (1.0/physics::Vec3_Mag(v1))*v1;
	v2 = (1.0/physics::Vec3_Mag(v2))*v2;
	v3 = (1.0/physics::Vec3_Mag(v3))*v3;
	v4 = (1.0/physics::Vec3_Mag(v4))*v4;
	/*std::cout<<std::endl <<"v1: "; 
	physics::Print_3Vec(v1);
	std::cout<<std::endl <<"v2: "; 
	physics::Print_3Vec(v2);
	std::cout<<std::endl <<"v3: "; 
	physics::Print_3Vec(v3);
	std::cout<<std::endl <<"v4: "; 
	physics::Print_3Vec(v4);*/
	delta_v = -v1;
	//std::cout<<std::endl <<"delta_v: "; 
	//physics::Print_3Vec(delta_v);
	beta_v = physics::Cross_Product(physics::Cross_Product(v2,v4),delta_v);
	//std::cout<<std::endl <<"Beta: "; 
	//physics::Print_3Vec(beta_v);
	gamma_v = physics::Cross_Product(delta_v,physics::Cross_Product(v3,v1));
	//std::cout<<std::endl <<"gamma: "; 
	//physics::Print_3Vec(gamma_v);
	
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
	if(alph < 0.0){
		alph = 360.0 + alph; 
	} 
	//std::cout<<std::endl <<"alpha: " <<alph;
	
	//std::cout<<"alpha post = " <<alph <<std::endl;
	//std::cout<<std::endl <<"aCos(-.5)" <<TMath::ACos(-0.5) <<std::endl;

	
	return alph; 
} //Alpha angle between scattering planes

float physics::alpha(int top, TLorentzVector k0, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, bool COM){
	float dotpro;
	float alph; 
	float sin, cos; 
	float theta_b, phi_b, phi_c; 
	TVector3 norm1, norm2, delta_v, gamma_v, beta_v, v1, v2, v3, v4, norm3;
	//if(!COM){
	//	physics::COM_gp(k0,p1,p2,p3,p4);
	//}
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
	
	//Make them all unit vectors. 
	v1 = (1.0/physics::Vec3_Mag(v1))*v1;
	v2 = (1.0/physics::Vec3_Mag(v2))*v2;
	v3 = (1.0/physics::Vec3_Mag(v3))*v3;
	v4 = (1.0/physics::Vec3_Mag(v4))*v4;
	/*std::cout<<std::endl <<"v1: "; 
	physics::Print_3Vec(v1);
	std::cout<<std::endl <<"v2: "; 
	physics::Print_3Vec(v2);
	std::cout<<std::endl <<"v3: "; 
	physics::Print_3Vec(v3);
	std::cout<<std::endl <<"v4: "; 
	physics::Print_3Vec(v4);*/
	delta_v = -v1;
	//std::cout<<std::endl <<"delta_v: "; 
	//physics::Print_3Vec(delta_v);
	beta_v = physics::Cross_Product(physics::Cross_Product(v2,v4),delta_v);
	//std::cout<<std::endl <<"Beta: "; 
	//physics::Print_3Vec(beta_v);
	gamma_v = physics::Cross_Product(delta_v,physics::Cross_Product(v3,v1));
	//std::cout<<std::endl <<"gamma: "; 
	//physics::Print_3Vec(gamma_v);
	
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
	if(alph < 0.0){
		alph = 360.0 + alph; 
	} 
	//std::cout<<std::endl <<"alpha: " <<alph;
	
	//std::cout<<"alpha post = " <<alph <<std::endl;
	//std::cout<<std::endl <<"aCos(-.5)" <<TMath::ACos(-0.5) <<std::endl;

	
	return alph; 
}

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

float physics::epsilon(TLorentzVector k0, float Energy, float Q_2){
	float event_energy, omega; 
	event_energy = k0[3];
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
float physics::Ev_Theta(int top, TLorentzVector k0, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, bool COM){
	float event_theta = NAN;
	//if(!COM){
	//	physics::COM_gp(k0,p1,p2,p3,p4);
	//}
	switch(top){
		case 0: event_theta = physics::get_theta(p4[2]/physics::Vec3_Mag(p4)); break;//looking for pim
		case 1: event_theta = physics::get_theta(p2[2]/physics::Vec3_Mag(p2)); break;//looking for pro
		case 2: event_theta = physics::get_theta(p3[2]/physics::Vec3_Mag(p3)); break;//looking for pip
	}
	return event_theta;
}

float physics::Ev_MM(int top, TLorentzVector k0, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, bool COM){
	float event_MM = NAN;
	//if(!COM){
	//	physics::COM_gp(k0,p1,p2,p3,p4);
	//}
	switch(top){
		case 0: event_MM = physics::MM_2(p2,p3); break;//looking for pim -> pro/pip
		case 1: event_MM = physics::MM_2(p3,p4); break;//looking for pro -> pip/pim
		case 2: event_MM = physics::MM_2(p2,p4); break;//looking for pip -> pro/pim
	}
	return event_MM;
}



