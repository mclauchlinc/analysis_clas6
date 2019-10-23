#include "event_selection.hpp"


bool Selection::Event_Selection(TLorentzVector k_mu_, Particle elec_, Particle h1_, Particle h2_){
	bool pass = false;
	float MM = physics::MM_event(0,k_mu_,elec_.Particle::Par_4Vec(),h1_.Particle::Par_4Vec(),h2_.Particle::Par_4Vec());
	if(h1_.Particle::Is_Pro() && h2_.Particle::Is_Pip()){ //Pim missing
		pass = cuts::MM_cut(2,MM);
	}else if(h1_.Particle::Is_Pro() && h2_.Particle::Is_Pim()){//Pip Missing
		pass = cuts::MM_cut(1,MM);
	}else if(h1_.Particle::Is_Pip() && h2_.Particle::Is_Pim()){//Pro Missing
		pass = cuts::MM_cut(0,MM);
	}
	return pass; 
}

bool Selection::Event_Selection(TLorentzVector k_mu_, Particle elec_, Particle h1_, Particle h2_, Particle h3_){
	bool pass = false;
	float MM = physics::MM_event(0,k_mu_,elec_.Particle::Par_4Vec(),h1_.Particle::Par_4Vec(),h2_.Particle::Par_4Vec(),h3_.Particle::Par_4Vec());
	if(h1_.Particle::Is_Pro() && h2_.Particle::Is_Pip() && h3_.Particle::Is_Pim()){ //Pim missing
		pass = cuts::MM_cut(2,MM);
	}
	return pass; 
}

/*
bool Selection::Topology(bool top_pos[4], std::shared_ptr<Histogram> hist_, std::shared_ptr<Branches> data, TLorentzVector k0, TLorentzVector k1, TLorentzVector k2, TLorentzVector k3, int idx1, int idx2, int idx3){//For all particles detected topologies
	bool pass = false; 
	if(top_pos[3]){
		float MM_z = physics::MM_event(0,0,k0,k1,k2,k3);
		float MM_z2 = physics::MM_event(0,1,k0,k1,k2,k3);
		_hists->Histogram::MM_Fill(3,MM_z,0,0);
		_hists->Histogram::MM_Fill(3,MM_z2,0,1);
		if((MM_z2 > (MM_zero_center2-MM_zero_sigma2))&&(MM_z2 < (MM_zero_center2+MM_zero_sigma2))){//Missing Mass Cut on Proton Mass
			_hists->Histogram::MM_Fill(3,MM_z,1,0);
			_hists->Histogram::MM_Fill(3,MM_z2,1,1);
			//_hists->Histogram::MM_Fill(4,MM_z,1,0);
			//_hists->Histogram::MM_Fill(4,MM_z2,1,1);
			topo[3]=true;
			std::cout<<std::endl <<"ZERO Missing topology passed";
			_hists->Histogram::WQ2_Fill(4,10,_W,_Q2);
			_hists->Histogram::Fid_Fill(4,theta[0],phi[0],0,10,0,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(4,data->Branches::p(0),data->Branches::etot(0),10,0,_W,sector[0]);
			_hists->Histogram::CC_Fill(4,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,0);
			for(int par = 0; par<3; par++){
				_hists->Histogram::Fid_Fill(4,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,_p[par]);
				_hists->Histogram::DT_Fill(4,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));
			}
		}else{
			_hists->Histogram::MM_Fill(3,MM_z,2,0);
			_hists->Histogram::MM_Fill(3,MM_z2,2,1);
			//_hists->Histogram::MM_Fill(4,MM_z,2,0);
			//_hists->Histogram::MM_Fill(4,MM_z2,2,1);
			_hists->Histogram::WQ2_Fill(4,10,_W,_Q2);
			_hists->Histogram::Fid_Fill(4,theta[0],phi[0],0,10,1,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(4,data->Branches::p(0),data->Branches::etot(0),10,1,_W,sector[0]);
			_hists->Histogram::CC_Fill(4,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,1);
			for(int par = 0; par<3; par++){
				_hists->Histogram::Fid_Fill(4,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,_p[par]);
				_hists->Histogram::DT_Fill(4,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));
			}
		}

	}
}

bool Selection::Topology(bool top_pos[4], std::shared_ptr<Histogram> hist_, std::shared_ptr<Branches> data, TLorentzVector k0, TLorentzVector k1, TLorentzVector k2, int idx1, int idx2){//For The particle Missing Topologies
	for(int i= 0; i<3; i++){
		if(top_pos[i]){

		}
	}
}*/