#include "event.hpp"

Event::Event(){

}
void Event::COM_4Vec(){
			  //physics::COM_gp(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
	_vec[0] = physics::COM_gp(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
	_vec[1] = physics::COM_gp(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
	_vec[2] = physics::COM_gp(2,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
	_vec[3] = physics::COM_gp(3,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
}

void Event::Vars(){
	for(int i = 0; i<3; i++){
		_alphab[i] = physics::alpha(i, _k1, _vec[0], _vec[1], _vec[2], _vec[3], true);
		_thetab[i] = physics::Ev_Theta(i, _k1, _vec[0], _vec[1], _vec[2], _vec[3], true);
		_MMb[i] = physics::Ev_MM(i, _k1, _vec[0], _vec[1], _vec[2], _vec[3], true);
	}
}

void Event::Fill_Event(std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, int top_, float W_, float Q2_, Particle &p1, Particle &p2, Particle &p3, int hel_){
	if(top_ == 3){
		std::cout<<"Incorrectly Filled Event: Not enough particles for topology" <<std::endl;
	}else{ 
		_W = W_;
		_Q2 = Q2_;
		_hel = hel_;
		std::cout<<"Fill event p1" <<std::endl;
		std::cout<<"		Inside Event pt1 Fill for event Combo: " <<std::endl;
		std::cout<<"			Elec " <<p1.Particle::Get_idx() <<" sim: " <<p1.Particle::Is_Elec() <<"| Pro " <<p2.Particle::Get_idx() <<" sim: " <<p2.Particle::Is_Pro() <<"| Pip " <<p3.Particle::Get_idx() <<" sim: " <<p3.Particle::Is_Pip() <<std::endl;
		//Check to make sure all particles were either simulated or not
		if(p1.Particle::Is_Sim() == p2.Particle::Is_Sim() && p2.Particle::Is_Sim() == p3.Particle::Is_Sim()){
			_sim = p1.Particle::Is_Sim();
			std::cout<<"Fill event p2" <<std::endl;
			std::cout<<"		Inside Event pt2 Fill for event Combo: " <<std::endl;
			std::cout<<"			Elec " <<p1.Particle::Get_idx() <<" sim: " <<p1.Particle::Is_Elec() <<"| Pro " <<p2.Particle::Get_idx() <<" sim: " <<p2.Particle::Is_Pro() <<"| Pip " <<p3.Particle::Get_idx() <<" sim: " <<p3.Particle::Is_Pip() <<std::endl;
			if(p1.Particle::Is_Thrown() == p2.Particle::Is_Thrown() && p2.Particle::Is_Thrown() == p3.Particle::Is_Thrown()){
				std::cout<<"Fill event p3" <<std::endl;
				std::cout<<"		Inside Event pt3a Fill for event Combo: " <<std::endl;
				std::cout<<"			Elec " <<p1.Particle::Get_idx() <<" sim: " <<p1.Particle::Is_Elec() <<"| Pro " <<p2.Particle::Get_idx() <<" sim: " <<p2.Particle::Is_Pro() <<"| Pip " <<p3.Particle::Get_idx() <<" sim: " <<p3.Particle::Is_Pip() <<std::endl;
				_thrown = p1.Particle::Is_Thrown();
				std::cout<<"		Inside Event pt3b Fill for event Combo: " <<std::endl;
				std::cout<<"			Elec " <<p1.Particle::Get_idx() <<" sim: " <<p1.Particle::Is_Elec() <<"| Pro " <<p2.Particle::Get_idx() <<" sim: " <<p2.Particle::Is_Pro() <<"| Pip " <<p3.Particle::Get_idx() <<" sim: " <<p3.Particle::Is_Pip() <<std::endl;
				_k1 = physics::Set_k_mu(p1.Particle::Get_set());
				std::cout<<"		Inside Event pt3c Fill for event Combo: " <<std::endl;
				std::cout<<"			Elec " <<p1.Particle::Get_idx() <<" sim: " <<p1.Particle::Is_Elec() <<"| Pro " <<p2.Particle::Get_idx() <<" sim: " <<p2.Particle::Is_Pro() <<"| Pip " <<p3.Particle::Get_idx() <<" sim: " <<p3.Particle::Is_Pip() <<std::endl;
				switch(top_){
					case 0: //Electron, Pip, Pim => Pro Missing
						if(p1.Particle::Is_Elec() && p2.Particle::Is_Pip() && p3.Particle::Is_Pim()){
							_filled_correctly = true; 
							_p_lab[0] = p1.Particle::Get_p();
							_p_lab[2] = p2.Particle::Get_p();
							_p_lab[3] = p3.Particle::Get_p();
							_theta_lab[0] = p1.Particle::Get_theta();
							_theta_lab[2] = p2.Particle::Get_theta();
							_theta_lab[3] = p3.Particle::Get_theta();
							_phi_lab[0] = p1.Particle::Get_phi();
							_phi_lab[2] = p2.Particle::Get_phi();
							_phi_lab[3] = p3.Particle::Get_phi();
							_vec_lab[0] = physics::Make_4Vector(true,_p_lab[0],_theta_lab[0],_phi_lab[0],me);
							_vec_lab[2] = physics::Make_4Vector(true,_p_lab[2],_theta_lab[2],_phi_lab[2],mpi);
							_vec_lab[3] = physics::Make_4Vector(true,_p_lab[3],_theta_lab[3],_phi_lab[3],mpi);
							_top[0] = true;
							hist_->Histogram::MM_Fill(envi_,0,physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[2],_vec_lab[3]),0,0,true);
							hist_->Histogram::MM_Fill(envi_,0,physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[2],_vec_lab[3]),0,1,true);
							_pass = Selection::Event_Selection(top_,_k1, _vec_lab[0], _vec_lab[2], _vec_lab[3]);
							p1.Particle::Fill_Par_Event(envi_,hist_,_W,top_,0,_pass); 
							p2.Particle::Fill_Par_Event(envi_,hist_,_W,top_,2,_pass); 
							p3.Particle::Fill_Par_Event(envi_,hist_,_W,top_,3,_pass); 
							if(_pass){
								_vec_lab[1] = _k1 + p_mu - _vec_lab[0] - _vec_lab[2] - _vec_lab[3];
								hist_->Histogram::MM_Fill(envi_,0,physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[2],_vec_lab[3]),1,0,true);
								hist_->Histogram::MM_Fill(envi_,0,physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[2],_vec_lab[3]),1,1,true);
							}else{
								hist_->Histogram::MM_Fill(envi_,0,physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[2],_vec_lab[3]),2,0,true);
								hist_->Histogram::MM_Fill(envi_,0,physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[2],_vec_lab[3]),2,1,true);
							}
						}
					break;
					case 1://Pip missing
						if(p1.Particle::Is_Elec() && p2.Particle::Is_Pro() && p3.Particle::Is_Pim()){
							_filled_correctly = true; 
							_p_lab[0] = p1.Particle::Get_p();
							_p_lab[1] = p2.Particle::Get_p();
							_p_lab[3] = p3.Particle::Get_p();
							_theta_lab[0] = p1.Particle::Get_theta();
							_theta_lab[1] = p2.Particle::Get_theta();
							_theta_lab[3] = p3.Particle::Get_theta();
							_phi_lab[0] = p1.Particle::Get_phi();
							_phi_lab[1] = p2.Particle::Get_phi();
							_phi_lab[3] = p3.Particle::Get_phi();
							_vec_lab[0] = physics::Make_4Vector(true,_p_lab[0],_theta_lab[0],_phi_lab[0],me);
							_vec_lab[1] = physics::Make_4Vector(true,_p_lab[1],_theta_lab[1],_phi_lab[1],mp);
							_vec_lab[3] = physics::Make_4Vector(true,_p_lab[3],_theta_lab[3],_phi_lab[3],mpi);
							_top[1] = true;
							hist_->Histogram::MM_Fill(envi_,1,physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[3]),0,0,true);
							hist_->Histogram::MM_Fill(envi_,1,physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[3]),0,1,true);
							_pass = Selection::Event_Selection(top_,_k1, _vec_lab[0], _vec_lab[1], _vec_lab[3]);
							p1.Particle::Fill_Par_Event(envi_,hist_,_W,top_,0,_pass); 
							p2.Particle::Fill_Par_Event(envi_,hist_,_W,top_,1,_pass); 
							p3.Particle::Fill_Par_Event(envi_,hist_,_W,top_,3,_pass);
							if(_pass){
								_vec_lab[2] = _k1 + p_mu - _vec_lab[0] - _vec_lab[1] - _vec_lab[3]; 
								hist_->Histogram::MM_Fill(envi_,1,physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[3]),1,0,true);
								hist_->Histogram::MM_Fill(envi_,1,physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[3]),1,1,true);
							}else{
								hist_->Histogram::MM_Fill(envi_,1,physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[3]),2,0,true);
								hist_->Histogram::MM_Fill(envi_,1,physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[3]),2,1,true);
							}
						}
					break;
					case 2://PIM missing
						std::cout<<"Fill event p4" <<std::endl;
						std::cout<<"		Inside Event pt1 Fill for event Combo: " <<std::endl;
							std::cout<<"			Elec " <<p1.Particle::Get_idx() <<" sim: " <<p1.Particle::Is_Elec() <<"| Pro " <<p2.Particle::Get_idx() <<" sim: " <<p2.Particle::Is_Pro() <<"| Pip " <<p3.Particle::Get_idx() <<" sim: " <<p3.Particle::Is_Pip() <<std::endl;
						if(p1.Particle::Is_Elec() && p2.Particle::Is_Pro() && p3.Particle::Is_Pip()){
							std::cout<<"		Inside Event pt5 Fill for event Combo: " <<std::endl;
							std::cout<<"			Elec " <<p1.Particle::Get_idx() <<" sim: " <<p1.Particle::Is_Elec() <<"| Pro " <<p2.Particle::Get_idx() <<" sim: " <<p2.Particle::Is_Pro() <<"| Pip " <<p3.Particle::Get_idx() <<" sim: " <<p3.Particle::Is_Pip() <<std::endl;
							_filled_correctly = true; 
							_p_lab[0] = p1.Particle::Get_p();
							_p_lab[1] = p2.Particle::Get_p();
							_p_lab[2] = p3.Particle::Get_p();
							_theta_lab[0] = p1.Particle::Get_theta();
							_theta_lab[1] = p2.Particle::Get_theta();
							_theta_lab[2] = p3.Particle::Get_theta();
							_phi_lab[0] = p1.Particle::Get_phi();
							_phi_lab[1] = p2.Particle::Get_phi();
							_phi_lab[2] = p3.Particle::Get_phi();
							_vec_lab[0] = physics::Make_4Vector(true,_p_lab[0],_theta_lab[0],_phi_lab[0],me);
							_vec_lab[1] = physics::Make_4Vector(true,_p_lab[1],_theta_lab[1],_phi_lab[1],mp);
							_vec_lab[2] = physics::Make_4Vector(true,_p_lab[2],_theta_lab[2],_phi_lab[2],mpi);
							_top[2] = true;
							hist_->Histogram::MM_Fill(envi_,2,physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2]),0,0,true);
							hist_->Histogram::MM_Fill(envi_,2,physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2]),0,1,true);
							std::cout<<"		Inside Event pt6 Fill for event Combo: " <<std::endl;
							std::cout<<"			Elec " <<p1.Particle::Get_idx() <<" sim: " <<p1.Particle::Is_Elec() <<"| Pro " <<p2.Particle::Get_idx() <<" sim: " <<p2.Particle::Is_Pro() <<"| Pip " <<p3.Particle::Get_idx() <<" sim: " <<p3.Particle::Is_Pip() <<std::endl;
							_pass = Selection::Event_Selection(top_,_k1, _vec_lab[0], _vec_lab[1], _vec_lab[2]);
							p1.Particle::Fill_Par_Event(envi_,hist_,_W,top_,0,_pass); 
							p2.Particle::Fill_Par_Event(envi_,hist_,_W,top_,1,_pass); 
							p3.Particle::Fill_Par_Event(envi_,hist_,_W,top_,2,_pass);
							std::cout<<"		Inside Event pt7 Fill for event Combo: " <<std::endl;
							std::cout<<"			Elec " <<p1.Particle::Get_idx() <<" sim: " <<p1.Particle::Is_Elec() <<"| Pro " <<p2.Particle::Get_idx() <<" sim: " <<p2.Particle::Is_Pro() <<"| Pip " <<p3.Particle::Get_idx() <<" sim: " <<p3.Particle::Is_Pip() <<std::endl;
							if(_pass){
								_vec_lab[3] = _k1 + p_mu - _vec_lab[0] - _vec_lab[1] - _vec_lab[2]; 
								hist_->Histogram::MM_Fill(envi_,2,physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2]),1,0,true);
								hist_->Histogram::MM_Fill(envi_,2,physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2]),1,1,true);
								std::cout<<"		Inside Event pt8 Fill for event Combo: " <<std::endl;
								std::cout<<"			Elec " <<p1.Particle::Get_idx() <<" sim: " <<p1.Particle::Is_Elec() <<"| Pro " <<p2.Particle::Get_idx() <<" sim: " <<p2.Particle::Is_Pro() <<"| Pip " <<p3.Particle::Get_idx() <<" sim: " <<p3.Particle::Is_Pip() <<std::endl;
							}else{
								hist_->Histogram::MM_Fill(envi_,2,physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2]),2,0,true);
								hist_->Histogram::MM_Fill(envi_,2,physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2]),2,1,true);
								std::cout<<"		Inside Event pt9 Fill for event Combo: " <<std::endl;
								std::cout<<"			Elec " <<p1.Particle::Get_idx() <<" sim: " <<p1.Particle::Is_Sim() <<"| Pro " <<p2.Particle::Get_idx() <<" sim: " <<p2.Particle::Is_Pro() <<"| Pip " <<p3.Particle::Get_idx() <<" sim: " <<p3.Particle::Is_Pip() <<std::endl;
							}
						}
					break;
				}
			}else{
				std::cout<<"			You have mixed thrown and unthrown events" <<std::endl;
			}
		}else{
			std::cout<<"			You have mixed simulated and experimental events" <<std::endl;
			std::cout<<"			Topology: " <<top_ <<std::endl;
			std::cout<<"			p1 sim status: " <<p1.Particle::Is_Sim() <<std::endl;
			std::cout<<"			p2 sim status: " <<p2.Particle::Is_Sim() <<std::endl;
			std::cout<<"			p3 sim status: " <<p3.Particle::Is_Sim() <<std::endl;
		}
	}
}

void Event::Fill_Event(std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, int top_, float W_, float Q2_, Particle &p1, Particle &p2, Particle &p3, Particle &p4, int hel_){
if(top_ != 3){
		std::cout<<"Incorrectly Filled Event: Too many particles for topology" <<std::endl;
	}else{ 
		_W = W_;
		_Q2 = Q2_;
		_hel = hel_;
		//Check to make sure all particles were either simulated or not
		if(p1.Particle::Is_Sim() == p2.Particle::Is_Sim() && p2.Particle::Is_Sim() == p3.Particle::Is_Sim()){
			_sim = p1.Particle::Is_Sim();
			if(p1.Particle::Is_Thrown() == p2.Particle::Is_Thrown() && p2.Particle::Is_Thrown() == p3.Particle::Is_Thrown()){
				_thrown = p1.Particle::Is_Thrown();
				_k1 = physics::Set_k_mu(p1.Particle::Get_set());
				if(p1.Particle::Is_Elec() && p2.Particle::Is_Pro() && p3.Particle::Is_Pip() && p4.Particle::Is_Pim()){
					_filled_correctly = true; 
					_p_lab[0] = p1.Particle::Get_p();
					_p_lab[1] = p2.Particle::Get_p();
					_p_lab[2] = p3.Particle::Get_p();
					_p_lab[3] = p4.Particle::Get_p();
					_theta_lab[0] = p1.Particle::Get_theta();
					_theta_lab[1] = p2.Particle::Get_theta();
					_theta_lab[2] = p3.Particle::Get_theta();
					_theta_lab[3] = p4.Particle::Get_theta();
					_phi_lab[0] = p1.Particle::Get_phi();
					_phi_lab[1] = p2.Particle::Get_phi();
					_phi_lab[2] = p3.Particle::Get_phi();
					_phi_lab[3] = p4.Particle::Get_phi();
					_vec_lab[0] = physics::Make_4Vector(true,_p_lab[0],_theta_lab[0],_phi_lab[0],me);
					_vec_lab[1] = physics::Make_4Vector(true,_p_lab[1],_theta_lab[1],_phi_lab[0],mp);
					_vec_lab[2] = physics::Make_4Vector(true,_p_lab[2],_theta_lab[2],_phi_lab[2],mpi);
					_vec_lab[3] = physics::Make_4Vector(true,_p_lab[3],_theta_lab[3],_phi_lab[3],mpi);
					_top[3] = true; 
					hist_->Histogram::MM_Fill(envi_,3,physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]),0,0,true);
					hist_->Histogram::MM_Fill(envi_,3,physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[2],_vec_lab[2],_vec_lab[3]),0,1,true);
					_pass = Selection::Event_Selection(top_,_k1, _vec_lab[0], _vec_lab[1], _vec_lab[2], _vec_lab[3]);
					p1.Particle::Fill_Par_Event(envi_,hist_,_W,top_,0,_pass); 
					p2.Particle::Fill_Par_Event(envi_,hist_,_W,top_,1,_pass); 
					p3.Particle::Fill_Par_Event(envi_,hist_,_W,top_,2,_pass); 
					p4.Particle::Fill_Par_Event(envi_,hist_,_W,top_,3,_pass); 
					if(_pass){
						hist_->Histogram::MM_Fill(envi_,3,physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]),1,0,true);
						hist_->Histogram::MM_Fill(envi_,3,physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]),1,1,true);
						_vec[0] = physics::COM_gp(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
						_vec[1] = physics::COM_gp(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
						_vec[2] = physics::COM_gp(2,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
						_vec[3] = physics::COM_gp(3,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);

					}else{
						hist_->Histogram::MM_Fill(envi_,3,physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]),2,0,true);
						hist_->Histogram::MM_Fill(envi_,3,physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]),2,1,true);
					}
				}
			}else{
				std::cout<<"		You have mixed thrown and unthrown events" <<std::endl;
			}
		}else{
			std::cout<<"		You have mixed simulated and experimental events" <<std::endl;
			std::cout<<"		Topology: " <<top_ <<std::endl;
			std::cout<<"		p1 sim status: " <<p1.Particle::Is_Sim() <<std::endl;
			std::cout<<"		p2 sim status: " <<p2.Particle::Is_Sim() <<std::endl;
			std::cout<<"		p3 sim status: " <<p3.Particle::Is_Sim() <<std::endl;
			std::cout<<"		p4 sim status: " <<p4.Particle::Is_Sim() <<std::endl;
		}
	}
}

void Fill_Event_Hists(std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_){
	//hist_->Histogram::WQ2_Fill(envi_, _top, _cut, _W, _Q2, _thrown)
	//hist_->Histogram::Friend_Fill(envi_, _top, _W, _Q2, _MM, _theta, _alpha, _phi , chan, _weight)
}

void Event::Assign_Weight(float weight_){
	_weight = weight_;
}

bool Event::Top(int i){
	return _top[i];
}
bool Event::Gevnt(){
	return _pass; 
}
float Event::Get_Weight(){
	return _weight; 
}
bool Event::Part_Present(int i){
	return _part[i];
}
float Event::Get_P(int i, bool COM_ ){
float _p_ = NAN;
	if(COM_){
		_p_ = _p[i];
	}else{
		_p_ = _p_lab[i];
	}
	return _p_;
}
float Event::Get_Theta(int i, bool COM_ ){
float _theta_ = NAN;
	if(COM_){
		_theta_ = _theta[i];
	}else{
		_theta_ = _theta_lab[i];
	}
	return _theta_;
}
float Event::Get_Phi(int i, bool COM_ ){
	float _phi_ = NAN;
	if(COM_){
		_phi_ = _phi[i];
	}else{
		_phi_ = _phi_lab[i];
	}
	return _phi_;
} 
TLorentzVector Event::Get_Beam(){
	return _k1; 
}
TLorentzVector Event::Get_4Vec(int i, bool COM_ ){
	TLorentzVector _boopers_;
	if(COM_){
		_boopers_ = _vec[i];
	} else{
		_boopers_ = _vec_lab[i];
	}	
	return _boopers_;
}
float Event::Get_W(){
	return _W;
}
float Event::Get_Q2(){
	return _Q2;
}
float Event::Get_MMb(int i){
	return _MMb[i];
}
float Event::Get_Thetab(int i){
	return _thetab[i];
}
float Event::Get_Alphab(int i){
	return _alphab[i];
}

//~Event()




