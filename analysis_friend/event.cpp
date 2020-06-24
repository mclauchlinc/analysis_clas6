#include "event.hpp"

Event::Event(){

}
void Event::COM_4Vec(){
			  //physics::COM_gp(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
	_vec[0] = physics::COM_gp(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
	_vec[1] = physics::COM_gp(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
	_vec[2] = physics::COM_gp(2,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
	_vec[3] = physics::COM_gp(3,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
	_k1 = physics::COM_gp(4,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
	_p1 = physics::COM_gp(5,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
	_COM = true;
}

void Event::Vars(){
	for(int i = 0; i<3; i++){
		if(_COM){
			_alphab[i] = physics::alpha(i, _k1, _vec[0], _vec[1], _vec[2], _vec[3], true);
			_thetab[i] = physics::Ev_Theta(i, _k1, _vec[0], _vec[1], _vec[2], _vec[3], true);
			_MMb[i] = physics::Ev_MM(i, _k1, _vec[0], _vec[1], _vec[2], _vec[3], true);
		}else{
			std::cout<<"You need to get your 4-Momenta into center of mass" <<std::endl;
		}
	}
}

void Event::Fill_Event(std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, int top_, float W_, float Q2_,  Particle p1,  Particle p2,  Particle p3, int hel_){
	if(top_ == 3){
		std::cout<<"Incorrectly Filled Event: Not enough particles for topology" <<std::endl;
	}else{ 
		_W = W_;
		_Q2 = Q2_;
		_hel = hel_;
		_set = p1.Particle::Get_set();
		//std::cout<<"	Fill event p1" <<std::endl;
		//std::cout<<"		Inside Event pt1 Fill for event Combo: " <<std::endl;
		
		//Check to make sure all particles were either simulated or not
		if(p1.Particle::Is_Sim() == p2.Particle::Is_Sim() && p2.Particle::Is_Sim() == p3.Particle::Is_Sim()){
			_sim = p1.Particle::Is_Sim();
			//std::cout<<"	Fill event p2" <<std::endl;
			//std::cout<<"		Inside Event pt2 Fill for event Combo: " <<std::endl;
			
			if(p1.Particle::Is_Thrown() == p2.Particle::Is_Thrown() && p2.Particle::Is_Thrown() == p3.Particle::Is_Thrown()){
				//std::cout<<"	Fill event p3" <<std::endl;
				//std::cout<<"		Inside Event pt3a Fill for event Combo: " <<std::endl;
				
				_thrown = p1.Particle::Is_Thrown();
				//std::cout<<"		Inside Event pt3b Fill for event Combo: " <<std::endl;
				
				_k1 = physics::Set_k_mu(p1.Particle::Get_set());
				//std::cout<<"		Inside Event pt3c Fill for event Combo: " <<std::endl;
				
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
							_MM = physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[2],_vec_lab[3]);
							_MM2 = physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[2],_vec_lab[3]);
							hist_->Histogram::MM_Fill(envi_,0,_MM,0,0,false);
							hist_->Histogram::MM_Fill(envi_,0,_MM,0,1,false);
							_pass = Selection::Event_Selection(top_,_MM);//_k1, _vec_lab[0], _vec_lab[2], _vec_lab[3]);
							p1.Particle::Fill_Par_Event(envi_,hist_,_W,top_,0,_pass); 
							p2.Particle::Fill_Par_Event(envi_,hist_,_W,top_,2,_pass); 
							p3.Particle::Fill_Par_Event(envi_,hist_,_W,top_,3,_pass); 
							if(_pass){
								_vec_lab[1] = _k1 + p_mu - _vec_lab[0] - _vec_lab[2] - _vec_lab[3];
								hist_->Histogram::MM_Fill(envi_,0,_MM,1,0,false);
								hist_->Histogram::MM_Fill(envi_,0,_MM2,1,1,false);
							}else{
								hist_->Histogram::MM_Fill(envi_,0,_MM,2,0,false);
								hist_->Histogram::MM_Fill(envi_,0,_MM2,2,1,false);
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
							_MM = physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[3]);
							_MM2 = physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[3]);
							hist_->Histogram::MM_Fill(envi_,1,_MM,0,0,false);
							hist_->Histogram::MM_Fill(envi_,1,_MM2,0,1,false);
							_pass = Selection::Event_Selection(top_,_MM);//_k1, _vec_lab[0], _vec_lab[1], _vec_lab[3]);
							p1.Particle::Fill_Par_Event(envi_,hist_,_W,top_,0,_pass); 
							p2.Particle::Fill_Par_Event(envi_,hist_,_W,top_,1,_pass); 
							p3.Particle::Fill_Par_Event(envi_,hist_,_W,top_,3,_pass);
							if(_pass){
								_vec_lab[2] = _k1 + p_mu - _vec_lab[0] - _vec_lab[1] - _vec_lab[3]; 
								hist_->Histogram::MM_Fill(envi_,1,_MM,1,0,false);
								hist_->Histogram::MM_Fill(envi_,1,_MM2,1,1,false);
							}else{
								hist_->Histogram::MM_Fill(envi_,1,_MM,2,0,false);
								hist_->Histogram::MM_Fill(envi_,1,_MM2,2,1,false);
							}
						}
					break;
					case 2://PIM missing
						//std::cout<<"	Fill event p4" <<std::endl;
						//std::cout<<"		Inside Event pt1 Fill for event Combo: " <<std::endl;
							
						if(p1.Particle::Is_Elec() && p2.Particle::Is_Pro() && p3.Particle::Is_Pip()){
						//	std::cout<<"		Inside Event pt5a Fill for event Combo: " <<std::endl;
							
							_filled_correctly = true; 
							_p_lab[0] = p1.Particle::Get_p();
							_p_lab[1] = p2.Particle::Get_p();
							_p_lab[2] = p3.Particle::Get_p();
							//std::cout<<"		Inside Event pt5b Fill for event Combo: " <<std::endl;
							
							_theta_lab[0] = p1.Particle::Get_theta();
							_theta_lab[1] = p2.Particle::Get_theta();
							_theta_lab[2] = p3.Particle::Get_theta();
							//std::cout<<"		Inside Event pt5c Fill for event Combo: " <<std::endl;
							
							_phi_lab[0] = p1.Particle::Get_phi();
							_phi_lab[1] = p2.Particle::Get_phi();
							_phi_lab[2] = p3.Particle::Get_phi();
							//std::cout<<"		Inside Event pt5d Fill for event Combo: " <<std::endl;
							
							_vec_lab[0] = physics::Make_4Vector(true,_p_lab[0],_theta_lab[0],_phi_lab[0],me);
							_vec_lab[1] = physics::Make_4Vector(true,_p_lab[1],_theta_lab[1],_phi_lab[1],mp);
							_vec_lab[2] = physics::Make_4Vector(true,_p_lab[2],_theta_lab[2],_phi_lab[2],mpi);
							//std::cout<<"		Inside Event pt5e Fill for event Combo: " <<std::endl;
							
							_top[2] = true; //For some reason this misplaces the index for electrons
							_MM = physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2]);
							_MM2 = physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2]);
							hist_->Histogram::MM_Fill(envi_,2,_MM,0,0,false);
							hist_->Histogram::MM_Fill(envi_,2,_MM2,0,1,false);
							//std::cout<<"		Inside Event pt6 Fill for event Combo: " <<std::endl;
							
							_pass = Selection::Event_Selection(top_,_MM);//_k1, _vec_lab[0], _vec_lab[1], _vec_lab[2]);
							p1.Particle::Fill_Par_Event(envi_,hist_,_W,top_,0,_pass); 
							p2.Particle::Fill_Par_Event(envi_,hist_,_W,top_,1,_pass); 
							p3.Particle::Fill_Par_Event(envi_,hist_,_W,top_,2,_pass);
							//std::cout<<"		Inside Event pt7 Fill for event Combo: " <<std::endl;
							
							if(_pass){
								_vec_lab[3] = _k1 + p_mu - _vec_lab[0] - _vec_lab[1] - _vec_lab[2]; 
								hist_->Histogram::MM_Fill(envi_,2,_MM,1,0,false);
								hist_->Histogram::MM_Fill(envi_,2,_MM2,1,1,false);
								//std::cout<<"		Inside Event pt8 Fill for event Combo: " <<std::endl;
								
							}else{
								hist_->Histogram::MM_Fill(envi_,2,_MM,2,0,false);
								hist_->Histogram::MM_Fill(envi_,2,_MM2,2,1,false);
								//std::cout<<"		Inside Event pt9 Fill for event Combo: " <<std::endl;
								//p1.Particle::Check_Particle();
								//p2.Particle::Check_Particle();
								//p3.Particle::Check_Particle();
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

void Event::Fill_Event(std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, int top_, float W_, float Q2_,  Particle p1,  Particle p2,  Particle p3,  Particle p4, int hel_){
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
					_MM = physics::MM_event(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
					_MM2 = physics::MM_event(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
					hist_->Histogram::MM_Fill(envi_,3,_MM,0,0,false);//The false refers to whether these will be used for fitting, which will only be done on events where there is one of each relevant particle measured
					hist_->Histogram::MM_Fill(envi_,3,_MM2,0,1,false);//The false refers to whether these will be used for fitting, which will only be done on events where there is one of each relevant particle measured
					_pass = Selection::Event_Selection(top_,_MM);//_k1, _vec_lab[0], _vec_lab[1], _vec_lab[2], _vec_lab[3]);
					p1.Particle::Fill_Par_Event(envi_,hist_,_W,top_,0,_pass); 
					p2.Particle::Fill_Par_Event(envi_,hist_,_W,top_,1,_pass); 
					p3.Particle::Fill_Par_Event(envi_,hist_,_W,top_,2,_pass); 
					p4.Particle::Fill_Par_Event(envi_,hist_,_W,top_,3,_pass); 
					if(_pass){
						hist_->Histogram::MM_Fill(envi_,3,_MM,1,0,false);//The false refers to whether these will be used for fitting, which will only be done on events where there is one of each relevant particle measured
						hist_->Histogram::MM_Fill(envi_,3,_MM2,1,1,false);//The false refers to whether these will be used for fitting, which will only be done on events where there is one of each relevant particle measured
						_vec[0] = physics::COM_gp(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
						_vec[1] = physics::COM_gp(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
						_vec[2] = physics::COM_gp(2,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
						_vec[3] = physics::COM_gp(3,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);

					}else{
						hist_->Histogram::MM_Fill(envi_,3,_MM,2,0,false);//The false refers to whether these will be used for fitting, which will only be done on events where there is one of each relevant particle measured
						hist_->Histogram::MM_Fill(envi_,3,_MM2,2,1,false);//The false refers to whether these will be used for fitting, which will only be done on events where there is one of each relevant particle measured
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

void Event::Fill_Event_Hists(std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, bool fit_){
	//std::cout<<"Event 1\n";
	for(int i = 0; i< 4; i++){//Over the different topologies
		if(_top[i] ){//&& _pass){
			//std::cout<<"Event 2\n";
			hist_->Histogram::WQ2_Fill(envi_, i+1, 10, _W, _Q2, _weight, _thrown);
			hist_->Histogram::WQ2_Fill(envi_, 5, 10, _W, _Q2, _weight, _thrown);
			//std::cout<<"Event 3\n";
			if(fit_ && !_thrown){
				//std::cout<<"Event 4\n";
				hist_->Histogram::MM_Fill(envi_,i,_MM,0,0,true);//The false refers to whether these will be used for fitting, which will only be done on events where there is one of each relevant particle measured
				//std::cout<<"Event 5\n";
				hist_->Histogram::MM_Fill(envi_,i,_MM2,0,1,true);
				//std::cout<<"Event 6\n";
				if(_pass){
					//std::cout<<"Event 7\n";
					hist_->Histogram::MM_Fill(envi_,i,_MM,1,0,true);//The false refers to whether these will be used for fitting, which will only be done on events where there is one of each relevant particle measured
					//std::cout<<"Event 8\n";
					hist_->Histogram::MM_Fill(envi_,i,_MM2,1,1,true);
					//std::cout<<"Event 9\n";
				}else{
					//std::cout<<"Event 10\n";
					hist_->Histogram::MM_Fill(envi_,i,_MM,2,0,true);//The false refers to whether these will be used for fitting, which will only be done on events where there is one of each relevant particle measured
					//std::cout<<"Event 11\n";
					hist_->Histogram::MM_Fill(envi_,i,_MM2,2,1,true);
					//std::cout<<"Event 12\n";
				}
			}
			for(int j = 0; j< 3; j++){//For the different kinematic expressions
				//std::cout<<"Event 13\n";
				hist_->Histogram::Friend_Fill(envi_, i+1, _W, _Q2, _MMb[j], _thetab[j], _alphab[j], _phi[j] , j, _weight);
				//std::cout<<"Event 14\n";
			}
		}
	}
}

void Event::Assign_Weight(float weight_){
	_weight = weight_;
}

bool Event::Top(int i){
	return _top[i];
}

int Event::Get_Top(){
	int out = -1; 
	for(int i = 0; i< 4; i++){
		if(_top[i]){
			out = i;
		}
	}
	return out;
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
	if(COM_ && _COM){
		_p_ = _p[i];
	}else{
		_p_ = _p_lab[i];
	}
	return _p_;
}
float Event::Get_Theta(int i, bool COM_ ){
float _theta_ = NAN;
	if(COM_ && _COM){
		_theta_ = _theta[i];
	}else{
		_theta_ = _theta_lab[i];
	}
	return _theta_;
}
float Event::Get_Phi(int i, bool COM_ ){
	float _phi_ = NAN;
	if(COM_ && _COM){
		_phi_ = _phi[i];
	}else{
		_phi_ = _phi_lab[i];
	}
	return _phi_;
} 
float Event::Get_Px(int i, bool COM_ ){
	float px = NAN; 
	if(COM_ && _COM){
		px = _vec[i][0];
	}else{
		px = _vec_lab[i][0];
	}
	return px;
}
float Event::Get_Py(int i, bool COM_ ){
	float py = NAN; 
	if(COM_ && _COM){
		py = _vec[i][1];
	}else{
		py = _vec_lab[i][1];
	}
	return py;
}
float Event::Get_Pz(int i, bool COM_ ){
	float pz = NAN; 
	if(COM_ && _COM){
		pz = _vec[i][2];
	}else{
		pz = _vec_lab[i][2];
	}
	return pz;
}
float Event::Get_P0(int i, bool COM_ ){
	float p0 = NAN; 
	if(COM_ && _COM){
		p0 = _vec[i][3];
	}else{
		p0 = _vec_lab[i][3];
	}
	return p0;
}
TLorentzVector Event::Get_Beam(bool COM_){
	TLorentzVector k1;
	if(COM_ && _COM){
		k1 = _k1;
	}else{
		k1 = _k1_lab;
	}
	return k1; 
}

TLorentzVector Event::Get_Target(bool COM_ ){
	TLorentzVector p1;
	if(COM_ && _COM){
		p1 = _p1;
	}else{
		p1 = p_mu;//Constants
	}
	return p1; 
}

float Event::Get_Beam_Comp(int i, bool COM_ ){
	float p_part = NAN;
	if(COM_ && _COM){
		p_part = _k1[i];
	}else{
		p_part = _k1_lab[i];
	}
	return p_part;
}

float Event::Get_Target_Comp(int i, bool COM_ ){
	float p_part = NAN;
	if(COM_ && _COM){
		p_part = _p1[i];
	}else{
		p_part = p_mu[i];//Constants
	}
	return p_part;
}

TLorentzVector Event::Get_4Vec(int i, bool COM_ ){
	TLorentzVector _boopers_;
	if(COM_ && _COM){
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

bool Event::Get_COM(){
	return _COM;
}

int Event::Get_PID(int i){
	int pid = 0;
	switch(i){
		case 0: pid = ELECTRON; break;
		case 1: pid = PROTON; break;
		case 2: pid = PION; break;
		case 3: pid = -PION; break;
	}
	return pid; 
}

int Event::Get_Set(){
	return _set; //{1,0} -> {e16,e1f}
}

int Event::Get_Hel(){
	return _hel;
}
//~Event()




